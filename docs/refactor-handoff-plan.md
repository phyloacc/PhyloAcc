# PhyloAcc Refactor Handoff Plan

This note summarizes the recommended next refactor sequence after the current
`PhyloAcc-common` extraction work. It is written for maintainers who know the
Python interface well but may not be comfortable with the C++ implementation.

The active C++ programs are:

- `src/PhyloAcc-ST/`: species-tree model executable.
- `src/PhyloAcc-GT/`: gene-tree model executable.
- `src/PhyloAcc-common/`: shared C++ support used by the active ST and GT
  binaries.

Ignore `src/PhyloAcc-ST-GBGC/` for this plan unless project scope changes.

## Current Architecture

At a high level, both ST and GT follow the same broad flow:

1. Parse config.
2. Load alignment/profile and tree inputs.
3. Construct one run-wide `BPP` object.
4. Resolve element IDs.
5. In an OpenMP element loop, construct one per-element `BPP_C` worker.
6. Run M0, M1, and M2.
7. Write rate, tree, likelihood, trace, and status outputs.

The main shared layer already covers the parts that are mostly mechanical:

- `config.*`: config defaults and parser support.
- `profile.*`: FASTA/BED profile loading.
- `newick.*`: species-tree and model-tree loading.
- `run.*`: program kind, model IDs, output bundle, run paths, element IDs,
  and per-element status rows.
- `bpp_constructor.*`: element layout, leaf encoding, missing-column filtering,
  subtree sets, missing-node helpers, and trace buffer allocation.
- `bpp_monitor.*`, `bpp_transition.*`, `bpp_output.*`,
  `bpp_active_output.*`, and `bpp_io.*`: smaller MCMC and output helpers.

The most important remaining ST/GT differences are still in the C++ model
state and sampler code:

- `BPP`: run-wide state. ST and GT share many fields, but GT also carries
  coalescent/theta state, gene-tree state, branch-move state, pi-sampling
  state, tree-output state, and extra likelihood summaries.
- `BPP_C`: per-element worker state. ST and GT share setup concepts, but GT
  owns gene-tree machinery, site/block handling, tree topology sampling, and
  extra traces.
- `src/PhyloAcc-GT/genetree.*`: GT-specific gene-tree model logic.
- `main.cpp`: separate executable orchestration. ST and GT should stay separate
  binaries even if more shared helpers are introduced.

The safe common-code style is therefore small, tested helper functions. Avoid a
large common base class or a single unified `BPP_C` until the model contracts
are much clearer.

## Recommended Order

### 1. Fix correctness hazards first

These should be handled before deeper architecture work because they can make
tests flaky or hide real behavior changes.

#### RNG and shared mutable state

The active seed contract is one public `SEED` parameter. ST and GT derive all
run-wide, worker, site-shuffle, and gene-tree-shuffle streams internally from
that seed plus stable coordinates such as program kind, chain, element ID,
block, and stream label.

This contract is intended to make repeated runs byte-identical for the same
seed, including `NUM_THREAD=1` versus `NUM_THREAD>1`. Deprecated GT parameters
`SEEDS` and `SEED2` are accepted but ignored with a warning.

#### Input and config validation

Known validation issues:

- Python `getOpt` references an undefined `optional` variable during `FILE` and
  `DIR` validation.
- Python `--labelmod` is populated from `args.labeltree` rather than
  `args.labelmod`.
- C++ profile sequence-length validation uses `assert`, which can disappear in
  release builds.
- C++ config parsing currently tolerates unknown keys and loose boolean parsing.

Recommended direction:

1. Fix the small Python parser bugs with focused tests.
2. Replace data/input validation asserts with explicit runtime checks.
3. Decide whether unknown config keys and invalid boolean tokens should become
   fatal errors. For scientific workflows, fail-fast behavior is safer unless a
   compatibility mode is explicitly needed.

#### Narrow RAII and ownership cleanup

There are still raw owning arrays and manual deletes in `BPP` and `BPP_C`.
Do not modernize everything at once.

Recommended first slice:

- Convert worker-local arrays such as `children2`, `parent2`, and `distances2`
  to `std::vector` or small fixed-shape containers where practical.
- Delete copy constructors and copy assignments for classes that cannot be
  safely copied yet.
- Leave the larger `BPP` ownership model for later unless it directly blocks a
  correctness fix.

### 2. Clarify missing-data semantics before changing encoding

Missing-data handling is an unresolved modeling decision, not just a code
cleanup.

Current behavior appears to differ between ST and GT:

- ST uses `MissingBasePolicy::GapOnly`.
- GT uses `MissingBasePolicy::GapNStar`.
- Both lowercase input before encoding.
- Canonical bases map to `0=a`, `1=c`, `2=g`, `3=t`.
- Ambiguity codes such as `r/y/k/m/s/w` map to unknown/nonstandard.
- ST treats only the configured gap character as explicit missing.
- GT treats the gap character, `n`, and `*` as explicit missing.
- Column filtering and per-species missing counts do not obviously use the same
  semantic boundary in all cases.

This may be intentional model behavior. Do not normalize it casually.

Before changing anything, decide and document how each model should treat:

- gaps
- `n`
- `*`
- IUPAC ambiguity codes
- unknown symbols
- per-column filtering
- per-species missing counts
- missing likelihood contributions

Only after that decision should encoding or filtering behavior change. Any
change here can alter scientific output and may require deliberate golden-file
updates.

### 3. Define output and artifact contracts

The output layer is partly shared already, and output writes appear to use
OpenMP critical sections in the most important paths. This makes output-row
builder duplication mostly a maintainability problem today, not the first
correctness blocker.

However, output contracts become critical before introducing model-level
execution.

Questions to answer before "run one model" support:

- What file proves M0 completed successfully?
- What file proves M1 completed successfully?
- What file proves M2 completed successfully?
- Should per-model status files exist, or should `*_elem_status.txt` be the
  single source of truth?
- How should partially failed elements be represented downstream?
- Should Python/Snakemake fail a batch when any element has `error`, or should
  it summarize and continue?
- What existing combined outputs must remain for backward compatibility?

Recommended direction:

1. Extract small row-building helpers only where output format is already
   stable.
2. Keep ST and GT model-specific payloads explicit.
3. Teach Python/post-processing about `*_elem_status.txt`.
4. Design model-level artifacts before changing executable behavior.

### 4. Consolidate tree parsing carefully

`src/PhyloAcc-common/newick.cpp` and `src/PhyloAcc-GT/newick2.cpp` are real
duplication.

The parser internals are copied:

- node and child structs
- custom allocation helpers
- recursive Newick parsing
- leaf counting
- postorder-ish node numbering
- DAG/name/distance extraction

The real differences are loader-level differences:

- `newick.cpp` reads `.mod` files with `BACKGROUND:`, `RATE_MAT:`, and `TREE:`.
- `newick2.cpp` reads a plain coalescent-unit Newick tree for GT.
- `newick.cpp` fills substitution model fields and a `PhyloTree`.
- `newick2.cpp` fills only names, topology, and distances in
  `PhyloTree_theta`.
- GT then derives theta values in `src/PhyloAcc-GT/main.cpp` by comparing the
  species tree to the coalescent tree.

Recommended abstraction:

- One shared Newick parser representation.
- Multiple thin loader adapters:
  - `.mod` species-tree loader.
  - plain Newick coalescent-tree loader.
  - tree-string loader used by GT gene-tree code.

Do not collapse this by forcing all callers into one semantic tree type in the
first patch. The safe abstraction is "one parser, multiple loaders."

### 5. Add tests before tree parser refactors

Tree parser refactors are dangerous because downstream arrays depend on node
IDs and order. The primary risk is not syntax parsing; it is silently changing
node numbering.

Before touching parser behavior, add tests for:

- tip order
- internal node order
- root ID
- `nodes_names` ordering
- branch lengths by node ID
- labeled internal-node behavior
- unlabeled internal-node behavior
- `LoadPhyloTree_text` behavior used by GT gene-tree code
- GT coalescent-tree loading
- GT species-tree/coalescent-tree topology matching

If a refactor changes node ordering, it must be treated as a behavioral change,
not a cosmetic cleanup.

### 6. Then introduce run-one-model functions

Only after the higher-priority contracts are stable should the project start
extracting model-level execution functions.

Recommended staged approach:

1. Add ST-internal helpers that run exactly one model but preserve the current
   M0, M1, M2 order and outputs.
2. Add GT-internal helpers that run exactly one model but preserve the current
   block iteration, tree output, and outputs.
3. Keep ST and GT helpers separate at first.
4. Add shared orchestration only around stable concepts such as `ModelId`,
   `ModelSpec`, output bundle access, and status reporting.
5. Only then expose selected-model or artifact-reuse behavior through config,
   CLI, and Python/Snakemake.

This sequencing avoids building model-level APIs on top of unresolved RNG,
validation, missing-data, and output-contract behavior.

## Test Coverage Needed

Current tests are useful but not sufficient for the next phase.

Minimum additions before major refactors:

- RNG reproducibility tests for repeated ST runs.
- RNG reproducibility tests for repeated GT runs.
- A GT test that exercises `NUM_THREAD > 1`.
- Tests for seed derivation.
- Tests for Python parser path validation and `--labelmod`.
- C++ profile-loader tests for malformed row lengths.
- Config parser tests for unknown keys and invalid booleans if behavior is
  changed.
- Missing-data encoding tests after a policy decision is made.
- Tree parser ordering and topology tests before consolidating `newick2`.
- Output/status tests that assert how errors, filtered elements, and model
  failures are represented.

Golden outputs should not be updated casually. If a golden changes, document:

- which behavior changed,
- why the change is intended,
- whether it is numerical, formatting-only, ordering-related, or a true model
  semantics change.

## Practical Sequence

A pragmatic next sequence would be:

1. Add test coverage around RNG determinism and tree parser ordering.
2. Fix RNG/shared mutable state across ST and GT.
3. Fix Python path parsing, `--labelmod`, and C++ profile validation.
4. Add narrow RAII protections around worker-owned state.
5. Decide and document missing-data semantics.
6. Define model artifact and status contracts.
7. Consolidate tree parsing behind one parser with separate loaders.
8. Extract ST and GT run-one-model helpers separately.
9. Introduce shared model orchestration only after the separate helpers are
   behaviorally locked down.

The guiding rule is: abstract stable mechanics, not unresolved model semantics.
