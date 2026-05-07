# Code Issues Handoff

This document tracks resolved, unfixed, and partially unresolved code-review
findings in the PhyloAcc repository. It is intended as a resume guide for a
future coding agent.

Current baseline commits already landed before this handoff:

- `a467add Share PhyloAcc run setup`
- `21d06c7 Record per-element run status`
- `56b6990 Add BPP constructor helpers`
- `a44232a Share BPP constructor setup`
- `560923e Add BPP_C leaf encoding helpers`
- `31227c3 Share BPP_C leaf encoding`
- `c7bbdc0 Share BPP_C post-encoding setup`
- `1b4bbe8 Implement one public RNG seed contract`

The last full verification after the one-public-seed RNG work passed:

- ST, GT, and C++ test binaries rebuilt with the project-local Pixi toolchain.
- `tests/test_cpp_unit.py`: `13 passed`
- `tests/test_rng_repro.py` with GT/RNG enabled: `6 passed`
- Full enabled suite with GT, optional test-data, and RNG tiers:
  `42 passed in 631.60s`

## Recently Resolved High Priority

### P1: GT shares one mutable RNG across workers

Status: fixed in the current refactor baseline and verified after `1b4bbe8`.

Primary locations:

- `src/PhyloAcc-GT/bpp_c.hpp`: each `BPP_C` worker allocates and owns its
  GSL RNG.
- `src/PhyloAcc-GT/bpp_c.hpp`: worker GSL, site shuffle, and gene-tree shuffle
  streams are derived from one public `SEED` plus stable stream coordinates.
- `src/PhyloAcc-GT/bpp.hpp`: the run-wide GSL RNG is also seeded through the
  derived `RunGsl` stream.
- `src/PhyloAcc-common/rng.*`: shared seed derivation helper and stream labels.

Resolution:

- GT workers no longer alias the parent `BPP::RNG`.
- `SEED` is the only active public seed. Deprecated GT `SEEDS` and `SEED2`
  are accepted but ignored with a warning.
- Repeated GT runs are covered for same-thread and cross-thread reproducibility.

Verification:

- C++ RNG tests cover stable derivation, stream uniqueness, `RunGsl`, and
  stable `MakeTwister` output.
- C++ config parser tests cover `SEED`, ignored `SEEDS`, and ignored `SEED2`.
- `tests/test_rng_repro.py` covers GT same-thread and `NUM_THREAD=1` versus
  `NUM_THREAD=2` reproducibility.

### P1: ST workers reuse the same seed and still call global `rand()`

Status: fixed in the current refactor baseline and verified after `1b4bbe8`.

Primary locations:

- `src/PhyloAcc-ST/bpp_c.hpp`: each `BPP_C` worker allocates and owns its
  GSL RNG.
- `src/PhyloAcc-ST/bpp_c.hpp`: worker GSL streams are derived from one public
  `SEED` plus stable stream coordinates.
- `src/PhyloAcc-ST/bpp.hpp`: the run-wide GSL RNG is seeded through the
  derived `RunGsl` stream.

Resolution:

- ST workers no longer start from the same GSL stream.
- Active ST sampling no longer depends on process-global `rand()`/`srand()`.
- Repeated ST runs are covered for same-thread and cross-thread reproducibility.

Verification:

- C++ RNG tests cover stable derivation, stream uniqueness, `RunGsl`, and
  stable `MakeTwister` output.
- `tests/test_rng_repro.py` covers ST same-thread and `NUM_THREAD=1` versus
  `NUM_THREAD=2` reproducibility.

## Correctness And Validation

### P2: Python option path validation references an undefined optional flag

Status: unfixed.

Primary location:

- `src/PhyloAcc-interface/phyloacc_lib/opt_parse.py`: `getOpt` checks `if not optional` for `FILE` and `DIR` validation, but `optional` is not defined.

Problem:

Bad required paths can raise `NameError` instead of the intended user-facing validation error. This makes invalid CLI/config input fail unclearly.

Recommended fix direction:

- Define the optional/required contract explicitly in `getOpt`, or pass it as a parameter if the caller needs to choose.
- Add focused Python tests for missing required file and directory options.

### P2: `--labelmod` parses the wrong argument

Status: unfixed.

Primary location:

- `src/PhyloAcc-interface/phyloacc_lib/opt_parse.py`: `globs['label-mod'] = getOpt(args.labeltree, "labelmod", ...)`.

Problem:

The `label-mod` option is populated from `args.labeltree` instead of `args.labelmod`. This means `--labelmod` can be ignored, while `--labeltree` can unintentionally toggle label-mod behavior.

Recommended fix direction:

- Change the source argument to `args.labelmod`.
- Add parser tests for `--labeltree`, `--labelmod`, and the two options together.

### P2: Runtime input validation depends on `assert`

Status: unfixed.

Primary location:

- `src/PhyloAcc-common/profile.cpp`: alignment row length is checked with `assert(wholeline.length() == prof.G)`.

Problem:

Assertions disappear from release builds compiled with `NDEBUG`. Malformed alignment lengths should fail explicitly in all builds, because this is input validation rather than an internal invariant.

Recommended fix direction:

- Replace the assertion with an explicit runtime check that reports the species/line and expected/observed lengths.
- Apply the same policy to other input/data validation asserts when encountered.

Testing needed:

- Add a profile-loading test with one malformed sequence row.
- Verify the error is visible and deterministic in both debug and release-like builds.

### P3: Config typos continue with defaults

Status: unfixed.

Primary locations:

- `src/PhyloAcc-common/config.cpp`: unknown parameters are printed and skipped.
- `src/PhyloAcc-common/config.cpp`: `StringToBool` returns `false` for anything except `true` or `1`.

Problem:

Unknown config keys and misspelled boolean values can silently change an analysis by leaving defaults in place. For scientific batch runs, config typos should fail fast unless there is an explicit compatibility mode.

Recommended fix direction:

- Turn unknown parameter keys into fatal config errors.
- Make boolean parsing accept documented values only and fail on invalid tokens.
- Consider a temporary compatibility flag only if old config files rely on loose parsing.

Testing needed:

- Unit-test unknown key handling.
- Unit-test valid and invalid boolean values.

## Architecture Blockers

### P2: Model runs are still hard-wired together

Status: unfixed. The shared run/config/output structures make this easier to address, but they did not decouple execution.

Primary locations:

- `src/PhyloAcc-ST/main.cpp`: each element still runs M0, M1, and M2 in sequence in one worker.
- `src/PhyloAcc-GT/main.cpp`: each element still runs null, restricted, and full models in sequence in one worker.
- `src/PhyloAcc-common/run_common.*`: current `ModelSpec` mapping and output bundle are useful foundations, but not a model artifact boundary.

Problem:

The current execution model recomputes null/full model work when users want to compare multiple target groups over the same alignment/tree. There is no way to run or reuse only M0, M1, or M2 artifacts.

Recommended fix direction:

- Introduce a model-level execution unit around `ModelId`.
- Define model-level artifacts before changing the Python/Snakemake interface.
- Preserve current output filenames and combined outputs until a migration plan exists.
- Add resume/reuse semantics only after model artifacts can be validated independently.

Suggested staged plan:

1. Add internal model-run functions that execute exactly one model but keep the same entrypoint order.
2. Write per-model status and metadata internally, while still producing current combined outputs.
3. Add config/CLI controls for selected models and artifact reuse.
4. Update the Python/Snakemake interface after the C++ model artifact contract is stable.

### P2: Generated workflow tracks only combined batch completion

Status: unfixed.

Primary location:

- `src/PhyloAcc-interface/phyloacc_lib/templates.py`: `rule all` checks batch completion through `*_elem_lik.txt`.

Problem:

The interface treats a PhyloAcc batch as complete when the combined likelihood file exists. It cannot skip or reuse only M0, M1, or M2 outputs.

Recommended fix direction:

- Do not change this before the C++ model-level artifacts are designed.
- Once model artifacts exist, expose separate Snakemake outputs for M0/M1/M2.
- Keep a backward-compatible combined output target for existing users.

### P2: Element status exists, but downstream policy is still incomplete

Status: partially fixed.

Primary locations:

- `src/PhyloAcc-common/run_common.*`: writes `*_elem_status.txt`.
- `src/PhyloAcc-ST/main.cpp` and `src/PhyloAcc-GT/main.cpp`: write `ok`, `filtered`, `model_failure`, or `error` rows.

Problem:

The silent per-element exception problem has been addressed at the C++ output level. However, the Python interface and post-processing still mostly reason from combined result files. Failed elements may need more visible summary reporting in user-facing workflows.

Recommended fix direction:

- Decide whether the interface should fail a batch when any element has `error`, or merely summarize and propagate failed element ids.
- Teach post-processing to read `*_elem_status.txt`.
- Add a concise run-level summary such as counts by status.

Testing needed:

- Add a test fixture with a forced element error once there is a stable way to trigger one.
- Verify status output does not break current post-processing.

## Ownership And Memory Safety

### P2: Raw owning arrays leak or copy unsafely

Status: unfixed.

Primary locations:

- `src/PhyloAcc-ST/bpp_c.hpp`: `children2`, `parent2`, and `distances2` are raw owning arrays.
- `src/PhyloAcc-GT/bpp_c.hpp`: `children2` and `gtree` use manual ownership.
- `src/PhyloAcc-ST/bpp.hpp` and `src/PhyloAcc-GT/bpp.hpp`: BPP objects still carry raw tree and model arrays.

Problem:

Some owned arrays are not cleaned up consistently, and copy operations are not explicitly deleted or made safe. That creates leak and accidental shallow-copy risk. It also makes larger refactors riskier because ownership is implicit.

Recommended fix direction:

- Convert narrow worker-local arrays first, preferably to `std::vector` or `std::array<int, 2>` where shapes are fixed.
- Delete copy constructors/assignments for classes that cannot yet be safely copied.
- Move toward RAII before deeper model-decoupling refactors.

Testing needed:

- Run the existing integration suite after each small ownership patch.
- If available, add an AddressSanitizer build target later; do not block small RAII patches on that infrastructure.

## ST/GT Semantics That Need Explicit Decisions

### P2: Missing-data semantics differ between ST and GT

Status: unfixed as a design question. Recent constructor helper tests preserve the current behavior.

Primary locations:

- `src/PhyloAcc-common/bpp_constructor.cpp`: shared leaf encoding and filtering helpers.
- `src/PhyloAcc-ST/bpp_c.hpp`: ST uses `MissingBasePolicy::GapOnly`.
- `src/PhyloAcc-GT/bpp_c.hpp`: GT uses `MissingBasePolicy::GapNStar`.
- `src/PhyloAcc-common/config.cpp`: ST and GT have different default `revgap` values.

Current behavior:

- Both ST and GT lowercase profile input before encoding.
- Canonical bases map to `0=a`, `1=c`, `2=g`, `3=t`.
- Ambiguity codes such as `r/y/k/m/s/w` map to `-1`.
- ST treats only the configured gap character as explicit missing (`Tg == 4`); `n` and `*` become unknown/nonstandard (`Tg == 5`).
- GT treats the gap character, `n`, and `*` as explicit missing (`Tg == 4`).
- Column filtering counts `Tg >= 4`, but per-species missing counts currently count only `Tg == 4`.
- ST missing likelihood uses `consToMis`/`nconsToMis` directly in `log_emission`.
- GT missing likelihood passes through gene-tree likelihood logic and uses different defaults.

Potential issue:

This may be intentional model behavior, but it is not documented as an explicit contract. If accidental, ST and GT are not analyzing `n`/`*` in the same way.

Recommended fix direction:

- Ask for a scientific/modeling decision before changing behavior.
- Document the intended treatment for gaps, `n`, `*`, ambiguity codes, and unknown symbols.
- Only then update encoding, filtering, and tests.

Testing needed:

- Extend `tests/cpp/test_bpp_constructor_main.cpp` to cover the final intended policy.
- Add at least one ST/GT integration case containing `n`, `*`, gaps, and ambiguity codes.

### P2: GT site-ordering behavior around simple/missing blocks is surprising

Status: unfixed and intentionally preserved during dedup.

Primary location:

- `src/PhyloAcc-GT/bpp_c.hpp`: GT shuffles sites with a derived site-shuffle RNG, then reorders simple or missing-like blocks and preserves `idblk_count`.

Problem:

The current logic derives a classification order from shuffled sites, but the final helper call uses the resulting positions as direct site indices. That behavior was preserved to avoid changing numerical output, but it deserves a focused audit.

Recommended fix direction:

- Add a small targeted test that locks down current ordering.
- Decide whether the current behavior is intended or a historical bug.
- If changing it, expect golden output changes and document them explicitly.

## Remaining Duplication

### P2: MCMC kernels have shared skeletons with model-specific payloads

Status: unfixed. Constructor setup was deduplicated, but core MCMC logic remains mostly separate.

Primary locations:

- `src/PhyloAcc-ST/bpp_c.cpp`: `MonitorChain`, `sample_rate`, and sampling helpers.
- `src/PhyloAcc-GT/bpp_c.cpp`: `MonitorChain`, `sample_rate`, and sampling helpers.
- `src/PhyloAcc-ST/bpp_c2.cpp`: `Output_init` and trace output.
- `src/PhyloAcc-GT/bpp_c2.cpp`: `Output_init`, tree output, and trace output.
- `src/PhyloAcc-common/bpp_monitor.*`: existing shared monitor helpers.

Current shape:

- `sample_transition` is already partly shared through common helper code.
- `MonitorChain` has a shared skeleton: copy traces, compute full likelihood, update maxima, emit optional diagnostics.
- GT adds pi priors, gene-tree topology records, `Max_pi`, and `Max_GT`.
- ST has a simpler species-tree likelihood path and simpler max-state records.
- `Output_init` and trace output share row/header structure but GT has extra tree/pi/indicator payloads.
- `sample_rate` has a common Metropolis-Hastings outline, but the payloads diverge substantially.

Recommended fix direction:

- Deduplicate output and monitor formatting before touching rate kernels.
- Extract small common row-building/data-reduction helpers first.
- Leave `sample_rate` until after narrow RAII fixes and more focused tests,
  because it mutates a lot of state and differs more deeply between ST and GT.

Testing needed:

- Keep golden output comparisons byte-identical for output-only refactors.
- Add unit tests around any common row builders before wiring them into ST/GT.

### P2: Newick parsing is duplicated

Status: unfixed.

Primary locations:

- `src/PhyloAcc-GT/newick2.cpp`
- `src/PhyloAcc-common/newick.*`

Problem:

GT carries a separate Newick parser that appears to duplicate common tree parsing and memory management. The common `PhyloTree` already has theta-related storage, so this may be collapsible into one parser API.

Recommended fix direction:

- First identify the exact GT-only fields or parse behavior that `newick2.cpp` still provides.
- Add tests covering theta/coalescent tree parsing before deleting or merging parsers.
- Collapse into common parser code only after the test contract is clear.

### P3: Interface carries two tree APIs

Status: unfixed.

Primary locations:

- `src/PhyloAcc-interface/phyloacc_lib/treeio.py`
- `src/PhyloAcc-interface/phyloacc.py`
- `src/PhyloAcc-interface/phyloacc_lib/cf.py`
- `src/PhyloAcc-interface/phyloacc_lib/plot.py`
- `src/PhyloAcc-interface/phyloacc_lib/html.py`

Problem:

The Python interface still branches between the class-based `Tree` API and legacy `tree_old` function-style API. The default is class-based, but the old branch shape is repeated across the interface.

Recommended fix direction:

- Confirm whether `tree-data-type = func` is still supported.
- If not supported, remove old branches in one focused interface cleanup.
- If still supported, isolate the compatibility layer so downstream code sees one tree API.

### P3: Legacy ST-GBGC tree duplicates active code

Status: unfixed.

Primary location:

- `src/PhyloAcc-ST-GBGC/`

Problem:

The ST-GBGC directory contains its own BPP/profile/newick/utils sources and build path. It duplicates active ST/common code and can mislead maintainers into thinking fixes were applied everywhere when they only touched active ST/GT.

Recommended fix direction:

- Decide whether ST-GBGC is supported.
- If supported, make it share common code or document why it cannot.
- If unsupported, mark it clearly as archived/legacy so future agents do not refactor it accidentally as active code.

### P3: HTML templates duplicate layout

Status: unfixed.

Primary location:

- `src/PhyloAcc-interface/phyloacc_lib/templates.py`

Problem:

Pre-run and post-run HTML templates duplicate page skeleton, navigation, run info, and styling assumptions. This is lower risk than C++ correctness work but will slow interface cleanup.

Recommended fix direction:

- Extract shared template fragments after higher-priority CLI/workflow issues are fixed.
- Add snapshot-style tests if practical, or at least preserve generated output for a known small config.

## Suggested Order Of Attack

1. Fix small parser/input-validation bugs: undefined `optional`, wrong `--labelmod` source, and runtime alignment validation.
2. Start RAII cleanup for narrow `BPP_C` worker-owned arrays and delete unsafe copy operations.
3. Decide and document missing-data semantics before changing any encoding behavior.
4. Introduce internal single-model execution functions as the first step toward model decoupling.
5. Teach interface/post-processing about per-element status and later model-level artifacts.
6. Deduplicate monitor/output formatting, then revisit heavier MCMC kernel duplication.
7. Clean up lower-risk interface/template/parser duplication once numerical behavior is protected by tests.

## General Test Guidance

Use the current baseline commands after production C++ changes:

```bash
pixi run build
pixi run cpp-tests
pixi run test-cpp
pixi run test-all
git diff --check
```

For small Python interface fixes, run the relevant Python tests plus `pixi run test-all` before committing.

For numerical or RNG changes, define the reproducibility contract first, then add tests that compare repeated runs with the same seed. Do not update golden outputs casually; if a golden output changes, explain why the numerical behavior changed and whether that change is intended.

## Already Addressed Or Partially Addressed

These older findings should not be restarted from scratch:

- The active ST/GT RNG contract now uses one public `SEED` with derived
  run-wide, worker, site-shuffle, and gene-tree-shuffle streams.
- Shared ST/GT run setup and output stream wiring were introduced in `src/PhyloAcc-common/run_common.*`.
- BPP and BPP_C constructor setup was partially deduplicated through `src/PhyloAcc-common/bpp_constructor.*`.
- Per-element status output now exists as `*_elem_status.txt`.

These changes are foundations, not complete architectural cleanup. Future agents should build on them rather than replacing them wholesale.
