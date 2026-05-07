# Duplicated Code Refactor Handoff

This report summarizes a review-only pass over duplicated code in the PhyloAcc repository. It is intended as a handoff for an agentic model or developer to continue from the review without repeating the initial survey.

No code changes were made during the review pass. This report is the only artifact produced from the findings.

## Review Context

- Local checkout reviewed: `/Users/tim/Work/PhyloAcc`
- Review type: duplicated-code discovery only
- Network access: not used
- Sub-agents: not used
- Environment note: `conda run -n phyloacc-test python --version` failed with `NoWritableEnvsDirError`, so the review used only direct read-only file inspection commands.

The codebase has two main duplication zones:

1. C++ model implementations, especially divergence between `src/PhyloAcc-ST`, `src/PhyloAcc-GT`, `src/PhyloAcc-ST-GBGC/SRC`, and `src/PhyloAcc-common`.
2. Python interface support code, especially active class-based tree logic versus legacy function-based tree logic, plus repeated run/report scaffolding.

## High-Level Refactor Order

Recommended order if the goal is maximum duplication reduction with manageable risk:

1. Move `ST-GBGC` onto the shared C++ runtime scaffolding in `src/PhyloAcc-common`.
2. Consolidate Newick parsing and theta/non-theta tree loading.
3. Extract duplicated C++ log-space math helpers into a common header.
4. Collapse Python `tree.py`/`tree_old.py` behavior behind one tree API or a compatibility adapter.
5. Remove or wrap legacy sCF implementations in `tree_old.py`.
6. Deduplicate Python params/option lifecycle helpers.
7. Deduplicate Python plot and HTML rendering primitives.

Run focused tests after each slice. Avoid batching all of these into one patch because the C++ runtime and Python tree/sCF paths have different risk profiles.

## Findings

### 1. ST-GBGC is still a fork of the shared C++ runtime

- Priority: P2
- Primary file: `src/PhyloAcc-ST-GBGC/SRC/main.cpp`
- Starting line: 84

`ST` and `GT` now route config/loading/output through `PhyloAcc-common`, but `ST-GBGC` still carries its own `LoadParams`, `DirectoryExists`, and run setup. This duplicates behavior already represented in `src/PhyloAcc-common/config.cpp` and `src/PhyloAcc-common/run.cpp`, and it makes option semantics likely to drift.

Suggested direction:

- Add a `ProgramKind::ST_GBGC` or equivalent shared config path.
- Move gBGC-specific options into `phyloacc::Config` with clear defaults.
- Reuse shared helpers for:
  - config loading
  - output directory validation
  - profile loading
  - species tree loading
  - run summary output
  - output bundle/path construction where possible

Useful nearby files:

- `src/PhyloAcc-common/config.cpp`
- `src/PhyloAcc-common/config.h`
- `src/PhyloAcc-common/run.cpp`
- `src/PhyloAcc-common/run.h`
- `src/PhyloAcc-ST/main.cpp`
- `src/PhyloAcc-GT/main.cpp`

### 2. GBGC BPP methods duplicate existing common helpers

- Priority: P2
- Primary file: `src/PhyloAcc-ST-GBGC/SRC/bpp.cpp`
- Starting line: 40

`ST-GBGC/SRC/bpp.cpp` locally implements tree array setup, profile/tree matching, MCMC storage initialization, proposal sampling, likelihood traversal, tree traversal, and init output. The active `ST` and `GT` implementations already delegate much of this work to `PhyloAcc-common`.

Existing common helpers to reuse:

- `phyloacc::InitializeTreeArrays`
- `phyloacc::MatchProfileToTree`
- `phyloacc::InitializeCommonMCMCStorage`
- `phyloacc::SampleProposal`
- `phyloacc::ComputeLogLikelihood`
- `phyloacc::CollectUpperTreeNodes`
- `phyloacc::CollectSubtreeNodes`
- `phyloacc::CollectSubtreeNodesUntilChildren`
- shared output helpers in `bpp_io.*`, `bpp_output.*`, and `bpp_active_output.*`

Suggested direction:

- First convert exact shared behavior in `ST-GBGC/SRC/bpp.cpp` to common helpers.
- Keep gBGC-specific behavior as explicit extension code rather than mixing it into copied general logic.
- After each conversion, compare output against existing ST-GBGC behavior if a fixture exists. If no fixture exists, add a minimal one before the migration.

### 3. Newick parsing exists in three near-copies

- Priority: P2
- Primary file: `src/PhyloAcc-GT/newick2.cpp`
- Starting line: 36

`src/PhyloAcc-common/newick.cpp`, `src/PhyloAcc-GT/newick2.cpp`, and `src/PhyloAcc-ST-GBGC/SRC/newick.cpp` all implement the same custom allocator, parser, traversal, and loader structure. Differences are mostly suffix/name/type changes and theta handling.

Suggested direction:

- Create one shared Newick parser representation in `PhyloAcc-common`.
- Support theta and non-theta extraction through loader/output policy functions rather than separate copied parsers.
- Replace `newick2.cpp` with a thin theta loader over the shared parser.
- Replace `ST-GBGC/SRC/newick.cpp` with the shared loader once `ST-GBGC` is using the common include path.

Important risk:

- Tree node ordering and label assignment are probably observable in output and tests. Preserve ordering exactly, or update tests intentionally with clear provenance.

Useful files:

- `src/PhyloAcc-common/newick.cpp`
- `src/PhyloAcc-common/newick.h`
- `src/PhyloAcc-GT/newick2.cpp`
- `src/PhyloAcc-GT/newick2.h`
- `src/PhyloAcc-ST-GBGC/SRC/newick.cpp`
- `src/PhyloAcc-ST-GBGC/SRC/newick.h`

### 4. Log-space math helpers are copied across all BPP headers

- Priority: P2
- Primary file: `src/PhyloAcc-ST/bpp.hpp`
- Starting line: 358

The following helpers are repeated in `ST`, `GT`, and `ST-GBGC` headers:

- `log_exp_multi`
- `log_multi`
- `log_multi2` or near equivalents
- `log_exp_colsum`
- `log_exp_sum`
- `log_sample`
- `log_sample_norm`
- `printZ`

Suggested direction:

- Add a small common header such as `src/PhyloAcc-common/bpp_math.h` or `src/PhyloAcc-common/log_math.h`.
- Move identical log-space functions there as inline functions in the `phyloacc` namespace.
- Keep compatibility wrappers in `BPP` only if too many call sites currently use `BPP::log_multi`.
- Normalize pass-by-value versus pass-by-reference differences carefully. Some current functions mutate their arguments, and some copies differ in signature.

Important risk:

- These helpers affect numerical stability and sampling probabilities. Use focused unit tests before and after extraction.

### 5. BPP_C tree/message updates share the same skeleton

- Priority: P2
- Primary file: `src/PhyloAcc-ST-GBGC/SRC/bpp_c.cpp`
- Starting line: 453

`getSubtree`, `getSubtree_missing`, changed-Z ancestor invalidation, `Update_Tg`, and `MonitorChain` recur across `ST`, `GT`, and `ST-GBGC`. Some ancestor/update code is already shared for `ST` and `GT`, but `ST-GBGC` still has local copies and model-specific variants.

Existing common helpers:

- `phyloacc::CollectSubtreeNodes`
- `phyloacc::MarkChangedZAncestors`
- `phyloacc::CopyTraceZ`
- `phyloacc::ComputeBaseFullLogLik`
- `phyloacc::UpdateMaxState`
- `phyloacc::SampleTransitionRates`

Suggested direction:

- First replace direct duplicates with existing helpers.
- For `Update_Tg`, consider a helper that owns traversal and message accumulation while accepting a transition-matrix selector callback or policy object.
- For `MonitorChain`, split common trace/max-state bookkeeping from model-specific extra terms such as gBGC trace fields or gene-tree recording.

Important risk:

- This is behaviorally sensitive MCMC code. Prefer small mechanical extractions with golden or deterministic tests after each step.

### 6. Python tree logic is split between class and legacy function APIs

- Priority: P2
- Primary file: `src/PhyloAcc-interface/phyloacc_lib/tree.py`
- Starting line: 458

`tree.py` and `tree_old.py` both implement branch categorization, clade/split handling, quartet sampling, branch-length rewriting, and tree labeling paths. The interface still has branching support for `tree-data-type`, but the default is class-based.

Active references to legacy `tree_old.py` remain in:

- `src/PhyloAcc-interface/phyloacc.py`
- `src/PhyloAcc-interface/phyloacc_lib/treeio.py`
- `src/PhyloAcc-interface/phyloacc_lib/cf.py`
- `src/PhyloAcc-interface/phyloacc_lib/plot.py`
- `src/PhyloAcc-interface/phyloacc_lib/batch.py`

Suggested direction:

- Decide whether the function-based tree API is still a supported mode.
- If not supported, remove `tree-data-type == "func"` branches and delete legacy implementations after tests pass.
- If still supported, build a compatibility adapter around `Tree` so callers use one semantic implementation.
- Preserve output labels and branch ordering in `treeio.writeCF` and plotting paths.

Useful tests:

- `tests/test_unit_tree_groups.py`
- `tests/test_unit_scf.py`
- interface tests that exercise plotting or tree labeling

### 7. sCF calculation has an active copy and a legacy copy

- Priority: P2
- Primary file: `src/PhyloAcc-interface/phyloacc_lib/tree_old.py`
- Starting line: 416

`src/PhyloAcc-interface/phyloacc_lib/cf.py` contains the current tested `locusSCF` and `scf` implementation. `tree_old.py` still carries older versions of the same algorithms.

Suggested direction:

- Remove `tree_old.locusSCF` and `tree_old.scf` if nothing active calls them.
- If compatibility is required, replace them with thin wrappers around `cf.py`.
- Keep all quartet scoring, skip-character behavior, and loop/zip site mode in one implementation.

Useful tests:

- `tests/test_unit_scf.py`

### 8. Python params and option parsing duplicate lifecycle scaffolding

- Priority: P3
- Primary file: `src/PhyloAcc-interface/phyloacc_lib/post_params.py`
- Starting line: 15

`params.py` and `post_params.py` duplicate `StrictDict`, metadata initialization, `info.yaml` loading, common runtime flags, log filename construction, and output/status keys. `opt_parse.py` and `post_opt_parse.py` also duplicate output directory creation, log setup, quiet/norun handling, and start-program logging.

Suggested direction:

- Create a shared `params_base.py` or similar module for:
  - `StrictDict`
  - start time/date metadata
  - Python version
  - `info.yaml` loading
  - common flags such as `quiet`, `norun`, `debug`, `nolog`, `overwrite`
- Create option/output helpers for:
  - path checks
  - output directory creation
  - log initialization
  - common start-program header lines
- Keep pre-run and post-run option differences declarative in their existing modules.

Risk:

- Low to medium. This is not numerical core code, but it can alter CLI behavior or error messages.

### 9. Plot and HTML generation repeat small rendering primitives

- Priority: P3
- Primary file: `src/PhyloAcc-interface/phyloacc_lib/plot.py`
- Starting line: 29

`genPlots` and `genPlotsPost` both set the same matplotlib theme and repeatedly hand-code histogram/scatter/bar save patterns. `writeHTML` and `writeHTMLPost` share the same template-fill lifecycle.

Suggested direction:

- Extract small helpers such as:
  - `apply_plot_theme`
  - `save_hist`
  - `save_scatter`
  - `save_bar`
  - `write_template_html`
- Keep report-specific labels and template arguments in the existing pre/post functions.

Risk:

- Low. This is mostly presentation code, but visual output and generated filenames should be checked.

## Additional Observations

The top-level `Makefile` already documents that shared implementations live in `src/PhyloAcc-common/` and that legacy local copies under `ST`/`GT` are excluded from active builds. That migration appears partially complete for `ST` and `GT`, while `ST-GBGC` remains a separate fork.

The tests documentation also acknowledges duplicate legacy `profile.*`, `newick.*`, and `utils.*` files for `ST`/`GT`. Do not treat those excluded legacy files as active implementation targets unless the build system changes.

## Validation Suggestions

After C++ common-helper refactors:

- Run C++ unit tests.
- Run ST minimal golden tests.
- Run GT golden tests when available/enabled.
- For ST-GBGC, add or run a dedicated fixture before migrating behavior if one is not already present.

After Python tree/sCF refactors:

- Run tree group tests.
- Run sCF unit tests.
- Run interface tests that cover tree reading, group parsing, output writing, plotting, and HTML generation.

After plotting/HTML refactors:

- Verify expected plot files are created.
- Verify HTML files render expected paths and summary values.

## Suggested Implementation Slices

### Slice A: Low-risk Python cleanup

1. Remove or wrap `tree_old.locusSCF` and `tree_old.scf`.
2. Confirm all sCF calls use `cf.py`.
3. Run `tests/test_unit_scf.py`.

### Slice B: C++ log math extraction

1. Add common log/math header.
2. Replace duplicated static methods with wrappers or direct calls.
3. Add narrow unit coverage for log-sum-exp helpers.
4. Run C++ unit tests and minimal ST/GT tests.

### Slice C: ST-GBGC shared runtime migration

1. Add shared config support for gBGC.
2. Move `main.cpp` onto `LoadConfigForProgram`-style flow.
3. Replace `BPP` helpers with common helper calls.
4. Add fixture or golden coverage before changing MCMC internals.

### Slice D: Newick parser consolidation

1. Write tests that lock current node ordering, labels, branch lengths, theta values, and topology checks.
2. Introduce shared parser representation.
3. Convert non-theta loader.
4. Convert theta loader.
5. Remove local parser forks after tests pass.

## Handoff Summary

The most meaningful duplicated-code refactor is not a cosmetic cleanup; it is completing the migration to `src/PhyloAcc-common`. `ST` and `GT` are already partly modernized, which gives clear examples for how `ST-GBGC` should be moved. The highest-risk areas are Newick parsing and MCMC message/trace code because they can change output semantics or numeric behavior. The lowest-risk starting points are Python sCF legacy wrappers and C++ log math extraction with focused tests.
