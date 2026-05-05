Test Suite (Pytest)

Purpose
- Protect refactors by locking down Python logic, C++ parsing/invariants, and end-to-end ST/GT behavior.
- Keep fast deterministic checks separate from heavier integration runs.
- Use golden files to catch numeric drift.

Quick Start
- Build binaries: `make`
- Build C++ unit binaries: `make cpp-tests`
- Run default suite: `pytest -q`
- Run GT and optional external-data coverage: `PHYLOACC_RUN_GT=1 PHYLOACC_RUN_TESTDATA=1 pytest -q`

Environment
- This repo expects the conda environment `phyloacc-test`.
- The C++ unit binaries in `tests/cpp/` link against libraries from `CONDA_PREFIX`.
- Regression validation for ST/GT binaries assumes the supported GCC14 wrapper/toolchain path, not an arbitrary local `g++`.
- Example supported rebuild:
  - `conda run -n phyloacc-test make -B CXX=/n/holylfs05/LABS/informatics/Lab/projects/gwct/phyloacc/PhyloAcc/dev/cc14-wrap.sh PREFIX=/n/home07/gthomas/miniconda3/envs/phyloacc-test PhyloAcc-ST PhyloAcc-GT`
- The duplicate `profile.*`, `newick.*`, and `utils.*` files still present under `src/PhyloAcc-ST/` and `src/PhyloAcc-GT/` are legacy copies excluded by `Makefile`; the active shared implementations live under `src/PhyloAcc-common/`.

Minimal Test Data
- `tests/data/minimal/` is hand-written synthetic data, not output copied from production analyses.
- `aln.fa` contains three 20 bp sequences arranged as two concatenated 10 bp loci so ST/GT smoke tests run quickly.
- `bed.bed` partitions that alignment into `locus1` and `locus2` using coordinates `0-10` and `10-20`.
- `model.mod` was written manually with the same three-species topology plus a non-singular 4x4 rate matrix so Armadillo/tree parsing works in the tiny test case.
- `tree_coal.tre` was written manually to mirror the same topology in coalescent units for GT tests.
- `tests/data/unit/` stores very small hand-written JSON/Newick fixtures for unit tests that should not embed their data inline in the test code.
- The optional larger integration tier uses files from `../PhyloAcc-test-data/` and is kept separate from the minimal synthetic fixtures.

Latest Suite Status
- Last pytest session: `2026-04-15`
- Command shape used here: `PHYLOACC_RUN_GT=1 /n/home07/gthomas/miniconda3/envs/phyloacc-test/bin/python -m pytest tests/test_st_gt.py tests/test_cpp_unit.py -q`
- Result: `12 passed in 204.00s`
- GT golden status: `tests/golden/minimal/gt_rate_postZ_M0.txt` exists

Status Legend
- `■ PASS`: last observed run passed
- `✖ FAIL`: last observed run failed
- `▲ SKIP`: last observed run skipped by design
- `◆ MIXED`: file had a mixture of passed and skipped tests in the same run
- `△ STALE`: file changed after its last observed run
- `? UNKNOWN`: file has no recorded run yet

Test Inventory

| Test | Type | Purpose | Data | Data origin | Status | Last run | Notes |
|---|---|---|---|---|---|---|---|
| `tests/test_unit_alignment.py` | Python unit | BED parsing, ID filtering, partitioning, alignment stats, low-quality detection, and label validation | `tests/data/minimal/aln.fa` and `tests/data/minimal/bed.bed` plus inline dict fixtures | Minimal files are hand-written; some edge cases are built directly in the test with inline literals. | `■ PASS` | `2026-04-15` | Fast deterministic Python coverage. |
| `tests/test_unit_scf.py` | Python unit | sCF counting, discordant-site accounting, skip handling, and zip vs loop consistency | `tests/data/unit/scf_cases.json` | The quartet alignments are stored as hand-written JSON fixtures so the expected site-pattern counts stay explicit and deterministic. | `■ PASS` | `2026-04-15` | Fast deterministic Python coverage. |
| `tests/test_unit_templates.py` | Python unit | Interface config template generation for ST and GT | Inline string literals only | No external data; template inputs are hand-written format arguments. | `■ PASS` | `2026-04-15` | Locks down emitted config structure. |
| `tests/test_unit_tree_groups.py` | Python unit | Propagation of target, conserved, and outgroup branch categories from tip assignments | `tests/data/unit/tree_groups.nwk` | The test tree is stored as a hand-written Newick fixture file. | `■ PASS` | `2026-04-15` | Covers internal branch categorization logic in tree.py. |
| `tests/test_unit_batch.py` | Python unit | Batch job-file generation for ST and GT configs, including ID files and model-specific options | Temporary files written under `tmp_path` from inline sequences/tree strings | Configs, model file, and coal tree are generated on the fly from hand-written literals. | `■ PASS` | `2026-04-15` | Covers config-writing logic in batch.py without running PhyloAcc. |
| `tests/test_interface.py` | Integration | Summarize-only interface execution on minimal synthetic data | `tests/data/minimal/aln.fa`, `tests/data/minimal/bed.bed`, `tests/data/minimal/model.mod` | These minimal files are hand-written synthetic fixtures in `tests/data/minimal/`. | `■ PASS` | `2026-04-15` | Exercises Python entrypoint without a full workflow run. |
| `tests/test_st_gt.py` | Integration + golden | ST execution on minimal data and GT execution against the ratite subset, both compared to goldens | ST: `tests/data/minimal/*`; GT: ratite subset from `../PhyloAcc-test-data/` plus `tests/golden/minimal/*` | Minimal ST fixtures are hand-written; GT uses a subset selected from the external test-data repo and recorded goldens. | `■ PASS` | `2026-04-15` | Requires PHYLOACC_RUN_GT=1 for GT coverage. |
| `tests/test_optional_testdata.py` | Optional integration + golden | Interface summarize and ST golden comparison using ../PhyloAcc-test-data | `../PhyloAcc-test-data/bioconda-test-data/*` plus `tests/golden/testdata/st_rate_postZ_M0.txt` | This data comes from the external `PhyloAcc-test-data` repo; the test subsets loci via `id-subset.txt`. | `■ PASS` | `2026-04-15` | Requires PHYLOACC_RUN_TESTDATA=1. |
| `tests/test_cpp_unit.py` | Pytest wrapper | Runs lightweight C++ unit binaries and failure-path wrappers | Minimal fixtures in `tests/data/minimal/` plus wrapper-generated temp malformed files | The wrappers use hand-written minimal fixtures and create malformed FASTA/BED inputs on the fly when needed. | `■ PASS` | `2026-04-15` | This is the Python harness for the C++ sanity binaries. |
| `tests/cpp/test_st_main.cpp` | C++ unit binary | ST tree/profile parse sanity on minimal data | `tests/data/minimal/model.mod`, `tests/data/minimal/aln.fa`, `tests/data/minimal/bed.bed` | All of these are hand-written synthetic fixtures under `tests/data/minimal/`. | `■ PASS` | `2026-04-15` | Observed via tests/test_cpp_unit.py. |
| `tests/cpp/test_gt_main.cpp` | C++ unit binary | GT coalescent-tree/profile parse sanity on minimal data | `tests/data/minimal/tree_coal.tre`, `tests/data/minimal/aln.fa`, `tests/data/minimal/bed.bed` | All of these are hand-written synthetic fixtures under `tests/data/minimal/`. | `■ PASS` | `2026-04-15` | Observed via tests/test_cpp_unit.py. |
| `tests/cpp/fail_st_main.cpp` | C++ failure-path binary | Intentional ST load failure to confirm transparent error behavior | No real data file; intentionally missing `.mod` path | The failure is generated by pointing the loader at a nonexistent hand-specified path. | `■ PASS` | `2026-04-15` | Observed via tests/test_cpp_unit.py. |
| `tests/cpp/fail_gt_main.cpp` | C++ failure-path binary | Intentional GT coalescent-tree load failure to confirm transparent error behavior | No real data file; intentionally missing `.tre` path | The failure is generated by pointing the loader at a nonexistent hand-specified path. | `■ PASS` | `2026-04-15` | Observed via tests/test_cpp_unit.py. |
| `tests/cpp/fail_st_profile_main.cpp` | C++ failure-path binary | Intentional ST profile-load failure for malformed FASTA input | Paths supplied by shell wrappers | The wrappers either point to missing files or write malformed FASTA/BED fixtures on the fly. | `■ PASS` | `2026-04-15` | Observed via tests/test_cpp_unit.py. |
| `tests/cpp/fail_gt_profile_main.cpp` | C++ failure-path binary | Intentional GT profile-load failure for malformed FASTA input | Paths supplied by shell wrappers | The wrappers either point to missing files or write malformed FASTA/BED fixtures on the fly. | `■ PASS` | `2026-04-15` | Observed via tests/test_cpp_unit.py. |
| `tests/cpp/wrap_fail_st.sh` | Shell wrapper | Validates the expected ST failure message | No real data file; missing `.mod` path | The wrapper simply invokes the helper against a nonexistent model path. | `■ PASS` | `2026-04-15` | Observed via tests/test_cpp_unit.py. |
| `tests/cpp/wrap_fail_gt.sh` | Shell wrapper | Validates the expected GT failure message | No real data file; missing `.tre` path | The wrapper simply invokes the helper against a nonexistent coalescent-tree path. | `■ PASS` | `2026-04-15` | Observed via tests/test_cpp_unit.py. |

Golden Files
- Record or refresh goldens: `PHYLOACC_RECORD_GOLDEN=1 pytest -q`
- Include optional test-data tier while recording: `PHYLOACC_RUN_TESTDATA=1 PHYLOACC_RECORD_GOLDEN=1 pytest -q`
- Record just the GT golden: `PHYLOACC_RUN_GT=1 PHYLOACC_RECORD_GOLDEN=1 pytest tests/test_st_gt.py::test_gt_minimal_run -q`
- Compare against existing goldens: `pytest -q`
- Tolerances: `atol=1e-6`, `rtol=1e-5`
- Locations:
  - `tests/golden/minimal/`
  - `tests/golden/testdata/`

GT Golden Provenance
- Metadata file: `tests/golden/minimal/gt_rate_postZ_M0.provenance.md`
- `tests/golden/minimal/gt_rate_postZ_M0.txt` should be treated as a conda-forge/bioconda build baseline, not a generic local-compiler baseline.
- The documented matching build path is: `v2.4.5` source plus the conda-forge GCC 14 toolchain.
- Earlier local rebuilds done with the system `g++ 8.5` path did not reproduce this GT golden and should not be used for regression baselines.

Environment Switches
- `PHYLOACC_RUN_GT=1`: enable GT integration coverage
- `PHYLOACC_RUN_TESTDATA=1`: enable optional `../PhyloAcc-test-data` integration tests
- `PHYLOACC_RECORD_GOLDEN=1`: record goldens instead of comparing

C++ Unit Tests
- Build: `make cpp-tests`
- Run binaries directly:
  - `./tests/cpp/phyloacc_cpp_tests_st`
  - `./tests/cpp/phyloacc_cpp_tests_gt`
- Run pytest wrapper:
  - `pytest tests/test_cpp_unit.py -q`
