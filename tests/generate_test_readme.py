import json
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TESTS_DIR = ROOT / "tests"
README_PATH = TESTS_DIR / "README.md"
STATUS_PATH = TESTS_DIR / ".last_test_status.json"
GT_GOLDEN = ROOT / "tests/golden/minimal/gt_rate_postZ_M0.txt"
GT_GOLDEN_PROVENANCE = ROOT / "tests/golden/minimal/gt_rate_postZ_M0.provenance.md"

MANIFEST = [
    {
        "path": "tests/test_unit_alignment.py",
        "type": "Python unit",
        "purpose": "BED parsing, ID filtering, partitioning, alignment stats, low-quality detection, and label validation",
        "notes": "Fast deterministic Python coverage.",
        "data": "`tests/data/minimal/aln.fa` and `tests/data/minimal/bed.bed` plus inline dict fixtures",
        "data_gen": "Minimal files are hand-written; some edge cases are built directly in the test with inline literals.",
    },
    {
        "path": "tests/test_unit_scf.py",
        "type": "Python unit",
        "purpose": "sCF counting, discordant-site accounting, skip handling, and zip vs loop consistency",
        "notes": "Fast deterministic Python coverage.",
        "data": "`tests/data/unit/scf_cases.json`",
        "data_gen": "The quartet alignments are stored as hand-written JSON fixtures so the expected site-pattern counts stay explicit and deterministic.",
    },
    {
        "path": "tests/test_unit_templates.py",
        "type": "Python unit",
        "purpose": "Interface config template generation for ST and GT",
        "notes": "Locks down emitted config structure.",
        "data": "Inline string literals only",
        "data_gen": "No external data; template inputs are hand-written format arguments.",
    },
    {
        "path": "tests/test_unit_tree_groups.py",
        "type": "Python unit",
        "purpose": "Propagation of target, conserved, and outgroup branch categories from tip assignments",
        "notes": "Covers internal branch categorization logic in tree.py.",
        "data": "`tests/data/unit/tree_groups.nwk`",
        "data_gen": "The test tree is stored as a hand-written Newick fixture file.",
    },
    {
        "path": "tests/test_unit_batch.py",
        "type": "Python unit",
        "purpose": "Batch job-file generation for ST and GT configs, including ID files and model-specific options",
        "notes": "Covers config-writing logic in batch.py without running PhyloAcc.",
        "data": "Temporary files written under `tmp_path` from inline sequences/tree strings",
        "data_gen": "Configs, model file, and coal tree are generated on the fly from hand-written literals.",
    },
    {
        "path": "tests/test_interface.py",
        "type": "Integration",
        "purpose": "Summarize-only interface execution on minimal synthetic data",
        "notes": "Exercises Python entrypoint without a full workflow run.",
        "data": "`tests/data/minimal/aln.fa`, `tests/data/minimal/bed.bed`, `tests/data/minimal/model.mod`",
        "data_gen": "These minimal files are hand-written synthetic fixtures in `tests/data/minimal/`.",
    },
    {
        "path": "tests/test_st_gt.py",
        "type": "Integration + golden",
        "purpose": "ST execution on minimal data and GT execution against the ratite subset, both compared to goldens",
        "notes": "Requires PHYLOACC_RUN_GT=1 for GT coverage.",
        "data": "ST: `tests/data/minimal/*`; GT: ratite subset from `../PhyloAcc-test-data/` plus `tests/golden/minimal/*`",
        "data_gen": "Minimal ST fixtures are hand-written; GT uses a subset selected from the external test-data repo and recorded goldens.",
    },
    {
        "path": "tests/test_optional_testdata.py",
        "type": "Optional integration + golden",
        "purpose": "Interface summarize and ST golden comparison using ../PhyloAcc-test-data",
        "notes": "Requires PHYLOACC_RUN_TESTDATA=1.",
        "data": "`../PhyloAcc-test-data/bioconda-test-data/*` plus `tests/golden/testdata/st_rate_postZ_M0.txt`",
        "data_gen": "This data comes from the external `PhyloAcc-test-data` repo; the test subsets loci via `id-subset.txt`.",
    },
    {
        "path": "tests/test_cpp_unit.py",
        "type": "Pytest wrapper",
        "purpose": "Runs lightweight C++ unit binaries and failure-path wrappers",
        "notes": "This is the Python harness for the C++ sanity binaries.",
        "data": "Minimal fixtures in `tests/data/minimal/` plus wrapper-generated temp malformed files",
        "data_gen": "The wrappers use hand-written minimal fixtures and create malformed FASTA/BED inputs on the fly when needed.",
    },
    {
        "path": "tests/cpp/test_st_main.cpp",
        "type": "C++ unit binary",
        "purpose": "ST tree/profile parse sanity on minimal data",
        "notes": "Observed via tests/test_cpp_unit.py.",
        "data": "`tests/data/minimal/model.mod`, `tests/data/minimal/aln.fa`, `tests/data/minimal/bed.bed`",
        "data_gen": "All of these are hand-written synthetic fixtures under `tests/data/minimal/`.",
        "status_key": "tests/test_cpp_unit.py",
        "watch_files": ["tests/cpp/test_st_main.cpp", "tests/test_cpp_unit.py"],
    },
    {
        "path": "tests/cpp/test_gt_main.cpp",
        "type": "C++ unit binary",
        "purpose": "GT coalescent-tree/profile parse sanity on minimal data",
        "notes": "Observed via tests/test_cpp_unit.py.",
        "data": "`tests/data/minimal/tree_coal.tre`, `tests/data/minimal/aln.fa`, `tests/data/minimal/bed.bed`",
        "data_gen": "All of these are hand-written synthetic fixtures under `tests/data/minimal/`.",
        "status_key": "tests/test_cpp_unit.py",
        "watch_files": ["tests/cpp/test_gt_main.cpp", "tests/test_cpp_unit.py"],
    },
    {
        "path": "tests/cpp/fail_st_main.cpp",
        "type": "C++ failure-path binary",
        "purpose": "Intentional ST load failure to confirm transparent error behavior",
        "notes": "Observed via tests/test_cpp_unit.py.",
        "data": "No real data file; intentionally missing `.mod` path",
        "data_gen": "The failure is generated by pointing the loader at a nonexistent hand-specified path.",
        "status_key": "tests/test_cpp_unit.py",
        "watch_files": ["tests/cpp/fail_st_main.cpp", "tests/test_cpp_unit.py", "tests/cpp/wrap_fail_st.sh"],
    },
    {
        "path": "tests/cpp/fail_gt_main.cpp",
        "type": "C++ failure-path binary",
        "purpose": "Intentional GT coalescent-tree load failure to confirm transparent error behavior",
        "notes": "Observed via tests/test_cpp_unit.py.",
        "data": "No real data file; intentionally missing `.tre` path",
        "data_gen": "The failure is generated by pointing the loader at a nonexistent hand-specified path.",
        "status_key": "tests/test_cpp_unit.py",
        "watch_files": ["tests/cpp/fail_gt_main.cpp", "tests/test_cpp_unit.py", "tests/cpp/wrap_fail_gt.sh"],
    },
    {
        "path": "tests/cpp/fail_st_profile_main.cpp",
        "type": "C++ failure-path binary",
        "purpose": "Intentional ST profile-load failure for malformed FASTA input",
        "notes": "Observed via tests/test_cpp_unit.py.",
        "data": "Paths supplied by shell wrappers",
        "data_gen": "The wrappers either point to missing files or write malformed FASTA/BED fixtures on the fly.",
        "status_key": "tests/test_cpp_unit.py",
        "watch_files": ["tests/cpp/fail_st_profile_main.cpp", "tests/test_cpp_unit.py", "tests/cpp/wrap_fail_st_profile.sh"],
    },
    {
        "path": "tests/cpp/fail_gt_profile_main.cpp",
        "type": "C++ failure-path binary",
        "purpose": "Intentional GT profile-load failure for malformed FASTA input",
        "notes": "Observed via tests/test_cpp_unit.py.",
        "data": "Paths supplied by shell wrappers",
        "data_gen": "The wrappers either point to missing files or write malformed FASTA/BED fixtures on the fly.",
        "status_key": "tests/test_cpp_unit.py",
        "watch_files": ["tests/cpp/fail_gt_profile_main.cpp", "tests/test_cpp_unit.py", "tests/cpp/wrap_fail_gt_profile.sh"],
    },
    {
        "path": "tests/cpp/wrap_fail_st.sh",
        "type": "Shell wrapper",
        "purpose": "Validates the expected ST failure message",
        "notes": "Observed via tests/test_cpp_unit.py.",
        "data": "No real data file; missing `.mod` path",
        "data_gen": "The wrapper simply invokes the helper against a nonexistent model path.",
        "status_key": "tests/test_cpp_unit.py",
        "watch_files": ["tests/cpp/wrap_fail_st.sh", "tests/test_cpp_unit.py"],
    },
    {
        "path": "tests/cpp/wrap_fail_gt.sh",
        "type": "Shell wrapper",
        "purpose": "Validates the expected GT failure message",
        "notes": "Observed via tests/test_cpp_unit.py.",
        "data": "No real data file; missing `.tre` path",
        "data_gen": "The wrapper simply invokes the helper against a nonexistent coalescent-tree path.",
        "status_key": "tests/test_cpp_unit.py",
        "watch_files": ["tests/cpp/wrap_fail_gt.sh", "tests/test_cpp_unit.py"],
    },
]

STATUS_DISPLAY = {
    "PASS": "■ PASS",
    "FAIL": "✖ FAIL",
    "SKIP": "▲ SKIP",
    "MIXED": "◆ MIXED",
    "UNKNOWN": "? UNKNOWN",
    "STALE": "△ STALE",
}


def load_status():
    if not STATUS_PATH.exists():
        return {"files": {}, "session": {}}
    return json.loads(STATUS_PATH.read_text(encoding="utf-8"))


def write_status(data):
    STATUS_PATH.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def format_last_run(last_run):
    if not last_run:
        return "never"
    return last_run.split("T", 1)[0]


def entry_status(entry, data):
    key = entry.get("status_key", entry["path"])
    record = data.get("files", {}).get(key)
    if not record:
        return STATUS_DISPLAY["UNKNOWN"], "never", entry["notes"]

    status = record.get("status", "UNKNOWN")
    last_run = record.get("last_run", "")
    note = entry["notes"]

    if last_run:
        last_dt = datetime.fromisoformat(last_run)
        for relpath in entry.get("watch_files", [entry["path"]]):
            abspath = ROOT / relpath
            if abspath.exists() and datetime.fromtimestamp(abspath.stat().st_mtime) > last_dt:
                status = "STALE"
                note = note + " Changed since the last observed run."
                break

    return STATUS_DISPLAY.get(status, STATUS_DISPLAY["UNKNOWN"]), format_last_run(last_run), note


def summary_lines(data):
    session = data.get("session", {})
    golden_state = "exists" if GT_GOLDEN.exists() else "is missing"
    if not session:
        return [
            "- Last pytest session: `never`",
            "- Result: `no recorded pytest session yet`",
            f"- GT golden status: `{GT_GOLDEN.relative_to(ROOT)}` {golden_state}",
        ]

    return [
        f"- Last pytest session: `{session.get('timestamp', 'unknown')[:10]}`",
        f"- Command shape used here: `{session.get('command', 'unknown')}`",
        f"- Result: `{session.get('summary', 'unknown')}`",
        f"- GT golden status: `{GT_GOLDEN.relative_to(ROOT)}` {golden_state}",
    ]


def render_readme(data=None):
    data = load_status() if data is None else data
    lines = [
        "Test Suite (Pytest)",
        "",
        "Purpose",
        "- Protect refactors by locking down Python logic, C++ parsing/invariants, and end-to-end ST/GT behavior.",
        "- Keep fast deterministic checks separate from heavier integration runs.",
        "- Use golden files to catch numeric drift.",
        "",
        "Quick Start",
        "- Build binaries: `make`",
        "- Build C++ unit binaries: `make cpp-tests`",
        "- Run default suite: `pytest -q`",
        "- Run GT and optional external-data coverage: `PHYLOACC_RUN_GT=1 PHYLOACC_RUN_TESTDATA=1 pytest -q`",
        "",
        "Environment",
        "- This repo expects the conda environment `phyloacc-test`.",
        "- The C++ unit binaries in `tests/cpp/` link against libraries from `CONDA_PREFIX`.",
        "- Regression validation for ST/GT binaries assumes the supported GCC14 wrapper/toolchain path, not an arbitrary local `g++`.",
        "- Example supported rebuild:",
        "  - `conda run -n phyloacc-test make -B CXX=/n/holylfs05/LABS/informatics/Lab/projects/gwct/phyloacc/PhyloAcc/dev/cc14-wrap.sh PREFIX=/n/home07/gthomas/miniconda3/envs/phyloacc-test PhyloAcc-ST PhyloAcc-GT`",
        "- The duplicate `profile.*`, `newick.*`, and `utils.*` files still present under `src/PhyloAcc-ST/` and `src/PhyloAcc-GT/` are legacy copies excluded by `Makefile`; the active shared implementations live under `src/PhyloAcc-common/`.",
        "",
        "Minimal Test Data",
        "- `tests/data/minimal/` is hand-written synthetic data, not output copied from production analyses.",
        "- `aln.fa` contains three 20 bp sequences arranged as two concatenated 10 bp loci so ST/GT smoke tests run quickly.",
        "- `bed.bed` partitions that alignment into `locus1` and `locus2` using coordinates `0-10` and `10-20`.",
        "- `model.mod` was written manually with the same three-species topology plus a non-singular 4x4 rate matrix so Armadillo/tree parsing works in the tiny test case.",
        "- `tree_coal.tre` was written manually to mirror the same topology in coalescent units for GT tests.",
        "- `tests/data/unit/` stores very small hand-written JSON/Newick fixtures for unit tests that should not embed their data inline in the test code.",
        "- The optional larger integration tier uses files from `../PhyloAcc-test-data/` and is kept separate from the minimal synthetic fixtures.",
        "",
        "Latest Suite Status",
    ]
    lines.extend(summary_lines(data))
    lines.extend([
        "",
        "Status Legend",
        "- `■ PASS`: last observed run passed",
        "- `✖ FAIL`: last observed run failed",
        "- `▲ SKIP`: last observed run skipped by design",
        "- `◆ MIXED`: file had a mixture of passed and skipped tests in the same run",
        "- `△ STALE`: file changed after its last observed run",
        "- `? UNKNOWN`: file has no recorded run yet",
        "",
        "Test Inventory",
        "",
        "| Test | Type | Purpose | Data | Data origin | Status | Last run | Notes |",
        "|---|---|---|---|---|---|---|---|",
    ])
    for entry in MANIFEST:
        status, last_run, note = entry_status(entry, data)
        lines.append(f"| `{entry['path']}` | {entry['type']} | {entry['purpose']} | {entry.get('data', 'NA')} | {entry.get('data_gen', 'NA')} | `{status}` | `{last_run}` | {note} |")
    lines.extend([
        "",
        "Golden Files",
        "- Record or refresh goldens: `PHYLOACC_RECORD_GOLDEN=1 pytest -q`",
        "- Include optional test-data tier while recording: `PHYLOACC_RUN_TESTDATA=1 PHYLOACC_RECORD_GOLDEN=1 pytest -q`",
        "- Record just the GT golden: `PHYLOACC_RUN_GT=1 PHYLOACC_RECORD_GOLDEN=1 pytest tests/test_st_gt.py::test_gt_minimal_run -q`",
        "- Compare against existing goldens: `pytest -q`",
        "- Tolerances: `atol=1e-6`, `rtol=1e-5`",
        "- Locations:",
        "  - `tests/golden/minimal/`",
        "  - `tests/golden/testdata/`",
        "",
        "GT Golden Provenance",
        f"- Metadata file: `{GT_GOLDEN_PROVENANCE.relative_to(ROOT)}`",
        "- `tests/golden/minimal/gt_rate_postZ_M0.txt` should be treated as a conda-forge/bioconda build baseline, not a generic local-compiler baseline.",
        "- The documented matching build path is: `v2.4.5` source plus the conda-forge GCC 14 toolchain.",
        "- Earlier local rebuilds done with the system `g++ 8.5` path did not reproduce this GT golden and should not be used for regression baselines.",
        "",
        "Environment Switches",
        "- `PHYLOACC_RUN_GT=1`: enable GT integration coverage",
        "- `PHYLOACC_RUN_TESTDATA=1`: enable optional `../PhyloAcc-test-data` integration tests",
        "- `PHYLOACC_RECORD_GOLDEN=1`: record goldens instead of comparing",
        "",
        "C++ Unit Tests",
        "- Build: `make cpp-tests`",
        "- Run binaries directly:",
        "  - `./tests/cpp/phyloacc_cpp_tests_st`",
        "  - `./tests/cpp/phyloacc_cpp_tests_gt`",
        "- Run pytest wrapper:",
        "  - `pytest tests/test_cpp_unit.py -q`",
    ])
    README_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_status(file_statuses, session_summary):
    data = load_status()
    data.setdefault("files", {})
    data["session"] = session_summary
    for path, record in file_statuses.items():
        data["files"][path] = record
    write_status(data)
    render_readme(data)


if __name__ == "__main__":
    render_readme()
