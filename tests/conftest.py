import os
import sys
import time
from pathlib import Path
import subprocess

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import generate_test_readme as gtr


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "tests" / "data" / "minimal"


@pytest.fixture()
def minimal_data():
    return {
        "aln": DATA_DIR / "aln.fa",
        "bed": DATA_DIR / "bed.bed",
        "mod": DATA_DIR / "model.mod",
        "coal": DATA_DIR / "tree_coal.tre",
    }


def _resolve_bin(env_var, default_name):
    if env_var in os.environ:
        candidate = Path(os.environ[env_var])
    else:
        candidate = ROOT / default_name
    if not candidate.exists():
        pytest.skip(f"{default_name} not found at {candidate}. Build it or set {env_var}.")
    return str(candidate)


def run_cmd(cmd, env=None, cwd=None):
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    proc = subprocess.run(
        cmd,
        env=merged_env,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if proc.returncode != 0:
        raise AssertionError(f"Command failed ({proc.returncode}): {cmd}\n{proc.stdout}")
    return proc.stdout


def _read_tsv_numeric(path):
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip() != ""]
    if not lines:
        raise AssertionError(f"Empty output file: {path}")
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) != len(header):
            raise AssertionError(f"Column mismatch in {path}: {len(parts)} vs {len(header)}")
        try:
            row = [float(p) for p in parts]
        except ValueError as exc:
            raise AssertionError(f"Non-numeric value in {path}: {exc}") from exc
        rows.append(row)
    if not rows:
        raise AssertionError(f"No data rows in {path}")
    return header, rows


def compare_or_record_golden(output_path, golden_path, atol=1e-6, rtol=1e-5):
    output_path = Path(output_path)
    golden_path = Path(golden_path)
    if os.environ.get("PHYLOACC_RECORD_GOLDEN", "0") == "1":
        golden_path.parent.mkdir(parents=True, exist_ok=True)
        golden_path.write_text(output_path.read_text(encoding="utf-8"), encoding="utf-8")
        return
    if not golden_path.exists():
        raise pytest.skip(f"Golden file missing: {golden_path} (set PHYLOACC_RECORD_GOLDEN=1 to create)")

    out_header, out_rows = _read_tsv_numeric(output_path)
    gold_header, gold_rows = _read_tsv_numeric(golden_path)

    if out_header != gold_header:
        raise AssertionError("Header mismatch between output and golden file.")
    if len(out_rows) != len(gold_rows):
        raise AssertionError("Row count mismatch between output and golden file.")

    for i, (out_row, gold_row) in enumerate(zip(out_rows, gold_rows)):
        if len(out_row) != len(gold_row):
            raise AssertionError(f"Column count mismatch in row {i}.")
        for j, (out_val, gold_val) in enumerate(zip(out_row, gold_row)):
            diff = abs(out_val - gold_val)
            tol = atol + rtol * abs(gold_val)
            if diff > tol:
                raise AssertionError(
                    f"Value mismatch at row {i}, col {j}: {out_val} vs {gold_val} "
                    f"(diff {diff} > tol {tol})"
                )


@pytest.fixture()
def phyloacc_st_bin():
    return _resolve_bin("PHYLOACC_ST_BIN", "PhyloAcc-ST")


@pytest.fixture()
def phyloacc_gt_bin():
    if os.environ.get("PHYLOACC_RUN_GT", "0") != "1":
        pytest.skip("GT tests disabled. Set PHYLOACC_RUN_GT=1 to enable.")
    return _resolve_bin("PHYLOACC_GT_BIN", "PhyloAcc-GT")


@pytest.fixture()
def phyloacc_py_env():
    env = {
        "PYTHONPATH": str(ROOT / "src" / "PhyloAcc-interface"),
    }
    return env


@pytest.fixture()
def phyloacc_gt_path():
    return _resolve_bin("PHYLOACC_GT_BIN", "PhyloAcc-GT")


_NODE_OUTCOMES = {}
_SESSION_START = None


def pytest_sessionstart(session):
    global _NODE_OUTCOMES, _SESSION_START
    _NODE_OUTCOMES = {}
    _SESSION_START = time.time()


def pytest_runtest_logreport(report):
    nodeid = report.nodeid
    if report.when == "call":
        _NODE_OUTCOMES[nodeid] = report.outcome
    elif report.when == "setup" and report.outcome in {"skipped", "failed"}:
        _NODE_OUTCOMES[nodeid] = report.outcome
    elif report.when == "teardown" and report.outcome == "failed":
        _NODE_OUTCOMES[nodeid] = "failed"


def _aggregate_file_status(outcomes):
    values = set(outcomes)
    if "failed" in values:
        return "FAIL"
    if values == {"skipped"}:
        return "SKIP"
    if "passed" in values and "skipped" in values:
        return "MIXED"
    if values == {"passed"}:
        return "PASS"
    return "UNKNOWN"


def pytest_sessionfinish(session, exitstatus):
    terminal = session.config.pluginmanager.get_plugin("terminalreporter")
    per_file = {}
    for nodeid, outcome in _NODE_OUTCOMES.items():
        path = nodeid.split("::", 1)[0]
        per_file.setdefault(path, []).append(outcome)

    timestamp = time.strftime("%Y-%m-%dT%H:%M:%S")
    duration = 0.0 if _SESSION_START is None else time.time() - _SESSION_START
    args = list(session.config.invocation_params.args)
    env_prefix = []
    for key in ("PHYLOACC_RUN_GT", "PHYLOACC_RUN_TESTDATA", "PHYLOACC_RECORD_GOLDEN"):
        if os.environ.get(key) == "1":
            env_prefix.append(f"{key}=1")
    command = " ".join(env_prefix + [sys.executable, "-m", "pytest", *args])

    stats = terminal.stats if terminal else {}
    summary_parts = []
    for key in ("passed", "failed", "skipped", "error", "xfailed", "xpassed"):
        count = len(stats.get(key, []))
        if count:
            summary_parts.append(f"{count} {key}")
    summary = ", ".join(summary_parts) if summary_parts else "no tests collected"
    summary = f"{summary} in {duration:.2f}s"

    file_statuses = {
        path: {"status": _aggregate_file_status(outcomes), "last_run": timestamp}
        for path, outcomes in per_file.items()
    }
    gtr.update_status(
        file_statuses,
        {
            "timestamp": timestamp,
            "command": command,
            "summary": summary,
            "exitstatus": exitstatus,
        },
    )
