import sys
from pathlib import Path

from conftest import run_cmd


CPP_DIR = Path("tests/cpp")
ENV_PREFIX = Path(sys.executable).resolve().parents[1]
CPP_ENV = {"LD_LIBRARY_PATH": str(ENV_PREFIX / "lib")}


def _cpp_bin(name):
    path = CPP_DIR / name
    if not path.exists():
        raise AssertionError(f"Missing C++ test binary: {path}. Build with `make cpp-tests`.")
    return str(path)


def test_cpp_st_sanity_binary_runs():
    output = run_cmd([_cpp_bin("phyloacc_cpp_tests_st")], env=CPP_ENV)
    assert "ST C++ unit tests passed." in output


def test_cpp_gt_sanity_binary_runs():
    output = run_cmd([_cpp_bin("phyloacc_cpp_tests_gt")], env=CPP_ENV)
    assert "GT C++ unit tests passed." in output


def test_cpp_st_missing_file_fails():
    output = run_cmd(["bash", str(CPP_DIR / "wrap_fail_st.sh")], env=CPP_ENV)
    assert "Observed expected ST failure." in output


def test_cpp_gt_missing_file_fails():
    output = run_cmd(["bash", str(CPP_DIR / "wrap_fail_gt.sh")], env=CPP_ENV)
    assert "Observed expected GT failure." in output


def test_cpp_st_malformed_profile_fails():
    output = run_cmd(["bash", str(CPP_DIR / "wrap_fail_st_profile.sh")], env=CPP_ENV)
    assert "Observed expected ST profile failure." in output


def test_cpp_gt_malformed_profile_fails():
    output = run_cmd(["bash", str(CPP_DIR / "wrap_fail_gt_profile.sh")], env=CPP_ENV)
    assert "Observed expected GT profile failure." in output


def test_cpp_st_missing_profile_fails():
    output = run_cmd(["bash", str(CPP_DIR / "wrap_fail_st_missing_profile.sh")], env=CPP_ENV)
    assert "Observed expected ST missing-profile failure." in output


def test_cpp_gt_missing_profile_fails():
    output = run_cmd(["bash", str(CPP_DIR / "wrap_fail_gt_missing_profile.sh")], env=CPP_ENV)
    assert "Observed expected GT missing-profile failure." in output


def test_cpp_st_missing_bed_fails():
    output = run_cmd(["bash", str(CPP_DIR / "wrap_fail_st_missing_bed.sh")], env=CPP_ENV)
    assert "Observed expected ST missing-bed failure." in output


def test_cpp_gt_missing_bed_fails():
    output = run_cmd(["bash", str(CPP_DIR / "wrap_fail_gt_missing_bed.sh")], env=CPP_ENV)
    assert "Observed expected GT missing-bed failure." in output
