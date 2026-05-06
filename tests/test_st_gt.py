import math
import os
import platform
from pathlib import Path

import pytest

from conftest import run_cmd, compare_or_record_golden, _read_tsv_numeric


GT_GOLDEN_PACKAGED = Path("tests/golden/minimal/gt_rate_postZ_M0.txt")
GT_GOLDEN_OSX_ARM64_PIXI = Path("tests/golden/minimal/gt_rate_postZ_M0.osx-arm64-clang-pixi.txt")


def _resolve_gt_integration_data(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    testdata_root = repo_root.parent / "PhyloAcc-test-data"
    testdata_dir = testdata_root / "bioconda-test-data"

    aln = testdata_dir / "simu_500_200_diffr_2-1.noanc.fa"
    bed = testdata_dir / "simu_500_200_diffr_2-1.bed"
    mod = testdata_dir / "ratite.mod"
    coal = testdata_root / "ratite.tre"
    ids = testdata_dir / "id-subset.txt"

    if not (aln.exists() and bed.exists() and mod.exists() and coal.exists() and ids.exists()):
        pytest.skip("GT integration test requires ../PhyloAcc-test-data with ratite inputs.")

    id_list = [line.strip() for line in ids.read_text(encoding="utf-8").splitlines() if line.strip()]
    id_set = set(id_list)
    subset_bed = tmp_path / "gt_subset.bed"
    gt_ids = tmp_path / "gt.ids"
    lines = []
    for line in bed.read_text(encoding="utf-8").splitlines():
        parts = line.split("	")
        if len(parts) >= 4 and parts[3] in id_set:
            lines.append(line)
    if not lines:
        pytest.skip("No matching GT subset IDs found in optional test data.")

    subset_bed.write_text("\n".join(lines) + "\n", encoding="utf-8")
    gt_ids.write_text("\n".join(str(int(locus_id) - 1) for locus_id in id_list) + "\n", encoding="utf-8")
    return {"aln": aln, "bed": subset_bed, "mod": mod, "coal": coal, "ids": gt_ids}


def _resolve_gt_golden():
    override = os.environ.get("PHYLOACC_GT_GOLDEN")
    if override:
        return Path(override)

    golden_id = os.environ.get("PHYLOACC_GT_GOLDEN_ID")
    if golden_id:
        if golden_id == "packaged-gcc14":
            return GT_GOLDEN_PACKAGED
        if golden_id == "osx-arm64-clang-pixi":
            return GT_GOLDEN_OSX_ARM64_PIXI
        raise AssertionError(f"Unknown PHYLOACC_GT_GOLDEN_ID: {golden_id}")

    is_osx_arm64_pixi = (
        platform.system() == "Darwin"
        and platform.machine() == "arm64"
        and os.environ.get("PIXI_PROJECT_NAME") == "phyloacc"
        and "clang++" in os.environ.get("CXX", "")
    )
    if is_osx_arm64_pixi:
        return GT_GOLDEN_OSX_ARM64_PIXI

    return GT_GOLDEN_PACKAGED


def _assert_gt_output_invariants(output_path):
    header, rows = _read_tsv_numeric(output_path)
    assert header[:6] == ["No.", "n_rate", "c_rate", "g_rate", "l_rate", "l2_rate"]
    assert len(rows) == 5

    for expected_id, row in enumerate(rows):
        assert row[0] == expected_id
        for value in row:
            assert math.isfinite(value)
        for value in row[1:6]:
            assert value >= 0.0
        for value in row[6:]:
            assert 0.0 <= value <= 1.0


def _assert_status_file(status_path, expected_rows, expected_mode):
    assert status_path.exists(), f"Expected status output not found: {status_path}"
    lines = [line.rstrip("\n") for line in status_path.read_text(encoding="utf-8").splitlines()]
    assert lines[0].split("\t") == [
        "chain",
        "No.",
        "element_name",
        "mode",
        "status",
        "completed_models",
        "message",
    ]
    rows = [line.split("\t") for line in lines[1:]]
    assert len(rows) == expected_rows
    assert {row[3] for row in rows} == {expected_mode}
    assert all(row[4] != "error" for row in rows)


def _write_st_cfg(cfg_path, minimal_data, out_dir):
    cfg_path.write_text(
        "\n".join(
            [
                f"PHYTREE_FILE {minimal_data['mod']}",
                f"ALIGN_FILE {minimal_data['aln']}",
                f"SEG_FILE {minimal_data['bed']}",
                f"RESULT_FOLDER {out_dir}",
                "PREFIX test",
                "BURNIN 5",
                "MCMC 10",
                "ADAPT_FREQ 1",
                "CHAIN 1",
                "TARGETSPECIES sp1",
                "OUTGROUP sp3",
                "CONSERVE sp2",
                "NUM_THREAD 1",
                "MIN_LEN 1",
            ]
        )
        + "\n"
    )


def _write_gt_cfg(cfg_path, minimal_data, out_dir, id_file):
    cfg_path.write_text(
        "\n".join(
            [
                f"PHYTREE_FILE {minimal_data['mod']}",
                f"TREE_IN_COALESCENT_UNIT {minimal_data['coal']}",
                f"ALIGN_FILE {minimal_data['aln']}",
                f"SEG_FILE {minimal_data['bed']}",
                f"ID_FILE {id_file}",
                f"RESULT_FOLDER {out_dir}",
                "PREFIX test",
                "BURNIN 5",
                "MCMC 10",
                "THIN 1",
                "CHAIN 1",
                "TARGETSPECIES sp1",
                "OUTGROUP sp3",
                "CONSERVE sp2",
                "NUM_THREAD 1",
                "MIN_LEN 1",
            ]
        )
        + "\n"
    )


def test_st_minimal_run(minimal_data, tmp_path, phyloacc_st_bin):
    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = tmp_path / "st.cfg"
    _write_st_cfg(cfg_path, minimal_data, out_dir)

    run_cmd([phyloacc_st_bin, str(cfg_path)])

    out_file = out_dir / "test_rate_postZ_M0.txt"
    assert out_file.exists(), f"Expected output not found: {out_file}"
    _assert_status_file(out_dir / "test_elem_status.txt", expected_rows=2, expected_mode="ST")
    compare_or_record_golden(
        out_file,
        Path("tests/golden/minimal/st_rate_postZ_M0.txt"),
        atol=1e-6,
        rtol=1e-5,
    )


def test_gt_minimal_run(minimal_data, tmp_path, phyloacc_gt_bin):
    gt_data = _resolve_gt_integration_data(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = tmp_path / "gt.cfg"
    _write_gt_cfg(cfg_path, gt_data, out_dir, gt_data["ids"])

    run_cmd([phyloacc_gt_bin, str(cfg_path)])

    out_file = out_dir / "test_rate_postZ_M0.txt"
    assert out_file.exists(), f"Expected output not found: {out_file}"
    _assert_gt_output_invariants(out_file)
    expected_status_rows = len([line for line in gt_data["ids"].read_text(encoding="utf-8").splitlines() if line.strip()])
    _assert_status_file(out_dir / "test_elem_status.txt", expected_status_rows, expected_mode="GT")
    compare_or_record_golden(
        out_file,
        _resolve_gt_golden(),
        atol=1e-6,
        rtol=1e-5,
    )
