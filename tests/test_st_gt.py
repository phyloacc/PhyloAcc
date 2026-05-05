from pathlib import Path

import pytest

from conftest import run_cmd, compare_or_record_golden


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
    compare_or_record_golden(
        out_file,
        Path("tests/golden/minimal/gt_rate_postZ_M0.txt"),
        atol=1e-6,
        rtol=1e-5,
    )
