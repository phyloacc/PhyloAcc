import os
import sys
from pathlib import Path

import pytest

from conftest import run_cmd, compare_or_record_golden


def test_optional_testdata_interface(tmp_path, phyloacc_py_env, phyloacc_st_bin, phyloacc_gt_path):
    if os.environ.get("PHYLOACC_RUN_TESTDATA", "0") != "1":
        pytest.skip("Optional test-data run disabled. Set PHYLOACC_RUN_TESTDATA=1 to enable.")

    repo_root = Path(__file__).resolve().parents[1]
    testdata_dir = repo_root.parent / "PhyloAcc-test-data" / "bioconda-test-data"

    aln = testdata_dir / "simu_500_200_diffr_2-1.noanc.fa"
    bed = testdata_dir / "simu_500_200_diffr_2-1.bed"
    mod = testdata_dir / "ratite.mod"
    ids = testdata_dir / "id-subset.txt"

    if not (aln.exists() and bed.exists() and mod.exists() and ids.exists()):
        pytest.skip("PhyloAcc-test-data not found. Skipping optional test-data run.")

    out_dir = tmp_path / "phyloacc-bioconda-test"
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg_path = tmp_path / "testdata-config.yaml"
    cfg_path.write_text(
        "\n".join(
            [
                f"aln_file: {aln}",
                f"bed_file: {bed}",
                f"id_file: {ids}",
                f"mod_file: {mod}",
                "targets: strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid",
                "outgroup: allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar",
                "run_mode: st",
                "num_procs: 1",
                f"out_dir: {out_dir}",
                "overwrite_flag: True",
                "summarize_flag: True",
                f"phyloacc_st_path: {phyloacc_st_bin}",
                f"phyloacc_gt_path: {phyloacc_gt_path}",
            ]
        )
        + "\n"
    )

    cmd = [
        sys.executable,
        "src/PhyloAcc-interface/phyloacc.py",
        "--config",
        str(cfg_path),
        "--local",
    ]
    run_cmd(cmd, env=phyloacc_py_env, cwd=str(repo_root))

    summary = out_dir / "phyloacc-pre-run-summary.html"
    assert summary.exists(), f"Expected summary file not found: {summary}"


def test_optional_testdata_st_golden(tmp_path, phyloacc_st_bin):
    if os.environ.get("PHYLOACC_RUN_TESTDATA", "0") != "1":
        pytest.skip("Optional test-data run disabled. Set PHYLOACC_RUN_TESTDATA=1 to enable.")

    repo_root = Path(__file__).resolve().parents[1]
    testdata_dir = repo_root.parent / "PhyloAcc-test-data" / "bioconda-test-data"

    aln = testdata_dir / "simu_500_200_diffr_2-1.noanc.fa"
    bed = testdata_dir / "simu_500_200_diffr_2-1.bed"
    mod = testdata_dir / "ratite.mod"
    ids = testdata_dir / "id-subset.txt"

    if not (aln.exists() and bed.exists() and mod.exists() and ids.exists()):
        pytest.skip("PhyloAcc-test-data not found. Skipping optional ST golden run.")

    subset_bed = tmp_path / "subset.bed"
    id_set = {line.strip() for line in ids.read_text(encoding="utf-8").splitlines() if line.strip()}
    lines = []
    for line in bed.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 4 and parts[3] in id_set:
            lines.append(line)
    if not lines:
        pytest.skip("No matching IDs found in bed file for subset.")
    subset_bed.write_text("\n".join(lines) + "\n", encoding="utf-8")

    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = tmp_path / "st.cfg"
    cfg_path.write_text(
        "\n".join(
            [
                f"PHYTREE_FILE {mod}",
                f"ALIGN_FILE {aln}",
                f"SEG_FILE {subset_bed}",
                f"RESULT_FOLDER {out_dir}",
                "PREFIX testdata",
                "BURNIN 5",
                "MCMC 10",
                "ADAPT_FREQ 1",
                "CHAIN 1",
                "TARGETSPECIES strCam;rhePen;rheAme;casCas;droNov;aptRow;aptHaa;aptOwe;anoDid",
                "OUTGROUP allMis;allSin;croPor;gavGan;chrPic;cheMyd;anoCar",
                "NUM_THREAD 1",
                "MIN_LEN 1",
            ]
        )
        + "\n"
    )

    run_cmd([phyloacc_st_bin, str(cfg_path)])

    out_file = out_dir / "testdata_rate_postZ_M0.txt"
    assert out_file.exists(), f"Expected output not found: {out_file}"
    compare_or_record_golden(
        out_file,
        Path("tests/golden/testdata/st_rate_postZ_M0.txt"),
        atol=1e-6,
        rtol=1e-5,
    )
