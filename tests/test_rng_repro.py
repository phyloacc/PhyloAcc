import os
from pathlib import Path

import pytest

from conftest import run_cmd
from test_st_gt import _resolve_gt_integration_data


pytestmark = pytest.mark.skipif(
    os.environ.get("PHYLOACC_RUN_RNG", "0") != "1",
    reason="RNG reproducibility tests disabled. Set PHYLOACC_RUN_RNG=1 to enable.",
)


def _write_st_cfg(cfg_path, minimal_data, out_dir, num_thread):
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
                f"NUM_THREAD {num_thread}",
                "MIN_LEN 1",
                "SEED 1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_gt_cfg(cfg_path, gt_data, out_dir, num_thread):
    cfg_path.write_text(
        "\n".join(
            [
                f"PHYTREE_FILE {gt_data['mod']}",
                f"TREE_IN_COALESCENT_UNIT {gt_data['coal']}",
                f"ALIGN_FILE {gt_data['aln']}",
                f"SEG_FILE {gt_data['bed']}",
                f"ID_FILE {gt_data['ids']}",
                f"RESULT_FOLDER {out_dir}",
                "PREFIX test",
                "BURNIN 5",
                "MCMC 10",
                "THIN 1",
                "CHAIN 1",
                "TARGETSPECIES sp1",
                "OUTGROUP sp3",
                "CONSERVE sp2",
                f"NUM_THREAD {num_thread}",
                "MIN_LEN 1",
                "SEED 1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _output_files(out_dir):
    files = {path.name: path for path in Path(out_dir).glob("*.txt")}
    if not files:
        raise AssertionError(f"No output text files found in {out_dir}")
    return files


def _assert_outputs_repeatable(first_dir, second_dir):
    first_files = _output_files(first_dir)
    second_files = _output_files(second_dir)
    assert first_files.keys() == second_files.keys()

    for name in sorted(first_files):
        assert first_files[name].read_bytes() == second_files[name].read_bytes(), name


@pytest.mark.parametrize("num_thread", [1, 2])
def test_st_rng_reproducible_same_thread_count(minimal_data, tmp_path, phyloacc_st_bin, num_thread):
    out1 = tmp_path / f"st_t{num_thread}_run1"
    out2 = tmp_path / f"st_t{num_thread}_run2"
    out1.mkdir()
    out2.mkdir()

    cfg1 = tmp_path / f"st_t{num_thread}_run1.cfg"
    cfg2 = tmp_path / f"st_t{num_thread}_run2.cfg"
    _write_st_cfg(cfg1, minimal_data, out1, num_thread)
    _write_st_cfg(cfg2, minimal_data, out2, num_thread)

    run_cmd([phyloacc_st_bin, str(cfg1)])
    run_cmd([phyloacc_st_bin, str(cfg2)])

    _assert_outputs_repeatable(out1, out2)


@pytest.mark.skipif(
    os.environ.get("PHYLOACC_RUN_GT", "0") != "1",
    reason="GT RNG reproducibility tests disabled. Set PHYLOACC_RUN_GT=1 to enable.",
)
@pytest.mark.parametrize("num_thread", [1, 2])
def test_gt_rng_reproducible_same_thread_count(tmp_path, phyloacc_gt_bin, num_thread):
    gt_data = _resolve_gt_integration_data(tmp_path)
    out1 = tmp_path / f"gt_t{num_thread}_run1"
    out2 = tmp_path / f"gt_t{num_thread}_run2"
    out1.mkdir()
    out2.mkdir()

    cfg1 = tmp_path / f"gt_t{num_thread}_run1.cfg"
    cfg2 = tmp_path / f"gt_t{num_thread}_run2.cfg"
    _write_gt_cfg(cfg1, gt_data, out1, num_thread)
    _write_gt_cfg(cfg2, gt_data, out2, num_thread)

    run_cmd([phyloacc_gt_bin, str(cfg1)])
    run_cmd([phyloacc_gt_bin, str(cfg2)])

    _assert_outputs_repeatable(out1, out2)


def test_st_rng_reproducible_across_thread_counts(minimal_data, tmp_path, phyloacc_st_bin):
    out1 = tmp_path / "st_t1"
    out2 = tmp_path / "st_t2"
    out1.mkdir()
    out2.mkdir()

    cfg1 = tmp_path / "st_t1.cfg"
    cfg2 = tmp_path / "st_t2.cfg"
    _write_st_cfg(cfg1, minimal_data, out1, 1)
    _write_st_cfg(cfg2, minimal_data, out2, 2)

    run_cmd([phyloacc_st_bin, str(cfg1)])
    run_cmd([phyloacc_st_bin, str(cfg2)])

    _assert_outputs_repeatable(out1, out2)


@pytest.mark.skipif(
    os.environ.get("PHYLOACC_RUN_GT", "0") != "1",
    reason="GT RNG reproducibility tests disabled. Set PHYLOACC_RUN_GT=1 to enable.",
)
def test_gt_rng_reproducible_across_thread_counts(tmp_path, phyloacc_gt_bin):
    gt_data = _resolve_gt_integration_data(tmp_path)
    out1 = tmp_path / "gt_t1"
    out2 = tmp_path / "gt_t2"
    out1.mkdir()
    out2.mkdir()

    cfg1 = tmp_path / "gt_t1.cfg"
    cfg2 = tmp_path / "gt_t2.cfg"
    _write_gt_cfg(cfg1, gt_data, out1, 1)
    _write_gt_cfg(cfg2, gt_data, out2, 2)

    run_cmd([phyloacc_gt_bin, str(cfg1)])
    run_cmd([phyloacc_gt_bin, str(cfg2)])

    _assert_outputs_repeatable(out1, out2)
