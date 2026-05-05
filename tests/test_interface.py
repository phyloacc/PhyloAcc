import sys

from conftest import run_cmd


def test_interface_summarize(minimal_data, tmp_path, phyloacc_py_env, phyloacc_st_bin, phyloacc_gt_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        "\n".join(
            [
                f"aln_file: {minimal_data['aln']}",
                f"bed_file: {minimal_data['bed']}",
                f"mod_file: {minimal_data['mod']}",
                "targets: sp1",
                "conserved: sp2",
                "outgroup: sp3",
                "run_mode: st",
                "num_procs: 1",
                f"out_dir: {out_dir}",
                f"phyloacc_st_path: {phyloacc_st_bin}",
                f"phyloacc_gt_path: {phyloacc_gt_path}",
                "overwrite_flag: True",
                "summarize_flag: True",
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
    run_cmd(cmd, env=phyloacc_py_env)

    summary = out_dir / "phyloacc-pre-run-summary.html"
    assert summary.exists(), f"Expected summary file not found: {summary}"
