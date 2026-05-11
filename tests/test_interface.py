import sys
from pathlib import Path

import pytest
from conftest import run_cmd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "PhyloAcc-interface"))

import phyloacc_lib.opt_parse as OP
import phyloacc_lib.params as PARAMS


def _test_globs():
    globs = PARAMS.init()
    globs["log-v"] = -1
    return globs


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


def test_getopt_missing_file_and_dir_fail_cleanly(tmp_path, capsys):
    globs = _test_globs()
    flags = {"mod_file": "-m", "aln_dir": "-d"}

    with pytest.raises(SystemExit):
        OP.getOpt(str(tmp_path / "missing.mod"), "mod_file", "FILE", False, {}, flags, globs)
    assert "does not exist or is not a file" in capsys.readouterr().out

    with pytest.raises(SystemExit):
        OP.getOpt(str(tmp_path / "missing-dir"), "aln_dir", "DIR", False, {}, flags, globs)
    assert "does not exist or is not a directory" in capsys.readouterr().out


def test_labeltree_and_labelmod_parse_independently(monkeypatch):
    globs = _test_globs()
    captured = {}

    def stop_after_label_options(parsed_globs, dep_check, dev_opt):
        captured["label-tree"] = parsed_globs["label-tree"]
        captured["label-mod"] = parsed_globs["label-mod"]
        raise RuntimeError("stop after label option parsing")

    monkeypatch.setattr(OP.PC, "execCheck", stop_after_label_options)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "phyloacc.py",
            "--quiet",
            "--depcheck",
            "--labelmod",
            "-st-path",
            "true",
            "-gt-path",
            "true",
        ],
    )

    with pytest.raises(RuntimeError, match="stop after label option parsing"):
        OP.optParse(globs)

    assert captured == {"label-tree": False, "label-mod": True}
