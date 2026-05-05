import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "PhyloAcc-interface"))

from phyloacc_lib import batch as BATCH
from phyloacc_lib import tree as TREE


def _base_globs(tmp_path):
    job_dir = tmp_path / "job-files"
    paths = {
        "job-dir": job_dir,
        "job-alns": job_dir / "alns",
        "job-cfgs": job_dir / "cfgs",
        "job-bed": job_dir / "bed",
        "job-out": job_dir / "phyloacc-output",
        "job-ids": job_dir / "ids",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)

    st = TREE.Tree("((sp1:0.1,sp2:0.1):0.1,sp3:0.2);")
    return {
        "theta": False,
        "batch": True,
        "alns": {
            "locus_st": {"sp1": "AAAA", "sp2": "AAAA", "sp3": "TTTT"},
            "locus_gt": {"sp1": "CCCC", "sp2": "CCCC", "sp3": "GGGG"},
        },
        "aln-stats": {
            "locus_st": {"batch-type": "st", "length": 4, "low-qual": False},
            "locus_gt": {"batch-type": "gt", "length": 4, "low-qual": False},
        },
        "no-inf-sites-loci": [],
        "filter-alns": False,
        "batch-size": 1,
        "st-batches": [],
        "gt-batches": [],
        "filtered-loci": 0,
        "st": st,
        "groups": {
            "targets": ["sp1"],
            "conserved": ["sp2"],
            "outgroup": ["sp3"],
        },
        "job-alns": str(paths["job-alns"]),
        "job-cfgs": str(paths["job-cfgs"]),
        "job-bed": str(paths["job-bed"]),
        "job-out": str(paths["job-out"]),
        "job-ids": str(paths["job-ids"]),
        "id-flag": True,
        "phyloacc-opts": ["SOME_OPT 7"],
        "dollo": True,
        "mod-file": str(tmp_path / "model.mod"),
        "coal-tree-file": str(tmp_path / "coal.tree"),
        "thin": "1",
        "burnin": "10",
        "mcmc": "20",
        "chain": "1",
        "procs-per-job": "2",
        "logfilename": str(tmp_path / "phyloacc.log"),
        "log-v": False,
        "test-cmd-flag": False,
    }


def test_gen_job_files_writes_expected_st_and_gt_configs(tmp_path, monkeypatch):
    globs = _base_globs(tmp_path)
    Path(globs["mod-file"]).write_text("dummy mod\n", encoding="utf-8")
    Path(globs["coal-tree-file"]).write_text("(sp1,sp2,sp3);\n", encoding="utf-8")

    monkeypatch.setattr(BATCH.PC, "report_step", lambda *args, **kwargs: False)
    monkeypatch.setattr(BATCH.PC, "printWrite", lambda *args, **kwargs: None)

    result = BATCH.genJobFiles(globs)

    assert result["st-batches"] == ["1"]
    assert result["gt-batches"] == ["2"]

    st_cfg = Path(result["job-cfgs"]) / "1-st.cfg"
    gt_cfg = Path(result["job-cfgs"]) / "2-gt.cfg"
    assert st_cfg.exists()
    assert gt_cfg.exists()

    st_text = st_cfg.read_text(encoding="utf-8")
    gt_text = gt_cfg.read_text(encoding="utf-8")

    assert "OUTGROUP sp3" in st_text
    assert "CONSERVE sp2" in st_text
    assert "ID_FILE" in st_text
    assert "TREE_IN_COALESCENT_UNIT" not in st_text
    assert "HYPER_LRATE2_A 0" in st_text
    assert "SOME_OPT 7" in st_text

    assert "TREE_IN_COALESCENT_UNIT" in gt_text
    assert "OUTGROUP sp3" in gt_text
    assert "ID_FILE" in gt_text

    st_ids = (Path(result["job-ids"]) / "1-st.id").read_text(encoding="utf-8")
    gt_ids = (Path(result["job-ids"]) / "2-gt.id").read_text(encoding="utf-8")
    assert st_ids == "0\n"
    assert gt_ids == "0\n"
