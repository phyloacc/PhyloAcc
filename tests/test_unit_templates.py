import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "PhyloAcc-interface"))

from phyloacc_lib import templates as TEMPLATES


def test_phyloacc_config_template_contains_required_fields():
    content = TEMPLATES.phyloaccConfig().format(
        mod_file="model.mod",
        bed_file="segments.bed",
        id_line="",
        aln_file="aln.fa",
        coal_tree_line="",
        outdir="out",
        batch="1",
        thin="1",
        burnin="10",
        mcmc="20",
        chain="1",
        targets="t1;t2",
        outgroup="",
        conserved="c1;c2",
        procs_per_job="1",
        phyloacc_opts="",
        dollo_str="",
    )

    required = [
        "PHYTREE_FILE",
        "SEG_FILE",
        "ALIGN_FILE",
        "RESULT_FOLDER",
        "PREFIX",
        "THIN",
        "BURNIN",
        "MCMC",
        "CHAIN",
        "TARGETSPECIES",
        "CONSERVE",
        "NUM_THREAD",
    ]
    for key in required:
        assert key in content
    assert "OUTGROUP" not in content


def test_phyloacc_config_template_includes_coal_tree_line():
    content = TEMPLATES.phyloaccConfig().format(
        mod_file="model.mod",
        bed_file="segments.bed",
        id_line="",
        aln_file="aln.fa",
        coal_tree_line="\nTREE_IN_COALESCENT_UNIT coal.tree",
        outdir="out",
        batch="1",
        thin="1",
        burnin="10",
        mcmc="20",
        chain="1",
        targets="t1;t2",
        outgroup="OUTGROUP o1;o2",
        conserved="c1;c2",
        procs_per_job="1",
        phyloacc_opts="",
        dollo_str="",
    )
    assert "TREE_IN_COALESCENT_UNIT coal.tree" in content
    assert "OUTGROUP o1;o2" in content
