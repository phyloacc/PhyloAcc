import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "PhyloAcc-interface"))

from phyloacc_lib import tree as TREE


def test_categorize_branches_propagates_tip_groups_to_internal_nodes():
    tree_str = (ROOT / "tests" / "data" / "unit" / "tree_groups.nwk").read_text(encoding="utf-8").strip()
    st = TREE.Tree(tree_str)
    globs = {
        "groups": {
            "targets": ["sp1", "sp2"],
            "conserved": ["sp3"],
            "outgroup": ["sp4"],
        }
    }

    targets, conserved, outgroups = TREE.categorizeBranches(globs, st)

    assert st.root in outgroups
    target_internal = next(node for node, desc in st.desc.items() if set(desc) == {"sp1", "sp2"})
    assert target_internal in targets
    parent_of_target = next(node for node, desc in st.desc.items() if target_internal in desc and "sp3" in desc)
    assert parent_of_target in conserved
