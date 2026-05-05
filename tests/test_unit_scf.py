import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "PhyloAcc-interface"))

from phyloacc_lib import cf as CF

SCAF_CASES = json.loads((ROOT / "tests" / "data" / "unit" / "scf_cases.json").read_text(encoding="utf-8"))


def _case(name):
    case = SCAF_CASES[name]
    quartet = tuple(tuple(side) for side in case["quartet"])
    return case["aln"], quartet


def test_locus_scf_simple_quartet():
    aln, quartet = _case("simple_quartet")
    quartets = {"N1": [quartet]}
    st_nodes = ["N1"]
    skip_chars = ["-", "N", "n"]

    node_scf, locus, quartet_scores = CF.locusSCF(
        ("locus1", aln, quartets, st_nodes, skip_chars, "zip")
    )

    assert locus == "locus1"
    assert len(node_scf) == 1
    assert node_scf[0] == 0.5

    scores = quartet_scores["N1"][quartet]
    assert scores["decisive-sites"] == 2
    assert scores["concordant-sites"] == 1
    assert scores["scf"] == 0.5


def test_locus_scf_counts_discordant_sites_and_skips_missing():
    aln, quartet = _case("discordant_and_skip")

    node_scf, _, quartet_scores = CF.locusSCF(
        ("locus2", aln, {"N1": [quartet]}, ["N1"], ["-", "N", "n"], "zip")
    )

    assert len(node_scf) == 1
    scores = quartet_scores["N1"][quartet]
    assert scores["variable-sites"] == 3
    assert scores["decisive-sites"] == 2
    assert scores["concordant-sites"] == 1
    assert scores["disco1-sites"] == 1
    assert scores["disco2-sites"] == 0
    assert scores["scf"] == pytest.approx(0.5)


def test_locus_scf_zip_and_loop_match():
    aln, quartet = _case("discordant_and_skip")
    args = ("locus3", aln, {"N1": [quartet]}, ["N1"], ["-", "N", "n"])

    zip_result = CF.locusSCF(args + ("zip",))
    loop_result = CF.locusSCF(args + ("loop",))

    assert zip_result[0] == loop_result[0]
    assert zip_result[2]["N1"][quartet] == loop_result[2]["N1"][quartet]


def test_locus_scf_returns_na_when_no_decisive_sites_exist():
    aln, quartet = _case("no_decisive")

    node_scf, _, quartet_scores = CF.locusSCF(
        ("locus4", aln, {"N1": [quartet]}, ["N1"], ["-", "N", "n"], "zip")
    )

    assert node_scf == []
    scores = quartet_scores["N1"][quartet]
    assert scores["variable-sites"] == 1
    assert scores["decisive-sites"] == 0
    assert scores["scf"] == "NA"


def test_locus_scf_skips_quartets_with_missing_species():
    aln, quartet = _case("missing_species")

    node_scf, _, quartet_scores = CF.locusSCF(
        ("locus5", aln, {"N1": [quartet]}, ["N1"], ["-", "N", "n"], "zip")
    )

    assert node_scf == []
    assert quartet_scores["N1"][quartet]["decisive-sites"] == 0
    assert quartet_scores["N1"][quartet]["scf"] == "NA"
