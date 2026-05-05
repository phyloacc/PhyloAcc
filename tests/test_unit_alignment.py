import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "PhyloAcc-interface"))

from phyloacc_lib import seq as SEQ


class _StubTree:
    def __init__(self, tips):
        self.tips = tips


def test_read_bed_without_id_column_assigns_sequential_ids(tmp_path):
    bed = tmp_path / "segments.bed"
    bed.write_text("chr1\t0\t2\nchr1\t2\t4\n", encoding="utf-8")

    result = SEQ.readBed(str(bed), {"bed-compression": "none", "id-file": False, "locus-ids": []})

    assert list(result.keys()) == ["1", "2"]
    assert result["1"] == {"scaff": "chr1", "start": 0, "end": 2}
    assert result["2"] == {"scaff": "chr1", "start": 2, "end": 4}


def test_read_bed_filters_by_id_file(tmp_path):
    bed = tmp_path / "segments.bed"
    bed.write_text(
        "chr1\t0\t2\tlocus1\nchr1\t2\t4\tlocus2\nchr1\t4\t6\tlocus3\n",
        encoding="utf-8",
    )

    result = SEQ.readBed(
        str(bed),
        {"bed-compression": "none", "id-file": True, "locus-ids": ["locus2", "locus3"]},
    )

    assert list(result.keys()) == ["locus2", "locus3"]
    assert result["locus2"]["start"] == 2
    assert result["locus3"]["end"] == 6


def test_partition_seqs_splits_concat_alignment():
    concat = {"sp1": "AACCGG", "sp2": "TTGGCC"}
    bed = {
        "l1": {"scaff": "chr1", "start": 0, "end": 2},
        "l2": {"scaff": "chr1", "start": 2, "end": 6},
    }

    result = SEQ.partitionSeqs(concat, bed)

    assert result["l1"] == {"sp1": "AA", "sp2": "TT"}
    assert result["l2"] == {"sp1": "CCGG", "sp2": "GGCC"}


def test_locus_aln_stats_counts():
    aln = {
        "sp1": "A-CG",
        "sp2": "ATCG",
        "sp3": "ATGG",
        "sp4": "ATGG",
    }
    locus, stats = SEQ.locusAlnStats(("l1", aln, ["-", "N", "n"]))

    assert locus == "l1"
    assert stats["num-seqs"] == 4
    assert stats["length"] == 4
    assert stats["variable-sites"] == 1
    assert stats["informative-sites"] == 1
    assert stats["num-sites-w-gap"] == 1
    assert stats["num-sites-half-gap"] == 0
    assert stats["num-seqs-half-missing"] == 0
    assert stats["num-seqs-all-missing"] == 0
    assert stats["low-qual"] is False
    assert stats["avg-nogap-seq-len"] == pytest.approx(3.75)


def test_locus_aln_stats_marks_low_quality_for_missing_sequences():
    aln = {
        "sp1": "AAAA",
        "sp2": "----",
        "sp3": "----",
        "sp4": "AAAA",
    }

    _, stats = SEQ.locusAlnStats(("l2", aln, ["-", "N", "n"]))

    assert stats["num-seqs-half-missing"] == 2
    assert stats["num-seqs-all-missing"] == 2
    assert stats["num-sites-half-gap"] == 4
    assert stats["low-qual"] is True


def test_check_aln_labels_mismatch_raises():
    globs = {
        "alns": {"l1": {"sp1": "AA", "sp2": "AA", "sp3": "AA"}},
        "aln-file": "dummy.fa",
        "st": _StubTree(["sp1", "sp2"]),
        "logfilename": "/tmp/phyloacc-test.log",
        "log-v": False,
        "endprog": False,
    }

    with pytest.raises(SystemExit):
        SEQ.checkAlnLabels(globs)


def test_check_aln_labels_unknown_header_raises():
    globs = {
        "alns": {"l1": {"sp1": "AA", "spX": "AA"}},
        "aln-file": False,
        "st": _StubTree(["sp1", "sp2"]),
        "logfilename": "/tmp/phyloacc-test.log",
        "log-v": False,
        "endprog": False,
    }

    with pytest.raises(SystemExit):
        SEQ.checkAlnLabels(globs)


def test_read_bed_invalid_coordinates_raise_value_error(tmp_path):
    bed = tmp_path / "bad_segments.bed"
    bed.write_text("chr1\tstart\t10\tlocus1\n", encoding="utf-8")

    with pytest.raises(ValueError):
        SEQ.readBed(str(bed), {"bed-compression": "none", "id-file": False, "locus-ids": []})


def test_partition_seqs_allows_empty_locus_slice():
    concat = {"sp1": "AACCGG", "sp2": "TTGGCC"}
    bed = {"empty": {"scaff": "chr1", "start": 3, "end": 3}}

    result = SEQ.partitionSeqs(concat, bed)

    assert result["empty"] == {"sp1": "", "sp2": ""}
