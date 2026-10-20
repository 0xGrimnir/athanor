"""BAM -> observations: strand-aware distances from the molecule ends."""
import numpy as np
import pytest

pysam = pytest.importorskip("pysam")

from athanor_core.adna.pileup import observations_from_bam  # noqa: E402


def _bam(tmp_path, reads):
    path = str(tmp_path / "t.bam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrM", "LN": 100}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for name, start, seq, rev, cigar in reads:
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence, a.reference_id, a.reference_start = name, seq, 0, start
            a.flag = 16 if rev else 0
            a.mapping_quality = 37
            a.cigarstring = cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            out.write(a)
    pysam.index(path)
    return path


def test_distances_follow_the_molecule(tmp_path):
    bam = _bam(tmp_path, [("f", 10, "ACGTA", False, "5M"), ("r", 20, "TTGCA", True, "5M"),
                          ("c", 40, "GGACGT", False, "2S4M")])
    o = observations_from_bam(bam, min_mapq=0)
    f = o.site < 15
    assert list(o.site[f]) == [10, 11, 12, 13, 14]
    assert list(o.dist5[f]) == [0, 1, 2, 3, 4] and list(o.dist3[f]) == [4, 3, 2, 1, 0]
    r = (o.site >= 20) & (o.site < 25)
    # reverse strand: the molecule's 5' end is the right-most reference base
    assert list(o.dist5[r]) == [4, 3, 2, 1, 0] and all(o.reverse[r])
    c = o.site >= 40
    # soft-clipped bases are not observations, but still count toward the distance
    assert list(o.site[c]) == [40, 41, 42, 43] and list(o.dist5[c]) == [2, 3, 4, 5]
    assert np.all(o.qual == 40)
