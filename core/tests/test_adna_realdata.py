"""Real-data helpers: circular folding, duplicates, indels, comparison, statistics."""
import numpy as np
import pytest

from athanor_core.adna.bamtools import apply_changes, majority_indels
from athanor_core.adna.compare import parse_cs
from athanor_core.adna.consensus import Observations, majority_consensus
from athanor_core.adna.damage import damage_from_observations, wilson_interval

def _write_bam(path, length, reads):
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "ref", "LN": length}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for d in reads:
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence = d["name"], d["seq"]
            a.reference_id, a.reference_start = 0, d["start"]
            a.flag = 16 if d.get("rev") else 0
            a.mapping_quality = d.get("mapq", 37)
            a.cigarstring = d.get("cigar", f"{len(d['seq'])}M")
            a.query_qualities = pysam.qualitystring_to_array(d.get("qual", "I" * len(d["seq"])))
            a.set_tags(d.get("tags", []))
            out.write(a)
    pysam.index(str(path))


def test_clean_bam_folds_rescues_and_deduplicates(tmp_path):
    pysam = pytest.importorskip("pysam")
    from athanor_core.adna.bamtools import clean_bam
    L = 1000
    reads = [
        # a read in the copied stretch, placed at L + 10 with its twin at 10: unique on the circle
        {"name": "0.1:f", "seq": "A" * 30, "start": L + 10, "mapq": 0,
         "tags": [("X0", 2), ("X1", 0), ("XA", "ref,+11,30M,0;")]},
        # same, but with another hit elsewhere: stays ambiguous
        {"name": "0.2:f", "seq": "C" * 30, "start": 50, "mapq": 0,
         "tags": [("X0", 2), ("X1", 0), ("XA", "ref,+501,30M,0;")]},
        # whole-molecule duplicates (same strand, start, end): the better copy is kept
        {"name": "0.3:f", "seq": "G" * 40, "start": 200, "qual": "5" * 40},
        {"name": "0.4:f", "seq": "G" * 40, "start": 200},
        # same coordinates but from another library: not a duplicate
        {"name": "1.5:f", "seq": "G" * 40, "start": 200},
        # partial reads (3' end not observed): duplicates by strand and 5' end
        {"name": "0.6:p", "seq": "T" * 60, "start": 300},
        {"name": "0.7:p", "seq": "T" * 55, "start": 300, "cigar": "55M"},
        # spans the origin: kept in elongated coordinates
        {"name": "0.8:f", "seq": "A" * 40, "start": L - 20},
        {"name": "0.9:f", "seq": "A" * 40, "start": 600, "mapq": 10},
    ]
    _write_bam(tmp_path / "in.bam", L + 300, sorted(reads, key=lambda r: r["start"]))
    st = clean_bam(tmp_path / "in.bam", tmp_path / "out.bam", 25, {"ref": L})
    assert st["rescued_circular"] == 1
    assert st["low_mapq"] == 2           # the ambiguous read and the MAPQ-10 read
    assert st["duplicates"] == 2
    with pysam.AlignmentFile(str(tmp_path / "out.bam")) as b:
        kept = {r.query_name: r for r in b.fetch(until_eof=True)}
    assert set(kept) == {"0.1:f", "0.4:f", "1.5:f", "0.6:p", "0.8:f"}
    assert kept["0.1:f"].reference_start == 10 and kept["0.1:f"].mapping_quality == 37
    assert kept["0.8:f"].reference_start == L - 20


def test_bam_to_fastq_restores_sequencing_orientation(tmp_path):
    pytest.importorskip("pysam")
    from athanor_core.adna.reads import bam_to_fastq
    _write_bam(tmp_path / "r.bam", 100, [
        {"name": "f", "seq": "AACG", "start": 5, "qual": "ABCD"},
        {"name": "r", "seq": "AACG", "start": 10, "rev": True, "qual": "ABCD"}])
    out = bam_to_fastq(tmp_path / "r.bam", tmp_path / "fq")
    lines = open(out["paths"]["single"]).read().split()
    rec = {lines[i][1:]: (lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)}
    assert rec["f"] == ("AACG", "ABCD")
    assert rec["r"] == ("CGTT", "DCBA")
    assert out["all_reads_mapped"]


def test_majority_indels_and_apply():
    from collections import Counter
    span = np.full(100, 10)
    dels = Counter({(20, 2): 7, (40, 1): 2})
    ins = Counter({(60, "GA"): 6, (61, "T"): 5})       # 5 of 10 is not a majority
    ev = majority_indels(dels, ins, span, min_reads=3)
    assert [(e["type"], e["pos"]) for e in ev] == [("deletion", 20), ("insertion", 60)]
    seq = list("ACGT" * 25)
    new, org = apply_changes(seq, list(range(100)), {5: "T"}, ev)
    assert len(new) == 100 - 2 + 2
    assert new[5] == "T"
    assert org[20] == 22                       # bases 20, 21 deleted
    i = 60 - 2                                 # insertion sits before original base 60
    assert new[i:i + 2] == ["G", "A"] and org[i:i + 2] == [-1, -1] and org[i + 2] == 60


def test_parse_cs_classifies_differences():
    r = parse_cs(":10*ag:5*na*an+tt:3-c:2", 0, 0)
    assert r["matches"] == 20
    assert r["substitutions"] == [{"theirs_pos": 10, "ours_pos": 10, "theirs": "A", "ours": "G"}]
    assert r["uncalled_theirs"] == 1 and r["uncalled_ours"] == 1
    assert r["insertions_in_ours"] == [{"theirs_pos": 18, "ours_pos": 18, "sequence": "TT"}]
    assert r["deletions_in_ours"] == [{"theirs_pos": 21, "ours_pos": 23, "sequence": "C"}]


def _obs(sites, bases, d5, d3, rev=None, full=None):
    n = len(sites)
    return Observations(site=np.array(sites), base=np.array(bases), qual=np.full(n, 30.0),
                        dist5=np.array(d5), dist3=np.array(d3),
                        reverse=np.zeros(n, bool) if rev is None else np.array(rev),
                        full=None if full is None else np.array(full))


def test_majority_consensus_trims_ends_and_needs_three_reads():
    # site 0: three reads say C, but one at a read end (trimmed) -> 2 reads -> N
    # site 1: A, A, A, G -> A ; site 2: A, A, G, G -> no majority -> N
    o = _obs([0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2], [1, 1, 1, 0, 0, 0, 2, 0, 0, 2, 2],
             [0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5], [5] * 11)
    seq, counts, depth = majority_consensus(o, 3)
    assert seq == "NAN"
    assert depth.tolist() == [2, 4, 4]


def test_damage_from_observations_ignores_cut_3p_ends():
    # true base C everywhere; 5' first base reads T in 1 of 4; 3' G->A measured only on whole molecules
    codes = np.array([1, 2])
    o = _obs([0, 0, 0, 0, 1, 1, 1], [3, 1, 1, 1, 0, 0, 2], [0, 0, 0, 0, 9, 9, 9], [9, 9, 9, 9, 0, 0, 0],
             full=[True] * 5 + [False, True])
    p = damage_from_observations(o, codes, max_pos=10, background_positions=2)
    assert p["raw_5p_CtoT"][0] == pytest.approx(0.25) and p["n_5p"][0] == 4
    assert p["raw_3p_GtoA"][0] == pytest.approx(0.5) and p["n_3p"][0] == 2


def test_intervals():
    from athanor_core.adna.contamination import clopper_pearson
    lo, hi = clopper_pearson(0, 100)
    assert lo == 0 and hi == pytest.approx(0.0362, abs=2e-4)       # exact binomial, 95%
    lo, hi = clopper_pearson(5, 100)
    assert lo == pytest.approx(0.0164, abs=2e-4) and hi == pytest.approx(0.1128, abs=2e-4)
    lo, hi = wilson_interval(5, 100)
    assert lo == pytest.approx(0.0215, abs=2e-4) and hi == pytest.approx(0.1118, abs=2e-4)
