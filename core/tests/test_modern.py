"""Modern ancestry and trait reports: estimators on simulated data with a
known answer, and the trait genotyper's strand handling."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from athanor_core.personal.modern import (admix, composition, decode_freq, encode_freq, paint, roh, roh_summary,
                                          segments)
from athanor_core.personal.rawfile import CODE, MISSING, RawGenotypes
from athanor_core.personal.traits import Genotyper, trait_cards

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "pgdata"))


def _panels(K, M, seed=0, spread=30):
    rng = np.random.default_rng(seed)
    anc = rng.beta(0.5, 0.5, M)
    return rng, np.clip(np.array([rng.beta(anc * spread + 0.1, (1 - anc) * spread + 0.1) for _ in range(K)]), 0.001, 0.999)


def test_frequency_encoding_error_is_small():
    p = np.linspace(0, 1, 10001)
    d = decode_freq(encode_freq(p))
    assert np.abs(d - np.clip(p, 1e-4, 1 - 1e-4)).max() < 0.0035
    # finer near the edges, where rare alleles live
    small = (p >= 0.002) & (p < 0.01)
    assert (np.abs(d[small] - p[small]) / p[small]).max() < 0.2
    assert np.abs(d[small] - p[small]).max() < 0.0007


def test_admix_recovers_simulated_proportions():
    rng, P = _panels(12, 120_000)
    q = np.zeros(12); q[[1, 4, 9]] = [0.55, 0.35, 0.10]
    g = rng.binomial(2, q @ P)
    est = admix(g, P)["q"]
    assert np.abs(est - q).max() < 0.02
    assert est[[0, 2, 3, 5, 6, 7, 8, 10, 11]].sum() < 0.03


def test_admix_haploid_x():
    rng, P = _panels(4, 40_000, seed=3)
    q = np.array([0.7, 0.3, 0, 0])
    g = rng.binomial(1, q @ P)
    est = admix(g, P, n=np.ones(len(g)))["q"]
    assert abs(est[0] - 0.7) < 0.04 and abs(est[1] - 0.3) < 0.04


def test_bootstrap_interval_contains_truth():
    rng, P = _panels(5, 60_000, seed=5)
    q = np.array([0.6, 0.4, 0, 0, 0])
    g = rng.binomial(2, q @ P).astype(np.int8)
    r = composition(g, P, np.arange(60_000) // 300, reps=30)
    lo, hi = np.percentile(r["boot"][:, 0], [2.5, 97.5])
    assert lo <= 0.6 <= hi and hi - lo < 0.05


def test_painting_follows_known_segments():
    rng, P = _panels(2, 30_000, seed=7, spread=8)
    chrom = np.repeat([1, 2], 15_000)
    cm = np.r_[np.linspace(0, 150, 15_000), np.linspace(0, 150, 15_000)]
    t = np.zeros((2, 30_000), int)
    for h in range(2):
        for c in (1, 2):
            idx = np.flatnonzero(chrom == c)
            cuts = np.sort(rng.uniform(0, 150, 10))
            src = rng.integers(0, 2, 11)
            t[h, idx] = src[np.searchsorted(cuts, cm[idx])]
    p0, p1 = P[t[0], np.arange(30_000)], P[t[1], np.arange(30_000)]
    g = (rng.random(30_000) < p0).astype(int) + (rng.random(30_000) < p1).astype(int)
    r = paint(g, P, chrom, cm, np.array([0.5, 0.5]), generations=8, temper=1.0)
    a, b = np.minimum(t[0], t[1]), np.maximum(t[0], t[1])
    acc = ((r["a"] == a) & (r["b"] == b)).mean()
    assert acc > 0.93
    assert abs(acc - r["posterior"].mean()) < 0.05
    segs = segments(chrom, np.arange(30_000), r["a"], r["b"], r["posterior"])
    assert 10 < len(segs) < 120


def test_roh_finds_inserted_run():
    rng = np.random.default_rng(11)
    n = 20_000
    chrom = np.ones(n, int)
    pos = np.cumsum(rng.integers(3000, 9000, n))
    het = rng.random(n) < 0.3
    s0, s1 = 8000, 9500                          # about 9 Mb of homozygous calls
    het[s0:s1] = False
    het[s0 + 700] = True                         # one heterozygous call inside, allowed by the window rule
    miss = rng.random(n) < 0.01
    segs = roh(chrom, pos, het & ~miss, miss)
    hit = [s for s in segs if s["start"] <= pos[s0 + 60] and s["end"] >= pos[s1 - 60]]
    assert len(hit) == 1
    summ = roh_summary(segs)
    assert summ["by_length_mb"]["5-10"] > 5
    # nothing long elsewhere
    assert sum(1 for s in segs if s["kb"] > 2000) == 1


def _raw(rows):
    """rows: (rsid, chrom, pos, a1, a2)"""
    return RawGenotypes("23andMe", "t.txt", "x", np.array([r[0] for r in rows], np.int64),
                        np.array([r[1] for r in rows], np.int8), np.array([r[2] for r in rows], np.int64),
                        np.array([CODE.get(r[3], MISSING) for r in rows], np.uint8),
                        np.array([CODE.get(r[4], MISSING) for r in rows], np.uint8))


def test_genotyper_strands():
    rec = {"rsid": "rs671", "chrom": 12, "pos37": 112241766, "ref": "G", "alts": ["A"], "effect": "A", "other": "G"}
    amb = {"rsid": "rs713598", "chrom": 7, "pos37": 141673345, "ref": "C", "alts": ["G"], "effect": "G", "other": "C"}
    g = _raw([(671, 12, 112241766, "A", "G"), (713598, 7, 141673345, "G", "G")])
    gt = Genotyper(g, g.pos, "the reference's + strand")
    assert gt.call(rec)[0] == 1
    assert gt.call(amb)[0] == 2
    # the same file reported on the other strand
    g2 = _raw([(671, 12, 112241766, "T", "C"), (713598, 7, 141673345, "C", "C")])
    gt2 = Genotyper(g2, g2.pos, "the opposite strand")
    assert gt2.call(rec)[0] == 1 and gt2.call(amb)[0] == 2
    # strand unknown: a non-ambiguous SNP is still readable, an A/T or C/G one is not
    gt3 = Genotyper(g2, g2.pos, "mixed or undetermined")
    assert gt3.call(rec)[0] == 1 and gt3.call(amb)[0] is None
    # not on the chip, and found by rs number when the position differs
    assert Genotyper(g, g.pos, "the reference's + strand").call({**rec, "pos37": 1, "chrom": 1})[0] == 1
    assert gt.call({**rec, "rsid": "rs1", "pos37": 5})[0] is None


def test_trait_cards_cover_every_variant():
    from traits_spec import REFS, VARIANTS
    rng = np.random.default_rng(1)
    variants, rows = {}, []
    for k, (rs, spec) in enumerate(VARIANTS.items()):
        e = spec["effect"] or "A"
        o = {"A": "G", "G": "A", "C": "T", "T": "C"}[e]
        variants[rs] = {"rsid": rs, "gene": spec["gene"], "label": spec["label"], "chrom": 1 + k % 22, "pos37": 1000 + k,
                        "ref": o, "alts": [e], "effect": e, "other": o, "freq": {"EUR": 0.3, "EAS": 0.1},
                        "gwas": [{"trait": "t", "first_author": "X", "p": 1e-9, "or": 1.4, "accession": "GCST1",
                                  "pubmed_id": "1"}] if spec.get("gwas") else []}
        n = int(rng.integers(0, 3))
        rows.append((int(rs[2:]), 1 + k % 22, 1000 + k, e if n >= 1 else o, e if n == 2 else o))
    traits = {"variants": variants, "references": {k: {"first_author": "A", "title": "T", "journal": "J", "year": "2000",
                                                       "url": "https://doi.org/x"} for k in REFS}}
    g = _raw(rows)

    def make(cid, title, group, grade, summary, result, limits, method, sources, ok=True):
        return {"id": cid, "group": group, "grade": grade, "summary": summary, "result": result, "ok": ok}

    cards = trait_cards(traits, Genotyper(g, g.pos, "the reference's + strand"), make)
    assert len(cards) == 18
    seen = {v["rsid"] for c in cards for v in c["result"]["variants"]}
    assert seen == set(VARIANTS)
    assert all(v["copies"] is not None for c in cards for v in c["result"]["variants"])
    assert all(c["summary"] and c["grade"] in "ABC" for c in cards)


def test_panel_size_does_not_decide_between_identical_panels():
    """Two panels drawn from one population, of different sizes, at SNPs that
    are mostly monomorphic in it (as most SNPs are within any one region): a
    fixed pseudo-count pulls a new member strongly towards the larger panel;
    pooled shrinkage does not."""
    from athanor_core.personal.modern import panel_freqs
    rng = np.random.default_rng(4)
    M = 80_000
    p = np.where(rng.random(M) < 0.6, 0.0, rng.beta(0.4, 0.4, M))     # 60% of SNPs monomorphic here
    a1, a2 = rng.binomial(22, p), rng.binomial(60, p)                   # panels of 11 and 30 people
    t1, t2 = np.full(M, 22), np.full(M, 60)
    g = rng.binomial(2, p)
    fixed = np.array([(a1 + 0.5) / (t1 + 1), (a2 + 0.5) / (t2 + 1)])
    shrunk = panel_freqs(np.array([a1, a2]), np.array([t1, t2]))
    q_fixed = admix(g, fixed)["q"]
    q_shrunk = admix(g, shrunk)["q"]
    assert q_fixed[1] > 0.8                      # the artefact the shrinkage removes
    assert 0.25 < q_shrunk[0] < 0.75
