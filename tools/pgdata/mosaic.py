"""Simulated admixed genomes with known ancestry segments, built from real
reference genomes. Used to calibrate and validate composition and painting.

Each of the two copies of every chromosome is cut at recombination points
(a Poisson process, `gens` switches per Morgan) and each piece is assigned a
source drawn from `props`. Where both copies come from the same source, that
source individual's own genotype is used (so linkage and heterozygosity are
real); where they differ, one allele is drawn from each source's genotype.
"""
from __future__ import annotations

import numpy as np


def tracks(chrom: np.ndarray, cm: np.ndarray, props: np.ndarray, gens: float, rng) -> np.ndarray:
    """(2, M) source index per copy per SNP."""
    out = np.zeros((2, len(cm)), np.int16)
    K = len(props)
    for h in range(2):
        for c in np.unique(chrom):
            idx = np.flatnonzero(chrom == c)
            x = cm[idx]
            lo, hi = float(x.min()), float(x.max())
            mean_gap = 100.0 / max(gens, 1e-9)
            n = int((hi - lo) / mean_gap * 3) + 10
            cuts = lo + np.cumsum(rng.exponential(mean_gap, n))
            cuts = cuts[cuts < hi]
            src = rng.choice(K, size=len(cuts) + 1, p=props)
            out[h, idx] = src[np.searchsorted(cuts, x, side="right")]
    return out


def build(sources: list[np.ndarray], props, chrom, cm, gens: float, seed: int):
    """sources: per-source alt dosage arrays (-1 missing). Returns
    (dosage, truth_a, truth_b) with truth_a <= truth_b source indices."""
    rng = np.random.default_rng(seed)
    props = np.asarray(props, float)
    t = tracks(chrom, cm, props / props.sum(), gens, rng)
    a, b = np.minimum(t[0], t[1]), np.maximum(t[0], t[1])
    M = len(cm)
    stack = np.vstack(sources)
    g = np.full(M, -1, np.int8)
    same = a == b
    g[same] = stack[a[same], np.flatnonzero(same)]

    def draw(d):
        return np.where(d < 0, -9, np.where(d == 2, 1, np.where(d == 0, 0, rng.integers(0, 2, len(d)))))
    diff = np.flatnonzero(~same)
    s = draw(stack[a[diff], diff]) + draw(stack[b[diff], diff])
    g[diff] = np.where(s < 0, -1, s).astype(np.int8)
    return g, a, b


def truth_props(a, b, k: int) -> np.ndarray:
    return np.array([((a == j).sum() + (b == j).sum()) / (2 * len(a)) for j in range(k)])
