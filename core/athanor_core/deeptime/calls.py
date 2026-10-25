"""Genotypes from ancient reads.

Pseudo-haploid calls. At low coverage a diploid genotype cannot be called, so
ancient-DNA studies take one allele at random from the reads covering a site
and use it twice (Skoglund et al., Science 344, 747, 2014). This keeps the
sample's allele frequencies unbiased, which is what the ancestry methods
need, and matches how the published panels of ancient individuals are built.
The random choice is seeded, so a rerun gives the same calls.

Diploid calls. Where coverage allows, genotypes are called from the base
qualities: for each of the ten genotypes, the probability of the observed
bases is computed with each read's own error probability, and the best is
taken when it is enough better than the next. This is the standard
genotype-likelihood model (Li, Bioinformatics 27, 2987, 2011), without
snpAD's position-dependent error profile (Prüfer, BMC Bioinformatics 19, 136,
2018), so a damaged library needs its read ends trimmed or its transitions
set aside instead.

Deamination and genotypes. A C->T or G->A difference can always be damage, so
every call here can be restricted to transversions, and bases near the ends
of molecules, where damage concentrates, can be dropped.
"""
from __future__ import annotations

import numpy as np

from .sites import BASES, CODE

TRANSITION = {(CODE["A"], CODE["G"]), (CODE["G"], CODE["A"]), (CODE["C"], CODE["T"]), (CODE["T"], CODE["C"])}


def _transition_mask(ref: np.ndarray, alt: np.ndarray) -> np.ndarray:
    return np.array([(int(r), int(a)) in TRANSITION for r, a in zip(ref, alt)])


def pseudo_haploid(p, ref: np.ndarray, alt: np.ndarray, seed: int = 1, min_qual: int = 20, trim: int = 0,
                   transversions_only: bool = False) -> dict:
    """One allele per site, drawn at random from the reads covering it.

    Returns the alternative-allele dosage used by the ancestry analyses: 0 or
    2 at called sites (the drawn allele counted twice), -1 where there is no
    call."""
    q = p.filtered(min_qual=min_qual, trim=trim)
    ref = np.asarray(ref, dtype=np.int8)
    alt = np.asarray(alt, dtype=np.int8)
    dos = np.full(p.n_sites, -1, dtype=np.int8)
    rng = np.random.default_rng(seed)
    # one random observation per site: shuffle, then keep the first of each site
    order = rng.permutation(len(q.site))
    site = q.site[order]
    base = q.base[order]
    first = np.full(p.n_sites, -1, dtype=np.int64)
    # later writes win, so reverse to keep the first of the shuffled order
    first[site[::-1]] = np.arange(len(site))[::-1]
    have = first >= 0
    idx = first[have]
    b = base[idx]
    s = np.nonzero(have)[0]
    is_ref, is_alt = b == ref[s], b == alt[s]
    keep = is_ref | is_alt
    if transversions_only:
        keep &= ~_transition_mask(ref[s], alt[s])
    dos[s[keep]] = np.where(is_alt[keep], 2, 0).astype(np.int8)
    counts = q.counts()
    depth = counts.sum(1)
    return {"dosage": dos, "sites_called": int((dos >= 0).sum()),
            "sites_covered": int((depth > 0).sum()),
            "sites_other_allele": int((have.sum() - keep.sum())),
            "mean_depth_at_covered": float(depth[depth > 0].mean()) if (depth > 0).any() else 0.0,
            "observations_used": int(q.reads_used),
            "filters": {"min_base_quality": min_qual, "trim_ends": trim,
                        "transversions_only": bool(transversions_only)},
            "method": "one allele drawn at random (seeded) from the reads covering each site and counted twice "
                      "(Skoglund et al. 2014)"}


def _genotype_loglik(counts: np.ndarray, err: float) -> tuple[np.ndarray, list]:
    """log P(observed bases | genotype) for the 10 genotypes, from per-site
    base counts and one error probability."""
    n = counts.shape[0]
    hom = np.log(1 - err)
    other = np.log(err / 3)
    het = np.log(0.5 * (1 - err) + 0.5 * err / 3)
    tot = counts.sum(1)
    out = np.empty((n, 10), dtype=np.float64)
    k = 0
    pairs = []
    for a in range(4):
        for b in range(a, 4):
            if a == b:
                out[:, k] = counts[:, a] * hom + (tot - counts[:, a]) * other
            else:
                m = counts[:, a] + counts[:, b]
                out[:, k] = m * het + (tot - m) * other
            pairs.append((a, b))
            k += 1
    return out, pairs


def diploid(counts: np.ndarray, min_depth: int = 8, max_depth: int | None = None, err: float = 0.005,
            min_quality: float = 30.0) -> dict:
    """Diploid genotypes from per-site base counts (n_sites x 4).

    Returns allele codes a1, a2 (-1 where no call) and the call's quality, the
    log-odds against the next-best genotype in decibans (as a Phred-scaled
    quality)."""
    counts = np.asarray(counts, dtype=np.int64)
    ll, pairs = _genotype_loglik(counts, err)
    depth = counts.sum(1)
    best = np.argmax(ll, axis=1)
    top = ll[np.arange(len(ll)), best]
    second = np.where(np.isneginf(top), -np.inf, np.partition(ll, -2, axis=1)[:, -2])
    qual = (top - second) * 10.0 / np.log(10.0)
    ok = (depth >= min_depth) & (qual >= min_quality)
    if max_depth:
        ok &= depth <= max_depth
    pa = np.array([p[0] for p in pairs], dtype=np.int8)
    pb = np.array([p[1] for p in pairs], dtype=np.int8)
    a1 = np.where(ok, pa[best], -1).astype(np.int8)
    a2 = np.where(ok, pb[best], -1).astype(np.int8)
    return {"a1": a1, "a2": a2, "quality": qual, "depth": depth, "called": ok,
            "method": f"genotype likelihoods with a flat prior, error probability {err}, depth at least "
                      f"{min_depth}{'' if not max_depth else f' and at most {max_depth}'}, and at least "
                      f"{min_quality} in favour of the called genotype over the next best"}


def panel_genotypes(p, ref: np.ndarray, alt: np.ndarray, min_qual: int = 20, trim: int = 0, **kw) -> dict:
    """Diploid alternative-allele dosage (0, 1, 2; -1 no call) at panel sites."""
    q = p.filtered(min_qual=min_qual, trim=trim)
    g = diploid(q.counts(), **kw)
    ref = np.asarray(ref, dtype=np.int8)
    alt = np.asarray(alt, dtype=np.int8)
    a1, a2 = g["a1"], g["a2"]
    known = g["called"] & np.isin(a1, [0, 1, 2, 3]) & ((a1 == ref) | (a1 == alt)) & ((a2 == ref) | (a2 == alt))
    dos = np.full(len(ref), -1, dtype=np.int8)
    dos[known] = ((a1[known] == alt[known]).astype(np.int8) + (a2[known] == alt[known]).astype(np.int8))
    return {"dosage": dos, "quality": g["quality"], "depth": g["depth"],
            "sites_called": int(known.sum()), "sites_third_allele": int((g["called"] & ~known).sum()),
            "method": g["method"]}


def haploid_calls(genome, chrom: str, positions: np.ndarray, min_depth: int = 2, min_fraction: float = 0.7,
                  min_qual: int = 20, trim: int = 0, min_mapq: int | None = None) -> tuple[dict, dict]:
    """Haploid bases (one copy only: Y or mitochondrion) at the given
    positions: position -> base, plus how many sites were covered and called."""
    from .sites import pileup
    p = pileup(genome, chrom, positions, min_mapq=min_mapq)
    q = p.filtered(min_qual=min_qual, trim=trim)
    counts = q.counts()
    depth = counts.sum(1)
    top = counts.max(1)
    which = counts.argmax(1)
    ok = (depth >= min_depth) & (top >= min_fraction * np.maximum(depth, 1))
    out = {int(pos): BASES[int(b)] for pos, b, good in zip(positions, which, ok) if good}
    return out, {"positions": len(positions), "covered": int((depth > 0).sum()), "called": len(out),
                 "mean_depth_at_covered": float(depth[depth > 0].mean()) if (depth > 0).any() else 0.0,
                 "rule": f"at least {min_depth} reads and at least {min_fraction:.0%} of them agreeing"}


def region_counts(genome, chrom: str, start: int, end: int, min_qual: int = 25,
                  min_mapq: int | None = None) -> np.ndarray:
    """Per-position base counts (A, C, G, T) over a region, 0-based half-open."""
    ref = genome.name_of(chrom)
    mq = genome.min_mapq if min_mapq is None else int(min_mapq)

    def keep(r):
        return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate) and r.mapping_quality >= mq

    cov = genome.bam.count_coverage(ref, start, end, quality_threshold=min_qual, read_callback=keep)
    return np.array(cov, dtype=np.int32).T          # positions x 4


def heterozygosity(counts: np.ndarray, callable_mask: np.ndarray | None = None, min_depth: int = 8,
                   max_depth: int | None = None, err: float = 0.005, min_quality: float = 30.0) -> dict:
    """Heterozygous sites per callable base, from per-position base counts.

    No reference sequence is needed: a site is heterozygous when the reads
    carry two alleles. Transitions are reported separately because C->T and
    G->A deamination mimics them, so the transversion rate is the one that
    survives damage; for an undamaged or enzyme-treated library the total is
    the comparable number."""
    g = diploid(counts, min_depth=min_depth, max_depth=max_depth, err=err, min_quality=min_quality)
    called = g["called"].copy()
    if callable_mask is not None:
        called &= np.asarray(callable_mask, dtype=bool)
    a1, a2 = g["a1"], g["a2"]
    het = called & (a1 != a2)
    ts = np.zeros(len(a1), dtype=bool)
    hi = np.nonzero(het)[0]
    if len(hi):
        ts[hi] = [(int(a1[i]), int(a2[i])) in TRANSITION for i in hi]
    n = int(called.sum())
    nh, nts = int(het.sum()), int((het & ts).sum())
    ntv = nh - nts
    return {"callable_sites": n, "heterozygous_sites": nh, "transitions": nts, "transversions": ntv,
            "heterozygosity": (nh / n) if n else None,
            "transversion_heterozygosity": (ntv / n) if n else None,
            "transition_share": (nts / nh) if nh else None, "method": g["method"]}


def bed_mask(path: str, chrom: str, start: int, end: int) -> np.ndarray:
    """Boolean mask over [start, end) from a BED file's intervals on ``chrom``."""
    import gzip
    from .bam import plain
    mask = np.zeros(end - start, dtype=bool)
    want = plain(chrom)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line or line[0] in "#tb":
                continue
            f = line.split()
            if len(f) < 3 or plain(f[0]) != want:
                continue
            a, b = int(f[1]), int(f[2])
            if b <= start or a >= end:
                continue
            mask[max(0, a - start):min(end - start, b - start)] = True
    return mask
