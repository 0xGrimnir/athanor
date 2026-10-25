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

A quality floor is also a depth floor, and a stricter one for homozygotes
than for heterozygotes: each read of a clean homozygote adds the same fixed
evidence, so at shallow depth homozygous sites fail the floor while
heterozygous ones pass it. Leaving the two inconsistent would drop
homozygotes from the denominator of every rate while keeping heterozygotes in
the numerator, so the depth floor is raised to where a homozygote can reach
the quality asked for.

Deamination and genotypes. A C->T or G->A difference can always be damage, so
every call here can be restricted to transversions, and bases near the ends
of molecules, where damage concentrates, can be dropped.
"""
from __future__ import annotations

import math

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
    is_ref_all = np.zeros(p.n_sites, dtype=bool); is_ref_all[s] = is_ref
    is_alt_all = np.zeros(p.n_sites, dtype=bool); is_alt_all[s] = is_alt
    if transversions_only:
        keep &= ~_transition_mask(ref[s], alt[s])
    keep_all = np.zeros(p.n_sites, dtype=bool); keep_all[s] = keep
    dos[s[keep]] = np.where(is_alt[keep], 2, 0).astype(np.int8)
    counts = q.counts()
    depth = counts.sum(1)
    return {"dosage": dos, "sites_called": int((dos >= 0).sum()),
            "sites_covered": int((depth > 0).sum()),
            "sites_with_an_unexpected_base": int((have & ~(is_ref_all | is_alt_all)).sum()),
            "sites_set_aside_as_transitions": int(((is_ref_all | is_alt_all) & ~keep_all).sum()),
            "mean_depth_at_covered": float(depth[depth > 0].mean()) if (depth > 0).any() else 0.0,
            "observations_used": int(q.reads_used),
            "filters": {"min_base_quality": min_qual, "trim_ends": trim,
                        "transversions_only": bool(transversions_only)},
            "method": "one allele drawn at random (seeded) from the reads covering each site and counted twice "
                      "(Skoglund et al. 2014)"}


GENOTYPES = [(a, b) for a in range(4) for b in range(a, 4)]


def _genotype_loglik(counts: np.ndarray, err: float) -> tuple[np.ndarray, list]:
    """log P(observed bases | genotype) for the 10 genotypes, from per-site
    base counts and one error probability."""
    n = counts.shape[0]
    hom = np.log(1 - err)
    other = np.log(err / 3)
    het = np.log(0.5 * (1 - err) + 0.5 * err / 3)
    tot = counts.sum(1)
    out = np.empty((n, 10), dtype=np.float64)
    for k, (a, b) in enumerate(GENOTYPES):
        if a == b:
            out[:, k] = counts[:, a] * hom + (tot - counts[:, a]) * other
        else:
            m = counts[:, a] + counts[:, b]
            out[:, k] = m * het + (tot - m) * other
    return out, GENOTYPES


def loglik_from_qualities(site: np.ndarray, base: np.ndarray, qual: np.ndarray, n_sites: int,
                          min_error: float = 0.005, max_error: float = 0.25) -> tuple[np.ndarray, np.ndarray]:
    """The same likelihoods, but with each read's own error probability.

    A base quality is the sequencer's estimate that the base is wrong, so a
    Q20 base is forty times as likely to be an error as a Q36 one and should
    not count the same. Observations are grouped by quality, which is a small
    set of integers, and each group's contribution is added in one pass.

    Taken at face value a high quality says a base is almost certainly the
    base that was in the molecule - and for ancient DNA that is false, because
    the commonest way a base is not what the genome says is deamination, which
    the sequencer cannot see and scores as a perfect read. A caller that
    believes the qualities therefore reads two damaged reads in thirty as a
    heterozygote. ``min_error`` is the floor on how wrong a base can be for
    reasons outside the sequencer's knowledge - damage, misalignment, a little
    contamination - and the quality is used only where it is worse than that."""
    ll = np.zeros((n_sites, 10), dtype=np.float64)
    depth = np.bincount(site, minlength=n_sites).astype(np.int64)
    for q in np.unique(qual):
        err = (min(max_error, max(min_error, 10.0 ** (-float(q) / 10.0))) if q >= 0 else max_error)
        hom, other = np.log(1 - err), np.log(err / 3)
        het = np.log(0.5 * (1 - err) + 0.5 * err / 3)
        m = qual == q
        c = np.zeros((n_sites, 4), dtype=np.int64)
        np.add.at(c, (site[m], base[m]), 1)
        tot = c.sum(1)
        for k, (a, b) in enumerate(GENOTYPES):
            if a == b:
                ll[:, k] += c[:, a] * hom + (tot - c[:, a]) * other
            else:
                pair = c[:, a] + c[:, b]
                ll[:, k] += pair * het + (tot - pair) * other
    return ll, depth


def minimum_depth_for_a_homozygote(err: float, min_quality: float) -> int:
    """How deep a site must be before a homozygote can clear ``min_quality``.

    Each read of a clean homozygote adds the same fixed evidence against the
    heterozygote that contains its allele, so a quality floor is also a depth
    floor — and it is a stricter one for homozygotes than for heterozygotes.
    Leaving the two inconsistent drops homozygous sites from the callable set
    while keeping heterozygous ones, which inflates every rate computed from
    them."""
    per_read = 10.0 * math.log10((1 - err) / (0.5 * (1 - err) + err / 6))
    return int(math.ceil(min_quality / per_read)) if per_read > 0 else 1


def diploid(counts: np.ndarray = None, min_depth: int = 8, max_depth: int | None = None, err: float = 0.005,
            min_quality: float = 30.0, pileup=None, n_sites: int | None = None) -> dict:
    """Diploid genotypes, from per-site base counts or from a pileup's own base
    qualities (``pileup=``, which is the better measurement).

    Returns allele codes a1, a2 (-1 where no call) and the call's quality: the
    evidence in decibans for the called genotype over the next best one, which
    is a likelihood ratio between two genotypes and not a posterior
    probability of being wrong. The depth
    floor is raised, if needed, to the depth at which a homozygote can reach
    the quality floor, so that callability does not depend on which genotype
    the site has."""
    if pileup is not None:
        n = int(n_sites if n_sites is not None else pileup.n_sites)
        ll, depth = loglik_from_qualities(pileup.site, pileup.base, pileup.qual, n, min_error=err)
        counts = None
        quality_source = (f"each read's own base quality, with a floor of {err} on how wrong a base can be for "
                          "reasons the sequencer cannot see, deamination above all")
    else:
        counts = np.asarray(counts, dtype=np.int64)
        ll, _ = _genotype_loglik(counts, err)
        depth = counts.sum(1)
        quality_source = f"one error probability of {err} for every read"
    floor = max(int(min_depth), minimum_depth_for_a_homozygote(err, min_quality))
    best = np.argmax(ll, axis=1)
    top = ll[np.arange(len(ll)), best]
    second = np.where(np.isneginf(top), -np.inf, np.partition(ll, -2, axis=1)[:, -2])
    qual = (top - second) * 10.0 / np.log(10.0)
    ok = (depth >= floor) & (qual >= min_quality)
    if max_depth:
        ok &= depth <= max_depth
    pa = np.array([p[0] for p in GENOTYPES], dtype=np.int8)
    pb = np.array([p[1] for p in GENOTYPES], dtype=np.int8)
    a1 = np.where(ok, pa[best], -1).astype(np.int8)
    a2 = np.where(ok, pb[best], -1).astype(np.int8)
    return {"a1": a1, "a2": a2, "quality": qual, "depth": depth, "called": ok, "depth_floor": floor,
            "method": f"genotype likelihoods with a flat prior and {quality_source}; a site is called at depth "
                      f"{floor} or more{'' if not max_depth else f' and at most {max_depth}'}, with at least "
                      f"{min_quality} decibans in favour of the called genotype over the next best. The depth "
                      f"floor is {floor} rather than the {min_depth} asked for because below it a homozygote "
                      f"cannot reach {min_quality} decibans, and dropping homozygotes but not heterozygotes "
                      "would inflate every rate computed from the calls"
                      if floor > int(min_depth) else
                      f"genotype likelihoods with a flat prior and {quality_source}; a site is called at depth "
                      f"{floor} or more{'' if not max_depth else f' and at most {max_depth}'}, with at least "
                      f"{min_quality} decibans in favour of the called genotype over the next best"}


def panel_genotypes(p, ref: np.ndarray, alt: np.ndarray, min_qual: int = 20, trim: int = 0, **kw) -> dict:
    """Diploid alternative-allele dosage (0, 1, 2; -1 no call) at panel sites."""
    q = p.filtered(min_qual=min_qual, trim=trim)
    g = diploid(pileup=q, **kw)
    ref = np.asarray(ref, dtype=np.int8)
    alt = np.asarray(alt, dtype=np.int8)
    a1, a2 = g["a1"], g["a2"]
    known = g["called"] & np.isin(a1, [0, 1, 2, 3]) & ((a1 == ref) | (a1 == alt)) & ((a2 == ref) | (a2 == alt))
    dos = np.full(len(ref), -1, dtype=np.int8)
    dos[known] = ((a1[known] == alt[known]).astype(np.int8) + (a2[known] == alt[known]).astype(np.int8))
    depth = g["depth"]
    return {"dosage": dos, "quality": g["quality"], "depth": depth,
            "sites_called": int(known.sum()), "sites_third_allele": int((g["called"] & ~known).sum()),
            "sites_covered": int((depth > 0).sum()), "depth_floor": g["depth_floor"],
            "mean_depth_at_covered": float(depth[depth > 0].mean()) if (depth > 0).any() else 0.0,
            "observations_used": int(q.reads_used), "method": g["method"]}


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
    """Per-position base counts (A, C, G, T) over a region, 0-based half-open.

    Counting straight off the alignment, as this does, counts both reads of an
    unmerged pair where they overlap; ``heterozygosity_over_region`` reads the
    same region molecule by molecule instead."""
    ref = genome.name_of(chrom)
    mq = genome.min_mapq if min_mapq is None else int(min_mapq)

    def keep(r):
        return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail
                    or r.is_duplicate) and r.mapping_quality >= mq

    cov = genome.bam.count_coverage(ref, start, end, quality_threshold=min_qual, read_callback=keep)
    return np.array(cov, dtype=np.int32).T          # positions x 4


def heterozygosity_over_region(genome, chrom: str, start: int, end: int, callable_mask=None, min_qual: int = 25,
                               trim: int = 0, min_depth: int = 8, max_depth: int | None = None,
                               min_quality: float = 30.0, err: float = 0.005, chunk: int = 100_000,
                               on_step=None) -> dict:
    """Heterozygous sites per callable base over a region, read in chunks.

    Every position is piled up as the panel's positions are: one observation
    per molecule (so the overlapping mates of a pair count once), each with its
    own base quality, and with the option to ignore the ends of molecules where
    deamination concentrates."""
    from .sites import pileup
    totals = {"callable_sites": 0, "heterozygous_sites": 0, "transitions": 0, "transversions": 0,
              "observations": 0, "depth_floor": 0}
    pos0 = start
    while pos0 < end:
        hi = min(end, pos0 + chunk)
        if on_step:
            on_step({"task": "deeptime", "stage": f"counting molecules on {chrom}:{pos0 + 1:,}-{hi:,}",
                     "fraction": 0.1 + 0.8 * (pos0 - start) / max(1, end - start)})
        positions = np.arange(pos0 + 1, hi + 1, dtype=np.int64)        # 1-based
        p = pileup(genome, chrom, positions, max_depth=max_depth or 0).filtered(min_qual=min_qual, trim=trim)
        mask = None if callable_mask is None else np.asarray(callable_mask, dtype=bool)[pos0 - start:hi - start]
        r = heterozygosity(None, callable_mask=mask, min_depth=min_depth, max_depth=max_depth, err=err,
                           min_quality=min_quality, pileup=p, n_sites=len(positions))
        for k in ("callable_sites", "heterozygous_sites", "transitions", "transversions"):
            totals[k] += r[k]
        totals["observations"] += int(len(p.site))
        totals["depth_floor"] = r["depth_floor"]
        totals["method"] = r["method"]
        pos0 = hi
    n, nh, ntv = totals["callable_sites"], totals["heterozygous_sites"], totals["transversions"]
    totals.update({"heterozygosity": (nh / n) if n else None, "heterozygosity_ci": _poisson_rate_ci(nh, n),
                   "transversion_heterozygosity": (ntv / n) if n else None,
                   "transversion_heterozygosity_ci": _poisson_rate_ci(ntv, n),
                   "transition_share": (totals["transitions"] / nh) if nh else None})
    return totals


def heterozygosity(counts: np.ndarray, callable_mask: np.ndarray | None = None, min_depth: int = 8,
                   max_depth: int | None = None, err: float = 0.005, min_quality: float = 30.0,
                   pileup=None, n_sites: int | None = None) -> dict:
    """Heterozygous sites per callable base, from per-position base counts.

    No reference sequence is needed: a site is heterozygous when the reads
    carry two alleles. Transitions are reported separately because C->T and
    G->A deamination mimics them, so the transversion rate is the one that
    survives damage; for an undamaged or enzyme-treated library the total is
    the comparable number.

    The denominator counts only sites deep enough for a homozygote to be
    called, which is what makes the rate a rate: a floor that lets
    heterozygotes through but not homozygotes would count the same site in the
    numerator and not in the denominator."""
    g = diploid(counts, min_depth=min_depth, max_depth=max_depth, err=err, min_quality=min_quality,
                pileup=pileup, n_sites=n_sites)
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
            "heterozygosity_ci": _poisson_rate_ci(nh, n),
            "transversion_heterozygosity": (ntv / n) if n else None,
            "transversion_heterozygosity_ci": _poisson_rate_ci(ntv, n),
            "transition_share": (nts / nh) if nh else None, "depth_floor": g["depth_floor"],
            "method": g["method"]}


def _poisson_rate_ci(k: int, n: int, z: float = 1.96) -> list[float] | None:
    """A rate this small is a count of rare events, so its uncertainty is the
    count's: an interval from one event is not a measurement."""
    if not n:
        return None
    lo = max(0.0, k - z * math.sqrt(max(k, 0)))
    hi = k + z * math.sqrt(max(k, 1))
    return [lo / n, hi / n]


def bed_mask(path: str, chrom: str, start: int, end: int) -> np.ndarray:
    """Boolean mask over [start, end) from a BED file's intervals on ``chrom``."""
    import gzip
    from .bam import plain
    mask = np.zeros(end - start, dtype=bool)
    want = plain(chrom)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            if len(f) < 3 or plain(f[0]) != want:
                continue
            a, b = int(f[1]), int(f[2])
            if b <= start or a >= end:
                continue
            mask[max(0, a - start):min(end - start, b - start)] = True
    return mask
