"""Deamination measured against the reference the reads were mapped to.

With the reference sequence at hand, every aligned base of every read can be
compared with it, as mapDamage2 does (Jonsson et al., Bioinformatics 29, 1682,
2013): C->T towards the molecule's 5' end and G->A towards its 3' end, counted
by distance from the end, with the rate in the middle of the molecules as the
background of real differences and sequencing errors.

Cytosine followed by guanine is reported separately. Methylated CpG cytosines
deaminate faster, so their rate runs several times higher than the rest, and a
figure quoted without saying which it is cannot be compared with a published
one. Panels of chosen SNPs are rich in CpG sites, which is why the rate
measured at a SNP panel's positions is not the rate a whole genome gives.
"""
from __future__ import annotations

import math

import numpy as np

from .bam import plain
from .sites import CODE, COMPLEMENT, MAXD, molecule_distances


def _empty():
    return {"hit": np.zeros(MAXD), "tot": np.zeros(MAXD), "hit_cpg": np.zeros(MAXD), "tot_cpg": np.zeros(MAXD),
            "bg_hit": 0.0, "bg_tot": 0.0, "bg_hit_cpg": 0.0, "bg_tot_cpg": 0.0}


MIN_OBSERVATIONS = 500      # below this a terminal rate is reported as an interval, not a number


def _wilson(k: float, n: float, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return (0.0, 1.0)
    p, d = k / n, 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(max(p * (1 - p) / n + z * z / (4 * n * n), 0.0)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def _difference_ci(k1, n1, k2, n2, z: float = 1.96) -> tuple[float, float] | None:
    """Interval for p1 - p2 (the terminal rate above its background)."""
    if n1 <= 0 or n2 <= 0:
        return None
    p1, p2 = k1 / n1, k2 / n2
    se = math.sqrt(max(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2, 0.0))
    return (round((p1 - p2) - z * se, 5), round((p1 - p2) + z * se, 5))


def _tally(a: dict, dist: np.ndarray, hit: np.ndarray, cpg: np.ndarray, inner: np.ndarray) -> None:
    """Add one read's observations: the profile by distance from the end, and
    the background from bases far from both ends."""
    near = dist < MAXD
    for sel, tot_k, hit_k in ((~cpg & near, "tot", "hit"), (cpg & near, "tot_cpg", "hit_cpg")):
        if sel.any():
            a[tot_k] += np.bincount(dist[sel], minlength=MAXD)[:MAXD]
            a[hit_k] += np.bincount(dist[sel], weights=hit[sel].astype(float), minlength=MAXD)[:MAXD]
    for sel, tot_k, hit_k in ((inner & ~cpg, "bg_tot", "bg_hit"), (inner & cpg, "bg_tot_cpg", "bg_hit_cpg")):
        a[tot_k] += float(sel.sum())
        a[hit_k] += float(hit[sel].sum())


def _finish(acc: dict, background_from: int = 25) -> dict:
    """Rates by distance from the molecule's end, with the background taken far
    from *both* ends: on a short molecule a base 25 from one end can still be
    near the other, and a single-stranded library deaminates at both."""
    out = {}
    for key, a in acc.items():
        d = {}
        for label, hit, tot, bh, bt in (("all", a["hit"] + a["hit_cpg"], a["tot"] + a["tot_cpg"],
                                         a["bg_hit"] + a["bg_hit_cpg"], a["bg_tot"] + a["bg_tot_cpg"]),
                                        ("cpg", a["hit_cpg"], a["tot_cpg"], a["bg_hit_cpg"], a["bg_tot_cpg"]),
                                        ("not_cpg", a["hit"], a["tot"], a["bg_hit"], a["bg_tot"])):
            with np.errstate(invalid="ignore", divide="ignore"):
                freq = np.where(tot > 0, hit / np.maximum(tot, 1), np.nan)
            bg_t, bg_h = float(bt), float(bh)
            bg = (bg_h / bg_t) if bg_t else float("nan")
            enough = tot[0] >= MIN_OBSERVATIONS and bg_t >= MIN_OBSERVATIONS
            d[label] = {"frequency": [None if np.isnan(x) else round(float(x), 5) for x in freq],
                        "observations": tot.astype(np.int64).tolist(),
                        "background": None if np.isnan(bg) else round(bg, 5),
                        "background_observations": int(bg_t),
                        "terminal": None if tot[0] == 0 else round(float(freq[0]), 5),
                        "terminal_ci": [round(x, 5) for x in _wilson(hit[0], tot[0])] if tot[0] else None,
                        "terminal_excess": (None if (tot[0] == 0 or np.isnan(bg))
                                            else round(float(freq[0] - bg), 5)),
                        "terminal_excess_ci": _difference_ci(hit[0], tot[0], bg_h, bg_t),
                        "enough_observations": bool(enough)}
        out[key] = {**d["all"], "cpg": d["cpg"], "not_cpg": d["not_cpg"]}
    first = [out.get(k, {}).get("not_cpg", {}).get("terminal_excess") for k in ("ct5", "ga3")]
    out["terminal_deamination"] = None if any(x is None for x in first) else round(float(np.mean(first)), 5)
    first_all = [out.get(k, {}).get("terminal_excess") for k in ("ct5", "ga3")]
    out["terminal_deamination_including_cpg"] = (None if any(x is None for x in first_all)
                                                 else round(float(np.mean(first_all)), 5))
    lows = [out.get(k, {}).get("not_cpg", {}).get("terminal_excess_ci") for k in ("ct5", "ga3")]
    out["terminal_deamination_ci"] = (None if any(x is None for x in lows) else
                                      [round(float(np.mean([x[0] for x in lows])), 5),
                                       round(float(np.mean([x[1] for x in lows])), 5)])
    out["enough_observations"] = all(out.get(k, {}).get("not_cpg", {}).get("enough_observations")
                                     for k in ("ct5", "ga3"))
    out["minimum_observations"] = MIN_OBSERVATIONS
    out["background_from"] = background_from
    return out


def from_reference(genome, fasta_path: str, chroms: list[str] | None = None, max_reads: int = 300_000,
                   min_qual: int = 20, background_from: int = 25, on_step=None) -> dict:
    """C->T and G->A profiles from every aligned base, against the reference.

    Each chromosome gets its own share of the read budget, so the sample is not
    the start of whichever one comes first, and the result records the
    chromosomes that were actually read."""
    import pysam
    fa = pysam.FastaFile(fasta_path)
    names = {n.replace("chr", "").upper(): n for n in fa.references}
    acc = {"ct5": _empty(), "ga3": _empty()}
    used = 0
    chroms = chroms or [c for c in (str(i) for i in range(1, 23)) if genome.name_of(c)]
    chroms = [c for c in chroms if c.replace("chr", "").upper() in names and genome.name_of(c)]
    if not chroms:
        raise ValueError("none of the alignment's chromosomes are in the reference FASTA "
                         f"(it has {', '.join(fa.references[:5])}…)")
    mismatched = [(c, genome.lengths[plain(c)], fa.get_reference_length(names[c.replace("chr", "").upper()]))
                  for c in chroms
                  if genome.lengths.get(plain(c)) != fa.get_reference_length(names[c.replace("chr", "").upper()])]
    if mismatched:
        c, a, b = mismatched[0]
        raise ValueError(f"the reference FASTA is not the one the reads were mapped to: chromosome {c} is {b:,} "
                         f"bases in it and {a:,} in the alignment")
    budget = max(1, max_reads // len(chroms))
    per_chrom, truncated = {}, False
    for k, chrom in enumerate(chroms):
        ref_name = names[chrom.replace("chr", "").upper()]
        if on_step:
            on_step({"stage": f"deamination on chromosome {chrom}", "fraction": k / max(1, len(chroms))})
        reads = 0
        for r in genome.bam.fetch(genome.name_of(chrom)):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_duplicate \
                    or r.mapping_quality < genome.min_mapq:
                continue
            reads += 1
            if reads > budget:
                truncated = True
                break
            seq = r.query_sequence
            if not seq:
                continue
            q = r.query_qualities
            pairs = r.get_aligned_pairs(matches_only=True)
            if not pairs:
                continue
            qpos = np.fromiter((p[0] for p in pairs), dtype=np.int64, count=len(pairs))
            rpos = np.fromiter((p[1] for p in pairs), dtype=np.int64, count=len(pairs))
            refseq = fa.fetch(ref_name, int(rpos[0]), int(rpos[-1]) + 2).upper()
            off = int(rpos[0])
            rb = np.fromiter((CODE.get(refseq[p - off], -1) if p - off < len(refseq) else -1 for p in rpos),
                             dtype=np.int8, count=len(rpos))
            nb = np.fromiter((CODE.get(refseq[p - off + 1], -1) if p - off + 1 < len(refseq) else -1 for p in rpos),
                             dtype=np.int8, count=len(rpos))
            pb = np.fromiter((CODE.get(refseq[p - off - 1], -1) if 0 <= p - off - 1 < len(refseq) else -1
                              for p in rpos), dtype=np.int8, count=len(rpos))
            ob = np.fromiter((CODE.get(seq[i], -1) for i in qpos), dtype=np.int8, count=len(qpos))
            if q is None:                      # no base qualities: the floor cannot be applied
                continue
            quals = np.fromiter((q[i] for i in qpos), dtype=np.int16, count=len(qpos))
            d5, d3, _ = molecule_distances(r, qpos)
            if r.is_reverse:
                refm = np.where(rb >= 0, COMPLEMENT[np.maximum(rb, 0)], -1)
                obsm = np.where(ob >= 0, COMPLEMENT[np.maximum(ob, 0)], -1)
                nxt = np.where(pb >= 0, COMPLEMENT[np.maximum(pb, 0)], -1)     # 3' on the molecule
                prv = np.where(nb >= 0, COMPLEMENT[np.maximum(nb, 0)], -1)
            else:
                refm, obsm, nxt, prv = rb, ob, nb, pb
            ok = (refm >= 0) & (obsm >= 0) & (quals >= min_qual)
            used += 1
            inner = (d5 >= background_from) & (d3 >= background_from)
            for key, base, prod, dist, ctx in (("ct5", CODE["C"], CODE["T"], d5, nxt == CODE["G"]),
                                               ("ga3", CODE["G"], CODE["A"], d3, prv == CODE["C"])):
                m = ok & (refm == base)
                if not m.any():
                    continue
                _tally(acc[key], dist[m], obsm[m] == prod, ctx[m], inner[m])
            per_chrom[chrom] = per_chrom.get(chrom, 0) + 1
    fa.close()
    out = _finish(acc, background_from)
    out["reads_used"] = used
    out["chromosomes"] = sorted(per_chrom, key=lambda c: chroms.index(c))
    out["reads_by_chromosome"] = per_chrom
    out["truncated"] = truncated
    out["method"] = ("every aligned base compared with the reference the reads were mapped to: C->T by distance "
                     "from the molecule's 5' end and G->A from its 3' end, on the molecule's own strand. The "
                     f"background of real differences and sequencing errors is the rate at bases {background_from} "
                     "or more from both ends of their molecule, and is subtracted. Cytosines followed by guanine "
                     "are counted separately, because methylated CpG deaminates several times faster. Each "
                     "chromosome contributes an equal share of the reads read.")
    return out


def accumulate_panel(p, ref: np.ndarray, alt: np.ndarray, acc: dict | None = None, min_qual: int = 20,
                     background_from: int = 25) -> dict:
    """Panel-site counts for one chromosome, to be pooled before finishing."""
    ref = np.asarray(ref, dtype=np.int8)
    alt = np.asarray(alt, dtype=np.int8)
    rb, ab = ref[p.site], alt[p.site]
    refm = np.where(p.reverse, COMPLEMENT[rb], rb)
    altm = np.where(p.reverse, COMPLEMENT[ab], ab)
    basem = np.where(p.reverse, COMPLEMENT[p.base], p.base)
    q = p.qual >= min_qual
    inner_all = (p.d5 >= background_from) & (p.d3 >= background_from)
    acc = acc or {"ct5": _empty(), "ga3": _empty()}
    for key, base, prod, dist, ctx in (("ct5", CODE["C"], CODE["T"], p.d5, p.nxt == CODE["G"]),
                                       ("ga3", CODE["G"], CODE["A"], p.d3, p.prev == CODE["C"])):
        m = q & (refm == base) & (altm != prod) & ((basem == base) | (basem == prod))
        _tally(acc[key], dist[m], basem[m] == prod, ctx[m], inner_all[m])
    return acc


def from_panel(p, ref: np.ndarray, alt: np.ndarray, min_qual: int = 20, background_from: int = 25) -> dict:
    """The same profiles without a reference sequence, using the panel's own
    alleles: only sites whose other allele is not the deamination product, so
    that a real allele cannot be counted as damage. The sites that survive that
    filter are not a random sample of the genome - it removes every C/T and G/A
    SNP, which are mostly CpG - so the rate it gives is not a genome-wide rate.
    The non-CpG figure is the closer of the two to a published one, but the
    reference sequence is what makes the measurement comparable."""
    return finish_panel(accumulate_panel(p, ref, alt, min_qual=min_qual, background_from=background_from),
                        background_from)


PANEL_METHOD = (
        "measured at the reference panel's own positions, the only places a reference base is known without the "
        "reference sequence: C->T against the molecule's 5' end and G->A against its 3' end, at sites whose other "
        "allele is not the deamination product, with the background taken from bases far from both ends of their "
        "molecule. Those sites are not a random sample of the genome: the filter removes the C/T and G/A SNPs, "
        "which are mostly CpG, while SNP panels are otherwise rich in CpG, where deamination is several times "
        "faster. The dinucleotide context here comes from the read itself rather than from the reference, so it "
        "too can carry damage. Give the reference FASTA for a measurement comparable with published rates.")


def finish_panel(acc: dict, background_from: int = 25) -> dict:
    out = _finish(acc, background_from)
    out["method"] = PANEL_METHOD
    return out
