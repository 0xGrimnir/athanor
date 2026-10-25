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

import numpy as np

from .sites import CODE, COMPLEMENT, MAXD


def _empty():
    return {"hit": np.zeros(MAXD), "tot": np.zeros(MAXD), "hit_cpg": np.zeros(MAXD), "tot_cpg": np.zeros(MAXD)}


def _finish(acc: dict, background_from: int = 15) -> dict:
    out = {}
    for key, a in acc.items():
        d = {}
        for label, hit, tot in (("all", a["hit"] + a["hit_cpg"], a["tot"] + a["tot_cpg"]),
                                ("cpg", a["hit_cpg"], a["tot_cpg"]), ("not_cpg", a["hit"], a["tot"])):
            with np.errstate(invalid="ignore", divide="ignore"):
                freq = np.where(tot > 0, hit / np.maximum(tot, 1), np.nan)
            bg_t, bg_h = tot[background_from:].sum(), hit[background_from:].sum()
            bg = float(bg_h / bg_t) if bg_t else float("nan")
            d[label] = {"frequency": [None if np.isnan(x) else round(float(x), 5) for x in freq],
                        "observations": tot.astype(np.int64).tolist(),
                        "background": None if np.isnan(bg) else round(bg, 5),
                        "terminal": None if tot[0] == 0 else round(float(freq[0]), 5),
                        "terminal_excess": (None if (tot[0] == 0 or np.isnan(bg))
                                            else round(float(freq[0] - bg), 5))}
        out[key] = {**d["all"], "cpg": d["cpg"], "not_cpg": d["not_cpg"]}
    first = [out.get(k, {}).get("not_cpg", {}).get("terminal_excess") for k in ("ct5", "ga3")]
    out["terminal_deamination"] = None if any(x is None for x in first) else round(float(np.mean(first)), 5)
    first_all = [out.get(k, {}).get("terminal_excess") for k in ("ct5", "ga3")]
    out["terminal_deamination_including_cpg"] = (None if any(x is None for x in first_all)
                                                 else round(float(np.mean(first_all)), 5))
    return out


def from_reference(genome, fasta_path: str, chroms: list[str] | None = None, max_reads: int = 300_000,
                   min_qual: int = 20, background_from: int = 15, on_step=None) -> dict:
    """C->T and G->A profiles from every aligned base, against the reference."""
    import pysam
    fa = pysam.FastaFile(fasta_path)
    names = {n.replace("chr", "").upper(): n for n in fa.references}
    acc = {"ct5": _empty(), "ga3": _empty()}
    reads = used = 0
    chroms = chroms or [c for c in (str(i) for i in range(1, 23)) if genome.name_of(c)]
    chroms = [c for c in chroms if c.replace("chr", "").upper() in names]
    if not chroms:
        raise ValueError("none of the alignment's chromosomes are in the reference FASTA "
                         f"(it has {', '.join(fa.references[:5])}…)")
    for k, chrom in enumerate(chroms):
        ref_name = names[chrom.replace("chr", "").upper()]
        if on_step:
            on_step({"stage": f"deamination on chromosome {chrom}", "fraction": k / max(1, len(chroms))})
        for r in genome.bam.fetch(genome.name_of(chrom)):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate \
                    or r.mapping_quality < genome.min_mapq:
                continue
            reads += 1
            if reads > max_reads:
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
            quals = (np.fromiter((q[i] for i in qpos), dtype=np.int16, count=len(qpos)) if q is not None
                     else np.full(len(qpos), 30, dtype=np.int16))
            qlen = r.query_length
            if r.is_reverse:
                d5, d3 = qlen - 1 - qpos, qpos
                refm = np.where(rb >= 0, COMPLEMENT[np.maximum(rb, 0)], -1)
                obsm = np.where(ob >= 0, COMPLEMENT[np.maximum(ob, 0)], -1)
                nxt = np.where(pb >= 0, COMPLEMENT[np.maximum(pb, 0)], -1)     # 3' on the molecule
                prv = np.where(nb >= 0, COMPLEMENT[np.maximum(nb, 0)], -1)
            else:
                d5, d3 = qpos, qlen - 1 - qpos
                refm, obsm, nxt, prv = rb, ob, nb, pb
            ok = (refm >= 0) & (obsm >= 0) & (quals >= min_qual)
            used += 1
            for key, base, prod, dist, ctx in (("ct5", CODE["C"], CODE["T"], d5, nxt == CODE["G"]),
                                               ("ga3", CODE["G"], CODE["A"], d3, prv == CODE["C"])):
                m = ok & (refm == base) & (dist < MAXD)
                if not m.any():
                    continue
                dd, hit, cpg = dist[m], (obsm[m] == prod), ctx[m]
                a = acc[key]
                a["tot"] += np.bincount(dd[~cpg], minlength=MAXD)[:MAXD]
                a["hit"] += np.bincount(dd[~cpg], weights=hit[~cpg].astype(float), minlength=MAXD)[:MAXD]
                a["tot_cpg"] += np.bincount(dd[cpg], minlength=MAXD)[:MAXD]
                a["hit_cpg"] += np.bincount(dd[cpg], weights=hit[cpg].astype(float), minlength=MAXD)[:MAXD]
        if reads > max_reads:
            break
    fa.close()
    out = _finish(acc, background_from)
    out["reads_used"] = used
    out["chromosomes"] = chroms
    out["truncated"] = reads > max_reads
    out["method"] = ("every aligned base compared with the reference the reads were mapped to: C->T by distance "
                     "from the molecule's 5' end and G->A from its 3' end, on the molecule's own strand. The rate "
                     f"at positions {background_from} and beyond is the background of real differences and "
                     "sequencing errors and is subtracted. Cytosines followed by guanine are counted separately, "
                     "because methylated CpG deaminates several times faster.")
    return out


def from_panel(p, ref: np.ndarray, alt: np.ndarray, min_qual: int = 20, background_from: int = 15) -> dict:
    """The same profiles without a reference sequence, using the panel's own
    alleles: only sites whose other allele is not the deamination product, so
    that a real allele cannot be counted as damage. A SNP panel is rich in CpG
    sites, so the figure to compare with a published one is the non-CpG rate."""
    ref = np.asarray(ref, dtype=np.int8)
    alt = np.asarray(alt, dtype=np.int8)
    rb, ab = ref[p.site], alt[p.site]
    refm = np.where(p.reverse, COMPLEMENT[rb], rb)
    altm = np.where(p.reverse, COMPLEMENT[ab], ab)
    basem = np.where(p.reverse, COMPLEMENT[p.base], p.base)
    q = p.qual >= min_qual
    acc = {"ct5": _empty(), "ga3": _empty()}
    for key, base, prod, dist, ctx in (("ct5", CODE["C"], CODE["T"], p.d5, p.nxt == CODE["G"]),
                                       ("ga3", CODE["G"], CODE["A"], p.d3, p.prev == CODE["C"])):
        m = q & (refm == base) & (altm != prod) & (dist < MAXD) & ((basem == base) | (basem == prod))
        dd, hit, cpg = dist[m], (basem[m] == prod), ctx[m]
        a = acc[key]
        a["tot"] += np.bincount(dd[~cpg], minlength=MAXD)[:MAXD]
        a["hit"] += np.bincount(dd[~cpg], weights=hit[~cpg].astype(float), minlength=MAXD)[:MAXD]
        a["tot_cpg"] += np.bincount(dd[cpg], minlength=MAXD)[:MAXD]
        a["hit_cpg"] += np.bincount(dd[cpg], weights=hit[cpg].astype(float), minlength=MAXD)[:MAXD]
    return finish_panel(acc, background_from)


def add(a: dict, b: dict) -> dict:
    """Combine two accumulators (used to pool chromosomes)."""
    out = {}
    for key in a:
        out[key] = {k: a[key][k] + b[key][k] for k in a[key]}
    return out


def accumulate_panel(p, ref: np.ndarray, alt: np.ndarray, acc: dict | None = None, min_qual: int = 20) -> dict:
    """Panel-site counts for one chromosome, to be pooled before finishing."""
    ref = np.asarray(ref, dtype=np.int8)
    alt = np.asarray(alt, dtype=np.int8)
    rb, ab = ref[p.site], alt[p.site]
    refm = np.where(p.reverse, COMPLEMENT[rb], rb)
    altm = np.where(p.reverse, COMPLEMENT[ab], ab)
    basem = np.where(p.reverse, COMPLEMENT[p.base], p.base)
    q = p.qual >= min_qual
    acc = acc or {"ct5": _empty(), "ga3": _empty()}
    for key, base, prod, dist, ctx in (("ct5", CODE["C"], CODE["T"], p.d5, p.nxt == CODE["G"]),
                                       ("ga3", CODE["G"], CODE["A"], p.d3, p.prev == CODE["C"])):
        m = q & (refm == base) & (altm != prod) & (dist < MAXD) & ((basem == base) | (basem == prod))
        dd, hit, cpg = dist[m], (basem[m] == prod), ctx[m]
        a = acc[key]
        a["tot"] += np.bincount(dd[~cpg], minlength=MAXD)[:MAXD]
        a["hit"] += np.bincount(dd[~cpg], weights=hit[~cpg].astype(float), minlength=MAXD)[:MAXD]
        a["tot_cpg"] += np.bincount(dd[cpg], minlength=MAXD)[:MAXD]
        a["hit_cpg"] += np.bincount(dd[cpg], weights=hit[cpg].astype(float), minlength=MAXD)[:MAXD]
    return acc


PANEL_METHOD = (
        "measured at the reference panel's own positions, the only places a reference base is known without the "
        "reference sequence: C->T against the molecule's 5' end and G->A against its 3' end, at sites whose other "
        "allele is not the deamination product. SNP panels are rich in CpG sites, where deamination is several "
        "times faster, so the non-CpG rate is the one to compare with a published genome-wide figure. Give the "
        "reference FASTA for the full measurement.")


def finish_panel(acc: dict, background_from: int = 15) -> dict:
    out = _finish(acc, background_from)
    out["method"] = PANEL_METHOD
    return out
