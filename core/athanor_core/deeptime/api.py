"""JSON-RPC methods for Deep Time (ancient nuclear genomes).

  deeptime_inspect {path}        -> what the alignment is: reference build, read counts, sex, molecule lengths
  deeptime_analyse {path, reference_path?, ...}
                                -> authentication, sex, genotypes and the deep-ancestry report cards
  deeptime_heterozygosity {path, chrom, start, end, mask_path?}
                                 -> heterozygous sites per callable base in a region (high coverage)
  deeptime_export {path, what}   -> write the last report (JSON) or its genotypes (EIGENSTRAT-style text)

Everything runs on the machine the reads are on. Nothing is uploaded.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np

from . import calls as C
from .analyse import analyse
from .bam import Genome, length_summary

PROGRESS = None
_LAST: dict = {}


def _step(d):
    if PROGRESS:
        PROGRESS(d)


def _path(p: dict, key: str = "path") -> str:
    v = (p or {}).get(key)
    if not v:
        raise ValueError(f"no {key} given")
    from ..server import to_linux_path
    return to_linux_path(v)


def m_inspect(p):
    """Read the header and index only: fast, whatever the file's size."""
    with Genome(_path(p), min_mapq=int((p or {}).get("min_mapq", 25))) as g:
        counts = g.counts()
        lp = g.lengths_and_damage_pass(max_reads=int((p or {}).get("max_reads", 200_000)))
        return {"file": os.path.basename(g.path), "bytes": os.path.getsize(g.path), "indexed": g.indexed,
                "build": g.build, "build_detail": g.build_detail, "mitochondrion_is_rCRS": g.mt_usable,
                "reads_by_chromosome": counts, "reads_mapped": int(sum(counts.values())),
                "read_groups": g.read_groups(), "programs": g.programs(), "sex": g.sex(),
                "lengths": length_summary(lp["length_histogram"]), "length_pass": lp}


def m_analyse(p):
    p = p or {}
    kw = {k: p[k] for k in ("mode", "seed", "min_mapq", "min_qual", "trim", "transversions_only", "min_depth",
                            "damage_chromosomes", "damage_max_reads") if k in p}
    if p.get("reference_path"):
        kw["reference_fasta"] = _path(p, "reference_path")
    if kw.get("mode") not in (None, "pseudo_haploid", "diploid"):
        raise ValueError("mode must be pseudo_haploid or diploid")
    res = analyse(_path(p), on_step=_step, **kw)
    _LAST.clear()
    _LAST.update(res)
    return res


def m_heterozygosity(p):
    p = p or {}
    chrom = str(p.get("chrom") or "")
    start, end = int(p.get("start", 0)), int(p.get("end", 0))
    if not chrom or end <= start:
        raise ValueError("give chrom, start and end (1-based start, inclusive end)")
    if end - start > 50_000_000:
        raise ValueError("the region is longer than 50 Mb; run it in parts")
    with Genome(_path(p), min_mapq=int(p.get("min_mapq", 25))) as g:
        if g.name_of(chrom) is None:
            raise ValueError(f"the alignment has no chromosome {chrom}")
        mask = None
        if p.get("mask_path"):
            mask = C.bed_mask(_path(p, "mask_path"), chrom, start - 1, end)
            if p.get("mask_is_excluded"):
                mask = ~mask
        kw = {k: p[k] for k in ("min_depth", "max_depth", "err", "min_quality", "trim") if k in p}
        res = C.heterozygosity_over_region(g, chrom, start - 1, end, callable_mask=mask,
                                           min_qual=int(p.get("min_qual", 25)), on_step=_step, **kw)
        _step({"task": "deeptime", "stage": "done", "fraction": 1.0})
        share = res.get("transition_share")
        nh = res.get("heterozygous_sites") or 0
        share_ci = None if not nh else [max(0.0, share - 1.96 * math.sqrt(share * (1 - share) / nh)),
                                        min(1.0, share + 1.96 * math.sqrt(share * (1 - share) / nh))]
        inflated = bool(share_ci and share_ci[0] > 0.72)   # above what human variation alone gives
        return {"region": {"chromosome": chrom, "start": start, "end": end, "bases": end - start + 1},
                "transitions_inflated_by_damage": inflated,
                "transition_share_ci": share_ci,
                "note": (f"Transitions are {100 * share:.0f}% of the heterozygous calls "
                         f"(95% {100 * share_ci[0]:.0f}-{100 * share_ci[1]:.0f}%), above the two thirds or so that "
                         "human variation alone gives, so deamination is being read as heterozygosity: take the "
                         "transversion rate as the sample's, or ignore the ends of the molecules."
                         if inflated else
                         ("Transitions are a normal share of the heterozygous calls, so the total is usable."
                          if nh else "No heterozygous sites were called.")),
                "build": g.build, "mask": ({"path": os.path.basename(_path(p, "mask_path")),
                                            "callable_bases": int(mask.sum()),
                                            "excluded_file": bool(p.get("mask_is_excluded"))} if mask is not None
                                           else None),
                "depth_mean": (res["observations"] / (end - start + 1)) if end > start else None, **res,
                "limits": [f"A site counts towards the total only when it has {res['depth_floor']} reads or more, "
                           "the depth at which a homozygote can be called as confidently as a heterozygote; "
                           "published values use callers built for ancient DNA (snpAD) and the study's own "
                           "filters, and a region is not the whole genome.",
                           "C->T and G->A deamination mimics heterozygous transitions, so the transversion rate is "
                           "the one that survives damage."]}


def m_export(p):
    p = p or {}
    what = p.get("what", "report")
    if not _LAST:
        raise ValueError("no Deep Time report to export yet")
    out = Path(_path(p))
    out.parent.mkdir(parents=True, exist_ok=True)
    if what == "report":
        out.write_text(json.dumps(_LAST, indent=1, default=_plain))
    else:
        raise ValueError("what must be report")
    return {"written": str(out), "bytes": out.stat().st_size}


def _plain(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


METHODS = {"deeptime_inspect": m_inspect, "deeptime_analyse": m_analyse,
           "deeptime_heterozygosity": m_heterozygosity, "deeptime_export": m_export}
