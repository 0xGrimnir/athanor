"""JSON-RPC methods for Deep Time (ancient nuclear genomes).

  deeptime_inspect {path}        -> what the alignment is: reference build, read counts, sex, molecule lengths
  deeptime_analyse {path, ...}   -> authentication, sex, genotypes and the deep-ancestry report cards
  deeptime_heterozygosity {path, chrom, start, end, mask_path?}
                                 -> heterozygous sites per callable base in a region (high coverage)
  deeptime_export {path, what}   -> write the last report (JSON) or its genotypes (EIGENSTRAT-style text)

Everything runs on the machine the reads are on. Nothing is uploaded.
"""
from __future__ import annotations

import json
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
    kw = {k: p[k] for k in ("mode", "seed", "min_mapq", "min_qual", "trim", "transversions_only", "min_depth")
          if k in p}
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
        _step({"task": "deeptime", "stage": f"counting bases on {chrom}", "fraction": 0.2})
        counts = C.region_counts(g, chrom, start - 1, end, min_qual=int(p.get("min_qual", 25)))
        mask = None
        if p.get("mask_path"):
            mask = C.bed_mask(_path(p, "mask_path"), chrom, start - 1, end)
            if p.get("mask_is_excluded"):
                mask = ~mask
        _step({"task": "deeptime", "stage": "calling genotypes", "fraction": 0.7})
        kw = {k: p[k] for k in ("min_depth", "max_depth", "err", "min_quality") if k in p}
        res = C.heterozygosity(counts, callable_mask=mask, **kw)
        _step({"task": "deeptime", "stage": "done", "fraction": 1.0})
        return {"region": {"chromosome": chrom, "start": start, "end": end, "bases": end - start + 1},
                "build": g.build, "mask": ({"path": os.path.basename(_path(p, "mask_path")),
                                            "callable_bases": int(mask.sum()),
                                            "excluded_file": bool(p.get("mask_is_excluded"))} if mask is not None
                                           else None),
                "depth_mean": float(counts.sum(1).mean()), **res,
                "limits": ["Heterozygosity here counts sites where the reads carry two alleles, so it needs enough "
                           "depth to see both; published values use callers built for ancient DNA (snpAD) and the "
                           "study's own filters, and a region is not the whole genome.",
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
