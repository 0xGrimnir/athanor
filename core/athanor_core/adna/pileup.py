"""Aligned reads (BAM) to per-base observations for the consensus model."""
from __future__ import annotations

import numpy as np

from .consensus import CODE, Observations


PARTIAL = ":p"   # read-name suffix: the read's 3' end is not the molecule's end (see reads.py)


def is_partial(name: str) -> bool:
    return name.endswith(PARTIAL)


def observations_from_bam(bam_path: str, contig: str | None = None, min_mapq: int = 25,
                          min_baseq: int = 0, circular_length: int | None = None) -> Observations:
    """One observation per aligned (non-gap, non-clipped) base of each primary,
    mapped read with mapping quality >= min_mapq.

    Distances are measured on the sequenced molecule: for reads on the reverse
    strand the 5' end is the end of the alignment in reference coordinates.
    Reads are expected to be full fragments (adapter-trimmed, pairs collapsed),
    so both ends carry the deamination signal, except reads named with the
    ``:p`` suffix, whose 3' end is where sequencing stopped.

    circular_length: for a circular genome mapped to a reference elongated
    by its own start (see bamtools.py), reference positions are folded back
    modulo this length."""
    import pysam

    site, base, qual, d5, d3, rev, full = [], [], [], [], [], [], []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam.fetch(contig) if contig else bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < min_mapq:
                continue
            seq = r.query_sequence
            q = r.query_qualities
            qlen = r.query_length
            pairs = np.array(r.get_aligned_pairs(matches_only=True), dtype=np.int64)
            if len(pairs) == 0:
                continue
            qpos, rpos = pairs[:, 0], pairs[:, 1]
            codes = np.array([CODE.get(seq[i], -1) for i in qpos])
            quals = np.array([q[i] for i in qpos]) if q is not None else np.full(len(qpos), 30)
            keep = (codes >= 0) & (quals >= min_baseq)
            if r.is_reverse:
                dist5, dist3 = qlen - 1 - qpos, qpos
            else:
                dist5, dist3 = qpos, qlen - 1 - qpos
            if circular_length:
                rpos = rpos % circular_length
            full.append(np.full(keep.sum(), not is_partial(r.query_name)))
            site.append(rpos[keep]); base.append(codes[keep]); qual.append(quals[keep])
            d5.append(dist5[keep]); d3.append(dist3[keep]); rev.append(np.full(keep.sum(), r.is_reverse))
    cat = lambda xs, dt: np.concatenate(xs).astype(dt) if xs else np.zeros(0, dt)  # noqa: E731
    return Observations(site=cat(site, np.int64), base=cat(base, np.int64), qual=cat(qual, np.float64),
                        dist5=cat(d5, np.int64), dist3=cat(d3, np.int64), reverse=cat(rev, bool),
                        full=cat(full, bool))
