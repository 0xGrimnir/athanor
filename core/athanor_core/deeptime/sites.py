"""Reads at chosen positions, and the deamination they carry.

A pileup at a list of positions: for every read overlapping one, the base it
shows, its base quality, which strand it came from, and how far the base is
from each end of the molecule. Everything downstream - genotypes, damage,
haplogroups - is counted from this.

Deamination. Cytosine loses an amino group with age and is read as thymine,
so ancient fragments carry C->T near their 5' ends and, as the complement,
G->A near their 3' ends (Briggs et al., PNAS 104, 14616, 2007). Measuring it
needs reference bases; here they come from the panel of known SNPs, using
only sites where the deamination product is not the site's other allele, so
that a real allele cannot be mistaken for damage. The rate far from the ends
is the background of real differences and sequencing errors, and the excess
over it at the ends is the damage.
"""
from __future__ import annotations

import numpy as np

BASES = "ACGT"
CODE = {b: i for i, b in enumerate(BASES)}
CODE.update({b.lower(): i for i, b in enumerate(BASES)})
MAXD = 25          # positions from each end kept for the damage profile
COMPLEMENT = np.array([3, 2, 1, 0], dtype=np.int8)      # A<->T, C<->G


class Pileup:
    """Observations at a set of positions on one chromosome.

    site       index of the position (into the position array given)
    base       0-3 (A, C, G, T)
    qual       base quality
    d5, d3     distance of the base from the molecule's 5' and 3' ends
    reverse    the read is on the reverse strand
    prev, nxt  the molecule's own neighbours of this base (-1 at the ends),
               which give the dinucleotide context: deamination is much faster
               at a cytosine followed by a guanine
    """

    __slots__ = ("site", "base", "qual", "d5", "d3", "reverse", "prev", "nxt", "n_sites", "reads", "reads_used")

    def __init__(self, site, base, qual, d5, d3, reverse, n_sites, reads, reads_used, prev=None, nxt=None):
        self.site, self.base, self.qual = site, base, qual
        self.d5, self.d3, self.reverse = d5, d3, reverse
        self.n_sites, self.reads, self.reads_used = n_sites, reads, reads_used
        z = np.full(len(site), -1, dtype=np.int8)
        self.prev = z if prev is None else prev
        self.nxt = z.copy() if nxt is None else nxt

    def counts(self) -> np.ndarray:
        """n_sites x 4 counts of each base."""
        out = np.zeros((self.n_sites, 4), dtype=np.int32)
        np.add.at(out, (self.site, self.base), 1)
        return out

    def filtered(self, min_qual: int = 20, trim: int = 0):
        keep = (self.qual >= min_qual) & (self.d5 >= trim) & (self.d3 >= trim)
        return Pileup(self.site[keep], self.base[keep], self.qual[keep], self.d5[keep], self.d3[keep],
                      self.reverse[keep], self.n_sites, self.reads, int(keep.sum()),
                      self.prev[keep], self.nxt[keep])


def pileup(genome, chrom: str, positions: np.ndarray, min_mapq: int | None = None, max_depth: int = 1000,
           drop_duplicates: bool = True) -> Pileup:
    """Reads over ``positions`` (1-based, sorted) on one chromosome."""
    ref = genome.name_of(chrom)
    positions = np.asarray(positions, dtype=np.int64)
    n = len(positions)
    if ref is None or n == 0:
        z = np.zeros(0, dtype=np.int32)
        return Pileup(z, z.astype(np.int8), z.astype(np.int16), z, z, np.zeros(0, bool), n, 0, 0)
    mq = genome.min_mapq if min_mapq is None else int(min_mapq)
    lo, hi = int(positions[0]), int(positions[-1])
    site_l, base_l, qual_l, d5_l, d3_l, rev_l, prev_l, next_l = [], [], [], [], [], [], [], []
    seen = np.zeros(n, dtype=np.int32)
    reads = used = 0
    for r in genome.bam.fetch(ref, max(0, lo - 1), hi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or (drop_duplicates and r.is_duplicate):
            continue
        reads += 1
        if r.mapping_quality < mq:
            continue
        seq = r.query_sequence
        if not seq:
            continue
        q = r.query_qualities
        pairs = r.get_aligned_pairs(matches_only=True)
        if not pairs:
            continue
        qpos = np.fromiter((p[0] for p in pairs), dtype=np.int64, count=len(pairs))
        rpos = np.fromiter((p[1] for p in pairs), dtype=np.int64, count=len(pairs)) + 1    # 1-based
        j = np.searchsorted(positions, rpos)
        on = (j < n) & (positions[np.minimum(j, n - 1)] == rpos)
        if not on.any():
            continue
        j, qpos = j[on], qpos[on]
        room = seen[j] < max_depth
        if not room.any():
            continue
        j, qpos = j[room], qpos[room]
        codes = np.fromiter((CODE.get(seq[i], -1) for i in qpos), dtype=np.int8, count=len(qpos))
        ok = codes >= 0
        if not ok.any():
            continue
        j, qpos, codes = j[ok], qpos[ok], codes[ok]
        quals = (np.fromiter((q[i] for i in qpos), dtype=np.int16, count=len(qpos)) if q is not None
                 else np.full(len(qpos), 30, dtype=np.int16))
        qlen = r.query_length
        if r.is_reverse:
            d5, d3 = qlen - 1 - qpos, qpos
        else:
            d5, d3 = qpos, qlen - 1 - qpos
        used += 1
        np.add.at(seen, j, 1)
        # the molecule's own neighbours: for a reverse read the stored sequence is
        # the reverse complement, so the molecule's 3' neighbour lies to the left
        s_arr = np.frombuffer(seq.encode(), dtype="S1")
        codes_all = np.full(len(seq), -1, dtype=np.int8)
        for b, c in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("a", 0), ("c", 1), ("g", 2), ("t", 3)):
            codes_all[s_arr == b.encode()] = c
        left = np.where(qpos > 0, codes_all[np.maximum(qpos - 1, 0)], -1)
        right = np.where(qpos < len(seq) - 1, codes_all[np.minimum(qpos + 1, len(seq) - 1)], -1)
        if r.is_reverse:
            prev_m = np.where(right >= 0, COMPLEMENT[np.maximum(right, 0)], -1)
            next_m = np.where(left >= 0, COMPLEMENT[np.maximum(left, 0)], -1)
        else:
            prev_m, next_m = left, right
        site_l.append(j.astype(np.int32)); base_l.append(codes); qual_l.append(quals)
        d5_l.append(d5.astype(np.int16)); d3_l.append(d3.astype(np.int16))
        rev_l.append(np.full(len(j), r.is_reverse, dtype=bool))
        prev_l.append(prev_m.astype(np.int8)); next_l.append(next_m.astype(np.int8))

    def cat(parts, dtype):
        return np.concatenate(parts).astype(dtype) if parts else np.zeros(0, dtype=dtype)

    return Pileup(cat(site_l, np.int32), cat(base_l, np.int8), cat(qual_l, np.int16), cat(d5_l, np.int16),
                  cat(d3_l, np.int16), cat(rev_l, bool), n, reads, used,
                  cat(prev_l, np.int8), cat(next_l, np.int8))
