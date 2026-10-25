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

import zlib

import numpy as np

BASES = "ACGT"
CODE = {b: i for i, b in enumerate(BASES)}
CODE.update({b.lower(): i for i, b in enumerate(BASES)})
MAXD = 40          # positions from each end kept for the damage profile
NO_QUALITY = -1    # a read that carries no base qualities: no quality filter can apply to it
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

    __slots__ = ("site", "base", "qual", "d5", "d3", "reverse", "prev", "nxt", "n_sites", "reads", "reads_used",
                 "reads_without_base_qualities")

    def __init__(self, site, base, qual, d5, d3, reverse, n_sites, reads, reads_used, prev=None, nxt=None):
        self.site, self.base, self.qual = site, base, qual
        self.d5, self.d3, self.reverse = d5, d3, reverse
        self.n_sites, self.reads, self.reads_used = n_sites, reads, reads_used
        self.reads_without_base_qualities = 0
        z = np.full(len(site), -1, dtype=np.int8)
        self.prev = z if prev is None else prev
        self.nxt = z.copy() if nxt is None else nxt

    def counts(self) -> np.ndarray:
        """n_sites x 4 counts of each base."""
        out = np.zeros((self.n_sites, 4), dtype=np.int32)
        np.add.at(out, (self.site, self.base), 1)
        return out

    def filtered(self, min_qual: int = 20, trim: int = 0):
        """Observations that pass a base-quality floor and are far enough from
        the molecule's ends. A read with no base qualities at all passes only
        when the floor is zero, so its bases are never taken for a quality they
        do not have."""
        keep = (self.qual >= min_qual) if min_qual > 0 else np.ones(len(self.qual), dtype=bool)
        keep = keep & (self.d5 >= trim) & (self.d3 >= trim)
        out = Pileup(self.site[keep], self.base[keep], self.qual[keep], self.d5[keep], self.d3[keep],
                     self.reverse[keep], self.n_sites, self.reads, int(keep.sum()),
                     self.prev[keep], self.nxt[keep])
        out.reads_without_base_qualities = self.reads_without_base_qualities
        return out


def molecule_distances(r, qpos: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """Distance of each aligned base from the molecule's 5' and 3' ends.

    The molecule is the whole sequenced fragment, which is longer than the
    stored sequence when the aligner hard-clipped it: those bases are gone from
    the record but were there in the molecule, so they still count towards the
    distances. Soft-clipped bases are in the stored sequence already."""
    cig = r.cigartuples or []
    lead_hard = cig[0][1] if cig and cig[0][0] == 5 else 0
    tail_hard = cig[-1][1] if len(cig) > 1 and cig[-1][0] == 5 else 0
    qlen = (r.query_length or 0) + lead_hard + tail_hard
    in_mol = qpos + lead_hard
    if r.is_reverse:                     # the stored sequence is the reverse complement
        return qlen - 1 - in_mol, in_mol, qlen
    return in_mol, qlen - 1 - in_mol, qlen


def _subsample(site: np.ndarray, n_sites: int, max_depth: int, rng) -> np.ndarray:
    """Keep at most ``max_depth`` observations per site, chosen uniformly.

    Taking the first ones in coordinate order would keep the reads that start
    earliest, which are exactly the ones whose bases sit near their 3' ends."""
    if max_depth is None or len(site) == 0:
        return np.ones(len(site), dtype=bool)
    depth = np.bincount(site, minlength=n_sites)
    if depth.max() <= max_depth:
        return np.ones(len(site), dtype=bool)
    order = np.lexsort((rng.random(len(site)), site))      # random order within each site
    rank = np.empty(len(site), dtype=np.int64)
    sorted_sites = site[order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_sites)) + 1]
    rank[order] = np.arange(len(site)) - np.repeat(starts, np.diff(np.r_[starts, len(site)]))
    return rank < max_depth


def _dedupe_mates(site, qual, molecule) -> np.ndarray:
    """One observation per molecule per site.

    The two reads of an unmerged pair are one molecule, so where they overlap
    they report the same base twice; counting both doubles the apparent depth
    and manufactures confidence. The higher-quality base of the pair is kept."""
    if len(site) == 0:
        return np.ones(0, dtype=bool)
    key = site.astype(np.int64) * (1 << 40) + molecule
    order = np.lexsort((-qual.astype(np.int64), key))
    first = np.r_[True, key[order][1:] != key[order][:-1]]
    keep = np.zeros(len(site), dtype=bool)
    keep[order[first]] = True
    return keep


def pileup(genome, chrom: str, positions: np.ndarray, min_mapq: int | None = None, max_depth: int = 1000,
           drop_duplicates: bool = True, seed: int = 0, merge_mates: bool = True) -> Pileup:
    """Reads over ``positions`` (1-based, sorted) on one chromosome."""
    ref = genome.name_of(chrom)
    positions = np.asarray(positions, dtype=np.int64)
    n = len(positions)
    if ref is None or n == 0:
        z = np.zeros(0, dtype=np.int32)
        return Pileup(z, z.astype(np.int8), z.astype(np.int16), z, z, np.zeros(0, bool), n, 0, 0)
    mq = genome.min_mapq if min_mapq is None else int(min_mapq)
    lo, hi = int(positions[0]), int(positions[-1])
    site_l, base_l, qual_l, d5_l, d3_l, rev_l, prev_l, next_l, mol_l = [], [], [], [], [], [], [], [], []
    reads = used = no_qual = 0
    # Everything collected is sampled down to max_depth afterwards, so an
    # ordinary site is sampled from all of its reads. This much higher limit
    # exists only so that a pathological pileup - a collapsed repeat at
    # thousands of times the genome's depth - cannot fill memory; a site that
    # reaches it keeps the reads that arrive first, which is a bias worth
    # having over running out of memory.
    collect_cap = (max_depth * 10) if max_depth else 0
    seen = np.zeros(n, dtype=np.int32)
    for r in genome.bam.fetch(ref, max(0, lo - 1), hi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail \
                or (drop_duplicates and r.is_duplicate):
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
        if collect_cap:
            room = seen[j] < collect_cap
            if not room.any():
                continue
            j, qpos = j[room], qpos[room]
            np.add.at(seen, j, 1)
        codes = np.fromiter((CODE.get(seq[i], -1) for i in qpos), dtype=np.int8, count=len(qpos))
        ok = codes >= 0
        if not ok.any():
            continue
        j, qpos, codes = j[ok], qpos[ok], codes[ok]
        if q is None:
            no_qual += 1
            quals = np.full(len(qpos), NO_QUALITY, dtype=np.int16)
        else:
            quals = np.fromiter((q[i] for i in qpos), dtype=np.int16, count=len(qpos))
        d5, d3, _ = molecule_distances(r, qpos)
        used += 1
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
        # the molecule a read belongs to: the pair share a name, so overlapping
        # mates of one fragment can be reduced to one observation per site
        mol_l.append(np.full(len(j), zlib.crc32(r.query_name.encode()), dtype=np.int64))

    def cat(parts, dtype):
        return np.concatenate(parts).astype(dtype) if parts else np.zeros(0, dtype=dtype)

    site, base = cat(site_l, np.int32), cat(base_l, np.int8)
    qual, d5, d3 = cat(qual_l, np.int16), cat(d5_l, np.int16), cat(d3_l, np.int16)
    rev, prev, nxt = cat(rev_l, bool), cat(prev_l, np.int8), cat(next_l, np.int8)
    mol = cat(mol_l, np.int64)
    if merge_mates and len(site):
        keep = _dedupe_mates(site, qual, mol)
        site, base, qual, d5, d3, rev, prev, nxt = (x[keep] for x in (site, base, qual, d5, d3, rev, prev, nxt))
    if max_depth and len(site):
        keep = _subsample(site, n, max_depth, np.random.default_rng(seed))
        site, base, qual, d5, d3, rev, prev, nxt = (x[keep] for x in (site, base, qual, d5, d3, rev, prev, nxt))
    p = Pileup(site, base, qual, d5, d3, rev, n, reads, used, prev, nxt)
    p.reads_without_base_qualities = no_qual
    return p
