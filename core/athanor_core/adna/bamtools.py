"""Mapping clean-up for real data: circular genomes, duplicates, indels.

Circular genomes. Mitochondrial genomes are circular, but BWA maps to linear
sequences, so molecules spanning the arbitrary start/end of the reference
would be lost. As in CircularMapper (Peltzer et al., Genome Biol. 17, 60,
2016), reads are mapped to the reference elongated by a copy of its first
ELONGATION bases, and positions are folded back modulo the genome length.
A read lying entirely within the copied stretch then has two equally good
placements and BWA gives it mapping quality 0; when its only other hit is
the same placement one genome length away (and BWA reports no suboptimal
hits), it is unique on the circle and receives the mapping quality BWA aln
gives a unique hit without suboptimal hits (37).

Duplicates. PCR copies of one molecule share both ends. Whole molecules are
duplicates when strand, start and end agree; reads cut off by the read
length when strand and 5' end agree (as in samtools rmdup -s). Duplicates
are only sought within one input, because separate libraries cannot
share PCR copies. The copy with the highest summed base quality is kept.

Indels. The consensus model calls bases at reference positions. Insertions
and deletions relative to the reference are called by majority: an event
supported by at least ``min_reads`` reads and by more than half of the reads
spanning the site is applied to the reference for the next mapping round.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np

from .pileup import is_partial

ELONGATION = 300   # longer than any molecule we expect to map


def write_fasta(path: Path, contigs: dict[str, str]) -> Path:
    with open(path, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n")
    return path


def elongate(seq: str, circular: bool) -> str:
    return seq + seq[:min(ELONGATION, len(seq))] if circular else seq


def _xa_hits(tag: str) -> list[tuple[str, int, bool]]:
    """BWA's XA tag: 'chr,+pos,CIGAR,NM;' entries (1-based positions, strand sign)."""
    out = []
    for e in tag.strip(";").split(";"):
        p = e.split(",")
        if len(p) >= 2 and p[1]:
            out.append((p[0], abs(int(p[1])) - 1, p[1].startswith("-")))
    return out


def clean_bam(in_bam: Path, out_bam: Path, min_mapq: int, circular: dict[str, int] | None = None,
              strict_mapq: bool = False) -> dict:
    """Fold circular placements, rescue reads in the elongated stretch,
    remove duplicates and reads below the mapping-quality threshold.
    ``circular`` maps contig name -> true genome length for elongated contigs.
    ``strict_mapq`` keeps reads with MAPQ > min_mapq instead of >=."""
    import pysam
    circular = circular or {}
    n = Counter()
    best: dict = {}
    with pysam.AlignmentFile(str(in_bam), "rb") as b:
        header = b.header.to_dict()
        for r in b.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            n["reads"] += 1
            if r.is_unmapped:
                n["unmapped"] += 1
                continue
            contig = r.reference_name
            L = circular.get(contig)
            if L and r.mapping_quality == 0 and r.has_tag("X0") and r.get_tag("X0") == 2 \
                    and r.has_tag("X1") and r.get_tag("X1") == 0 and r.has_tag("XA"):
                alts = _xa_hits(r.get_tag("XA"))
                if len(alts) == 1 and alts[0][0] == contig and alts[0][2] == r.is_reverse \
                        and abs(alts[0][1] - r.reference_start) == L:
                    r.mapping_quality = 37
                    n["rescued_circular"] += 1
            if L and r.reference_start >= L:
                r.reference_start -= L          # lies wholly in the copied stretch
            mq = r.mapping_quality
            if (mq <= min_mapq) if strict_mapq else (mq < min_mapq):
                n["low_mapq"] += 1
                continue
            n["mapped_mapq"] += 1
            src = r.query_name.split(".", 1)[0]
            span = (L or 10 ** 12)
            start, end = r.reference_start % span, r.reference_end % span
            if is_partial(r.query_name):
                key = (src, contig, r.is_reverse, end if r.is_reverse else start, "p")
            else:
                key = (src, contig, r.is_reverse, start, end)
            score = sum(r.query_qualities) if r.query_qualities is not None else 0
            prev = best.get(key)
            if prev is None or score > prev[0]:
                best[key] = (score, r.to_string())
    n["duplicates"] = n["mapped_mapq"] - len(best)
    n["kept"] = len(best)
    tmp = Path(str(out_bam) + ".unsorted.bam")
    with pysam.AlignmentFile(str(tmp), "wb", header=header) as out:
        hdr = out.header
        for _, line in best.values():
            out.write(pysam.AlignedSegment.fromstring(line, hdr))
    pysam.sort("-o", str(out_bam), str(tmp))
    pysam.index(str(out_bam))
    tmp.unlink()
    stats = dict(n)
    stats["duplication_rate"] = n["duplicates"] / n["mapped_mapq"] if n["mapped_mapq"] else None
    return stats


def coverage_span(bam: Path, length: int, circular: bool) -> np.ndarray:
    """Reads spanning each position (aligned bases and deletions)."""
    import pysam
    cov = np.zeros(length + 1, dtype=np.int64)
    with pysam.AlignmentFile(str(bam), "rb") as b:
        for r in b.fetch(until_eof=True):
            if r.is_unmapped:
                continue
            s, e = r.reference_start, r.reference_end
            if circular and e > length:
                cov[s:length] += 1
                cov[0:e - length] += 1
            else:
                cov[s:min(e, length)] += 1
    return cov[:length]


def indel_events(bam: Path, length: int, circular: bool) -> tuple[Counter, Counter]:
    """Counts of deletions {(pos, n): reads} and insertions {(pos, seq): reads};
    an insertion at pos sits immediately before reference position pos."""
    import pysam
    dels, ins = Counter(), Counter()
    with pysam.AlignmentFile(str(bam), "rb") as b:
        for r in b.fetch(until_eof=True):
            if r.is_unmapped or not r.cigartuples:
                continue
            ref, q = r.reference_start, 0
            seq = r.query_sequence
            for op, ln in r.cigartuples:
                if op in (0, 7, 8):
                    ref += ln; q += ln
                elif op == 1:
                    ins[(ref % length if circular else ref, seq[q:q + ln])] += 1
                    q += ln
                elif op == 2:
                    dels[(ref % length if circular else ref, ln)] += 1
                    ref += ln
                elif op == 4:
                    q += ln
    return dels, ins


def majority_indels(dels: Counter, ins: Counter, span: np.ndarray, min_reads: int = 3,
                    majority: float = 0.5) -> list[dict]:
    """Indels supported by >= min_reads and by more than ``majority`` of the
    spanning reads; at most one event per position (the best supported)."""
    cand = []
    for (pos, ln), c in dels.items():
        cand.append({"type": "deletion", "pos": int(pos), "length": int(ln), "reads": int(c)})
    for (pos, s), c in ins.items():
        cand.append({"type": "insertion", "pos": int(pos), "sequence": s, "reads": int(c)})
    out, used = [], set()
    for e in sorted(cand, key=lambda e: -e["reads"]):
        d = int(span[e["pos"] % len(span)]) if len(span) else 0
        e["spanning_reads"] = d
        if e["reads"] >= min_reads and e["reads"] > majority * d and e["pos"] not in used:
            out.append(e)
            used.add(e["pos"])
    return sorted(out, key=lambda e: e["pos"])


def apply_changes(seq: list[str], origin: list[int], subs: dict[int, str], indels: list[dict]) -> tuple[list, list]:
    """New sequence (and, per base, its index in the original reference or -1
    for inserted bases) after substitutions and indels at current positions."""
    dele = {}
    for e in indels:
        if e["type"] == "deletion":
            dele[e["pos"]] = e["length"]
    insb = {e["pos"]: e["sequence"] for e in indels if e["type"] == "insertion"}
    new, org = [], []
    skip = 0
    for i, (b, o) in enumerate(zip(seq, origin)):
        if i in insb and skip == 0:
            new.extend(insb[i]); org.extend([-1] * len(insb[i]))
        if skip:
            skip -= 1
            continue
        if i in dele:
            skip = dele[i] - 1
            continue
        new.append(subs.get(i, b)); org.append(o)
    return new, org
