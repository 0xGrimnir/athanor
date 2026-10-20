"""Map reference coordinates to the true sample genome (for scoring).

The truth genome is aligned to the mapping reference with minimap2
(asm20 preset, --cs); each reference position gets the truth base, '-' if
the truth has a deletion there, or None where no primary alignment covers it.
"""
from __future__ import annotations

import re

CS_OP = re.compile(r"(:[0-9]+|\*[a-z][a-z]|[+-][a-z]+)")


def truth_by_reference_position(paf_text: str, ref_len: int, truth_seq: str) -> list:
    out: list = [None] * ref_len
    truth_seq = truth_seq.upper()
    for line in paf_text.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        tags = dict(t.split(":", 2)[0::2] for t in f[12:] if t.count(":") >= 2)
        if tags.get("tp") not in (None, "P"):
            continue
        cs = next((t[5:] for t in f[12:] if t.startswith("cs:Z:")), None)
        if cs is None:
            continue
        qstart, strand, tstart = int(f[2]), f[4], int(f[7])
        if strand != "+":
            continue  # mitogenomes of related species align on the forward strand
        q, t = qstart, tstart
        for op in CS_OP.findall(cs):
            if op[0] == ":":
                n = int(op[1:])
                for k in range(n):
                    out[t + k] = truth_seq[q + k]
                q += n; t += n
            elif op[0] == "*":
                out[t] = op[2].upper(); q += 1; t += 1
            elif op[0] == "+":
                q += len(op) - 1          # insertion in the truth: no reference position
            elif op[0] == "-":
                for k in range(len(op) - 1):
                    out[t + k] = "-"
                t += len(op) - 1
    return out


def indel_sites(paf_text: str) -> list:
    """Reference positions adjacent to an insertion or deletion of the truth
    relative to the reference (from the same minimap2 alignment)."""
    sites = []
    for line in paf_text.splitlines():
        f = line.split("\t")
        if len(f) < 12 or f[4] != "+":
            continue
        tags = dict(t.split(":", 2)[0::2] for t in f[12:] if t.count(":") >= 2)
        if tags.get("tp") not in (None, "P"):
            continue
        cs = next((t[5:] for t in f[12:] if t.startswith("cs:Z:")), None)
        if cs is None:
            continue
        t = int(f[7])
        for op in CS_OP.findall(cs):
            if op[0] == ":":
                t += int(op[1:])
            elif op[0] == "*":
                t += 1
            elif op[0] == "+":
                sites.append(t)
            elif op[0] == "-":
                sites.extend(range(t, t + len(op) - 1)); t += len(op) - 1
    return sorted(set(sites))
