"""Compare a reconstructed sequence with another (e.g. a deposited GenBank record).

minimap2 aligns ours (query) to theirs (target) and reports every
difference through its cs tag (Li, Bioinformatics 34, 3094, 2018). Sites
where either sequence is N or another ambiguity code are counted
separately: they are not disagreements, only sites one side left uncalled.
"""
from __future__ import annotations

import re
from pathlib import Path

CS = re.compile(r":(\d+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)")
ACGT = set("acgt")


def parse_cs(cs: str, t_start: int, q_start: int, q_reverse: bool = False) -> dict:
    """Walk a short-form cs string. Returns matches, substitutions between
    called bases, uncalled sites and indels, with target/query positions."""
    t, q = t_start, q_start
    out = {"matches": 0, "substitutions": [], "uncalled_ours": 0, "uncalled_theirs": 0,
           "insertions_in_ours": [], "deletions_in_ours": []}
    for m in CS.finditer(cs):
        if m.group(1):
            n = int(m.group(1)); out["matches"] += n; t += n; q += n
        elif m.group(2):
            tb, qb = m.group(2), m.group(3)
            if tb not in ACGT:
                out["uncalled_theirs"] += 1
            elif qb not in ACGT:
                out["uncalled_ours"] += 1
            else:
                out["substitutions"].append({"theirs_pos": t, "ours_pos": q, "theirs": tb.upper(), "ours": qb.upper()})
            t += 1; q += 1
        elif m.group(4):
            s = m.group(4)
            out["insertions_in_ours"].append({"theirs_pos": t, "ours_pos": q, "sequence": s.upper()}); q += len(s)
        elif m.group(5):
            s = m.group(5)
            out["deletions_in_ours"].append({"theirs_pos": t, "ours_pos": q, "sequence": s.upper()}); t += len(s)
    return out


def compare_sequences(step, work: Path, ours: str, theirs: str, theirs_name: str = "deposited") -> dict:
    from .bamtools import write_fasta
    from .pipeline import BIO
    d = work / "compare"
    d.mkdir(parents=True, exist_ok=True)
    write_fasta(d / "ours.fa", {"athanor": ours})
    write_fasta(d / "theirs.fa", {theirs_name: theirs})
    paf = step.run("compare with deposited sequence (minimap2)",
                   [BIO / "bin/minimap2", "-cx", "asm5", "--cs", d / "theirs.fa", d / "ours.fa"], d)
    total = {"matches": 0, "substitutions": [], "uncalled_ours": 0, "uncalled_theirs": 0,
             "insertions_in_ours": [], "deletions_in_ours": []}
    alignments = []
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        tags = {x[:2]: x[5:] for x in f[12:]}
        if tags.get("tp", "P") != "P" or "cs" not in tags:
            continue
        q_start, strand, t_start = int(f[2]), f[4], int(f[7])
        alignments.append({"ours": [q_start, int(f[3])], "theirs": [t_start, int(f[8])], "strand": strand})
        r = parse_cs(tags["cs"], t_start, q_start)
        for k, v in r.items():
            total[k] = total[k] + v
    called_both = total["matches"] + len(total["substitutions"])
    aligned_theirs = sum(a["theirs"][1] - a["theirs"][0] for a in alignments)
    return {"theirs_name": theirs_name, "theirs_length": len(theirs), "ours_length": len(ours),
            "alignments": alignments,
            "theirs_aligned_fraction": aligned_theirs / len(theirs) if theirs else None,
            "sites_called_by_both": called_both, "matches": total["matches"],
            "substitutions": total["substitutions"],
            "agreement": total["matches"] / called_both if called_both else None,
            "uncalled_ours": total["uncalled_ours"], "uncalled_theirs": total["uncalled_theirs"],
            "theirs_ambiguity_codes": sum(1 for b in theirs.upper() if b not in "ACGT"),
            "ours_N": ours.upper().count("N"),
            "insertions_in_ours": total["insertions_in_ours"], "deletions_in_ours": total["deletions_in_ours"]}
