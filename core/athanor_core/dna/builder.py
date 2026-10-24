"""Build a DNA duplex atom by atom from base-pair and base-pair-step parameters.

The route is the one 3DNA's rebuild takes (Lu & Olson 2003): base-pair frames are laid down one after
another from the step parameters; each pair's two base frames follow from its base-pair parameters; the
standard bases are placed in their frames; and each nucleotide's sugar-phosphate backbone is attached in
its base's frame from a template measured on high-resolution crystal structures of the same helix form.
Because the backbone rides on the base, the O3'-P links between nucleotides are only approximately
1.6 A; the builder reports how far they are off.

Parameter sets (data/dna/steps.json):
  B   sequence-dependent: the average step parameters of each dimer step in protein-DNA crystal complexes
      (Olson et al., PNAS 95, 11163, 1998, Table 1)
  A, Z  uniform values measured on crystal structures of those forms (recorded with their sources)
"""
from __future__ import annotations

import json

import numpy as np

from .geometry import next_frame, pair_frames
from .structure import DATA, standard_bases

COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def _data(name):
    return json.loads((DATA / name).read_text())


def complement_step(p):
    """Parameters of the same step read along the other strand: shift and tilt change sign."""
    p = list(p)
    return [-p[0], p[1], p[2], -p[3], p[4], p[5]]


def step_params(seq: str, form: str) -> list[list[float]]:
    d = _data("steps.json")[form]
    out = []
    for k in range(len(seq) - 1):
        if d.get("uniform"):
            alt = d["uniform"] if not isinstance(d["uniform"], dict) else d["uniform"][seq[k:k + 2]]
            out.append(list(alt))
            continue
        st = seq[k:k + 2]
        if st in d["dimers"]:
            out.append(list(d["dimers"][st]))
        else:
            rc = COMP[st[1]] + COMP[st[0]]
            out.append(complement_step(d["dimers"][rc]))
    return out


def pair_params(seq: str, form: str) -> list[list[float]]:
    d = _data("steps.json")[form]
    return [list(d["pair"]) for _ in seq]


def build(seq: str, form: str = "B", steps=None, pairs=None) -> dict:
    """Atoms of the duplex: {"atoms": [(chain, resnum, resname, atom, xyz)], "frames": [...],
    "o3p": O3'-P distances}. steps/pairs override the parameter tables."""
    seq = "".join(seq.split()).upper()
    if form == "Z" and any(seq[i:i + 2] not in ("CG", "GC") for i in range(len(seq) - 1)):
        raise ValueError("Z-DNA is built only for alternating C and G")
    std = standard_bases()
    tmpl = _data("templates.json")[form]
    steps = steps if steps is not None else step_params(seq, form)
    pairs = pairs if pairs is not None else pair_params(seq, form)
    R, o = np.eye(3), np.zeros(3)
    frames = [(R, o)]
    for p in steps:
        R, o = next_frame(R, o, p)
        frames.append((R, o))
    n = len(seq)
    atoms = []
    strand2 = []
    for k, (b, (Rm, om), pp) in enumerate(zip(seq, frames, pairs)):
        (Ri, oi), (Rii, oii) = pair_frames(Rm, om, pp)
        for base, Rb, ob, chain, num in ((b, Ri, oi, "A", k + 1), (COMP[b], Rii, oii, "B", 2 * n - k)):
            t = tmpl.get(_key(form, base, k if chain == "A" else n - 1 - k, n)) or tmpl[base]
            res = []
            for name, xyz in t.items():
                if name in std[base]:
                    xyz = std[base][name]
                res.append((chain, num, "D" + base, name, ob + Rb @ np.asarray(xyz)))
            (atoms if chain == "A" else strand2).append(res)
    ordered = [a for r in atoms for a in r] + [a for r in reversed(strand2) for a in r]
    o3p = []
    for res_list in (atoms, list(reversed(strand2))):
        for r1, r2 in zip(res_list, res_list[1:]):
            o3 = next((a[4] for a in r1 if a[3] == "O3'"), None)
            p = next((a[4] for a in r2 if a[3] == "P"), None)
            if o3 is not None and p is not None:
                o3p.append(float(np.linalg.norm(o3 - p)))
    return {"sequence": seq, "form": form, "atoms": ordered, "steps": steps, "pairs": pairs,
            "frames": [(R.tolist(), o.tolist()) for R, o in frames], "o3p": o3p}


def _key(form, base, pos, n):
    """Z-DNA alternates two nucleotide conformations (C anti, G syn); templates are keyed by base."""
    return base


def to_pdb(model: dict) -> str:
    lines = []
    for k, (chain, num, res, name, xyz) in enumerate(model["atoms"], 1):
        el = name.lstrip("0123456789")[0]
        nm = f" {name:<3s}" if len(name) < 4 else name
        lines.append(f"ATOM  {k:5d} {nm} {res:>3s} {chain}{num:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                     f"  1.00  0.00          {el:>2s}")
    lines.append("END")
    return "\n".join(lines) + "\n"
