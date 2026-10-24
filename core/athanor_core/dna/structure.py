"""Nucleic-acid structures: reading mmCIF/PDB files, fitting base reference frames, and analysing a
duplex into base-pair, step and groove parameters.

Base frames are fitted the way 3DNA fits them (Lu & Olson 2003): the ring atoms of the standard base
(placed in its reference frame, Olson et al. 2001) are superposed on the observed ring atoms by least
squares, and the fitted rotation and translation are the base's frame.
"""
from __future__ import annotations

import json
import shlex
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .geometry import FLIP, bz_check, kabsch, step_parameters

PURINE_RING = ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"]
PYRIMIDINE_RING = ["N1", "C2", "N3", "C4", "C5", "C6"]
BASE_OF = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U", "A": "A", "C": "C", "G": "G", "U": "U", "T": "T"}
DATA = Path(__file__).resolve().parents[1] / "data" / "dna"


def ring(base: str) -> list[str]:
    return PURINE_RING if base in "AG" else PYRIMIDINE_RING


@dataclass
class Residue:
    chain: str
    number: int
    insertion: str
    name: str
    atoms: dict = field(default_factory=dict)        # atom name -> xyz

    @property
    def base(self) -> str | None:
        return BASE_OF.get(self.name)


def _float(s):
    try:
        return float(s)
    except ValueError:
        return None


def read_cif(text: str, model: int | None = None) -> list[Residue]:
    """Nucleotide residues from an mmCIF file (first model, first alternate location)."""
    lines = text.splitlines()
    cols, rows, i = [], [], 0
    while i < len(lines):
        if lines[i].strip() == "loop_" and i + 1 < len(lines) and lines[i + 1].startswith("_atom_site."):
            i += 1
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].split(".", 1)[1].strip())
                i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(("_", "loop_", "#")):
                rows.append(shlex.split(lines[i], posix=True))
                i += 1
            break
        i += 1
    if not cols:
        raise ValueError("no atom_site table in the file")
    c = {k: n for n, k in enumerate(cols)}

    def get(r, *names):
        for nme in names:
            if nme in c and c[nme] < len(r):
                v = r[c[nme]]
                return None if v in ("?", ".") else v
        return None

    out: dict = {}
    first_model = None
    for r in rows:
        m = get(r, "pdbx_PDB_model_num")
        if first_model is None:
            first_model = m
        if (model is None and m != first_model) or (model is not None and m != str(model)):
            continue
        alt = get(r, "label_alt_id")
        if alt not in (None, "A", "1"):
            continue
        name = get(r, "auth_comp_id", "label_comp_id")
        if name not in BASE_OF:
            continue
        chain = get(r, "auth_asym_id", "label_asym_id")
        num = int(get(r, "auth_seq_id", "label_seq_id"))
        ins = get(r, "pdbx_PDB_ins_code") or ""
        atom = get(r, "auth_atom_id", "label_atom_id").replace("*", "'")
        xyz = np.array([_float(get(r, "Cartn_x")), _float(get(r, "Cartn_y")), _float(get(r, "Cartn_z"))])
        key = (chain, num, ins)
        if key not in out:
            out[key] = Residue(chain, num, ins, name)
        out[key].atoms.setdefault(atom, xyz)
    return list(out.values())


def read_pdb(text: str) -> list[Residue]:
    """Nucleotide residues from a PDB-format file (first model, first alternate location)."""
    out: dict = {}
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        alt = line[16].strip()
        if alt not in ("", "A", "1"):
            continue
        name = line[17:20].strip()
        if name not in BASE_OF:
            continue
        chain, num, ins = line[21].strip(), int(line[22:26]), line[26].strip()
        atom = line[12:16].strip().replace("*", "'")
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        key = (chain, num, ins)
        if key not in out:
            out[key] = Residue(chain, num, ins, name)
        out[key].atoms.setdefault(atom, xyz)
    return list(out.values())


def standard_bases() -> dict:
    """Standard base coordinates in their reference frames, {base: {atom: xyz}}."""
    d = json.loads((DATA / "standard_bases.json").read_text())
    return {b: {a: np.array(x) for a, x in atoms.items()} for b, atoms in d["bases"].items()}


def base_frame(res: Residue, std: dict):
    """(R, origin, rmsd) of the residue's base, or None if ring atoms are missing."""
    b = res.base
    names = ring(b)
    if not all(n in res.atoms for n in names):
        return None
    P = np.array([std[b][n] for n in names])
    Q = np.array([res.atoms[n] for n in names])
    R, t, rmsd = kabsch(P, Q)
    return R, t, rmsd


def chains(residues: list[Residue]) -> dict:
    by: dict = {}
    for r in residues:
        by.setdefault(r.chain, []).append(r)
    return by


def find_pairs(residues: list[Residue], std: dict, max_origin: float = 4.0):
    """Watson-Crick-like pairs of a two-stranded duplex: for each base of the first strand, the base of
    the second strand whose frame origin is nearest, provided the two z axes point in roughly opposite
    directions (antiparallel strands) and the origins are within max_origin A."""
    by = chains(residues)
    if len(by) < 2:
        raise ValueError("a duplex needs two strands")
    names = list(by)[:2]
    s1, s2 = by[names[0]], by[names[1]]
    f1 = [base_frame(r, std) for r in s1]
    f2 = [base_frame(r, std) for r in s2]
    pairs = []
    for i, a in enumerate(f1):
        if a is None:
            continue
        best = None
        for j, b in enumerate(f2):
            if b is None:
                continue
            d = np.linalg.norm(a[1] - b[1])
            if d < max_origin and np.dot(a[0][:, 2], b[0][:, 2]) < -0.5 and (best is None or d < best[0]):
                best = (d, j)
        if best:
            pairs.append((i, best[1]))
    return s1, s2, f1, f2, pairs


COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}


TURN = np.diag([-1.0, 1.0, -1.0])


def analyse(residues: list[Residue], std: dict | None = None) -> dict:
    """Base-pair and step parameters, and minor/major groove widths, for a duplex.

    As in 3DNA (analyze, ref_frames): each strand II base frame is reversed (y and z negated) to
    run with strand I; and if every Watson-Crick base's z axis points against the direction to the next
    pair - a left-handed, Z-DNA structure - all base frames are turned 180 degrees about y (x and z
    negated) before any parameter is computed. Steps between Watson-Crick pairs then get 3DNA's B-Z
    junction check."""
    std = std or standard_bases()
    s1, s2, f1, f2, pairs = find_pairs(residues, std)
    RI = [f1[i][0] for i, _ in pairs]
    oI = [f1[i][1] for i, _ in pairs]
    RII = [f2[j][0] @ FLIP for _, j in pairs]
    oII = [f2[j][1] for _, j in pairs]
    wc = [COMPLEMENT.get(s1[i].base) == s2[j].base or (s1[i].base, s2[j].base) == ("A", "U") for i, j in pairs]
    n = len(pairs)
    left = total = 0
    for Rs, os_ in ((RI, oI), (RII, oII)):
        d = None
        for k in range(n):
            if k < n - 1:
                d = os_[k + 1] - os_[k]
            if d is not None and wc[k]:
                total += 1
                left += float(np.dot(d, Rs[k][:, 2])) < 0
    z_dna = total > 0 and left == total
    if z_dna:
        RI = [R @ TURN for R in RI]
        RII = [R @ TURN for R in RII]
    bps = []
    for k, (i, j) in enumerate(pairs):
        p, Rm, om = step_parameters(RII[k], oII[k], RI[k], oI[k])
        bps.append({"i": i, "j": j, "pair": f"{s1[i].base}-{s2[j].base}",
                    "residues": [f"{s1[i].chain}{s1[i].number}{s1[i].insertion}",
                                 f"{s2[j].chain}{s2[j].number}{s2[j].insertion}"],
                    "watson_crick": wc[k], "params": p.tolist(), "frame": (Rm.tolist(), om.tolist()),
                    "fit_rmsd": [f1[i][2], f2[j][2]]})
    steps = []
    for a, b in zip(bps, bps[1:]):
        if b["i"] != a["i"] + 1:
            continue
        Ra, oa = (np.array(x) for x in a["frame"])
        Rb, ob = (np.array(x) for x in b["frame"])
        if a["watson_crick"] and b["watson_crick"]:
            Ra, Rb = bz_check(Ra, oa, Rb, ob)
        p, _, _ = step_parameters(Ra, oa, Rb, ob)
        steps.append({"step": f"{a['pair'][0]}{b['pair'][0]}/{b['pair'][2]}{a['pair'][2]}", "params": p.tolist()})
    return {"strands": [s1[0].chain if s1 else None, s2[0].chain if s2 else None],
            "sequence": "".join(s1[a["i"]].base for a in bps), "left_handed": z_dna,
            "pairs": bps, "steps": steps, "grooves": grooves(s1, s2, pairs)}


def grooves(s1, s2, pairs) -> list[dict]:
    """Minor and major groove widths per step from cross-strand phosphorus distances, as 3DNA reports
    them (direct P-P distances; El Hassan & Calladine, J. Mol. Biol. 282, 331, 1998). With base pairs
    numbered from 1 and P_I(k), P_II(k) the phosphorus of each strand's nucleotide in pair k, step s
    (between pairs s and s+1) has
        minor = (|P_I(s+2) - P_II(s-2)| + |P_I(s+3) - P_II(s-1)|) / 2,   major = |P_I(s-1) - P_II(s+2)|
    defined for 3 <= s <= n-3. Subtract 5.8 A (the phosphate groups' van der Waals radii) for the
    width of the groove itself; both are returned."""
    n = len(pairs)
    p1 = [s1[i].atoms.get("P") for i, _ in pairs]
    p2 = [s2[j].atoms.get("P") for _, j in pairs]

    def P(arr, k):                 # 1-based pair index
        return arr[k - 1] if 1 <= k <= n else None

    def dist(a, b):
        return float(np.linalg.norm(a - b)) if a is not None and b is not None else None

    out = []
    for st in range(1, n):
        row = {"step": st - 1, "minor_pp": None, "major_pp": None, "minor": None, "major": None}
        if 3 <= st <= n - 3:
            d1, d2 = dist(P(p1, st + 2), P(p2, st - 2)), dist(P(p1, st + 3), P(p2, st - 1))
            if d1 is not None and d2 is not None:
                row["minor_pp"] = (d1 + d2) / 2
                row["minor"] = row["minor_pp"] - 5.8
            d = dist(P(p1, st - 1), P(p2, st + 2))
            if d is not None:
                row["major_pp"] = d
                row["major"] = d - 5.8
        out.append(row)
    return out
