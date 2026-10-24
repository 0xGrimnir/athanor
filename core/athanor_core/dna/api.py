"""DNA Structure methods for the app: build a helix, melt it, save it."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import melting
from .builder import build, to_pdb
from .structure import DATA, Residue, analyse

ELEMENT = {"P": "P", "O": "O", "N": "N", "C": "C", "H": "H"}


def _element(name: str) -> str:
    return ELEMENT.get(name.lstrip("0123456789")[0], "C")


def _residues(model: dict) -> list[Residue]:
    out: dict = {}
    for chain, num, resname, name, xyz in model["atoms"]:
        out.setdefault((chain, num), Residue(chain, num, "", resname)).atoms[name] = np.asarray(xyz)
    return list(out.values())


def _bonds(atoms) -> list[list[int]]:
    """Covalent bonds: atoms of one residue closer than 1.9 A (2.0 for P), and O3'(i)-P(i+1)."""
    bonds = []
    by_res: dict = {}
    for i, a in enumerate(atoms):
        by_res.setdefault((a[0], a[1]), []).append(i)
    for idx in by_res.values():
        for x in range(len(idx)):
            for y in range(x + 1, len(idx)):
                i, j = idx[x], idx[y]
                lim = 2.0 if "P" in (atoms[i][3], atoms[j][3]) else 1.9
                if np.linalg.norm(atoms[i][4] - atoms[j][4]) < lim:
                    bonds.append([i, j])
    order = list(by_res)
    for (c1, n1), (c2, n2) in zip(order, order[1:]):
        if c1 != c2:
            continue
        o3 = next((i for i in by_res[(c1, n1)] if atoms[i][3] == "O3'"), None)
        p = next((i for i in by_res[(c2, n2)] if atoms[i][3] == "P"), None)
        if o3 is not None and p is not None:
            bonds.append([o3, p])
    return bonds


def validation() -> dict:
    f = DATA / "validation.json"
    return json.loads(f.read_text()) if f.exists() else {}


def m_dna_build(p: dict) -> dict:
    seq = "".join(str(p.get("sequence", "")).split()).upper()
    if not 2 <= len(seq) <= 200:
        raise ValueError("give 2 to 200 base pairs")
    form = p.get("form", "B")
    model = build(seq, form)
    an = analyse(_residues(model))
    info = json.loads((DATA / "steps.json").read_text())[form]
    v = validation()
    atoms = [[c, n, r, a, _element(a), float(x[0]), float(x[1]), float(x[2])] for c, n, r, a, x in model["atoms"]]
    return {"sequence": seq, "form": form, "atoms": atoms, "bonds": _bonds(model["atoms"]),
            "frames": model["frames"], "o3p": model["o3p"],
            "analysis": {"pairs": [{k: x[k] for k in ("pair", "params")} for x in an["pairs"]],
                         "steps": an["steps"], "grooves": an["grooves"]},
            "method": info.get("method", ""), "sources": info.get("sources", []),
            "validation": v.get("summary", {}).get(form)}


def m_dna_melt(p: dict) -> dict:
    seq = p["sequence"]
    ct, na, t = float(p.get("ct", 1e-6)), float(p.get("na", 0.15)), float(p.get("temp", 25))
    r = melting.tm(seq, ct, na)
    z = melting.zipper(seq, t, ct, na)
    v = validation().get("melting")
    return {**r, "temp": t, "fraction_duplex": z["fraction_duplex"], "closed": z["closed"], "checked": v}


def m_dna_export(p: dict) -> dict:
    model = build(p["sequence"], p.get("form", "B"))
    Path(p["path"]).write_text(to_pdb(model))
    return {"path": p["path"], "atoms": len(model["atoms"])}


METHODS = {"dna_build": m_dna_build, "dna_melt": m_dna_melt, "dna_export": m_dna_export}
