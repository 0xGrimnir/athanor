"""DNA Structure methods for the app: build a helix, melt it, save it."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import md, melting
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


def _load(name) -> dict:
    f = DATA / name
    return json.loads(f.read_text()) if f.exists() else {}


def validation(form: str) -> str | None:
    """What CI measured for this form's builder (data/dna/validation.json, written by tools/dnadata)."""
    v = _load("validation.json")
    r = v.get("summary", {}).get("rebuild", {}).get(form)
    if not r or not r.get("n"):
        return None
    how = ("built from sequence alone with these parameters" if form == "B"
           else "each built from parameters measured on the others")
    text = (f"{r['n']} crystal {form}-DNA duplexes from the PDB, {how}, are matched to a median all-atom RMSD of "
            f"{r['generic_median']:.2f} \u00c5 (90th percentile {r['generic_p90']:.2f} \u00c5).")
    if all(x["pass"] for x in v.get("results", [])):
        text += " Base-pair, step and groove measurements equal 3DNA's on every comparable crystal."
    return text


def melting_check() -> dict | None:
    m = _load("melting_validation.json").get("errors", {}).get("owczarzy2004", {})
    a = m.get("all")
    if not a:
        return None
    return {"mean_abs": round(a["mean_abs"], 2), "mean": round(a["mean"], 2), "n": a["n"],
            "salts": len([k for k in m if k != "all"]), "duplexes": a["n"] // max(1, len([k for k in m if k != "all"]))}


def _seq(p) -> str:
    seq = "".join(str(p.get("sequence", "")).split()).upper()
    if not 2 <= len(seq) <= 200:
        raise ValueError("give 2 to 200 base pairs")
    return seq


def _response(seq, form, atoms, frames, o3p, extra=None) -> dict:
    """What the app needs to draw and describe a duplex given as Athanor atoms."""
    an = analyse(_residues({"atoms": atoms}))
    info = json.loads((DATA / "steps.json").read_text())[form]
    out = [[c, n, r, a, _element(a), float(x[0]), float(x[1]), float(x[2])] for c, n, r, a, x in atoms]
    return {"sequence": seq, "form": form, "atoms": out, "bonds": _bonds(atoms),
            "frames": frames, "o3p": o3p,
            "analysis": {"pairs": [{k: x[k] for k in ("pair", "params")} for x in an["pairs"]],
                         "steps": an["steps"], "grooves": an["grooves"]},
            "method": info.get("method", ""), "sources": info.get("sources", []),
            "validation": validation(form), **(extra or {})}


def m_dna_build(p: dict) -> dict:
    seq, form = _seq(p), p.get("form", "B")
    model = build(seq, form)
    return _response(seq, form, model["atoms"], model["frames"], model["o3p"],
                     {"relaxed": False, "md_available": md.available()})


def m_dna_relax(p: dict) -> dict:
    """The built helix with its backbone closed by an OpenMM energy minimisation (bases restrained)."""
    seq, form = _seq(p), p.get("form", "B")
    model = build(seq, form)
    r = md.relax(model["atoms"], salt=float(p.get("salt", 0.15)))
    return _response(seq, form, r["atoms"], model["frames"], r["o3p"],
                     {"relaxed": True, "relax": {k: r[k] for k in ("energy_kcal", "base_shift_max", "base_shift_rms",
                                                                    "seconds", "platform", "method")},
                      "md_available": True})


_MD: dict = {}


def m_dna_md_start(p: dict) -> dict:
    """Start molecular dynamics of the built helix; returns the starting model and a run id."""
    seq, form = _seq(p), p.get("form", "B")
    model = build(seq, form)
    _MD.clear()                                    # one simulation at a time
    d = md.Dynamics(model["atoms"], salt=float(p.get("salt", 0.15)), temperature=float(p.get("temperature", 300.0)),
                    prefer=p.get("platform"))
    f = d.run(float(p.get("ps", 1.0)))
    rid = f"md{id(d)}"
    _MD[rid] = {"dyn": d, "sum": None, "sq": None, "n": 0, "seq": seq, "form": form}
    resp = _response(seq, form, f["atoms"], model["frames"], md._o3p(f["atoms"]))
    resp.update({"run": rid, "frame": _frame(rid, f), "method": d.describe()})
    return resp


def _frame(rid, f) -> dict:
    st = _MD[rid]
    an = analyse(_residues({"atoms": f["atoms"]}))
    steps = np.array([s["params"] for s in an["steps"]]) if an["steps"] else np.zeros((0, 6))
    grooves = np.array([g["minor"] if g["minor"] is not None else np.nan for g in an["grooves"]], dtype=float)
    vals = np.concatenate([steps.ravel(), grooves])
    if st["sum"] is None or len(st["sum"]) != len(vals):
        st["sum"], st["sq"], st["cnt"], st["n"] = np.zeros_like(vals), np.zeros_like(vals), np.zeros_like(vals), 0
    ok = ~np.isnan(vals)
    st["sum"][ok] += vals[ok]
    st["sq"][ok] += vals[ok] ** 2
    st["cnt"][ok] += 1
    st["n"] += 1
    cnt = np.maximum(st["cnt"], 1)
    mean = st["sum"] / cnt
    sd = np.sqrt(np.maximum(st["sq"] / cnt - mean ** 2, 0.0))
    k = steps.size
    return {"coords": [[round(float(v), 3) for v in a[4]] for a in f["atoms"]], "time_ps": f["time_ps"],
            "temperature": f["temperature"], "potential_kcal": f["potential_kcal"], "ns_per_day": f["ns_per_day"],
            "platform": f["platform"], "n_frames": st["n"],
            "steps": steps.tolist(), "minor": [None if np.isnan(g) else float(g) for g in grooves],
            "mean_steps": mean[:k].reshape(-1, 6).tolist(), "sd_steps": sd[:k].reshape(-1, 6).tolist(),
            "mean_minor": [None if c == 0 else float(m) for m, c in zip(mean[k:], st["cnt"][k:])],
            "sd_minor": [None if c == 0 else float(x) for x, c in zip(sd[k:], st["cnt"][k:])]}


def m_dna_md_step(p: dict) -> dict:
    rid = p["run"]
    if rid not in _MD:
        raise ValueError("that simulation has ended")
    f = _MD[rid]["dyn"].run(float(p.get("ps", 5.0)))
    return _frame(rid, f)


def m_dna_md_stop(p: dict) -> dict:
    _MD.pop(p.get("run"), None)
    return {"stopped": True}


def m_dna_variants(p: dict) -> dict:
    """The variants the app can show in a helix (reference sequence either side, data/dna/flanks.json)."""
    f = _load("flanks.json")
    return {"variants": sorted(f.get("variants", {})), "flank": f.get("flank")}


def m_dna_variant(p: dict) -> dict:
    """A B-DNA helix of the GRCh37 forward-strand sequence around a variant, with ``allele`` at its centre."""
    f = _load("flanks.json")
    v = f.get("variants", {}).get(p.get("rsid"))
    if not v:
        raise ValueError(f"no reference sequence for {p.get('rsid')}")
    allele = str(p.get("allele") or v["ref"]).upper()
    if allele not in "ACGT" or len(allele) != 1:
        raise ValueError("a single base is needed")
    seq = v["left"] + allele + v["right"]
    r = m_dna_build({"sequence": seq, "form": "B"})
    r["variant"] = {"rsid": p["rsid"], "chrom": v["chrom"], "pos37": v["pos37"], "ref": v["ref"], "alts": v["alts"],
                    "allele": allele, "index": len(v["left"]), "source": f.get("source")}
    return r


def m_dna_melt(p: dict) -> dict:
    seq = p["sequence"]
    ct, na, t = float(p.get("ct", 1e-6)), float(p.get("na", 0.15)), float(p.get("temp", 25))
    r = melting.tm(seq, ct, na)
    z = melting.zipper(seq, t, ct, na)
    return {**r, "temp": t, "fraction_duplex": z["fraction_duplex"], "closed": z["closed"], "checked": melting_check()}


def m_dna_export(p: dict) -> dict:
    model = build(p["sequence"], p.get("form", "B"))
    Path(p["path"]).write_text(to_pdb(model))
    return {"path": p["path"], "atoms": len(model["atoms"])}


METHODS = {"dna_build": m_dna_build, "dna_relax": m_dna_relax, "dna_md_start": m_dna_md_start,
           "dna_md_step": m_dna_md_step, "dna_md_stop": m_dna_md_stop, "dna_melt": m_dna_melt, "dna_export": m_dna_export,
           "dna_variants": m_dna_variants, "dna_variant": m_dna_variant}
