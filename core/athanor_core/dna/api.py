"""DNA Structure methods for the app: build a helix, melt it, save it."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import damage, energy, md, melting, shape
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
    """Save the structure on screen (atoms as the app holds them), or a fresh build of a sequence."""
    if p.get("atoms"):
        model = {"atoms": [(a[0], int(a[1]), a[2], a[3], np.array(a[5:8], dtype=float)) for a in p["atoms"]]}
    else:
        model = build(p["sequence"], p.get("form", "B"))
    Path(p["path"]).write_text(to_pdb(model))
    return {"path": p["path"], "atoms": len(model["atoms"])}


# ------------------------------------------------------------------ shape, energies, damage, scales
def m_dna_shape(p: dict) -> dict:
    """Deep DNAshape's predicted shape of the sequence as free DNA in solution."""
    seq = _seq(p)
    feats = p.get("features") or ["MGW", "ProT", "Roll", "HelT", "EP", "MGW-FL"]
    layer = int(p.get("layer", 4))
    out = {f: [round(float(x), 3) for x in shape.predict(f, seq, layer)] for f in feats if f in shape.available()}
    v = _load("dnashape/validation.json")
    return {"sequence": seq, "layer": layer, "shape": out, "units": {f: shape.UNITS.get(f.replace("-FL", ""), "deg") for f in out},
            "names": {f: shape.NAMES.get(f.replace("-FL", ""), f) + (" (fluctuation)" if f.endswith("-FL") else "") for f in out},
            "per": {f: "step" if f.replace("-FL", "") in shape.INTER else "pair" for f in out},
            "method": "Deep DNAshape (Li, Chiu & Rohs, Nature Commun. 15, 1243, 2024), trained on Monte Carlo simulations "
                      "of free DNA; run by Athanor's NumPy implementation of the published network and weights "
                      f"(BSD-3-Clause), {layer} bp of flanking sequence considered.",
            "checked": "identical to the original TensorFlow model in CI" if v.get("equal_to_original") else None}


_RELAXED: dict = {}
PROGRESS = None            # set by the server: sends progress notifications to the app


def _relaxed(seq, form, salt=0.15):
    key = (seq, form, round(float(salt), 3))
    if key not in _RELAXED:
        _RELAXED.clear()
        _RELAXED[key] = md.relax(build(seq, form)["atoms"], salt=salt)
    return _RELAXED[key]


def m_dna_energies(p: dict) -> dict:
    """Pairing energy of every base pair and stacking energy of every step (GFN2-xTB, gas phase)."""
    seq, form = _seq(p), p.get("form", "B")
    if len(seq) > 60:
        raise ValueError("energies are computed for up to 60 base pairs")
    r = _relaxed(seq, form)
    h = energy.helix_energies(r["atoms"], r["hydrogens"], len(seq), "gfn2-xtb", on_step=PROGRESS)
    return {"sequence": seq, "form": form, "pairs": [round(x, 2) for x in h["pairs"]], "steps": [round(x, 2) for x in h["steps"]],
            "method": "GFN2-xTB interaction energies of the capped bases of the relaxed helix, in vacuum", "checked": _load("energy_validation.json")}


def m_dna_energy_one(p: dict) -> dict:
    """One pairing (kind 'pair', index k) or stacking (kind 'step', index k) energy with DFT."""
    seq, form = _seq(p), p.get("form", "B")
    r = _relaxed(seq, form)
    fr = energy.fragments(r["atoms"], r["hydrogens"], len(seq))
    k = int(p["index"])
    if p.get("kind") == "step":
        a, b = list(fr[k]), list(fr[k + 1])
    else:
        a, b = [fr[k][0]], [fr[k][1]]
    e = energy.interaction_dft(a, b)
    return {"kind": p.get("kind", "pair"), "index": k, "energy": round(e, 2),
            "method": "B3LYP-D3(BJ)/def2-TZVP with counterpoise correction (Psi4), capped bases in vacuum"}


def m_dna_damage(p: dict) -> dict:
    """One ancient molecule: overhangs and uracils drawn from the Briggs model or a measured profile."""
    seq = _seq(p)
    rng = np.random.default_rng(int(p["seed"]) if p.get("seed") is not None else None)
    lam, dd, ds = float(p.get("lambda", 0.4)), float(p.get("delta_d", 0.01)), float(p.get("delta_s", 0.3))
    prof = p.get("profile")
    m = damage.sample(seq, lam, dd, ds, rng, profile=prof)
    resp = _response(seq, "B", m["atoms"], m["frames"], md._o3p(m["atoms"]) if hasattr(md, "_o3p") else [])
    resp.update({"damage": {"overhang5": m["overhang5"], "overhang3": m["overhang3"], "uracils": m["uracils"],
                            "p_uracil": m["p_uracil"], "model": m["model"],
                            "expected_ct": damage.p_uracil(np.arange(min(25, len(seq))), lam, dd, ds).round(4).tolist(),
                            "parameters": {"lambda": lam, "delta_d": dd, "delta_s": ds}}})
    return resp


def _scale(name):
    return _load(f"scales/{name}.json")


def m_dna_scales(p: dict) -> dict:
    """What the zoom through the scales offers, with the checks CI ran."""
    v = _load("scales/validation.json")
    cx = _scale("complexes") or []
    return {"available": bool(v), "validation": v,
            "complexes": [{k: c.get(k) for k in ("pdb", "name", "citation", "published_bend", "at_least", "bend", "note",
                                                  "largest_roll", "n_pairs", "sequence")} for c in cx]}


def m_dna_nucleosome(p: dict) -> dict:
    """A sequence threaded onto the 1KX5 nucleosome's DNA path, with the histone octamer."""
    nuc = _scale("nucleosome")
    if not nuc:
        raise ValueError("nucleosome data are not installed")
    ref = nuc["sequence"]
    n = len(ref)
    seq = "".join(str(p.get("sequence") or "").split()).upper()
    # the user's sequence goes in the middle, centred on the dyad; the rest keeps the crystal's sequence
    if seq:
        seq = seq[:n]
        start = max(0, nuc["dyad_index"] - len(seq) // 2)
        start = min(start, n - len(seq))
        full = ref[:start] + seq + ref[start + len(seq):]
    else:
        start, full = 0, ref
    m = build(full, "B", steps=nuc["steps"], pairs=nuc["pairs"])
    R0, o0 = np.array(nuc["frame0"][0]), np.array(nuc["frame0"][1])
    atoms = [(c, k, r, a, R0 @ np.asarray(x) + o0) for c, k, r, a, x in m["atoms"]]
    resp = _response(full, "B", atoms, m["frames"], m["o3p"])
    resp.update({"scale": "nucleosome", "histones": nuc["histones"], "superhelix": nuc["superhelix"],
                 "superhelix_axis": nuc.get("superhelix_axis"),
                 "user": {"start": start, "length": len(seq)} if seq else None, "dyad_index": nuc["dyad_index"],
                 "pdb": nuc["pdb"], "title": nuc["title"], "citation": nuc["citation"],
                 "method": "The 147 base-pair steps and pairs of the 1KX5 crystal structure, measured as in 3DNA, with your "
                           "sequence's bases placed in them (3DNA's 'mutate' approach): the DNA keeps the crystal's path, "
                           "so this shows where your sequence would sit, not how well it positions a nucleosome."})
    return resp


def m_dna_complex(p: dict) -> dict:
    cx = {c["pdb"]: c for c in (_scale("complexes") or [])}
    c = cx.get(p.get("pdb"))
    if not c:
        raise ValueError("unknown complex")
    atoms = [(a[0], a[1], a[2], a[3], np.array(a[4:7])) for a in c["dna"]]
    seq = c["sequence"]
    resp = _response(seq, "B", atoms, [], [])
    resp.update({"scale": "complex", "protein": c["protein"], **{k: c.get(k) for k in ("pdb", "name", "citation", "published_bend",
                                                                                      "at_least", "bend", "note", "largest_roll")}})
    return resp


def m_dna_fibre(p: dict) -> dict:
    f = _scale("fibre")
    if not f:
        raise ValueError("fibre data are not installed")
    return {"scale": "fibre", **f}


def m_dna_nucleus(p: dict) -> dict:
    f = _scale("nucleus")
    if not f:
        raise ValueError("nucleus data are not installed")
    out = {"scale": "nucleus", **f}
    chrom, pos = p.get("chrom"), p.get("pos")
    if chrom and pos:
        b = int(pos) // f["bin"]
        hits = []
        for hom, v in f["homologs"].items():
            if hom.split("(")[0] == str(chrom).replace("chr", "") and b in v["bins"]:
                hits.append({"homolog": hom, "xyz": v["xyz"][v["bins"].index(b)]})
        out["highlight"] = {"chrom": chrom, "pos": pos, "bin": b, "hits": hits}
    return out


METHODS = {"dna_build": m_dna_build, "dna_relax": m_dna_relax, "dna_md_start": m_dna_md_start,
           "dna_md_step": m_dna_md_step, "dna_md_stop": m_dna_md_stop, "dna_shape": m_dna_shape,
           "dna_energies": m_dna_energies, "dna_energy_one": m_dna_energy_one, "dna_damage": m_dna_damage,
           "dna_scales": m_dna_scales, "dna_nucleosome": m_dna_nucleosome, "dna_complex": m_dna_complex,
           "dna_fibre": m_dna_fibre, "dna_nucleus": m_dna_nucleus, "dna_melt": m_dna_melt, "dna_export": m_dna_export,
           "dna_variants": m_dna_variants, "dna_variant": m_dna_variant}
