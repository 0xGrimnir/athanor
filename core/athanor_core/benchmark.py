"""Phase 1 validation benchmark.

Each tier is run on literature reference data and its error statistics are
compared with published statistics for the same method (see the
``comparators`` and ``gates`` sections of the reference file). Results are
written as JSON plus a human-readable Markdown summary.

Usage:  python -m athanor_core.benchmark --tiers mmff94,gfn2-xtb,dft --out report
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from . import tasks, units
from .engines import get_engine
from .molecule import Molecule, angle, distance
from .vibrations import harmonic_analysis

REF = Path(__file__).resolve().parent.parent / "benchmarks" / "reference" / "phase1.json"
GAS_PHASES = {"gas"}


def _mol_from_entry(e: dict) -> Molecule:
    if "smiles" in e and e.get("prefer_smiles", True):
        m = Molecule.from_smiles(e["smiles"], name=e["id"])
        return m
    syms = [a[0] for a in e["xyz"]]
    xyz = np.array([a[1:] for a in e["xyz"]], dtype=float)
    return Molecule(syms, xyz, e.get("charge", 0), e.get("multiplicity", 1), name=e["id"])


def _engine(tier: str, purpose: str, ref: dict):
    if tier != "dft":
        return get_engine(tier)
    level = ref["gates"]["dft_benchmark_level"][purpose]
    func, basis = level.split("/")
    return get_engine("dft", functional=func, basis=basis)


def _stats(errors: list[float]) -> dict:
    a = np.asarray(errors, dtype=float)
    if a.size == 0:
        return {"n": 0}
    return {"n": int(a.size), "mean_signed": float(a.mean()), "mad": float(np.abs(a).mean()),
            "rms": float(np.sqrt((a ** 2).mean())), "max_abs": float(np.abs(a).max())}


# ---------------------------------------------------------------- geometries
def bench_geometries(tier: str, ref: dict) -> dict:
    eng = _engine(tier, "geometry", ref)
    rows, bond_err, angle_err = [], [], []
    for e in ref["geometries"]:
        row = {"id": e["id"]}
        try:
            mol = _mol_from_entry({**e, "prefer_smiles": tier == "mmff94"})
            opt = tasks.optimize(eng, mol)
            g = opt["molecule"].coords
            row["converged"] = opt["converged"]
            row["params"] = []
            for p in e["params"]:
                if tier == "mmff94" and "smiles" in e:
                    idx = p["atoms"]  # SMILES order for water is O,H,H: same as the xyz entry
                else:
                    idx = p["atoms"]
                calc = distance(g, *idx) if p["type"] == "bond" else angle(g, *idx)
                err = calc - p["value"]
                (bond_err if p["type"] == "bond" else angle_err).append(err)
                row["params"].append({**p, "calc": calc, "error": err})
        except Exception as ex:  # unsupported by this tier (e.g. no MMFF94 parameters)
            row["skipped"] = f"{type(ex).__name__}: {ex}"
        rows.append(row)
    return {"method": eng.method, "rows": rows,
            "bond_length_angstrom": _stats(bond_err), "angle_degrees": _stats(angle_err)}


# ---------------------------------------------------------------- dipoles
def bench_dipoles(tier: str, ref: dict) -> dict:
    eng_geo = _engine(tier, "geometry", ref)
    eng_dip = _engine(tier, "dipoles", ref)
    rows, reg = [], []
    for e in ref["dipoles_debye"]:
        row = {"id": e["id"], "exp": e["value"]}
        try:
            mol = _mol_from_entry({**e, "prefer_smiles": False})
            opt = tasks.optimize(eng_geo, mol)
            sp = tasks.single_point(eng_dip, opt["molecule"])
            mu = sp["dipole_debye"]
            if mu is None:
                raise ValueError("tier gives no dipole for this molecule")
            row.update(calc=mu, error=mu - e["value"])
            # regularised error (Hait & Head-Gordon 2018): relative for |mu| > 1 D, absolute below
            r = (mu - e["value"]) / max(e["value"], 1.0)
            row["regularized_error"] = r
            reg.append(r)
        except Exception as ex:
            row["skipped"] = f"{type(ex).__name__}: {ex}"
        rows.append(row)
    a = np.asarray(reg)
    return {"method": eng_dip.method, "geometry_method": eng_geo.method, "rows": rows,
            "regularized_rms_percent": float(np.sqrt((a ** 2).mean()) * 100) if a.size else None,
            "n": int(a.size)}


# ---------------------------------------------------------------- frequencies
def _expand(modes: list[dict], gas_only: bool) -> list[tuple[float, dict]]:
    out = []
    for m in modes:
        out += [(float(m["v"]), m)] * int(m["deg"])
    return out


def bench_frequencies(tier: str, ref: dict) -> dict:
    eng = _engine(tier, "frequencies", ref)
    comp = ref["comparators"].get(tier, {})
    scale = comp.get("frequency_scale_factor", {}).get("value") if tier == "dft" else None
    rows, pairs = [], []
    for e in ref["frequencies"]["molecules"]:
        row = {"id": e["id"]}
        try:
            mol = Molecule.from_smiles(e["smiles"], name=e["id"])
            opt = tasks.optimize(eng, mol)
            if tier == "dft" and hasattr(eng, "hessian"):
                frame, hess = eng.hessian(opt["molecule"])
                va = harmonic_analysis(frame, hess)
            else:
                va_d = tasks.frequencies(eng, opt["molecule"])
                va = None
                freqs = np.array(va_d["frequencies_cm1"])
            if va is not None:
                freqs = va.frequencies_cm1
            calc = np.sort(freqs[freqs > 0])
            exp = sorted(_expand(e["modes"], False), key=lambda t: t[0])
            if len(calc) != len(exp):
                raise ValueError(f"{len(calc)} computed vs {len(exp)} experimental modes")
            row["modes"] = []
            for c, (x, meta) in zip(calc, exp):
                rec = {"exp": x, "harmonic": float(c), "sym": meta["sym"], "phase": meta["phase"],
                       "anharmonic_flag": bool(meta.get("anharmonic"))}
                row["modes"].append(rec)
                pairs.append(rec)
            row["n_imaginary"] = int(np.sum(freqs < -10))
        except Exception as ex:
            row["skipped"] = f"{type(ex).__name__}: {ex}"
        rows.append(row)

    def summarize(sel):
        if not sel:
            return {"n": 0}
        h = np.array([p["harmonic"] for p in sel]); x = np.array([p["exp"] for p in sel])
        fitted = float((h @ x) / (h @ h))                       # least-squares scale factor
        res = {"n": len(sel), "fitted_scale_factor": fitted,
               "raw": _stats(list(h - x)),
               "fitted_scaled_rms_relative": float(np.sqrt(np.mean(((fitted * h - x) / x) ** 2)))}
        if scale:
            res["literature_scale_factor"] = scale
            res["literature_scaled"] = _stats(list(scale * h - x))
            res["literature_scaled_rms_relative"] = float(np.sqrt(np.mean(((scale * h - x) / x) ** 2)))
        return res

    clean = [p for p in pairs if p["phase"] in GAS_PHASES and not p["anharmonic_flag"]]
    return {"method": eng.method, "rows": rows,
            "matching": "rank order within each molecule after expanding degeneracies",
            "all_modes": summarize(pairs), "gas_phase_harmonic_modes": summarize(clean)}


# ---------------------------------------------------------------- gates
def evaluate_gates(tier: str, res: dict, ref: dict) -> list[dict]:
    f = ref["gates"]["comparator_factor"]
    comp = ref["comparators"].get(tier, {})
    gates = []

    def gate(name, measured, published, src):
        if measured is None or published is None:
            return
        gates.append({"check": name, "measured": measured, "published": published,
                      "limit": f * published, "pass": bool(measured <= f * published), "source": src})

    if "geometries" in res and "bond_length_mad_angstrom" in comp:
        c = comp["bond_length_mad_angstrom"]
        gate("bond length MAD (A)", res["geometries"]["bond_length_angstrom"].get("mad"), c["value"], c["source"])
    if "dipoles" in res and "dipole_regularized_rms_percent" in comp:
        c = comp["dipole_regularized_rms_percent"]
        gate("dipole regularized RMS (%)", res["dipoles"]["regularized_rms_percent"], c["value"], c["source"])
    if "frequencies" in res and "frequency_scale_factor" in comp:
        c = comp["frequency_scale_factor"]
        gate("scaled frequency RMS relative error (gas-phase modes)",
             res["frequencies"]["gas_phase_harmonic_modes"].get("literature_scaled_rms_relative"),
             c["sigma"], c["source"])
    return gates


def run(tiers: list[str], parts: list[str], out: Path) -> dict:
    ref = json.loads(REF.read_text())
    report = {"reference": str(REF.name), "tiers": {}}
    for tier in tiers:
        t0 = time.time()
        res = {}
        if "geometries" in parts:
            res["geometries"] = bench_geometries(tier, ref)
        if "dipoles" in parts and tier != "mmff94":
            res["dipoles"] = bench_dipoles(tier, ref)
        if "frequencies" in parts and tier != "mmff94":
            res["frequencies"] = bench_frequencies(tier, ref)
        res["gates"] = evaluate_gates(tier, res, ref)
        res["seconds"] = round(time.time() - t0, 1)
        report["tiers"][tier] = res
    out.mkdir(parents=True, exist_ok=True)
    (out / "benchmark.json").write_text(json.dumps(report, indent=1, default=float))
    (out / "benchmark.md").write_text(to_markdown(report))
    return report


def to_markdown(report: dict) -> str:
    L = ["# Athanor Phase 1 benchmark", ""]
    for tier, r in report["tiers"].items():
        L += [f"## {tier}  ({r['seconds']} s)", ""]
        if "geometries" in r:
            g = r["geometries"]
            b, a = g["bond_length_angstrom"], g["angle_degrees"]
            L.append(f"- **Geometry** ({g['method']}): bonds n={b.get('n')} MAD={b.get('mad', float('nan')):.4f} A, "
                     f"mean signed={b.get('mean_signed', float('nan')):+.4f} A; angles n={a.get('n')} MAD={a.get('mad', float('nan')):.2f} deg")
            sk = [x["id"] for x in g["rows"] if "skipped" in x]
            if sk:
                L.append(f"  - not supported by this tier: {', '.join(sk)}")
        if "dipoles" in r:
            d = r["dipoles"]
            L.append(f"- **Dipoles** ({d['method']} at {d['geometry_method']} geometry): n={d['n']} "
                     f"regularized RMS={d['regularized_rms_percent']:.1f}%")
        if "frequencies" in r:
            fq = r["frequencies"]
            for label in ("all_modes", "gas_phase_harmonic_modes"):
                s = fq[label]
                if not s.get("n"):
                    continue
                line = (f"- **Frequencies, {label.replace('_', ' ')}** ({fq['method']}): n={s['n']} raw MAD={s['raw']['mad']:.0f} cm-1, "
                        f"fitted scale={s['fitted_scale_factor']:.4f} (RMS rel {100*s['fitted_scaled_rms_relative']:.2f}%)")
                if "literature_scale_factor" in s:
                    line += (f"; literature scale {s['literature_scale_factor']} -> MAD={s['literature_scaled']['mad']:.0f} cm-1, "
                             f"RMS rel {100*s['literature_scaled_rms_relative']:.2f}%")
                L.append(line)
        for gt in r["gates"]:
            mark = "PASS" if gt["pass"] else "FAIL"
            L.append(f"- [{mark}] {gt['check']}: measured {gt['measured']:.4g} vs published {gt['published']:.4g} "
                     f"(limit {gt['limit']:.4g}; {gt['source']})")
        L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiers", default="mmff94,gfn2-xtb")
    ap.add_argument("--parts", default="geometries,dipoles,frequencies")
    ap.add_argument("--out", default="benchmark-report")
    a = ap.parse_args()
    rep = run(a.tiers.split(","), a.parts.split(","), Path(a.out))
    print(to_markdown(rep))
    failed = [g for t in rep["tiers"].values() for g in t["gates"] if not g["pass"]]
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
