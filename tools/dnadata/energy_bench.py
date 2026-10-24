"""Interaction-energy benchmarks for DNA Structure's pairing and stacking energies.

1. S22 (Jurecka et al., Phys. Chem. Chem. Phys. 8, 1985, 2006): 22 noncovalent complexes - hydrogen-
   bonded, dispersion-bound and mixed - including the adenine-thymine Watson-Crick pair (S22-7) and
   the stacked adenine-thymine pair (S22-15). Reference interaction energies are CCSD(T) at the
   complete-basis-set limit, as revised by Marshall et al. (J. Chem. Phys. 135, 194102, 2011; S22B);
   geometries and reference values are read from Psi4's S22 database file. Athanor's two methods are
   run on the same geometries: GFN2-xTB, and B3LYP-D3(BJ)/def2-TZVP with counterpoise correction.
2. In the helix: every base pair of a relaxed built duplex with GFN2-xTB; G-C pairs must come out more
   strongly bound than A-T pairs (three hydrogen bonds against two), and a G-C and an A-T pair are also
   computed with DFT to show how far the fast method is from the slow one.

Usage: python energy_bench.py WORK_DIR [OUT_JSON] [--no-dft]
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.dna import energy, md  # noqa: E402
from athanor_core.dna.builder import build  # noqa: E402

WORK = Path(sys.argv[1])
WORK.mkdir(parents=True, exist_ok=True)
DFT = "--no-dft" not in sys.argv
results, summary = [], {}


def gate(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ------------------------------------------------------------------ S22 from Psi4's database file
db = Path(sys.prefix) / "share" / "psi4" / "databases" / "S22.py"
text = db.read_text()
ref = {int(k): float(v) for k, v in re.findall(r"BIND_S22B\[\s*'%s-%s'\s*%\s*\(\s*dbse\s*,\s*'?(\d+)'?\s*\)\s*\]\s*=\s*(-?[\d.]+)", text)}
names = dict(re.findall(r"TAGL\[\s*'%s-%s'\s*%\s*\(\s*dbse\s*,\s*'?(\d+)'?\s*\)\s*\]\s*=\s*'([^']*)'", text))
geos = {}
for num, body in re.findall(r"GEOS\[\s*'%s-%s-%s'\s*%\s*\(\s*dbse\s*,\s*'?(\d+)'?\s*,\s*'dimer'\s*\)\s*\]\s*=\s*qcdb\.Molecule\(\s*\"\"\"(.*?)\"\"\"", text, re.S):
    frags, cur = [], []
    for line in body.strip().splitlines():
        line = line.strip()
        if line == "--":
            frags.append(cur)
            cur = []
            continue
        p = line.split()
        if len(p) == 4 and re.fullmatch(r"[A-Za-z]{1,2}", p[0]):
            cur.append((p[0].capitalize(), [float(x) for x in p[1:]]))
    frags.append(cur)
    geos[int(num)] = [([a for a, _ in f], np.array([x for _, x in f])) for f in frags if f]
print(f"S22 from {db}: {len(geos)} geometries, {len(ref)} S22B reference energies")
if not (len(geos) == 22 and len(ref) == 22):
    i = text.find("GEOS[")
    print(text[i:i + 1500])
    i = text.find("BIND_S22B[")
    print(text[i:i + 400])
    sys.exit(2)

rows = {}
for k in sorted(geos):
    a, b = geos[k]
    t = time.time()
    g = energy.interaction_gfn2([a], [b])
    row = {"name": names.get(str(k), str(k)), "ref": ref[k], "gfn2": g, "atoms": len(a[0]) + len(b[0])}
    if DFT:
        row["dft"] = energy.interaction_dft([a], [b])
    rows[k] = row
    print(f"   S22-{k:2d} {row['name'][:38]:38s} ref {ref[k]:7.3f}  GFN2-xTB {g:7.3f}" +
          (f"  DFT {row['dft']:7.3f}" if DFT else "") + f"  ({time.time() - t:.0f} s)", flush=True)
summary["s22"] = rows


def mad(key):
    return float(np.mean([abs(r[key] - r["ref"]) for r in rows.values()]))


summary["mad"] = {"gfn2": round(mad("gfn2"), 3)}
gate("GFN2-xTB on S22 against CCSD(T)/CBS (S22B): mean absolute error <= 1.5 kcal/mol",
     mad("gfn2") <= 1.5, f"MAD {mad('gfn2'):.2f} kcal/mol; adenine-thymine Watson-Crick {rows[7]['gfn2']:.2f} vs {rows[7]['ref']:.2f}, "
     f"stacked {rows[15]['gfn2']:.2f} vs {rows[15]['ref']:.2f}")
if DFT:
    summary["mad"]["dft"] = round(mad("dft"), 3)
    worst_at = max(abs(rows[k]["dft"] - rows[k]["ref"]) for k in (7, 15))
    gate("B3LYP-D3(BJ)/def2-TZVP (counterpoise) on S22 against CCSD(T)/CBS: mean absolute error <= 0.6 kcal/mol, "
         "and both adenine-thymine complexes within 1.0 kcal/mol", mad("dft") <= 0.6 and worst_at <= 1.0,
         f"MAD {mad('dft'):.2f} kcal/mol; adenine-thymine Watson-Crick {rows[7]['dft']:.2f} vs {rows[7]['ref']:.2f}, "
         f"stacked {rows[15]['dft']:.2f} vs {rows[15]['ref']:.2f}")

# ------------------------------------------------------------------ in the helix
SEQ = "CGCGAATTCGCG"
r = md.relax(build(SEQ, "B")["atoms"])
h = energy.helix_energies(r["atoms"], r["hydrogens"], len(SEQ))
gc = [e for b, e in zip(SEQ, h["pairs"]) if b in "GC"]
at = [e for b, e in zip(SEQ, h["pairs"]) if b in "AT"]
summary["helix"] = {"sequence": SEQ, "pairs_gfn2": [round(x, 2) for x in h["pairs"]], "steps_gfn2": [round(x, 2) for x in h["steps"]]}
print("helix pairs (GFN2-xTB):", summary["helix"]["pairs_gfn2"])
print("helix stacking (GFN2-xTB):", summary["helix"]["steps_gfn2"])
gate("in a built helix, G-C pairs are more strongly bound than A-T pairs (GFN2-xTB; mean difference >= 5 kcal/mol)",
     np.mean(gc) - np.mean(at) <= -5, f"G-C {np.mean(gc):.1f}, A-T {np.mean(at):.1f} kcal/mol")
if DFT:
    fr = energy.fragments(r["atoms"], r["hydrogens"], len(SEQ))
    d_gc = energy.interaction_dft([fr[0][0]], [fr[0][1]])
    d_at = energy.interaction_dft([fr[5][0]], [fr[5][1]])
    summary["helix"]["dft_pairs"] = {"C1-G24": round(d_gc, 2), "T6-A19": round(d_at, 2)}
    print(f"helix DFT: C-G {d_gc:.2f} (GFN2 {h['pairs'][0]:.2f}), A-T {d_at:.2f} (GFN2 {h['pairs'][5]:.2f})")
(WORK / "energy_report.json").write_text(json.dumps({"results": results, "summary": summary}, indent=1, default=float))
# what the app shows (two decimals, so it regenerates identically)
out = {"s22_mad_kcal": {k: round(v, 2) for k, v in summary["mad"].items()},
       "s22_at": {m: {"watson_crick": round(rows[7][m], 2), "stacked": round(rows[15][m], 2)} for m in ("gfn2", "dft") if m in rows[7]},
       "s22_at_ref": {"watson_crick": rows[7]["ref"], "stacked": rows[15]["ref"]},
       "passed": all(x["pass"] for x in results)}
if len(sys.argv) > 2 and not sys.argv[2].startswith("--"):
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
sys.exit(0 if all(x["pass"] for x in results) else 1)
