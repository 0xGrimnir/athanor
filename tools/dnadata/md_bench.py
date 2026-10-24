"""DNA dynamics benchmarks: Athanor's OpenMM set-up against AMBER, and a simulated dodecamer against
its solution structure.

1. Same energies as AMBER. The Dickerson-Drew dodecamer (CGCGAATTCGCG) is built and relaxed by
   Athanor; the same coordinates are parameterised by AmberTools' tleap (leaprc.DNA.OL15, mbondi3
   radii for GBn2) and evaluated by OpenMM from that prmtop and by AMBER's own sander (igb=8, 150 mM
   salt). Gate: Athanor's force-field set-up gives the same energy as both.
2. A dodecamer in motion. Langevin dynamics as the app runs it (GBn2, 150 mM, 300 K, 4 fs with hydrogen
   mass repartitioning). After equilibration, the central eight base pairs are compared with the NMR
   solution structure of the same sequence refined against residual dipolar couplings (PDB 1NAJ; Tjandra
   et al., J. Am. Chem. Soc. 122, 6190, 2000) and with the 1.4 A crystal structure (PDB 355D; Shui et al.,
   Biochemistry 37, 8341, 1998). Gates: the duplex stays B-DNA close to the solution structure.

Usage: python md_bench.py WORK_DIR [production_ps]
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
import openmm as mm  # noqa: E402
import openmm.app as app  # noqa: E402
import openmm.unit as u  # noqa: E402

from athanor_core.dna import md  # noqa: E402
from athanor_core.dna.api import _residues  # noqa: E402
from athanor_core.dna.builder import build  # noqa: E402
from athanor_core.dna.geometry import kabsch  # noqa: E402
from athanor_core.dna.structure import analyse, read_pdb  # noqa: E402

WORK = Path(sys.argv[1])
WORK.mkdir(parents=True, exist_ok=True)
PROD = float(sys.argv[2]) if len(sys.argv) > 2 else 500.0
SEQ = "CGCGAATTCGCG"
results, summary = [], {}


def gate(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ------------------------------------------------------------------ 1. energies against AMBER
r = md.relax(build(SEQ, "B")["atoms"])
# hydrogens are re-added by Modeller in the same order as in relax(); no constraints, so every bond term is
# there to compare (the app constrains bonds to hydrogen during dynamics)
s = md.System(r["atoms"], constraints=None)
ctx = mm.Context(s.system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
# put the relaxed hydrogens back exactly: rebuild the full coordinate set in the System's atom order
full = {(c, int(n), name): xyz for c, n, _, name, xyz in r["atoms"] + r["hydrogens"]}
pos = []
for (c, n, name), p0 in zip(s.keys, s.positions.value_in_unit(u.angstrom)):
    pos.append(full.get((c, n, name), np.array(p0)))
pos = np.array(pos)
ctx.setPositions(pos * 0.1)


def energy(context):
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)


e_ath = energy(ctx)
pdb = WORK / "ddd.pdb"
with open(pdb, "w") as fh:
    app.PDBFile.writeFile(s.topology, pos * 0.1 * u.nanometer, fh, keepIds=True)
# AMBER's terminal residue names: 5' residues D?5, 3' residues D?3
lines, res_order = [], []
for line in pdb.read_text().splitlines():
    if line.startswith(("ATOM", "HETATM")):
        key = (line[21], int(line[22:26]))
        if key not in res_order:
            res_order.append(key)
first = {}
last = {}
for c, n in res_order:
    first.setdefault(c, n)
    last[c] = n
for line in pdb.read_text().splitlines():
    if line.startswith(("ATOM", "HETATM")):
        c, n = line[21], int(line[22:26])
        rn = line[17:20].strip()
        if n == first[c]:
            rn = rn + "5"
        elif n == last[c]:
            rn = rn + "3"
        line = line[:17] + f"{rn:>3s}" + line[20:]
    lines.append(line)
(WORK / "ddd_amber.pdb").write_text("\n".join(lines) + "\n")
(WORK / "leap.in").write_text("source leaprc.DNA.OL15\nset default PBRadii mbondi3\nm = loadpdb ddd_amber.pdb\n"
                              "saveamberparm m ddd.prmtop ddd.inpcrd\nquit\n")
lp = subprocess.run(["tleap", "-f", "leap.in"], cwd=WORK, capture_output=True, text=True)
print(lp.stdout[-1500:])
prm = app.AmberPrmtopFile(str(WORK / "ddd.prmtop"))
crd = app.AmberInpcrdFile(str(WORK / "ddd.inpcrd"))
sp = prm.createSystem(nonbondedMethod=app.NoCutoff, implicitSolvent=app.GBn2, implicitSolventSaltConc=0.15 * u.molar,
                      temperature=300 * u.kelvin, sasaMethod=None)       # no surface-area term, as sander's gbsa=0
c2 = mm.Context(sp, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
c2.setPositions(crd.positions)
e_prm = energy(c2)
(WORK / "sp.in").write_text("single point\n &cntrl\n  imin=1, maxcyc=1, ncyc=1, ntb=0, igb=8, saltcon=0.15, cut=999.0, rgbmax=999.0, gbsa=0, ntpr=1,\n /\n")
sd = subprocess.run(["sander", "-O", "-i", "sp.in", "-p", "ddd.prmtop", "-c", "ddd.inpcrd", "-o", "sp.out"], cwd=WORK,
                    capture_output=True, text=True)
out = (WORK / "sp.out").read_text() if (WORK / "sp.out").exists() else sd.stderr
terms = dict(re.findall(r"(BOND|ANGLE|DIHED|VDWAALS|EEL|EGB|1-4 VDW|1-4 EEL)\s*=\s*(-?[\d.]+)", out))
terms = {k: float(v) for k, v in terms.items()}
print("sander terms:", terms)


def by_term(system, positions):
    """Energy of each part of an OpenMM system, in sander's groupings (kcal/mol)."""
    kinds = {"HarmonicBondForce": "bonds", "HarmonicAngleForce": "angles", "PeriodicTorsionForce": "torsions",
             "NonbondedForce": "coulomb+vdW", "CustomGBForce": "GB solvation", "GBSAOBCForce": "GB solvation"}
    out = {}
    for k, f in enumerate(system.getForces()):
        f.setForceGroup(k)
    c = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    c.setPositions(positions)
    for k, f in enumerate(system.getForces()):
        name = kinds.get(type(f).__name__)
        if name:
            e = c.getState(getEnergy=True, groups={k}).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
            out[name] = out.get(name, 0.0) + e
    return out


# evaluate Athanor's system at exactly the coordinates tleap read (rounded to the PDB file's 0.001 A), matching
# atoms by residue and name since tleap orders atoms by its own templates
crd_nm = np.asarray(crd.positions.value_in_unit(u.nanometer))
amber_xyz = {(a.residue.index, a.name): crd_nm[a.index] for a in prm.topology.atoms()}
same_pos = []
for a in s.topology.atoms():
    x = amber_xyz.get((a.residue.index, a.name))
    same_pos.append(x if x is not None else pos[a.index] * 0.1)
missing = sum(1 for a in s.topology.atoms() if (a.residue.index, a.name) not in amber_xyz)
print(f"atoms matched to tleap's: {s.topology.getNumAtoms() - missing} of {s.topology.getNumAtoms()}")
ours = by_term(s.system, np.array(same_pos))
prm_terms = by_term(sp, crd.positions)
sander = {"bonds": terms.get("BOND"), "angles": terms.get("ANGLE"), "torsions": terms.get("DIHED"),
          "coulomb+vdW": sum(terms.get(k, 0.0) for k in ("VDWAALS", "EEL", "1-4 VDW", "1-4 EEL")), "GB solvation": terms.get("EGB")}
print(f"{'term':14s} {'Athanor':>12s} {'prmtop':>12s} {'sander':>12s}")
for t in sander:
    print(f"{t:14s} {ours.get(t, float('nan')):12.4f} {prm_terms.get(t, float('nan')):12.4f} {sander[t]:12.4f}")
summary["energy"] = {"athanor": round(e_ath, 2), "openmm_prmtop": round(e_prm, 2), "sander": round(sum(sander.values()), 2),
                     "terms": {t: {"athanor": round(ours.get(t, 0), 3), "sander": sander[t]} for t in sander}}
# the same force field from two routes (OpenMM's amber14 files, tleap's prmtop) must give the same numbers
# term by term; against sander, a different program, the terms may differ by its constants (AMBER's
# Coulomb constant 332.0522 vs OpenMM's 332.0637) and its own GBn2 code, so a relative tolerance
same = all(abs(ours[t] - prm_terms[t]) <= 0.01 for t in sander)
rel = {t: abs(ours[t] - sander[t]) / max(abs(sander[t]), 1.0) for t in sander}
gate("Athanor's OpenMM set-up gives AMBER's energies for the same coordinates (OL15 + GBn2, 150 mM): every term "
     "equal to the same force field read from tleap's prmtop (0.01 kcal/mol), and within 0.05% of sander's",
     same and max(rel.values()) <= 5e-4,
     f"totals: Athanor {e_ath:.2f}, prmtop {e_prm:.2f}, sander {sum(sander.values()):.2f} kcal/mol; largest term difference "
     f"from prmtop {max(abs(ours[t] - prm_terms[t]) for t in sander):.4f} kcal/mol, from sander {100 * max(rel.values()):.3f}%")
if PROD <= 0:
    (WORK / "md_report.json").write_text(json.dumps({"results": results, "summary": summary}, indent=1))
    sys.exit(0 if all(x["pass"] for x in results) else 1)

# ------------------------------------------------------------------ 2. dynamics against the solution structure
def fetch(pid):
    f = WORK / f"{pid}.pdb"
    if not f.exists():
        with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=120) as resp:
            f.write_bytes(resp.read())
    return [r for r in read_pdb(f.read_text())]


def heavy(res):
    for r in res:
        r.atoms = {k: v for k, v in r.atoms.items() if not k.lstrip("0123456789").startswith("H")}
    return res


refs = {pid: heavy(fetch(pid)) for pid in ("1NAJ", "355D")}
ref_an = {pid: analyse(res) for pid, res in refs.items()}
for pid, a in ref_an.items():
    print(pid, a["sequence"], len(a["pairs"]), "pairs")
CENTRAL = range(2, 9)        # steps 3-9: the central eight base pairs


def central_steps(an, n=12):
    """Steps 3-9 by position along strand I (NaN where an end has frayed and a pair is missing)."""
    out = np.full((n - 1, 6), np.nan)
    pi = [p["i"] for p in an["pairs"]]
    k = 0
    for a, b in zip(pi, pi[1:]):
        if b == a + 1:
            out[a] = an["steps"][k]["params"]
            k += 1
    return out[list(CENTRAL)]


def core_atoms(atoms_by_key, res_list, n=12):
    """Matched heavy atoms of base pairs 3-10 (both strands) between our model and a reference."""
    by = {}
    for r in res_list:
        by.setdefault(r.chain, []).append(r)
    ch = list(by)
    s1, s2 = by[ch[0]], by[ch[1]]
    P, Q = [], []
    for k in range(2, 10):
        for (c, num), cr in ((("A", k + 1), s1[k]), (("B", 2 * n - k), s2[n - 1 - k])):
            o = atoms_by_key.get((c, num), {})
            for nm, xyz in cr.atoms.items():
                if nm in o:
                    P.append(o[nm]); Q.append(xyz)
    return np.array(P), np.array(Q)


dyn = md.Dynamics(build(SEQ, "B")["atoms"], seed=11)
dyn.run(50.0)                                   # equilibration
frames = []
while dyn.time_ps < 50.0 + PROD - 1e-6:
    f = dyn.run(5.0)
    frames.append(f["atoms"])
    if len(frames) % 20 == 0:
        print(f"   {dyn.time_ps:.0f} ps, {f['temperature']:.0f} K, {f['ns_per_day']:.1f} ns/day", flush=True)
md_steps, rms = [], {pid: [] for pid in refs}
for atoms in frames:
    an = analyse(_residues({"atoms": atoms}))
    md_steps.append(central_steps(an))
    by = {}
    for c, n, _, name, xyz in atoms:
        by.setdefault((c, n), {})[name] = xyz
    for pid, res in refs.items():
        P, Q = core_atoms(by, res)
        rms[pid].append(kabsch(P, Q)[2])
md_mean = np.nanmean(np.nanmean(md_steps, axis=0), axis=0)
summary["md"] = {"production_ps": PROD, "frames": len(frames), "mean_twist": round(float(md_mean[5]), 2),
                 "mean_rise": round(float(md_mean[2]), 3), "mean_roll": round(float(md_mean[4]), 2),
                 "mean_slide": round(float(md_mean[1]), 3)}
for pid in refs:
    ref_mean = np.nanmean(central_steps(ref_an[pid]), axis=0)
    summary[pid] = {"mean_twist": round(float(ref_mean[5]), 2), "mean_rise": round(float(ref_mean[2]), 3),
                    "mean_roll": round(float(ref_mean[4]), 2), "mean_slide": round(float(ref_mean[1]), 3),
                    "rmsd_median": round(float(np.median(rms[pid])), 2), "rmsd_p90": round(float(np.percentile(rms[pid], 90)), 2)}
    print(pid, summary[pid])
print("MD", summary["md"])
nmr = summary["1NAJ"]
gate("in dynamics the dodecamer stays close to its NMR solution structure (central 8 bp, heavy atoms; median "
     "RMSD over the production run <= 2.5 A)", nmr["rmsd_median"] <= 2.5,
     f"median {nmr['rmsd_median']} A, 90th percentile {nmr['rmsd_p90']} A over {len(frames)} frames "
     f"({PROD:.0f} ps); to the crystal (355D): median {summary['355D']['rmsd_median']} A")
dtw = abs(summary["md"]["mean_twist"] - nmr["mean_twist"])
dri = abs(summary["md"]["mean_rise"] - nmr["mean_rise"])
gate("average helix twist and rise of the central 8 bp match the NMR structure (within 2.5 deg and 0.2 A)",
     dtw <= 2.5 and dri <= 0.2,
     f"twist {summary['md']['mean_twist']} vs {nmr['mean_twist']} deg; rise {summary['md']['mean_rise']} vs "
     f"{nmr['mean_rise']} A (crystal 355D: {summary['355D']['mean_twist']} deg, {summary['355D']['mean_rise']} A)")
(WORK / "md_report.json").write_text(json.dumps({"results": results, "summary": summary}, indent=1))
sys.exit(0 if all(x["pass"] for x in results) else 1)
