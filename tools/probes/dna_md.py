"""Probe: OpenMM with amber14 (DNA.OL15) and GBn2 implicit solvent on an Athanor-built helix."""
import io, time, sys
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
from athanor_core.dna.builder import build, to_pdb

print("openmm", mm.__version__, [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])
import os
d = os.path.join(os.path.dirname(app.__file__), "data")
print([f for f in os.listdir(os.path.join(d, "amber14")) if "DNA" in f or "RNA" in f], os.listdir(os.path.join(d, "implicit")))
for seq, form in (("CGCGAATTCGCG", "B"), ("CGCGCGCGCGCG", "Z")):
    m = build(seq, form)
    # AMBER terminal residues: 5' ends have no phosphate
    n = len(seq)
    atoms = [a for a in m["atoms"] if not ((a[0] == "A" and a[1] == 1 or a[0] == "B" and a[1] == n + 1) and a[3] in ("P", "OP1", "OP2"))]
    pdb = app.PDBFile(io.StringIO(to_pdb({"atoms": atoms})))
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    mod = app.Modeller(pdb.topology, pdb.positions)
    try:
        mod.addHydrogens(ff)
    except Exception as e:
        print("addHydrogens failed:", e); continue
    print(form, "atoms with H", mod.topology.getNumAtoms(), [r.name for r in list(mod.topology.residues())[:3]], [r.name for r in list(mod.topology.residues())[-3:]])
    sysm = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds, implicitSolventSaltConc=0.15 * u.molar)
    integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 1 / u.picosecond, 0.002 * u.picoseconds)
    sim = app.Simulation(mod.topology, sysm, integ, mm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(mod.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=500)
    e1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.angstrom)
    names = [(a.residue.chain.index, a.residue.index, a.name) for a in mod.topology.atoms()]
    idx = {k: i for i, k in enumerate(names)}
    o3p = []
    res = list(mod.topology.residues())
    for r1, r2 in zip(res, res[1:]):
        if r1.chain.index != r2.chain.index: continue
        o3 = idx.get((r1.chain.index, r1.index, "O3'")); p = idx.get((r2.chain.index, r2.index, "P"))
        if o3 is not None and p is not None: o3p.append(np.linalg.norm(pos[o3] - pos[p]))
    print(form, "energy", e0, "->", e1, "O3'-P after minimise", round(min(o3p), 3), round(max(o3p), 3))
    t = time.time(); sim.step(2000); dt = time.time() - t
    print(form, f"2000 steps (4 ps) in {dt:.1f} s -> {4e-3 / dt * 86400:.1f} ns/day on CPU")
