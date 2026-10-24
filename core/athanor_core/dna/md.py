"""Relaxation and molecular dynamics of a DNA duplex with OpenMM.

Force field: AMBER OL15 for DNA (Zgarbova et al., J. Chem. Theory Comput. 11, 5723, 2015; amber14 set
in OpenMM), which was refined for both B- and Z-DNA. Solvent: the GBn2 generalized-Born implicit model
(Nguyen, Roe & Simmerling, J. Chem. Theory Comput. 9, 2020, 2013) with a Debye-Hueckel salt screening
term. Implicit solvent is fast enough to run on a desktop without a GPU, at the cost of no explicit
water or counterions; it is a model of DNA in solution, not of the crystal.

* :func:`relax` closes the backbone of a built helix: hydrogens are added, the bases are held by
  harmonic restraints to where the builder put them, and the energy is minimised. The O3'-P links
  then take their force-field length (1.61 A) instead of whatever the rigid nucleotide templates gave.
* :class:`Dynamics` runs Langevin dynamics in chunks, so the app can show frames while it runs and
  stop at any time. It uses the fastest OpenMM platform available (CUDA, OpenCL, then CPU).
"""
from __future__ import annotations

import io
import math
import time

import numpy as np

KCAL = 4.184  # kJ per kcal

try:  # OpenMM ships with the app's runtime; the analysis and builder work without it
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u
except ImportError:  # pragma: no cover
    mm = app = u = None

FF_FILES = ("amber14-all.xml", "implicit/gbn2.xml")
BASE_ATOMS = {"N9", "C8", "N7", "C5", "C6", "O6", "N6", "N1", "C2", "N2", "N3", "C4", "O2", "N4", "O4", "C7", "C1'"}


def available() -> bool:
    return mm is not None


def version() -> str | None:
    return mm.__version__ if mm else None


def kappa(salt_molar: float, temperature_k: float = 300.0, dielectric: float = 78.5) -> float:
    """Debye-Hueckel screening parameter in 1/nm (as OpenMM's AMBER loader computes it)."""
    return 50.33355 * math.sqrt(salt_molar / dielectric / temperature_k) * 10.0


def _pdb(atoms) -> str:
    lines, k, last = [], 0, None
    for chain, num, res, name, xyz in atoms:
        if last is not None and chain != last:
            lines.append("TER")
        last = chain
        k += 1
        nm = f" {name:<3s}" if len(name) < 4 else name
        el = name.lstrip("0123456789")[0]
        lines.append(f"ATOM  {k:5d} {nm} {res:>3s} {chain}{num:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {el:>2s}")
    lines += ["TER", "END"]
    return "\n".join(lines) + "\n"


def _five_prime(atoms):
    """The first residue of each chain, in the order the atoms are listed (5' to 3')."""
    first = {}
    for c, n, *_ in atoms:
        first.setdefault(c, n)
    return first


def strip_5p_phosphate(atoms):
    """Synthetic and crystallised oligonucleotides have a 5'-OH: no phosphate on the first residue."""
    first = _five_prime(atoms)
    return [a for a in atoms if not (first[a[0]] == a[1] and a[3] in ("P", "OP1", "OP2"))]


class System:
    """An OpenMM system for a duplex given as Athanor atoms [(chain, num, resname, name, xyz)]."""

    def __init__(self, atoms, salt: float = 0.15, temperature: float = 300.0, hmr: bool = False):
        if mm is None:
            raise RuntimeError("OpenMM is not installed in this runtime")
        atoms = strip_5p_phosphate(atoms)
        pdb = app.PDBFile(io.StringIO(_pdb(atoms)))
        self.ff = app.ForceField(*FF_FILES)
        mod = app.Modeller(pdb.topology, pdb.positions)
        mod.addHydrogens(self.ff)
        self.topology, self.positions = mod.topology, mod.positions
        self.salt, self.temperature = salt, temperature
        # hydrogen mass repartitioning (1.5 amu hydrogens) allows a 4 fs time step (Hopkins et al., JCTC 11, 1864, 2015)
        extra = {"hydrogenMass": 1.5 * u.amu} if hmr else {}
        self.system = self.ff.createSystem(self.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds,
                                           implicitSolventKappa=kappa(salt, temperature) / u.nanometer, **extra)
        self.keys = [(a.residue.chain.id, int(a.residue.id), a.name) for a in self.topology.atoms()]
        self.heavy = [i for i, a in enumerate(self.topology.atoms()) if a.element is not None and a.element.symbol != "H"]
        self.order = atoms

    def restrain(self, indices, k_kcal_per_a2: float):
        f = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        f.addGlobalParameter("k", k_kcal_per_a2 * KCAL * 100.0)  # kcal/mol/A^2 -> kJ/mol/nm^2
        for p in ("x0", "y0", "z0"):
            f.addPerParticleParameter(p)
        pos = self.positions.value_in_unit(u.nanometer)
        for i in indices:
            f.addParticle(i, list(pos[i]))
        self.system.addForce(f)
        return f

    def atoms_from(self, positions_nm) -> list:
        """Heavy atoms back in Athanor's format and order (positions in nm -> A)."""
        idx = {k: i for i, k in enumerate(self.keys)}
        out = []
        for c, n, r, name, _ in self.order:
            i = idx.get((c, int(n), name))
            if i is not None:
                out.append((c, n, r, name, np.asarray(positions_nm[i]) * 10.0))
        return out

    def hydrogens_from(self, positions_nm) -> list:
        out = []
        for i, a in enumerate(self.topology.atoms()):
            if a.element is not None and a.element.symbol == "H":
                out.append((a.residue.chain.id, int(a.residue.id), "D" + a.residue.name[1:2], a.name, np.asarray(positions_nm[i]) * 10.0))
        return out


def platform(prefer: str | None = None):
    """The fastest OpenMM platform that works here: CUDA, OpenCL, then CPU."""
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    order = [prefer] if prefer else []
    order += ["CUDA", "OpenCL", "CPU", "Reference"]
    for n in order:
        if n in names:
            try:
                p = mm.Platform.getPlatformByName(n)
                # a platform can be listed and still fail on this machine's driver: try a tiny context
                s = mm.System()
                s.addParticle(1.0)
                mm.Context(s, mm.VerletIntegrator(0.001), p)
                return p
            except Exception:  # noqa: BLE001
                continue
    raise RuntimeError("no OpenMM platform works on this machine")


def _o3p(atoms) -> list[float]:
    by = {}
    for c, n, _, name, xyz in atoms:
        by.setdefault((c, n), {})[name] = xyz
    out = []
    keys = list(by)
    for (c1, n1), (c2, n2) in zip(keys, keys[1:]):
        if c1 == c2 and "O3'" in by[(c1, n1)] and "P" in by[(c2, n2)]:
            out.append(float(np.linalg.norm(by[(c1, n1)]["O3'"] - by[(c2, n2)]["P"])))
    return out


def relax(atoms, salt: float = 0.15, restraint: float = 50.0, tolerance: float = 1.0, max_iter: int = 0) -> dict:
    """Close the backbone: minimise with the bases restrained (kcal/mol/A^2) to their built positions.

    Returns heavy atoms in Athanor's format, hydrogens separately, the O3'-P distances, and how far the
    restrained base atoms moved."""
    s = System(atoms, salt)
    base = [i for i in s.heavy if s.keys[i][2] in BASE_ATOMS]
    s.restrain(base, restraint)
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(s.system, integ, platform())
    ctx.setPositions(s.positions)
    e0 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole) / KCAL
    t = time.time()
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance * KCAL * 10.0, max_iter)
    st = ctx.getState(getEnergy=True, getPositions=True)
    pos = st.getPositions(asNumpy=True).value_in_unit(u.nanometer)
    e1 = st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole) / KCAL
    p0 = np.asarray(s.positions.value_in_unit(u.nanometer))
    moved = np.linalg.norm(pos[base] - p0[base], axis=1) * 10.0
    heavy = s.atoms_from(pos)
    return {"atoms": heavy, "hydrogens": s.hydrogens_from(pos), "o3p": _o3p(heavy),
            "energy_kcal": [e0, e1], "base_shift_max": float(moved.max()), "base_shift_rms": float(np.sqrt((moved ** 2).mean())),
            "seconds": time.time() - t, "platform": ctx.getPlatform().getName(),
            "method": f"OpenMM {mm.__version__}: AMBER OL15 DNA force field, GBn2 implicit solvent, {salt * 1000:.0f} mM salt; "
                      f"energy minimised with the bases restrained ({restraint:g} kcal/mol/A^2) to the built structure"}


class Dynamics:
    """Langevin dynamics of a duplex in GBn2 implicit solvent, run in chunks.

    Starts from a relaxed structure: a short restrained minimisation, then free dynamics at the given
    temperature. Each :meth:`run` advances the simulation and returns the heavy-atom coordinates at
    the end of the chunk."""

    def __init__(self, atoms, salt: float = 0.15, temperature: float = 300.0, timestep_fs: float = 4.0,
                 friction_per_ps: float = 1.0, seed: int = 0, prefer: str | None = None):
        self.s = System(atoms, salt, temperature, hmr=timestep_fs > 2.5)
        integ = mm.LangevinMiddleIntegrator(temperature * u.kelvin, friction_per_ps / u.picosecond, timestep_fs * u.femtoseconds)
        if seed:
            integ.setRandomNumberSeed(seed)
        self.dt_ps = timestep_fs / 1000.0
        self.plat = platform(prefer)
        self.sim = app.Simulation(self.s.topology, self.s.system, integ, self.plat)
        self.sim.context.setPositions(self.s.positions)
        self.sim.minimizeEnergy(tolerance=10.0 * KCAL * 10.0 * u.kilojoule_per_mole / u.nanometer, maxIterations=2000)
        self.sim.context.setVelocitiesToTemperature(temperature * u.kelvin, seed or 1)
        self.time_ps = 0.0
        self.temperature = temperature

    def describe(self) -> str:
        return (f"OpenMM {mm.__version__} on {self.plat.getName()}: AMBER OL15 DNA force field, GBn2 implicit solvent "
                f"with {self.s.salt * 1000:.0f} mM salt, Langevin dynamics at {self.temperature:.0f} K, "
                f"{self.dt_ps * 1000:.0f} fs steps with hydrogen mass repartitioning")

    def run(self, ps: float) -> dict:
        n = max(1, int(round(ps / self.dt_ps)))
        t = time.time()
        self.sim.step(n)
        wall = time.time() - t
        self.time_ps += n * self.dt_ps
        st = self.sim.context.getState(getPositions=True, getEnergy=True)
        pos = st.getPositions(asNumpy=True).value_in_unit(u.nanometer)
        ke = st.getKineticEnergy().value_in_unit(u.kilojoule_per_mole)
        ndof = 3 * self.s.system.getNumParticles() - self.s.system.getNumConstraints()
        temp = 2 * ke / (ndof * 0.00831446261815324)
        return {"atoms": self.s.atoms_from(pos), "time_ps": self.time_ps, "temperature": temp,
                "potential_kcal": st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole) / KCAL,
                "ns_per_day": (n * self.dt_ps / 1000.0) / max(wall, 1e-9) * 86400.0, "platform": self.plat.getName()}
