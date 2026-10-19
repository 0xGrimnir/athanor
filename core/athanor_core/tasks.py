"""High-level tasks. Every function returns plain data plus a provenance record."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from . import units
from .engines import Engine
from .molecule import Molecule
from .provenance import Provenance, hash_input
from .vibrations import finite_difference_hessian, harmonic_analysis

# Gaussian's default ("normal") convergence thresholds on forces, Eh/bohr.
# Reference: Gaussian 16 documentation, Opt keyword.
MAX_FORCE = 4.5e-4
RMS_FORCE = 3.0e-4


def _prov(engine: Engine, mol: Molecule, task: str, extra: dict | None = None) -> Provenance:
    params = dict(engine.parameters())
    if extra:
        params.update(extra)
    return Provenance(
        tier=engine.tier, engine=engine.engine_name, engine_version=engine.version(),
        method=engine.method, task=task, parameters=params,
        input_hash=hash_input(mol.symbols, mol.coords, mol.charge, mol.multiplicity, params),
    )


def single_point(engine: Engine, mol: Molecule) -> dict:
    prov = _prov(engine, mol, "energy")
    r = engine.compute(mol)
    out = {
        "energy_hartree": r.energy,
        "dipole_debye": None if r.dipole is None else float(np.linalg.norm(r.dipole) * units.AU_TO_DEBYE),
        "dipole_vector_debye": None if r.dipole is None else (np.asarray(r.dipole) * units.AU_TO_DEBYE).tolist(),
        "dipole_kind": engine.dipole_kind,
        "charges": None if r.charges is None else [float(q) for q in r.charges],
        "charge_model": r.charge_model,
        "max_force": float(np.abs(r.gradient).max()),
        "extras": r.extras,
    }
    out["provenance"] = prov.finish().to_dict()
    return out


def optimize(engine: Engine, mol: Molecule, max_steps: int = 500) -> dict:
    """Quasi-Newton (L-BFGS) minimisation on analytic gradients in Cartesian
    coordinates. Converged when max and RMS force meet Gaussian's defaults."""
    prov = _prov(engine, mol, "optimize", {"optimizer": "scipy L-BFGS-B",
                                           "max_force": MAX_FORCE, "rms_force": RMS_FORCE})
    trajectory: list[dict] = []
    cache: dict = {}

    def fg(x):
        key = x.tobytes()
        if key not in cache:
            m = mol.copy(x.reshape(-1, 3) * units.BOHR_TO_ANGSTROM)
            r = engine.compute(m)
            cache.clear()
            cache[key] = (r.energy, r.gradient.ravel().copy())
            trajectory.append({"energy": r.energy, "coords": m.coords.round(6).tolist(),
                               "max_force": float(np.abs(r.gradient).max())})
        return cache[key]

    x0 = (mol.coords * units.ANGSTROM_TO_BOHR).ravel()
    res = minimize(lambda x: fg(x)[0], x0, jac=lambda x: fg(x)[1], method="L-BFGS-B",
                   options={"maxiter": max_steps, "gtol": MAX_FORCE / 10, "ftol": 1e-14, "maxcor": 30})
    final = mol.copy(res.x.reshape(-1, 3) * units.BOHR_TO_ANGSTROM)
    rf = engine.compute(final)
    g = rf.gradient.ravel()
    max_f, rms_f = float(np.abs(g).max()), float(np.sqrt(np.mean(g ** 2)))
    converged = max_f < MAX_FORCE and rms_f < RMS_FORCE
    if not converged:
        prov.notes.append(f"not converged: max force {max_f:.2e}, rms {rms_f:.2e} Eh/bohr")
    return {
        "molecule": final,
        "energy_hartree": rf.energy,
        "converged": converged,
        "max_force": max_f,
        "rms_force": rms_f,
        "steps": len(trajectory),
        "trajectory": trajectory,
        "provenance": prov.finish().to_dict(),
    }


def frequencies(engine: Engine, mol: Molecule, step_bohr: float = 0.005) -> dict:
    """Harmonic frequencies and IR intensities at the given geometry. Run on an
    optimized structure: frequencies are only meaningful at a stationary point."""
    prov = _prov(engine, mol, "frequencies", {"hessian": "central finite differences of analytic gradients",
                                              "step_bohr": step_bohr})
    r0 = engine.compute(mol)
    max_f = float(np.abs(r0.gradient).max())
    hess, dmu = finite_difference_hessian(engine.compute, mol, step_bohr)
    va = harmonic_analysis(mol, hess, dmu)
    if max_f > MAX_FORCE:
        va.notes.append(f"geometry is not stationary (max force {max_f:.2e} Eh/bohr); frequencies are unreliable")
    prov.notes.extend(va.notes)
    out = va.to_dict()
    out["provenance"] = prov.finish().to_dict()
    return out
