"""High-level tasks. Every function returns plain data plus a provenance record."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from . import units
from .engines import Engine
from .molecule import Molecule
from .provenance import Provenance, hash_input
from .vibrations import finite_difference_hessian, harmonic_analysis

# Gaussian's default ("normal") optimization convergence criteria.
# Reference: Gaussian 16 documentation, Opt keyword. Forces in Eh/bohr,
# displacements (change between successive steps) in bohr. As in Gaussian, the
# geometry is converged when all four are met, or when the forces alone are
# below 1/100 of their thresholds.
MAX_FORCE = 4.5e-4
RMS_FORCE = 3.0e-4
MAX_DISP = 1.8e-3
RMS_DISP = 1.2e-3


def _gaussian_converged(g: np.ndarray, dx: np.ndarray | None) -> bool:
    max_f, rms_f = float(np.abs(g).max()), float(np.sqrt(np.mean(g ** 2)))
    if max_f < MAX_FORCE / 100 and rms_f < RMS_FORCE / 100:
        return True
    if not (max_f < MAX_FORCE and rms_f < RMS_FORCE) or dx is None:
        return False
    return float(np.abs(dx).max()) < MAX_DISP and float(np.sqrt(np.mean(dx ** 2))) < RMS_DISP


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


def optimize(engine: Engine, mol: Molecule, max_steps: int = 500, on_step=None,
             pre_optimization: dict | None = None) -> dict:
    """Quasi-Newton (L-BFGS) minimisation on analytic gradients in Cartesian
    coordinates, stopped by Gaussian's default criteria (forces and step size)."""
    extra = {"optimizer": "scipy L-BFGS-B", "max_force": MAX_FORCE, "rms_force": RMS_FORCE,
             "max_displacement": MAX_DISP, "rms_displacement": RMS_DISP}
    if pre_optimization:
        extra["pre_optimization"] = pre_optimization
    prov = _prov(engine, mol, "optimize", extra)
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
            if on_step:
                on_step({"task": "optimize", "step": len(trajectory), "energy": r.energy,
                         "max_force": trajectory[-1]["max_force"], "target_max_force": MAX_FORCE})
        return cache[key]

    x0 = (mol.coords * units.ANGSTROM_TO_BOHR).ravel()
    state = {"prev": x0.copy(), "x": None}

    def check(xk):
        key = xk.tobytes()
        g = cache[key][1] if key in cache else fg(xk)[1]
        dx = xk - state["prev"]
        state["prev"] = xk.copy()
        if _gaussian_converged(g, dx):
            state["x"] = xk.copy()
            raise StopIteration

    res = minimize(lambda x: fg(x)[0], x0, jac=lambda x: fg(x)[1], method="L-BFGS-B", callback=check,
                   options={"maxiter": max_steps, "gtol": MAX_FORCE / 100, "ftol": 1e-15, "maxcor": 30})
    xf = state["x"] if state["x"] is not None else res.x
    final = mol.copy(xf.reshape(-1, 3) * units.BOHR_TO_ANGSTROM)
    key = xf.tobytes()
    if key in cache:
        ef, g = cache[key][0], cache[key][1]
    else:
        rf = engine.compute(final)
        ef, g = rf.energy, rf.gradient.ravel()
    max_f, rms_f = float(np.abs(g).max()), float(np.sqrt(np.mean(g ** 2)))
    converged = state["x"] is not None or (max_f < MAX_FORCE and rms_f < RMS_FORCE)
    if not converged:
        prov.notes.append(f"not converged: max force {max_f:.2e}, rms {rms_f:.2e} Eh/bohr")
    return {
        "molecule": final,
        "energy_hartree": ef,
        "converged": converged,
        "max_force": max_f,
        "rms_force": rms_f,
        "steps": len(trajectory),
        "trajectory": trajectory,
        "provenance": prov.finish().to_dict(),
    }


def frequencies(engine: Engine, mol: Molecule, step_bohr: float = 0.005, on_step=None) -> dict:
    """Harmonic frequencies and IR intensities at the given geometry. Run on an
    optimized structure: frequencies are only meaningful at a stationary point."""
    analytic = getattr(engine, "supports_analytic_hessian", lambda m: False)(mol)
    if analytic:
        prov = _prov(engine, mol, "frequencies", {
            "hessian": "analytic (coupled-perturbed SCF)",
            "dipole_derivatives": "analytic atomic polar tensors"})
    else:
        prov = _prov(engine, mol, "frequencies", {"hessian": "central finite differences of analytic gradients",
                                                  "step_bohr": step_bohr})
    r0 = engine.compute(mol)
    max_f = float(np.abs(r0.gradient).max())
    if analytic:
        if on_step:
            on_step({"task": "frequencies", "step": 0, "total": 1, "analytic": True})
        hess, dmu = engine.second_derivatives(mol)
        va = harmonic_analysis(mol, hess, dmu)
        if max_f > MAX_FORCE:
            va.notes.append(f"geometry is not stationary (max force {max_f:.2e} Eh/bohr); frequencies are unreliable")
        prov.notes.extend(va.notes)
        out = va.to_dict()
        out["provenance"] = prov.finish().to_dict()
        return out
    total = 6 * len(mol.symbols)
    done = [0]

    def compute(m):
        r = engine.compute(m)
        done[0] += 1
        if on_step:
            on_step({"task": "frequencies", "step": done[0], "total": total})
        return r

    hess, dmu = finite_difference_hessian(compute, mol, step_bohr)
    va = harmonic_analysis(mol, hess, dmu)
    if max_f > MAX_FORCE:
        va.notes.append(f"geometry is not stationary (max force {max_f:.2e} Eh/bohr); frequencies are unreliable")
    prov.notes.extend(va.notes)
    out = va.to_dict()
    out["provenance"] = prov.finish().to_dict()
    return out
