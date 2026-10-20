"""Acceptance test for the GPU tier, run on the user's own machine.

CI has no GPU, so the GPU code path cannot be verified there. Before the GPU
tier is offered, this test computes the same molecules with GPU4PySCF and with
CPU PySCF at identical settings (functional, basis, density fitting, grid,
convergence) and requires agreement well inside chemical significance. The
CPU path is itself checked against Psi4 on every CI build.
"""
from __future__ import annotations

import time

import numpy as np

from .engines.pyscf_dft import PySCFDFT, gpu_available
from .molecule import Molecule

# Tolerances for GPU vs CPU at identical settings. Both converge the SCF to
# 1e-10 Eh; differences come from floating-point summation order only.
TOL = {"energy_hartree": 1e-6, "gradient_hartree_per_bohr": 1e-5, "dipole_au": 1e-4, "charges_e": 1e-4}

CASES = [
    ("water", "O", "B3LYP", "def2-TZVP"),
    ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "B3LYP", "def2-SVP"),
]


def run(quick: bool = False) -> dict:
    ok, detail = gpu_available()
    report = {"gpu": detail, "tolerances": TOL, "cases": [], "passed": False}
    if not ok:
        report["error"] = f"GPU not usable: {detail}"
        return report
    passed = True
    for name, smiles, xc, basis in (CASES[:1] if quick else CASES):
        try:
            report["cases"].append(_case(name, smiles, xc, basis))
        except Exception as e:  # noqa: BLE001 - report, don't crash the core
            import traceback
            tb = traceback.extract_tb(e.__traceback__)[-1]
            report["cases"].append({"molecule": name, "passed": False,
                                    "error": f"{type(e).__name__}: {e} ({tb.filename.rsplit('/', 1)[-1]}:{tb.lineno})"})
        passed &= report["cases"][-1]["passed"]
    report["passed"] = passed
    # Analytic Hessian + IR on the GPU vs finite differences of GPU gradients.
    # Failing this does not disable the GPU tier; frequencies then use
    # finite differences, as before.
    if passed and not quick:
        try:
            report["analytic_hessian"] = _analytic_case()
        except Exception as e:  # noqa: BLE001
            report["analytic_hessian"] = {"passed": False, "error": f"{type(e).__name__}: {e}"}
    if not passed and not report.get("error"):
        errs = [c["error"] for c in report["cases"] if c.get("error")]
        if errs:
            report["error"] = errs[0]
    return report


def _case(name: str, smiles: str, xc: str, basis: str) -> dict:
    mol = Molecule.from_smiles(smiles, name=name)
    PySCFDFT(xc, basis, device="gpu", density_fitting=True).compute(mol)  # warm-up: CUDA context, kernels
    # Fresh instances, so neither run starts from a converged guess.
    gpu = PySCFDFT(xc, basis, device="gpu", density_fitting=True)
    cpu = PySCFDFT(xc, basis, device="cpu", density_fitting=True)
    t0 = time.perf_counter(); rg = gpu.compute(mol); tg = time.perf_counter() - t0
    t0 = time.perf_counter(); rc = cpu.compute(mol); tc = time.perf_counter() - t0
    diff = {
        "energy_hartree": abs(rg.energy - rc.energy),
        "gradient_hartree_per_bohr": float(np.abs(rg.gradient - rc.gradient).max()),
        "dipole_au": float(np.abs(rg.dipole - rc.dipole).max()),
    }
    diff["charges_e"] = float(np.abs(rg.charges - rc.charges).max())
    case_ok = all(diff[k] <= TOL[k] for k in TOL)
    return {
        "molecule": name, "atoms": len(mol.symbols), "method": gpu.method,
        "energy_gpu": rg.energy, "energy_cpu": rc.energy, "differences": diff,
        "seconds_gpu": round(tg, 2), "seconds_cpu": round(tc, 2),
        "speedup": round(tc / tg, 1) if tg > 0 else None, "passed": case_ok,
    }


HESS_TOL = {"frequency_cm1": 2.0, "ir_relative": 0.03, "ir_absolute_km_mol": 0.5}


def _analytic_case() -> dict:
    from .tasks import frequencies, optimize
    from .vibrations import finite_difference_hessian, harmonic_analysis
    fd = PySCFDFT("B3LYP", "def2-SVP", device="gpu", density_fitting=True)
    an = PySCFDFT("B3LYP", "def2-SVP", device="gpu", density_fitting=True, hessian="analytic")
    mol = optimize(fd, Molecule.from_smiles("C=O", name="formaldehyde"))["molecule"]
    t0 = time.perf_counter(); h, dmu = finite_difference_hessian(fd.compute, mol, 0.005); t_fd = time.perf_counter() - t0
    t0 = time.perf_counter(); fa = frequencies(an, mol); t_an = time.perf_counter() - t0
    ff = harmonic_analysis(mol, h, dmu)
    f1, f2 = np.array(ff.frequencies_cm1), np.array(fa["frequencies_cm1"])
    i1, i2 = np.array(ff.ir_intensities_km_mol), np.array(fa["ir_intensities_km_mol"])
    dfreq = float(np.abs(f1 - f2).max())
    ok_ir = bool(np.all(np.abs(i1 - i2) <= HESS_TOL["ir_absolute_km_mol"] + HESS_TOL["ir_relative"] * np.abs(i1)))
    return {"molecule": "formaldehyde", "method": an.method, "max_frequency_difference_cm1": round(dfreq, 3),
            "ir_finite_difference": [round(x, 2) for x in i1], "ir_analytic": [round(x, 2) for x in i2],
            "seconds_finite_difference": round(t_fd, 1), "seconds_analytic": round(t_an, 1),
            "tolerances": HESS_TOL, "passed": dfreq <= HESS_TOL["frequency_cm1"] and ok_ir}


if __name__ == "__main__":
    import json
    import sys
    rep = run(quick="--quick" in sys.argv)
    print(json.dumps(rep, indent=2))
    sys.exit(0 if rep["passed"] else 1)
