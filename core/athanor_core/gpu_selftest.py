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


if __name__ == "__main__":
    import json
    import sys
    rep = run(quick="--quick" in sys.argv)
    print(json.dumps(rep, indent=2))
    sys.exit(0 if rep["passed"] else 1)
