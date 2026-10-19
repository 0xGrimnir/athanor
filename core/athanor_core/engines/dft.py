"""Tier 3: Kohn-Sham DFT and Hartree-Fock via Psi4.

Reference: D. G. A. Smith et al., "PSI4 1.4: Open-source software for
high-throughput quantum chemistry", J. Chem. Phys. 152, 184108 (2020).
DOI: 10.1063/5.0006002
"""
from __future__ import annotations

import os
import tempfile

import numpy as np

from ..molecule import Molecule
from .base import Engine, EngineResult


class Psi4DFT(Engine):
    tier = "dft"
    engine_name = "Psi4"
    dipole_kind = "SCF electronic + nuclear dipole"

    def __init__(self, functional: str = "B3LYP", basis: str = "def2-TZVP",
                 memory: str = "2 GB", threads: int | None = None,
                 scf_convergence: float = 1e-8, d_convergence: float = 1e-8):
        self.functional = functional
        self.basis = basis
        self.method = f"{functional}/{basis}"
        self.memory = memory
        self.threads = threads or max(1, (os.cpu_count() or 2) - 1)
        self.e_conv = scf_convergence
        self.d_conv = d_convergence

    def version(self) -> str:
        import psi4
        return psi4.__version__

    def parameters(self) -> dict:
        return {"functional": self.functional, "basis": self.basis,
                "e_convergence": self.e_conv, "d_convergence": self.d_conv,
                "dft_grid": "Psi4 default (75 radial, 302 spherical)"}

    def _geometry(self, mol: Molecule) -> str:
        lines = [f"{mol.charge} {mol.multiplicity}"]
        lines += [f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in zip(mol.symbols, mol.coords)]
        lines += ["units angstrom", "symmetry c1", "no_reorient", "no_com"]
        return "\n".join(lines)

    def compute(self, mol: Molecule) -> EngineResult:
        import psi4
        quiet_psi4(psi4)
        psi4.set_memory(self.memory)
        psi4.set_num_threads(self.threads)
        psi4.core.clean_options()
        psi4.set_options({
            "basis": self.basis,
            "reference": "rks" if mol.multiplicity == 1 else "uks",
            "e_convergence": self.e_conv,
            "d_convergence": self.d_conv,
        })
        geom = psi4.geometry(self._geometry(mol))
        name = "scf" if self.functional.upper() == "HF" else self.functional
        grad, wfn = psi4.gradient(name, molecule=geom, return_wfn=True)
        psi4.oeprop(wfn, "DIPOLE", "MULLIKEN_CHARGES")
        dip = np.array(wfn.variable("SCF DIPOLE")).ravel()        # e*bohr
        charges = np.array(wfn.atomic_point_charges())
        result = EngineResult(
            energy=float(wfn.energy()),
            gradient=np.array(grad.np),
            dipole=dip,
            charges=charges,
            charge_model="Mulliken",
            extras={"homo_lumo": self._frontier(wfn)},
        )
        psi4.core.clean()
        return result

    def hessian(self, mol: Molecule):
        """Hessian from Psi4's symmetry-adapted finite differences of analytic
        gradients. Symmetry is enabled here only to reduce the number of
        displacements; Psi4 may reorient the molecule, so the geometry in Psi4's
        frame is returned with the Hessian (Eh/bohr^2)."""
        import psi4
        from .. import units
        quiet_psi4(psi4)
        psi4.set_memory(self.memory)
        psi4.set_num_threads(self.threads)
        psi4.core.clean_options()
        psi4.set_options({"basis": self.basis,
                          "reference": "rks" if mol.multiplicity == 1 else "uks",
                          "e_convergence": self.e_conv, "d_convergence": self.d_conv})
        lines = [f"{mol.charge} {mol.multiplicity}"]
        lines += [f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in zip(mol.symbols, mol.coords)]
        lines += ["units angstrom"]
        geom = psi4.geometry("\n".join(lines))
        name = "scf" if self.functional.upper() == "HF" else self.functional
        h, wfn = psi4.hessian(name, molecule=geom, return_wfn=True, dertype=1)
        g = np.array(wfn.molecule().geometry().np) * units.BOHR_TO_ANGSTROM
        symbols = [wfn.molecule().symbol(i).capitalize() for i in range(wfn.molecule().natom())]
        frame_mol = Molecule(symbols, g, mol.charge, mol.multiplicity, name=mol.name)
        psi4.core.clean()
        return frame_mol, np.array(h.np)

    @staticmethod
    def _frontier(wfn) -> dict:
        eps = np.array(wfn.epsilon_a().nph[0]) if wfn.epsilon_a().nirrep() else np.array([])
        nocc = wfn.nalpha()
        if len(eps) <= nocc:
            return {}
        return {"homo_eh": float(eps[nocc - 1]), "lumo_eh": float(eps[nocc])}


_SCRATCH: str | None = None


def quiet_psi4(psi4) -> None:
    """Discard Psi4's text output and use one scratch folder per process.

    ``psi4.core.be_quiet()`` opens the literal path /dev/null, which fails on
    Windows, so the platform's null device is used instead."""
    global _SCRATCH
    psi4.core.set_output_file(os.devnull, False)
    if _SCRATCH is None:
        base = os.environ.get("PSI_SCRATCH") or None
        _SCRATCH = tempfile.mkdtemp(prefix="athanor-psi4-", dir=base)
    psi4.core.IOManager.shared_object().set_default_path(_SCRATCH)
