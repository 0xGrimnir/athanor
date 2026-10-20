"""Kohn-Sham DFT via PySCF, on the CPU or on an NVIDIA GPU via GPU4PySCF.

The GPU tier is what the optional compute pack runs (Linux only, so on Windows
it runs inside WSL2). The CPU path is the same code with ``device="cpu"`` and
is what CI checks against Psi4, because CI machines have no GPU.

References:
  Q. Sun et al., "Recent developments in the PySCF program package",
  J. Chem. Phys. 153, 024109 (2020). DOI: 10.1063/5.0006074
  R. Li, Q. Sun, X. Zhang, G. K.-L. Chan, "Introducing GPU acceleration into
  the Python-based simulations of chemistry framework", J. Phys. Chem. A 129,
  1459-1468 (2025). DOI: 10.1021/acs.jpca.4c05876
"""
from __future__ import annotations

import numpy as np

from ..molecule import Molecule
from .base import Engine, EngineResult


def _np(x) -> np.ndarray:
    """CuPy arrays (GPU) to NumPy; NumPy passes through."""
    return np.asarray(x.get() if hasattr(x, "get") else x)


class PySCFDFT(Engine):
    dipole_kind = "SCF electronic + nuclear dipole"

    def __init__(self, functional: str = "B3LYP", basis: str = "def2-TZVP", device: str = "gpu",
                 density_fitting: bool | None = None, auxbasis: str = "def2-universal-jkfit",
                 radial_points: int = 75, angular_points: int = 302,
                 conv_tol: float = 1e-10, max_memory_mb: int = 8000):
        if device not in ("gpu", "cpu"):
            raise ValueError("device must be 'gpu' or 'cpu'")
        self.functional = functional
        self.basis = basis
        self.device = device
        # Density fitting (RI-J/K) is GPU4PySCF's fast path; on CPU it is off by
        # default so the CPU path can be compared directly with Psi4's exact integrals.
        self.density_fitting = (device == "gpu") if density_fitting is None else density_fitting
        self.auxbasis = auxbasis
        self.radial_points = radial_points
        self.angular_points = angular_points
        self.conv_tol = conv_tol
        self.max_memory_mb = max_memory_mb
        self.method = f"{functional}/{basis}" + ("/RI" if self.density_fitting else "")

    @property
    def tier(self) -> str:  # type: ignore[override]
        return "dft-gpu" if self.device == "gpu" else "dft-pyscf"

    @property
    def engine_name(self) -> str:  # type: ignore[override]
        return "GPU4PySCF" if self.device == "gpu" else "PySCF"

    def version(self) -> str:
        import pyscf
        if self.device == "gpu":
            import gpu4pyscf
            return f"{gpu4pyscf.__version__} (PySCF {pyscf.__version__})"
        return pyscf.__version__

    def parameters(self) -> dict:
        p = {"functional": self.functional, "basis": self.basis, "device": self.device,
             "density_fitting": self.density_fitting,
             "grid": f"{self.radial_points} radial x {self.angular_points} angular (Lebedev), "
                     "Treutler-Ahlrichs radii, Becke partition",
             "conv_tol": self.conv_tol}
        if self.density_fitting:
            p["auxbasis"] = self.auxbasis
        return p

    def _scf(self, mol: Molecule):
        from pyscf import dft, gto
        pmol = gto.M(atom=[(s, tuple(map(float, x))) for s, x in zip(mol.symbols, mol.coords)],
                     unit="Angstrom", charge=mol.charge, spin=mol.multiplicity - 1,
                     basis=self.basis, verbose=0, max_memory=self.max_memory_mb)
        mf = dft.RKS(pmol) if mol.multiplicity == 1 else dft.UKS(pmol)
        mf.xc = self.functional
        mf.grids.atom_grid = (self.radial_points, self.angular_points)
        mf.conv_tol = self.conv_tol
        mf.max_cycle = 200
        if self.density_fitting:
            mf = mf.density_fit(auxbasis=self.auxbasis)
        if self.device == "gpu":
            mf = mf.to_gpu()
        return mf

    def compute(self, mol: Molecule) -> EngineResult:
        mf = self._scf(mol)
        energy = float(mf.kernel())
        if not mf.converged:
            raise RuntimeError(f"SCF did not converge in {mf.max_cycle} cycles")
        grad = _np(mf.nuc_grad_method().kernel())
        dip = _np(mf.dip_moment(unit="AU", verbose=0))  # about the coordinate origin
        charges = mulliken_charges(mf)
        mo_e = _np(mf.mo_energy)
        occ = _np(mf.mo_occ)
        extras = {}
        if mo_e.ndim == 1 and (occ > 0).any() and (occ == 0).any():
            extras["homo_lumo"] = {"homo_eh": float(mo_e[occ > 0].max()), "lumo_eh": float(mo_e[occ == 0].min())}
        return EngineResult(energy=energy, gradient=grad.reshape(-1, 3), dipole=dip.ravel(),
                            charges=charges.ravel(), charge_model="Mulliken", extras=extras)


def mulliken_charges(mf) -> np.ndarray:
    """Mulliken charges from the SCF density, computed by PySCF's CPU routine
    on both paths. GPU4PySCF's SCF objects deliberately leave
    ``mulliken_pop`` unimplemented (``NotImplemented``), so the method is
    never called on them."""
    from pyscf.scf import hf as hf_cpu
    dm = _np(mf.make_rdm1())
    s = mf.mol.intor_symmetric("int1e_ovlp")
    return np.asarray(hf_cpu.mulliken_pop(mf.mol, dm, s, verbose=0)[1]).ravel()


def gpu_available() -> tuple[bool, str]:
    """Whether GPU4PySCF can run here, with the reason when it cannot."""
    try:
        import cupy
        n = cupy.cuda.runtime.getDeviceCount()
        if n < 1:
            return False, "no CUDA device visible"
        props = cupy.cuda.runtime.getDeviceProperties(0)
        name = props["name"].decode() if isinstance(props["name"], bytes) else str(props["name"])
        import gpu4pyscf  # noqa: F401
        return True, f"{name} (compute capability {props['major']}.{props['minor']})"
    except Exception as e:  # noqa: BLE001 - report any failure as unavailable
        msg = f"{type(e).__name__}: {e}"
        if "InsufficientDriver" in msg or "driver version is insufficient" in msg:
            msg += (" -- the Windows NVIDIA driver is older than this pack's CUDA 13 runtime; "
                    "update the NVIDIA driver (580 or newer) and run the self-test again")
        return False, msg
