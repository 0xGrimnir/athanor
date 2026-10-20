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
                 conv_tol: float = 1e-10, max_memory_mb: int = 8000,
                 hessian: str = "finite-difference"):
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
        if hessian not in ("finite-difference", "analytic"):
            raise ValueError("hessian must be 'finite-difference' or 'analytic'")
        self.hessian = hessian
        self.method = f"{functional}/{basis}" + ("/RI" if self.density_fitting else "")
        # Converged density of the previous call, reused as the initial guess
        # for the next geometry (same atoms, charge and spin). This changes
        # only how fast the SCF converges, not what it converges to.
        self._guess: tuple | None = None

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
             "conv_tol": self.conv_tol, "hessian": self.hessian}
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
        key = (tuple(mol.symbols), mol.charge, mol.multiplicity)
        dm0 = self._guess[1] if self._guess and self._guess[0] == key else None
        energy = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
        if not mf.converged:
            raise RuntimeError(f"SCF did not converge in {mf.max_cycle} cycles")
        grad = _np(mf.nuc_grad_method().kernel())
        dip = _np(mf.dip_moment(unit="AU", verbose=0))  # about the coordinate origin
        charges = mulliken_charges(mf)
        self._guess = (key, np.array(_np(mf.make_rdm1()), dtype=float))  # plain array: no stale MO tags
        mo_e = _np(mf.mo_energy)
        occ = _np(mf.mo_occ)
        extras = {}
        if mo_e.ndim == 1 and (occ > 0).any() and (occ == 0).any():
            extras["homo_lumo"] = {"homo_eh": float(mo_e[occ > 0].max()), "lumo_eh": float(mo_e[occ == 0].min())}
        return EngineResult(energy=energy, gradient=grad.reshape(-1, 3), dipole=dip.ravel(),
                            charges=charges.ravel(), charge_model="Mulliken", extras=extras)


    # ------------------------------------------------------------ second derivatives
    def supports_analytic_hessian(self, mol: Molecule) -> bool:
        return self.hessian == "analytic" and mol.multiplicity == 1

    def second_derivatives(self, mol: Molecule) -> tuple[np.ndarray, np.ndarray]:
        """Analytic Hessian (3N, 3N; Eh/bohr^2) and dipole derivatives (3, 3N; e)
        at ``mol``, from one SCF plus the coupled-perturbed (CPHF/CPKS) equations."""
        mf = self._scf(mol)
        mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF did not converge")
        return analytic_hessian_and_apt(mf, gpu=(self.device == "gpu"))


def analytic_hessian_and_apt(mf, gpu: bool) -> tuple[np.ndarray, np.ndarray]:
    """Closed-shell analytic Hessian and atomic polar tensors (dipole derivatives).

    The polar tensors follow GPU4PySCF's infrared routine
    (gpu4pyscf/properties/ir.py, Apache-2.0): the first-order density from the
    CPKS response orbitals contracted with the dipole integrals, plus the
    derivative of the dipole integrals and the nuclear term. GPU4PySCF stores
    the response orbitals in the MO basis, PySCF in the AO basis; that is the
    only difference between the two paths. CI checks the CPU path against
    finite differences, and the GPU self-test checks the GPU path the same way.
    """
    mol = mf.mol
    natm, nao = mol.natm, mol.nao
    hobj = mf.Hessian()
    h = _np(hobj.kernel())                                   # (natm, natm, 3, 3)
    H = h.transpose(0, 2, 1, 3).reshape(3 * natm, 3 * natm)
    if gpu:
        import cupy as xp
    else:
        xp = np
    C = xp.asarray(mf.mo_coeff)
    occ = xp.asarray(mf.mo_occ)
    eps = xp.asarray(mf.mo_energy)
    atmlst = range(natm)
    h1 = hobj.make_h1(C, occ, None, atmlst)
    fx = hobj.gen_vind(C, occ) if hasattr(hobj, "gen_vind") and callable(getattr(hobj, "gen_vind")) else None
    mo1, _ = hobj.solve_mo1(eps, C, occ, h1, fx, atmlst, hobj.max_memory, 0)
    mo1 = _np(xp.asarray(mo1))  # as in GPU4PySCF's IR routine
    C, occ = _np(C), _np(occ)
    mocc = C[:, occ > 0]
    dm0 = _np(mf.make_rdm1())
    with mol.with_common_orig((0, 0, 0)):
        hmu = mol.intor("int1e_r")                                 # (3, nao, nao)
        hmu11 = -mol.intor("int1e_irp").reshape(3, 3, nao, nao)
    apt = np.zeros((3, 3, natm))
    for ia, (_s0, _s1, p0, p1) in enumerate(mol.aoslice_by_atom()):
        h11 = np.zeros((3, 3, nao, nao))
        h11[:, :, :, p0:p1] += hmu11[:, :, :, p0:p1]
        h11[:, :, p0:p1] += hmu11[:, :, :, p0:p1].transpose(0, 1, 3, 2)
        if gpu:   # mo1[ia]: (3, nmo, nocc), MO basis
            dm1 = np.stack([C @ mo1[ia, y] @ mocc.T for y in range(3)])
        else:     # mo1[ia]: (3, nao, nocc), AO basis
            dm1 = np.stack([mo1[ia, y] @ mocc.T for y in range(3)])
        apt[:, :, ia] = (-4.0 * np.einsum("xuv,yuv->xy", hmu, dm1)
                         - np.einsum("xyuv,vu->xy", h11, dm0)
                         + mol.atom_charge(ia) * np.eye(3))
    dmu = apt.transpose(0, 2, 1).reshape(3, 3 * natm)             # dmu[x, 3*atom + y]
    return 0.5 * (H + H.T), dmu


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
