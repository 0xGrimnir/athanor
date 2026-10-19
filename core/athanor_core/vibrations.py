"""Harmonic vibrational analysis shared by every tier.

Method (double-harmonic approximation), following Wilson, Decius & Cross,
"Molecular Vibrations" (1955), and the standard projection scheme used by
quantum-chemistry codes (e.g. the Gaussian white paper "Vibrational Analysis
in Gaussian", Ochterski 1999):

1. Cartesian Hessian H (Eh/bohr^2) and dipole derivatives D = d mu / d x
   (e*bohr/bohr) by central finite differences of analytic gradients and
   dipoles.
2. Mass-weight: H_mw = M^-1/2 H M^-1/2, masses in amu (most abundant isotope).
3. Remove the 6 (5 if linear) translation/rotation coordinates exactly by
   diagonalising H_mw in the orthonormal complement of that subspace.
4. Frequencies from eigenvalues; IR intensities from |D . M^-1/2 L_k|^2.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import units
from .molecule import Molecule


@dataclass
class VibrationalAnalysis:
    frequencies_cm1: np.ndarray        # signed; negative = imaginary
    ir_intensities_km_mol: np.ndarray  # may be all-NaN if no dipole derivatives
    reduced_masses_amu: np.ndarray
    normal_modes: np.ndarray           # (nmodes, N, 3) normalised Cartesian displacements
    zpe_hartree: float
    linear: bool
    n_imaginary: int
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "frequencies_cm1": [float(v) for v in self.frequencies_cm1],
            "ir_intensities_km_mol": [None if np.isnan(v) else float(v) for v in self.ir_intensities_km_mol],
            "reduced_masses_amu": [float(v) for v in self.reduced_masses_amu],
            "normal_modes": self.normal_modes.round(6).tolist(),
            "zpe_hartree": self.zpe_hartree,
            "zpe_kcal_mol": self.zpe_hartree * units.HARTREE_TO_KCALMOL,
            "linear": self.linear,
            "n_imaginary": self.n_imaginary,
            "notes": self.notes,
        }


def is_linear(coords_ang: np.ndarray, masses: np.ndarray, tol: float = 1e-3) -> bool:
    """Linear if the smallest principal moment of inertia vanishes."""
    if len(masses) == 1:
        return False
    com = (masses[:, None] * coords_ang).sum(0) / masses.sum()
    r = coords_ang - com
    inertia = np.zeros((3, 3))
    for m, v in zip(masses, r):
        inertia += m * (np.dot(v, v) * np.eye(3) - np.outer(v, v))
    w = np.linalg.eigvalsh(inertia)
    return bool(w[0] < tol * max(w[-1], 1e-12))


def _tr_basis(coords_bohr: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Orthonormal basis (3N x k) of mass-weighted translations and rotations."""
    n = len(masses)
    sq = np.sqrt(masses)
    com = (masses[:, None] * coords_bohr).sum(0) / masses.sum()
    r = coords_bohr - com
    vecs = []
    for a in range(3):                                   # translations
        v = np.zeros((n, 3)); v[:, a] = sq
        vecs.append(v.ravel())
    for a in range(3):                                   # rotations about x, y, z
        axis = np.zeros(3); axis[a] = 1.0
        v = np.cross(axis, r) * sq[:, None]
        vecs.append(v.ravel())
    b = np.array(vecs).T                                 # (3N, 6)
    u, s, _ = np.linalg.svd(b, full_matrices=False)
    keep = s > 1e-6 * s.max()
    return u[:, keep]


def harmonic_analysis(mol: Molecule, hessian: np.ndarray,
                      dipole_derivatives: np.ndarray | None = None) -> VibrationalAnalysis:
    """hessian: (3N, 3N) in Eh/bohr^2. dipole_derivatives: (3, 3N) in e (au)."""
    masses = mol.masses
    n = mol.natoms
    coords_bohr = mol.coords * units.ANGSTROM_TO_BOHR
    h = 0.5 * (hessian + hessian.T)
    inv_sqrt_m = np.repeat(1.0 / np.sqrt(masses), 3)
    h_mw = h * inv_sqrt_m[:, None] * inv_sqrt_m[None, :]          # Eh/(bohr^2 amu)

    tr = _tr_basis(coords_bohr, masses)
    # orthonormal complement of the translation/rotation space
    q, _ = np.linalg.qr(np.hstack([tr, np.eye(3 * n)]))
    comp = q[:, tr.shape[1]:3 * n]
    h_int = comp.T @ h_mw @ comp
    evals, evecs = np.linalg.eigh(h_int)
    l_mw = comp @ evecs                                           # (3N, nmodes), orthonormal

    # eigenvalue units Eh/(bohr^2 amu) -> Eh/(bohr^2 me), then to wavenumbers
    lam_me = evals / units.AMU_TO_ME
    freqs = np.sign(lam_me) * np.sqrt(np.abs(lam_me)) * units.AU_SQRT_HESSIAN_TO_CM1

    cart = (l_mw * inv_sqrt_m[:, None])                           # Cartesian displacement per unit Q
    norms = np.linalg.norm(cart, axis=0)
    reduced_masses = 1.0 / norms ** 2                             # amu, Gaussian convention
    modes = (cart / norms).T.reshape(-1, n, 3)

    if dipole_derivatives is not None:
        dmu_dq = dipole_derivatives @ cart                        # (3, nmodes) in e/amu^1/2
        ir = (dmu_dq ** 2).sum(0) * units.IR_INTENSITY_E2_PER_AMU_TO_KM_PER_MOL
    else:
        ir = np.full(len(freqs), np.nan)

    real = freqs[freqs > 0]
    # ZPE = 1/2 sum h c nu~ ; with nu~ in cm^-1 -> J -> Eh
    from scipy import constants as c
    zpe = 0.5 * float(np.sum(real)) * 100.0 * c.h * c.c / units.HARTREE_J

    notes = []
    n_imag = int(np.sum(freqs < -10.0))  # |nu| < 10 cm^-1 treated as numerical noise
    if n_imag:
        notes.append(f"{n_imag} imaginary frequency(ies): structure is not a minimum")
    return VibrationalAnalysis(freqs, ir, reduced_masses, modes, zpe,
                               is_linear(mol.coords, masses), n_imag, notes)


def finite_difference_hessian(compute, mol: Molecule, step_bohr: float = 0.005):
    """Central-difference Hessian and dipole derivatives from analytic gradients.

    ``compute(mol)`` must return an object with ``gradient`` (N,3; Eh/bohr) and
    ``dipole`` (3,; e*bohr, or None). Cost: 6N gradient evaluations. The default
    step of 0.005 bohr matches the xtb program's numerical Hessian.
    """
    n = mol.natoms
    x0 = mol.coords * units.ANGSTROM_TO_BOHR
    hess = np.zeros((3 * n, 3 * n))
    dmu = np.zeros((3, 3 * n))
    have_dipole = True
    for i in range(3 * n):
        out = []
        for sgn in (+1, -1):
            x = x0.copy().ravel(); x[i] += sgn * step_bohr
            m = mol.copy(x.reshape(-1, 3) * units.BOHR_TO_ANGSTROM)
            out.append(compute(m))
        hess[i] = (out[0].gradient.ravel() - out[1].gradient.ravel()) / (2 * step_bohr)
        if out[0].dipole is None or out[1].dipole is None:
            have_dipole = False
        else:
            dmu[:, i] = (np.asarray(out[0].dipole) - np.asarray(out[1].dipole)) / (2 * step_bohr)
    return 0.5 * (hess + hess.T), (dmu if have_dipole else None)
