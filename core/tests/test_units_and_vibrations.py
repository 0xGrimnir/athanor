"""Unit tests for constants and the shared vibrational analysis."""
import numpy as np
import pytest

from athanor_core import units
from athanor_core.molecule import Molecule
from athanor_core.vibrations import harmonic_analysis


def test_known_conversion_factors():
    # Literature values these derived constants must reproduce.
    assert units.BOHR_TO_ANGSTROM == pytest.approx(0.529177210903, rel=1e-9)
    assert units.AU_TO_DEBYE == pytest.approx(2.541746, rel=1e-5)
    assert units.HARTREE_TO_KCALMOL == pytest.approx(627.5095, rel=1e-6)
    # Commonly tabulated IR conversion: 974.88 km/mol per (e^2/amu)
    assert units.IR_INTENSITY_E2_PER_AMU_TO_KM_PER_MOL == pytest.approx(974.88, rel=1e-4)
    # sqrt(Eh/(a0^2 me))/(2 pi c) = 5140.487 cm^-1 * sqrt(me/amu) relation:
    assert units.AU_SQRT_HESSIAN_TO_CM1 / np.sqrt(units.AMU_TO_ME) == pytest.approx(5140.487, rel=1e-5)


def test_diatomic_harmonic_oscillator_is_exact():
    """For a diatomic with a pure quadratic bond potential k/2 (r-re)^2 the
    single vibrational frequency is analytic: nu = sqrt(k/mu)/(2 pi c)."""
    mol = Molecule(["H", "Cl"], np.array([[0, 0, 0], [0, 0, 1.2746]]))
    k = 0.32  # Eh/bohr^2
    # Hessian of k/2 (z2 - z1)^2 in Cartesians
    h = np.zeros((6, 6))
    h[2, 2] = h[5, 5] = k
    h[2, 5] = h[5, 2] = -k
    va = harmonic_analysis(mol, h)
    assert va.linear and len(va.frequencies_cm1) == 1
    m1, m2 = mol.masses
    mu = m1 * m2 / (m1 + m2) * units.AMU_TO_ME
    expected = np.sqrt(k / mu) * units.AU_SQRT_HESSIAN_TO_CM1
    assert va.frequencies_cm1[0] == pytest.approx(expected, rel=1e-9)
    # Reduced mass follows the Gaussian/Psi4 convention 1/sum_i(l_i^2/m_i) for the
    # normalised mass-weighted mode l (Ochterski 1999), not the classical mu:
    # for a diatomic it equals m1*m2*(m1+m2)/(m1^2+m2^2) (H2 -> 1.008 amu).
    assert va.reduced_masses_amu[0] == pytest.approx(m1 * m2 * (m1 + m2) / (m1**2 + m2**2), rel=1e-9)


def test_translation_rotation_removed_for_nonlinear():
    mol = Molecule(["O", "H", "H"], np.array([[0, 0, 0], [0.958, 0, 0], [-0.24, 0.927, 0]]))
    va = harmonic_analysis(mol, np.zeros((9, 9)) + 0.0)
    assert not va.linear and len(va.frequencies_cm1) == 3
