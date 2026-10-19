"""Physical constants and unit conversions.

All constants come from scipy.constants, which implements the CODATA
recommended values shipped with the installed SciPy version. Nothing here is
hand-typed except derived combinations, each of which is computed from those
constants rather than copied from a table.
"""
from __future__ import annotations

import math

from scipy import constants as _c

# --- base CODATA values (SI) ---
HARTREE_J = _c.physical_constants["Hartree energy"][0]          # J
BOHR_M = _c.physical_constants["Bohr radius"][0]                # m
ELECTRON_MASS_KG = _c.m_e                                       # kg
AMU_KG = _c.physical_constants["atomic mass constant"][0]       # kg
ELEMENTARY_CHARGE_C = _c.e                                      # C
SPEED_OF_LIGHT = _c.c                                           # m/s
AVOGADRO = _c.N_A                                               # 1/mol
EPSILON_0 = _c.epsilon_0                                        # F/m
HBAR = _c.hbar                                                  # J s

# --- derived conversions ---
BOHR_TO_ANGSTROM = BOHR_M * 1e10
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_KCALMOL = HARTREE_J * AVOGADRO / 4184.0
HARTREE_TO_KJMOL = HARTREE_J * AVOGADRO / 1000.0
HARTREE_TO_EV = HARTREE_J / _c.e
AMU_TO_ME = AMU_KG / ELECTRON_MASS_KG

# Debye is defined as 1e-21 / c  C*m (statcoulomb-centimetre definition).
DEBYE_CM = 1e-21 / SPEED_OF_LIGHT
AU_DIPOLE_CM = ELEMENTARY_CHARGE_C * BOHR_M                    # e*a0 in C*m
AU_TO_DEBYE = AU_DIPOLE_CM / DEBYE_CM

# sqrt(Eh / (a0^2 * me)) is the atomic unit of angular frequency for a
# mass-weighted Hessian in Eh/(bohr^2 me). Wavenumber = omega / (2 pi c).
AU_SQRT_HESSIAN_TO_CM1 = (
    math.sqrt(HARTREE_J / (BOHR_M ** 2 * ELECTRON_MASS_KG))
    / (2.0 * math.pi * SPEED_OF_LIGHT) / 100.0
)

# IR integrated absorption coefficient, double-harmonic approximation:
#   A = N_A / (12 eps0 c^2) * |d mu / d Q|^2        (SI, m/mol)
# with Q the mass-weighted normal coordinate (kg^1/2 m). For dmu/dQ given in
# e / amu^1/2 (i.e. (e*a0)/(a0*amu^1/2)), convert to SI first.
_DMUDQ_AU_TO_SI = ELEMENTARY_CHARGE_C / math.sqrt(AMU_KG)
IR_INTENSITY_E2_PER_AMU_TO_KM_PER_MOL = (
    AVOGADRO / (12.0 * EPSILON_0 * SPEED_OF_LIGHT ** 2) * _DMUDQ_AU_TO_SI ** 2 / 1000.0
)
