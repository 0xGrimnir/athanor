"""Correctness cross-checks against independent implementations.

* GFN2-xTB: Athanor (GFN2 gradients + Athanor vibrational analysis) vs the xtb
  program's own single point and numerical Hessian; and the tblite library vs
  the xtb program as interchangeable backends.
* DFT: Athanor's finite-difference Hessian from Psi4 gradients vs Psi4's own
  frequency analysis.
* Vibrational analysis: vs ASE's independent implementation.

Each check skips when its reference program is not installed; CI installs all.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from athanor_core import tasks, units
from athanor_core.engines import get_engine
from athanor_core.engines.dft import quiet_psi4
from athanor_core.molecule import Molecule
from athanor_core.vibrations import harmonic_analysis

HAS_XTB = shutil.which("xtb") is not None
try:
    import psi4  # noqa: F401
    HAS_PSI4 = True
except Exception:
    HAS_PSI4 = False
try:
    import ase  # noqa: F401
    HAS_ASE = True
except Exception:
    HAS_ASE = False


def _xtb_run(mol: Molecule, *args) -> tuple[str, str]:
    d = tempfile.mkdtemp(prefix="athanor-xtb-")
    open(os.path.join(d, "in.xyz"), "w").write(mol.to_xyz())
    env = dict(os.environ, OMP_NUM_THREADS="1")
    r = subprocess.run(["xtb", "in.xyz", "--gfn", "2", "--chrg", str(mol.charge), *args],
                       cwd=d, capture_output=True, text=True, env=env)
    return r.stdout + r.stderr, d


def _xtb_vibspectrum(d: str) -> np.ndarray:
    """Turbomole-format $vibrational spectrum written by xtb 6.7:
    mode  [symmetry]  wavenumber  IR-intensity  IR-selection
    The symmetry column is blank for translations/rotations, so the
    wavenumber is read as the 3rd field from the end."""
    freqs = []
    for line in open(os.path.join(d, "vibspectrum"), encoding="utf-8", errors="replace"):
        p = line.split()
        if len(p) >= 4 and p[0].isdigit():
            try:
                freqs.append(float(p[-3]))
            except ValueError:
                pass
    return np.array([f for f in freqs if abs(f) > 1.0])


@pytest.mark.skipif(not HAS_XTB, reason="xtb program not installed")
@pytest.mark.parametrize("smiles", ["O", "C=O", "N"])
def test_gfn2_matches_xtb_program(smiles):
    eng = get_engine("gfn2-xtb")
    mol = tasks.optimize(eng, Molecule.from_smiles(smiles))["molecule"]
    mine_e = eng.compute(mol).energy
    out, d = _xtb_run(mol, "--hess")
    m = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", out)
    assert m, out[-2000:]
    assert mine_e == pytest.approx(float(m.group(1)), abs=1e-5)
    theirs = np.sort(_xtb_vibspectrum(d))
    mine = np.sort(np.array(tasks.frequencies(eng, mol)["frequencies_cm1"]))
    mine = mine[np.abs(mine) > 1.0]
    assert len(mine) == len(theirs)
    np.testing.assert_allclose(mine, theirs, atol=3.0)  # cm-1; two independent GFN2 codes


try:
    import tblite.interface  # noqa: F401
    HAS_TBLITE = sys.platform != "win32"  # conda-forge tblite crashes on Windows; see engines/xtb.py
except Exception:
    HAS_TBLITE = False


@pytest.mark.skipif(not (HAS_XTB and HAS_TBLITE), reason="needs both tblite and the xtb program")
@pytest.mark.parametrize("smiles", ["O", "C=O", "[NH4+]", "[O]", "c1ccccc1"])
def test_gfn2_backends_agree(smiles):
    """The tblite library and the xtb program are independent GFN2-xTB codes;
    the app uses whichever runs on the platform, so they must agree."""
    mol = Molecule.from_smiles(smiles)
    a = get_engine("gfn2-xtb", backend="tblite").compute(mol)
    b = get_engine("gfn2-xtb", backend="xtb-program").compute(mol)
    assert a.energy == pytest.approx(b.energy, abs=1e-6)
    np.testing.assert_allclose(a.gradient, b.gradient, atol=1e-5)
    np.testing.assert_allclose(a.dipole, b.dipole, atol=1e-3)
    np.testing.assert_allclose(a.charges, b.charges, atol=1e-3)


@pytest.mark.skipif(not HAS_PSI4, reason="Psi4 not installed")
def test_fd_hessian_matches_psi4_frequency_analysis():
    import psi4
    eng = get_engine("dft", functional="B3LYP", basis="cc-pVDZ", memory="1 GB")
    mol = tasks.optimize(eng, Molecule.from_smiles("O"))["molecule"]
    mine = np.sort(np.array(tasks.frequencies(eng, mol)["frequencies_cm1"]))
    frame, hess = eng.hessian(mol)
    via_psi4_hessian = np.sort(harmonic_analysis(frame, hess).frequencies_cm1)
    # Psi4's own vibrational analysis on the same level of theory
    quiet_psi4(psi4); psi4.core.clean_options()
    psi4.set_options({"basis": "cc-pVDZ", "e_convergence": 1e-8, "d_convergence": 1e-8})
    geom = psi4.geometry("\n".join(["0 1"] + [f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in zip(mol.symbols, mol.coords)] + ["units angstrom"]))
    e, wfn = psi4.frequency("b3lyp", molecule=geom, return_wfn=True, dertype=1)
    theirs = np.sort(np.array(wfn.frequencies().np))
    np.testing.assert_allclose(via_psi4_hessian, theirs, atol=0.5)
    np.testing.assert_allclose(mine, theirs, atol=2.0)


@pytest.mark.skipif(not HAS_ASE, reason="ASE not installed")
def test_vibrations_match_ase():
    from ase import Atoms
    from ase.vibrations import Infrared
    from tblite.ase import TBLite
    eng = get_engine("gfn2-xtb")
    mol = tasks.optimize(eng, Molecule.from_smiles("C=O"))["molecule"]
    mine = tasks.frequencies(eng, mol)
    atoms = Atoms(mol.symbols, positions=mol.coords)
    atoms.calc = TBLite(method="GFN2-xTB", verbosity=0)
    cwd = os.getcwd(); d = tempfile.mkdtemp(); os.chdir(d)
    try:
        ir = Infrared(atoms, delta=0.005 * units.BOHR_TO_ANGSTROM, nfree=2)
        ir.run()
        f = ir.get_frequencies().real
        inten = ir.intensities * 42.255  # (D/A)^2/amu -> km/mol
    finally:
        os.chdir(cwd)
    keep = np.argsort(np.abs(f))[-len(mine["frequencies_cm1"]):]
    order = np.argsort(f[keep])
    np.testing.assert_allclose(np.sort(mine["frequencies_cm1"]), f[keep][order], atol=3.0)
    np.testing.assert_allclose(np.array(mine["ir_intensities_km_mol"])[np.argsort(mine["frequencies_cm1"])],
                               inten[keep][order], rtol=0.05, atol=0.5)
