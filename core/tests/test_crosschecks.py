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
    import pyscf  # noqa: F401
    HAS_PYSCF = True
except Exception:
    HAS_PYSCF = False
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
    """Turbomole-format $vibrational spectrum written by xtb:
    mode  [symmetry]  wavenumber  IR-intensity  selection-rule column(s)
    The symmetry label is blank for translations/rotations and the number of
    selection-rule columns varies between xtb versions, so the wavenumber is
    the first numeric field after the mode index."""
    freqs = []
    for line in open(os.path.join(d, "vibspectrum"), encoding="utf-8", errors="replace"):
        p = line.split()
        if len(p) >= 3 and p[0].isdigit():
            for tok in p[1:3]:
                try:
                    freqs.append(float(tok))
                    break
                except ValueError:
                    continue
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
    HAS_TBLITE = True  # CI installs the PyPI wheel on Linux; conda-forge builds crash (see engines/xtb.py)
except Exception:
    HAS_TBLITE = False


@pytest.mark.skipif(not (HAS_XTB and HAS_TBLITE), reason="needs both tblite and the xtb program")
@pytest.mark.parametrize("smiles", ["O", "C=O", "[NH4+]", "[O]", "c1ccccc1"])
def test_gfn2_backends_agree(smiles):
    """The tblite library and the xtb program are independent GFN2-xTB codes;
    the app runs the xtb program, and tblite must reproduce it."""
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


@pytest.mark.skipif(not (HAS_PSI4 and HAS_PYSCF), reason="needs Psi4 and PySCF")
def test_pyscf_dft_matches_psi4():
    """The GPU tier runs PySCF's code path (via GPU4PySCF); on CPU it must
    reproduce Psi4 at the same functional, basis and grid size. Residual
    differences come from the two programs' grid partitioning schemes."""
    mol = tasks.optimize(get_engine("gfn2-xtb"), Molecule.from_smiles("C=O"))["molecule"]
    psi = get_engine("dft", functional="B3LYP", basis="def2-SVP", memory="1 GB",
                     scf_convergence=1e-10, d_convergence=1e-10)
    pys = get_engine("dft-pyscf", functional="B3LYP", basis="def2-SVP", density_fitting=False)
    a, b = psi.compute(mol), pys.compute(mol)
    assert a.energy == pytest.approx(b.energy, abs=2e-5)
    np.testing.assert_allclose(a.gradient, b.gradient, atol=1e-4)
    np.testing.assert_allclose(a.dipole, b.dipole, atol=1e-3)
    # Harmonic frequencies at each program's own optimized geometry.
    fa = tasks.frequencies(psi, tasks.optimize(psi, mol)["molecule"])["frequencies_cm1"]
    fb = tasks.frequencies(pys, tasks.optimize(pys, mol)["molecule"])["frequencies_cm1"]
    np.testing.assert_allclose(np.sort(fa), np.sort(fb), atol=1.5)


class _AthanorASECalculator:
    """Minimal ASE calculator around an Athanor engine, so ASE's own finite
    differences and normal-mode analysis can be compared with Athanor's."""

    def __new__(cls, engine, symbols, charge=0, multiplicity=1):
        from ase.calculators.calculator import Calculator, all_changes
        from ase.units import Bohr, Hartree

        class _Calc(Calculator):
            implemented_properties = ["energy", "forces", "dipole"]

            def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
                super().calculate(atoms, properties, system_changes)
                mol = Molecule(symbols=list(symbols), coords=self.atoms.positions.copy(),
                               charge=charge, multiplicity=multiplicity)
                r = engine.compute(mol)
                self.results = {"energy": r.energy * Hartree,
                                "forces": -np.asarray(r.gradient) * Hartree / Bohr,
                                "dipole": np.asarray(r.dipole) * Bohr}

        return _Calc()


@pytest.mark.skipif(not (HAS_ASE and HAS_XTB), reason="needs ASE and the xtb program")
def test_vibrations_match_ase():
    """Same GFN2-xTB gradients and dipoles, independent vibrational analysis (ASE)."""
    from ase import Atoms
    from ase.vibrations import Infrared
    eng = get_engine("gfn2-xtb")
    mol = tasks.optimize(eng, Molecule.from_smiles("C=O"))["molecule"]
    mine = tasks.frequencies(eng, mol)
    atoms = Atoms(mol.symbols, positions=mol.coords)
    atoms.calc = _AthanorASECalculator(eng, mol.symbols)
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


@pytest.mark.skipif(not HAS_PYSCF, reason="PySCF not installed")
def test_pyscf_engine_avoids_methods_gpu4pyscf_lacks(monkeypatch):
    """GPU4PySCF's SCF objects set mulliken_pop to NotImplemented; the engine
    must compute charges without it (the 0.2.0 GPU self-test failed on this).
    Charges must match PySCF's own Mulliken analysis."""
    from pyscf import dft
    eng = get_engine("dft-pyscf", functional="B3LYP", basis="def2-SVP", density_fitting=True)
    mol = Molecule.from_smiles("C=O")
    ref = eng._scf(mol)
    ref.kernel()
    expected = ref.mulliken_pop(verbose=0)[1]
    for cls in (dft.rks.RKS, dft.uks.UKS):
        monkeypatch.setattr(cls, "mulliken_pop", NotImplemented, raising=False)
    r = eng.compute(mol)
    np.testing.assert_allclose(r.charges, expected, atol=1e-8)
    assert abs(r.charges.sum()) < 1e-6
    # open shell (UKS)
    rad = eng.compute(Molecule.from_smiles("[CH3]"))
    assert abs(rad.charges.sum()) < 1e-6
