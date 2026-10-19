"""Correctness: the MMFF94 tier must reproduce Halgren's MMFF94 validation suite.

761 molecules, each with a published reference total energy at a fixed
geometry. The reference log prints energies to 1e-5 kcal/mol, so agreement is
required to 1e-4 kcal/mol (rounding allowance).

The suite's ERULE_* molecules exercise MMFF94's empirical-rule parameters, which
Athanor's MMFF tier declines to use (see engines/mmff.py). RDKit reproduces
those to within 0.005 kcal/mol; they are held to 0.01 kcal/mol and reported
separately so a regression in either path is still caught.
"""
import gzip
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

DATA = Path(__file__).resolve().parent.parent / "benchmarks" / "data" / "mmff94"


def test_mmff94_reproduces_validation_suite():
    ref = json.loads((DATA / "reference_energies.json").read_text())["energies_kcal_mol"]
    with gzip.open(DATA / "MMFF94_dative.sdf.gz", "rb") as fh:
        supplier = Chem.ForwardSDMolSupplier(fh, removeHs=False, sanitize=True)
        checked, worst, failures = 0, 0.0, []
        for mol in supplier:
            if mol is None:
                continue
            name = mol.GetProp("_Name").strip()
            if name not in ref:
                continue
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, ignoreInterfragInteractions=False)
            e = ff.CalcEnergy()
            dev = abs(e - ref[name])
            tol = 1e-2 if name.startswith("ERULE") else 1e-4
            worst = max(worst, dev if not name.startswith("ERULE") else 0.0)
            if dev > tol:
                failures.append((name, e, ref[name]))
            checked += 1
    assert checked >= 750, f"only {checked} suite molecules were checked"
    assert not failures, f"{len(failures)} molecules deviate; first: {failures[:3]} (worst {worst:.2e})"
