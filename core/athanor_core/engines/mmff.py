"""Tier 1: MMFF94 force field via RDKit.

Reference: T. A. Halgren, "Merck molecular force field. I-V",
J. Comput. Chem. 17, 490-641 (1996). DOI: 10.1002/(SICI)1096-987X(199604)17:5/6<490::AID-JCC1>3.0.CO;2-P
"""
from __future__ import annotations

import numpy as np
import rdkit
from rdkit.Chem import AllChem

from .. import units
from ..molecule import Molecule
from .base import Engine, EngineResult

_KCALMOL_PER_ANG_TO_HARTREE_PER_BOHR = units.BOHR_TO_ANGSTROM / units.HARTREE_TO_KCALMOL

import json as _json
from pathlib import Path as _Path

_BOND_KEYS = {tuple(k) for k in _json.loads(
    (_Path(__file__).resolve().parent.parent / "data" / "mmff94_explicit_params.json").read_text())["bond_keys"]}


class MMFF94(Engine):
    tier = "mmff94"
    engine_name = "RDKit"
    dipole_kind = "classical point-charge dipole from MMFF94 partial charges (not validated)"

    def __init__(self, variant: str = "MMFF94"):
        if variant not in ("MMFF94", "MMFF94s"):
            raise ValueError("variant must be MMFF94 or MMFF94s")
        self.variant = variant
        self.method = variant

    def version(self) -> str:
        return rdkit.__version__

    def parameters(self) -> dict:
        return {"variant": self.variant}

    def _setup(self, mol: Molecule):
        rd = mol.to_rdkit()
        if any(a.GetNumRadicalElectrons() for a in rd.GetAtoms()):
            raise ValueError("MMFF94 is defined for closed-shell molecules only; perceived bonding leaves "
                             "unpaired electrons (supply explicit connectivity, e.g. SMILES)")
        props = AllChem.MMFFGetMoleculeProperties(rd, mmffVariant=self.variant)
        if props is None:
            raise ValueError("MMFF94 has no parameters for this molecule (unsupported atom type)")
        if not AllChem.MMFFHasAllMoleculeParams(rd):
            raise ValueError("MMFF94 lacks parameters for part of this molecule")
        # RDKit silently fills missing bond-stretch parameters from MMFF94's empirical
        # rules (MMFFBNDK), which are far less accurate than the fitted table. Check
        # every bond against the fitted MMFFBOND.PAR table and refuse rather than guess.
        missing = []
        for b in rd.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            btype = props.GetMMFFBondStretchParams(rd, i, j)[0]
            ti, tj = sorted((props.GetMMFFAtomType(i), props.GetMMFFAtomType(j)))
            if (btype, ti, tj) not in _BOND_KEYS:
                missing.append(f"{rd.GetAtomWithIdx(i).GetSymbol()}{i}-{rd.GetAtomWithIdx(j).GetSymbol()}{j}"
                               f" (types {ti}-{tj})")
        if missing:
            raise ValueError("outside MMFF94's fitted parameter set (bond stretch would come from empirical "
                             "rules): " + ", ".join(missing))
        ff = AllChem.MMFFGetMoleculeForceField(rd, props, ignoreInterfragInteractions=False)
        return rd, props, ff

    def compute(self, mol: Molecule) -> EngineResult:
        rd, props, ff = self._setup(mol)
        pos = mol.coords.ravel().tolist()
        e_kcal = ff.CalcEnergy(pos)
        g_kcal_ang = np.array(ff.CalcGrad(pos)).reshape(-1, 3)
        q = np.array([props.GetMMFFPartialCharge(i) for i in range(mol.natoms)])
        com_free = mol.coords * units.ANGSTROM_TO_BOHR
        dipole = (q[:, None] * com_free).sum(0) if abs(q.sum()) < 1e-6 else None
        return EngineResult(
            energy=e_kcal / units.HARTREE_TO_KCALMOL,
            gradient=g_kcal_ang * _KCALMOL_PER_ANG_TO_HARTREE_PER_BOHR,
            dipole=dipole,
            charges=q,
            charge_model="MMFF94 partial charges",
            extras={"energy_kcal_mol": e_kcal},
        )
