"""Molecule representation and conversions."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

_PT = Chem.GetPeriodicTable()


@dataclass
class Molecule:
    """Atoms, Cartesian coordinates in angstrom, total charge, spin multiplicity.

    ``bonds`` (atom index pairs with bond order) is optional connectivity used
    only for display and for the force-field tier; quantum tiers ignore it.
    """

    symbols: list[str]
    coords: np.ndarray                     # shape (N, 3), angstrom
    charge: int = 0
    multiplicity: int = 1
    bonds: list[tuple[int, int, float]] = field(default_factory=list)
    name: str = ""
    formal_charges: list[int] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.coords = np.asarray(self.coords, dtype=float).reshape(-1, 3)
        if len(self.symbols) != len(self.coords):
            raise ValueError("symbols and coords differ in length")
        n_electrons = sum(self.atomic_numbers) - self.charge
        if (n_electrons + self.multiplicity - 1) % 2 != 0:
            raise ValueError(
                f"charge {self.charge} and multiplicity {self.multiplicity} are "
                f"inconsistent with {n_electrons} electrons"
            )

    # --- basic properties ---
    @property
    def natoms(self) -> int:
        return len(self.symbols)

    @property
    def atomic_numbers(self) -> list[int]:
        return [_PT.GetAtomicNumber(s) for s in self.symbols]

    @property
    def masses(self) -> np.ndarray:
        """Mass of the most abundant isotope of each atom, in amu (standard for
        harmonic frequency calculations)."""
        return np.array([_PT.GetMostCommonIsotopeMass(s) for s in self.symbols])

    @property
    def unpaired_electrons(self) -> int:
        return self.multiplicity - 1

    def formula(self) -> str:
        counts: dict[str, int] = {}
        for s in self.symbols:
            counts[s] = counts.get(s, 0) + 1
        order = (["C", "H"] if "C" in counts else []) + sorted(k for k in counts if k not in ("C", "H") or "C" not in counts)
        return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}" for e in order)

    def copy(self, coords: Optional[np.ndarray] = None) -> "Molecule":
        return Molecule(list(self.symbols), np.array(self.coords if coords is None else coords),
                        self.charge, self.multiplicity, list(self.bonds), self.name,
                        list(self.formal_charges))

    # --- I/O ---
    def to_xyz(self, comment: str = "") -> str:
        lines = [str(self.natoms), comment or self.name]
        for s, (x, y, z) in zip(self.symbols, self.coords):
            lines.append(f"{s:<2s} {x:16.10f} {y:16.10f} {z:16.10f}")
        return "\n".join(lines) + "\n"

    @classmethod
    def from_xyz(cls, text: str, charge: int = 0, multiplicity: int = 1) -> "Molecule":
        lines = [ln for ln in text.strip().splitlines()]
        n = int(lines[0].split()[0])
        symbols, coords = [], []
        for ln in lines[2:2 + n]:
            parts = ln.split()
            symbols.append(parts[0].capitalize() if len(parts[0]) < 3 else parts[0])
            coords.append([float(v) for v in parts[1:4]])
        return cls(symbols, np.array(coords), charge, multiplicity, name=lines[1].strip())

    @classmethod
    def from_smiles(cls, smiles: str, seed: int = 0xA7A, name: str = "") -> "Molecule":
        """Build 3D coordinates with RDKit ETKDGv3 (deterministic seed). The
        geometry is only a starting point; optimize with any tier afterwards."""
        rd = Chem.MolFromSmiles(smiles)
        if rd is None:
            raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
        rd = Chem.AddHs(rd)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        if AllChem.EmbedMolecule(rd, params) != 0:
            raise RuntimeError(f"3D embedding failed for {smiles!r}")
        return cls.from_rdkit(rd, name=name or smiles)

    @classmethod
    def from_rdkit(cls, rd: Chem.Mol, conf_id: int = -1, name: str = "") -> "Molecule":
        conf = rd.GetConformer(conf_id)
        symbols = [a.GetSymbol() for a in rd.GetAtoms()]
        coords = np.array(conf.GetPositions())
        charge = Chem.GetFormalCharge(rd)
        radicals = sum(a.GetNumRadicalElectrons() for a in rd.GetAtoms())
        bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()) for b in rd.GetBonds()]
        fcs = [a.GetFormalCharge() for a in rd.GetAtoms()]
        return cls(symbols, coords, charge, radicals + 1, bonds, name, fcs)

    def to_rdkit(self) -> Chem.Mol:
        """RDKit molecule with explicit connectivity. If no bonds are stored,
        connectivity and bond orders are perceived from the geometry."""
        if self.bonds:
            rw = Chem.RWMol()
            fcs = self.formal_charges or [0] * self.natoms
            for s, fc in zip(self.symbols, fcs):
                at = Chem.Atom(s)
                at.SetFormalCharge(int(fc))
                at.SetNoImplicit(True)
                rw.AddAtom(at)
            btype = {1.0: Chem.BondType.SINGLE, 2.0: Chem.BondType.DOUBLE,
                     3.0: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC}
            for i, j, o in self.bonds:
                rw.AddBond(int(i), int(j), btype.get(float(o), Chem.BondType.SINGLE))
            mol = rw.GetMol()
            Chem.SanitizeMol(mol)
        else:
            from rdkit.Chem import rdDetermineBonds
            mol = Chem.MolFromXYZBlock(self.to_xyz())
            rdDetermineBonds.DetermineBonds(mol, charge=self.charge)
        conf = Chem.Conformer(self.natoms)
        for i, (x, y, z) in enumerate(self.coords):
            conf.SetAtomPosition(i, (float(x), float(y), float(z)))
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        return mol


# --- geometry helpers ---
def distance(coords: np.ndarray, i: int, j: int) -> float:
    return float(np.linalg.norm(coords[i] - coords[j]))


def angle(coords: np.ndarray, i: int, j: int, k: int) -> float:
    """Angle i-j-k in degrees."""
    a, b = coords[i] - coords[j], coords[k] - coords[j]
    c = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))
