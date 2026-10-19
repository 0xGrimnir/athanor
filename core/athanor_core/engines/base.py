"""Common engine interface. Every tier returns the same result shape in atomic units."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ..molecule import Molecule


@dataclass
class EngineResult:
    energy: float                              # Eh
    gradient: np.ndarray                       # (N, 3) Eh/bohr
    dipole: Optional[np.ndarray] = None        # (3,) e*bohr, about the origin of the input frame
    charges: Optional[np.ndarray] = None       # (N,) e, scheme named in ``charge_model``
    charge_model: str = ""
    extras: dict = field(default_factory=dict)


class Engine:
    tier: str = ""
    engine_name: str = ""
    method: str = ""
    #: what the dipole from this engine means, shown next to the value in the UI
    dipole_kind: str = ""

    def version(self) -> str:
        raise NotImplementedError

    def compute(self, mol: Molecule) -> EngineResult:
        raise NotImplementedError

    def parameters(self) -> dict:
        return {}

    def __call__(self, mol: Molecule) -> EngineResult:
        return self.compute(mol)
