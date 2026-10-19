"""Tier 2: GFN2-xTB via tblite.

Reference: C. Bannwarth, S. Ehlert, S. Grimme, "GFN2-xTB - An Accurate and
Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method
with Multipole Electrostatics and Density-Dependent Dispersion Contributions",
J. Chem. Theory Comput. 15, 1652-1671 (2019). DOI: 10.1021/acs.jctc.8b01176
"""
from __future__ import annotations

import numpy as np

from .. import units
from ..molecule import Molecule
from .base import Engine, EngineResult


class GFN2xTB(Engine):
    tier = "gfn2-xtb"
    engine_name = "tblite"
    method = "GFN2-xTB"
    dipole_kind = "GFN2-xTB electronic dipole (charges + atomic dipoles)"

    def __init__(self, accuracy: float = 1.0, max_iterations: int = 250,
                 electronic_temperature: float = 300.0):
        self.accuracy = accuracy
        self.max_iterations = max_iterations
        self.electronic_temperature = electronic_temperature

    def version(self) -> str:
        from tblite.library import get_version
        return ".".join(str(v) for v in get_version())

    def parameters(self) -> dict:
        return {"accuracy": self.accuracy, "max_iterations": self.max_iterations,
                "electronic_temperature_K": self.electronic_temperature}

    def compute(self, mol: Molecule) -> EngineResult:
        from tblite.interface import Calculator
        calc = Calculator(self.method, np.array(mol.atomic_numbers),
                          mol.coords * units.ANGSTROM_TO_BOHR,
                          charge=float(mol.charge), uhf=mol.unpaired_electrons)
        calc.set("verbosity", 0)
        calc.set("accuracy", self.accuracy)
        calc.set("max-iter", self.max_iterations)
        calc.set("temperature", self.electronic_temperature * 3.166811563e-6)  # K -> Eh (k_B)
        res = calc.singlepoint()
        return EngineResult(
            energy=float(res.get("energy")),
            gradient=np.array(res.get("gradient")),
            dipole=np.array(res.get("dipole")),
            charges=np.array(res.get("charges")),
            charge_model="GFN2-xTB Mulliken-type charges",
        )
