"""Tier 2: GFN2-xTB.

Reference: C. Bannwarth, S. Ehlert, S. Grimme, "GFN2-xTB - An Accurate and
Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method
with Multipole Electrostatics and Density-Dependent Dispersion Contributions",
J. Chem. Theory Comput. 15, 1652-1671 (2019). DOI: 10.1021/acs.jctc.8b01176

Two interchangeable backends compute the same method:

* ``tblite``: the tblite library, called in-process.
* ``xtb-program``: the Grimme lab's ``xtb`` executable, run once per
  energy/gradient in a scratch folder.

The conda-forge tblite build crashes with an access violation inside its
first single point on Windows, so ``auto`` uses the xtb program there and
tblite elsewhere. The CI cross-check holds the two backends to 1e-5 Eh.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from .. import units
from ..molecule import Molecule
from .base import Engine, EngineResult

K_TO_HARTREE = 3.166811563e-6  # Boltzmann constant in Eh/K


def default_backend() -> str:
    forced = os.environ.get("ATHANOR_GFN2_BACKEND")
    if forced:
        return forced
    return "xtb-program" if sys.platform == "win32" else "tblite"


class GFN2xTB(Engine):
    tier = "gfn2-xtb"
    method = "GFN2-xTB"
    dipole_kind = "GFN2-xTB electronic dipole (charges + atomic dipoles)"

    def __init__(self, accuracy: float = 1.0, max_iterations: int = 250,
                 electronic_temperature: float = 300.0, backend: str = "auto"):
        self.accuracy = accuracy
        self.max_iterations = max_iterations
        self.electronic_temperature = electronic_temperature
        self.backend = default_backend() if backend == "auto" else backend
        if self.backend not in ("tblite", "xtb-program"):
            raise ValueError(f"unknown GFN2-xTB backend {self.backend!r}")
        self._version = None

    @property
    def engine_name(self) -> str:  # type: ignore[override]
        return "tblite" if self.backend == "tblite" else "xtb"

    # ---------------------------------------------------------------- version
    def version(self) -> str:
        if self._version is None:
            if self.backend == "tblite":
                from tblite.library import get_version
                self._version = ".".join(str(v) for v in get_version())
            else:
                out = subprocess.run([self._exe(), "--version"], capture_output=True, text=True,
                                     **_no_window()).stdout
                tok = [l for l in out.splitlines() if "xtb version" in l]
                self._version = tok[0].split("version")[1].split()[0] if tok else "unknown"
        return self._version

    def parameters(self) -> dict:
        return {"accuracy": self.accuracy, "max_iterations": self.max_iterations,
                "electronic_temperature_K": self.electronic_temperature, "backend": self.backend}

    # ---------------------------------------------------------------- compute
    def compute(self, mol: Molecule) -> EngineResult:
        if self.backend == "tblite":
            return self._compute_tblite(mol)
        return self._compute_program(mol)

    def _compute_tblite(self, mol: Molecule) -> EngineResult:
        from tblite.interface import Calculator
        calc = Calculator(self.method, np.array(mol.atomic_numbers),
                          mol.coords * units.ANGSTROM_TO_BOHR,
                          charge=float(mol.charge), uhf=mol.unpaired_electrons)
        calc.set("verbosity", 0)
        calc.set("accuracy", self.accuracy)
        calc.set("max-iter", self.max_iterations)
        calc.set("temperature", self.electronic_temperature * K_TO_HARTREE)
        res = calc.singlepoint()
        return EngineResult(
            energy=float(res.get("energy")),
            gradient=np.array(res.get("gradient")),
            dipole=np.array(res.get("dipole")),
            charges=np.array(res.get("charges")),
            charge_model="GFN2-xTB Mulliken-type charges",
        )

    @staticmethod
    def _exe() -> str:
        exe = shutil.which("xtb")
        if not exe:
            raise RuntimeError("the xtb program is not installed in this environment")
        return exe

    def _compute_program(self, mol: Molecule) -> EngineResult:
        n = len(mol.symbols)
        with tempfile.TemporaryDirectory(prefix="athanor-xtb-") as d:
            Path(d, "in.xyz").write_text(mol.to_xyz())
            cmd = [self._exe(), "in.xyz", "--gfn", "2", "--grad", "--json",
                   "--chrg", str(mol.charge), "--uhf", str(mol.unpaired_electrons),
                   "--acc", repr(self.accuracy), "--iterations", str(self.max_iterations),
                   "--etemp", repr(self.electronic_temperature)]
            env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "1"),
                       OMP_STACKSIZE=os.environ.get("OMP_STACKSIZE", "500M"))
            r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, env=env, **_no_window())
            out = r.stdout + r.stderr
            if r.returncode != 0 or "normal termination" not in out:
                tail = "\n".join(out.strip().splitlines()[-15:])
                raise RuntimeError(f"xtb did not terminate normally (exit {r.returncode}):\n{tail}")
            data = _read_xtb_json(Path(d, "xtbout.json"))
            missing = {"total energy", "dipole / a.u.", "partial charges"} - data.keys()
            if missing:
                raise RuntimeError(f"xtb {data.get('xtb version', '?')} output lacks {sorted(missing)}; "
                                   "Athanor needs xtb 6.7 or newer")
            grad = _read_turbomole_gradient(Path(d, "gradient"), n)
        return EngineResult(
            energy=float(data["total energy"]),
            gradient=grad,
            dipole=np.array(data["dipole / a.u."], dtype=float),
            charges=np.array(data["partial charges"], dtype=float),
            charge_model="GFN2-xTB Mulliken-type charges",
        )


def _read_xtb_json(path: Path) -> dict:
    """xtbout.json as written by xtb 6.7. On Windows its "program call" entry
    holds raw backslashes from the executable path, which is invalid JSON.
    That entry's backslashes become forward slashes; any other stray backslash
    is escaped while valid escape sequences are kept."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    lines = [l.replace("\\", "/") if l.lstrip().startswith('"program call"') else l for l in lines]
    text = re.sub(r'\\(["\\/bfnrt]|u[0-9a-fA-F]{4})?',
                  lambda m: m.group(0) if m.group(1) else "\\\\", "\n".join(lines))
    return json.loads(text)


def _read_turbomole_gradient(path: Path, n: int) -> np.ndarray:
    """Last cycle of a Turbomole ``$grad`` file: n coordinate lines, then n gradient lines (Eh/bohr)."""
    lines = path.read_text().splitlines()
    starts = [i for i, l in enumerate(lines) if l.strip().startswith("cycle")]
    if not starts:
        raise RuntimeError("xtb gradient file has no cycle block")
    block = lines[starts[-1] + 1 + n: starts[-1] + 1 + 2 * n]
    g = np.array([[float(x.replace("D", "E").replace("d", "e")) for x in l.split()[:3]] for l in block])
    if g.shape != (n, 3):
        raise RuntimeError(f"xtb gradient file: expected {n} gradient lines, read {g.shape}")
    return g


def _no_window() -> dict:
    """Keep console windows from flashing up when the desktop app runs xtb on Windows."""
    if sys.platform == "win32":
        return {"creationflags": 0x08000000}  # CREATE_NO_WINDOW
    return {}
