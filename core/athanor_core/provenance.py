"""Provenance records: what produced a number, so any result can be reproduced."""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import platform
import sys
from dataclasses import asdict, dataclass, field
from typing import Any

from . import __version__


@dataclass
class Provenance:
    tier: str                         # "mmff94" | "gfn2-xtb" | "dft"
    engine: str                       # e.g. "RDKit", "tblite", "Psi4"
    engine_version: str
    method: str                       # e.g. "MMFF94", "GFN2-xTB", "B3LYP/def2-TZVP"
    task: str                         # "energy" | "optimize" | "frequencies" | ...
    parameters: dict[str, Any] = field(default_factory=dict)
    input_hash: str = ""
    athanor_version: str = __version__
    python: str = field(default_factory=lambda: sys.version.split()[0])
    platform: str = field(default_factory=lambda: f"{platform.system()} {platform.machine()}")
    started_utc: str = field(default_factory=lambda: _dt.datetime.now(_dt.timezone.utc).isoformat())
    finished_utc: str = ""
    notes: list[str] = field(default_factory=list)

    def finish(self) -> "Provenance":
        self.finished_utc = _dt.datetime.now(_dt.timezone.utc).isoformat()
        return self

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def hash_input(symbols, coords, charge, multiplicity, extra: dict | None = None) -> str:
    """Stable SHA-256 of the physical input (coordinates rounded to 1e-8 A)."""
    payload = {
        "symbols": list(symbols),
        "coords": [[round(float(v), 8) for v in row] for row in coords],
        "charge": int(charge),
        "multiplicity": int(multiplicity),
        "extra": extra or {},
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
