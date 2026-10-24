"""DNA Structure data files from 3DNA v2.4 (BSD-3-Clause, Lu & Olson): the standard bases that define
the base reference frames (config/Atomic_?.pdb, Olson et al. 2001) and the nucleotide geometries 3DNA
rebuilds B- and A-DNA with (config/atomic/BDNA_?.pdb, ADNA_?.pdb: Clowney et al. 1996 bases, Arnott
fiber sugar-phosphate backbones with chi -110 and -160 degrees).

Usage: python make_data.py X3DNA_DIR OUT_DIR
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

X3, OUT = Path(sys.argv[1]), Path(sys.argv[2])
NOTICE = ("From 3DNA v2.4 (https://github.com/x3dna-dssr/x3dna-v2.4), Copyright (c) 1999-2026, Xiang-Jun Lu "
          "and Wilma K. Olson, BSD 3-Clause License. Cite: Lu & Olson, Nucleic Acids Res. 31, 5108 (2003); "
          "Nature Protocols 3, 1213 (2008).")
RENAME = {"O1P": "OP1", "O2P": "OP2", "C7": "C5M"}


def read(path: Path) -> dict:
    atoms = {}
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip().replace("*", "'")
            atoms[RENAME.get(name, name)] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    return atoms


OUT.mkdir(parents=True, exist_ok=True)
bases = {b: read(X3 / "config" / f"Atomic_{b}.pdb") for b in "ACGTU"}
(OUT / "standard_bases.json").write_text(json.dumps({"notice": NOTICE, "source": "config/Atomic_?.pdb",
                                                    "bases": bases}, indent=1))
templates = {}
for form, prefix in (("B", "BDNA"), ("A", "ADNA")):
    templates[form] = {b: read(X3 / "config" / "atomic" / f"{prefix}_{b}.pdb") for b in "ACGT"}
(OUT / "templates.json").write_text(json.dumps({"notice": NOTICE, "source": "config/atomic/BDNA_?.pdb, ADNA_?.pdb",
                                               **templates}, indent=1))
print({b: len(a) for b, a in bases.items()}, {f: {b: len(a) for b, a in t.items()} for f, t in templates.items()})
