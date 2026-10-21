"""JSON-RPC 2.0 over stdin/stdout, one JSON object per line.

The desktop app spawns ``python -m athanor_core.server`` and talks to it
through pipes only: no network socket is opened. Methods:

  ping                          -> {"version": ...}
  engines                       -> tiers available on this machine, with versions
  build {smiles}                -> molecule (ETKDGv3 3D coordinates)
  parse_xyz {xyz, charge?, multiplicity?}
  single_point {molecule, tier, options?}
  optimize {molecule, tier, options?}
  frequencies {molecule, tier, options?}
  validation                    -> published benchmark comparators for the UI's error bars
  adna_catalog / adna_run       -> Ancient DNA Lab, simulated data (Linux engines pack)
  adna_real_run {inputs, ...}   -> Ancient DNA Lab, the user's own reads (Linux engines pack)
  adna_export {path, what}      -> save the last real-data consensus or report
"""
from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

import numpy as np

from . import __version__, tasks
from .engines import get_engine
from .molecule import Molecule

REF = Path(__file__).resolve().parent.parent / "benchmarks" / "reference" / "phase1.json"


def mol_to_json(m: Molecule) -> dict:
    return {"symbols": m.symbols, "coords": np.round(m.coords, 8).tolist(), "charge": m.charge,
            "multiplicity": m.multiplicity, "bonds": [[int(i), int(j), float(o)] for i, j, o in m.bonds],
            "formal_charges": m.formal_charges, "name": m.name, "formula": m.formula()}


def mol_from_json(d: dict) -> Molecule:
    return Molecule(d["symbols"], np.array(d["coords"], dtype=float), int(d.get("charge", 0)),
                    int(d.get("multiplicity", 1)), [tuple(b) for b in d.get("bonds", [])],
                    d.get("name", ""), list(d.get("formal_charges", [])))


def _engine(params: dict):
    opts = params.get("options") or {}
    return get_engine(params["tier"], **opts)


def m_ping(_):
    return {"version": __version__}


def m_engines(_):
    out = {}
    for tier, opts in (("mmff94", {}), ("gfn2-xtb", {}), ("dft", {})):
        try:
            e = get_engine(tier, **opts)
            out[tier] = {"available": True, "engine": e.engine_name, "version": e.version(), "method": e.method}
        except Exception as ex:
            out[tier] = {"available": False, "reason": f"{type(ex).__name__}: {ex}"}
    return out


def m_gpu_status(_):
    from .engines.pyscf_dft import gpu_available
    ok, detail = gpu_available()
    return {"available": ok, "detail": detail}


def m_gpu_selftest(p):
    from .gpu_selftest import run
    return run(quick=bool(p.get("quick", False)))


def m_adna_catalog(_):
    from .adna.pipeline import reference_catalog, tool_versions
    return {"references": reference_catalog(), "tools": tool_versions()}


def m_adna_run(p):
    """Simulate ancient reads from a mitogenome and reconstruct it against references.
    Runs only inside the Linux engines pack (the tools are Linux programs)."""
    import os
    from .adna.pipeline import run
    allowed = {"endogenous", "references", "coverage", "damage", "contamination", "length_loc",
               "length_scale", "min_length", "read_length", "mapdamage_stats", "min_mapq", "iterations"}
    kw = {k: v for k, v in (p or {}).items() if k in allowed}
    if "references" in kw:
        kw["references"] = tuple(kw["references"])
    if "damage" in kw:
        kw["damage"] = tuple(float(x) for x in kw["damage"])
    work = os.path.join(os.path.expanduser("~"), "athanor-adna", "latest")
    return run(work, on_step=_progress({"tier": "adna"}), **kw)


_LAST_REAL: dict = {}


def to_linux_path(p: str) -> str:
    """Windows paths from the desktop app's file dialogs, as seen inside the
    WSL2 engines pack (drives are mounted under /mnt)."""
    import re
    s = str(p)
    for pre in ("\\\\?\\",):
        if s.startswith(pre):
            s = s[len(pre):]
    m = re.match(r"^([A-Za-z]):[\\/](.*)$", s)
    if m:
        return f"/mnt/{m.group(1).lower()}/" + m.group(2).replace("\\", "/")
    m = re.match(r"^[\\/]{2}wsl(?:\.localhost|\$)[\\/][^\\/]+[\\/](.*)$", s, re.IGNORECASE)
    if m:
        return "/" + m.group(1).replace("\\", "/")
    return s


def m_adna_real_run(p):
    """Reconstruct a mitogenome from the user's own reads (FASTQ or BAM).
    Files are read where they are; nothing leaves the computer."""
    import os
    from .adna.pipeline import read_fasta
    from .adna.realdata import run_real
    p = p or {}
    inputs = []
    for inp in p.get("inputs", []):
        d = {"kind": inp["kind"], "path": to_linux_path(inp["path"])}
        if inp.get("path2"):
            d["path2"] = to_linux_path(inp["path2"])
        for k in ("trimmed", "read_length"):
            if k in inp:
                d[k] = inp[k]
        inputs.append(d)
    if not inputs:
        raise ValueError("no read files given")
    kw = {k: p[k] for k in ("iterations", "min_mapq", "min_length", "min_depth", "min_confidence",
                            "mapdamage_stats", "circular", "sample") if k in p}
    if p.get("reference_path"):
        kw["reference_fasta"] = to_linux_path(p["reference_path"])
    else:
        kw["reference"] = p.get("reference", "NC_007596")
    kw["contaminant"] = "NC_012920" if p.get("contamination_check", True) else None
    if p.get("compare_path"):
        name, seq = read_fasta(to_linux_path(p["compare_path"]))
        kw["compare_with"] = (name, seq)
    work = os.path.join(os.path.expanduser("~"), "athanor-adna", "real")
    res = run_real(work, inputs, on_step=_progress({"tier": "adna"}), **kw)
    _LAST_REAL.clear()
    _LAST_REAL.update(res)
    return res


def m_adna_export(p):
    """Write the last real-data result: consensus FASTA, majority-rule FASTA or full report (JSON)."""
    if not _LAST_REAL:
        raise ValueError("no reconstruction to export yet")
    path = to_linux_path(p["path"])
    what = p.get("what", "fasta")
    R = _LAST_REAL
    if what == "fasta":
        text = R["consensus"]["fasta"]
    elif what == "majority_fasta":
        seq = R["majority_rule"]["sequence"]
        text = (f">{R['sample']} majority-rule consensus ({R['majority_rule']['rule']}; reference "
                f"{R['reference'].get('accession', '')})\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n")
    elif what == "json":
        text = json.dumps(R, indent=1)
    else:
        raise ValueError(f"unknown export {what!r}")
    with open(path, "w", newline="\n") as fh:
        fh.write(text)
    return {"path": p["path"], "bytes": len(text)}


def m_build(p):
    return mol_to_json(Molecule.from_smiles(p["smiles"], name=p.get("name", p["smiles"])))


def m_parse_xyz(p):
    return mol_to_json(Molecule.from_xyz(p["xyz"], int(p.get("charge", 0)), int(p.get("multiplicity", 1))))


def m_parse_molblock(p):
    from rdkit import Chem
    rd = Chem.MolFromMolBlock(p["molblock"], removeHs=False)
    if rd is None:
        raise ValueError("could not parse MOL/SDF block")
    if rd.GetNumConformers() == 0 or not rd.GetConformer().Is3D():
        raise ValueError("MOL/SDF file has no 3D coordinates")
    return mol_to_json(Molecule.from_rdkit(rd, name=rd.GetProp("_Name") if rd.HasProp("_Name") else ""))


def m_single_point(p):
    return tasks.single_point(_engine(p), mol_from_json(p["molecule"]))


_PROTO = None


def _notify(method: str, params: dict) -> None:
    """JSON-RPC notification (no id) on the protocol channel, e.g. progress."""
    if _PROTO is not None:
        _PROTO.write(json.dumps({"jsonrpc": "2.0", "method": method, "params": params}, default=float) + "\n")
        _PROTO.flush()


def _progress(p):
    tier = p.get("tier")
    return lambda d: _notify("progress", {**d, "tier": tier})


def m_optimize(p):
    r = tasks.optimize(_engine(p), mol_from_json(p["molecule"]), on_step=_progress(p),
                       pre_optimization=p.get("pre_optimization"))
    r["molecule"] = mol_to_json(r["molecule"])
    return r


def m_frequencies(p):
    return tasks.frequencies(_engine(p), mol_from_json(p["molecule"]), on_step=_progress(p))


def m_validation(_):
    ref = json.loads(REF.read_text())
    return {"comparators": ref["comparators"], "sources": ref["sources"]}


METHODS = {"ping": m_ping, "engines": m_engines, "build": m_build, "parse_xyz": m_parse_xyz,
           "parse_molblock": m_parse_molblock,
           "single_point": m_single_point, "optimize": m_optimize, "frequencies": m_frequencies,
           "validation": m_validation, "gpu_status": m_gpu_status, "gpu_selftest": m_gpu_selftest,
           "adna_catalog": m_adna_catalog, "adna_run": m_adna_run,
           "adna_real_run": m_adna_real_run, "adna_export": m_adna_export}


def handle(line: str) -> dict:
    try:
        req = json.loads(line)
    except json.JSONDecodeError as ex:
        return {"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": f"parse error: {ex}"}}
    rid = req.get("id")
    fn = METHODS.get(req.get("method"))
    if fn is None:
        return {"jsonrpc": "2.0", "id": rid, "error": {"code": -32601, "message": f"unknown method {req.get('method')!r}"}}
    try:
        return {"jsonrpc": "2.0", "id": rid, "result": fn(req.get("params") or {})}
    except Exception as ex:
        return {"jsonrpc": "2.0", "id": rid, "error": {"code": -32000, "message": f"{type(ex).__name__}: {ex}",
                                                       "data": traceback.format_exc(limit=4)}}


def main():
    global _PROTO
    # Engines print to stdout in places; keep the protocol channel clean.
    proto = _PROTO = sys.stdout
    sys.stdout = sys.stderr
    proto.write(json.dumps({"jsonrpc": "2.0", "method": "ready", "params": {"version": __version__}}) + "\n")
    proto.flush()
    for line in sys.stdin:
        if not line.strip():
            continue
        resp = handle(line)
        proto.write(json.dumps(resp, default=float) + "\n")
        proto.flush()


if __name__ == "__main__":
    main()
