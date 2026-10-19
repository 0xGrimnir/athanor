"""Smoke-test a packed runtime exactly as the installed app launches it.

Extracts env.tar.gz into a temporary folder, runs conda-unpack, then starts
`python -m athanor_core.server` with the same environment variables the Rust
shell sets, and performs real calculations on every tier.
Usage: python smoke_runtime.py <runtime_dir>
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

runtime = Path(sys.argv[1]).resolve()
work = Path(tempfile.mkdtemp(prefix="athanor-smoke-"))
env_dir = work / "env"
env_dir.mkdir()
t0 = time.time()
# Same extraction the installed app performs (system tar).
subprocess.run(["tar", "-xzf", str(runtime / "env.tar.gz"), "-C", str(env_dir)], check=True)
print(f"extracted in {time.time() - t0:.0f} s")

win = os.name == "nt"
py = env_dir / ("python.exe" if win else "bin/python")
paths = ([env_dir, env_dir / "Library/mingw-w64/bin", env_dir / "Library/usr/bin", env_dir / "Library/bin",
          env_dir / "Scripts", env_dir / "bin"] if win else [env_dir / "bin"])
env = {k: v for k, v in os.environ.items() if not k.startswith(("CONDA", "PYTHON", "MAMBA"))}
env["PATH"] = os.pathsep.join(str(p) for p in paths) + os.pathsep + os.environ.get("PATH", "")
env.update(CONDA_PREFIX=str(env_dir), PYTHONPATH=str(runtime / "core"), PYTHONNOUSERSITE="1",
           PYTHONUTF8="1", PSI_SCRATCH=str(work))

unpack_exe = env_dir / "Scripts" / "conda-unpack.exe"
unpack_py = env_dir / ("Scripts/conda-unpack-script.py" if win else "bin/conda-unpack")
if win and unpack_exe.exists():
    subprocess.run([str(unpack_exe)], env=env, check=True)
elif unpack_py.exists():
    subprocess.run([str(py), str(unpack_py)], env=env, check=True)
print("conda-unpack done")

proc = subprocess.Popen([str(py), "-m", "athanor_core.server"], cwd=str(work), env=env, text=True,
                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=sys.stderr)
print("server:", proc.stdout.readline().strip())
ids = iter(range(1, 1000))


def call(method, params=None):
    proc.stdin.write(json.dumps({"jsonrpc": "2.0", "id": next(ids), "method": method, "params": params or {}}) + "\n")
    proc.stdin.flush()
    resp = json.loads(proc.stdout.readline())
    if "error" in resp:
        raise SystemExit(f"{method} failed: {resp['error']}")
    return resp["result"]


engines = call("engines")
print("engines:", json.dumps(engines))
assert all(v["available"] for v in engines.values()), "an engine is unavailable in the packed runtime"
water = call("build", {"smiles": "O"})
summary = {}
for tier, opts in (("mmff94", {}), ("gfn2-xtb", {}), ("dft", {"functional": "B3LYP", "basis": "def2-SVP"})):
    t = time.time()
    opt = call("optimize", {"molecule": water, "tier": tier, "options": opts})
    sp = call("single_point", {"molecule": opt["molecule"], "tier": tier, "options": opts})
    fq = call("frequencies", {"molecule": opt["molecule"], "tier": tier, "options": opts})
    summary[tier] = {"converged": opt["converged"], "E": round(sp["energy_hartree"], 8),
                     "dipole_D": round(sp["dipole_debye"] or 0, 3),
                     "freqs": [round(f, 1) for f in fq["frequencies_cm1"]], "seconds": round(time.time() - t, 1)}
    assert opt["converged"] and fq["n_imaginary"] == 0, f"{tier}: {summary[tier]}"
print("SMOKE_OK", json.dumps(summary))
proc.stdin.close()
proc.wait(timeout=30)
