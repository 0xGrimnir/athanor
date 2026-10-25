"""Smoke-test a packed runtime exactly as the installed app launches it.

Extracts env.tar.gz into a temporary folder, relocates it with unpack_env.py, then starts
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

# Same relocation step the installed app runs.
t0 = time.time()
r = subprocess.run([str(py), "-I", "-S", str(runtime / "unpack_env.py"), str(env_dir), str(runtime / "relocate.json")],
                   env=env, capture_output=True, text=True)
print(r.stdout.strip().splitlines()[-1] if r.stdout.strip() else "", r.stderr.strip())
print(f"relocated in {time.time() - t0:.1f} s")
r.check_returncode()
print("relocation done")

proc = subprocess.Popen([str(py), "-m", "athanor_core.server"], cwd=str(work), env=env, text=True,
                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=sys.stderr)
print("server:", proc.stdout.readline().strip())
ids = iter(range(1, 1000))


def call(method, params=None, expect_error=False):
    proc.stdin.write(json.dumps({"jsonrpc": "2.0", "id": next(ids), "method": method, "params": params or {}}) + "\n")
    proc.stdin.flush()
    while True:  # skip notifications (progress) until the response with our id
        resp = json.loads(proc.stdout.readline())
        if "id" in resp:
            break
    if "error" in resp:
        if expect_error:
            return resp["error"]
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
# DNA Structure: the shipped parameter files load and build each form
for seq, form in (("CGCGAATTCGCG", "B"), ("GGGGCCCC", "A"), ("CGCGCGCGCGCG", "Z")):
    h = call("dna_build", {"sequence": seq, "form": form})
    tw = sum(s["params"][5] for s in h["analysis"]["steps"]) / len(h["analysis"]["steps"])
    summary[f"dna_{form}"] = {"atoms": len(h["atoms"]), "twist": round(tw, 1)}
    assert len(h["atoms"]) > 20 * len(seq) and (tw < 0) == (form == "Z"), summary[f"dna_{form}"]
v = call("dna_variants")["variants"]
assert v and call("dna_variant", {"rsid": v[0]})["variant"]["rsid"] == v[0]
mt = call("dna_melt", {"sequence": "CGCGAATTCGCG"})
summary["dna_melt_tm"] = round(mt["tm"], 1)
rx = call("dna_relax", {"sequence": "CGCGAATTCGCG"})
summary["dna_relax_o3p"] = [round(min(rx["o3p"]), 3), round(max(rx["o3p"]), 3), rx["relax"]["platform"]]
assert 1.55 < min(rx["o3p"]) and max(rx["o3p"]) < 1.67, summary["dna_relax_o3p"]
sh = call("dna_shape", {"sequence": "CGCGAATTCGCG"})
summary["dna_shape_mgw"] = sh["shape"]["MGW"][:3]
assert abs(sh["shape"]["MGW"][0] - 5.4325) < 0.01
en = call("dna_energies", {"sequence": "GCAT"})
summary["dna_energies"] = en["pairs"]
assert en["pairs"][0] < en["pairs"][2] < 0                      # G-C binds more strongly than A-T
md0 = call("dna_md_start", {"sequence": "CGCG", "ps": 0.2})
md1 = call("dna_md_step", {"run": md0["run"], "ps": 0.2})
call("dna_md_stop", {"run": md0["run"]})
summary["dna_md"] = [md1["platform"], round(md1["ns_per_day"], 1)]
for name in ("dna_nucleosome", "dna_fibre", "dna_nucleus"):
    r = call(name, {"sequence": "CGCGAATTCGCG"} if name == "dna_nucleosome" else {})
    summary[name] = True
summary["dna_complexes"] = len(call("dna_scales")["complexes"])
# Deep Time runs in the Linux engines pack (it needs pysam); the Windows runtime
# must list the methods and say so clearly rather than failing obscurely.
err = call("deeptime_inspect", {"path": "none.bam"}, expect_error=True)
msg = (err.get("message", "") if isinstance(err, dict) else str(err))
summary["deeptime"] = msg[:140]
# either the pack is missing (the usual case on Windows) or the file is; both
# must say so plainly rather than failing obscurely
assert any(x in msg.lower() for x in ("linux engines pack", "no such file", "none.bam")), summary["deeptime"]
print("SMOKE_OK", json.dumps(summary))
proc.stdin.close()
proc.wait(timeout=30)
