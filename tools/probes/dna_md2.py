"""Probe: OpenMM dynamics speed in CI, and AmberTools (tleap/sander) for an energy cross-check."""
import os, subprocess, time
import numpy as np
import openmm as mm, openmm.app as app, openmm.unit as u
from athanor_core.dna.builder import build
from athanor_core.dna import md
for name in ("tleap", "sander", "cpptraj"):
    print(name, subprocess.run(["which", name], capture_output=True, text=True).stdout.strip())
d = md.Dynamics(build("CGCGAATTCGCG", "B")["atoms"])
for k in range(4):
    f = d.run(10.0)
    print(round(f["time_ps"], 1), "ps", round(f["temperature"]), "K", round(f["ns_per_day"], 1), "ns/day", f["platform"], os.cpu_count())
leap = os.path.join(os.environ["CONDA_PREFIX"], "dat", "leap", "cmd")
print(sorted(x for x in os.listdir(leap) if "DNA" in x or "OL" in x)[:20])
