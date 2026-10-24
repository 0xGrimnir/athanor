"""Probe: gargammel deamSim's Briggs damage on known fragments - per-position C->T and G->A."""
import random, subprocess, shutil, os
import numpy as np
print(shutil.which("deamSim"))
print(subprocess.run(["deamSim"], capture_output=True, text=True).stderr[:3000])
random.seed(1)
N, L = 100000, 60
with open("frags.fa", "w") as f:
    for k in range(N):
        f.write(f">f{k}\n" + "".join(random.choice("ACGT") for _ in range(L)) + "\n")
nu, lam, dd, ds = 0.03, 0.4, 0.01, 0.3
r = subprocess.run(["deamSim", "-damage", f"{nu},{lam},{dd},{ds}", "-name", "-o", "out.fa", "frags.fa"], capture_output=True, text=True)
print("rc", r.returncode, r.stderr[:800])
if not os.path.exists("out.fa"):
    r = subprocess.run(f"deamSim -damage {nu},{lam},{dd},{ds} frags.fa > out.fa", shell=True, capture_output=True, text=True)
    print("rc2", r.returncode, r.stderr[:800])
orig = [l.strip() for l in open("frags.fa") if not l.startswith(">")]
out_lines = open("out.fa").read().splitlines()
print(out_lines[:6])
dmg = [l.strip() for l in out_lines if not l.startswith(">")]
print(len(orig), len(dmg))
ct5 = np.zeros(L); c5 = np.zeros(L); ga3 = np.zeros(L); g3 = np.zeros(L)
lower = 0
for o, d in zip(orig, dmg):
    lower += sum(ch.islower() for ch in d)
    d = d.upper()
    for i in range(L):
        if o[i] == "C":
            c5[i] += 1; ct5[i] += d[i] == "T"
        j = L - 1 - i
        if o[j] == "G":
            g3[i] += 1; ga3[i] += d[j] == "A"
p5, p3 = ct5 / c5, ga3 / g3
print("lowercase letters:", lower)
print("C->T from 5':", np.round(p5[:15], 4).tolist(), "middle", round(float(p5[25:35].mean()), 4))
print("G->A from 3':", np.round(p3[:15], 4).tolist(), "middle", round(float(p3[25:35].mean()), 4))
i = np.arange(L)
for name, f in (("(1-l)^(i+1)", (1 - lam) ** (i + 1)), ("l^(i+1)", lam ** (i + 1)), ("(1-l)^i", (1 - lam) ** i), ("l^i", lam ** i),
                ("l(1-l)^i...geom P(len>i)", (1 - lam) ** (i + 1))):
    pred = ds * f + dd * (1 - f)
    print(name, np.round(pred[:8], 4).tolist(), "max |diff| 5'", round(float(np.max(np.abs(pred[:20] - p5[:20]))), 4))
