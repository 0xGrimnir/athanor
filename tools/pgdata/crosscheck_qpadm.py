"""Check Athanor's qpAdm against ADMIXTOOLS qpAdm on a well-powered simulation.

This is the correctness test for the estimator itself, kept separate from the
reference-pack validation: it simulates a population history whose admixture
proportions are known, with enough samples and SNPs that the weights are
pinned down (standard errors well under 0.02). Two correct implementations
must then land on the same answer, and both must recover the truth.

Usage: python crosscheck_qpadm.py ADMIXTOOLS_BIN [WORKDIR]
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import msprime
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.personal.fstats import blocks_from_cm, qpadm  # noqa: E402

BIN = Path(sys.argv[1])
WORK = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(tempfile.mkdtemp())
WORK.mkdir(parents=True, exist_ok=True)
TRUE = [0.5, 0.3, 0.2]
NPOP = 20
fails = []


def gate(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)
    if not ok:
        fails.append(name)


def hashit(s: str) -> int:
    h = 0
    for ch in s.encode():
        h = (h * 23 + ch) & 0xFFFFFFFF
    return h


def hasharr(names) -> int:
    h = 0
    for n in names:
        h = ((h * 17) & 0xFFFFFFFF) ^ hashit(n)
    return h


d = msprime.Demography()
for p in ["O", "R1", "R2", "R3", "R4", "R5", "S1", "S2", "S3", "T", "A12", "A123", "AN", "ROOT", "A45"]:
    d.add_population(name=p, initial_size=10000)
d.add_admixture(time=40, derived="T", ancestral=["S1", "S2", "S3"], proportions=TRUE)
d.add_population_split(time=600, derived=["S1"], ancestral="R1")
d.add_population_split(time=600, derived=["S2"], ancestral="R2")
d.add_population_split(time=600, derived=["S3"], ancestral="R3")
d.add_population_split(time=1200, derived=["R1", "R2"], ancestral="A12")
d.add_population_split(time=1600, derived=["A12", "R3"], ancestral="A123")
d.add_population_split(time=1400, derived=["R4", "R5"], ancestral="A45")
d.add_population_split(time=2000, derived=["A123", "A45"], ancestral="AN")
d.add_population_split(time=4000, derived=["AN", "O"], ancestral="ROOT")
d.sort_events()
pops = ["O", "R1", "R2", "R3", "R4", "R5", "S1", "S2", "S3", "T"]
samples = {p: NPOP for p in pops}
rng = np.random.default_rng(3)
G, POS, CH = [], [], []
for c in range(1, 11):
    ts = msprime.sim_ancestry(samples=samples, demography=d, sequence_length=20_000_000,
                              recombination_rate=1e-8, random_seed=c)
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=100 + c, model=msprime.BinaryMutationModel())
    gm = ts.genotype_matrix()
    pos = ts.tables.sites.position.astype(np.int64)
    keep = rng.random(len(pos)) < 0.25
    G.append((gm[keep][:, 0::2] + gm[keep][:, 1::2]).astype(np.int8))
    POS.append(pos[keep] + 1)
    CH.append(np.full(int(keep.sum()), c))
    del gm, ts
G = np.vstack(G)
POS = np.concatenate(POS)
CH = np.concatenate(CH)
nsnp = G.shape[0]
names = [(f"{p}_{i}", p) for p in pops for i in range(NPOP)]
snp_names = [f"rs{i}" for i in range(nsnp)]
print(f"simulated {nsnp} SNPs, {len(names)} individuals", flush=True)

pre = WORK / "cc"
with open(str(pre) + ".ind", "w") as f:
    for n, p in names:
        f.write(f"{n} U {p}\n")
with open(str(pre) + ".snp", "w") as f:
    for i in range(nsnp):
        f.write(f"{snp_names[i]} {CH[i]} {POS[i] * 1e-8:.6f} {POS[i]} A G\n")
rlen = (nsnp + 3) // 4
with open(str(pre) + ".geno", "wb") as f:
    f.write(f"TGENO {len(names)} {nsnp} {hasharr([n for n, _ in names]):x} {hasharr(snp_names):x}"
            .encode().ljust(48, b"\x00"))
    pad = np.full((len(names), rlen * 4), 3, np.uint8)
    pad[:, :nsnp] = (2 - G).T                      # copies of allele 1
    q = pad.reshape(len(names), rlen, 4)
    f.write(((q[:, :, 0] << 6) | (q[:, :, 1] << 4) | (q[:, :, 2] << 2) | q[:, :, 3]).tobytes())

(WORK / "left").write_text("T\nS1\nS2\nS3\n")
(WORK / "right").write_text("O\nR1\nR2\nR3\nR4\nR5\n")
(WORK / "par").write_text(f"genotypename: {pre}.geno\nsnpname: {pre}.snp\nindivname: {pre}.ind\n"
                          f"popleft: {WORK / 'left'}\npopright: {WORK / 'right'}\ndetails: YES\nallsnps: NO\n")
out = subprocess.run([str(BIN / "qpAdm"), "-p", str(WORK / "par")], capture_output=True, text=True,
                     timeout=7200).stdout
line = [l for l in out.splitlines() if l.strip().startswith("best coefficients:")]
err = [l for l in out.splitlines() if l.strip().startswith("std. errors:")]
if not line:
    print(out[-2000:])
    raise SystemExit("ADMIXTOOLS qpAdm produced no coefficients")
w_ad = np.array([float(x) for x in line[0].split(":")[1].split()])
se_ad = np.array([float(x) for x in err[0].split(":")[1].split()]) if err else np.full(3, np.nan)

lab = np.array([p for _, p in names])
F = np.vstack([G[:, lab == p].mean(1) / 2 for p in pops])
blk = blocks_from_cm(CH, POS * 1e-6, 5.0)
mine = qpadm(F, blk, 9, [6, 7, 8], [0, 1, 2, 3, 4, 5])
w_me, se_me = np.array(mine["weights"]), np.array(mine["se"])
print(f"   ADMIXTOOLS {np.round(w_ad, 4)} (se {np.round(se_ad, 4)})")
print(f"   Athanor    {np.round(w_me, 4)} (se {np.round(se_me, 4)}), p {mine['p']:.3f}, "
      f"{mine['snps']} SNPs, {mine['blocks']} blocks")
gate("the simulation is well powered (every standard error below 0.02)", float(np.max(se_me)) < 0.02,
     f"largest standard error {np.max(se_me):.4f}")
gate("Athanor's qpAdm weights match ADMIXTOOLS within 0.005", float(np.max(np.abs(w_me - w_ad))) <= 0.005,
     f"largest difference {np.max(np.abs(w_me - w_ad)):.5f}")
gate("Athanor's standard errors match ADMIXTOOLS within 0.005", float(np.max(np.abs(se_me - se_ad))) <= 0.005,
     f"largest difference {np.nanmax(np.abs(se_me - se_ad)):.5f}")
# The sources are proxies: they drift on after the admixture, which biases the
# fitted proportions slightly. Both implementations carry the same bias, so this
# gate allows for it while still catching anything larger.
gate("both recover the simulated proportions within 0.03",
     float(max(np.max(np.abs(w_me - TRUE)), np.max(np.abs(w_ad - TRUE)))) <= 0.03,
     f"Athanor {np.round(np.abs(w_me - TRUE), 4)}, ADMIXTOOLS {np.round(np.abs(w_ad - TRUE), 4)}")
summary = (f"Athanor's qpAdm reproduces ADMIXTOOLS qpAdm {np.round(w_ad, 3).tolist()} with "
           f"{np.round(w_me, 3).tolist()} on a well-powered simulation: every weight agrees to within "
           f"{np.max(np.abs(w_me - w_ad)):.4f} and every standard error to within "
           f"{np.nanmax(np.abs(se_me - se_ad)):.4f}, and both recover the simulated proportions "
           f"{TRUE} to within {max(np.max(np.abs(w_me - TRUE)), np.max(np.abs(w_ad - TRUE))):.3f}.")
(WORK / "crosscheck.json").write_text(json.dumps({"ancient_ancestry": summary, "fails": fails}, indent=1))
print(f"{len(fails)} gate(s) failed: {fails}")
sys.exit(1 if fails else 0)
