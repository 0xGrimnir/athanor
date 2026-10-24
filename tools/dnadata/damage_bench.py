"""Ancient-damage formulas against gargammel's deamSim (the simulator the Ancient DNA Lab uses).

Random fragments of known sequence are damaged by deamSim with the Briggs model (-damage nu,lambda,
delta_d,delta_s), and the C->T frequency at each distance from the 5' end and G->A from the 3' end are
counted. Athanor's formulas (core/athanor_core/dna/damage.py) must predict both curves within sampling
error, for two parameter sets. This pins down what each parameter means in gargammel 1.1.4.

Usage: python damage_bench.py WORK_DIR
"""
from __future__ import annotations

import gzip
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.dna import damage  # noqa: E402

WORK = Path(sys.argv[1])
WORK.mkdir(parents=True, exist_ok=True)
results, summary = [], {}
N, L, NPOS = 200_000, 80, 25


def gate(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


rng = np.random.default_rng(5)
frags = ["".join(rng.choice(list("ACGT"), L)) for _ in range(N)]
fa = WORK / "frags.fa"
fa.write_text("".join(f">f{k}\n{s}\n" for k, s in enumerate(frags)))
worst = 0.0
for nu, lam, dd, ds in ((0.03, 0.4, 0.01, 0.3), (0.01, 0.25, 0.005, 0.5)):
    out = WORK / f"dmg_{lam}.fa"
    r = subprocess.run(["deamSim", "-damage", f"{nu},{lam},{dd},{ds}", "-o", str(out), str(fa)], capture_output=True, text=True)
    raw = out.read_bytes()
    lines = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode().splitlines()
    dmg = [x.strip().upper() for x in lines if not x.startswith(">")]
    assert len(dmg) == N, (len(dmg), r.stderr[:300])
    O = np.frombuffer("".join(frags).encode(), dtype="S1").reshape(N, L)
    D = np.frombuffer("".join(dmg).encode(), dtype="S1").reshape(N, L)
    isC, isG = O == b"C", O == b"G"
    ct = ((D == b"T") & isC).sum(0) / isC.sum(0)
    ga = ((D == b"A") & isG).sum(0) / isG.sum(0)
    ga3 = ga[::-1]                                     # by distance from the 3' end
    i = np.arange(L, dtype=float)
    pred5 = damage.read_profile(L, nu, lam, dd, ds)
    j = np.arange(L, dtype=float)
    pss3 = damage.p_single(j, lam)
    nicked = 1.0 - (1.0 - nu) ** ((L - 1 - j) + 1.0)   # a nick at or before that base, counted from the 5' end
    pred3 = ds * pss3 + dd * (1.0 - pss3) * nicked
    se5 = np.sqrt(pred5 * (1 - pred5) / isC.sum(0))
    se3 = np.sqrt(pred3 * (1 - pred3) / isG.sum(0)[::-1])
    z5 = np.abs(ct[:NPOS] - pred5[:NPOS]) / se5[:NPOS]
    z3 = np.abs(ga3[:NPOS] - pred3[:NPOS]) / se3[:NPOS]
    key = f"nu={nu},lambda={lam},dd={dd},ds={ds}"
    summary[key] = {"ct5_obs": np.round(ct[:10], 4).tolist(), "ct5_pred": np.round(pred5[:10], 4).tolist(),
                    "ga3_obs": np.round(ga3[:10], 4).tolist(), "ga3_pred": np.round(pred3[:10], 4).tolist(),
                    "max_z5": round(float(z5.max()), 2), "max_z3": round(float(z3.max()), 2)}
    print(key, json.dumps(summary[key]))
    worst = max(worst, float(z5.max()), float(z3.max()))
gate("Briggs-model formulas predict deamSim's C->T (5') and G->A (3') at every one of the first 25 positions, "
     "for two parameter sets (largest deviation <= 4.5 standard errors)", worst <= 4.5, f"largest deviation {worst:.2f} SE")
(WORK / "damage_report.json").write_text(json.dumps({"results": results, "summary": summary}, indent=1))
sys.exit(0 if all(x["pass"] for x in results) else 1)
