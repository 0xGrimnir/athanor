"""Ancient-damage formulas against gargammel's deamSim (the simulator the Ancient DNA Lab uses).

Random fragments of known sequence are damaged by deamSim with the Briggs model (-damage nu,lambda,
delta_d,delta_s), and the C->T frequency at each distance from the 5' end and G->A from the 3' end are
counted. Athanor's formulas (core/athanor_core/dna/damage.py) must predict both curves within sampling
error. The Briggs model as published is checked against gargammel 1.1.5's corrected -damage (built from
source at a fixed commit); gargammel 1.1.4's model - which the Ancient DNA Lab's simulations use - against
1.1.4 itself and 1.1.5's -damagelegacy.

Usage: python damage_bench.py WORK_DIR [DEAMSIM_1.1.5_OR_LATER]
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
O = np.frombuffer("".join(frags).encode(), dtype="S1").reshape(N, L)
isC, isG = O == b"C", O == b"G"
NEW = sys.argv[2] if len(sys.argv) > 2 else None        # deamSim built from gargammel >= 1.1.5
# (program, option, Athanor model): bioconda's 1.1.4 is what the Ancient DNA Lab runs
RUNS = [("deamSim", "-damage", "legacy", "gargammel 1.1.4")]
if NEW:
    RUNS += [(NEW, "-damage", "briggs", "gargammel 1.1.5"), (NEW, "-damagelegacy", "legacy", "gargammel 1.1.5")]
worst = {}
for prog, opt, model, label in RUNS:
    for nu, lam, dd, ds in ((0.024, 0.36, 0.0097, 0.68), (0.03, 0.4, 0.01, 0.3)):
        out = WORK / "dmg.fa"
        out.unlink(missing_ok=True)
        r = subprocess.run([prog, opt, f"{nu},{lam},{dd},{ds}", "-o", str(out), str(fa)], capture_output=True, text=True)
        raw = out.read_bytes()
        lines = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode().splitlines()
        dmg = [x.strip().upper() for x in lines if not x.startswith(">")]
        assert len(dmg) == N, (len(dmg), r.stderr[:300])
        D = np.frombuffer("".join(dmg).encode(), dtype="S1").reshape(N, L)
        ct = ((D == b"T") & isC).sum(0) / isC.sum(0)
        ga3 = (((D == b"A") & isG).sum(0) / isG.sum(0))[::-1]
        pred5, pred3 = damage.read_profile(L, nu, lam, dd, ds, model=model)
        se5 = np.sqrt(pred5 * (1 - pred5) / isC.sum(0))
        se3 = np.sqrt(pred3 * (1 - pred3) / isG.sum(0)[::-1])
        z = max(float((np.abs(ct[:NPOS] - pred5[:NPOS]) / se5[:NPOS]).max()), float((np.abs(ga3[:NPOS] - pred3[:NPOS]) / se3[:NPOS]).max()))
        key = f"{label} {opt} ({model}) nu={nu},lambda={lam},dd={dd},ds={ds}"
        summary[key] = {"ct5_obs": np.round(ct[:8], 4).tolist(), "ct5_pred": np.round(pred5[:8], 4).tolist(),
                        "ga3_obs": np.round(ga3[:8], 4).tolist(), "ga3_pred": np.round(pred3[:8], 4).tolist(), "max_z": round(z, 2)}
        print(key, json.dumps(summary[key]), flush=True)
        worst[f"{label} {opt}"] = max(worst.get(f"{label} {opt}", 0.0), z)
gate("Athanor's damage formulas predict deamSim's C->T (5') and G->A (3') at every one of the first 25 positions: the "
     "Briggs (2007) model against gargammel 1.1.5's -damage, gargammel 1.1.4's model against its -damage and 1.1.5's "
     "-damagelegacy; two parameter sets each (largest deviation <= 4.5 standard errors)",
     NEW is not None and max(worst.values()) <= 4.5, "; ".join(f"{k}: {v:.2f} SE" for k, v in worst.items()))
(WORK / "damage_report.json").write_text(json.dumps({"results": results, "summary": summary}, indent=1))
sys.exit(0 if all(x["pass"] for x in results) else 1)
