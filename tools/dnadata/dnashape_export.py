"""Export Deep DNAshape's trained weights for Athanor's NumPy implementation, and check it.

1. Clone Deep DNAshape (BSD-3-Clause) at a fixed commit and read every model checkpoint with
   TensorFlow's checkpoint reader; write the arrays as raw little-endian float32 files with a JSON
   index (byte-for-byte reproducible, unlike .npz archives, which carry timestamps).
2. Gate: Athanor's NumPy network gives the same predictions as the original TensorFlow predictor, for
   every exported feature, on random sequences of several lengths and at flanking layers 2, 4 and 7.
3. Report: agreement of the predictions with the experimental and simulation data the authors
   published (Zenodo 10403307, CC BY 4.0). Those data were used to build the model, so this shows the
   port reproduces the published model's behaviour; it is not an independent test of Deep DNAshape.

Usage: python dnashape_export.py OUT_DIR WORK_DIR
"""
from __future__ import annotations

import io
import json
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

import numpy as np

COMMIT = "fe12e2444cb74d10bc870251345936d8adb65cde"
FEATURES = ["MGW", "ProT", "Roll", "HelT", "EP", "Shear", "Stretch", "Stagger", "Buckle", "Opening",
            "Shift", "Slide", "Rise", "Tilt", "MGW-FL", "ProT-FL", "Roll-FL", "HelT-FL"]
OUT, WORK = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
results = []


def gate(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


src = WORK / "deepDNAshape"
if not src.exists():
    subprocess.run(["git", "clone", "-q", "https://github.com/JinsenLi/deepDNAshape.git", str(src)], check=True)
subprocess.run(["git", "-C", str(src), "checkout", "-q", COMMIT], check=True)

import tensorflow as tf  # noqa: E402

index = {"source": "Deep DNAshape (Li, Chiu & Rohs, Nature Commun. 15, 1243, 2024), "
                   f"https://github.com/JinsenLi/deepDNAshape commit {COMMIT}, BSD-3-Clause",
         "params": json.loads((src / "deepDNAshape" / "params.json").read_text()), "features": {}}
for f in FEATURES:
    r = tf.train.load_checkpoint(str(src / "deepDNAshape" / "models" / f))
    names = sorted(k for k in r.get_variable_to_shape_map() if k.endswith("/.ATTRIBUTES/VARIABLE_VALUE")
                   and "optimizer" not in k.lower())
    blob, arrays, off = [], {}, 0
    for k in names:
        a = np.asarray(r.get_tensor(k), dtype="<f4")
        key = k.replace("/.ATTRIBUTES/VARIABLE_VALUE", "")
        arrays[key] = [off, list(a.shape)]
        blob.append(a.ravel())
        off += a.size
    fn = f"{f}.f32"
    np.concatenate(blob).astype("<f4").tofile(OUT / fn)
    index["features"][f] = {"file": fn, "arrays": arrays}
(OUT / "index.json").write_text(json.dumps(index, indent=1, sort_keys=True) + "\n")
shutil.copy(src / "LICENSE", OUT / "LICENSE-DeepDNAshape.txt")
print("exported", len(FEATURES), "features", flush=True)

# ------------------------------------------------------------------ 2. same predictions as the original
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.dna import shape  # noqa: E402

shape.SHAPE = OUT
shape.index.cache_clear()
shape.weights.cache_clear()
sys.path.insert(0, str(src))
from deepDNAshape import predictor  # noqa: E402

P = predictor.predictor()
rng = np.random.default_rng(7)
seqs = ["".join(rng.choice(list("ACGT"), n)) for n in (5, 8, 12, 20, 33, 60) for _ in range(4)] + ["CGCGAATTCGCG", "AAAAAAAAAA", "GCGCGCGC"]
worst = {}
for f in FEATURES:
    d = 0.0
    for s in seqs:
        for layer in (2, 4, 7):
            a = np.asarray(P.predict(f, s, layer), dtype=float)
            b = shape.predict(f, s, layer)
            scale = index["params"][f].get("percentile_range") or index["params"][f].get("std") or 1.0
            d = max(d, float(np.max(np.abs(a - b))) / scale)
    worst[f] = d
gate("Athanor's NumPy Deep DNAshape equals the original TensorFlow model (largest difference, in units of each "
     "feature's normalisation range, <= 1e-4)", max(worst.values()) <= 1e-4,
     ", ".join(f"{k} {v:.1e}" for k, v in worst.items()))

# ------------------------------------------------------------------ 3. against the published data
report = {}
try:
    with urllib.request.urlopen("https://zenodo.org/api/records/10403307", timeout=120) as r:
        meta = json.loads(r.read())
    with urllib.request.urlopen(meta["files"][0]["links"]["self"], timeout=600) as r:
        z = zipfile.ZipFile(io.BytesIO(r.read()))
    for kind in ("Expt", "MD"):
        for f in ("MGW", "ProT", "Roll", "HelT"):
            name = next((n for n in z.namelist() if f"/{kind}/{f}_aggre" in n), None)
            if not name:
                continue
            xs, ys = [], []
            for line in io.TextIOWrapper(z.open(name)):
                parts = line.split()
                if len(parts) < 2 or any(b not in "ACGT" for b in parts[0]):
                    continue
                vals = np.array([float(v) for v in parts[1:]])
                pred = shape.predict(f, parts[0], 4)
                if len(pred) != len(vals):
                    continue
                ok = ~np.isnan(vals)
                xs += list(pred[ok]); ys += list(vals[ok])
            if len(xs) > 10:
                rr = float(np.corrcoef(xs, ys)[0, 1])
                report[f"{kind}:{f}"] = {"n": len(xs), "pearson": round(rr, 3),
                                         "mae": round(float(np.mean(np.abs(np.array(xs) - np.array(ys)))), 3)}
                print(f"   {kind} {f}: Pearson r {rr:.3f}, mean |error| {report[f'{kind}:{f}']['mae']} over {len(xs)} positions")
except Exception as e:  # noqa: BLE001
    print("published data not checked:", repr(e)[:200])
(WORK / "dnashape_report.json").write_text(json.dumps({"results": results, "agreement": report}, indent=1))
# the committed summary holds only rounded figures, so it regenerates byte for byte
# (the agreement with the authors' processed data stays in the CI report: the alignment of those files'
# values to positions is not documented well enough to quote the figures)
v = {"equal_to_original": all(r["pass"] for r in results)}
(OUT / "validation.json").write_text(json.dumps(v, indent=1, sort_keys=True) + "\n")
sys.exit(0 if all(r["pass"] for r in results) else 1)
