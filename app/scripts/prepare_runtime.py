"""Assemble app/src-tauri/resources/runtime for the installer.

  runtime/env.tar.gz   conda-pack archive of the science environment (made by CI)
  runtime/env_id.txt   hash of the environment lock; changes trigger re-extraction
  runtime/core/        athanor_core package plus the benchmark reference file

Also writes app/ui/data/benchmark.json from benchmark reports, so the
Validation view shows the numbers measured for this build.
Usage: python prepare_runtime.py --lock env-lock.txt [--bench a.json b.json ...]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNTIME = ROOT / "app" / "src-tauri" / "resources" / "runtime"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lock", required=True, help="explicit environment lock file")
    ap.add_argument("--bench", nargs="*", default=[], help="benchmark.json files to merge")
    a = ap.parse_args()

    RUNTIME.mkdir(parents=True, exist_ok=True)
    lock = Path(a.lock).read_bytes()
    env_id = hashlib.sha256(lock).hexdigest()[:16]
    (RUNTIME / "env_id.txt").write_text(env_id)

    core_dst = RUNTIME / "core"
    shutil.rmtree(core_dst, ignore_errors=True)
    shutil.copytree(ROOT / "core" / "athanor_core", core_dst / "athanor_core",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    (core_dst / "benchmarks" / "reference").mkdir(parents=True)
    shutil.copy2(ROOT / "core" / "benchmarks" / "reference" / "phase1.json",
                 core_dst / "benchmarks" / "reference" / "phase1.json")

    merged = {"tiers": {}, "sources": []}
    for b in a.bench:
        p = Path(b)
        if p.exists():
            rep = json.loads(p.read_text())
            merged["tiers"].update(rep.get("tiers", {}))
            merged["sources"].append(p.parent.name)
    out = ROOT / "app" / "ui" / "data" / "benchmark.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    if merged["tiers"]:
        out.write_text(json.dumps(merged))
    print(f"runtime prepared: env_id={env_id}, tiers in benchmark: {sorted(merged['tiers'])}")


if __name__ == "__main__":
    main()
