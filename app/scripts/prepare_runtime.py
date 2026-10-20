"""Assemble app/src-tauri/resources/runtime for the installer.

  runtime/env.tar.gz   conda-pack archive of the science environment (made by CI)
  runtime/env_id.txt   hash of the environment lock; changes trigger re-extraction
  runtime/core/        athanor_core package plus the benchmark reference file
  runtime/relocate.json  the prefix records whose files actually change on
                       relocation, found by testing every recorded file in the
                       archive here, so the user's machine opens only those

Also writes app/ui/data/benchmark.json from benchmark reports, so the
Validation view shows the numbers measured for this build.
Usage: python prepare_runtime.py --lock env-lock.txt [--bench a.json b.json ...]
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import shutil
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNTIME = ROOT / "app" / "src-tauri" / "resources" / "runtime"


def relocation_list(archive: Path) -> tuple[list, int]:
    """Run conda-pack's own replacement code over every recorded file in the
    archive and keep the records whose contents would change. File contents
    are identical on every machine, so the result holds for all installs."""
    norm = lambda n: n[2:] if n.startswith("./") else n  # noqa: E731
    wanted = {"Scripts/conda-unpack-script.py", "bin/conda-unpack",
              "Scripts/conda_unpack_progress.py", "bin/conda_unpack_progress.py"}
    found: dict = {}
    with tarfile.open(archive, "r:gz") as tf:
        for m in tf:
            if norm(m.name) in wanted and m.isfile():
                found[norm(m.name)] = tf.extractfile(m).read().decode("utf-8")
    script = found.get("Scripts/conda-unpack-script.py") or found.get("bin/conda-unpack")
    progress = found.get("Scripts/conda_unpack_progress.py") or found.get("bin/conda_unpack_progress.py") or ""
    if script is None:
        raise RuntimeError("no conda-unpack script in the archive")
    import sys
    import types
    sys.modules["conda_unpack_progress"] = mod = types.ModuleType("conda_unpack_progress")
    exec(compile(progress, "conda_unpack_progress.py", "exec"), mod.__dict__)
    if script.startswith("#!"):
        script = script.split("\n", 1)[1]
    ns: dict = {"__name__": "athanor_relocation_scan"}
    exec(compile(script, "conda-unpack-script.py", "exec"), ns)
    replace_prefix, records = ns["replace_prefix"], ns["_prefix_records"]
    by_path = {r[0].replace("\\", "/"): r for r in records}
    probe = "C:/__athanor_relocation_probe__"
    keep = []
    with tarfile.open(archive, "r:gz") as tf:
        for m in tf:
            rec = by_path.get(norm(m.name))
            if rec is None or not m.isfile():
                continue
            data = tf.extractfile(m).read()
            if replace_prefix(data, rec[2], rec[1], probe) != data:
                keep.append(list(rec))
    return keep, len(records)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lock", required=True, help="explicit environment lock file")
    ap.add_argument("--bench", nargs="*", default=[], help="benchmark.json files to merge")
    a = ap.parse_args()

    RUNTIME.mkdir(parents=True, exist_ok=True)
    lock = Path(a.lock).read_bytes()
    unpacker = (ROOT / "app" / "runtime" / "unpack_env.py").read_bytes()
    shutil.copy2(ROOT / "app" / "runtime" / "unpack_env.py", RUNTIME / "unpack_env.py")
    # Without a trustworthy list the app falls back to scanning every record.
    reloc = b""
    (RUNTIME / "relocate.json").unlink(missing_ok=True)
    try:
        keep, total = relocation_list(RUNTIME / "env.tar.gz")
        if total == 0:
            raise RuntimeError("archive has no prefix records")
        reloc = json.dumps(keep).encode()
        (RUNTIME / "relocate.json").write_bytes(reloc)
        print(f"relocation: {len(keep)} of {total} recorded files change")
    except Exception as e:  # noqa: BLE001
        print(f"relocation list not built ({e}); the app will scan every record")
    # A new environment or a new relocation step both force a fresh extraction.
    env_id = hashlib.sha256(lock + b"\0" + unpacker + b"\0" + reloc).hexdigest()[:16]
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
