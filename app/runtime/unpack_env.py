"""Relocate an extracted conda-pack environment (replaces ``conda-unpack``).

conda-pack's own ``conda-unpack`` opens every recorded file read-write, even
those it leaves unchanged. On Windows that fails for DLLs (.pyd) its own Python
process has loaded, e.g. ``DLLs/_bz2.pyd``. This runner reuses conda-pack's
prefix records and replacement functions from ``conda-unpack-script.py`` but
reads each file first and rewrites only those whose contents change. Run it
with the environment's own interpreter: ``python -I -S unpack_env.py <env>``.

Prefix-replacement logic is conda-pack's (BSD-3-Clause), executed unmodified.
"""
from __future__ import annotations

import os
import sys


def relocate(env: str) -> tuple[int, int]:
    env = os.path.abspath(env)
    win = sys.platform == "win32"
    scripts = os.path.join(env, "Scripts" if win else "bin")
    script = os.path.join(scripts, "conda-unpack-script.py" if win else "conda-unpack")
    if not os.path.exists(script):
        return 0, 0
    sys.path.insert(0, scripts)  # conda_unpack_progress lives beside the script
    ns: dict = {"__name__": "athanor_conda_unpack", "__file__": script}
    with open(script, encoding="utf-8") as fh:
        src = fh.read()
    if src.startswith("#!"):
        src = src.split("\n", 1)[1]
    exec(compile(src, script, "exec"), ns)
    replace_prefix = ns["replace_prefix"]
    records = ns["_prefix_records"]

    changed = 0
    for path, placeholder, mode in records:
        full = os.path.join(env, path)
        if win:
            full = full.replace("\\", "/")
        new_prefix = env.replace("\\", "/") if (win and mode == "text") else env
        with open(full, "rb") as fh:
            data = fh.read()
        new = replace_prefix(data, mode, placeholder, new_prefix)
        if new != data:
            with open(full, "wb") as fh:
                fh.write(new)
            changed += 1
    return len(records), changed


if __name__ == "__main__":
    n, changed = relocate(sys.argv[1])
    print(f"relocated {changed} of {n} recorded files")
