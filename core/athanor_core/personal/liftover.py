"""Coordinate conversion between genome builds with UCSC chain files.

A chain lists aligned blocks between two assemblies; a position inside a block
maps by a constant offset, a position in a gap does not map. Only chains to
the same chromosome on the + strand are used, which covers the markers on
consumer genotyping chips; anything else is reported as unmapped.
Chain format: https://genome.ucsc.edu/goldenPath/help/chain.html
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np


class Liftover:
    def __init__(self, chain_path: str | Path):
        opener = gzip.open if str(chain_path).endswith(".gz") else open
        blocks: dict = {}
        with opener(chain_path, "rt") as fh:
            cur = None
            for line in fh:
                f = line.split()
                if not f:
                    cur = None
                    continue
                if f[0] == "chain":
                    t_name, t_strand, t_start = f[2], f[4], int(f[5])
                    q_name, q_strand, q_start = f[7], f[9], int(f[10])
                    usable = t_strand == "+" and q_strand == "+" and _norm(t_name) == _norm(q_name)
                    cur = [_norm(t_name), t_start, q_start] if usable else None
                    continue
                if cur is None:
                    continue
                size = int(f[0])
                chrom, t, q = cur
                blocks.setdefault(chrom, []).append((t, t + size, q - t))
                if len(f) == 3:
                    cur[1] = t + size + int(f[1])
                    cur[2] = q + size + int(f[2])
        self.blocks = {}
        for c, b in blocks.items():
            a = np.array(sorted(b), dtype=np.int64)
            self.blocks[c] = (a[:, 0], a[:, 1], a[:, 2])

    def map(self, chrom: str, pos1: np.ndarray) -> np.ndarray:
        """1-based positions on ``chrom`` -> 1-based positions in the target
        build, -1 where the position does not map."""
        out = np.full(len(pos1), -1, dtype=np.int64)
        b = self.blocks.get(_norm(chrom))
        if b is None:
            return out
        starts, ends, off = b
        p0 = np.asarray(pos1, dtype=np.int64) - 1
        i = np.searchsorted(starts, p0, side="right") - 1
        ok = (i >= 0)
        ok[ok] = p0[ok] < ends[i[ok]]
        out[ok] = p0[ok] + off[i[ok]] + 1
        return out


def _norm(name: str) -> str:
    n = name.lower().replace("chr", "")
    return {"m": "mt"}.get(n, n)
