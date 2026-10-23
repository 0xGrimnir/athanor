"""Reading the AADR 1240K release (EIGENSTRAT .snp/.ind, TGENO .geno, .anno).

Shared by the reference-pack build scripts. Genotypes are stored as copies of
allele 1 (the reference allele in the .snp file), 3 = missing; this module
returns alternative-allele dosages (2 - g) with -1 for missing.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
LUT = np.array([[(b >> (6 - 2 * k)) & 3 for k in range(4)] for b in range(256)], dtype=np.uint8)


class AADR:
    def __init__(self, prefix: str | Path, anno: str | Path):
        prefix = Path(prefix)
        lines = prefix.with_suffix(".snp").read_text().splitlines()
        n = len(lines)
        self.chrom = np.zeros(n, np.int16)
        self.pos = np.zeros(n, np.int64)
        self.cm = np.zeros(n)
        self.a1 = np.full(n, 255, np.uint8)
        self.a2 = np.full(n, 255, np.uint8)
        self.rsid = np.full(n, -1, np.int64)
        for i, l in enumerate(lines):
            f = l.split()
            self.chrom[i] = int(f[1]); self.cm[i] = float(f[2]) * 100; self.pos[i] = int(f[3])
            self.a1[i] = BASE.get(f[4], 255); self.a2[i] = BASE.get(f[5], 255)
            if f[0].startswith("rs") and f[0][2:].isdigit():
                self.rsid[i] = int(f[0][2:])
        self.nsnp = n
        ind = [l.split() for l in prefix.with_suffix(".ind").read_text().splitlines()]
        self.ids = [x[0] for x in ind]
        self.group = np.array([x[2] for x in ind])
        self.nind = len(ind)
        rows = list(csv.reader(open(anno, encoding="utf-8", errors="replace"), delimiter="\t"))
        hdr = [h.lower() for h in rows[0]]

        def col(pred):
            return next(i for i, h in enumerate(hdr) if pred(h))

        self.C = {"iid": 2, "date": col(lambda h: h.startswith("date mean in bp")),
                  "snps": col(lambda h: "1240k snpset" in h), "sex": col(lambda h: h.startswith("molecular sex")),
                  "country": col(lambda h: h.startswith("political entity")),
                  "lat": col(lambda h: h == "latitude"), "lon": col(lambda h: h == "longitude")}
        self.anno = {r[0]: r for r in rows[1:] if r}
        self.date_bp = np.array([self._num(self.get(g, "date")) if self.get(g, "date") is not None else np.nan
                                 for g in self.ids])
        self.snps_hit = np.array([self._num(self.get(g, "snps")) or 0 for g in self.ids])
        geno = prefix.with_suffix(".geno")
        with open(geno, "rb") as fh:
            head = fh.read(48).split(b"\x00")[0].decode().split()
        assert head[0] == "TGENO" and int(head[1]) == self.nind and int(head[2]) == n, head
        self.rlen = (n + 3) // 4
        assert geno.stat().st_size == 48 + self.nind * self.rlen, "unexpected TGENO layout"
        self.mm = np.memmap(geno, dtype=np.uint8, mode="r", offset=48, shape=(self.nind, self.rlen))

    @staticmethod
    def _num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    def get(self, gid, key):
        r = self.anno.get(gid)
        if not r or self.C[key] >= len(r):
            return None
        v = r[self.C[key]].strip()
        return None if v in ("", "..", "n/a") else v

    def sex(self, i: int) -> str:
        return (self.get(self.ids[i], "sex") or "U").upper()[:1]

    def rows(self, idx: list[int]) -> np.ndarray:
        """Packed genotype rows (one contiguous read per individual)."""
        out = np.empty((len(idx), self.rlen), np.uint8)
        for k, i in enumerate(idx):
            out[k] = self.mm[i]
        return out

    def dosage(self, packed: np.ndarray, snp_idx: np.ndarray) -> np.ndarray:
        """(rows, len(snp_idx)) alt-allele dosage, -1 missing, from packed rows."""
        out = np.empty((packed.shape[0], len(snp_idx)), np.int8)
        step = 200_000
        for c0 in range(0, len(snp_idx), step):
            s = snp_idx[c0:c0 + step]
            b = packed[:, s // 4]
            g = (b >> (6 - 2 * (s % 4)).astype(np.uint8)) & 3
            out[:, c0:c0 + len(s)] = np.where(g == 3, -1, 2 - g.astype(np.int8))
        return out

    def best_per_individual(self, idx: list[int]) -> list[int]:
        best = {}
        for i in idx:
            iid = self.get(self.ids[i], "iid") or self.ids[i]
            if iid not in best or self.snps_hit[i] > self.snps_hit[best[iid]]:
                best[iid] = i
        return sorted(best.values())
