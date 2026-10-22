"""The Personal Genome reference data pack, and alignment of a raw file to it.

The pack is built in CI (tools/pgdata/build.py) from public data with
recorded versions and checksums, published to the release feed, and
downloaded once by the app. It holds:

  snps.npz        GRCh37 SNP table of the AADR 1240K set (chromosome,
                  position, reference/alternative alleles, rs number,
                  5 cM jackknife block)
  freqs.npz       alternative-allele counts for the populations used in
                  f-statistics (sources, outgroups, archaic genomes)
  pca_*.npz       principal-component loadings and reference coordinates
  ancients.json.gz  ancient individuals: dates, places, haplogroups,
                  publications, PCA coordinates
  archaic.npz     alleles shared with Neanderthals but nearly absent in Africa
  reference_stats.json  archaic-ancestry statistics of present-day groups
  mt_tree.json.gz, y_tree.json.gz  haplogroup trees
  hg38ToHg19.over.chain.gz  UCSC chain for GRCh38 files
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
from functools import cached_property
from pathlib import Path

import numpy as np

from .rawfile import CODE, MISSING, RawGenotypes

PACK_NAME = "pgdata-1"
COMP = np.full(256, 255, dtype=np.uint8)
COMP[:6] = [3, 2, 1, 0, 4, 5]                                 # A<->T, C<->G; I, D unchanged


def default_dir() -> Path:
    env = os.environ.get("ATHANOR_PGDATA")
    if env:
        return Path(env)
    base = os.environ.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), ".local", "share")
    return Path(base) / "Athanor" / PACK_NAME


class RefData:
    def __init__(self, root: str | Path | None = None):
        self.root = Path(root) if root else default_dir()
        if not (self.root / "manifest.json").exists():
            raise FileNotFoundError(f"Personal Genome reference data not found in {self.root}")
        self.manifest = json.loads((self.root / "manifest.json").read_text())

    def verify(self) -> list[str]:
        """Files whose SHA-256 differs from the manifest."""
        bad = []
        for f in self.manifest["files"]:
            h = hashlib.sha256()
            with open(self.root / f["name"], "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 22), b""):
                    h.update(chunk)
            if h.hexdigest() != f["sha256"]:
                bad.append(f["name"])
        return bad

    @cached_property
    def snps(self) -> dict:
        return dict(np.load(self.root / "snps.npz"))

    @cached_property
    def freqs(self) -> dict:
        d = np.load(self.root / "freqs.npz", allow_pickle=False)
        names = [str(x) for x in d["names"]]
        return {"names": names, "alt": d["alt"], "total": d["total"]}

    def pca(self, name: str) -> dict:
        d = np.load(self.root / f"pca_{name}.npz", allow_pickle=False)
        return {k: d[k] for k in d.files}

    @cached_property
    def ancients(self) -> list[dict]:
        with gzip.open(self.root / "ancients.json.gz", "rt") as fh:
            return json.load(fh)

    @cached_property
    def archaic(self) -> dict:
        d = np.load(self.root / "archaic.npz", allow_pickle=False)
        return {k: d[k] for k in d.files}

    @cached_property
    def reference_stats(self) -> dict:
        return json.loads((self.root / "reference_stats.json").read_text())

    def tree(self, which: str) -> dict:
        with gzip.open(self.root / f"{which}_tree.json.gz", "rt") as fh:
            return json.load(fh)

    def chain(self) -> Path:
        return self.root / "hg38ToHg19.over.chain.gz"


# ------------------------------------------------------------------ alignment
def _keys(chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    return chrom.astype(np.int64) * 1_000_000_000 + pos.astype(np.int64)


def to_grch37(g: RawGenotypes, ref: RefData) -> tuple[np.ndarray, str, dict]:
    """GRCh37 positions for the file's markers (-1 = unmapped) and the build
    it was found to be on, judged by how many positions agree with the
    reference table (rs numbers where present, otherwise position matches)."""
    s = ref.snps
    keys37 = np.sort(_keys(s["chrom"], s["pos"]))
    rs_ok = s["rsid"] > 0
    rs_sorted = np.argsort(s["rsid"][rs_ok])
    rs_vals = s["rsid"][rs_ok][rs_sorted]
    rs_pos = s["pos"][rs_ok][rs_sorted]

    def agreement(pos: np.ndarray) -> float:
        auto = (g.chrom >= 1) & (g.chrom <= 22) & (pos > 0)
        has_rs = auto & (g.rsid > 0)
        if has_rs.sum() > 1000:
            i = np.searchsorted(rs_vals, g.rsid[has_rs])
            i = np.clip(i, 0, len(rs_vals) - 1)
            found = rs_vals[i] == g.rsid[has_rs]
            if found.sum() > 500:
                return float((rs_pos[i][found] == pos[has_rs][found]).mean())
        k = _keys(g.chrom[auto], pos[auto])
        i = np.clip(np.searchsorted(keys37, k), 0, len(keys37) - 1)
        return float((keys37[i] == k).mean())

    a37 = agreement(g.pos)
    info = {"agreement_as_grch37": a37}
    if a37 >= 0.9 or g.build_declared == "37":
        return g.pos.copy(), "GRCh37", info
    from .liftover import Liftover
    lo = Liftover(ref.chain())
    lifted = np.full(len(g.pos), -1, dtype=np.int64)
    from .rawfile import CHROM_NAME
    for c in np.unique(g.chrom):
        m = g.chrom == c
        lifted[m] = lo.map(CHROM_NAME[int(c)], g.pos[m])
    a38 = agreement(lifted)
    info["agreement_after_grch38_liftover"] = a38
    if a38 > a37 and a38 >= 0.8:
        info["unmapped_after_liftover"] = int((lifted < 0).sum())
        return lifted, "GRCh38 (converted to GRCh37)", info
    if a37 >= 0.5:
        return g.pos.copy(), "GRCh37 (partial agreement)", info
    raise ValueError("the file's positions match neither GRCh37 nor GRCh38; it may use an older build (NCBI 36)")


def align(g: RawGenotypes, ref: RefData) -> dict:
    """Alternative-allele dosage (0, 1, 2; -1 missing) of the user at every
    reference SNP. Markers whose alleles match neither the reference strand
    nor its complement are dropped; so are A/T and C/G SNPs, whose strand
    cannot be checked from the alleles."""
    pos37, build, info = to_grch37(g, ref)
    s = ref.snps
    order = np.argsort(_keys(s["chrom"], s["pos"]))
    keys = _keys(s["chrom"], s["pos"])[order]
    auto_or_x = (g.chrom >= 1) & (g.chrom <= 23) & (pos37 > 0) & (g.a1 != MISSING)
    uk = _keys(g.chrom[auto_or_x], pos37[auto_or_x])
    i = np.clip(np.searchsorted(keys, uk), 0, len(keys) - 1)
    hit = keys[i] == uk
    sidx = order[i[hit]]
    a1 = g.a1[auto_or_x][hit]
    a2 = g.a2[auto_or_x][hit]
    r, a = s["ref"][sidx], s["alt"][sidx]
    ambiguous = (COMP[r] == a)        # A/T and C/G: the alleles alone cannot show the strand
    direct = np.isin(a1, [0, 1, 2, 3]) & (((a1 == r) | (a1 == a)) & ((a2 == r) | (a2 == a)))
    c1, c2 = COMP[a1], COMP[a2]
    flipped = ~direct & (((c1 == r) | (c1 == a)) & ((c2 == r) | (c2 == a)))
    # Decide the file's strand from the markers that do show it. The testing
    # companies report the reference's + strand, so almost all of those should
    # match directly; when they do, the ambiguous markers can be used as they
    # are instead of being thrown away (about one marker in six).
    d_n, f_n = int((direct & ~ambiguous).sum()), int((flipped & ~ambiguous).sum())
    decided = d_n + f_n
    same_strand = decided >= 1000 and d_n / decided > 0.99
    other_strand = decided >= 1000 and f_n / decided > 0.99
    amb_ok = np.isin(a1, [0, 1, 2, 3]) & ambiguous & (
        (((a1 == r) | (a1 == a)) & ((a2 == r) | (a2 == a))) if same_strand else
        ((((c1 == r) | (c1 == a)) & ((c2 == r) | (c2 == a))) if other_strand else False))
    use = ((direct | flipped) & ~ambiguous) | amb_ok
    flip_it = np.where(ambiguous, other_strand, flipped)
    b1 = np.where(flip_it, c1, a1)
    b2 = np.where(flip_it, c2, a2)
    dos = np.full(len(s["pos"]), -1, dtype=np.int8)
    dos[sidx[use]] = (b1[use] == a[use]).astype(np.int8) + (b2[use] == a[use]).astype(np.int8)
    strand = ("the reference's + strand" if same_strand else
              "the opposite strand" if other_strand else "mixed or undetermined")
    stats = {"build": build, **info, "markers_on_reference": int(hit.sum()),
             "used": int(use.sum()), "strand": strand,
             "strand_flipped": int((flip_it & use).sum()),
             "ambiguous_used": int(amb_ok.sum()),
             "ambiguous_dropped": int((ambiguous & ~amb_ok & (direct | flipped | ambiguous)).sum()),
             "allele_mismatch_dropped": int((~direct & ~flipped & ~ambiguous).sum())}
    return {"dosage": dos, "positions37": pos37, "stats": stats}


def chromosome_calls(g: RawGenotypes, chrom: int, positions37: np.ndarray | None = None) -> dict[int, str]:
    """Haploid calls on chromosome 24 (Y, GRCh37) or 26 (MT, rCRS): position -> base.
    Heterozygous calls (heteroplasmy, errors) are skipped."""
    m = (g.chrom == chrom) & (g.a1 != MISSING) & (g.a1 == g.a2) & (g.a1 < 4)
    pos = (positions37 if (positions37 is not None and chrom == 24) else g.pos)[m]
    bases = g.a1[m]
    return {int(p): "ACGT"[int(b)] for p, b in zip(pos, bases) if p > 0}
