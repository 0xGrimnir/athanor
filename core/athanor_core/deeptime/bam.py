"""Reading an aligned ancient genome: which reference it is on, how much of it
there is, how long the molecules are, and which sex chromosomes are present.

Reference build. A BAM says only the name and length of each sequence it was
mapped to. Chromosome lengths differ between GRCh37 and GRCh38 by millions of
bases, so the build is read off the lengths, and the naming style ("1" or
"chr1") off the names. hg19's mitochondrial sequence (chrM, 16,571 bp) is not
the rCRS (16,569 bp) used by the mtDNA tree, and its coordinates differ past
position 3106, so mitochondrial calls are refused on that reference rather
than silently misplaced.

Coverage and sex. Read counts per chromosome come from the BAM index
(samtools idxstats), so no reads have to be read. Genetic sex follows
Skoglund et al. (Nat. Commun. 4, 2013) R_y = n_Y / (n_X + n_Y) with their
thresholds, and Mittnik et al. (Sci. Rep. 6, 31326, 2016) R_x, the X-to-
autosome coverage ratio (about 1 with two X chromosomes, about 0.5 with one).
Both are reported with binomial confidence intervals, and a sample is called
only when the interval falls wholly inside one range.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np

AUTOSOMES = tuple(str(i) for i in range(1, 23))

GRCH37 = {"1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260, "6": 171115067,
          "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
          "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753, "17": 81195210, "18": 78077248,
          "19": 59128983, "20": 63025520, "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566}
GRCH38 = {"1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979,
          "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
          "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285,
          "19": 58617616, "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415}
RCRS = 16569          # NC_012920, the mtDNA tree's reference
HG19_CHRM = 16571     # NC_001807, hg19's chrM: different coordinates past 3106

# Skoglund et al. 2013, table 1
RY_FEMALE, RY_MALE = 0.016, 0.077


def plain(name: str) -> str:
    """Chromosome name without the build's decoration: chr1 -> 1, chrM/chrMT -> MT."""
    n = re.sub(r"^chr", "", str(name), flags=re.I)
    return "MT" if n.upper() in ("M", "MT") else n.upper() if n.upper() in ("X", "Y") else n


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if n <= 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


class Genome:
    """An aligned ancient genome: the BAM, its reference, and its index counts."""

    def __init__(self, path: str | Path, min_mapq: int = 25):
        import pysam
        self.path = str(path)
        self.min_mapq = int(min_mapq)
        self.bam = pysam.AlignmentFile(self.path, "rb")
        self.lengths = {plain(n): int(L) for n, L in zip(self.bam.references, self.bam.lengths)}
        self.names = {plain(n): n for n in self.bam.references}
        self.build, self.build_detail = self._build()
        self.indexed = self.bam.has_index()

    def close(self):
        self.bam.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def _build(self) -> tuple[str, dict]:
        best, detail = None, {}
        for name, ref in (("GRCh37", GRCH37), ("GRCh38", GRCH38)):
            shared = [c for c in ref if c in self.lengths]
            agree = [c for c in shared if self.lengths[c] == ref[c]]
            detail[name] = {"chromosomes_shared": len(shared), "lengths_agree": len(agree)}
            if len(shared) >= 20 and len(agree) == len(shared):
                best = name
        mt = self.lengths.get("MT")
        style = "chr-prefixed" if any(str(n).lower().startswith("chr") for n in self.bam.references) else "plain"
        d = {"naming": style, "mitochondrion": ("rCRS" if mt == RCRS else "hg19 chrM" if mt == HG19_CHRM else
                                                None if mt is None else f"{mt} bp, unrecognised"), **detail}
        return (best or "unrecognised"), d

    @property
    def mt_usable(self) -> bool:
        return self.lengths.get("MT") == RCRS

    def name_of(self, chrom: str) -> str | None:
        return self.names.get(plain(chrom))

    def counts(self) -> dict:
        """Mapped reads per chromosome, from the index."""
        out = {}
        for row in self.bam.get_index_statistics():
            out[plain(row.contig)] = int(row.mapped)
        return out

    def read_groups(self) -> list[dict]:
        rgs = self.bam.header.to_dict().get("RG", [])
        return [{k: v for k, v in rg.items() if k in ("ID", "LB", "SM", "PL", "DS")} for rg in rgs]

    def programs(self) -> list[str]:
        return [pg.get("PN") or pg.get("ID", "") for pg in self.bam.header.to_dict().get("PG", [])]

    def _sex_counts(self, max_reads: int = 2_000_000) -> tuple[dict, bool]:
        """Reads on each chromosome that pass the same filters the rest of the
        module uses. The index's own counts include duplicates, secondary and
        supplementary alignments and multi-mapping reads at mapping quality 0,
        and the published R_y thresholds were set on filtered unique
        alignments, so the index cannot answer this question."""
        want = [*AUTOSOMES, "X", "Y"]
        out, truncated, seen = {}, False, 0
        for chrom in want:
            ref = self.name_of(chrom)
            if ref is None:
                continue
            n = 0
            for r in self.bam.fetch(ref):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_duplicate \
                        or r.mapping_quality < self.min_mapq:
                    continue
                n += 1
                seen += 1
                if seen >= max_reads:
                    truncated = True
                    break
            out[chrom] = n
            if truncated:
                break
        return out, truncated

    def sex(self, filtered: bool = False) -> dict:
        """Genetic sex from the numbers of reads on X and Y (R_y) and from the
        X-to-autosome coverage ratio (R_x).

        These counts are whole chromosomes, so they answer only for shotgun
        sequencing. Capture data has as many reads on a chromosome as the
        experiment put targets there, and both ratios shift; ``sex_from_panel``
        uses the panel's own positions instead and holds either way. The
        thresholds were set on filtered unique alignments, so a call is made
        only from counts that pass the same filters (``filtered=True``);
        otherwise the ratios are reported without one."""
        c, truncated = (self._sex_counts() if filtered else (self.counts(), False))
        nx, ny = c.get("X", 0), c.get("Y", 0)
        auto = {k: c.get(k, 0) for k in AUTOSOMES}
        res = {"reads_X": nx, "reads_Y": ny, "reads_autosomal": int(sum(auto.values())),
               "counts_filtered": bool(filtered), "counts_truncated": truncated,
               "basis": ("reads over whole chromosomes passing the mapping-quality and flag filters" if filtered
                         else "read counts from the index, which include duplicates and multi-mapping reads")}
        if nx + ny < 100:
            return {**res, "call": None, "r_y": None, "r_x": None,
                    "reason": f"only {nx + ny} reads on the sex chromosomes"}
        ry = ny / (nx + ny)
        lo, hi = wilson(ny, nx + ny)
        call = ("XX" if hi < RY_FEMALE else "XY" if lo > RY_MALE else None) if filtered else None
        lx = self.lengths.get("X")
        rx = rx_ci = None
        if lx and res["reads_autosomal"] > 0:
            rx, rx_ci = _ratio_of_rates(nx, lx, [(auto[k], self.lengths[k]) for k in AUTOSOMES
                                                 if auto.get(k, 0) > 0 and self.lengths.get(k)])
        reason = None
        if not filtered:
            reason = ("not called from these counts: they come from the index, which includes duplicates and "
                      "reads of mapping quality 0, while the thresholds were set on filtered unique alignments")
        elif not call:
            reason = (f"R_y = {ry:.3f}, interval {lo:.3f}-{hi:.3f}, spans the thresholds "
                      f"{RY_FEMALE} (XX) and {RY_MALE} (XY)")
        return {**res, "r_y": ry, "r_y_ci": [lo, hi], "r_x": rx, "r_x_ci": rx_ci, "call": call, "reason": reason,
                "method": f"R_y = n_Y/(n_X+n_Y) (Skoglund et al. 2013): below {RY_FEMALE} XX, above {RY_MALE} XY; "
                          "R_x is the X-to-autosome coverage ratio (Mittnik et al. 2016), about 1 with two X "
                          "chromosomes and 0.5 with one. Whole chromosomes include regions where the X maps less "
                          "well than the autosomes, so R_x reads low here; the panel-position R_x does not have "
                          "that problem. The thresholds hold for shotgun sequencing only."}


    def lengths_and_damage_pass(self, chroms: list[str] | None = None, max_reads: int = 2_000_000,
                                max_length: int = 200) -> dict:
        """One pass over reads: length distribution, strand balance, mapping
        qualities. Fragment length is the read's own length for collapsed
        single-end ancient reads; for pairs it is the template length."""
        hist = np.zeros(max_length + 2, dtype=np.int64)
        tlen = np.zeros(max_length + 2, dtype=np.int64)
        n = nrev = npaired = nq = 0
        mapqs = np.zeros(256, dtype=np.int64)
        chroms = chroms or [c for c in AUTOSOMES if c in self.names]
        for chrom in chroms:
            ref = self.name_of(chrom)
            if ref is None:
                continue
            for r in self.bam.fetch(ref):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                    continue
                n += 1
                mapqs[min(255, r.mapping_quality)] += 1
                if r.mapping_quality < self.min_mapq:
                    continue
                nq += 1
                nrev += bool(r.is_reverse)
                hist[min(max_length + 1, r.query_length or 0)] += 1
                if r.is_paired:
                    npaired += 1
                    t = abs(r.template_length or 0)
                    if t:
                        tlen[min(max_length + 1, t)] += 1
                if n >= max_reads:
                    break
            if n >= max_reads:
                break
        return {"reads_seen": n, "reads_passing_mapq": nq, "reverse_fraction": (nrev / nq) if nq else None,
                "paired_fraction": (npaired / nq) if nq else None, "length_histogram": hist.tolist(),
                "template_length_histogram": tlen.tolist() if tlen.sum() else None,
                "mapq_histogram": mapqs.tolist(), "max_length_binned": max_length,
                "chromosomes": chroms, "truncated": n >= max_reads}


# Mittnik et al. 2016: the 95% interval of R_x decides, and nothing is called when it spans both
RX_MALE, RX_FEMALE = 0.6, 0.8
MIN_X_OBSERVATIONS = 400      # below this the counting noise alone is over 5%


def _ratio_of_rates(x_count: float, x_denom: float, autosomes: list[tuple[float, float]],
                    z: float = 1.96) -> tuple[float | None, list[float] | None]:
    """R_x and an interval that carries both sources of uncertainty.

    R_x is the X's reads per unit against each autosome's, averaged. The X's
    own count appears in every one of those ratios, so the spread between them
    says nothing about how well that count is known: at low coverage its
    Poisson noise is the whole uncertainty, and an interval built from the
    spread alone misses the truth most of the time. The two are combined here:
    the relative variance of the X count and of the pooled autosomal count
    (Poisson), plus the variance between the per-autosome ratios."""
    autosomes = [(n, L) for n, L in autosomes if n > 0 and L]
    if not autosomes or x_count <= 0 or x_denom <= 0:
        return None, None
    ratios = [(x_count / x_denom) / (n / L) for n, L in autosomes]
    rx = float(np.mean(ratios))
    n_auto = float(sum(n for n, _ in autosomes))
    var_rel = 1.0 / x_count + 1.0 / n_auto                      # counting noise, both sides
    if len(ratios) > 1:                                          # how much the autosomes disagree
        var_rel += float(np.var(ratios, ddof=1)) / (len(ratios) * rx * rx)
    half = z * rx * math.sqrt(var_rel)
    return rx, [round(max(0.0, rx - half), 4), round(rx + half, 4)]


def sex_from_panel(depth_by_chrom: dict, sites_by_chrom: dict) -> dict:
    """Genetic sex from the depth at the reference panel's own positions.

    The panel is a fixed set of positions, so the comparison is between the
    same targets on the X and on each autosome whether the reads come from
    shotgun sequencing or from a capture experiment aimed at those very
    positions. R_x is the X's depth per site against each autosome's, averaged
    over the autosomes (Mittnik et al., Sci. Rep. 6, 31326, 2016): about 1 with
    two X chromosomes, about 0.5 with one."""
    autos = [c for c in AUTOSOMES if sites_by_chrom.get(c) and depth_by_chrom.get(c)]
    x_sites, x_depth = sites_by_chrom.get("X", 0), depth_by_chrom.get("X", 0)
    res = {"reads_at_X_sites": int(x_depth), "X_sites": int(x_sites),
           "reads_at_autosomal_sites": int(sum(depth_by_chrom.get(c, 0) for c in autos)),
           "autosomal_sites": int(sum(sites_by_chrom.get(c, 0) for c in autos)),
           "basis": "depth at the reference panel's positions"}
    if not autos or x_sites < 100 or x_depth < MIN_X_OBSERVATIONS:
        return {**res, "r_x": None, "call": None,
                "reason": f"too little coverage at the panel's X positions ({int(x_depth)} observations; at least "
                          f"{MIN_X_OBSERVATIONS} are needed before the ratio means anything)"}
    rx, ci = _ratio_of_rates(x_depth, x_sites, [(depth_by_chrom[c], sites_by_chrom[c]) for c in autos])
    lo, hi = ci
    call = "XX" if lo > RX_FEMALE else "XY" if hi < RX_MALE else None
    return {**res, "r_x": rx, "r_x_ci": [lo, hi], "autosomes_used": len(autos), "call": call,
            "reason": None if call else f"R_x = {rx:.2f}, interval {lo:.2f}-{hi:.2f}, which does not sit wholly "
                                        f"below {RX_MALE} (XY) or above {RX_FEMALE} (XX)",
            "method": "R_x, the depth per panel position on the X against the same on each autosome, averaged over "
                      f"the autosomes (Mittnik et al. 2016): wholly below {RX_MALE} is XY, wholly above "
                      f"{RX_FEMALE} is XX. The 95% interval carries the counting noise in the X's own observations "
                      "and in the autosomes', as well as the disagreement between the autosomes. Measured at the "
                      "panel's positions, so capture data and shotgun data are read the same way."}


def length_summary(hist: list[int] | np.ndarray, max_length: int = 200) -> dict:
    h = np.asarray(hist, dtype=np.float64)[: max_length + 1]
    x = np.arange(len(h))
    tot = h.sum()
    if tot == 0:
        return {"reads": 0}
    cdf = np.cumsum(h) / tot
    q = {p: float(np.searchsorted(cdf, p)) for p in (0.25, 0.5, 0.75)}
    return {"reads": int(tot), "mean": float((x * h).sum() / tot), "median": q[0.5],
            "quartiles": [q[0.25], q[0.75]], "shortest_5pc": float(np.searchsorted(cdf, 0.05)),
            "longest_5pc": float(np.searchsorted(cdf, 0.95)),
            "at_or_above_binned_maximum": int(np.asarray(hist)[max_length + 1]) if len(np.asarray(hist)) > max_length + 1 else 0}
