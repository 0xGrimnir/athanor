"""Present-day human contamination of non-human ancient DNA by competitive mapping.

Method (Feuerborn et al., "Competitive mapping allows for the identification
and exclusion of human DNA contamination in ancient faunal genomic
datasets", BMC Genomics 21, 844, 2020): map all reads with BWA aln to a
concatenation of the target species' reference and the human reference;
keep reads with mapping quality above 30; the contamination estimate is
the fraction of those reads assigned to the human sequence. Reads that fit
both equally well get mapping quality 0 and are discarded, so the estimate
is a lower bound. Here the references are the mitochondrial genomes (the
target species' and the revised Cambridge Reference Sequence).

As in that study, the human-assigned reads are also compared with the
target reads for read length and terminal deamination: genuine
present-day contamination is longer and lacks the damage signal.

The interval is the exact binomial (Clopper-Pearson) 95% interval.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .bamtools import clean_bam, write_fasta
from .consensus import CODE
from .damage import damage_from_observations, terminal_signal
from .pileup import observations_from_bam

MIN_MAPQ = 30       # reads must exceed this (Feuerborn et al. 2020)


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    from scipy.stats import beta
    if n == 0:
        return (0.0, 1.0)
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return (lo, hi)


def competitive_mapping(step, work: Path, fastq: Path, target: tuple[str, str], human: tuple[str, str],
                        threads: int = 4, prefiltered: bool = False) -> dict:
    """``target`` and ``human`` are (accession, sequence)."""
    from .pipeline import BIO
    d = work / "contamination"
    d.mkdir(parents=True, exist_ok=True)
    ref = write_fasta(d / "competitive.fa", {target[0]: target[1], human[0]: human[1]})
    step.run("index competitive reference (bwa)", [BIO / "bin/bwa", "index", ref], d)
    step.run("competitive mapping (bwa aln)", [BIO / "bin/bwa", "aln", "-l", 1024, "-n", 0.04, "-o", 2,
                                               "-t", threads, ref, fastq], d, stdout_path=d / "reads.sai")
    step.run("competitive mapping (bwa samse)", [BIO / "bin/bwa", "samse", ref, d / "reads.sai", fastq],
             d, stdout_path=d / "reads.sam")
    step.run("sort competitive BAM (samtools)", [BIO / "bin/samtools", "sort", "-o", d / "raw.bam", d / "reads.sam"], d)
    (d / "reads.sam").unlink()
    clean = step.python("assign reads (MAPQ > 30, duplicates removed)",
                        lambda: clean_bam(d / "raw.bam", d / "assigned.bam", MIN_MAPQ, strict_mapq=True))

    def group(acc: str, seq: str) -> dict:
        obs = observations_from_bam(str(d / "assigned.bam"), contig=acc, min_mapq=0)
        reads = _reads_on(d / "assigned.bam", acc)
        codes = np.array([CODE.get(b, -1) for b in seq])
        prof = damage_from_observations(obs, codes)
        lengths = np.array([ln for ln, partial in reads if not partial])
        return {"reads": len(reads), "whole_molecules": int(lengths.size),
                "median_length": float(np.median(lengths)) if lengths.size else None,
                "damage": terminal_signal(prof)}

    groups = step.python("contamination statistics", lambda: {"target": group(*target), "human": group(*human)})
    k, n = groups["human"]["reads"], groups["human"]["reads"] + groups["target"]["reads"]
    lo, hi = clopper_pearson(k, n)
    notes = ["Lower bound: reads matching both genomes equally well are discarded (Feuerborn et al. 2020)."]
    if prefiltered:
        notes.append("The input BAM held only reads already mapped to a reference before deposit, so most "
                     "contaminant reads were removed beforehand; this estimate cannot detect them.")
    return {"method": "competitive mapping (Feuerborn et al., BMC Genomics 21, 844, 2020)",
            "references": {"target": target[0], "human": human[0]}, "min_mapq_exclusive": MIN_MAPQ,
            "mapping": clean, "human_reads": k, "assigned_reads": n,
            "fraction": k / n if n else None, "ci95": (lo, hi),
            "target_reads": groups["target"], "human_reads_detail": groups["human"], "notes": notes}


def _reads_on(bam: Path, contig: str) -> list[tuple[int, bool]]:
    import pysam
    from .pileup import is_partial
    with pysam.AlignmentFile(str(bam), "rb") as b:
        return [(r.query_length, is_partial(r.query_name)) for r in b.fetch(contig)]
