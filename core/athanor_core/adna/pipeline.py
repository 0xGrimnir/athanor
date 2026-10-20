"""Ancient-DNA pipeline: simulate -> trim/merge -> map -> damage -> consensus -> score.

Runs the published tools from the Linux engines pack (bioconda environment at
/opt/athanor/bio). Each step records the exact command and tool version.

  simulate   gargammel (Renaud, Hanghøj, Willerslev & Orlando, Bioinformatics 33, 577, 2017):
             fragSim fragment lengths, deamSim Briggs-model damage, ART Illumina reads
  merge      AdapterRemoval v3 (Schubert, Lindgreen & Orlando, BMC Res Notes 9, 88, 2016)
  map        BWA aln with ancient-DNA settings -l 1024 -n 0.04 -o 2 (Li & Durbin 2009;
             settings per Schubert et al., BMC Genomics 13, 178, 2012)
  damage     mapDamage2 (Jónsson et al., Bioinformatics 29, 1682, 2013)
  consensus  athanor_core.adna.consensus (damage-aware posteriors)
  score      minimap2 (Li, Bioinformatics 34, 3094, 2018) aligns the truth to the reference
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from .consensus import BASES, calibration, call_consensus
from .damage import bayesian_estimates, damage_profiles
from .pileup import observations_from_bam
from .truth import indel_sites, truth_by_reference_position

BIO = Path(os.environ.get("ATHANOR_BIO", "/opt/athanor/bio"))
REFS = Path(os.environ.get("ATHANOR_REFS", "/opt/athanor/refs"))


class Step:
    def __init__(self, log: list, on_step=None):
        self.log, self.on_step = log, on_step

    def run(self, name: str, cmd: list, cwd: Path, stdout_path: Path | None = None) -> str:
        env = dict(os.environ, PATH=f"{BIO / 'bin'}:{os.environ.get('PATH', '')}")
        t0 = time.time()
        if self.on_step:
            self.on_step({"task": "adna", "stage": name})
        out = open(stdout_path, "wb") if stdout_path else subprocess.PIPE
        try:
            r = subprocess.run([str(c) for c in cmd], cwd=cwd, env=env, stdout=out, stderr=subprocess.PIPE)
        finally:
            if stdout_path:
                out.close()
        self.log.append({"step": name, "command": " ".join(map(str, cmd)), "seconds": round(time.time() - t0, 1),
                         "exit": r.returncode})
        if r.returncode != 0:
            tail = r.stderr.decode(errors="replace")[-1500:]
            raise RuntimeError(f"{name} failed (exit {r.returncode}): {tail}")
        return (r.stdout or b"").decode(errors="replace") if not stdout_path else ""


def tool_versions() -> dict:
    env = dict(os.environ, PATH=f"{BIO / 'bin'}:{os.environ.get('PATH', '')}")

    def first(cmd, pat=None):
        try:
            r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=60)
            text = (r.stdout + r.stderr).strip().splitlines()
            lines = [l for l in text if (pat in l if pat else l.strip())]
            return lines[0].strip() if lines else "unknown"
        except Exception:  # noqa: BLE001
            return "unavailable"
    return {"bwa": first(["bwa"], "Version"), "samtools": first(["samtools", "--version"]),
            "adapterremoval": first(["adapterremoval3", "--version"]), "mapDamage": first(["mapDamage", "--version"]),
            "minimap2": first(["minimap2", "--version"]), "gargammel": "gargammel (bioconda; see pack lock)"}


def reference_catalog() -> dict:
    return json.loads((REFS / "manifest.json").read_text())


def read_fasta(path: Path) -> tuple[str, str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = fh.read().splitlines()
    return lines[0][1:].split()[0], "".join(l.strip() for l in lines[1:] if not l.startswith(">")).upper()


def _map_once(step, rdir: Path, label: str, ref_seq: str, merged: str, stats: bool, min_mapq: int, on_step):
    """Map merged reads to ``ref_seq`` and estimate damage; returns (obs, profile, bayes, flagstat)."""
    rdir.mkdir(parents=True, exist_ok=True)
    ref_fa = rdir / "ref.fa"
    ref_fa.write_text(">ref\n" + "\n".join(ref_seq[i:i + 70] for i in range(0, len(ref_seq), 70)) + "\n")
    step.run(f"index {label} (bwa)", [BIO / "bin/bwa", "index", ref_fa], rdir)
    step.run(f"map to {label} (bwa aln)", [BIO / "bin/bwa", "aln", "-l", 1024, "-n", 0.04, "-o", 2, "-t", 4,
                                           ref_fa, merged], rdir, stdout_path=rdir / "reads.sai")
    step.run(f"map to {label} (bwa samse)", [BIO / "bin/bwa", "samse", ref_fa, rdir / "reads.sai", merged],
             rdir, stdout_path=rdir / "reads.sam")
    step.run(f"sort {label} (samtools)", [BIO / "bin/samtools", "sort", "-o", rdir / "reads.bam", rdir / "reads.sam"], rdir)
    step.run(f"index BAM {label} (samtools)", [BIO / "bin/samtools", "index", rdir / "reads.bam"], rdir)
    (rdir / "reads.sam").unlink()
    flag = step.run(f"flagstat {label}", [BIO / "bin/samtools", "flagstat", rdir / "reads.bam"], rdir)
    md = [BIO / "bin/mapDamage", "-i", rdir / "reads.bam", "-r", ref_fa, "-d", rdir / "mapdamage", "-q"]
    if not stats:
        md.append("--no-stats")
    step.run(f"damage on {label} (mapDamage2)", md, rdir)
    prof = damage_profiles(rdir / "mapdamage")
    bayes = bayesian_estimates(rdir / "mapdamage") if stats else {}
    if on_step:
        on_step({"task": "adna", "stage": f"consensus on {label}"})
    return observations_from_bam(str(rdir / "reads.bam"), min_mapq=min_mapq), prof, bayes, flag


def _reconstruct(step, work: Path, acc: str, cat: dict, merged: str, truth_seq: str, endogenous: str,
                 iterations: int, mapdamage_stats: bool, min_mapq: int, on_step) -> dict:
    base = work / f"ref_{acc}"
    base.mkdir()
    ref_id, ref_seq = read_fasta(REFS / f"{acc}.fa")
    paf = step.run(f"align truth to {acc} (minimap2)", [BIO / "bin/minimap2", "-cx", "asm20", "--cs",
                                                        REFS / f"{acc}.fa", REFS / f"{endogenous}.fa"], base)
    tmap = truth_by_reference_position(paf, len(ref_seq), truth_seq)
    tb = np.array([t if t is not None else "N" for t in tmap])
    known = np.array([t is not None and t in BASES for t in tmap])
    ref_arr = np.array(list(ref_seq))
    current = ref_seq
    history = []
    for it in range(1, max(1, iterations) + 1):
        last = it == max(1, iterations)
        label = acc if iterations == 1 else f"{acc} iteration {it}"
        obs, prof, bayes, flag = _map_once(step, base / f"it{it}", label, current, merged,
                                           mapdamage_stats and last, min_mapq, on_step)
        # calls use the divergence prior relative to the ORIGINAL reference
        cons = call_consensus(ref_seq, obs, prof["delta5"], prof["delta3"], min_depth=1)
        calls = np.array(list(cons.call))
        scorable = known & (cons.depth > 0)
        correct = calls[scorable] == tb[scorable]
        c_it = calibration(cons.confidence[scorable], correct)
        history.append({"iteration": it, "mapped_reads": _flagstat(flag, " mapped ("),
                        "predicted_error": c_it["predicted_error_rate"], "observed_error": c_it["observed_error_rate"],
                        "ece": c_it["expected_calibration_error"],
                        "covered_fraction": float((cons.depth > 0).mean()), "mean_depth": float(cons.depth.mean()),
                        "accuracy": float(correct.mean()) if correct.size else None,
                        "estimated_divergence_theta": cons.theta})
        if not last:
            upd = (cons.confidence >= 0.99) & (cons.depth > 0)
            current = "".join(c if u else r for c, u, r in zip(cons.call, upd, current))
    cal = calibration(cons.confidence[scorable], correct)
    err = np.flatnonzero(scorable)[~correct]
    return {
        "reference": cat[acc], "reference_id": ref_id, "length": len(ref_seq),
        "total_reads": _flagstat(flag, " in total"), "mapped_reads": _flagstat(flag, " mapped ("),
        "flagstat": flag.strip(), "observations": len(obs),
        "mean_depth": float(cons.depth.mean()), "covered_fraction": float((cons.depth > 0).mean()),
        "truth_vs_reference_divergence": float((ref_arr[known] != tb[known]).mean()) if known.any() else None,
        "estimated_divergence_theta": cons.theta,
        "truth_gap_positions": int(sum(t == "-" for t in tmap)),
        "unaligned_positions": int(sum(t is None for t in tmap)),
        "scored_sites": int(scorable.sum()), "accuracy": float(correct.mean()) if correct.size else None,
        "error_positions": err.tolist(),
        "error_diagnostics": _diagnose(err, calls, tb, ref_arr, cons, indel_sites(paf)),
        "iterations": history,
        "calibration": cal, "damage_profile": {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in prof.items()},
        "mapdamage_bayesian": bayes,
        "consensus": cons.call, "confidence": cons.confidence.round(6).tolist(), "depth": cons.depth.tolist(),
    }


def _flagstat(text: str, key: str):
    """First count on the samtools flagstat line containing ``key``."""
    for line in (text or "").splitlines():
        if key in line:
            return int(line.split()[0])
    return None


def _diagnose(err_pos, calls, truth, ref, cons, indels) -> dict:
    """Why are calls wrong? Reference bias shows up as wrong calls equal to the
    reference base where the truth differs from it."""
    if len(err_pos) == 0:
        return {"n": 0}
    ind = np.array(indels) if indels else np.array([-10**9])
    dist = np.array([int(np.min(np.abs(ind - p))) for p in err_pos])
    called_ref = calls[err_pos] == ref[err_pos]
    truth_differs = truth[err_pos] != ref[err_pos]
    return {
        "n": int(len(err_pos)),
        "called_reference_where_truth_differs": int((called_ref & truth_differs).sum()),
        "within_10bp_of_truth_indel": int((dist <= 10).sum()),
        "median_depth_at_errors": float(np.median(cons.depth[err_pos])),
        "median_depth_overall": float(np.median(cons.depth[cons.depth > 0])) if (cons.depth > 0).any() else 0.0,
        "mean_confidence_at_errors": float(cons.confidence[err_pos].mean()),
        "examples": [{"pos": int(p), "call": str(calls[p]), "truth": str(truth[p]), "ref": str(ref[p]),
                      "depth": int(cons.depth[p]), "confidence": float(cons.confidence[p]), "dist_to_indel": int(d)}
                     for p, d in list(zip(err_pos, dist))[:15]],
    }


def run(work: str | Path, endogenous: str = "NC_007596", references: tuple = ("NC_007596", "NC_000934", "NC_005129"),
        coverage: float = 10.0, damage: tuple = (0.03, 0.4, 0.01, 0.3), contamination: float = 0.0,
        length_loc: float = 3.7, length_scale: float = 0.35, min_length: int = 25, read_length: int = 75,
        mapdamage_stats: bool = True, min_mapq: int = 25, iterations: int = 1, on_step=None) -> dict:
    """Simulate ancient reads from ``endogenous`` and reconstruct it against each reference.

    damage = Briggs et al. (2007) parameters (nick frequency, overhang geometric
    parameter, double-stranded deamination, single-stranded deamination).
    Fragment lengths are log-normal with log-location ``length_loc`` and
    log-scale ``length_scale`` (exp(3.7) ≈ 40 bp median).

    iterations > 1 re-maps the reads to a reference updated with the confident
    consensus bases (iterative mapping, as used for the Neanderthal
    mitogenome; Briggs et al., Science 325, 318, 2009), which reduces
    reference bias. Final calls keep the divergence prior relative to the
    original reference, so the reads are not counted twice."""
    work = Path(work)
    shutil.rmtree(work, ignore_errors=True)
    for d in ("input/endo", "input/cont", "input/bact", "sim", "merge"):
        (work / d).mkdir(parents=True, exist_ok=True)
    log: list = []
    step = Step(log, on_step)
    cat = {g["accession"]: g for g in reference_catalog()["genomes"]}
    shutil.copy(REFS / f"{endogenous}.fa", work / "input/endo/endo.fa")
    shutil.copy(REFS / "NC_012920.fa", work / "input/cont/cont.fa")   # present-day human contaminant
    shutil.copy(REFS / "NC_012920.fa", work / "input/bact/bact.fa")   # placeholder; bacterial fraction is 0
    truth_id, truth_seq = read_fasta(REFS / f"{endogenous}.fa")

    # 1. simulate
    step.run("simulate (gargammel)", [
        BIO / "bin/gargammel", "-c", coverage, "--comp", f"0,{contamination},{1 - contamination}",
        "-damagee", ",".join(map(str, damage)), "--loc", length_loc, "--scale", length_scale,
        "--minsize", min_length, "--maxsize", 200, "-rl", read_length, "-ss", "HS25",
        "-o", work / "sim/sim", work / "input"], work)
    s1 = sorted(glob.glob(str(work / "sim/*_s1.fq.gz")))
    s2 = sorted(glob.glob(str(work / "sim/*_s2.fq.gz")))
    if not s1 or not s2:
        raise RuntimeError(f"gargammel produced no paired reads: {os.listdir(work / 'sim')}")

    # 2. trim adapters and merge overlapping pairs into full fragments
    step.run("trim and merge (AdapterRemoval)", [
        BIO / "bin/adapterremoval3", "--in-file1", s1[-1], "--in-file2", s2[-1], "--out-prefix", work / "merge/m",
        "--merge", "--threads", "4"], work)
    merged = next(iter(sorted(glob.glob(str(work / "merge/m*merged*")))), None)
    if merged is None:
        raise RuntimeError(f"AdapterRemoval produced no merged reads: {os.listdir(work / 'merge')}")

    results = []
    for acc in references:
        results.append(_reconstruct(step, work, acc, cat, merged, truth_seq, endogenous, iterations,
                                    mapdamage_stats, min_mapq, on_step))
    return {"simulation": {"endogenous": cat[endogenous], "coverage": coverage,
                           "damage_briggs_v_l_d_s": list(damage), "contamination": contamination,
                           "fragment_length_lognormal": [length_loc, length_scale], "min_length": min_length,
                           "read_length": read_length, "platform": "Illumina HiSeq 2500 profile (ART)",
                           "iterations": iterations},
            "tools": tool_versions(), "steps": log, "references": results,
            "reference_catalog": reference_catalog()}
