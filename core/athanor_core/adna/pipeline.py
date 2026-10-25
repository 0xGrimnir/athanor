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
import math
import gzip
import json
import os
import re
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
GARGAMMEL = Path(os.environ.get("ATHANOR_GARGAMMEL", "/opt/athanor/gargammel"))


# Expected seconds per stage on a reference machine (CI, 10x mammoth demo);
# the progress model rescales them by how fast this machine actually runs.
# (measured in CI: every tool except mapDamage2's Bayesian estimate takes well
# under a second on mitogenome data; the full MCMC takes about 100 s).
EXPECTED = {"simulate": 1.0, "merge": 0.5, "truth": 0.5, "index": 0.3, "map": 0.3, "sort": 0.3, "flagstat": 0.3,
            "damage_fast": 2.0, "damage_stats": 100.0, "consensus": 0.5,
            # real-data stages
            "read_input": 1.0, "trim": 0.5, "clean": 0.5, "count": 0.5, "compare": 0.3}


class Progress:
    """Continuous percentage and ETA over a known plan of stages.

    Each stage has an expected duration; within a running stage, progress
    advances with elapsed time (capped at 95% of the stage until it finishes).
    The ETA scales the remaining expected time by the ratio of actual to
    expected time over the stages completed so far."""

    def __init__(self, plan: list[tuple[str, float]], on_step=None):
        self.plan, self.on_step = plan, on_step
        self.total = sum(e for _, e in plan) or 1.0
        self.idx, self.exp_done, self.act_done = 0, 0.0, 0.0
        self.kind_exp: dict = {}
        self.kind_act: dict = {}
        self.t_stage = time.time()
        self.t0 = time.time()
        self.last_fraction = 0.0

    def expect(self, remaining: list[tuple[str, float]]):
        """Replace the stages still to come (from the current one on) once
        the work is known better, e.g. after the inputs have been read."""
        self.plan = self.plan[:self.idx] + list(remaining)
        self.total = sum(e for _, e in self.plan) or 1.0

    def _ratio(self, kind: str | None = None) -> float:
        """Actual/expected speed. Per kind of stage once one has completed
        (BWA speed says nothing about R's MCMC speed), otherwise over all
        completed stages; shrunk toward 1 so one odd stage cannot swing it."""
        if kind is not None and self.kind_exp.get(kind, 0) > 0:
            e, a = self.kind_exp[kind], self.kind_act[kind]
            k = 0.2 * e
        else:
            e, a = self.exp_done, self.act_done
            k = min(30.0, 0.1 * self.total)
        return max(0.25, min(8.0, (a + k) / (e + k)))

    @staticmethod
    def _stage_fraction(el: float, expected: float) -> float:
        """Linear to 80% at 80% of the expected time, then an exponential tail:
        always increasing, never reaching 100% before the stage really ends."""
        if expected <= 0:
            return 0.0
        knee = 0.8 * expected
        if el <= knee:
            return el / expected
        return 0.8 + 0.2 * (1 - math.exp(-(el - knee) / (0.5 * expected)))

    def emit(self, label: str):
        if not self.on_step:
            return
        kind_cur = self.plan[self.idx][0] if self.idx < len(self.plan) else None
        exp_cur = self.plan[self.idx][1] if self.idx < len(self.plan) else 0.0
        ratio = self._ratio(kind_cur)
        el = time.time() - self.t_stage
        frac_cur = self._stage_fraction(el, ratio * exp_cur)
        done = self.exp_done + frac_cur * exp_cur
        over = bool(exp_cur) and el > ratio * exp_cur
        # remaining: the rest of the current stage (at least a little while it runs over) plus later stages
        rest_cur = max(ratio * exp_cur - el, 0.15 * ratio * exp_cur if over else 0.0)
        later = sum(self._ratio(k) * e for k, e in self.plan[self.idx + 1:])
        remaining = rest_cur + later
        # never move backwards, even when the plan is revised
        self.last_fraction = max(self.last_fraction, min(0.999, done / self.total))
        self.on_step({"task": "adna", "stage": label, "step": self.idx + 1, "steps": len(self.plan),
                      "fraction": round(self.last_fraction, 4),
                      "eta_seconds": round(max(0.0, remaining)), "elapsed_seconds": round(time.time() - self.t0),
                      "overrunning": bool(over)})

    def finish_stage(self):
        if self.idx < len(self.plan):
            kind, exp = self.plan[self.idx]
            act = time.time() - self.t_stage
            self.exp_done += exp
            self.act_done += act
            self.kind_exp[kind] = self.kind_exp.get(kind, 0.0) + exp
            self.kind_act[kind] = self.kind_act.get(kind, 0.0) + act
            self.idx += 1
        self.t_stage = time.time()


class Step:
    def __init__(self, log: list, progress: Progress | None = None):
        self.log, self.progress = log, progress

    def run(self, name: str, cmd: list, cwd: Path, stdout_path: Path | None = None) -> str:
        env = dict(os.environ, PATH=f"{BIO / 'bin'}:{os.environ.get('PATH', '')}")
        t0 = time.time()
        out_f = open(stdout_path, "wb") if stdout_path else subprocess.PIPE
        try:
            proc = subprocess.Popen([str(c) for c in cmd], cwd=cwd, env=env, stdout=out_f, stderr=subprocess.PIPE)
            # drain stderr in the background so a chatty tool cannot block
            import threading
            err_chunks: list = []
            t = threading.Thread(target=lambda: err_chunks.append(proc.stderr.read()), daemon=True)
            t.start()
            out_chunks: list = []
            if not stdout_path:
                to = threading.Thread(target=lambda: out_chunks.append(proc.stdout.read()), daemon=True)
                to.start()
            last_emit = 0.0
            while proc.poll() is None:
                if self.progress and time.time() - last_emit >= 1.0:
                    self.progress.emit(name)
                    last_emit = time.time()
                time.sleep(0.2)
            t.join()
            if not stdout_path:
                to.join()
        finally:
            if stdout_path:
                out_f.close()
        if self.progress:
            self.progress.finish_stage()
        self.log.append({"step": name, "command": " ".join(map(str, cmd)), "seconds": round(time.time() - t0, 1),
                         "exit": proc.returncode})
        if proc.returncode != 0:
            tail = (err_chunks[0] if err_chunks else b"").decode(errors="replace")[-1500:]
            raise RuntimeError(f"{name} failed (exit {proc.returncode}): {tail}")
        return (out_chunks[0] if out_chunks else b"").decode(errors="replace")

    def python(self, name: str, fn):
        """Time an in-process stage (BAM reading, consensus) as part of the plan."""
        t0 = time.time()
        if self.progress:
            self.progress.emit(name)
        result = fn()
        if self.progress:
            self.progress.finish_stage()
        self.log.append({"step": name, "command": "(in-process)", "seconds": round(time.time() - t0, 1), "exit": 0})
        return result


def _plan(references, iterations: int, stats: bool) -> list[tuple[str, float]]:
    E = EXPECTED
    plan = [("simulate", E["simulate"]), ("merge", E["merge"])]
    for _ in references:
        plan.append(("truth", E["truth"]))
        for it in range(1, iterations + 1):
            last = it == iterations
            plan += [("index", E["index"]), ("map aln", E["map"]), ("map samse", 1), ("sort", E["sort"]),
                     ("index bam", 1), ("flagstat", E["flagstat"]),
                     ("damage stats" if (stats and last) else "damage", E["damage_stats"] if (stats and last) else E["damage_fast"]),
                     ("consensus", E["consensus"])]
    return plan


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
            "minimap2": first(["minimap2", "--version"]), "gargammel": _gargammel_version()}


def _gargammel_version() -> str:
    try:
        commit = (GARGAMMEL / "COMMIT").read_text().strip()[:7]
        notes = (GARGAMMEL / "RELEASE_NOTES.md").read_text().splitlines()
        ver = next((l.strip() for l in notes if re.fullmatch(r"\d+\.\d+\.\d+", l.strip())), "?")
        return f"gargammel {ver} (built from source, commit {commit})"
    except OSError:
        return "unavailable"


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
    obs = observations_from_bam(str(rdir / "reads.bam"), min_mapq=min_mapq)
    return obs, prof, bayes, flag


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
        cons = step.python(f"consensus on {label}",
                           lambda: call_consensus(ref_seq, obs, prof["delta5"], prof["delta3"], min_depth=1))
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


def simulate(step, work: Path, endogenous: str = "NC_007596", coverage: float = 10.0,
             damage: tuple = (0.03, 0.4, 0.01, 0.3), contamination: float = 0.0, length_loc: float = 3.7,
             length_scale: float = 0.35, min_length: int = 25, read_length: int = 75,
             contaminant: str = "NC_012920", seed: int | None = None) -> dict:
    """gargammel paired-end reads from ``endogenous`` with a fraction
    ``contamination`` of present-day ``contaminant`` molecules (no damage),
    then AdapterRemoval trimming and merging. Returns the read files."""
    work = Path(work)
    for d in ("input/endo", "input/cont", "input/bact", "sim", "merge"):
        (work / d).mkdir(parents=True, exist_ok=True)
    shutil.copy(REFS / f"{endogenous}.fa", work / "input/endo/endo.fa")
    shutil.copy(REFS / f"{contaminant}.fa", work / "input/cont/cont.fa")
    shutil.copy(REFS / f"{contaminant}.fa", work / "input/bact/bact.fa")   # placeholder; bacterial fraction is 0
    step.run("simulate (gargammel)", [
        BIO / "bin/perl", GARGAMMEL / "gargammel.pl", *(["--seed", seed] if seed is not None else []),
        "-c", coverage, "--comp", f"0,{contamination},{1 - contamination}",
        "-damagee", ",".join(map(str, damage)), "--loc", length_loc, "--scale", length_scale,
        "--minsize", min_length, "--maxsize", 200, "-rl", read_length, "-ss", "HS25",
        "-o", work / "sim/sim", work / "input"], work)
    s1 = sorted(glob.glob(str(work / "sim/*_s1.fq.gz")))
    s2 = sorted(glob.glob(str(work / "sim/*_s2.fq.gz")))
    if not s1 or not s2:
        raise RuntimeError(f"gargammel produced no paired reads: {os.listdir(work / 'sim')}")
    # trim adapters and merge overlapping pairs into full fragments
    step.run("trim and merge (AdapterRemoval)", [
        BIO / "bin/adapterremoval3", "--in-file1", s1[-1], "--in-file2", s2[-1], "--out-prefix", work / "merge/m",
        "--merge", "--threads", "4"], work)
    merged = next(iter(sorted(glob.glob(str(work / "merge/m*merged*")))), None)
    if merged is None:
        raise RuntimeError(f"AdapterRemoval produced no merged reads: {os.listdir(work / 'merge')}")
    return {"r1": s1[-1], "r2": s2[-1], "merged": merged}


def run(work: str | Path, endogenous: str = "NC_007596", references: tuple = ("NC_007596", "NC_000934", "NC_005129"),
        coverage: float = 10.0, damage: tuple = (0.03, 0.4, 0.01, 0.3), contamination: float = 0.0,
        length_loc: float = 3.7, length_scale: float = 0.35, min_length: int = 25, read_length: int = 75,
        mapdamage_stats: bool = True, min_mapq: int = 25, iterations: int = 3, on_step=None,
        seed: int | None = None) -> dict:
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
    log: list = []
    progress = Progress(_plan(references, max(1, iterations), mapdamage_stats), on_step)
    step = Step(log, progress)
    cat = {g["accession"]: g for g in reference_catalog()["genomes"]}
    truth_id, truth_seq = read_fasta(REFS / f"{endogenous}.fa")
    sim = simulate(step, work, endogenous, coverage, damage, contamination, length_loc, length_scale,
                   min_length, read_length, seed=seed)
    merged = sim["merged"]

    results = []
    for acc in references:
        results.append(_reconstruct(step, work, acc, cat, merged, truth_seq, endogenous, iterations,
                                    mapdamage_stats, min_mapq, on_step))
    return {"simulation": {"endogenous": cat[endogenous], "coverage": coverage,
                           "damage_briggs_v_l_d_s": list(damage), "contamination": contamination,
                           "fragment_length_lognormal": [length_loc, length_scale], "min_length": min_length,
                           "read_length": read_length, "platform": "Illumina HiSeq 2500 profile (ART)", "seed": seed,
                           "iterations": iterations},
            "tools": tool_versions(), "steps": log, "references": results,
            "reference_catalog": reference_catalog()}
