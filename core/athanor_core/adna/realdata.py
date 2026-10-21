"""Ancient-DNA pipeline for real sequencing data (mitochondrial scale).

  reads           FASTQ or BAM in, molecules out (reads.py)
  reconstruction  iterative mapping with BWA aln (ancient-DNA settings) to a
                  circular reference (bamtools.py), duplicate removal, damage
                  measured against the sample's own consensus, damage-aware
                  posterior calls (consensus.py), majority-rule indels;
                  repeated until the sequence stops changing
  authentication  terminal deamination with intervals, fragment lengths,
                  mapDamage2 profiles of the final mapping
  outputs         consensus FASTA with an N mask (depth and posterior
                  thresholds), per-base posteriors and depths; the
                  majority-rule consensus (>= 3 reads, > 50%, 2 bases trimmed
                  from read ends) as a conventional comparison
  contamination   competitive mapping against the reconstructed sequence and
                  the human mitogenome (contamination.py; Feuerborn et al.,
                  BMC Genomics 21, 844, 2020)

Everything runs locally in the Linux engines pack. Nothing is uploaded.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np

from .bamtools import (apply_changes, clean_bam, coverage_span, elongate, indel_events, majority_indels,
                       write_fasta)
from .compare import compare_sequences
from .consensus import BASES, CODE, call_consensus, majority_consensus
from .contamination import competitive_mapping
from .damage import damage_from_observations, damage_profiles, terminal_signal
from .pileup import observations_from_bam
from .pipeline import BIO, EXPECTED, REFS, Progress, Step, read_fasta, reference_catalog, tool_versions
from .reads import prepare_reads, sha256


def _map_stages() -> list:
    E = EXPECTED
    return [("index", E["index"]), ("map aln", E["map"]), ("map samse", 1), ("sort", E["sort"]),
            ("clean", E["clean"]), ("consensus", E["consensus"])]


def _prep_stages(inputs) -> list:
    out = []
    for inp in inputs:
        if inp["kind"] == "bam":
            out.append(("read_input", EXPECTED["read_input"]))
            if not inp.get("trimmed", True):
                out.append(("trim", EXPECTED["trim"]))
        elif inp.get("path2") or not inp.get("trimmed", False):
            out.append(("trim", EXPECTED["trim"]))
    return out


def _contam_stages() -> list:
    E = EXPECTED
    return [("index", E["index"]), ("map aln", E["map"]), ("map samse", 1), ("sort", E["sort"]),
            ("clean", E["clean"]), ("count", E["count"])]


def _final_stages(stats: bool, compare: int) -> list:
    return [("damage stats" if stats else "damage", EXPECTED["damage_stats"] if stats else EXPECTED["damage_fast"])] \
        + [("compare", EXPECTED["compare"])] * compare


def _site_observations(obs, p: int, limit: int = 60) -> list[dict]:
    """Every read base at site p (reference orientation), with its strand and
    distance from each molecule end, so a difference can be traced to
    damage at read ends or to genuine disagreement among reads."""
    idx = np.flatnonzero(obs.site == p)[:limit]
    full = obs.full if obs.full is not None else np.ones(len(obs), bool)
    return [{"base": BASES[int(obs.base[i])], "strand": "-" if obs.reverse[i] else "+",
             "from_5p": int(obs.dist5[i]), "from_3p": int(obs.dist3[i]) if full[i] else None,
             "quality": int(obs.qual[i])} for i in idx]


def run_real(work: str | Path, inputs: list[dict], reference: str = "NC_007596", reference_fasta: str | None = None,
             circular: bool = True, contaminant: str | None = "NC_012920", iterations: int = 4,
             min_mapq: int = 25, min_length: int = 30, min_depth: int = 3, min_confidence: float = 0.99,
             mapdamage_stats: bool = False, threads: int = 4, sample: str = "sample",
             compare_with: tuple[str, str] | None = None, on_step=None) -> dict:
    """Reconstruct a mitogenome from real reads.

    inputs: [{"kind": "fastq", "path": ..., "path2": optional mate file,
              "trimmed": bool, "read_length": optional},
             {"kind": "bam", "path": ..., "trimmed": True}]
    reference: accession of a genome in the pack, or ``reference_fasta``.
    compare_with: (name, sequence) of a published sequence to compare with."""
    work = Path(work)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    log: list = []
    iterations = max(1, int(iterations))
    n_compare = 2 if compare_with else 0
    plan = _prep_stages(inputs) + _map_stages() * iterations + (_contam_stages() if contaminant else []) \
        + _final_stages(mapdamage_stats, n_compare)
    progress = Progress(plan, on_step)
    step = Step(log, progress)

    cat = {g["accession"]: g for g in reference_catalog()["genomes"]} if (REFS / "manifest.json").exists() else {}
    ref_path = Path(reference_fasta) if reference_fasta else REFS / f"{reference}.fa"
    ref_id, ref_seq = read_fasta(ref_path)
    ref_info = cat.get(reference, {"accession": ref_id, "file": ref_path.name})
    ref_info = dict(ref_info, sha256=sha256(ref_path), length=len(ref_seq))

    # 1. reads
    prepared = prepare_reads(step, work, inputs, min_length=min_length, threads=threads)
    fastq = prepared["fastq"]
    prefiltered = any(s.get("all_reads_mapped") for s in prepared["inputs"])
    progress.expect(_map_stages() * iterations + (_contam_stages() if contaminant else [])
                    + _final_stages(mapdamage_stats, n_compare))

    # 3. iterative reconstruction
    cur, origin = list(ref_seq), list(range(len(ref_seq)))
    history, applied = [], []
    converged = False
    oscillating = False
    seen = {hash("".join(cur))}
    for it in range(1, iterations + 1):
        L = len(cur)
        seq = "".join(cur)
        d = work / "map" / f"it{it}"
        d.mkdir(parents=True)
        fa = write_fasta(d / "ref.fa", {"ref": elongate(seq, circular)})
        label = f"iteration {it}"
        step.run(f"index reference, {label} (bwa)", [BIO / "bin/bwa", "index", fa], d)
        step.run(f"map, {label} (bwa aln)", [BIO / "bin/bwa", "aln", "-l", 1024, "-n", 0.04, "-o", 2,
                                             "-t", threads, fa, fastq], d, stdout_path=d / "reads.sai")
        step.run(f"map, {label} (bwa samse)", [BIO / "bin/bwa", "samse", fa, d / "reads.sai", fastq],
                 d, stdout_path=d / "reads.sam")
        step.run(f"sort, {label} (samtools)", [BIO / "bin/samtools", "sort", "-o", d / "raw.bam", d / "reads.sam"], d)
        (d / "reads.sam").unlink()
        mapping = step.python(f"fold circle, remove duplicates, {label}",
                              lambda: clean_bam(d / "raw.bam", d / "clean.bam", min_mapq,
                                                {"ref": L} if circular else None))
        prior_ref = "".join(ref_seq[o] if o >= 0 else "N" for o in origin)

        def analyse():
            obs = observations_from_bam(str(d / "clean.bam"), min_mapq=0, circular_length=L if circular else None)
            # damage measured first against the current reference, then against the sample's consensus
            prof0 = damage_from_observations(obs, np.array([CODE.get(b, -1) for b in seq]))
            c0 = call_consensus(prior_ref, obs, prof0["delta5"], prof0["delta3"], min_depth=1)
            trusted = (c0.confidence >= 0.99) & (c0.depth >= min_depth)
            prof = damage_from_observations(obs, c0.posterior.argmax(axis=1), trusted)
            cons = call_consensus(prior_ref, obs, prof["delta5"], prof["delta3"], min_depth=1)
            span = coverage_span(d / "clean.bam", L, circular)
            dels, ins = indel_events(d / "clean.bam", L, circular)
            return obs, prof, cons, majority_indels(dels, ins, span, min_reads=min_depth)

        obs, prof, cons, indels = step.python(f"consensus, {label}", analyse)
        subs = {i: cons.call[i] for i in range(L)
                if cons.confidence[i] >= 0.99 and cons.depth[i] > 0 and cons.call[i] != cur[i]}
        history.append({"iteration": it, "reference_length": L, "mapping": mapping,
                        "covered_fraction": float((cons.depth > 0).mean()), "mean_depth": float(cons.depth.mean()),
                        "estimated_divergence_theta": cons.theta,
                        "substitutions_found": len(subs), "indels_found": len(indels)})
        if not subs and not indels:
            converged = True
            progress.expect((_contam_stages() if contaminant else []) + _final_stages(mapdamage_stats, n_compare))
            break
        if it < iterations:
            new, new_origin = apply_changes(cur, origin, subs, indels)
            key = hash("".join(new))
            if key in seen:
                # the changes lead back to a sequence already mapped against:
                # the calls alternate between states, so iterating further cannot settle them
                oscillating = True
                progress.expect((_contam_stages() if contaminant else []) + _final_stages(mapdamage_stats, n_compare))
                break
            seen.add(key)
            applied.append({"after_iteration": it, "substitutions": len(subs), "indels": indels})
            cur, origin = new, new_origin

    # 4. outputs on the final mapping
    L = len(cur)
    conf = cons.confidence
    called = (cons.depth >= min_depth) & (conf >= min_confidence)
    athanor_seq = "".join(c if k else "N" for c, k in zip(cons.call, called))
    majority_seq, counts, mdepth = majority_consensus(obs, L, min_depth=min_depth, majority=0.5, trim=2)

    # 5. contamination: competitive mapping against the sample's own sequence
    # (uncalled sites filled from the mapping reference) and the human
    # mitogenome. Using the sample's sequence rather than the species
    # reference keeps endogenous reads from being lost to sample-reference
    # differences, which would inflate the estimate.
    contamination = None
    if contaminant:
        h_id, h_seq = read_fasta(REFS / f"{contaminant}.fa")
        target_seq = "".join(c if c in "ACGT" else r for c, r in zip(cons.call, cur))
        contamination = competitive_mapping(step, work, fastq, (f"{sample}_consensus", target_seq), (h_id, h_seq),
                                            threads=threads, prefiltered=prefiltered)

    md = [BIO / "bin/mapDamage", "-i", d / "clean.bam", "-r", fa, "-d", d / "mapdamage", "-q"]
    if not mapdamage_stats:
        md.append("--no-stats")
    step.run("damage profile (mapDamage2)", md, d)
    try:
        md_prof = damage_profiles(d / "mapdamage")
        md_prof = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in md_prof.items()}
    except Exception as e:  # noqa: BLE001
        md_prof = {"error": str(e)}

    comparison = None
    if compare_with:
        comparison = {
            "athanor": compare_sequences(step, work / "cmp_athanor", athanor_seq, compare_with[1], compare_with[0]),
            "majority_rule": compare_sequences(step, work / "cmp_majority", majority_seq, compare_with[1], compare_with[0]),
        }
        for c in comparison.values():
            for s in c["substitutions"]:
                p = s["ours_pos"]
                s["depth"] = int(cons.depth[p])
                s["posterior"] = float(conf[p])
                s["base_counts_trimmed"] = dict(zip(BASES, map(int, counts[p])))
                s["observations"] = _site_observations(obs, p)

    fasta = (f">{sample} Athanor consensus (reference {ref_id}; N where depth < {min_depth} "
             f"or posterior < {min_confidence})\n" + "\n".join(athanor_seq[i:i + 70] for i in range(0, L, 70)) + "\n")
    depth = cons.depth
    pending = [] if converged else [{"substitutions": sorted(int(i) for i in subs), "indels": indels}]
    return {
        "mode": "real data",
        "sample": sample, "inputs": prepared["inputs"], "fragment_lengths": prepared["fragment_lengths"],
        "reference": ref_info, "circular": circular,
        "settings": {"iterations": iterations, "min_mapq": min_mapq, "min_length": min_length,
                     "min_depth": min_depth, "min_posterior": min_confidence, "mapdamage_stats": mapdamage_stats,
                     "bwa_aln": "-l 1024 -n 0.04 -o 2"},
        "contamination": contamination,
        "iterations": history, "changes_applied": applied, "converged": converged, "oscillating": oscillating,
        "pending_changes": pending,
        "mapping": history[-1]["mapping"],
        "coverage": {"length": L, "mean_depth": float(depth.mean()), "median_depth": float(np.median(depth)),
                     "fraction_ge1": float((depth >= 1).mean()), "fraction_ge3": float((depth >= 3).mean()),
                     "fraction_ge10": float((depth >= 10).mean())},
        "authentication": {"damage_vs_consensus": {k: (v.tolist() if hasattr(v, "tolist") else v)
                                                   for k, v in prof.items()},
                           "terminal_signal": terminal_signal(prof), "mapdamage": md_prof,
                           "fragment_lengths": prepared["fragment_lengths"]},
        "consensus": {"sequence": athanor_seq, "called_sites": int(called.sum()),
                      "estimated_divergence_theta": cons.theta,
                      "predicted_errors_in_called_sites": float((1 - conf[called]).sum()),
                      "confidence": conf.round(6).tolist(), "depth": depth.tolist(), "fasta": fasta},
        "majority_rule": {"sequence": majority_seq, "called_sites": int(sum(b != "N" for b in majority_seq)),
                          "rule": f">= {min_depth} reads after trimming 2 bases from each read end; > 50% agreement"},
        "comparison": comparison,
        "tools": tool_versions(), "steps": log,
    }
