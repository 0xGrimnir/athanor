"""Validate the real-data ancient-DNA pipeline (runs inside the pack image).

A. Simulations with a known truth
   A1  contamination: mammoth reads with 0%, 2% and 10% present-day human
       molecules; the competitive-mapping estimate must match the true
       fraction of human molecules among the reads analysed.
   A2  reconstruction across a divergent reference: mammoth reads mapped
       iteratively starting from the African elephant mitogenome (about 9%
       divergent, with indels); the consensus must agree with the mammoth
       truth at the sites it calls.
B. Reproduction of published mitogenomes (Hohle Fels mammoth ivory; reads
   ENA PRJEB103048, sequences GenBank PZ244049-PZ244060): the deposited
   reads are run through the pipeline and the result is compared with the
   deposited sequence.

Gates are fixed here before the results are seen; each prints PASS or FAIL.
Usage: python ci_reproduce.py REPRO_DIR OUT_JSON
"""
import gzip
import json
import sys
import time
from pathlib import Path

from athanor_core.adna.pipeline import REFS, Step, read_fasta, simulate
from athanor_core.adna.realdata import run_real

repro, out_json = Path(sys.argv[1]), Path(sys.argv[2])
fails, report = [], {"simulations": {}, "reproduction": []}


def gate(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)
    if not ok:
        fails.append(name)


def slim(res: dict) -> dict:
    r = json.loads(json.dumps(res))
    r["consensus"].pop("confidence", None)
    r["consensus"].pop("depth", None)
    r["consensus"].pop("fasta", None)
    return r


def true_fraction(merged_fq: Path, min_length: int, contaminant_id: str) -> tuple[int, int]:
    """gargammel names each read after the sequence it came from."""
    cont = total = 0
    with gzip.open(merged_fq, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:]
            elif i % 4 == 1 and len(line.strip()) >= min_length:
                total += 1
                cont += name.startswith(contaminant_id)
    return cont, total


t0 = time.time()
# ---- A1: contamination ------------------------------------------------------
for c in (0.0, 0.02, 0.10):
    w = Path(f"/tmp/contam_{c}")
    sim = simulate(Step([]), w / "sim", "NC_007596", coverage=6.0, contamination=c)
    res = run_real(w / "run", [{"kind": "fastq", "path": sim["r1"], "path2": sim["r2"]}],
                   reference="NC_007596", iterations=2, sample=f"sim-contamination-{c}")
    k, n = true_fraction(w / "run/reads/in0/ar.merged.fastq.gz", 30, read_fasta(REFS / "NC_012920.fa")[0])
    truth = k / n
    ct = res["contamination"]
    est, (lo, hi) = ct["fraction"], ct["ci95"]
    print(f"contamination {c}: true {truth:.4f} ({k}/{n}) estimate {est:.4f} [{lo:.4f}, {hi:.4f}] "
          f"human-read 5' C>T {ct['human_reads_detail']['damage']['CtoT_5p_position1']['frequency']:.3f} "
          f"target 5' C>T {ct['target_reads']['damage']['CtoT_5p_position1']['frequency']:.3f}", flush=True)
    gate(f"contamination {c}: estimate within 20% + 0.5 points of the truth",
         abs(est - truth) <= 0.2 * truth + 0.005, f"{est:.4f} vs {truth:.4f}")
    report["simulations"][f"contamination_{c}"] = {"true": truth, "estimate": est, "ci95": [lo, hi],
                                                   "detail": ct}

# ---- A2: reconstruction from a divergent reference ---------------------------
_, mammoth = read_fasta(REFS / "NC_007596.fa")
w = Path("/tmp/divergent")
sim = simulate(Step([]), w / "sim", "NC_007596", coverage=15.0, contamination=0.02)
res = run_real(w / "run", [{"kind": "fastq", "path": sim["r1"], "path2": sim["r2"]}],
               reference="NC_000934", iterations=8, sample="sim-mammoth-from-elephant",
               compare_with=("NC_007596 (truth)", mammoth))
cmp = res["comparison"]["athanor"]
print("divergent reference: iterations", [(h["iteration"], h["mapping"].get("kept"), h["substitutions_found"],
                                           h["indels_found"]) for h in res["iterations"]],
      "converged", res["converged"], flush=True)
print(f"  agreement {cmp['agreement']:.6f} on {cmp['sites_called_by_both']} sites; truth aligned "
      f"{cmp['theirs_aligned_fraction']:.4f}; uncalled {cmp['uncalled_ours']}; substitutions {len(cmp['substitutions'])}; "
      f"indel differences {len(cmp['insertions_in_ours']) + len(cmp['deletions_in_ours'])}; "
      f"predicted errors {res['consensus']['predicted_errors_in_called_sites']:.2f}", flush=True)
for s in cmp["substitutions"][:20]:
    print("   ", s)
gate("divergent reference: agreement with truth >= 99.9% of called sites", cmp["agreement"] >= 0.999,
     f"{cmp['agreement']:.6f}")
gate("divergent reference: >= 95% of the truth called", cmp["sites_called_by_both"] >= 0.95 * len(mammoth),
     f"{cmp['sites_called_by_both']} of {len(mammoth)}")
report["simulations"]["divergent_reference"] = slim(res)

# ---- B: reproduction of published mitogenomes ---------------------------------
manifest = json.loads((repro / "manifest.json").read_text())
for s in manifest["samples"]:
    _, deposited = read_fasta(repro / s["isolate"] / "deposited.fa")
    inputs = [{"kind": "bam", "path": str(repro / f["path"]), "trimmed": True} for f in s["files"]]
    t1 = time.time()
    res = run_real(Path("/tmp/repro") / s["isolate"], inputs, reference="NC_007596", iterations=6,
                   sample=s["isolate"], compare_with=(s["genbank"], deposited))
    a, m = res["comparison"]["athanor"], res["comparison"]["majority_rule"]
    ts = res["authentication"]["terminal_signal"]
    fl = res["fragment_lengths"]
    cov = res["coverage"]
    print(f"\n== {s['isolate']} -> {s['genbank']} ({len(deposited)} bp; {len(inputs)} BAMs; "
          f"{sum(f['reads'] for f in s['files'])} reads) {time.time() - t1:.0f} s", flush=True)
    print(f"  mapping {res['mapping']}; coverage mean {cov['mean_depth']:.1f}, >=3 reads {cov['fraction_ge3']:.4f}")
    print(f"  iterations {[(h['iteration'], h['substitutions_found'], h['indels_found']) for h in res['iterations']]}"
          f" converged {res['converged']} oscillating {res['oscillating']}; pending {res['pending_changes']}; "
          f"theta vs NC_007596 {res['consensus']['estimated_divergence_theta']:.5f}")
    print(f"  fragments: whole {fl['whole_molecules']} median {fl['median']}; cut by read length {fl['cut_by_read_length']}")
    print(f"  5' C>T pos1 {ts['CtoT_5p_position1']['frequency']:.4f} pos2 {ts['CtoT_5p_position2']['frequency']:.4f} "
          f"interior {ts['CtoT_interior']['frequency']:.4f}; 3' G>A pos1 {ts['GtoA_3p_position1']['frequency']:.4f}")
    ctm = res["contamination"]
    print(f"  contamination {ctm['fraction']:.4f} [{ctm['ci95'][0]:.4f}, {ctm['ci95'][1]:.4f}] ({ctm['human_reads']}/{ctm['assigned_reads']})")
    for name, c in (("athanor", a), ("majority rule", m)):
        print(f"  {name}: agreement {c['agreement']:.6f} on {c['sites_called_by_both']} sites; aligned "
              f"{c['theirs_aligned_fraction']:.4f}; uncalled ours {c['uncalled_ours']} theirs {c['uncalled_theirs']}; "
              f"substitutions {len(c['substitutions'])}; insertions {len(c['insertions_in_ours'])}; "
              f"deletions {len(c['deletions_in_ours'])}")
        for x in c["substitutions"][:12]:
            print("     ", {k: v for k, v in x.items() if k != "observations"})
            for o in x.get("observations", []):
                print("         read base", o)
        for x in (c["insertions_in_ours"] + c["deletions_in_ours"])[:12]:
            print("      indel", x)
    print(f"  predicted errors among called sites: {res['consensus']['predicted_errors_in_called_sites']:.3f}")
    gate(f"{s['isolate']}: terminal deamination excess (authentic ancient DNA)", ts["terminal_excess"],
         f"{ts['CtoT_5p_position1']['frequency']:.4f} vs interior {ts['CtoT_interior']['frequency']:.4f}")
    gate(f"{s['isolate']}: Athanor consensus agrees with {s['genbank']} at >= 99.9% of sites both call",
         a["agreement"] is not None and a["agreement"] >= 0.999, f"{a['agreement']}")
    gate(f"{s['isolate']}: majority-rule consensus agrees with {s['genbank']} at >= 99.9% of sites both call",
         m["agreement"] is not None and m["agreement"] >= 0.999, f"{m['agreement']}")
    report["reproduction"].append({"sample": s, "result": slim(res)})

print(f"\ntotal {time.time() - t0:.0f} s; {len(fails)} gate(s) failed: {fails}")
out_json.write_text(json.dumps(report))
sys.exit(1 if fails else 0)
