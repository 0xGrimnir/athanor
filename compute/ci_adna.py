"""Run the ancient-DNA pipeline on the mammoth demo inside the pack image and
report the validation numbers (CI has no GPU; nothing here needs one)."""
import json
import sys
import time

from athanor_core.adna.pipeline import run

t0 = time.time()
res = run("/tmp/adna", coverage=10.0, on_step=lambda d: print("stage:", d.get("stage"), flush=True))
out = sys.argv[1] if len(sys.argv) > 1 else "/tmp/adna.json"
json.dump(res, open(out, "w"))
sim = res["simulation"]
print("tools:", json.dumps(res["tools"]))
print("simulated damage (v,l,d,s):", sim["damage_briggs_v_l_d_s"])
for r in res["references"]:
    b = r["mapdamage_bayesian"]
    est = {k: round(v.get("Mean", float("nan")), 4) for k, v in b.items() if isinstance(v, dict)}
    c = r["calibration"]
    print(f"[{r['reference']['accession']}] {r['reference']['organism']}: mapped={r['mapped_reads']} "
          f"depth={r['mean_depth']:.1f} covered={r['covered_fraction']:.3f} "
          f"truth-ref divergence={r['truth_vs_reference_divergence']:.4f} theta={r['estimated_divergence_theta']:.4f} "
          f"accuracy={r['accuracy']:.5f} predicted_err={c['predicted_error_rate']:.5f} observed_err={c['observed_error_rate']:.5f} "
          f"ECE={c['expected_calibration_error']:.5f} gaps={r['truth_gap_positions']} unaligned={r['unaligned_positions']}")
    print("   mapDamage Bayesian means:", est)
    dg = r["error_diagnostics"]
    if dg.get("n"):
        print(f"   errors: n={dg['n']} called-reference-where-truth-differs={dg['called_reference_where_truth_differs']} "
              f"within-10bp-of-indel={dg['within_10bp_of_truth_indel']} median-depth errors/all={dg['median_depth_at_errors']}/{dg['median_depth_overall']} "
              f"mean-confidence-at-errors={dg['mean_confidence_at_errors']:.4f}")
        for ex in dg["examples"][:8]:
            print("     ", ex)
    print("   5' C>T raw:", [round(x, 4) for x in r["damage_profile"]["raw_5p_CtoT"][:6]],
          " background:", round(r["damage_profile"]["background_CtoT"], 4))
print(f"total {time.time() - t0:.0f} s")
