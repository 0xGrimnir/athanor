"""Run the ancient-DNA pipeline on the mammoth demo inside the pack image and
report the validation numbers (CI has no GPU; nothing here needs one)."""
import json
import sys
import time

from athanor_core.adna.pipeline import run

t0 = time.time()
events = []


def on_step(d):
    events.append(dict(d, t=time.time() - t0))
    if len(events) % 20 == 1:
        print(f"progress {100 * d['fraction']:.1f}% eta {d['eta_seconds']}s step {d['step']}/{d['steps']} {d['stage']}", flush=True)


res = run("/tmp/adna", coverage=10.0, iterations=3, on_step=on_step, seed=11)   # reproducible: gargammel --seed
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
    for h in r["iterations"]:
        print(f"   iteration {h['iteration']}: mapped={h['mapped_reads']} covered={h['covered_fraction']:.3f} "
              f"accuracy={h['accuracy']:.5f} predicted={h['predicted_error']:.5f} observed={h['observed_error']:.5f} "
              f"ECE={h['ece']:.5f} theta={h['estimated_divergence_theta']:.4f}")
    dg = r["error_diagnostics"]
    if dg.get("n"):
        print(f"   errors: n={dg['n']} called-reference-where-truth-differs={dg['called_reference_where_truth_differs']} "
              f"within-10bp-of-indel={dg['within_10bp_of_truth_indel']} median-depth errors/all={dg['median_depth_at_errors']}/{dg['median_depth_overall']} "
              f"mean-confidence-at-errors={dg['mean_confidence_at_errors']:.4f}")
        for ex in dg["examples"][:8]:
            print("     ", ex)
    print("   5' C>T raw:", [round(x, 4) for x in r["damage_profile"]["raw_5p_CtoT"][:6]],
          " background:", round(r["damage_profile"]["background_CtoT"], 4))
total = time.time() - t0
print(f"total {total:.0f} s")
# How good was the ETA? Compare predicted remaining time with the actual remaining time.
for mark in (0.1, 0.25, 0.5, 0.75):
    e = next((x for x in events if x["fraction"] >= mark), None)
    if e:
        actual = total - e["t"]
        print(f"[INFO] ETA at {int(mark * 100)}%: predicted {e['eta_seconds']} s, actual {actual:.0f} s")
fr = [x["fraction"] for x in events]
print(f"[INFO] progress events {len(events)}, monotone {fr == sorted(fr)}, steps planned {events[-1]['steps'] if events else '?'}, "
      f"steps run {len(res['steps'])}")

# ---- validation gates -------------------------------------------------------
# Principled checks on the simulated truth; each prints PASS/FAIL.
fails = []

def gate(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        fails.append(name)

v, lam, dd, ds = sim["damage_briggs_v_l_d_s"]
endo = sim["endogenous"]["accession"]
for r in res["references"]:
    acc = r["reference"]["accession"]
    b = r["mapdamage_bayesian"]
    if acc == endo:
        gate("self-reference accuracy >= 99.9%", r["accuracy"] >= 0.999, f"{r['accuracy']:.5f}")
        gate("mapDamage recovers deltaD within 30%", abs(b["DeltaD"]["Mean"] - dd) <= 0.3 * dd, f"{b['DeltaD']['Mean']:.4f} vs {dd}")
        gate("mapDamage recovers lambda within 20%", abs(b["Lambda"]["Mean"] - lam) <= 0.2 * lam, f"{b['Lambda']['Mean']:.4f} vs {lam}")
        # gargammel 1.1.5 implements the Briggs model as published, which mapDamage2 fits (1.1.4 drew an
        # overhang at every end, and deltaS came out about twice the setting)
        gate("mapDamage recovers deltaS within 35% (gargammel 1.1.5's Briggs model)", abs(b["DeltaS"]["Mean"] - ds) <= 0.35 * ds,
             f"{b['DeltaS']['Mean']:.4f} vs {ds}")
        continue
    h = r["iterations"]
    gate(f"{acc}: iterative mapping does not lose accuracy", h[-1]["accuracy"] >= h[0]["accuracy"],
         f"{h[0]['accuracy']:.5f} -> {h[-1]['accuracy']:.5f}")
    true_div = r["truth_vs_reference_divergence"]
    gate(f"{acc}: iterative mapping moves divergence estimate toward the truth",
         abs(h[-1]["estimated_divergence_theta"] - true_div) <= abs(h[0]["estimated_divergence_theta"] - true_div),
         f"{h[0]['estimated_divergence_theta']:.4f} -> {h[-1]['estimated_divergence_theta']:.4f} (true {true_div:.4f})")
    gate(f"{acc}: final consensus >= 99.9% correct", r["accuracy"] >= 0.999, f"{r['accuracy']:.5f}")
    gate(f"{acc}: expected calibration error <= 0.002", r["calibration"]["expected_calibration_error"] <= 0.002,
         f"{r['calibration']['expected_calibration_error']:.5f}")
    gate(f"{acc}: mapDamage deltaD within 50% after iteration", abs(b["DeltaD"]["Mean"] - dd) <= 0.5 * dd,
         f"{b['DeltaD']['Mean']:.4f} vs {dd}")
gate("progress plan matches the steps actually run", bool(events) and events[-1]["steps"] == len(res["steps"]),
     f"planned {events[-1]['steps'] if events else 0}, ran {len(res['steps'])}")
gate("progress never goes backwards", fr == sorted(fr), "monotone" if fr == sorted(fr) else "decreased")
sys.exit(1 if fails else 0)
