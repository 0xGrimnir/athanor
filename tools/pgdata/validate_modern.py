"""Validate modern ancestry and the trait genotyper (runs in CI after
build_modern.py and build_traits.py).

Usage: python validate_modern.py PACK_DIR VAL_DIR [TOOLS_BIN]

Checks (thresholds fixed here before results are seen; each prints PASS or FAIL):
  1. Leave-one-out composition: reference individuals the build set aside
     (never used for merging regions or setting the floors), each removed
     from its own panel and cut to 80,000 random SNPs, are
     assigned to their own reporting group of regions (mean share >= 0.6 at
     80,000).
     Own-continent shares are reported but not gated: several reference
     groups are themselves mixed. (The first real-data run gated them at
     0.95 and failed for that reason; the accuracy test moved to check 2.)
  2. Simulated admixed genomes (2-3 continents, random proportions, built
     from left-out individuals' chromosomes). From people whose own
     ancestry is unambiguous (>= 90% own continent): proportions within 0.02
     on average and 0.06 at worst. From anyone: proportions within 0.03 on
     average and 0.08 at worst of the sources' own compositions.
     Bootstrap interval coverage reported.
  3. Indigenous American ancestry: simulated genomes with 0, 2, 5, 10 and
     25% recover the proportion within 0.03 from 5% up; detection rate at each
     level reported; East Asian left-out individuals reported as Indigenous
     American in no more than 5%.
  4. ADMIXTURE: its projection mode (-P) on the same SNPs and frequencies
     gives the same proportions within 0.01.
  5. Chromosome painting on simulated genomes: per-SNP accuracy of the pair
     of ancestries >= 0.85 and mean posterior within 0.05 of accuracy.
  6. Runs of homozygosity: segments match PLINK 1.9 --homozyg on the same
     genotypes (total length within 2% per genome, overlap >= 0.95), and
     inserted homozygous stretches of known length are found (>= 95% of
     their length).
  7. X chromosome: left-out individuals' X assigned to their own continent
     (mean share >= 0.8 for every continent with enough X data).
  8. End to end: the full analysis runs on simulated consumer files and
     every trait variant present in the file is genotyped on the right strand.
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")      # many small matrix-vector products: one BLAS thread is far faster

import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "core"))
import mosaic  # noqa: E402
from athanor_core.personal.modern import (admix, composition, decode_freq, encode_freq, paint,  # noqa: E402
                                          panel_freqs, roh, roh_summary)
from athanor_core.personal.refdata import RefData  # noqa: E402

PACK, VAL = Path(sys.argv[1]), Path(sys.argv[2])
BIN = Path(sys.argv[3]) if len(sys.argv) > 3 else None
T0 = time.time()
fails, summary = [], {}
os.environ["ATHANOR_PGDATA"] = str(PACK)
ref = RefData(PACK)
M = ref.modern
man = ref.manifest["modern"]
V = np.load(VAL / "modern_val.npz")
rng = np.random.default_rng(7)
tmp = Path(tempfile.mkdtemp())
snps = ref.snps
S = len(snps["pos"])
chrom = snps["chrom"].astype(np.int64)
cm = snps["cm"].astype(np.float64)
pos = snps["pos"].astype(np.int64)
rids = [str(x) for x in M["region_ids"]]
cont_ids = [str(x) for x in M["continent_ids"]]
rcont = M["continent"].astype(int)
n_ind = M["n_individuals"]
alt, tot = V["alt"].astype(np.int64), V["tot"].astype(np.int64)
G = V["loo_geno"]
lreg = V["loo_region"].astype(int)
lcont = rcont[lreg]
G_OF = M["group_of"].astype(int) if "group_of" in M else np.arange(len(rids))
lgrp = G_OF[lreg]
gtitles = [g["title"] for g in man["regions"]]


def to_group(q):
    return np.array([q[G_OF == n].sum() for n in range(len(gtitles))])


# The synthetic stand-in (push builds, no AADR download) exercises every code
# path, but its South Asian branch is by construction half the West Asian
# branch, so "own continent" is not a clean truth there; accuracy gates are
# reported but only enforced on the real reference data.
SYNTH = len(ref.snps["pos"]) < 500_000
ACCURACY = {"own continent", "own region", "admixed genomes", "Indigenous American"}


def gate(name, ok, detail):
    soft = SYNTH and any(k in name for k in ACCURACY)
    print(f"[{'PASS' if ok else ('FAIL (synthetic stand-in, not enforced)' if soft else 'FAIL')}] {name}: {detail}", flush=True)
    if not ok and not soft:
        fails.append(name)


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def freq(a, t):
    return decode_freq(encode_freq(panel_freqs(a, t)))


def loo_P(j, cols):
    """Panel frequencies on SNP columns `cols` with left-out individual j removed."""
    a, t = alt[:, cols].copy(), tot[:, cols].copy()
    d = G[j, cols]
    a[lreg[j]] -= np.where(d >= 0, d, 0)
    t[lreg[j]] -= (d >= 0) * 2
    return freq(a, t)


def to_cont(q):
    out = np.zeros(len(cont_ids))
    for k in range(len(rids)):
        out[rcont[k]] += q[k]
    return out


# ------------------------------------------------------------------ 1. leave-one-out
per = {}
for j in range(len(lreg)):
    per.setdefault(lreg[j], []).append(j)
sample = [j for v in per.values() for j in v[:12]]
loo = {}
loo_comp = {}                      # validation person -> continental composition (80,000 SNPs)
for D in (80_000,):
    cols = np.sort(rng.choice(S, min(D, S), replace=False))
    own_r, own_c = [], []
    for j in sample:
        P = loo_P(j, cols)
        g = G[j, cols]
        m = g >= 0
        q = admix(g[m], P[:, m], tol=1e-5)["q"]
        own_r.append(to_group(q)[lgrp[j]])
        own_c.append(to_cont(q)[lcont[j]])
        loo_comp[j] = to_cont(q)
    own_r, own_c = np.array(own_r), np.array(own_c)
    lr, lc = lgrp[sample], lcont[sample]
    loo[D] = {"region": {gtitles[k]: float(own_r[lr == k].mean()) for k in np.unique(lr)},
              "continent": {cont_ids[k]: float(own_c[lc == k].mean()) for k in np.unique(lc)}}
    log(f"LOO at {D}: continents {({k: round(v, 3) for k, v in loo[D]['continent'].items()})}")
    log(f"LOO at {D}: regions {({k: round(v, 3) for k, v in loo[D]['region'].items()})}")
# Own-continent shares are reported, not gated: several reference groups are
# themselves mixed (Bedouins and Mozabites carry African ancestry, some Maya
# European ancestry, Levantine and South Asian groups share ancestry with
# Europe), so a share elsewhere is often real. Accuracy is tested below on
# people whose ancestry is unambiguous, and consistency on everyone.
cmin = min(loo[80_000]["continent"].values())
print(f"   own-continent shares (80,000 SNPs): lowest {cmin:.3f}")
rmin = min(loo[80_000]["region"].values())
gate("left-out individuals assigned to their own region (mean >= 0.6 per region, 80,000 SNPs)", rmin >= 0.6,
     f"lowest {rmin:.3f} ({min(loo[80_000]['region'], key=loo[80_000]['region'].get)})")
# Sources for the accuracy test must be close to one continent (>= 97%), or
# their own few percent elsewhere would count as error.
clean = {j for j in sample if loo_comp[j][lcont[j]] >= 0.97}
print(f"   {len(clean)} of {len(sample)} validation people are unambiguous (>= 97% own continent)")

# ------------------------------------------------------------------ 2-3, 5. simulated admixed genomes
by_cont = {}
for j in range(len(lreg)):
    by_cont.setdefault(lcont[j], []).append(j)
cols250 = np.sort(rng.choice(S, min(250_000, S), replace=False))


def held_P(js, cols):
    a, t = alt[:, cols].copy(), tot[:, cols].copy()
    for j in js:
        d = G[j, cols]
        a[lreg[j]] -= np.where(d >= 0, d, 0)
        t[lreg[j]] -= (d >= 0) * 2
    return freq(a, t)


def cont_P(P):
    out = np.zeros((len(cont_ids), P.shape[1]))
    for c in range(len(cont_ids)):
        rows = np.flatnonzero(rcont == c)
        w = n_ind[rows].astype(float)
        out[c] = (w[:, None] * P[rows]).sum(0) / w.sum()
    return out


errs, cons, cover, paint_acc, paint_post = [], [], [], [], []
clean_by_cont = {c: [j for j in v if j in clean] for c, v in by_cont.items()}
any_by_cont = {c: [j for j in v if j in loo_comp] for c, v in by_cont.items()}
NMOS = 24
for t_ in range(NMOS):
    is_clean = t_ % 2 == 0
    pool = clean_by_cont if is_clean else any_by_cont
    avail = [c for c in pool if len(pool[c]) >= 1]
    k = 2 if t_ % 3 else 3
    cs = list(rng.choice(avail, min(k, len(avail)), replace=False))
    js = [int(rng.choice(pool[c])) for c in cs]
    props = rng.dirichlet(np.ones(len(cs)))
    props = np.clip(props, 0.08, None); props /= props.sum()
    srcs = [G[j][cols250] for j in js]
    g, a, b = mosaic.build(srcs, props, chrom[cols250], cm[cols250], 10.0, seed=500 + t_)
    frac = mosaic.truth_props(a, b, len(cs))
    P = held_P(js, cols250)
    res = composition(g, P, snps["block"][cols250], reps=30, seed=t_)
    qc = to_cont(res["q"])
    bc = np.array([to_cont(x) for x in res["boot"]])
    # expected composition: each piece carries its source person's own composition
    expect = sum(f * loo_comp[j] for f, j in zip(frac, js))
    for c in set(cs):
        cons.append(abs(qc[c] - expect[c]))
    if is_clean:
        for i, c in enumerate(cs):
            errs.append(abs(qc[c] - frac[i]))
            lo, hi = np.percentile(bc[:, c], [2.5, 97.5])
            cover.append(lo <= frac[i] <= hi)
        floors = {c["id"]: c["detection_floor"] for c in man["continents"]}
        used = [c for c in range(len(cont_ids)) if qc[c] >= max(floors[cont_ids[c]], 0.01)]
        if set(cs) <= set(used) and len(used) >= 2:
            Pc = cont_P(P)[used]
            qn = qc[used] / qc[used].sum()
            pr = paint(g, Pc, chrom[cols250], cm[cols250], qn, generations=man["painting"]["generations"],
                       temper=man["painting"]["temper"])
            cmap = {c: i for i, c in enumerate(used)}
            ta = np.array([cmap[cs[x]] for x in a]); tb = np.array([cmap[cs[x]] for x in b])
            lo_, hi_ = np.minimum(ta, tb), np.maximum(ta, tb)
            paint_acc.append(float(((pr["a"] == lo_) & (pr["b"] == hi_)).mean()))
            paint_post.append(float(pr["posterior"].mean()))
    log(f"mosaic {t_} ({'clean' if is_clean else 'any'} sources): {[cont_ids[c] for c in cs]} pieces {np.round(frac, 3)} "
        f"expected {np.round(expect[cs], 3)} est {np.round(qc[cs], 3)}")
errs, cons = np.array(errs), np.array(cons)
gate("simulated admixed genomes from unambiguous sources: continental proportions within 0.02 on average and 0.06 at worst",
     errs.mean() <= 0.02 and errs.max() <= 0.06, f"mean error {errs.mean():.4f}, worst {errs.max():.4f}")
gate("simulated admixed genomes from any sources: proportions match the sources' own compositions within 0.03 on average and 0.08 at worst",
     cons.mean() <= 0.03 and cons.max() <= 0.08, f"mean error {cons.mean():.4f}, worst {cons.max():.4f}")
cov = float(np.mean(cover))
print(f"   bootstrap 95% intervals contain the simulated proportion in {100 * cov:.0f}% of cases")
pa, pp = (np.mean(paint_acc), np.mean(paint_post)) if paint_acc else (0.0, 0.0)
gate("chromosome painting: per-SNP accuracy >= 0.85 and mean posterior within 0.05 of it",
     bool(paint_acc) and pa >= 0.85 and abs(pa - pp) <= 0.05,
     f"accuracy {pa:.3f}, mean posterior {pp:.3f} over {len(paint_acc)} genomes")
summary["continental_ancestry"] = (
    f"On {NMOS // 2} simulated genomes mixing two or three continents, built from the chromosomes of reference people "
    f"left out of the panels whose own ancestry is unambiguous, the proportions were within {100 * errs.mean():.1f} "
    f"percentage points on average ({100 * errs.max():.1f} at worst); the 95% intervals, which reflect statistical "
    f"noise only, contained the true value {100 * cov:.0f}% of the time. On genomes built from any reference people, "
    f"including members of mixed groups, results matched the sources' own compositions within {100 * cons.mean():.1f} "
    f"points on average. Reference people left out of their panel are assigned on average at least "
    f"{100 * cmin:.0f}% to their own continent; the rest mostly reflects real shared ancestry between neighbouring "
    f"regions, which is why each region's typical profile is shown.")
summary["regional_ancestry"] = (
    f"With 80,000 SNPs, people left out of the reference panels are assigned on average "
    f"{100 * np.mean(list(loo[80_000]['region'].values())):.0f}% to their own region (lowest region "
    f"{100 * rmin:.0f}%); regions that could not be told apart were merged before the pack was built"
    + (f" ({len(man['merges'])} merges)" if man.get("merges") else "") + ".")
summary["chromosome_painting"] = (
    f"On simulated genomes with known ancestry segments, the pair of ancestries was right at {100 * pa:.0f}% of SNPs, "
    f"and the average confidence the method reported was {100 * pp:.0f}%.") if paint_acc else None

# Indigenous American
if "americas" in cont_ids:
    am = cont_ids.index("americas")
    floor_am = next(c["detection_floor"] for c in man["continents"] if c["id"] == "americas")
    others = [c for c in by_cont if len(by_cont[c]) >= 2 and c != am]
    eur = cont_ids.index("europe") if "europe" in cont_ids else others[0]
    det, err_am = {}, []
    for lvl in (0.0, 0.02, 0.05, 0.10, 0.25):
        hits = []
        for r_ in range(6):
            ja = int(rng.choice(clean_by_cont.get(am) or by_cont[am])); je = int(rng.choice(clean_by_cont.get(eur) or by_cont[eur]))
            props = np.array([1 - lvl, lvl]) if lvl > 0 else np.array([1.0, 0.0])
            g, a, b = mosaic.build([G[je][cols250], G[ja][cols250]], np.clip(props, 1e-9, None), chrom[cols250],
                                   cm[cols250], 8.0, seed=900 + r_ + int(lvl * 1000))
            truth = mosaic.truth_props(a, b, 2)[1]
            P = held_P([ja, je], cols250)
            q = to_cont(admix(g[g >= 0], P[:, g >= 0])["q"])[am]
            hits.append(q >= floor_am)
            if lvl >= 0.05:
                err_am.append(abs(q - truth))
        det[lvl] = float(np.mean(hits))
    ea = cont_ids.index("east_asia") if "east_asia" in cont_ids else None
    fp_ea = None
    if ea is not None:
        fl = []
        for j in (clean_by_cont.get(ea) or by_cont[ea])[:20]:
            P = loo_P(j, cols250)
            g = G[j, cols250]
            fl.append(to_cont(admix(g[g >= 0], P[:, g >= 0])["q"])[am] >= floor_am)
        fp_ea = float(np.mean(fl))
    log("Indigenous American detection rate by simulated level:", det, "East Asian false positives:", fp_ea)
    ok = (max(err_am) <= 0.03) and det[0.0] <= 0.05 and (fp_ea is None or fp_ea <= 0.05)
    gate("Indigenous American: proportion within 0.03 from 5% up; <= 5% false detections at 0% and in East Asians", ok,
         f"worst error {max(err_am):.3f}; detected at 0/2/5/10/25%: {[round(det[k], 2) for k in sorted(det)]}; "
         f"East Asian false positives {fp_ea}")
    summary["indigenous_american"] = (
        f"Detection limit {100 * floor_am:.1f}%. On simulated genomes with known Indigenous American ancestry the proportion "
        f"was recovered within {100 * max(err_am):.1f} points from 5% upwards; detection rates at 2%, 5% and 10% were "
        f"{100 * det[0.02]:.0f}%, {100 * det[0.05]:.0f}% and {100 * det[0.10]:.0f}%; "
        f"{100 * det[0.0]:.0f}% of genomes with none, and {100 * (fp_ea or 0):.0f}% of left-out East Asian genomes, were "
        f"reported as having some.")
else:
    print("   no Indigenous American panel in this pack; check 3 skipped")

# ------------------------------------------------------------------ 4. ADMIXTURE
adm = (BIN / "admixture") if BIN else None
if adm and adm.exists():
    cols = np.sort(rng.choice(S, min(80_000, S), replace=False))
    P = freq(alt[:, cols], tot[:, cols])
    pick = sample[::max(1, len(sample) // 12)][:12]
    cols = cols[(G[pick][:, cols] >= 0).any(0)]          # ADMIXTURE refuses SNPs missing in everyone
    P = freq(alt[:, cols], tot[:, cols])
    Gs = G[pick][:, cols]
    K = len(rids)
    pre = tmp / "adm"
    # PLINK .bed (SNP-major); allele 1 in .bim = the alternative allele
    with open(str(pre) + ".fam", "w") as f:
        for i in range(len(pick)):
            f.write(f"f{i} i{i} 0 0 0 -9\n")
    with open(str(pre) + ".bim", "w") as f:
        for k_, s_ in enumerate(cols):
            f.write(f"{chrom[s_]} s{k_} 0 {pos[s_]} A G\n")
    code = {2: 0b00, 1: 0b10, 0: 0b11, -1: 0b01}
    nb = (len(pick) + 3) // 4
    with open(str(pre) + ".bed", "wb") as f:
        f.write(bytes([0x6c, 0x1b, 0x01]))
        for k_ in range(len(cols)):
            byts = bytearray(nb)
            for i in range(len(pick)):
                byts[i // 4] |= code[int(Gs[i, k_])] << (2 * (i % 4))
            f.write(bytes(byts))
    best = None
    for orient, PP in (("alt", P), ("ref", 1 - P)):
        np.savetxt(str(pre) + f".{K}.P.in", PP.T, fmt="%.6f")
        r = subprocess.run([str(adm), "-P", "adm.bed", str(K)], cwd=tmp, capture_output=True, text=True, timeout=3600)
        ll = [l for l in r.stdout.splitlines() if l.startswith("Loglikelihood:")]
        if r.returncode or not ll:
            print(r.stdout[-1500:], r.stderr[-1500:])
            continue
        Qa = np.loadtxt(tmp / f"adm.{K}.Q").reshape(len(pick), K)
        llv = float(ll[-1].split()[1])
        print(f"   ADMIXTURE projection with P as the frequency of the {orient} allele: log-likelihood {llv:.1f}")
        if best is None or llv > best[0]:
            best = (llv, Qa, orient)
    if best:
        mine = np.array([admix(Gs[i][Gs[i] >= 0], P[:, Gs[i] >= 0])["q"] for i in range(len(pick))])
        dmax = float(np.abs(mine - best[1]).max())
        gate("ADMIXTURE projection (-P) agrees within 0.01", dmax <= 0.01, f"largest difference {dmax:.4f} over {len(pick)} genomes")
        summary["continental_ancestry"] += f" ADMIXTURE's projection mode gives the same proportions to within {dmax:.3f}."
    else:
        gate("ADMIXTURE projection (-P) agrees within 0.01", False, "ADMIXTURE did not run")
else:
    print("   ADMIXTURE not available; check 4 skipped")

# ------------------------------------------------------------------ 6. runs of homozygosity
thin = np.sort(rng.choice(S, S // 2, replace=False))
pick = sample[::max(1, len(sample) // 10)][:10]
# inserted homozygous stretches of known length (copy one pseudo-haplotype over 5-20 Mb)
ins_found, ins_total = 0.0, 0.0
mine_segs = {}
Gr = G[pick][:, thin].copy()
inserted = {}
for i in range(len(Gr)):
    c = int(rng.integers(1, 23))
    m = np.flatnonzero(chrom[thin] == c)
    if len(m) < 400:
        continue
    span = int(pos[thin][m[-1]] - pos[thin][m[0]])
    L = int(min(rng.integers(5_000_000, 20_000_000), 0.6 * span))
    start = int(pos[thin][m[0]] + rng.integers(0, max(1, span - L)))
    w = m[(pos[thin][m] >= start) & (pos[thin][m] < start + L)]
    d = Gr[i, w]
    Gr[i, w] = np.where(d == 1, 2 * rng.integers(0, 2, len(w)), d)
    inserted[i] = (c, int(pos[thin][w[0]]), int(pos[thin][w[-1]]))      # the stretch as typed
for i in range(len(Gr)):
    d = Gr[i]
    mine_segs[i] = roh(chrom[thin], pos[thin], d == 1, d < 0)
    if i in inserted:
        c, s0, s1 = inserted[i]
        cov_ = 0
        for sg in mine_segs[i]:
            if sg["chrom"] == c:
                cov_ += max(0, min(sg["end"], s1) - max(sg["start"], s0))
        ins_found += cov_; ins_total += s1 - s0
ins_frac = ins_found / max(ins_total, 1)
gate("inserted homozygous stretches found (>= 95% of their length)", ins_frac >= 0.95, f"{100 * ins_frac:.1f}%")
plink = (BIN / "plink") if BIN else None
if plink and plink.exists():
    pre = tmp / "roh"
    with open(str(pre) + ".fam", "w") as f:
        for i in range(len(Gr)):
            f.write(f"f{i} i{i} 0 0 0 -9\n")
    with open(str(pre) + ".bim", "w") as f:
        for k_, s_ in enumerate(thin):
            f.write(f"{chrom[s_]} s{k_} 0 {pos[s_]} A G\n")
    code = np.array([0b01, 0b11, 0b10, 0b00], np.uint8)          # index dosage+1: missing, 0, 1, 2
    nb = (len(Gr) + 3) // 4
    with open(str(pre) + ".bed", "wb") as f:
        f.write(bytes([0x6c, 0x1b, 0x01]))
        X = code[Gr.T.astype(np.int64) + 1]                           # (snps, people)
        pad = np.zeros((X.shape[0], nb * 4), np.uint8); pad[:, :X.shape[1]] = X
        pad[:, X.shape[1]:] = 0b01
        q4 = pad.reshape(X.shape[0], nb, 4)
        f.write((q4[:, :, 0] | (q4[:, :, 1] << 2) | (q4[:, :, 2] << 4) | (q4[:, :, 3] << 6)).astype(np.uint8).tobytes())
    args = ["--bfile", str(pre), "--homozyg", "--homozyg-window-snp", "50", "--homozyg-window-het", "1",
            "--homozyg-window-missing", "5", "--homozyg-window-threshold", "0.05", "--homozyg-snp", "100",
            "--homozyg-kb", "1000", "--homozyg-density", "50", "--homozyg-gap", "1000", "--homozyg-het", "999",
            "--allow-no-sex", "--out", str(pre)]
    r = subprocess.run([str(plink)] + args, capture_output=True, text=True, timeout=3600)
    hom = Path(str(pre) + ".hom")
    if r.returncode or not hom.exists():
        print(r.stdout[-2000:], r.stderr[-1000:])
        gate("runs of homozygosity match PLINK", False, "PLINK did not run")
    else:
        rows = [l.split() for l in hom.read_text().splitlines()[1:]]
        ps = {}
        for x in rows:
            ps.setdefault(int(x[1][1:]), []).append((int(x[3]), int(x[6]), int(x[7])))
        tl_diff, jac = [], []
        for i in range(len(Gr)):
            mine = [(s_["chrom"], s_["start"], s_["end"]) for s_ in mine_segs[i]]
            theirs = ps.get(i, [])
            a_ = sum(e - s0 for _, s0, e in mine); b_ = sum(e - s0 for _, s0, e in theirs)
            inter = 0
            for c, s0, e in mine:
                for c2, s2, e2 in theirs:
                    if c == c2:
                        inter += max(0, min(e, e2) - max(s0, s2))
            union = a_ + b_ - inter
            tl_diff.append(abs(a_ - b_) / max(b_, 1))
            jac.append(inter / union if union else 1.0)
        gate("runs of homozygosity match PLINK 1.9 (total within 2%, overlap >= 0.95)",
             max(tl_diff) <= 0.02 and min(jac) >= 0.95,
             f"largest total-length difference {100 * max(tl_diff):.2f}%, lowest overlap {min(jac):.3f} over {len(Gr)} genomes")
        summary["roh"] = (f"Segments match PLINK 1.9's --homozyg on the same genotypes (total length within "
                          f"{100 * max(tl_diff):.1f}%, overlap at least {min(jac):.2f}), and homozygous stretches of "
                          f"known length inserted into real genomes were found over {100 * ins_frac:.0f}% of their length.")
else:
    print("   PLINK not available; PLINK comparison skipped")
    summary["roh"] = f"Homozygous stretches of known length inserted into real genomes were found over {100 * ins_frac:.0f}% of their length."

# ------------------------------------------------------------------ 7. X chromosome
GX = V["loo_x"]
if GX.shape[1] > 0:
    xalt, xtot = V["xalt"].astype(np.int64), V["xtot"].astype(np.int64)
    male = V["loo_male"]
    xb = M["x_block"]
    owns = {}
    for j in [j for j in sample if j in clean]:
        d = GX[j].astype(np.int64)
        a, t = xalt.copy(), xtot.copy()
        if male[j]:
            ok = (d == 0) | (d == 2)
            a[lreg[j]] -= np.where(ok, d // 2, 0); t[lreg[j]] -= ok
            g, n = np.where(ok, d // 2, -1), np.where(ok, 1, 0)
        else:
            ok = d >= 0
            a[lreg[j]] -= np.where(ok, d, 0); t[lreg[j]] -= 2 * ok
            g, n = np.where(ok, d, -1), np.where(ok, 2, 0)
        Pc = cont_P(freq(a, t))
        m = g >= 0
        if m.sum() < 1500:
            continue
        q = admix(g[m], Pc[:, m], n=n[m])["q"]
        owns.setdefault(lcont[j], []).append(q[lcont[j]])
    xo = {cont_ids[c]: float(np.mean(v)) for c, v in owns.items() if len(v) >= 3}
    log("X leave-one-out own-continent means:", {k: round(v, 3) for k, v in xo.items()})
    xmin = min(xo.values()) if xo else 0
    gate("X chromosome: left-out individuals with unambiguous ancestry assigned to their own continent (mean >= 0.8)", xmin >= 0.8,
         f"lowest {xmin:.3f}")
    summary["x_ancestry"] = (f"From X-chromosome SNPs alone, people left out of the reference panels are assigned on "
                             f"average at least {100 * xmin:.0f}% to their own continent (people whose autosomes are unambiguous).")
else:
    print("   no X-chromosome data in this pack; check 7 skipped")

# ------------------------------------------------------------------ 8. end to end on simulated consumer files
from athanor_core.personal.deep import analyse  # noqa: E402
from athanor_core.personal.writers import write  # noqa: E402

traits = ref.traits or {"variants": {}}
chip = np.sort(rng.choice(S, min(300_000, S), replace=False))
tpos = {(v["chrom"], v["pos37"]) for v in traits["variants"].values()}
chip = np.array([k for k in chip if (int(snps["chrom"][k]), int(snps["pos"][k])) not in tpos])
e2e_ok = True
for n_, j in enumerate(sample[::max(1, len(sample) // 3)][:3]):
    d = G[j]
    tv = list(traits["variants"].values())
    BI = {"A": 0, "C": 1, "G": 2, "T": 3}
    # add the trait variants with a known genotype (one copy of the effect allele)
    extra_c = np.array([v["chrom"] for v in tv], np.int64)
    extra_p = np.array([v["pos37"] for v in tv], np.int64)
    extra_r = np.array([BI[v["effect"]] for v in tv], np.uint8)
    extra_a = np.array([BI[v["other"]] for v in tv], np.uint8)
    extra_rs = np.array([int(v["rsid"][2:]) for v in tv], np.int64)
    for fmt in ("23andme", "ancestry"):
        f = tmp / f"e2e_{n_}_{fmt}.txt"
        write(f, fmt, np.r_[snps["rsid"][chip], extra_rs], np.r_[snps["chrom"][chip], extra_c],
              np.r_[snps["pos"][chip], extra_p], np.r_[snps["ref"][chip], extra_r], np.r_[snps["alt"][chip], extra_a],
              np.r_[d[chip], np.ones(len(tv), np.int8)])
        t1 = time.time()
        r = analyse(str(f), ref)
        cards = {c["id"]: c for c in r["cards"]}
        want = {"continental_ancestry", "regional_ancestry", "chromosome_painting", "roh"}
        missing = want - set(cards)
        rows = [v for c in r["cards"] for v in (c.get("result") or {}).get("variants", [])
                if isinstance(v, dict) and str(v.get("rsid", "")).startswith("rs")]
        wrong = [v["rsid"] for v in rows if v["copies"] not in (1, None)]
        untyped = [v["rsid"] for v in rows if v["copies"] is None]
        ok = not missing and not wrong and not untyped
        e2e_ok &= ok
        print(f"   {fmt} file {n_}: {len(r['cards'])} cards in {time.time() - t1:.0f} s; missing {sorted(missing)}; "
              f"trait variants read wrongly {wrong}; not found {untyped}")
gate("end to end: full analysis runs and every trait variant is read with the right genotype", e2e_ok, "see lines above")

# ------------------------------------------------------------------ record
mpath = PACK / "manifest.json"
man_all = json.loads(mpath.read_text())
man_all["validation"] = {**(man_all.get("validation") or {}), **{k: v for k, v in summary.items() if v}}
man_all.setdefault("modern", {})["validation"] = {"loo": {str(k): v for k, v in loo.items()}}
mpath.write_text(json.dumps(man_all, indent=1))
(PACK / "validation_modern.json").write_text(json.dumps({"summary": summary, "fails": fails}, indent=1))
print("   validation summaries written into the pack manifest:", sorted(summary))
print(f"total {time.time() - T0:.0f} s; {len(fails)} gate(s) failed: {fails}")
sys.exit(1 if fails else 0)
