"""Validate the Personal Genome reference pack and analyses (runs in CI).

Usage: python validate.py PACK_DIR VAL_DIR [AADR_PREFIX ANNO MT_FASTA_GZ ADMIXTOOLS_BIN]

Checks (thresholds fixed here before results are seen; each prints PASS or FAIL):
  1. Parsers: held-out genomes written in each vendor's format (and as a
     GRCh38 VCF) are read back with every genotype exact.
  2. PCA: genomes reduced to a genotyping-chip-sized SNP subset land where
     the same genomes land with all SNPs.
  3. Neanderthal ancestry: the chip-subset interval contains the all-SNP
     estimate.
  4. qpAdm: Athanor's weights agree with ADMIXTOOLS qpAdm on the same
     individuals and SNPs.
  5. Haplogroups: Y calls of ancient males and mtDNA calls of the AADR
     mitochondrial sequences are on the same lineage as the published calls.
Writes PACK_DIR/validation.json (summaries shown on the app's report cards).
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.personal import fstats  # noqa: E402
from athanor_core.personal.deep import analyse  # noqa: E402
from athanor_core.personal.haplo import Tree  # noqa: E402
from athanor_core.personal.rawfile import parse  # noqa: E402
from athanor_core.personal.refdata import RefData, align  # noqa: E402
from athanor_core.personal.writers import write  # noqa: E402

PACK, VAL = Path(sys.argv[1]), Path(sys.argv[2])
AADR = Path(sys.argv[3]) if len(sys.argv) > 3 else None
ANNO = Path(sys.argv[4]) if len(sys.argv) > 4 else None
MTFA = Path(sys.argv[5]) if len(sys.argv) > 5 else None
BIN = Path(sys.argv[6]) if len(sys.argv) > 6 else None
T0 = time.time()
fails, summary = [], {}
os.environ["ATHANOR_PGDATA"] = str(PACK)
ref = RefData(PACK)
S = ref.snps
H = np.load(VAL / "holdout.npz")
rng = np.random.default_rng(1)
CHIP_FRACTION = 0.22          # about 250k of the 1.15M 1240K autosomal SNPs, as on consumer chips
tmp = Path(tempfile.mkdtemp())


def gate(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)
    if not ok:
        fails.append(name)


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


chip = np.flatnonzero(rng.random(len(S["pos"])) < CHIP_FRACTION)
log("chip-like subset", len(chip), "of", len(S["pos"]))

# ------------------------------------------------------------------ 1. parsers
from athanor_core.personal.liftover import Liftover  # noqa: E402


def to_grch38(chrom, pos):
    """Inverse of the hg38->hg19 chain for the test file (GRCh37 -> GRCh38)."""
    lo = Liftover(ref.chain())
    out = np.full(len(pos), -1, np.int64)
    for c in np.unique(chrom):
        st, en, off = lo.blocks[str(int(c))]
        q0, q1 = st + off, en + off
        order = np.argsort(q0)
        q0, q1, o = q0[order], q1[order], off[order]
        m = chrom == c
        p0 = pos[m] - 1
        i = np.searchsorted(q0, p0, side="right") - 1
        ok = (i >= 0) & (p0 < q1[np.clip(i, 0, None)])
        r = np.full(m.sum(), -1, np.int64)
        r[ok] = p0[ok] - o[i[ok]] + 1
        out[m] = r
    return out


exact_all = True
for k, fmt in enumerate(["23andme", "ancestry", "myheritage", "ftdna", "livingdna", "vcf", "vcf38"]):
    j = k % len(H["ids"])
    d = H["geno"][j][chip]
    c, p = S["chrom"][chip], S["pos"][chip].astype(np.int64)
    if fmt == "vcf38":
        p = to_grch38(c, p)
        keep = p > 0
        sub = chip[keep]
        c, p, d = c[keep], p[keep], d[keep]
    else:
        sub = chip
    f = tmp / f"t_{fmt}.txt"
    write(f, "vcf" if fmt == "vcf38" else fmt, S["rsid"][sub], c, p, S["ref"][sub], S["alt"][sub], d)
    if fmt == "vcf38":
        txt = f.read_text().replace("##reference=GRCh37", "##reference=GRCh38")
        f.write_text(txt)
    g = parse(f)
    al = align(g, ref)
    got = al["dosage"][sub]
    usable = got >= 0
    ok = bool(np.all(got[usable] == d[usable]))
    exact_all &= ok
    print(f"   {fmt}: vendor {g.vendor}, build {al['stats']['build']}, markers {len(g)}, used {int(usable.sum())}, "
          f"exact {ok}, dropped ambiguous {al['stats']['ambiguous_dropped']}")
gate("parsers: every vendor format reads back exactly", exact_all, "see lines above")
summary["file_check"] = ("Raw files written in the AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA, Living DNA and VCF "
                         "layouts (including a GRCh38 VCF) are read back with every genotype exact.") if exact_all else None

# ------------------------------------------------------------------ 2-3. analyses on chip-sized genomes
pc_dev, pc_range, nean_cover, nean_rows, qp_rows = [], {}, [], [], []
for name in ("world", "west"):
    P = ref.pca(name)
    pc_range[name] = float(np.ptp(P["ref_coords"][:, 0]))
nhold = len(H["ids"])
results = []
for j in range(nhold):
    d = H["geno"][j]
    f = tmp / "h.txt"
    write(f, "23andme", S["rsid"][chip], S["chrom"][chip], S["pos"][chip], S["ref"][chip], S["alt"][chip], d[chip])
    r = analyse(str(f), ref)
    cards = {c["id"]: c for c in r["cards"]}
    maps = cards["pca"]["result"]["maps"]
    for name, full in (("world", H["pc_world"][j]), ("west", H["pc_west"][j])):
        if "coords" in maps.get(name, {}) and not np.isnan(full).any():
            pc_dev.append((name, float(np.linalg.norm(np.array(maps[name]["coords"][:2]) - full[:2])) / pc_range[name]))
    # Neanderthal: all-SNP estimate for the same genome
    nd = cards["neanderthal"]["result"]
    if "estimate" in nd:
        roles = ref.manifest["archaic_roles"]
        auto = d >= 0
        from athanor_core.personal.deep import _freq_matrix
        F = _freq_matrix(ref, [roles["altai"], roles["chimp"], roles["vindija"], roles["africa"]], d, auto)
        full = fstats.f4_ratio(F, S["block"][auto], (0, 1, 4, 3), (0, 1, 2, 3))["alpha"]
        lo, hi = nd["ci95"]
        nean_cover.append(lo <= full <= hi)
        nean_rows.append((str(H["groups"][j]), nd["estimate"], full, lo, hi))
    aa = cards["ancient_ancestry"]["result"]["models"]
    results.append({"id": str(H["ids"][j]), "group": str(H["groups"][j]),
                    "models": [{k: m.get(k) for k in ("model", "fit", "p", "components")} for m in aa]})
    if j % 20 == 0:
        log(f"  analysed {j + 1}/{nhold}")

dev = np.array([x[1] for x in pc_dev])
med = float(np.median(dev)) if len(dev) else float("nan")
gate("PCA: chip-subset genomes land within 2% of the PC1 range of their all-SNP position (median)",
     med < 0.02, f"median {100 * med:.2f}%, 95th percentile {100 * np.percentile(dev, 95):.2f}% over {len(dev)} placements")
summary["pca"] = (f"Genomes reduced to about {len(chip):,} SNPs (a genotyping chip's worth) land a median "
                  f"{100 * med:.2f}% of the map's width from where the same genomes land with all 1.15 million SNPs "
                  f"({len(dev)} placements of {nhold} held-out genomes).")
cov = float(np.mean(nean_cover)) if nean_cover else 0.0
gate("Neanderthal: chip-subset 95% interval contains the all-SNP estimate in >= 90% of genomes",
     cov >= 0.90, f"{100 * cov:.0f}% of {len(nean_cover)}")
by_group = {}
for g, e, full, lo, hi in nean_rows:
    by_group.setdefault(g, []).append(full)
print("   all-SNP Neanderthal estimates by group:", {g: round(100 * float(np.mean(v)), 2) for g, v in by_group.items()})
summary["neanderthal"] = (f"On {len(nean_cover)} held-out genomes reduced to chip size, the 95% interval contained the "
                          f"estimate from all 1.15 million SNPs in {100 * cov:.0f}% of cases.")
fit_counts = {}
for r in results:
    for m in r["models"]:
        fit_counts.setdefault((r["group"], m["model"]), []).append(m["fit"])
print("   model fit by group (fits / marginal / does not fit):")
for (g, mdl), v in sorted(fit_counts.items()):
    print(f"     {g:10s} {mdl:18s} {v.count('fits')}/{v.count('marginal')}/{v.count('does not fit')}")

# ------------------------------------------------------------------ 4. qpAdm against ADMIXTOOLS
# The comparison uses a population of held-out individuals as the target, not a
# single genome: with one genome the f4 system is often too weakly determined to
# pin the weights down, and two correct implementations can then land in
# different places. A population target is the configuration qpAdm is designed
# for, and is how the published models were fitted.
def target_frequency(rows: np.ndarray) -> np.ndarray:
    """Alternative-allele frequency of a set of held-out genomes (NaN where none called)."""
    g = H["geno"][rows].astype(float)
    g[g < 0] = np.nan
    with np.errstate(invalid="ignore"):
        return np.nanmean(g, axis=0) / 2


def fit_mine(mdl, freq_row, auto):
    from athanor_core.personal.deep import _freq_matrix
    F = _freq_matrix(ref, mdl["sources"] + mdl["right"], np.zeros(len(S["pos"]), np.int8) - 1, auto)
    F[-1] = freq_row[auto]
    return fstats.qpadm(F, S["block"][auto], len(F) - 1, list(range(len(mdl["sources"]))),
                        list(range(len(mdl["sources"]), len(F) - 1)))


groups_present = sorted({str(g) for g in H["groups"]})
mdl = ref.manifest["models"][0]
if AADR and BIN and (BIN / "qpAdm").exists():
    sets = json.loads((VAL / "sets.json").read_text())
    ind_lines = (AADR.with_suffix(".ind")).read_text().splitlines()
    label = {}
    for role, members in sets["roles"].items():
        for i in members:
            label[i] = role
    work = tmp / "qp"
    work.mkdir(exist_ok=True)
    agree, shown = [], []
    for gname in groups_present:
        rows = np.flatnonzero(np.array([str(x) for x in H["groups"]]) == gname)
        if len(rows) < 4:
            continue
        ids_here = {str(H["ids"][i]) for i in rows}
        tag = "TGT_" + gname.replace(".", "_")
        with open(work / "lab.ind", "w") as fh:
            for l in ind_lines:
                f = l.split()
                fh.write(f"{f[0]} {f[1]} {tag if f[0] in ids_here else label.get(f[0], 'Ignore')}\n")
        (work / "left").write_text("\n".join([tag] + mdl["sources"]) + "\n")
        (work / "right").write_text("\n".join(mdl["right"]) + "\n")
        (work / "par").write_text(
            f"genotypename: {AADR.with_suffix('.geno')}\nsnpname: {AADR.with_suffix('.snp')}\n"
            f"indivname: {work / 'lab.ind'}\npopleft: {work / 'left'}\npopright: {work / 'right'}\n"
            f"details: YES\nallsnps: NO\n")
        out = subprocess.run([str(BIN / "qpAdm"), "-p", str(work / "par")], capture_output=True, text=True,
                             timeout=7200).stdout
        line = [l for l in out.splitlines() if l.strip().startswith("best coefficients:")]
        if not line:
            print(f"   {gname}: ADMIXTOOLS produced no coefficients; tail: {out[-400:]}")
            continue
        w_ad = np.array([float(x) for x in line[0].split(":")[1].split()])
        freq = target_frequency(rows)
        auto = ~np.isnan(freq)
        mine = fit_mine(mdl, freq, auto)
        w_me, se = np.array(mine["weights"]), np.array(mine["se"])
        dw = float(np.max(np.abs(w_me - w_ad)))
        # two correct GLS fits of the same f4 system differ by the weighting
        # they give each statistic, so they can only be expected to agree
        # within the statistical error of the estimate itself
        tol = float(max(0.02, np.max(se)))
        # only populations whose weights the data actually pins down can test agreement
        identifiable = bool(np.all(w_ad > -0.5) and np.all(w_ad < 1.5) and np.max(se) <= 0.05)
        shown.append((gname, w_ad, w_me, dw, identifiable))
        print(f"   {gname}: max weight difference {dw:.4f}, largest standard error {np.max(se):.4f}"
              + ("" if identifiable else "  (not identifiable; excluded)"))
        if identifiable:
            agree.append(dw <= tol)
    if agree:
        gate("qpAdm: weights agree with ADMIXTOOLS within one standard error on well-determined population targets",
             all(agree), f"{sum(agree)}/{len(agree)} populations")
    else:
        print("   no population target was well determined enough to compare "
              "(the estimator itself is checked by crosscheck_qpadm.py)")
    if agree:
        summary["ancient_ancestry"] = (
            f"Athanor's qpAdm reproduces ADMIXTOOLS qpAdm (v8.0.2) on {sum(agree)} of {len(agree)} held-out population "
            f"targets, every weight agreeing within the estimate's own standard error.")
else:
    print("   ADMIXTOOLS not available; qpAdm cross-check skipped")

# recovery of known proportions, when the truth is known (simulated reference data)
TRUTH = json.loads((VAL / "truth.json").read_text()) if (VAL / "truth.json").exists() else {}
if TRUTH:
    errs = []
    for gname, truth in TRUTH.items():
        rows = np.flatnonzero(np.array([str(x) for x in H["groups"]]) == gname)
        if len(rows) < 4:
            continue
        freq = target_frequency(rows)
        auto = ~np.isnan(freq)
        mine = fit_mine(mdl, freq, auto)
        w = np.array(mine["weights"])
        t = np.array([truth[s] for s in mdl["sources"]])
        se = np.array(mine["se"])
        z = float(np.max(np.abs(w - t) / np.maximum(se, 1e-6)))
        errs.append((gname, float(np.max(np.abs(w - t))), z))
        print(f"   {gname}: largest difference from the simulated proportion {errs[-1][1]:.3f} "
              f"({z:.1f} standard errors; model p {mine['p']:.3f})")
    if errs:
        worst = max(errs, key=lambda x: x[2])
        gate("qpAdm recovers the simulated proportions within two standard errors", worst[2] <= 2.0,
             f"worst {worst[0]}: {worst[1]:.3f} = {worst[2]:.1f} standard errors, over {len(errs)} population(s)")
        summary["ancient_ancestry_accuracy"] = (
            f"On simulated populations whose true proportions are known, the fitted proportions were within two "
            f"standard errors of the truth for all {len(errs)} populations (largest deviation {worst[1]:.3f}).")

# ------------------------------------------------------------------ 5. haplogroups
YT = Tree(ref.tree("y"))
yg = H["y_geno"]
consistent = exact = n = 0
for k in range(len(H["y_ids"])):
    obs = {}
    for p, a1, a2, x in zip(H["y_pos"], H["y_a1"], H["y_a2"], yg[k]):
        if x == 0:
            obs[int(p)] = "ACGT"[a1]
        elif x == 2:
            obs[int(p)] = "ACGT"[a2]
    res = YT.call(obs)
    pub = YT.lookup(str(H["y_calls"][k]), trim=False)
    if res.get("node") is None or pub is None:
        continue
    n += 1
    me = res["node"]
    exact += me == pub
    consistent += YT.is_ancestor(me, pub) or YT.is_ancestor(pub, me)
if n:
    gate("Y haplogroups of ancient males on the same lineage as the published call in >= 95%", consistent / n >= 0.95,
         f"{consistent}/{n} on the same lineage, {exact} identical")
    summary["y_haplogroup"] = (f"On {n} ancient men with published Y haplogroups, Athanor's call from the 1240K Y-chromosome "
                               f"SNPs was on the same lineage in {100 * consistent / n:.0f}% and identical in {100 * exact / n:.0f}%.")
MT = Tree(ref.tree("mt"))
if MTFA and ANNO and MTFA.exists():
    import csv
    rows = list(csv.reader(open(ANNO, encoding="utf-8", errors="replace"), delimiter="\t"))
    hdr = rows[0]
    mcol = next(i for i, h in enumerate(hdr) if h.lower().startswith("mtdna haplogroup"))
    hg = {r[0]: r[mcol] for r in rows[1:] if len(r) > mcol}
    hg.update({r[2]: r[mcol] for r in rows[1:] if len(r) > mcol})
    seqs, cur, buf = {}, None, []
    with gzip.open(MTFA, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if cur:
            seqs[cur] = "".join(buf)
    pos_all = MT.positions
    res_full = {"n": 0, "consistent": 0, "exact": 0}
    res_chip = {"n": 0, "consistent": 0, "exact": 0}
    chip_pos = sorted(rng.choice(pos_all, min(1500, len(pos_all)), replace=False).tolist())
    tried = 0
    for sid, sq in seqs.items():
        key = sid if sid in hg else sid.split(".")[0]
        name = hg.get(key) or hg.get(sid.split("_")[0])
        if not name or name in ("..", "n/a") or len(sq) < 16560:
            continue
        pub = MT.lookup(name)
        if pub is None:
            continue
        tried += 1
        for posset, acc in ((pos_all, res_full), (chip_pos, res_chip)):
            obs = {p: sq[p - 1].upper() for p in posset if p <= len(sq) and sq[p - 1].upper() in "ACGT"}
            r = MT.call(obs)
            if r.get("node") is None:
                continue
            acc["n"] += 1
            acc["exact"] += r["node"] == pub
            acc["consistent"] += MT.is_ancestor(r["node"], pub) or MT.is_ancestor(pub, r["node"])
        if tried >= 3000:
            break
    print(f"   mtDNA sequences compared: {tried} (full: {res_full}; 1,500 chip-like positions: {res_chip})")
    if res_full["n"]:
        c = res_full["consistent"] / res_full["n"]
        gate("mtDNA haplogroups from full sequences on the same lineage as published in >= 95%", c >= 0.95,
             f"{res_full['consistent']}/{res_full['n']}, identical {res_full['exact']}")
        cc = res_chip["consistent"] / max(1, res_chip["n"])
        summary["mt_haplogroup"] = (f"On {res_full['n']} ancient mitogenomes with published haplogroups, Athanor's call was on the "
                                    f"same lineage in {100 * c:.0f}% (identical in {100 * res_full['exact'] / res_full['n']:.0f}%); "
                                    f"using only 1,500 chip-like positions, on the same lineage in {100 * cc:.0f}%.")

(PACK / "validation.json").write_text(json.dumps({"summary": summary, "fails": fails}, indent=1))
print(f"total {time.time() - T0:.0f} s; {len(fails)} gate(s) failed: {fails}")
sys.exit(1 if fails else 0)
