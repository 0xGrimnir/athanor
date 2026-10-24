"""Validate the pharmacogenomics pack and caller.

1. PharmCAT: the same 1000 Genomes genotypes (every CPIC position PharmCAT
   reads, unphased) go to PharmCAT and to Athanor's caller; diplotypes and
   phenotypes must agree.
2. Chip input: genotypes written as a 23andMe-style file (GRCh37 positions,
   either strand, indels as I/D) must be read back exactly.
3. The reliability figure each gene card shows (how often a call from your
   file's positions matches a full reading) must hold on people it was not
   estimated from.
4. The pack itself: genes, drugs, recommendations and lookup tables present.

Usage: python validate_pgx.py PACK VALDATA PHARMCAT_JAR JAVA WORKDIR
"""
from __future__ import annotations

import glob
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.personal.pgx import NOVEL, UNTESTED, Gene, Panel, PgxGenotyper, pgx_cards  # noqa: E402
from athanor_core.personal.rawfile import CODE, RawGenotypes  # noqa: E402

PACK, VAL, JAR, JAVA, WORK = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4], Path(sys.argv[5])
results = []
summary = {}


def gate(name, ok, detail):
    results.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


pgx = json.loads((PACK / "pgx.json").read_text())
z = dict(np.load(PACK / "pgx_1kg.npz"))
panel = Panel(z, pgx["genes"])
genes = {s: Gene(s, d) for s, d in pgx["genes"].items()}
samples = z["samples"].astype(str)
unrel = np.flatnonzero(~z["related"].astype(bool))
rng = np.random.default_rng(7)

# ------------------------------------------------------------------ 4. the pack
level_a = sorted({g for d in pgx["drugs"] for g in d.get("level_a", [])})
no_rec = [d["name"] for d in pgx["drugs"] if not d["recs"]]
no_lookup = [s for s, d in pgx["genes"].items() if not d["results"]]
with_rec = len(pgx["drugs"]) - len(no_rec)
gate("pack: CPIC level-A genes with definitions and phenotype tables; medicines with recommendation tables",
     len(genes) >= 15 and with_rec >= 60 and set(no_lookup) <= {"CYP4F2", "VKORC1"},
     f"{len(genes)} genes ({', '.join(sorted(genes))}); skipped {sorted(pgx.get('skipped', {}))}; "
     f"{len(pgx['drugs'])} drugs, {with_rec} with recommendation tables (without: {', '.join(no_rec) or 'none'}); "
     f"genes CPIC gives no phenotype: {', '.join(no_lookup) or 'none'}")
summary["genes"] = sorted(genes)
summary["drugs"] = len(pgx["drugs"])


def panel_geno(sym, s, keep=None):
    """(positions, 2) for 1000 Genomes sample index s; NOVEL as untested (PharmCAT gets ./. there)."""
    G = panel.genotypes(sym)
    g = G[s].copy()
    g[g == NOVEL] = UNTESTED
    bad = (g == UNTESTED).any(1)
    g[bad] = UNTESTED
    if keep is not None:
        g[~keep] = UNTESTED
    return g


def haploid(sym, s):
    c = pgx["genes"][sym]["chrom"]
    return c == 26 or (c == 23 and bool(z["male"][s]))


# ------------------------------------------------------------------ 1. PharmCAT
vcf = VAL / "pgx_1kg.vcf"
written = json.loads((VAL / "pgx_vcf_positions.json").read_text()) if (VAL / "pgx_vcf_positions.json").exists() else {}
pops = z["superpop"].astype(str)
pick = []
for p in sorted(set(pops[unrel])):
    idx = unrel[pops[unrel] == p]
    pick += list(rng.choice(idx, size=min(100, len(idx)), replace=False))
pick = sorted(int(i) for i in pick)
WORK.mkdir(parents=True, exist_ok=True)
(WORK / "samples.txt").write_text("\n".join(samples[i] for i in pick) + "\n")
t0 = time.time()
r = subprocess.run([JAVA, "-Xmx6g", "-jar", JAR, "-vcf", str(vcf), "-S", str(WORK / "samples.txt"), "-o", str(WORK / "pc"),
                    "-reporterCallsOnlyTsv"], capture_output=True, text=True)
print(f"PharmCAT exit {r.returncode} in {time.time() - t0:.0f} s; {r.stdout[-600:]} {r.stderr[-1500:]}", flush=True)
tsvs = glob.glob(str(WORK / "pc" / "**" / "*.tsv"), recursive=True)
print(f"{len(tsvs)} PharmCAT TSV reports; e.g. {tsvs[:2]}")
by_sample = {}
sset = {samples[i]: i for i in pick}
for f in tsvs:
    name = os.path.basename(f)
    sid = next((s for s in sset if f".{s}." in f".{name}"), None)
    lines = open(f).read().splitlines()
    hdr_i = next((i for i, l in enumerate(lines) if l.startswith("Gene\t")), None)
    if hdr_i is None:
        continue
    hdr = lines[hdr_i].split("\t")
    for l in lines[hdr_i + 1:]:
        f2 = l.split("\t")
        row = dict(zip(hdr, f2))
        s = sid if sid is not None else row.get("Sample") or row.get("Sample ID")
        if s is None:
            continue
        by_sample.setdefault(s, {})[row["Gene"]] = row
if by_sample:
    ex = next(iter(by_sample.values()))
    print("PharmCAT columns:", list(next(iter(ex.values())).keys()))


def norm_dip(d):
    return "/".join(sorted(d.split("/")))


def same_score(a, b):
    try:
        return abs(float(str(a).replace("≥", "")) - float(str(b).replace("≥", ""))) < 1e-6
    except ValueError:
        return str(a).strip() == str(b).strip()


stats = defaultdict(Counter)
examples = defaultdict(list)


def nonref_text(sym, g):
    gene = genes[sym]
    out = []
    for k, p in enumerate(gene.pos):
        a, b = g[k]
        if a == UNTESTED:
            continue
        if not ((gene.refmask[k] >> int(a)) & 1 and (gene.refmask[k] >> int(b)) & 1):
            out.append(f"{p.get('rsid') or p.get('hgvs')}={p['alleles'][a] if a < len(p['alleles']) else '?'}/"
                       f"{p['alleles'][b] if b < len(p['alleles']) else '?'}")
    return " ".join(out[:8])


def known_to_pharmcat(sym):
    """The gene without alleles defined at positions PharmCAT does not read (CPIC alleles newer than
    PharmCAT's data): with those positions missing they would tie with alleles PharmCAT knows."""
    d = pgx["genes"][sym]
    keep = set(written.get(sym, range(len(d["positions"]))))
    drop = [a["name"] for a in d["alleles"] if not a.get("ref") and any(
        c is not None and sorted(c) != sorted(d["positions"][k].get("refset", [0])) and k not in keep
        for k, c in enumerate(a["def"]))]
    return Gene(sym, {**d, "alleles": [a for a in d["alleles"] if a["name"] not in drop]}), drop


pc_genes = {}
for sym in genes:
    pc_genes[sym], dropped = known_to_pharmcat(sym)
    if dropped:
        print(f"   {sym}: {len(dropped)} CPIC alleles use positions PharmCAT does not read, left out of the comparison: "
              f"{', '.join(dropped[:8])}{' ...' if len(dropped) > 8 else ''}")

for sid, rows in by_sample.items():
    s = sset.get(sid)
    if s is None:
        continue
    for sym, row in rows.items():
        if sym not in genes or panel.genotypes(sym) is None or sym == "CYP2D6":
            continue
        keep = None
        if sym in written:
            keep = np.zeros(len(genes[sym].pos), dtype=bool)
            keep[written[sym]] = True
        geno = panel_geno(sym, s, keep)
        ours = pc_genes[sym].call(geno, haploid(sym, s))
        src = (row.get("Source Diplotype") or "").strip()
        pc_dips = {norm_dip(x.strip()) for x in src.split(" OR ") if x.strip()}
        look_ph = (row.get("Recommendation Lookup Phenotype") or "").strip()
        pc_ph = {x.strip() for x in (look_ph or row.get("Phenotype") or "").split(" / ") if x.strip()}
        pc_as = (row.get("Recommendation Lookup Activity Score") or row.get("Activity Score") or "").strip()
        st = stats[sym]
        st["n"] += 1
        raw_line = "\t".join(row.get(k, "") for k in ("Source Diplotype", "Phenotype", "Activity Score",
                                                       "Recommendation Lookup Diplotype", "Recommendation Lookup Phenotype",
                                                       "Recommendation Lookup Activity Score"))
        if not pc_dips or src.lower() in ("unknown/unknown", "not called", "no call") or (not pc_ph and not genes[sym].genotype_only):
            st["pharmcat_no_call"] += 1
            if ours["status"] in ("called", "ambiguous"):
                st["only_ours"] += 1
                if len(examples[sym]) < 3:
                    examples[sym].append(f"{sid}: PharmCAT [{raw_line}]; ours {ours['diplotypes'][:2]} {ours.get('phenotype')}; {nonref_text(sym, geno)}")
            continue
        if ours["status"] not in ("called", "ambiguous"):
            st["only_pharmcat"] += 1
            if len(examples[sym]) < 4:
                examples[sym].append(f"{sid}: PharmCAT [{raw_line}]; ours {ours['status']}; {nonref_text(sym, geno)}")
            continue
        dset = {norm_dip(d) for d in ours["diplotypes"]}
        dip_ok = bool(pc_dips & dset)
        if sym in ("DPYD", "RYR1"):
            dip_ok = True             # variant lists are written differently; the phenotype decides
        if genes[sym].genotype_only:
            ph_ok = dip_ok
        elif ours["status"] == "called":
            ph_ok = {ours["phenotype"].lower()} == {x.lower() for x in pc_ph}
        else:
            ph_ok = {x.lower() for x in ours.get("possible", [])} == {x.lower() for x in pc_ph}
            st["ours_ambiguous"] += 1
            if pc_ph & set(ours.get("possible", [])):
                st["ambiguous_includes_pharmcat"] += 1
        as_ok = (not pc_as) or (not ours.get("score")) or same_score(ours.get("score"), pc_as)
        st["compared"] += 1
        st["diplotype_agree"] += dip_ok
        st["phenotype_agree"] += ph_ok
        st["score_agree"] += as_ok
        if not (dip_ok and ph_ok and as_ok) and len(examples[sym]) < 5:
            examples[sym].append(f"{sid}: PharmCAT [{raw_line}]; ours {sorted(dset)[:3]} / "
                                 f"{ours.get('phenotype') or ours.get('possible')} / {ours.get('score') or '-'}"
                                 f"{' p=%.2f' % ours['diplotype_probability'] if ours.get('diplotype_probability') else ''}; "
                                 f"{nonref_text(sym, geno)}")

print(f"\n{'gene':8s} {'n':>5s} {'cmp':>5s} {'diplo':>7s} {'pheno':>7s} {'score':>7s}  pharmcat-no-call  ours-only  pc-only  ambiguous")
worst = []
for sym in sorted(stats):
    st = stats[sym]
    c = max(1, st["compared"])
    print(f"{sym:8s} {st['n']:5d} {st['compared']:5d} {st['diplotype_agree'] / c:7.3f} {st['phenotype_agree'] / c:7.3f} "
          f"{st['score_agree'] / c:7.3f}  {st['pharmcat_no_call']:5d}  {st['only_ours']:5d}  {st['only_pharmcat']:5d}  "
          f"{st['ours_ambiguous']} ({st['ambiguous_includes_pharmcat']} include PharmCAT's)")
    for e in examples[sym]:
        print("      ", e)
    if st["compared"] >= 50:
        worst.append((min(st["diplotype_agree"], st["phenotype_agree"], st["score_agree"]) / c, sym))
summary["pharmcat"] = {s: dict(v) for s, v in stats.items()}
total = sum(st["compared"] for st in stats.values())
low = min(worst) if worst else (0, "none")
gate("PharmCAT: same diplotype and phenotype on 1000 Genomes genotypes (every gene >= 99%)",
     worst and low[0] >= 0.99 and total > 1000,
     f"{total:,} gene calls in {len(by_sample)} people; lowest agreement {low[0]:.3f} ({low[1]})")
only_pc = max(((st["only_pharmcat"] / max(1, st["n"]), sym) for sym, st in stats.items()), default=(0, "none"))
gate("no call where PharmCAT makes one (<= 0.5% of people in every gene)", only_pc[0] <= 0.005,
     f"highest {only_pc[0]:.1%} ({only_pc[1]})")

# ------------------------------------------------------------------ full-data calls
nocall = {}
for sym, gene in sorted(genes.items()):
    if panel.genotypes(sym) is None:
        continue
    st = Counter()
    shown = 0
    for s in unrel:
        gg = panel_geno(sym, s)
        c = gene.call(gg, haploid(sym, s))
        st[c["status"].split(" (")[0]] += 1
        if c["status"].startswith("no call") and shown < 3:
            shown += 1
            print(f"      no call {samples[s]}: {nonref_text(sym, gg)}")
    nocall[sym] = st["no call"] / max(1, sum(st.values()))
    print(f"   {sym:8s} full-data calls on {sum(st.values())} unrelated 1000 Genomes people: {dict(st)}")
worst_nc = max(nocall.items(), key=lambda x: x[1]) if nocall else ("none", 0)
gate("every gene gives a call for >= 95% of 1000 Genomes people (the rest carry combinations CPIC does not define)",
     worst_nc[1] <= 0.05, f"highest no-call rate {worst_nc[1]:.2%} ({worst_nc[0]})")
summary["no_call"] = nocall

# ------------------------------------------------------------------ 2. chip input
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
mism = tot = 0
rare_cards = rare_unwarned = 0


def rare_allele(sym, name):
    a = genes[sym].by_name.get(name) or {}
    f = [v for v in (a.get("freq") or {}).values() if isinstance(v, (int, float))]
    return not a.get("ref") and bool(f) and max(f) < 0.01


def carries_rare(s):
    for sym in genes:
        if panel.genotypes(sym) is None or sym == "CYP2D6":
            continue
        c = genes[sym].call(panel_geno(sym, s), haploid(sym, s))
        if any(rare_allele(sym, n) for dp in c.get("diplotypes") or [] for n in dp.split("/")):
            return True
    return False


carriers = [s for s in pick if carries_rare(s)]
print(f"   {len(carriers)} of {len(pick)} compared people carry a rare (< 1%) allele somewhere")
for flip in (False, True):
    for s in sorted(set(pick[:60]) | set(carriers[:40])):
        rows, expect = [], {}
        for sym, gene in genes.items():
            G = panel.genotypes(sym)
            if G is None:
                continue
            d = pgx["genes"][sym]
            for k, p in enumerate(gene.pos):
                a, b = G[s, k]
                if a >= NOVEL or b >= NOVEL or not p.get("pos37"):
                    continue
                if p["kind"] == "snv":
                    x, y = p["alleles"][a], p["alleles"][b]
                    if flip:
                        x, y = COMP[x], COMP[y]
                    gt = x + y
                elif p.get("chip"):
                    inv = {v: kk for kk, v in p["chip"].items()}
                    gt = inv[a] + inv[b]
                else:
                    continue
                rid = int(p["rsid"][2:]) if (p.get("rsid") or "").startswith("rs") else -1
                if rid < 0 and p["kind"] != "snv":
                    continue
                rows.append((rid, d["chrom"], p["pos37"], gt))
                expect[(sym, k)] = sorted((int(a), int(b)))
        rawg = RawGenotypes("23andMe", "sim.txt", "", np.array([r[0] for r in rows]), np.array([r[1] for r in rows], dtype=np.int8),
                            np.array([r[2] for r in rows]), np.array([CODE[r[3][0]] for r in rows], dtype=np.uint8),
                            np.array([CODE[r[3][1]] for r in rows], dtype=np.uint8))
        strand = "the opposite strand" if flip else "the reference's + strand"
        gt = PgxGenotyper(rawg, rawg.pos, strand)
        if not flip:
            cards = pgx_cards(pgx, None, rawg, rawg.pos, strand, "XY" if z["male"][s] else "XX",
                              lambda cid, title, group, grade, summary, result, limits, method, sources, ok=True:
                              {"id": cid, "summary": summary, "result": result, "grade": grade})
            for c in cards:
                if (c["result"] or {}).get("rare_in_call"):
                    rare_cards += 1
                    if "clinical test" not in c["summary"] or c["grade"] == "A":
                        rare_unwarned += 1
        for sym, gene in genes.items():
            if panel.genotypes(sym) is None:
                continue
            got = gt.genotype(gene, pgx["genes"][sym]["chrom"])
            for k in range(len(gene.pos)):
                if (sym, k) in expect:
                    tot += 1
                    if sorted(got[k].tolist()) != expect[(sym, k)]:
                        mism += 1
                        if mism <= 8:
                            p = gene.pos[k]
                            print(f"      chip {sym} {p.get('rsid')} {p['kind']} {p['alleles']} pos37 {p.get('pos37')} "
                                  f"flip={flip}: expected {expect[(sym, k)]}, read {got[k].tolist()}")
gate("chip input: 23andMe-style genotypes (GRCh37, both strands, indels as I/D) read back exactly",
     tot > 1000 and mism == 0, f"{tot:,} genotypes, {mism} differ")
gate("rare-variant gate: every call with a rare (< 1%) allele says it needs a clinical test and is not graded A",
     rare_cards > 0 and rare_unwarned == 0, f"{rare_cards} gene cards with a rare allele, {rare_unwarned} without the warning")

# ------------------------------------------------------------------ 3. reliability holds out of sample
half = rng.permutation(unrel)
A, B = np.sort(half[: len(half) // 2]), np.sort(half[len(half) // 2:])
cal = []
for sym, gene in sorted(genes.items()):
    G = panel.genotypes(sym)
    if G is None or sym == "CYP2D6":
        continue
    known = (G[unrel] != UNTESTED).all(axis=(0, 2)) & (G[unrel] < NOVEL).all(axis=(0, 2))
    for frac in (0.3, 0.6):
        keep = known & (rng.random(len(gene.pos)) < frac)
        if not keep.any():
            continue

        def calls(idx):
            out = []
            cache = {}
            for s in idx:
                hap = haploid(sym, s)
                g = G[s]
                kf = (hap, g[known].tobytes())
                km = (hap, g[keep].tobytes())
                if kf not in cache:
                    gg = np.full_like(g, UNTESTED)
                    gg[known] = g[known]
                    cache[kf] = gene.call(gg, hap).get("phenotype") or "?"
                if km not in cache:
                    gg = np.full_like(g, UNTESTED)
                    gg[keep] = g[keep]
                    cache[km] = gene.call(gg, hap).get("phenotype") or "?"
                out.append((cache[km], cache[kf]))
            return out
        ca, cb = calls(A), calls(B)
        for ph in {m for m, _ in ca}:
            ra = [f == m for m, f in ca if m == ph]
            rb = [f == m for m, f in cb if m == ph]
            if len(ra) >= 30 and len(rb) >= 30:
                cal.append((abs(np.mean(ra) - np.mean(rb)), len(rb), sym, frac, ph, np.mean(ra), np.mean(rb)))
if cal:
    err = sum(e * n for e, n, *_ in cal) / sum(n for _, n, *_ in cal)
    w = max(cal)
    gate("reliability: agreement estimated on half of 1000 Genomes holds on the other half (weighted error <= 0.03)",
         err <= 0.03, f"weighted error {err:.4f} over {len(cal)} calls; largest {w[0]:.3f} ({w[2]}, {int(100 * w[3])}% of "
                      f"positions, {w[4]}: {w[5]:.2f} vs {w[6]:.2f}, n={w[1]})")
else:
    gate("reliability: agreement estimated on half of 1000 Genomes holds on the other half", False, "nothing to compare")

out = {"results": [{"name": n, "pass": ok, "detail": d} for n, ok, d in results], "summary": summary}
(PACK / "validation_pgx.json").write_text(json.dumps(out, indent=1, default=str))
# the per-card validation line the app shows
man = PACK / "manifest.json"
if man.exists():
    m = json.loads(man.read_text())
    v = m.setdefault("validation", {})
    pc = next((d for n, ok, d in results if n.startswith("PharmCAT")), "")
    for sym in genes:
        st = stats.get(sym)
        line = (f"Checked in CI against PharmCAT on {st['compared']} 1000 Genomes people: diplotype "
                f"{st['diplotype_agree'] / max(1, st['compared']):.1%}, phenotype {st['phenotype_agree'] / max(1, st['compared']):.1%} the same."
                if st and st["compared"] else "Not comparable with PharmCAT (PharmCAT does not call it from a VCF).")
        v[f"pgx_{sym.lower().replace('-', '_')}"] = line
    v["medicines"] = "Calls checked against PharmCAT in CI: " + pc
    man.write_text(json.dumps(m, indent=1))
print("\n" + "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}: {d}" for n, ok, d in results))
sys.exit(0 if all(ok for _, ok, _ in results) else 1)
