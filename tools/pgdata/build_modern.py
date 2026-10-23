"""Modern-ancestry part of the Personal Genome reference pack.

Run after build.py, on the same AADR release. Writes into the pack:
  modern.npz          allele frequencies of the regional reference panels on
                      the pack's autosomal SNPs and on the X chromosome
                      (arcsine-encoded uint8), panel sizes, continents
  roh_reference.json  runs of homozygosity in present-day reference groups
and adds a "modern" section to manifest.json. The validation bundle gets the
genotypes and per-region allele counts that validate_modern.py needs.

How the panels are made:
  1. Present-day individuals from the regional groups below (1000 Genomes,
     HGDP and other present-day genomes in the AADR), one version each.
  2. Outliers are removed: each individual's continental composition is
     estimated leaving that individual out of the reference; anyone whose
     own-continent share is more than 0.15 below their group's median is
     dropped. This removes mislabelled and recently admixed
     individuals without removing a group's shared history.
  3. Left-out individuals of each region are split once into a calibration
     set and a validation set (the latter is used only by validate_modern.py).
     Each calibration individual's regional composition is estimated with
     that individual removed from its own panel, on 80,000 random SNPs (the
     low end of what consumer files share with the reference). A region whose
     members are assigned on average under 65% to their own region is merged
     into the region on the same continent that absorbs most of the rest;
     repeat until every region passes.
  4. Detection floors: for each region and continent, the 99th percentile of
     the proportion assigned to it among calibration individuals from other
     continents whose own ancestry is unambiguous (at least 90% assigned to
     their own continent), and at least 1%. The average result of each
     region's calibration individuals ships as that region's typical profile. Smaller proportions are reported as not distinguishable
     from noise.
  5. Chromosome painting's down-weighting of correlated SNPs is fitted on
     simulated genomes built from left-out individuals' chromosomes.

Usage: python build_modern.py AADR_PREFIX ANNO PACK_DIR VALDIR
"""
from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")      # many small matrix-vector products: one BLAS thread is far faster

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "core"))
from aadr import AADR  # noqa: E402
import mosaic  # noqa: E402
from athanor_core.personal.fstats import blocks_from_cm  # noqa: E402
from athanor_core.personal.modern import (admix, decode_freq, encode_freq, paint, panel_freqs, roh,  # noqa: E402
                                          roh_summary)

PREFIX, ANNO, PACK, VALDIR = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
T0 = time.time()
rng = np.random.default_rng(20261023)


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


CONTINENTS = [("africa", "Sub-Saharan Africa"), ("europe", "Europe"), ("west_asia", "West Asia and North Africa"),
              ("south_asia", "South Asia"), ("east_asia", "East and Southeast Asia"),
              ("north_asia", "Siberia"), ("oceania", "Oceania"), ("americas", "Indigenous Americas")]
# (id, title, continent, AADR group labels). Labels with "-o" (outliers flagged by
# the AADR) and duplicate "Discovery" versions are not listed and so not used.
REGIONS = [
    ("west_africa", "West Africa", "africa", ["YRI", "ESN", "GWD", "MSL", "Yoruba", "Mandenka", "Mende", "Gambian", "Igbo"]),
    ("east_africa", "East Africa (Bantu-speaking)", "africa", ["LWK", "Luhya", "BantuKenya"]),
    ("south_central_africa", "Central and Southern Africa (Bantu-speaking)", "africa",
     ["BantuSA", "Namibia_Bantu_Herero", "BantuSA_Ovambo", "Kongo"]),
    ("rainforest_hg", "Central African rainforest hunter-gatherers", "africa", ["Mbuti", "MbutiPygmy", "Biaka"]),
    ("san", "Southern African San", "africa", ["Ju_hoan_North", "Khomani_San"]),
    ("nw_europe", "Northwestern Europe", "europe", ["CEU", "GBR", "English", "Scottish", "Orcadian", "Icelandic", "Norwegian"]),
    ("france", "France", "europe", ["French"]),
    ("finland", "Finland", "europe", ["FIN", "Finnish"]),
    ("east_europe", "Eastern Europe", "europe", ["Russian", "Polish", "Czech", "Hungarian", "Estonian", "Mari"]),
    ("iberia", "Iberia", "europe", ["IBS", "Spanish", "Spanish_North"]),
    ("basque", "Basque Country", "europe", ["Basque"]),
    ("italy", "Italy", "europe", ["TSI", "Italian_North"]),
    ("sardinia", "Sardinia", "europe", ["Sardinian"]),
    ("balkans", "Balkans and Greece", "europe", ["Greek", "Greek_Crete", "Greek_1", "Albanian", "Bulgarian"]),
    ("levant", "Levant", "west_asia", ["Palestinian", "Jordanian", "Samaritan"]),
    ("druze", "Druze", "west_asia", ["Druze"]),
    ("arabia", "Arabian Peninsula and Bedouin", "west_asia", ["BedouinA", "BedouinB", "BedouinB_1", "BedouinB_2", "Jew_Yemenite"]),
    ("caucasus", "Caucasus", "west_asia", ["Adygei", "Georgian", "Armenian", "Abkhasian", "Lezgin", "Chechen", "Ossetian", "Avar"]),
    ("north_africa", "North Africa", "west_asia", ["Mozabite", "Saharawi"]),
    ("nw_south_asia", "Northwestern South Asia", "south_asia",
     ["PJL", "Pakistani", "Punjabi", "Pathan", "Sindhi_Pakistan", "Balochi", "Brahui", "Burusho"]),
    ("gujarat", "Gujarat", "south_asia", ["GIH"]),
    ("bengal", "Bengal", "south_asia", ["BEB"]),
    ("south_india", "South India and Sri Lanka", "south_asia", ["ITU", "STU"]),
    ("north_china", "Northern and Northeast China", "east_asia",
     ["CHB", "Han", "Tu", "Hezhen", "Oroqen", "Daur", "Xibo", "Mongola"]),
    ("south_china", "Southern China", "east_asia", ["CHS", "She", "Tujia", "Miao", "Yi", "Naxi"]),
    ("japan_korea", "Japan and Korea", "east_asia", ["JPT", "Japanese", "Korean"]),
    ("se_asia", "Mainland Southeast Asia", "east_asia",
     ["CDX", "KHV", "Dai", "Kinh_Vietnamese", "Cambodian", "Thai", "Burmese", "China_Lahu"]),
    ("north_asia", "Siberia", "north_asia",
     ["Yakut", "Yakut_1", "Yakut_2", "Even", "Buryat", "Ulchi", "Koryak", "Chukchi", "Itelmen",
      "Eskimo_ChaplinSireniki", "Eskimo_Naukan"]),
    ("papua", "Papua New Guinea, Melanesia and Australia", "oceania", ["Papuan", "Nasioi", "Australian"]),
    ("americas", "Indigenous Americas", "americas",
     ["Karitiana", "Surui", "Piapoco", "Pima", "Mayan", "Mixe", "Quechua", "Aymara", "Huichol", "Nahua", "Mexico_Zapotec",
      "Mixtec_1", "Mixtec_2", "Tsimshian", "Yukpa"]),
]
GROUP_NAMES = {  # 1000 Genomes population codes
    "CEU": "Utah residents with Northern and Western European ancestry", "GBR": "British in England and Scotland",
    "FIN": "Finnish in Finland", "IBS": "Iberian populations in Spain", "TSI": "Toscani in Italy",
    "YRI": "Yoruba in Ibadan, Nigeria", "ESN": "Esan in Nigeria", "GWD": "Gambian in Western Division, The Gambia",
    "MSL": "Mende in Sierra Leone", "LWK": "Luhya in Webuye, Kenya", "CHB": "Han Chinese in Beijing",
    "CHS": "Southern Han Chinese", "CDX": "Chinese Dai in Xishuangbanna", "JPT": "Japanese in Tokyo",
    "KHV": "Kinh in Ho Chi Minh City, Vietnam", "GIH": "Gujarati Indians in Houston, Texas",
    "PJL": "Punjabi in Lahore, Pakistan", "BEB": "Bengali in Bangladesh", "STU": "Sri Lankan Tamil in the UK",
    "ITU": "Indian Telugu in the UK", "ASW": "African Ancestry in Southwest USA", "ACB": "African Caribbean in Barbados",
    "MXL": "Mexican Ancestry in Los Angeles", "PUR": "Puerto Rican in Puerto Rico", "CLM": "Colombian in Medellín",
    "PEL": "Peruvian in Lima"}
PAR = [(60001, 2699520), (154931044, 155260560)]           # GRCh37 pseudo-autosomal regions of X

A = AADR(PREFIX, ANNO)
synthetic = A.nsnp < 500_000
MIN_REGION = 8 if synthetic else 15      # smaller panels estimate frequencies too noisily to be told apart
LOO_SNPS = 80_000
CAL_PER_REGION = 16
VAL_PER_REGION = 12
OWN_MIN = 0.65         # with margin: the validation set must still reach 0.6
log("AADR", A.nind, "individuals", A.nsnp, "SNPs", "(synthetic stand-in)" if synthetic else "")

keep = (A.chrom >= 1) & (A.chrom <= 22) & (A.a1 < 4) & (A.a2 < 4)
S_idx = np.flatnonzero(keep)
snps = np.load(PACK / "snps.npz")
assert len(S_idx) == len(snps["pos"]) and np.array_equal(A.pos[S_idx], snps["pos"]), "pack SNP table differs"
S = len(S_idx)
par = np.zeros(A.nsnp, bool)
for lo, hi in PAR:
    par |= (A.pos >= lo) & (A.pos <= hi)
X_idx = np.flatnonzero((A.chrom == 23) & (A.a1 < 4) & (A.a2 < 4) & ~par)
log("autosomal SNPs", S, "X SNPs (outside the pseudo-autosomal regions)", len(X_idx))

present = (A.date_bp == 0)
reg_members = {}
for rid, title, cont, groups in REGIONS:
    idx = A.best_per_individual([i for i in range(A.nind) if present[i] and A.group[i] in groups])
    if len(idx) < MIN_REGION:
        log(f"region {rid}: only {len(idx)} individuals; left out")
        continue
    reg_members[rid] = idx
regions = [r for r in REGIONS if r[0] in reg_members]
panel = sorted({i for v in reg_members.values() for i in v})
log("regions", len(regions), "panel individuals", len(panel))

rows_all = sorted(panel)
PK = A.rows(rows_all)
row_of = {i: k for k, i in enumerate(rows_all)}
log(f"genotype rows in memory: {PK.nbytes / 1e9:.2f} GB")

# a fixed random subset of SNPs for the leave-one-out runs (the low end of consumer overlap)
sub = np.sort(rng.choice(S, min(LOO_SNPS, S), replace=False))
DSUB = A.dosage(PK, S_idx[sub])                      # (rows_all, sub)
log("subset dosages", DSUB.shape)

# per-individual allele counts on the subset, and per-region totals on all SNPs
ind_alt_sub = np.where(DSUB >= 0, DSUB, 0).astype(np.int8)
ind_tot_sub = np.where(DSUB >= 0, 2, 0).astype(np.int8)


def region_counts(members: dict[str, list[int]], snp_idx: np.ndarray, chunk=100_000):
    rids = list(members)
    alt = np.zeros((len(rids), len(snp_idx)), np.int32)
    tot = np.zeros((len(rids), len(snp_idx)), np.int32)
    for c0 in range(0, len(snp_idx), chunk):
        s = snp_idx[c0:c0 + chunk]
        for k, r in enumerate(rids):
            D = A.dosage(PK[[row_of[i] for i in members[r]]], s)
            alt[k, c0:c0 + len(s)] = np.where(D >= 0, D, 0).sum(0)
            tot[k, c0:c0 + len(s)] = (D >= 0).sum(0) * 2
    return alt, tot


def freq(alt, tot):
    """Panel frequencies exactly as the app will see them (shrunk, then encoded)."""
    return decode_freq(encode_freq(panel_freqs(alt, tot)))


def loo_q(members: dict[str, list[int]], group_of_region: dict[str, int], ngroups: int, people: list[int]):
    """Leave-one-out proportions on the SNP subset for the listed people,
    with references = members grouped by group_of_region (region -> column)."""
    alt = np.zeros((ngroups, len(sub)), np.int64)
    tot = np.zeros((ngroups, len(sub)), np.int64)
    home = {}
    for r, idx in members.items():
        k = group_of_region[r]
        rr = [row_of[i] for i in idx]
        alt[k] += ind_alt_sub[rr].sum(0, dtype=np.int64)
        tot[k] += ind_tot_sub[rr].sum(0, dtype=np.int64)
        for i in idx:
            home[i] = k
    out = np.zeros((len(people), ngroups))
    for j, i in enumerate(people):
        a, t = alt.copy(), tot.copy()
        if i in home:
            a[home[i]] -= ind_alt_sub[row_of[i]]
            t[home[i]] -= ind_tot_sub[row_of[i]]
        P = freq(a, t)
        g = DSUB[row_of[i]]
        m = g >= 0
        out[j] = admix(g[m], P[:, m], tol=1e-4, max_iter=400)["q"]
    return out


cont_ids = [c for c, _ in CONTINENTS if any(r[2] == c for r in regions)]
cont_index = {c: k for k, c in enumerate(cont_ids)}

# ---------------------------------------------------------------- 2. outliers
reg_cont = {r[0]: cont_index[r[2]] for r in regions}
q = loo_q(reg_members, reg_cont, len(cont_ids), panel)
own = np.array([q[j, reg_cont[next(r for r, v in reg_members.items() if i in v)]] for j, i in enumerate(panel)])
by_group = {}
for j, i in enumerate(panel):
    by_group.setdefault(A.group[i], []).append(own[j])
med = {g: float(np.median(v)) for g, v in by_group.items()}
dropped = [(A.ids[i], str(A.group[i]), round(float(own[j]), 3)) for j, i in enumerate(panel)
           if own[j] < med[A.group[i]] - 0.15]
drop_ids = {d[0] for d in dropped}
for r in list(reg_members):
    reg_members[r] = [i for i in reg_members[r] if A.ids[i] not in drop_ids]
    if len(reg_members[r]) < MIN_REGION:
        log(f"region {r}: {len(reg_members[r])} left after outlier removal; left out")
        del reg_members[r]
regions = [r for r in regions if r[0] in reg_members]
log(f"outliers removed: {len(dropped)}", dropped[:40])

# ---------------------------------------------------------------- 3. merge regions that cannot be told apart
# Left-out individuals are split once: a calibration set (merging, the
# small-panel correction and detection floors) and a validation set that the
# build never looks at (validate_modern.py measures accuracy on it).
merges = []
titles: dict[str, list[str]] = {}
cal_sample, val_sample = {}, {}
for r, v in reg_members.items():
    perm = [int(x) for x in rng.permutation(v)]
    nc = min(CAL_PER_REGION, int(np.ceil(0.6 * len(perm))))
    cal_sample[r] = sorted(perm[:nc])
    val_sample[r] = sorted(perm[nc:nc + VAL_PER_REGION])


while True:
    rids = list(reg_members)
    col = {r: k for k, r in enumerate(rids)}
    people = [i for r in rids for i in cal_sample[r]]
    who = np.array([col[r] for r in rids for _ in cal_sample[r]])
    Q = loo_q(reg_members, col, len(rids), people)
    own = np.array([Q[who == k, k].mean() for k in range(len(rids))])
    conf = np.array([Q[who == k].mean(0) for k in range(len(rids))])
    reg_cont_now = {r[0]: r[2] for r in regions}
    log("own-region means:", {r: round(float(own[col[r]]), 2) for r in rids})
    worst = sorted([r for r in rids if own[col[r]] < OWN_MIN], key=lambda r: own[col[r]])
    cand = None
    for r in worst:
        sibs = [x for x in rids if x != r and reg_cont_now[x] == reg_cont_now[r]]
        if sibs:
            cand = (r, max(sibs, key=lambda x: conf[col[r], col[x]]))
            break
    if cand is None:
        break
    r, s_ = cand
    # merge the smaller into the larger; the merged region keeps the larger one's id and a combined title
    small, big = sorted((r, s_), key=lambda x: len(reg_members[x]))
    rb = next(x for x in regions if x[0] == big)
    rs_ = next(x for x in regions if x[0] == small)
    parts = titles.get(big, [rb[1]]) + titles.get(small, [rs_[1]])
    titles[big] = parts
    name = (" and ".join(parts) if len(parts) == 2 else ", ".join(parts[:-1]) + " and " + parts[-1]) if len(parts) <= 3 \
        else f"{dict(CONTINENTS)[rb[2]]} ({len(parts)} regions combined)"
    merged = (big, name, rb[2], rb[3] + rs_[3])
    regions = [merged if x[0] == big else x for x in regions if x[0] != small]
    reg_members[big] = sorted(reg_members[big] + reg_members[small])
    cal_sample[big] = sorted(cal_sample[big] + cal_sample[small])
    val_sample[big] = sorted(val_sample[big] + val_sample[small])
    del reg_members[small], cal_sample[small], val_sample[small]
    merges.append({"merged": small, "into": big, "own_share": round(float(own[col[r]]), 3)})
    log(f"merged {small} into {big} (own share of {r} was {own[col[r]]:.2f})")

rids = list(reg_members)
col = {r: k for k, r in enumerate(rids)}
reg_cont_now = {r[0]: r[2] for r in regions}
rc = np.array([cont_index[reg_cont_now[r]] for r in rids])
Qc = np.zeros((len(people), len(cont_ids)))
for k in range(len(rids)):
    Qc[:, rc[k]] += Q[:, k]
# "Clean" calibration people: assigned at least 90% to their own continent.
# Many reference groups are themselves mixed (Bedouins and Mozabites carry
# African ancestry, some Maya European ancestry), so people from those groups
# legitimately get shares elsewhere; the detection floors measure the noise
# seen in people whose ancestry is unambiguous.
clean = Qc[np.arange(len(people)), rc[who]] >= 0.9
floors = {}
for r in rids:
    k = col[r]
    others = Q[clean & (rc[who] != rc[k]), k]
    floors[r] = max(0.01, float(np.percentile(others, 99))) if len(others) else 0.01
cfloors = {}
for c in cont_ids:
    k = cont_index[c]
    others = Qc[clean & (rc[who] != k), k]
    cfloors[c] = max(0.01, float(np.percentile(others, 99))) if len(others) else 0.01
# typical result of reference people from each region (left out of their own panel)
PROFILE = np.array([Q[who == k].mean(0) for k in range(len(rids))])
log("detection floors:", {r: round(v, 3) for r, v in floors.items()}, {c: round(v, 3) for c, v in cfloors.items()})

# ---------------------------------------------------------------- final panel frequencies
alt, tot = region_counts(reg_members, S_idx)
P = panel_freqs(alt, tot)
n_ind = np.array([len(reg_members[r]) for r in rids])
log("panel frequencies done")

# X chromosome: men carry one copy (the AADR records it as homozygous)
male = {i: A.sex(i) == "M" for i in rows_all}
xalt = np.zeros((len(rids), len(X_idx)), np.int32)
xtot = np.zeros((len(rids), len(X_idx)), np.int32)
DX = A.dosage(PK, X_idx) if len(X_idx) else np.zeros((len(rows_all), 0), np.int8)
for k, r in enumerate(rids):
    for i in reg_members[r]:
        d = DX[row_of[i]]
        if male[i]:
            ok = (d == 0) | (d == 2)
            xalt[k] += np.where(ok, d // 2, 0)
            xtot[k] += ok
        else:
            ok = d >= 0
            xalt[k] += np.where(ok, d, 0)
            xtot[k] += 2 * ok
log("X frequencies done")

# ---------------------------------------------------------------- 5. painting calibration on simulated genomes
cm = snps["cm"].astype(np.float64)
chrom = snps["chrom"].astype(np.int64)


# hold out the sources from the continental frequencies used for calibration
cal_people = {}
for c in cont_ids:
    rs_ = [r for r in rids if reg_cont_now[r] == c]
    pool = [i for r in rs_ for i in cal_sample[r] if clean[people.index(i)]] or [i for r in rs_ for i in cal_sample[r]]
    cal_people[c] = [int(x) for x in rng.choice(pool, min(4, len(pool)), replace=False)]
hold = {i for v in cal_people.values() for i in v}
alt_c = np.zeros((len(cont_ids), S)); tot_c = np.zeros((len(cont_ids), S))
for k, r in enumerate(rids):
    c = cont_index[reg_cont_now[r]]
    alt_c[c] += alt[k]; tot_c[c] += tot[k]
full_cal = A.dosage(PK[[row_of[i] for i in sorted(hold)]], S_idx)
cal_geno = {i: full_cal[k] for k, i in enumerate(sorted(hold))}
for i, d in cal_geno.items():
    c = next(cc for cc, v in cal_people.items() if i in v)
    alt_c[cont_index[c]] -= np.where(d >= 0, d, 0)
    tot_c[cont_index[c]] -= (d >= 0) * 2
Pc_cal = freq(alt_c, tot_c)
pairs = [(a_, b_) for a_ in range(len(cont_ids)) for b_ in range(a_ + 1, len(cont_ids))]
rng.shuffle(pairs)
pairs = pairs[:6]
GENS = 10.0
grid = [0.25, 0.35, 0.5, 0.7, 1.0]
csub = np.sort(rng.choice(S, min(250_000, S), replace=False))
cal = {t: [] for t in grid}
for n_, (ca, cb) in enumerate(pairs):
    srcs = [cal_geno[cal_people[cont_ids[ca]][0]], cal_geno[cal_people[cont_ids[cb]][0]]]
    g, ta, tb = mosaic.build(srcs, [0.5, 0.5], chrom, cm, GENS, seed=100 + n_)
    for t in grid:
        r_ = paint(g[csub], Pc_cal[[ca, cb]][:, csub], chrom[csub], cm[csub], np.array([0.5, 0.5]),
                   generations=GENS, temper=t)
        acc = float(((r_["a"] == ta[csub]) & (r_["b"] == tb[csub])).mean())
        cal[t].append((acc, float(r_["posterior"].mean())))
summary_cal = {t: {"accuracy": float(np.mean([x[0] for x in v])), "mean_posterior": float(np.mean([x[1] for x in v]))}
               for t, v in cal.items()}
best_t = min(grid, key=lambda t: abs(summary_cal[t]["accuracy"] - summary_cal[t]["mean_posterior"])
             + max(0.0, max(s["accuracy"] for s in summary_cal.values()) - summary_cal[t]["accuracy"]))
log("painting calibration (temper -> accuracy, mean posterior):",
    {t: (round(v["accuracy"], 3), round(v["mean_posterior"], 3)) for t, v in summary_cal.items()}, "chosen", best_t)

# ---------------------------------------------------------------- runs of homozygosity in reference groups
roh_groups = {}
for i in range(A.nind):
    if present[i] and A.snps_hit[i] > 0:
        roh_groups.setdefault(A.group[i], []).append(i)
roh_groups = {g: A.best_per_individual(v) for g, v in roh_groups.items()}
roh_groups = {g: v for g, v in roh_groups.items() if len(v) >= 8 and "-o" not in g and "Discovery" not in g}
thin = np.sort(rng.choice(S, S // 2, replace=False))
pos = snps["pos"].astype(np.int64)
roh_ref = {}
for gname, v in sorted(roh_groups.items()):
    pick = sorted(int(x) for x in rng.choice(v, min(20, len(v)), replace=False))
    D = A.dosage(A.rows(pick), S_idx[thin])
    res = []
    for d in D:
        segs = roh(chrom[thin], pos[thin], d == 1, d < 0)
        res.append(roh_summary(segs))
    fr = np.array([x["froh"] for x in res])
    lm = np.array([x["by_length_mb"]["10+"] for x in res])
    reg = next((r[1] for r in regions if gname in r[3]), None)
    roh_ref[gname] = {"n": len(res), "region": reg, "label": GROUP_NAMES.get(gname, gname),
                      "froh_median": float(np.median(fr)), "froh_p90": float(np.percentile(fr, 90)),
                      "long_mb_median": float(np.median(lm)), "long_mb_p90": float(np.percentile(lm, 90))}
(PACK / "roh_reference.json").write_text(json.dumps({"snps_used": int(len(thin)), "groups": roh_ref}, indent=1))
log("ROH reference groups", len(roh_ref))

# ---------------------------------------------------------------- write
xs = X_idx
xblock = blocks_from_cm(A.chrom[xs], A.cm[xs], 5.0) if len(xs) else np.zeros(0, np.int32)
np.savez_compressed(PACK / "modern.npz", region_ids=np.array(rids), P=encode_freq(P), n_individuals=n_ind,
                    profile=PROFILE.astype(np.float32),
                    continent=np.array([cont_index[reg_cont_now[r]] for r in rids], np.int8),
                    continent_ids=np.array(cont_ids),
                    x_pos=A.pos[xs].astype(np.int32), x_ref=A.a1[xs], x_alt=A.a2[xs], x_rsid=A.rsid[xs],
                    x_cm=A.cm[xs].astype(np.float32), x_block=xblock.astype(np.int32),
                    Px=encode_freq(panel_freqs(xalt, xtot)))
manifest = json.loads((PACK / "manifest.json").read_text())
manifest["modern"] = {
    "regions": [{"id": r, "title": next(x[1] for x in regions if x[0] == r), "continent": reg_cont_now[r],
                 "groups": sorted({A.group[i] for i in reg_members[r]}),
                 "group_names": {g: GROUP_NAMES[g] for g in sorted({A.group[i] for i in reg_members[r]}) if g in GROUP_NAMES},
                 "individuals": len(reg_members[r]), "detection_floor": floors[r],
                 "own_share_loo": float(own[col[r]])}
                for r in rids],
    "continents": [{"id": c, "title": dict(CONTINENTS)[c], "detection_floor": cfloors[c]} for c in cont_ids],
    "outliers_removed": len(dropped), "merges": merges,
    "loo": {"snps": int(len(sub)), "calibration_per_region": CAL_PER_REGION, "validation_per_region": VAL_PER_REGION,
            "own_share_threshold": OWN_MIN},
    "painting": {"generations": GENS, "temper": best_t, "calibration": summary_cal},
    "x_snps": int(len(xs)), "not_covered": [
        "Central Asia (its populations are themselves mixtures of West and East Eurasian ancestry)",
        "Iran (too few reference genomes)", "Northeast Africa and the Horn of Africa (too few reference genomes)",
        "Polynesia and Island Southeast Asia (too few reference genomes)",
        "Indigenous North America (few reference genomes; most Indigenous American references are from Central and South America)",
        "Jewish diaspora groups (not in the reference set)"]}
(PACK / "manifest.json").write_text(json.dumps(manifest, indent=1))

VALDIR.mkdir(parents=True, exist_ok=True)
loo_people = [i for r in rids for i in val_sample[r]]
G_loo = A.dosage(PK[[row_of[i] for i in loo_people]], S_idx)
np.savez_compressed(VALDIR / "modern_val.npz",
                    loo_ids=np.array([A.ids[i] for i in loo_people]), loo_region=np.array([col[r] for r in rids for _ in val_sample[r]]),
                    loo_group=np.array([A.group[i] for i in loo_people]),
                    loo_male=np.array([male[i] for i in loo_people]),
                    loo_geno=G_loo, loo_x=DX[[row_of[i] for i in loo_people]],
                    alt=alt.astype(np.int32), tot=tot.astype(np.int32), xalt=xalt, xtot=xtot)
log("done")
