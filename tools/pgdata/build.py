"""Build the Personal Genome reference data pack (runs in CI).

Inputs (downloaded by the workflow, checksums verified there):
  AADR v66.p1 1240K: .anno .ind .snp .geno (TGENO, individual-major)
  Yleaf (YFull tree and GRCh37 SNP positions), PhyloTree 17 RSRS (Haplogrep)
  UCSC hg38ToHg19 chain

Usage: python build.py AADR_PREFIX ANNO YLEAF_DIR PHYLOTREE_DIR CHAIN OUT_DIR [VALIDATION_DIR]

Writes the pack files described in core/athanor_core/personal/refdata.py,
plus (optionally) a validation bundle of held-out genotypes used by
validate.py. All statistics use the same code as the app
(athanor_core.personal).
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.personal.fstats import blocks_from_cm  # noqa: E402
from athanor_core.personal.haplo import Tree, build_mt_tree, build_y_tree  # noqa: E402

PREFIX, ANNO, YLEAF, PHYLO, CHAIN, OUT = [Path(a) for a in sys.argv[1:7]]
VALDIR = Path(sys.argv[7]) if len(sys.argv) > 7 else None
OUT.mkdir(parents=True, exist_ok=True)
T0 = time.time()
K = 8                       # principal components kept
CHUNK = 20000               # SNPs per pass chunk


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ------------------------------------------------------------------ SNPs and individuals
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
snp_lines = (PREFIX.with_suffix(".snp")).read_text().splitlines()
nsnp_all = len(snp_lines)
chrom = np.zeros(nsnp_all, np.int16); pos = np.zeros(nsnp_all, np.int64); cm = np.zeros(nsnp_all)
a1 = np.full(nsnp_all, 255, np.uint8); a2 = np.full(nsnp_all, 255, np.uint8); rsid = np.full(nsnp_all, -1, np.int64)
for i, l in enumerate(snp_lines):
    f = l.split()
    chrom[i] = int(f[1]); cm[i] = float(f[2]) * 100; pos[i] = int(f[3])
    a1[i] = BASE.get(f[4], 255); a2[i] = BASE.get(f[5], 255)
    if f[0].startswith("rs") and f[0][2:].isdigit():
        rsid[i] = int(f[0][2:])
keep = (chrom >= 1) & (chrom <= 22) & (a1 < 4) & (a2 < 4)
ysnp = (chrom == 24) & (a1 < 4) & (a2 < 4)
log("SNPs", nsnp_all, "autosomal biallelic", int(keep.sum()), "Y", int(ysnp.sum()))

ind = [l.split() for l in (PREFIX.with_suffix(".ind")).read_text().splitlines()]
nind = len(ind)
ids = [x[0] for x in ind]
group = np.array([x[2] for x in ind])
rows = list(csv.reader(open(ANNO, encoding="utf-8", errors="replace"), delimiter="\t"))
hdr = rows[0]


def col(pred):
    return next(i for i, h in enumerate(hdr) if pred(h.lower()))


C = {"gid": 0, "iid": 2, "pub": 6, "doi": 7, "date": col(lambda h: h.startswith("date mean in bp")),
     "sd": col(lambda h: h.startswith("date standard deviation")), "group": col(lambda h: h.startswith("group id")),
     "loc": col(lambda h: h == "locality"), "country": col(lambda h: h.startswith("political entity")),
     "lat": col(lambda h: h == "latitude"), "lon": col(lambda h: h == "longitude"),
     "snps": col(lambda h: "1240k snpset" in h), "sex": col(lambda h: h.startswith("molecular sex")),
     "y": col(lambda h: h.startswith("y haplogroup in terminal")), "y_isogg": col(lambda h: h.startswith("y haplogroup  in isogg") or h.startswith("y haplogroup in isogg")),
     "y_manual": col(lambda h: h.startswith("y haplogroup manually")), "mt": col(lambda h: h.startswith("mtdna haplogroup")),
     "assess": col(lambda h: h == "assessment"), "full_date": col(lambda h: h.startswith("full date"))}
anno = {r[0]: r for r in rows[1:] if r}
log("individuals", nind, "annotated", len(anno))


def a_get(gid, key):
    r = anno.get(gid)
    if not r or C[key] >= len(r):
        return None
    v = r[C[key]].strip()
    return None if v in ("", "..", "n/a") else v


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


snps_hit = np.array([num(a_get(g, "snps")) or 0 for g in ids])
date_bp = np.array([num(a_get(g, "date")) if a_get(g, "date") is not None else np.nan for g in ids])

# ------------------------------------------------------------------ genotype access (TGENO)
geno_path = PREFIX.with_suffix(".geno")
with open(geno_path, "rb") as fh:
    head = fh.read(48).split(b"\x00")[0].decode().split()
assert head[0] == "TGENO" and int(head[1]) == nind and int(head[2]) == nsnp_all, head
RLEN = (nsnp_all + 3) // 4
assert geno_path.stat().st_size == 48 + nind * RLEN, "unexpected TGENO layout"
MM = np.memmap(geno_path, dtype=np.uint8, mode="r", offset=48, shape=(nind, RLEN))
LUT = np.array([[(b >> (6 - 2 * k)) & 3 for k in range(4)] for b in range(256)], dtype=np.uint8)


def decode(rows_idx: np.ndarray, s0: int, s1: int, source=None) -> np.ndarray:
    """(len(rows), s1 - s0) genotypes = copies of allele 1, 3 = missing."""
    b0, b1 = s0 // 4, (s1 + 3) // 4
    raw = np.asarray((MM if source is None else source)[rows_idx, b0:b1])
    g = LUT[raw].reshape(len(rows_idx), -1)
    off = s0 - b0 * 4
    return g[:, off:off + (s1 - s0)]


def load_rows(rows: list) -> tuple[np.ndarray, dict]:
    """Copy whole genotype rows into memory (each row is contiguous on disk).
    Pass 1 touches the same individuals for every SNP chunk, so reading them
    once this way replaces millions of scattered reads."""
    arr = np.empty((len(rows), RLEN), np.uint8)
    for k, i in enumerate(rows):
        arr[k] = MM[i]
    return arr, {i: k for k, i in enumerate(rows)}


# orientation check: the hg19 reference individual must be homozygous for allele 1
ref_ind = [i for i, gname in enumerate(group) if gname == "hg19ref"]
if ref_ind:
    gtest = decode(np.array(ref_ind[:1]), 0, 200000)[0]
    called = gtest != 3
    frac2 = float((gtest[called] == 2).mean())
    log("hg19ref homozygous for allele 1:", round(frac2, 4))
    assert frac2 > 0.95, "genotype coding is not 'copies of allele 1'"

# ------------------------------------------------------------------ population roles
def members(groups: list[str], best_per_individual=True, min_snps=0) -> list[int]:
    idx = [i for i, gname in enumerate(group) if gname in groups and snps_hit[i] >= min_snps]
    if best_per_individual:          # one version per individual (highest coverage)
        best = {}
        for i in idx:
            iid = a_get(ids[i], "iid") or ids[i]
            if iid not in best or snps_hit[i] > snps_hit[best[iid]]:
                best[iid] = i
        idx = sorted(best.values())
    return idx


def single_best(groups):
    idx = members(groups)
    return [max(idx, key=lambda i: snps_hit[i])] if idx else []


def regex_groups(pat):
    return sorted({gname for gname in set(group) if re.search(pat, gname) and not re.search(r"[-_]o($|[A-Z_-])|contam|_d$|lc$", gname)})


ROLES = {
    "Mbuti": members(["Mbuti", "MbutiPygmy"]),
    "YRI": members(["YRI"]),
    "Chimp": single_best(["Chimp"]),
    "Altai": single_best(["Altai_Neanderthal"]),
    "Vindija": single_best(["VindijaCave_Neanderthal"]),
    "Chagyrskaya": single_best(["ChagyrskayaCave_Neanderthal"]),
    "Denisova": single_best(["Russia_Denisova"]),
    # Western hunter-gatherers: the western European Mesolithic groups. The
    # eastern (Karelia, Latvia, Ukraine) and Scandinavian groups carry eastern
    # hunter-gatherer ancestry as well, and the Iron Gates groups carry early
    # farmer ancestry, so none of those belong in this source.
    "WHG": members(["Luxembourg_Loschbour_Mesolithic", "France_Mesolithic", "Spain_Mesolithic",
                    "Belgium_Mesolithic", "England_Mesolithic", "Netherlands_Mesolithic",
                    "Germany_GrossFredenwalde_Mesolithic"], min_snps=30000),
    "Anatolia_N": members(["Turkey_N"]),
    "Yamnaya": members(["Russia_Samara_EBA_Yamnaya"]),
    "Iran_N": members(["Iran_GanjDareh_N"]),
    "CHG": members(["Georgia_KotiasKlde_Mesolithic", "Georgia_Satsurblia_LateUP"]),
    "Steppe_MLBA": members(["Russia_Chelyabinsk_MLBA_Sintashta"]),
    "Onge": members(["ONG"]),
    "UstIshim": members(["Russia_UstIshim_IUP"]),
    "Kostenki": members(["Russia_Kostenki_UP"]),
    "MA1": members(["Russia_Malta_UP"]),
    "Han": members(["Han"]),
    "Papuan": members(["Papuan"]),
    "Chukchi": members(["Chukchi"]),
    "Karitiana": members(["Karitiana"]),
    "Natufian": members(["Israel_Natufian"]),
    "Mota": members(["Ethiopia_MotaCave_4500BP"]),
}
MIN_MEMBERS = {"WHG": 15, "Anatolia_N": 15, "Yamnaya": 10, "Iran_N": 5, "YRI": 20, "Mbuti": 5}
for k, v in ROLES.items():
    log(f"role {k}: {len(v)} individuals: {sorted({group[i] for i in v})} e.g. {[ids[i] for i in v[:3]]}")
    assert v, f"no individuals for {k}"
    assert len(v) >= MIN_MEMBERS.get(k, 1), f"{k}: only {len(v)} individuals, expected at least {MIN_MEMBERS[k]}"
    assert len(v) <= 127, f"too many individuals for uint8 counts in {k}"

MODELS = [
    {"id": "europe", "title": "European three-way model",
     "description": "Western hunter-gatherers, Anatolian early farmers and Yamnaya steppe herders (Haak et al. 2015).",
     "sources": ["WHG", "Anatolia_N", "Yamnaya"],
     "labels": {"WHG": "Western hunter-gatherer", "Anatolia_N": "Anatolian early farmer", "Yamnaya": "Steppe herder (Yamnaya)"},
     "right": ["Mbuti", "UstIshim", "Kostenki", "MA1", "Han", "Papuan", "Onge", "Chukchi", "Karitiana"]},
    {"id": "europe_west_asia", "title": "European and West Asian four-way model",
     "description": "Adds early Iranian farmers for ancestry from the Caucasus and Near East (southern Europe, West Asia).",
     "sources": ["WHG", "Anatolia_N", "Yamnaya", "Iran_N"],
     "labels": {"WHG": "Western hunter-gatherer", "Anatolia_N": "Anatolian early farmer", "Yamnaya": "Steppe herder (Yamnaya)",
                "Iran_N": "Iranian early farmer"},
     "right": ["Mbuti", "UstIshim", "Kostenki", "MA1", "Han", "Papuan", "Onge", "Chukchi", "Karitiana", "Natufian", "Mota"]},
    {"id": "south_asia", "title": "South Asian three-way model",
     "description": "Iranian-related farmers, Middle-to-Late Bronze Age steppe herders and an ancestral South Indian proxy (Narasimhan et al. 2019).",
     "sources": ["Iran_N", "Steppe_MLBA", "Onge"],
     "labels": {"Iran_N": "Iranian-related farmer", "Steppe_MLBA": "Steppe herder (Bronze Age)", "Onge": "Ancestral South Indian (Onge proxy)"},
     "right": ["Mbuti", "UstIshim", "Kostenki", "MA1", "Han", "Papuan", "Karitiana", "Natufian", "Anatolia_N", "WHG", "Mota"]},
]
ARCH = {"altai": "Altai", "chimp": "Chimp", "vindija": "Vindija", "africa": "YRI", "denisova": "Denisova"}

# ------------------------------------------------------------------ PCA reference sets
REGION = {}
for p in ("CEU", "GBR", "FIN", "IBS", "TSI", "French", "Sardinian", "Basque", "Italian_North", "Orcadian", "Russian",
          "English", "Spanish", "Scottish", "Finnish", "Icelandic", "Norwegian", "Polish", "Czech", "Hungarian",
          "Estonian", "Bulgarian", "Greek", "Albanian", "Saami"):
    REGION[p] = "Europe"
for p in ("Adygei", "Druze", "Palestinian", "BedouinA", "BedouinB", "Jordanian", "Iranian", "Armenian", "Georgian",
          "Turkish", "Abkhasian", "Lezgin", "Chechen", "Samaritan", "Jew_Yemenite", "Jew_Iraqi", "Ossetian"):
    REGION[p] = "West Asia"
for p in ("Mozabite", "Saharawi"):
    REGION[p] = "North Africa"
for p in ("YRI", "LWK", "GWD", "MSL", "ESN", "Yoruba", "Mandenka", "Biaka", "Mbuti", "BantuKenya", "Luhya",
          "Mende", "Gambian", "Ju_hoan_North", "BantuSA", "Dinka", "Masai", "Luo", "Igbo", "Khomani_San", "Lemande"):
    REGION[p] = "Africa"
for p in ("GIH", "PJL", "BEB", "STU", "ITU", "Pakistani", "Brahui", "Balochi", "Makrani", "Sindhi_Pakistan",
          "Pathan", "Burusho", "Kalash", "Punjabi", "Brahmin", "Irula", "Kapu", "Madiga", "Mala", "Yadava", "Kusunda"):
    REGION[p] = "South Asia"
for p in ("Hazara", "Uyghur", "Tajik", "Uzbek", "Kyrgyz_Kyrgyzstan", "Kazakh", "Turkmen", "Karakalpak"):
    REGION[p] = "Central Asia"
for p in ("CHB", "CHS", "CDX", "JPT", "KHV", "Han", "Japanese", "Dai", "Miao", "She", "Tu", "Yi", "Naxi", "Tujia",
          "Cambodian", "Kinh_Vietnamese", "Korean", "Thai", "Burmese", "Ami", "Atayal", "Lahu", "China_Lahu"):
    REGION[p] = "East Asia"
for p in ("Yakut", "Hezhen", "Mongola", "Oroqen", "Daur", "Xibo", "Buryat", "Even", "Chukchi", "Koryak", "Itelmen",
          "Ulchi", "Aleut", "Altaian", "Tubalar", "Mansi"):
    REGION[p] = "Siberia and North Asia"
for p in ("Papuan", "Nasioi", "Australian", "Dusun", "Igorot", "Hawaiian", "Maori"):
    REGION[p] = "Oceania and Island Southeast Asia"
for p in ("Karitiana", "Surui", "Pima", "Mayan", "Mixe", "Piapoco", "Quechua", "Aymara", "Huichol", "Mexico_Zapotec",
          "Nahua", "Tsimshian", "Yukpa"):
    REGION[p] = "Americas"
ADMIXED = {"ASW", "ACB", "PUR", "CLM", "MXL", "PEL"}

present = [i for i in range(nind) if date_bp[i] == 0 and group[i] in REGION
           and not re.search(r"Discovery", group[i])]
rng = np.random.default_rng(20261022)


def capped(idx, cap):
    by = {}
    for i in idx:
        by.setdefault(group[i], []).append(i)
    out = []
    for gname, v in sorted(by.items()):
        v = sorted(v)
        out += list(rng.choice(v, cap, replace=False)) if len(v) > cap else v
    return sorted(int(x) for x in out)


# held-out individuals for validation (not used to build the PCA or population frequencies)
VAL_POPS = ["GBR", "CEU", "FIN", "IBS", "TSI", "YRI", "CHB", "JPT", "GIH", "PJL", "PEL", "MXL", "PUR", "Sardinian",
            "Druze", "Papuan", "Karitiana"]
holdout = []
for p in VAL_POPS:
    v = [i for i in range(nind) if group[i] == p and date_bp[i] == 0]
    used = set(i for r in ROLES.values() for i in r)
    v = [i for i in v if i not in used]
    holdout += sorted(int(x) for x in rng.choice(v, min(16, len(v)), replace=False)) if v else []
hold = set(holdout)
world_ref = capped([i for i in present if i not in hold], 25)
west_ref = capped([i for i in present if i not in hold and REGION[group[i]] in ("Europe", "West Asia", "North Africa")], 40)
log("PCA reference individuals: world", len(world_ref), "west", len(west_ref), "held out", len(holdout))

# ancient individuals to describe and project
assess = {i: (a_get(ids[i], "assess") or "") for i in range(nind)}
anc_all = [i for i in range(nind) if not np.isnan(date_bp[i]) and date_bp[i] > 0 and snps_hit[i] >= 20000
           and assess[i].upper().startswith("PASS") and group[i] not in ("Chimp", "Gorilla", "Ancestor")
           and not re.search(r"Neanderthal|Denisova", group[i])]          # archaic hominins are not placed on modern-human PCAs
best = {}
for i in anc_all:
    iid = a_get(ids[i], "iid") or ids[i]
    if iid not in best or snps_hit[i] > snps_hit[best[iid]]:
        best[iid] = i
ancients = sorted(best.values())
log("ancient individuals (PASS, >= 20k SNPs, one version each)", len(ancients))

# Y validation set: ancient males with Y calls and good coverage
ymales = [i for i in ancients if (a_get(ids[i], "sex") or "") == "M" and a_get(ids[i], "y") and snps_hit[i] >= 500000][:600]

# ------------------------------------------------------------------ pass 1
S_idx = np.flatnonzero(keep)
S = len(S_idx)
role_names = list(ROLES)
alt_cnt = np.zeros((len(role_names), S), np.uint8)
tot_cnt = np.zeros((len(role_names), S), np.uint8)
SUPER = {"EUR": ["CEU", "GBR", "FIN", "IBS", "TSI"], "EAS": ["CHB", "CHS", "CDX", "JPT", "KHV"],
         "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"], "AMR": ["MXL", "PUR", "CLM", "PEL"],
         "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN"], "AFR_all": ["YRI", "LWK", "GWD", "MSL", "ESN", "Mbuti", "Yoruba",
                                                              "Mandenka", "Biaka", "Luhya", "Mende", "Gambian", "Ju_hoan_North",
                                                              "BantuKenya", "BantuSA", "Dinka", "Masai", "Luo", "Igbo"]}
super_idx = {k: [i for i in range(nind) if group[i] in v and date_bp[i] == 0] for k, v in SUPER.items()}
super_alt = {k: np.zeros(S, np.int32) for k in SUPER}
super_tot = {k: np.zeros(S, np.int32) for k in SUPER}

pca_sets = {"world": np.array(world_ref), "west": np.array(west_ref)}
pca_stats = {k: {"p": np.zeros(S, np.float32), "ok": np.zeros(S, bool)} for k in pca_sets}
grm = {k: np.zeros((len(v), len(v))) for k, v in pca_sets.items()}

present_all = [i for i in range(nind) if date_bp[i] == 0 and group[i] in REGION]
blocks_all = blocks_from_cm(chrom[S_idx], cm[S_idx], 5.0)
NB = int(blocks_all.max()) + 1
indstat = {k: np.zeros((len(present_all), NB)) for k in ("nean_num", "nean_den", "den_num", "den_den")}

need = sorted(set(int(i) for r in ROLES.values() for i in r) | set(world_ref) | set(west_ref)
              | set(i for v in super_idx.values() for i in v) | set(present_all))
NEED, pos_in_need = load_rows(need)
need_arr = np.arange(len(need))
log(f"held {len(need)} genotype rows in memory ({NEED.nbytes / 1e9:.2f} GB)")


def freq_of(G, idx):
    sub = G[[pos_in_need[i] for i in idx]]
    called = sub != 3
    alt = np.where(called, 2 - sub, 0).sum(0)
    tot = 2 * called.sum(0)
    return alt, tot


log("pass 1: frequencies, archaic statistics, PCA relationship matrices")
for c0 in range(0, S, CHUNK):
    c1 = min(S, c0 + CHUNK)
    s_lo, s_hi = int(S_idx[c0]), int(S_idx[c1 - 1]) + 1
    G = decode(need_arr, s_lo, s_hi, NEED)[:, S_idx[c0:c1] - s_lo]    # (len(need), chunk)
    for k, r in enumerate(role_names):
        a, t = freq_of(G, ROLES[r])
        alt_cnt[k, c0:c1] = a; tot_cnt[k, c0:c1] = t
    for k, v in super_idx.items():
        a, t = freq_of(G, v)
        super_alt[k][c0:c1] = a; super_tot[k][c0:c1] = t
    # per-individual archaic statistics (present-day)
    fr = {}
    for r in ("Altai", "Chimp", "Vindija", "YRI", "Denisova"):
        k = role_names.index(r)
        t = tot_cnt[k, c0:c1].astype(float)
        fr[r] = np.where(t > 0, alt_cnt[k, c0:c1] / np.maximum(t, 1), np.nan)
    sub = G[[pos_in_need[i] for i in present_all]]
    x = np.where(sub != 3, (2 - sub) / 2.0, np.nan)                 # (individuals, chunk)
    blk = blocks_all[c0:c1]
    A_, C_, V_, Y_, D_ = fr["Altai"], fr["Chimp"], fr["Vindija"], fr["YRI"], fr["Denisova"]
    base_ok = ~(np.isnan(A_) | np.isnan(C_) | np.isnan(V_) | np.isnan(Y_) | np.isnan(D_))
    ok = base_ok[None, :] & ~np.isnan(x)
    xn = np.where(ok, x, 0.0)
    A0, C0, V0, Y0, D0 = (np.nan_to_num(v) for v in (A_, C_, V_, Y_, D_))
    terms = {"nean_num": (A0 - C0)[None, :] * (xn - Y0[None, :]),
             "nean_den": np.broadcast_to(((A0 - C0) * (V0 - Y0))[None, :], xn.shape),
             "den_num": (Y0[None, :] - xn) * (D0 - V0)[None, :],
             "den_den": (Y0[None, :] + xn - 2 * Y0[None, :] * xn) * (D0 + V0 - 2 * D0 * V0)[None, :]}
    starts = np.flatnonzero(np.r_[True, blk[1:] != blk[:-1]])
    ub = blk[starts]
    for key, arr in terms.items():
        indstat[key][:, ub] += np.add.reduceat(np.where(ok, arr, 0.0), starts, axis=1)
    # PCA: allele frequency within the reference set, then the GRM
    for name, members_ in pca_sets.items():
        sub = G[[pos_in_need[i] for i in members_]].astype(np.float32)
        called = sub != 3
        n = called.sum(0)
        p = np.where(n > 0, np.where(called, 2 - sub, 0).sum(0) / np.maximum(2 * n, 1), 0)
        ok = (n >= 0.95 * len(members_)) & (p >= 0.01) & (p <= 0.99)
        pca_stats[name]["p"][c0:c1] = p; pca_stats[name]["ok"][c0:c1] = ok
        Z = np.where(called, ((2 - sub) - 2 * p) / np.sqrt(np.maximum(p * (1 - p), 1e-6)), 0)[:, ok]
        grm[name] += Z @ Z.T
    if (c0 // CHUNK) % 10 == 0:
        log(f"  {c1}/{S}")

eig = {}
for name in pca_sets:
    w, U = np.linalg.eigh(grm[name])
    o = np.argsort(w)[::-1][:K]
    eig[name] = (w[o], U[:, o])
    log(f"PCA {name}: top eigenvalues {np.round(w[o] / w.sum() * 100, 2)} % of variance")

# ------------------------------------------------------------------ pass 2: loadings and projections
loadings = {k: np.zeros((S, K), np.float32) for k in pca_sets}
log("pass 2a: SNP loadings")
for name, members_ in pca_sets.items():
    lam, U = eig[name]
    ok = pca_stats[name]["ok"]
    p_all = pca_stats[name]["p"]
    sc_all = np.sqrt(np.maximum(p_all * (1 - p_all), 1e-6))
    acc = np.zeros((S, K), np.float64)
    for c0 in range(0, S, CHUNK):
        c1 = min(S, c0 + CHUNK)
        s_lo, s_hi = int(S_idx[c0]), int(S_idx[c1 - 1]) + 1
        G = decode(np.array([pos_in_need[i] for i in members_]), s_lo, s_hi, NEED)[:, S_idx[c0:c1] - s_lo]
        G = G.astype(np.float32)
        p, sc = p_all[c0:c1], sc_all[c0:c1]
        Z = np.where(G != 3, ((2 - G) - 2 * p) / sc, 0)
        acc[c0:c1] = (Z.T @ U) / np.sqrt(lam)[None, :]
    acc[~ok] = 0
    loadings[name] = acc.astype(np.float32)
    log(f"  {name} done")

proj_rows = sorted(set(ancients) | set(holdout) | set(present_all))
AtA = {k: np.zeros((len(proj_rows), K * K)) for k in pca_sets}
Atb = {k: np.zeros((len(proj_rows), K)) for k in pca_sets}
log("pass 2b: projecting", len(proj_rows), "individuals")
BATCH = 64
SUB = 200000
for p0 in range(0, len(proj_rows), BATCH):
    rows = proj_rows[p0:p0 + BATCH]
    block_rows = np.empty((len(rows), RLEN), np.uint8)
    for k, i in enumerate(rows):
        block_rows[k] = MM[i]                       # one contiguous read per individual
    idx = np.arange(len(rows))
    for name in pca_sets:
        ok = pca_stats[name]["ok"]
        p_all = pca_stats[name]["p"]
        sc_all = np.sqrt(np.maximum(p_all * (1 - p_all), 1e-6))
        V = loadings[name]
        for c0 in range(0, S, SUB):
            c1 = min(S, c0 + SUB)
            s_lo, s_hi = int(S_idx[c0]), int(S_idx[c1 - 1]) + 1
            G = decode(idx, s_lo, s_hi, block_rows)[:, S_idx[c0:c1] - s_lo]
            m = (G != 3) & ok[c0:c1][None, :]
            Zp = np.where(m, ((2 - G.astype(np.float32)) - 2 * p_all[c0:c1]) / sc_all[c0:c1], 0)
            Vc = V[c0:c1]
            Atb[name][p0:p0 + len(rows)] += Zp @ Vc
            VV = (Vc[:, :, None] * Vc[:, None, :]).reshape(len(Vc), K * K)
            AtA[name][p0:p0 + len(rows)] += m.astype(np.float32) @ VV
    if (p0 // BATCH) % 20 == 0:
        log(f"  {p0 + len(rows)}/{len(proj_rows)}")

proj = {}
for name in pca_sets:
    P = np.full((len(proj_rows), K), np.nan)
    for j in range(len(proj_rows)):
        A = AtA[name][j].reshape(K, K)
        if A[0, 0] > 0 and np.linalg.cond(A) < 1e10:
            P[j] = np.linalg.solve(A, Atb[name][j])
    proj[name] = {i: P[k] for k, i in enumerate(proj_rows)}

# ------------------------------------------------------------------ trees
seq = "".join(l.strip() for l in open(PHYLO / "src/rsrs.fasta") if not l.startswith(">"))
mt = build_mt_tree((PHYLO / "src/tree.xml").read_text(), seq, (PHYLO / "src/weights.txt").read_text())
yt = build_y_tree(json.load(open(YLEAF / "yleaf/data/hg_prediction_tables/tree.json")),
                  (YLEAF / "yleaf/data/hg19/new_positions.txt").read_text())
MT, YT = Tree(mt), Tree(yt)
log("trees: mt", len(mt["nodes"]), "y", len(yt["nodes"]))

# ------------------------------------------------------------------ write the pack
s_sub = S_idx
np.savez_compressed(OUT / "snps.npz", chrom=chrom[s_sub].astype(np.int8), pos=pos[s_sub].astype(np.int32),
                    ref=a1[s_sub], alt=a2[s_sub], rsid=rsid[s_sub], cm=cm[s_sub].astype(np.float32),
                    block=blocks_all.astype(np.int32))
np.savez_compressed(OUT / "freqs.npz", names=np.array(role_names), alt=alt_cnt, total=tot_cnt)
for name, members_ in pca_sets.items():
    lam, U = eig[name]
    ok = pca_stats[name]["ok"]
    snp_idx = np.flatnonzero(ok).astype(np.int32)
    p = pca_stats[name]["p"][ok]
    np.savez_compressed(OUT / f"pca_{name}.npz", snp_idx=snp_idx, p=p.astype(np.float32),
                        scale=np.sqrt(np.maximum(p * (1 - p), 1e-6)).astype(np.float32),
                        loadings=loadings[name][ok].astype(np.float16), eigvals=lam,
                        ref_coords=(U * np.sqrt(lam)[None, :]).astype(np.float32),
                        ref_group=np.array([group[i] for i in members_]),
                        ref_region=np.array([REGION[group[i]] for i in members_]),
                        ref_id=np.array([ids[i] for i in members_]))
    log(f"pca_{name}: {len(snp_idx)} SNPs")

anc_out, mt_matched, y_matched = [], 0, 0
for i in ancients:
    gid = ids[i]
    mt_name = a_get(gid, "mt")
    y_name = a_get(gid, "y_manual") or a_get(gid, "y")
    mt_node = MT.lookup(mt_name) if mt_name else None
    y_node = YT.lookup(y_name, trim=False) if y_name else None
    mt_matched += mt_node is not None
    y_matched += y_node is not None
    anc_out.append({"id": gid, "group": group[i], "date_bp": int(date_bp[i]), "date_sd": num(a_get(gid, "sd")),
                    "full_date": a_get(gid, "full_date"), "locality": a_get(gid, "loc"), "country": a_get(gid, "country"),
                    "lat": num(a_get(gid, "lat")), "lon": num(a_get(gid, "lon")), "sex": a_get(gid, "sex"),
                    "mt": mt_name, "y": y_name, "y_isogg": a_get(gid, "y_isogg"), "mt_node": mt_node, "y_node": y_node,
                    "publication": a_get(gid, "pub"), "doi": a_get(gid, "doi"), "snps": int(snps_hit[i]),
                    "pc": {name: (None if np.isnan(proj[name][i]).any() else np.round(proj[name][i], 5).tolist())
                           for name in pca_sets}})
log("ancients written", len(anc_out), "mt matched", mt_matched, "y matched", y_matched)
with gzip.open(OUT / "ancients.json.gz", "wt") as fh:
    json.dump(anc_out, fh)
with gzip.open(OUT / "mt_tree.json.gz", "wt") as fh:
    json.dump(mt, fh)
with gzip.open(OUT / "y_tree.json.gz", "wt") as fh:
    json.dump(yt, fh)
shutil.copy(CHAIN, OUT / "hg38ToHg19.over.chain.gz")

# archaic alleles: all three Neanderthals carry it (homozygous where called; at least two called),
# under 1% in sub-Saharan Africans (>= 100 chromosomes), and present outside Africa
nea = [role_names.index(r) for r in ("Altai", "Vindija", "Chagyrskaya")]
na_alt = alt_cnt[nea].astype(int); na_tot = tot_cnt[nea].astype(int)
afr_f = super_alt["AFR_all"] / np.maximum(super_tot["AFR_all"], 1)
out_alt = sum(super_alt[k] for k in ("EUR", "EAS", "SAS")); out_tot = sum(super_tot[k] for k in ("EUR", "EAS", "SAS"))
called = (na_tot > 0).sum(0)
hom_alt = ((na_alt == na_tot) | (na_tot == 0)).all(0) & (called >= 2)
hom_ref = ((na_alt == 0) | (na_tot == 0)).all(0) & (called >= 2)
nonafr_f = out_alt / np.maximum(out_tot, 1)
arch_alt = hom_alt & (afr_f < 0.01) & (super_tot["AFR_all"] >= 100) & (nonafr_f >= 0.01)
arch_ref = hom_ref & ((1 - afr_f) < 0.01) & (super_tot["AFR_all"] >= 100) & ((1 - nonafr_f) >= 0.01)
sel = np.flatnonzero(arch_alt | arch_ref)
allele = np.where(arch_alt[sel], 1, 0).astype(np.uint8)
freqs_out = {}
for k in ("EUR", "EAS", "SAS", "AMR", "AFR"):
    f = super_alt[k][sel] / np.maximum(super_tot[k][sel], 1)
    freqs_out[f"freq_{k}"] = np.where(allele == 1, f, 1 - f).astype(np.float32)
np.savez_compressed(OUT / "archaic.npz", snp_idx=sel.astype(np.int32), allele=allele, **freqs_out)
log("archaic alleles", len(sel))

# per-group summaries of the archaic statistics in present-day people
def summ(num_key, den_key, ratio=True):
    out = {}
    for j, i in enumerate(present_all):
        n, d = indstat[num_key][j].sum(), indstat[den_key][j].sum()
        if d == 0:
            continue
        out.setdefault(group[i], []).append(n / d)
    res = []
    for gname, v in sorted(out.items()):
        v = np.array(v)
        res.append({"group": gname, "region": REGION[gname], "n": int(len(v)), "mean": float(v.mean()),
                    "sd": float(v.std(ddof=1)) if len(v) > 1 else None})
    by_region = {}
    for j, i in enumerate(present_all):
        n, d = indstat[num_key][j].sum(), indstat[den_key][j].sum()
        if d:
            by_region.setdefault(REGION[group[i]], []).append(n / d)
    reg = [{"region": r, "n": len(v), "mean": float(np.mean(v)), "sd": float(np.std(v, ddof=1)) if len(v) > 1 else None}
           for r, v in sorted(by_region.items())]
    return {"groups": res, "regions": reg}


refstats = {"neanderthal": summ("nean_num", "nean_den"), "denisovan": summ("den_num", "den_den")}
(OUT / "reference_stats.json").write_text(json.dumps(refstats, indent=1))

# ------------------------------------------------------------------ validation bundle (not shipped)
if VALDIR:
    VALDIR.mkdir(parents=True, exist_ok=True)
    gh = np.zeros((len(holdout), S), np.int8)
    for c0 in range(0, S, 100000):
        c1 = min(S, c0 + 100000)
        s_lo, s_hi = int(S_idx[c0]), int(S_idx[c1 - 1]) + 1
        G = decode(np.array(holdout), s_lo, s_hi)[:, S_idx[c0:c1] - s_lo]
        gh[:, c0:c1] = np.where(G == 3, -1, 2 - G.astype(np.int8))
    yidx = np.flatnonzero(ysnp)
    have_y = bool(ymales) and len(yidx) > 0
    gy = (decode(np.array(ymales), int(yidx.min()), int(yidx.max()) + 1)[:, yidx - int(yidx.min())]
          if have_y else np.zeros((0, 0), np.uint8))
    np.savez_compressed(VALDIR / "holdout.npz", geno=gh, ids=np.array([ids[i] for i in holdout]),
                        groups=np.array([group[i] for i in holdout]),
                        pc_world=np.array([proj["world"][i] for i in holdout]),
                        pc_west=np.array([proj["west"][i] for i in holdout]),
                        y_geno=np.where(gy == 3, -1, 2 - gy.astype(np.int8)) if have_y else np.zeros((0, 0), np.int8),
                        y_pos=pos[yidx], y_a1=a1[yidx], y_a2=a2[yidx],
                        y_ids=np.array([ids[i] for i in ymales] if have_y else [], dtype=object).astype(str),
                        y_calls=np.array([a_get(ids[i], "y") or "" for i in ymales] if have_y else [], dtype=object).astype(str))
    # population-level frequencies of the validation groups for ADMIXTOOLS cross-checks
    json.dump({"holdout": [ids[i] for i in holdout], "roles": {k: [ids[i] for i in v] for k, v in ROLES.items()},
               "world_ref": [ids[i] for i in world_ref], "west_ref": [ids[i] for i in west_ref]},
              open(VALDIR / "sets.json", "w"))
    t = PREFIX.with_suffix(".truth.json")
    if t.exists():
        shutil.copy(t, VALDIR / "truth.json")
    log("validation bundle written")

manifest = {"name": "pgdata-1", "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "sources": json.loads(Path(OUT / "sources.json").read_text()) if (OUT / "sources.json").exists() else {},
            "roles": {k: {"groups": sorted({group[i] for i in v}), "individuals": len(v)} for k, v in ROLES.items()},
            "models": MODELS, "archaic_roles": ARCH,
            "archaic_definition": "All of the Altai, Vindija and Chagyrskaya Neanderthals that are called carry the allele (at least two called); under 1% in sub-Saharan African reference genomes (at least 100 chromosomes); at least 1% in European, East Asian and South Asian 1000 Genomes samples combined.",
            "pca": [{"id": "world", "title": "Worldwide", "reference_individuals": len(world_ref)},
                    {"id": "west", "title": "West Eurasia", "reference_individuals": len(west_ref)}],
            "snps": int(S), "ancients": len(anc_out), "mt_matched": mt_matched, "y_matched": y_matched}
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=1))
log("done")
