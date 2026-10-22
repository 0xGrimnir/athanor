"""A small synthetic stand-in for the AADR files (same formats and labels),
for exercising build.py and validate.py without the 7 GB download.

Usage: python synthetic.py OUT_PREFIX
Writes OUT_PREFIX.{snp,ind,geno (TGENO),anno}.
"""
import sys

import msprime
import numpy as np

out = sys.argv[1]
rng = np.random.default_rng(7)
# (label, present-day?, individuals, branch, sampling time in generations before present)
G0 = 29                      # years per generation, for reading the dates below
POPS = [("Chimp", False, 1, "chimp", 0), ("Altai_Neanderthal", False, 1, "nea", 3600),
        ("VindijaCave_Neanderthal", False, 1, "nea", 3000), ("ChagyrskayaCave_Neanderthal", False, 1, "nea", 3400),
        ("Russia_Denisova", False, 2, "den", 3600),
        ("Mbuti", True, 6, "afr", 0), ("YRI", True, 30, "afr", 0), ("LWK", True, 6, "afr", 0),
        ("Yoruba", True, 6, "afr", 0),
        ("Luxembourg_Loschbour_Mesolithic", False, 16, "whg", 280),
        ("Turkey_N", False, 20, "anf_pop", 300),
        ("Russia_Samara_EBA_Yamnaya", False, 20, "ste", 170),
        ("Iran_GanjDareh_N", False, 6, "irn", 330),
        ("Georgia_KotiasKlde_Mesolithic", False, 4, "irn", 350),
        ("Russia_Chelyabinsk_MLBA_Sintashta", False, 6, "ste", 135),
        ("ONG", True, 4, "asi", 0), ("Russia_UstIshim_IUP", False, 1, "basal", 2500),
        ("Russia_Kostenki_UP", False, 1, "basal", 2200), ("Russia_Malta_UP", False, 3, "ane", 1600),
        ("Han", True, 8, "eas", 0), ("CHB", True, 12, "eas", 0), ("JPT", True, 12, "eas", 0),
        ("Papuan", True, 6, "pap", 0), ("Chukchi", True, 1, "eas", 0), ("Karitiana", True, 6, "amr", 0),
        ("Israel_Natufian", False, 4, "nat", 420), ("Ethiopia_MotaCave_4500BP", False, 2, "afr", 155),
        ("GBR", True, 36, "eur", 0), ("CEU", True, 36, "eur", 0), ("FIN", True, 28, "eur", 0),
        ("IBS", True, 28, "eur", 0), ("TSI", True, 28, "eur", 0), ("French", True, 8, "eur", 0),
        ("Sardinian", True, 10, "anf_pop", 0), ("Druze", True, 10, "irn", 0),
        ("GIH", True, 12, "sas", 0), ("PJL", True, 12, "sas", 0), ("PEL", True, 10, "amr", 0),
        ("MXL", True, 10, "amr", 0), ("PUR", True, 10, "eur", 0),
        ("Germany_BellBeaker", False, 6, "eur", 100), ("hg19ref", True, 1, "chimp", 0)]
BR = ["chimp", "nea", "den", "afr", "basal", "ane", "whg", "ehg", "anf_pop", "nat", "irn", "ste", "eur",
      "eas", "pap", "amr", "sas", "asi"]
INTERNAL = ["root", "anc_all", "anc_arch", "anc_hum", "anc_ooa", "ooa", "anc_main", "anc_west", "anc_east",
            "anc_wf", "anc_sw", "anc_irn", "anf_base", "eas_anc", "nat_base"]
d = msprime.Demography()
for b_ in BR + INTERNAL:
    d.add_population(name=b_, initial_size=2000)
# Every sampled population is a leaf, and every internal node is its own
# population: msprime deactivates a population once it has split, so a branch
# that is sampled can never also be a split's ancestor.
d.add_admixture(time=120, derived="eur", ancestral=["whg", "anf_pop", "ste"], proportions=[0.15, 0.45, 0.40])
d.add_admixture(time=160, derived="sas", ancestral=["irn", "ste", "asi"], proportions=[0.5, 0.2, 0.3])
d.add_admixture(time=250, derived="ste", ancestral=["ehg", "irn"], proportions=[0.5, 0.5])
d.add_admixture(time=350, derived="ehg", ancestral=["whg", "ane"], proportions=[0.6, 0.4])
d.add_admixture(time=500, derived="anf_pop", ancestral=["anf_base", "nat"], proportions=[0.75, 0.25])
d.add_population_split(time=500, derived=["eas", "amr"], ancestral="eas_anc")
# Natufians (and through them the Anatolian farmers) carry Basal Eurasian
# ancestry, which is what lets the African and early-Eurasian outgroups tell
# the farmer source apart from the hunter-gatherer source.
d.add_admixture(time=1300, derived="nat", ancestral=["nat_base", "basal"], proportions=[0.6, 0.4])
d.add_population_split(time=1400, derived=["irn", "nat_base"], ancestral="anc_irn")
d.add_population_split(time=1600, derived=["anc_irn", "anf_base"], ancestral="anc_sw")
d.add_population_split(time=1800, derived=["eas_anc", "pap", "asi"], ancestral="anc_east")
d.add_population_split(time=1900, derived=["whg", "anc_sw"], ancestral="anc_wf")
d.add_population_split(time=2200, derived=["anc_wf", "ane"], ancestral="anc_west")
d.add_population_split(time=2600, derived=["anc_west", "anc_east"], ancestral="anc_main")
d.add_population_split(time=2800, derived=["anc_main", "basal"], ancestral="ooa")
d.add_admixture(time=3400, derived="ooa", ancestral=["anc_ooa", "nea"], proportions=[0.98, 0.02])
d.add_population_split(time=4000, derived=["anc_ooa", "afr"], ancestral="anc_hum")
d.add_population_split(time=15000, derived=["nea", "den"], ancestral="anc_arch")
d.add_population_split(time=25000, derived=["anc_hum", "anc_arch"], ancestral="anc_all")
d.add_population_split(time=60000, derived=["anc_all", "chimp"], ancestral="root")
d.sort_events()
groups = {}
for lab, present, n, br, t in POPS:
    groups.setdefault((br, t), 0)
    groups[(br, t)] += n
samples = [msprime.SampleSet(n, population=br, time=t, ploidy=2) for (br, t), n in groups.items()]
G, POS, CH = [], [], []
for c in range(1, 23):
    ts = msprime.sim_ancestry(samples=samples, demography=d, sequence_length=7_000_000, recombination_rate=1e-8,
                              random_seed=c)
    ts = msprime.sim_mutations(ts, rate=2e-8, random_seed=50 + c, model=msprime.BinaryMutationModel())
    gm = ts.genotype_matrix()
    p = ts.tables.sites.position.astype(np.int64)
    keep = rng.random(len(p)) < 0.3
    gm = gm[keep]
    # collapse to per-individual allele counts straight away to keep memory down
    G.append((gm[:, 0::2] + gm[:, 1::2]).astype(np.int8))
    POS.append(p[keep] + 1); CH.append(np.full(keep.sum(), c))
    del gm, ts
G = np.vstack(G); POS = np.concatenate(POS); CH = np.concatenate(CH)
# haplotype columns grouped by branch in sample order
cols, k = {}, 0
for key, n in groups.items():
    cols[key] = list(range(k, k + 2 * n)); k += 2 * n
used = {key: 0 for key in groups}
inds = []
for lab, present, n, br, t in POPS:
    key = (br, t)
    for j in range(n):
        h = cols[key][used[key]:used[key] + 2]; used[key] += 2
        inds.append((f"{lab}_{j}", lab, present, h, t))
nsnp, nind = G.shape[0], len(inds)
geno = np.zeros((nind, nsnp), np.uint8)
for i, (_, lab, present, h, t) in enumerate(inds):
    alt = G[:, h[0] // 2].astype(np.int64)
    if lab == "hg19ref":
        alt = np.zeros(nsnp, int)
    elif not present and lab not in ("Altai_Neanderthal", "VindijaCave_Neanderthal", "ChagyrskayaCave_Neanderthal", "Chimp"):
        # ancient: pseudo-haploid with missing data
        # pseudo-haploid: one allele drawn at random at heterozygous sites,
        # as ancient genomes are represented in the AADR
        cnt = G[:, h[0] // 2]
        hap = np.where(cnt == 2, 2, np.where(cnt == 1, rng.integers(0, 2, nsnp) * 2, 0))
        alt = np.where(rng.random(nsnp) < 0.25, -1, hap)
    geno[i] = np.where(alt < 0, 3, 2 - alt)           # copies of allele 1 (reference)
with open(out + ".snp", "w") as f:
    for i in range(nsnp):
        f.write(f"rs{100000 + i} {CH[i]} {POS[i] * 1e-8:.6f} {POS[i]} {'ACGT'[i % 4]} {'CGTA'[i % 4]}\n")
with open(out + ".ind", "w") as f:
    for name, lab, present, _, t in inds:
        f.write(f"{name} U {lab}\n")
def _hashit(s):
    h = 0
    for ch in s.encode():
        h = (h * 23 + ch) & 0xFFFFFFFF
    return h


def _hasharr(names):
    h = 0
    for n in names:
        h = ((h * 17) & 0xFFFFFFFF) ^ _hashit(n)
    return h


rlen = (nsnp + 3) // 4
ihash = _hasharr([name for name, *_ in inds])
shash = _hasharr([f"rs{100000 + i}" for i in range(nsnp)])
with open(out + ".geno", "wb") as f:
    f.write(f"TGENO {nind} {nsnp} {ihash:x} {shash:x}".encode().ljust(48, b"\x00"))
    pad = np.full((nind, rlen * 4), 3, np.uint8); pad[:, :nsnp] = geno
    q = pad.reshape(nind, rlen, 4)
    f.write(((q[:, :, 0] << 6) | (q[:, :, 1] << 4) | (q[:, :, 2] << 2) | q[:, :, 3]).astype(np.uint8).tobytes())
hdr = ["Genetic ID", "Persistent Genetic ID", "Individual ID", "Skeletal code", "Skeletal element",
       "First publication", "Publication abbreviation", "doi for publication of this representation of the data",
       "Link", "Method for Determining Date", "Date mean in BP in years before 1950 CE", "Date standard deviation in BP",
       "Full Date One of two formats", "Age at death", "Group ID", "Locality", "Political Entity", "Latitude", "Longitude",
       "Pulldown", "Suffices", "Data type", "No. Libraries", "Mean coverage", "Mean coverage non-targeted",
       "SNPs hit on autosomal targets (Computed using easystats on enhance 2M capture subset)",
       "SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)", "SNPs HO", "SNPs compat", "SNPs compat HO",
       "Molecular Sex", "Family relations", "ROH1", "ROH2",
       "Y haplogroup in terminal mutation notation automatically called based on Y-full 12.03",
       "Y haplogroup  in ISOGG notation automatically called", "Y haplogroup manually called if different from automatic",
       "mtDNA coverage", "mtDNA haplogroup if >2x or published", "mtDNA match", "Damage", "Sex ratio", "ANGSD", "hapConX",
       "Library type", "Libraries", "endogenous", "ASSESSMENT", "ASSESSMENT WARNINGS"]
MT = ["H1c3", "U5b1a", "K1a", "T2b", "J1c", "H", "N9a", "B4", "L3e"]
YS = ["R-M269", "I-M253", "G-M201", "J-M172", "R-M417", "I-V6473", "N-Y6503"]
with open(out + ".anno", "w") as f:
    f.write("\t".join(hdr) + "\n")
    for i, (name, lab, present, _, t) in enumerate(inds):
        row = [".."] * len(hdr)
        row[0], row[1], row[2] = name, str(i), name
        row[6], row[7] = "SyntheticTest", "https://doi.org/10.0000/synthetic"
        row[10] = "0" if present else str(int(t * G0))
        row[11] = "0" if present else "100"
        row[14], row[15], row[16] = lab, "Somewhere", "Nowhere"
        row[17], row[18] = f"{rng.uniform(30, 60):.2f}", f"{rng.uniform(-10, 60):.2f}"
        row[26] = str(int(nsnp * (1.0 if present else 0.6)))
        row[30] = "M" if i % 2 == 0 else "F"
        row[34] = YS[i % len(YS)] if i % 2 == 0 else ".."
        row[38] = MT[i % len(MT)]
        row[47] = "PASS"
        f.write("\t".join(row) + "\n")
# the proportions the European populations were simulated with, for validate.py
truth = {lab: {"WHG": 0.15, "Anatolia_N": 0.45, "Yamnaya": 0.40}
         for lab, present, n, br, t in POPS if br == "eur" and present}
import json as _json
open(out + ".truth.json", "w").write(_json.dumps(truth, indent=1))
print("SNPs", nsnp, "individuals", nind, "truth groups", len(truth))
