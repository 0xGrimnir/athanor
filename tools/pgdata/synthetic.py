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
# (label, present-day?, individuals, branch)
POPS = [("Chimp", False, 1, "chimp"), ("Altai_Neanderthal", False, 1, "nea"), ("VindijaCave_Neanderthal", False, 1, "nea"),
        ("ChagyrskayaCave_Neanderthal", False, 1, "nea"), ("Russia_Denisova", False, 2, "den"),
        ("Mbuti", True, 6, "afr"), ("YRI", True, 30, "afr"), ("LWK", True, 6, "afr"), ("Yoruba", True, 6, "afr"),
        ("Luxembourg_Loschbour_Mesolithic", False, 2, "whg"), ("Turkey_N", False, 8, "eef"),
        ("Russia_Samara_EBA_Yamnaya", False, 8, "ste"), ("Iran_GanjDareh_N", False, 6, "irn"),
        ("Georgia_KotiasKlde_Mesolithic", False, 2, "irn"), ("Russia_Chelyabinsk_MLBA_Sintashta", False, 6, "ste"),
        ("ONG", True, 4, "asi"), ("Russia_UstIshim_IUP", False, 1, "ooa"), ("Russia_Kostenki_UP", False, 1, "ooa"),
        ("Russia_Malta_UP", False, 1, "ane"), ("Han", True, 8, "eas"), ("CHB", True, 12, "eas"), ("JPT", True, 12, "eas"),
        ("Papuan", True, 6, "pap"), ("Chukchi", True, 1, "eas"), ("Karitiana", True, 6, "amr"),
        ("Israel_Natufian", False, 4, "eef"), ("Ethiopia_MotaCave_4500BP", False, 2, "afr"),
        ("GBR", True, 20, "eur"), ("CEU", True, 20, "eur"), ("FIN", True, 14, "eur"), ("IBS", True, 14, "eur"),
        ("TSI", True, 14, "eur"), ("French", True, 8, "eur"), ("Sardinian", True, 10, "eef"), ("Druze", True, 10, "irn"),
        ("GIH", True, 12, "sas"), ("PJL", True, 12, "sas"), ("PEL", True, 10, "amr"), ("MXL", True, 10, "amr"),
        ("PUR", True, 10, "eur"), ("Germany_BellBeaker", False, 6, "eur"), ("hg19ref", True, 1, "chimp")]
BR = ["chimp", "nea", "den", "afr", "ooa", "ane", "whg", "eef", "irn", "ste", "eur", "eas", "pap", "amr", "sas", "asi"]
d = msprime.Demography()
for b in BR + ["anc_hum", "anc_arch", "anc_all", "anc_ooa", "anc_west", "anc_east", "anc_wf", "root"]:
    d.add_population(name=b, initial_size=5000)
d.add_admixture(time=150, derived="eur", ancestral=["whg", "eef", "ste"], proportions=[0.15, 0.45, 0.40])
d.add_admixture(time=160, derived="sas", ancestral=["irn", "ste", "asi"], proportions=[0.5, 0.2, 0.3])
d.add_admixture(time=1700, derived="ooa", ancestral=["anc_ooa", "nea"], proportions=[0.98, 0.02])
d.add_population_split(time=200, derived=["ste"], ancestral="irn")
d.add_population_split(time=400, derived=["eef"], ancestral="irn")
d.add_population_split(time=600, derived=["irn", "whg"], ancestral="anc_wf")
d.add_population_split(time=700, derived=["anc_wf", "ane"], ancestral="anc_west")
d.add_population_split(time=500, derived=["amr"], ancestral="eas")
d.add_population_split(time=900, derived=["eas", "pap", "asi"], ancestral="anc_east")
d.add_population_split(time=1500, derived=["anc_west", "anc_east"], ancestral="ooa")
d.add_population_split(time=2000, derived=["anc_ooa", "afr"], ancestral="anc_hum")
d.add_population_split(time=15000, derived=["nea", "den"], ancestral="anc_arch")
d.add_population_split(time=25000, derived=["anc_hum", "anc_arch"], ancestral="anc_all")
d.add_population_split(time=60000, derived=["anc_all", "chimp"], ancestral="root")
d.sort_events()
counts = {}
for lab, _, n, br in POPS:
    counts[br] = counts.get(br, 0) + n
samples = [msprime.SampleSet(n, population=br, ploidy=2) for br, n in counts.items()]
G, POS, CH = [], [], []
for c in range(1, 23):
    ts = msprime.sim_ancestry(samples=samples, demography=d, sequence_length=3_000_000, recombination_rate=1e-8,
                              random_seed=c)
    ts = msprime.sim_mutations(ts, rate=2e-8, random_seed=50 + c, model=msprime.BinaryMutationModel())
    gm = ts.genotype_matrix()
    p = ts.tables.sites.position.astype(np.int64)
    keep = rng.random(len(p)) < 0.3
    G.append(gm[keep]); POS.append(p[keep] + 1); CH.append(np.full(keep.sum(), c))
G = np.vstack(G); POS = np.concatenate(POS); CH = np.concatenate(CH)
# haplotype columns grouped by branch in sample order
offset, cols = {}, {}
k = 0
for br, n in counts.items():
    cols[br] = list(range(k, k + 2 * n)); k += 2 * n
used = {br: 0 for br in counts}
inds = []
for lab, present, n, br in POPS:
    for j in range(n):
        h = cols[br][used[br]:used[br] + 2]; used[br] += 2
        inds.append((f"{lab}_{j}", lab, present, h))
nsnp, nind = G.shape[0], len(inds)
geno = np.zeros((nind, nsnp), np.uint8)
for i, (_, lab, present, h) in enumerate(inds):
    alt = G[:, h[0]] + G[:, h[1]]
    if lab == "hg19ref":
        alt = np.zeros(nsnp, int)
    elif not present and lab not in ("Altai_Neanderthal", "VindijaCave_Neanderthal", "ChagyrskayaCave_Neanderthal", "Chimp"):
        # ancient: pseudo-haploid with missing data
        hap = G[:, h[0]] * 2
        alt = np.where(rng.random(nsnp) < 0.4, -1, hap)
    geno[i] = np.where(alt < 0, 3, 2 - alt)           # copies of allele 1 (reference)
with open(out + ".snp", "w") as f:
    for i in range(nsnp):
        f.write(f"rs{100000 + i} {CH[i]} {POS[i] * 1e-8:.6f} {POS[i]} {'ACGT'[i % 4]} {'CGTA'[i % 4]}\n")
with open(out + ".ind", "w") as f:
    for name, lab, present, _ in inds:
        f.write(f"{name} U {lab}\n")
rlen = (nsnp + 3) // 4
with open(out + ".geno", "wb") as f:
    f.write(f"TGENO {nind} {nsnp} 0 0".encode().ljust(48, b"\x00"))
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
    for i, (name, lab, present, _) in enumerate(inds):
        row = [".."] * len(hdr)
        row[0], row[1], row[2] = name, str(i), name
        row[6], row[7] = "SyntheticTest", "https://doi.org/10.0000/synthetic"
        row[10] = "0" if present else str(int(rng.integers(3000, 45000)))
        row[11] = "0" if present else "100"
        row[14], row[15], row[16] = lab, "Somewhere", "Nowhere"
        row[17], row[18] = f"{rng.uniform(30, 60):.2f}", f"{rng.uniform(-10, 60):.2f}"
        row[26] = str(int(nsnp * (1.0 if present else 0.6)))
        row[30] = "M" if i % 2 == 0 else "F"
        row[34] = YS[i % len(YS)] if i % 2 == 0 else ".."
        row[38] = MT[i % len(MT)]
        row[47] = "PASS"
        f.write("\t".join(row) + "\n")
print("SNPs", nsnp, "individuals", nind)
