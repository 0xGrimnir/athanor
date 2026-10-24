"""Data and checks for DNA Structure's zoom through the scales.

1. Nucleosome: the 1.9 A nucleosome core particle (PDB 1KX5; Davey et al., J. Mol. Biol. 319, 1097, 2002),
   147 bp on a histone octamer. Its DNA is analysed into base-pair and step parameters, so any sequence
   can be threaded onto the same path. Checks: rebuilt from its own parameters, the crystal's bases come
   back; the fitted superhelix has the radius and pitch Richmond & Davey reported (Nature 423, 145,
   2003: 41.9 A and 25.9 A) and about 1.65 turns (Luger et al., Nature 389, 251, 1997).
2. Protein-DNA complexes: DNA bend angles measured on crystal structures, checked against the bends the
   original papers report (TBP 80 deg, CAP 90 deg, IHF more than 160 deg).
3. Chromatin fibre: a six-nucleosome array with linker histone H1 (PDB 6HKT; Garcia-Saez et al., Mol.
   Cell 72, 902, 2018). Check: it is a two-start zigzag - each nucleosome is closer to the one two
   along than to its neighbour on the DNA.
4. Nucleus: single-cell 3D genome structures of GM12878 cells from Dip-C (Tan et al., Science 361, 924,
   2018; GEO GSE117876), coarse-grained to 1 Mb. Check: gene-rich chromosomes sit nearer the nuclear
   centre (Croft et al., J. Cell Biol. 145, 1119, 1999: chromosome 19 interior, 18 peripheral).

Usage: python scales.py OUT_DIR WORK_DIR
"""
from __future__ import annotations

import gzip
import json
import sys
import time
import urllib.request
from collections import OrderedDict
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.dna.builder import build  # noqa: E402
from athanor_core.dna.structure import Residue, analyse, read_pdb  # noqa: E402

OUT, WORK = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
results, summary = [], {}
NUC = {"DA", "DC", "DG", "DT", "DU"}


def gate(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def get(url, tries=4):
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "athanor-ci", "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(3 * (k + 1))
    raise RuntimeError(f"{url}: {err!r}")


def pdb_text(pid):
    f = WORK / f"{pid}.pdb"
    if not f.exists():
        f.write_bytes(get(f"https://files.rcsb.org/download/{pid}.pdb"))
    return f.read_text()


def protein_ca(text, chains=None):
    """C-alpha trace per protein chain (first model, first altloc): {chain: [[x,y,z], ...]}."""
    out = OrderedDict()
    seen = set()
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[12:16].strip() != "CA" or line[17:20].strip() in NUC:
            continue
        if line[16] not in (" ", "A"):
            continue
        c = line[21]
        if chains and c not in chains:
            continue
        key = (c, line[22:27])
        if key in seen:
            continue
        seen.add(key)
        out.setdefault(c, []).append([round(float(line[30:38]), 2), round(float(line[38:46]), 2), round(float(line[46:54]), 2)])
    return out


def strands(residues, std=None):
    """Two strands from any number of DNA chains: chains that pair with each other go on opposite
    strands; chains on the same strand are joined in the order that keeps the duplex continuous."""
    from athanor_core.dna.structure import base_frame, chains as by_chain, standard_bases
    std = std or standard_bases()
    by = by_chain(residues)
    names = list(by)
    fr = {c: [base_frame(r, std) for r in by[c]] for c in names}
    link = {}
    for a in names:
        for b in names:
            if a >= b:
                continue
            n = 0
            for x in fr[a]:
                if x is None:
                    continue
                for y in fr[b]:
                    if y is not None and np.linalg.norm(x[1] - y[1]) < 4.0 and np.dot(x[0][:, 2], y[0][:, 2]) < -0.5:
                        n += 1
                        break
            if n:
                link[(a, b)] = link[(b, a)] = n
    side = {names[0]: 0}
    changed = True
    while changed:
        changed = False
        for (a, b) in link:
            if a in side and b not in side:
                side[b] = 1 - side[a]
                changed = True
    paired = [c for c in names if c in side]
    s1 = [c for c in paired if side[c] == 0]
    s2 = [c for c in paired if side[c] == 1]

    def order(group, other):
        # position of each chain along the partner strand: the mean index of the partner bases it pairs with
        partner = [r for c in other for r in by[c]]
        pf = [base_frame(r, std) for r in partner]
        key = {}
        for c in group:
            idx = []
            for x in fr[c]:
                if x is None:
                    continue
                for j, y in enumerate(pf):
                    if y is not None and np.linalg.norm(x[1] - y[1]) < 4.0:
                        idx.append(j)
            key[c] = -np.mean(idx) if idx else 0.0          # antiparallel: later partner index comes first
        return sorted(group, key=lambda c: key[c])
    s2 = order(s2, s1) if len(s2) > 1 else s2
    s1 = order(s1, s2) if len(s1) > 1 else s1
    out = []
    for name, group in (("A", s1), ("B", s2)):
        k = 0
        for c in group:
            for r in by[c]:
                k += 1
                out.append(Residue(name, k, "", r.name, r.atoms))
    return out, s1, s2


def helix_bend(an, m=4):
    """Angle (deg) between the mean base-pair z axes of the first and last m base pairs."""
    R = [np.array(p["frame"][0]) for p in an["pairs"]]
    a = np.mean([r[:, 2] for r in R[:m]], axis=0)
    b = np.mean([r[:, 2] for r in R[-m:]], axis=0)
    c = np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def superhelix(points):
    """Fit a helix to points: axis, radius, pitch and number of turns."""
    P = np.asarray(points)
    c0 = P.mean(0)
    w, v = np.linalg.eigh(np.cov((P - c0).T))
    axis0 = v[:, 0]                       # flat superhelix: least variance along its axis

    def unpack(x):
        ax = x[:3] / np.linalg.norm(x[:3])
        return ax, x[3:6]

    def resid(x):
        ax, c = unpack(x)
        d = P - c
        h = d @ ax
        rvec = d - np.outer(h, ax)
        r = np.linalg.norm(rvec, axis=1)
        # angle about the axis, unwrapped along the path
        e1 = rvec[0] / np.linalg.norm(rvec[0])
        e2 = np.cross(ax, e1)
        th = np.unwrap(np.arctan2(rvec @ e2, rvec @ e1))
        A = np.vstack([th, np.ones_like(th)]).T
        coef, *_ = np.linalg.lstsq(A, h, rcond=None)
        return np.concatenate([r - r.mean(), (h - A @ coef)])
    sol = least_squares(resid, np.concatenate([axis0, c0]))
    ax, c = unpack(sol.x)
    d = P - c
    h = d @ ax
    rvec = d - np.outer(h, ax)
    r = np.linalg.norm(rvec, axis=1)
    e1 = rvec[0] / np.linalg.norm(rvec[0])
    e2 = np.cross(ax, e1)
    th = np.unwrap(np.arctan2(rvec @ e2, rvec @ e1))
    slope = np.polyfit(th, h, 1)[0]
    turns = abs(th[-1] - th[0]) / (2 * np.pi)
    # handedness: a left-handed superhelix advances along the axis while turning clockwise
    hand = "left" if np.sign(slope) * np.sign(th[-1] - th[0]) < 0 else "right"
    return {"radius": float(r.mean()), "pitch": float(abs(slope) * 2 * np.pi), "turns": float(turns), "hand": hand,
            "axis": ax.tolist(), "centre": c.tolist()}


# ------------------------------------------------------------------ 1. nucleosome
txt = pdb_text("1KX5")
res = [r for r in read_pdb(txt) if r.chain in ("I", "J")]
an = analyse(res)
n = len(an["pairs"])
seq = an["sequence"]
print(f"1KX5: {n} base pairs, {len(an['steps'])} steps, sequence {seq[:20]}...")
steps = [s["params"] for s in an["steps"]]
pairs = [p["params"] for p in an["pairs"]]
R0, o0 = (np.array(x) for x in an["pairs"][0]["frame"])
m = build(seq, "B", steps=steps, pairs=pairs)
# place the built DNA in the crystal's frame: the builder starts from the identity frame
built = {(a[0], a[1], a[3]): R0 @ np.asarray(a[4]) + o0 for a in m["atoms"]}
s1 = [r for r in res if r.chain == "I"]
s2 = [r for r in res if r.chain == "J"]
P, Q, PB, QB = [], [], [], []
BACK = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"}
for k, pinfo in enumerate(an["pairs"]):
    for chain, num, cr in (("A", k + 1, s1[pinfo["i"]]), ("B", 2 * n - k, s2[pinfo["j"]])):
        for nm, xyz in cr.atoms.items():
            key = (chain, num, nm)
            if key in built:
                P.append(built[key]); Q.append(xyz)
                if nm not in BACK:
                    PB.append(built[key]); QB.append(xyz)
rb = float(np.sqrt(np.mean(np.sum((np.array(PB) - np.array(QB)) ** 2, axis=1))))
ra = float(np.sqrt(np.mean(np.sum((np.array(P) - np.array(Q)) ** 2, axis=1))))
gate("the 1KX5 nucleosome's DNA, rebuilt from its own parameters in place, puts the crystal's bases back "
     "(base-atom RMSD <= 0.1 A, no superposition)", rb <= 0.1, f"bases {rb:.3f} A, all atoms {ra:.2f} A over {len(P)} atoms")
origins = [np.array(p["frame"][1]) for p in an["pairs"]]
sh = superhelix(origins)
print("superhelix:", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in sh.items()})
gate("the 1KX5 DNA follows a left-handed superhelix of the published radius, pitch and length (Richmond & Davey "
     "2003: 41.9 A and 25.9 A; about 1.65 turns, Luger et al. 1997): within 1.5 A, 3 A and 1.55-1.80 turns",
     sh["hand"] == "left" and abs(sh["radius"] - 41.9) <= 1.5 and abs(sh["pitch"] - 25.9) <= 3.0 and 1.55 <= sh["turns"] <= 1.80,
     f"{sh['hand']}-handed, radius {sh['radius']:.1f} A, pitch {sh['pitch']:.1f} A, {sh['turns']:.2f} turns")
HISTONE = {"A": "H3", "E": "H3", "B": "H4", "F": "H4", "C": "H2A", "G": "H2A", "D": "H2B", "H": "H2B"}
nuc = {"pdb": "1KX5", "title": "Nucleosome core particle, 147 bp on a Xenopus histone octamer (1.9 A)",
       "citation": "Davey, Sargent, Luger, Maeder & Richmond, J. Mol. Biol. 319, 1097 (2002)",
       "sequence": seq, "steps": np.round(steps, 3).tolist(), "pairs": np.round(pairs, 3).tolist(),
       "frame0": [np.round(R0, 6).tolist(), np.round(o0, 3).tolist()],
       "histones": [{"chain": c, "type": HISTONE.get(c, "?"), "ca": ca} for c, ca in protein_ca(txt, set(HISTONE)).items()],
       "superhelix": {k: (round(v, 2) if isinstance(v, float) else v) for k, v in sh.items() if k != "axis" and k != "centre"},
       "superhelix_axis": np.round(sh["axis"], 4).tolist(), "superhelix_centre": np.round(sh["centre"], 2).tolist(),
       "dyad_index": n // 2}
(OUT / "nucleosome.json").write_text(json.dumps(nuc, separators=(",", ":")) + "\n")

# ------------------------------------------------------------------ 2. protein-DNA complexes
COMPLEXES = [
    {"pdb": "1CDW", "name": "TATA-binding protein on the TATA box (human TBP, adenovirus major late promoter)",
     "citation": "Nikolov et al., Proc. Natl. Acad. Sci. USA 93, 4862 (1996)", "published_bend": 80, "tolerance": 15,
     "note": "TBP saddles the minor groove and kinks the DNA at both ends of the TATA box"},
    {"pdb": "1CGP", "name": "Catabolite activator protein (CAP) dimer with cAMP",
     "citation": "Schultz, Shields & Steitz, Science 253, 1001 (1991)", "published_bend": 90, "tolerance": 15,
     "note": "two kinks of about 40 deg at the TG/CA steps add to an overall bend of about 90 deg"},
    {"pdb": "1IHF", "name": "Integration host factor (IHF)",
     "citation": "Rice, Yang, Mizuuchi & Nash, Cell 87, 1295 (1996)", "published_bend": 160, "at_least": True,
     "note": "the DNA wraps around the protein, bent by more than 160 deg"},
    {"pdb": "1AAY", "name": "Zif268 zinc-finger domains", "citation": "Elrod-Erickson et al., Structure 4, 1171 (1996)",
     "note": "three zinc fingers read the major groove of nearly straight B-DNA"},
    {"pdb": "1TSR", "name": "p53 tumour suppressor DNA-binding domain", "citation": "Cho, Gorina, Jeffrey & Pavletich, Science 265, 346 (1994)",
     "note": "the most frequently mutated residues in human cancers contact this DNA"},
    {"pdb": "1EFA", "name": "Lac repressor bound to its operator", "citation": "Bell & Lewis, Nat. Struct. Biol. 7, 209 (2000)",
     "note": "the hinge helices insert into the minor groove and open it"},
]
cx_out = []
for c in COMPLEXES:
    t = pdb_text(c["pdb"])
    dna = [r for r in read_pdb(t)]
    merged, g1, g2 = strands(dna)
    a = analyse(merged)
    bend = helix_bend(a, 4)
    c2 = {k: v for k, v in c.items()}
    c2.update({"strands": [g1, g2], "sequence": a["sequence"], "n_pairs": len(a["pairs"]), "bend": round(bend, 1),
               "steps": [{"step": s["step"], "params": np.round(s["params"], 2).tolist()} for s in a["steps"]],
               "dna": [[r.chain, r.number, r.name, nm, *np.round(x, 2).tolist()] for r in merged for nm, x in r.atoms.items()],
               "protein": [{"chain": ch, "ca": ca} for ch, ca in protein_ca(t).items()]})
    # largest roll: where the protein kinks the DNA
    if a["steps"]:
        k = int(np.argmax([abs(s["params"][4]) for s in a["steps"]]))
        c2["largest_roll"] = {"step": k + 1, "dimer": a["steps"][k]["step"], "roll": round(a["steps"][k]["params"][4], 1)}
    cx_out.append(c2)
    print(f"{c['pdb']}: strands {g1} / {g2}, {len(a['pairs'])} bp, bend {bend:.1f} deg (published {c.get('published_bend')})")
checked = [c for c in cx_out if c.get("published_bend")]
ok = all((c["bend"] >= c["published_bend"] - 15) if c.get("at_least") else abs(c["bend"] - c["published_bend"]) <= c["tolerance"]
         for c in checked)
gate("DNA bend angles measured on the crystal structures match the published bends (TBP 80, CAP 90 deg within 15 deg; "
     "IHF more than 160 deg, measured >= 145 deg)", ok, "; ".join(f"{c['pdb']} {c['bend']:.0f} deg" for c in checked))
(OUT / "complexes.json").write_text(json.dumps(cx_out, separators=(",", ":")) + "\n")

# ------------------------------------------------------------------ 3. chromatin fibre
t = pdb_text("6HKT")
dna_res = [r for r in read_pdb(t) if r.chain in ("I", "J")]
ca = protein_ca(t)
# group histone chains into octamers by proximity of their centres
cents = {c: np.mean(v, axis=0) for c, v in ca.items() if len(v) > 50}
octs, used = [], set()
for c in sorted(cents):
    if c in used:
        continue
    free = [d for d in cents if d not in used]
    grp = sorted(free, key=lambda d: np.linalg.norm(cents[d] - cents[c]))[:8]
    if max(np.linalg.norm(cents[d] - cents[c]) for d in grp) > 75.0:
        continue                                       # not a whole octamer (e.g. linker histone)
    used.update(grp)
    octs.append(grp)
centres = [np.mean([cents[c] for c in g], axis=0) for g in octs]
# order nucleosomes along the DNA: by the strand-I index of the nearest phosphate
s1 = [r for r in dna_res if r.chain == "I"]
pP = np.array([r.atoms.get("P", r.atoms.get("C1'")) for r in s1 if ("P" in r.atoms or "C1'" in r.atoms)])
order_idx = [int(np.argmin(np.linalg.norm(pP - cc, axis=1))) for cc in centres]
order = np.argsort(order_idx)
centres = [centres[i] for i in order]
octs = [octs[i] for i in order]
d1 = [np.linalg.norm(centres[i] - centres[i + 1]) for i in range(len(centres) - 1)]
d2 = [np.linalg.norm(centres[i] - centres[i + 2]) for i in range(len(centres) - 2)]
print(f"6HKT: {len(octs)} octamers of {[len(g) for g in octs]} chains; i,i+1 distances {np.round(d1, 0).tolist()}; i,i+2 {np.round(d2, 0).tolist()}")
gate("the 6HKT array is a two-start zigzag: six nucleosomes, each nearer the nucleosome two along the DNA than its "
     "neighbour (mean centre distances)", len(octs) == 6 and np.mean(d2) < np.mean(d1),
     f"{len(octs)} nucleosomes; neighbours {np.mean(d1):.0f} A apart, next-but-one {np.mean(d2):.0f} A")
pt = {"I": [], "J": []}
for r in dna_res:
    x = r.atoms.get("P", r.atoms.get("C1'"))
    if x is not None:
        pt[r.chain].append(np.round(x, 1).tolist())
fibre = {"pdb": "6HKT", "title": "Six-nucleosome array with linker histone H1 (9.7 A)",
         "citation": "Garcia-Saez, Menoni, Boopathi, Shukla, Soueidan, Marcoux-Roussel, Nikolov, Elsasser, Bhatt, Dimitrov "
                     "& Petosa, Mol. Cell 72, 902 (2018)",
         "dna": pt, "nucleosomes": [{"chains": g, "centre": np.round(cc, 1).tolist(),
                                    "ca": {c: ca[c] for c in g}} for g, cc in zip(octs, centres)],
         "other": {c: v for c, v in ca.items() if c not in used},
         "distances": {"neighbour": round(float(np.mean(d1)), 1), "next_but_one": round(float(np.mean(d2)), 1)}}
(OUT / "fibre.json").write_text(json.dumps(fibre, separators=(",", ":")) + "\n")

# ------------------------------------------------------------------ 4. nucleus (Dip-C)
BIN = 1_000_000
filelist = get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117876/suppl/filelist.txt").decode().splitlines()
CELLS = []
for line in filelist:
    f = line.split("\t")[1] if "\t" in line else ""
    if f.endswith(".impute3.round4.clean.3dg.txt.gz") and "_gm12878_" in f and "_rep" not in f:
        gsm, rest = f.split("_", 1)
        CELLS.append((gsm, rest.replace(".impute3.round4.clean.3dg.txt.gz", "")))
CELLS = CELLS[:6]
print("Dip-C cells:", CELLS)


def load_3dg(gsm, name):
    f = WORK / f"{gsm}_{name}.3dg.txt.gz"
    if not f.exists():
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{gsm}_{name}.impute3.round4.clean.3dg.txt.gz"
        f.write_bytes(get(url))
    acc = {}
    for line in gzip.decompress(f.read_bytes()).decode().splitlines():
        p = line.split()
        if len(p) < 5:
            continue
        hom, pos = p[0], int(p[1])
        acc.setdefault((hom, pos // BIN), []).append([float(x) for x in p[2:5]])
    out = {}
    for (hom, b), xs in acc.items():
        out.setdefault(hom, []).append((b, np.mean(xs, axis=0)))
    return {h: sorted(v, key=lambda t: t[0]) for h, v in out.items()}


cells = []
for gsm, name in CELLS:
    try:
        cells.append((gsm, name, load_3dg(gsm, name)))
        print(f"{gsm} {name}: {len(cells[-1][2])} homologs")
    except Exception as e:  # noqa: BLE001
        print(f"{gsm} {name}: not loaded ({e!r})"[:200])
# gene density per chromosome: protein-coding genes per Mb (Ensembl GRCh37), counted in 5 Mb windows
assembly = json.loads(get("https://grch37.rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json"))
lengths = {r["name"]: r["length"] for r in assembly["top_level_region"] if r["coord_system"] == "chromosome"}
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
density = {}
for ch in CHROMS:
    L = lengths[ch]
    n_genes = 0
    for s in range(1, L + 1, 5_000_000):
        e = min(L, s + 5_000_000 - 1)
        js = json.loads(get(f"https://grch37.rest.ensembl.org/overlap/region/human/{ch}:{s}-{e}?feature=gene;biotype=protein_coding;content-type=application/json"))
        n_genes += sum(1 for g in js if s <= g["start"] <= e)
        time.sleep(0.07)
    density[ch] = n_genes / (L / 1e6)
print("protein-coding genes per Mb:", {k: round(v, 1) for k, v in density.items()})


def chrom_of(hom):
    return hom.split("(")[0]


radial = {ch: [] for ch in CHROMS}
c19_inside = 0
for gsm, name, cell in cells:
    allp = np.vstack([x for v in cell.values() for _, x in v])
    centre = allp.mean(0)
    rms = np.sqrt(np.mean(np.sum((allp - centre) ** 2, axis=1)))
    per = {}
    for hom, v in cell.items():
        ch = chrom_of(hom)
        if ch in radial:
            per.setdefault(ch, []).extend(list(np.linalg.norm(np.array([x for _, x in v]) - centre, axis=1) / rms))
    for ch, vals in per.items():
        radial[ch].append(float(np.mean(vals)))
    if "19" in per and "18" in per and np.mean(per["19"]) < np.mean(per["18"]):
        c19_inside += 1
mean_rad = {ch: float(np.mean(v)) for ch, v in radial.items() if v}
rho = float(spearmanr([density[c] for c in mean_rad], [mean_rad[c] for c in mean_rad]).correlation) if mean_rad else 0.0
print("mean normalised radial position:", {k: round(v, 3) for k, v in mean_rad.items()})
gate("in the Dip-C single-cell structures, gene-rich chromosomes sit nearer the centre of the nucleus: chromosome 19 is "
     "more central than 18 in at least 5 of 6 cells, and radial position falls with gene density (Spearman rho <= -0.5)",
     len(cells) == 6 and c19_inside >= 5 and rho <= -0.5, f"19 inside 18 in {c19_inside} of {len(cells)} cells; rho {rho:.2f}")
gsm, name, cell = cells[0]
nucleus = {"source": f"Dip-C, GEO {gsm} ({name}); Tan, Xing, Chang, Li & Xie, Science 361, 924 (2018)",
           "assembly": "GRCh37 (hg19)", "bin": BIN,
           "homologs": {hom: {"bins": [b for b, _ in v], "xyz": [np.round(x, 3).tolist() for _, x in v]} for hom, v in cell.items()},
           "gene_density": {k: round(v, 2) for k, v in density.items()},
           "radial": {k: round(v, 3) for k, v in mean_rad.items()}, "rho_density_radial": round(rho, 2)}
(OUT / "nucleus.json").write_text(json.dumps(nucleus, separators=(",", ":")) + "\n")

v = {"results": [{"name": r["name"], "pass": r["pass"]} for r in results],
     "nucleosome": {"superhelix": nuc["superhelix"]},
     "complexes": {c["pdb"]: {"bend": round(c["bend"]), "published": c.get("published_bend")} for c in cx_out},
     "fibre": {"neighbour": round(np.mean(d1)), "next_but_one": round(np.mean(d2))},
     "nucleus": {"chr19_inside_18": c19_inside, "cells": len(cells), "rho": round(rho, 2)}}
(OUT / "validation.json").write_text(json.dumps(v, indent=1, sort_keys=True, default=float) + "\n")
(WORK / "scales_report.json").write_text(json.dumps({"results": results, "summary": summary}, indent=1, default=float))
sys.exit(0 if all(r["pass"] for r in results) else 1)
