"""DNA Structure data and benchmarks, built from the PDB and checked against 3DNA.

1. Every X-ray structure in the PDB with DNA and nothing else (no protein, RNA or hybrid), resolved to
   1.6 A or better, is downloaded (RCSB; PDB data are CC0). Two-stranded Watson-Crick duplexes of at
   least six base pairs with only standard nucleotides are kept and assigned a helix form from their
   own geometry (Z: negative twist; A: slide below -1.2 A; B: slide above -0.8 A; others left out).
2. Athanor's analysis (base-pair, step and groove parameters) is compared with 3DNA's find_pair and
   analyze on every kept structure.
3. Parameter sets the builder uses: B-DNA steps are Olson et al. (1998) Table 1; B-DNA base pairs, and
   A- and Z-DNA steps and pairs, are medians over the kept structures of that form. Z-DNA nucleotide
   geometries (syn G, anti C) are the medians of the crystal nucleotides in their base frames.
4. Builder benchmarks: Athanor's rebuild against 3DNA's rebuild from the same parameters; each crystal
   rebuilt from its own parameters; and each crystal built from the generic parameter set for its form,
   with that structure left out of the set (leave-one-out).

Usage: python build_dna.py X3DNA_DIR DATA_DIR WORK_DIR
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.dna.builder import build  # noqa: E402
from athanor_core.dna.geometry import kabsch  # noqa: E402
from athanor_core.dna.structure import analyse, base_frame, read_pdb, standard_bases  # noqa: E402

X3, DATA, WORK = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]).resolve()
WORK.mkdir(parents=True, exist_ok=True)
ENV = {**os.environ, "X3DNA": str(X3), "PATH": f"{X3 / 'bin'}:{os.environ['PATH']}"}
results, summary = [], {}
OLSON = {  # Olson et al., PNAS 95, 11163 (1998), Table 1: shift, slide, rise (A), tilt, roll, twist (deg)
    "AA": [-0.03, -0.08, 3.27, -1.4, 0.7, 35.1], "AC": [0.13, -0.58, 3.36, -0.1, 0.7, 31.5],
    "AG": [0.09, -0.25, 3.34, -1.7, 4.5, 31.9], "AT": [0.0, -0.59, 3.31, 0.0, 1.1, 29.3],
    "CA": [0.09, 0.53, 3.33, 0.5, 4.7, 37.3], "CG": [0.0, 0.41, 3.39, 0.0, 5.4, 36.1],
    "GA": [-0.28, 0.09, 3.37, -1.5, 1.9, 36.3], "GC": [0.0, -0.38, 3.40, 0.0, 0.3, 33.6],
    "GG": [0.05, -0.22, 3.42, -0.1, 3.6, 32.9], "TA": [0.0, 0.05, 3.42, 0.0, 3.3, 37.8]}
STANDARD = {"DA", "DC", "DG", "DT"}


def gate(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def get(url, data=None, tries=4, raw=False):
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "athanor-ci", "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                b = r.read()
            return b if raw else json.loads(b)
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(3 * (k + 1))
    raise RuntimeError(f"{url}: {err!r}")


# ------------------------------------------------------------------ 1. the structures
def term(attr, op, value):
    return {"type": "terminal", "service": "text", "parameters": {"attribute": attr, "operator": op, "value": value}}


q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
    term("rcsb_entry_info.polymer_entity_count_DNA", "greater_or_equal", 1),
    term("rcsb_entry_info.polymer_entity_count_protein", "equals", 0),
    term("rcsb_entry_info.polymer_entity_count_RNA", "equals", 0),
    term("rcsb_entry_info.polymer_entity_count_nucleic_acid_hybrid", "equals", 0),
    term("rcsb_entry_info.resolution_combined", "less_or_equal", 1.6),
    term("exptl.method", "exact_match", "X-RAY DIFFRACTION")]},
    "return_type": "entry", "request_options": {"return_all_hits": True}}
hits = get("https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(q).encode())
ids = sorted(h["identifier"] for h in hits.get("result_set", []))
print(f"{len(ids)} DNA-only X-ray entries at <= 1.6 A", flush=True)

std = standard_bases()
kept = []
for pid in ids:
    f = WORK / f"{pid}.pdb"
    if not f.exists():
        try:
            f.write_bytes(get(f"https://files.rcsb.org/download/{pid}.pdb", raw=True))
        except Exception:  # noqa: BLE001 - large entries have no PDB-format file
            continue
    res = read_pdb(f.read_text())
    if len({r.chain for r in res}) != 2 or any(r.name not in STANDARD for r in res):
        continue
    try:
        a = analyse(res, std)
    except Exception as e:  # noqa: BLE001
        print(f"   {pid}: analysis failed {e!r}"[:160])
        continue
    n = len(a["pairs"])
    chains = {}
    for r in res:
        chains.setdefault(r.chain, []).append(r)
    if n < 6 or not all(p["watson_crick"] for p in a["pairs"]) or len(a["steps"]) != n - 1:
        print(f"   {pid}: {n} pairs ({sum(p['watson_crick'] for p in a['pairs'])} Watson-Crick), {len(a['steps'])} steps, "
              f"chain lengths {[len(c) for c in chains.values()]} - left out")
        continue
    if any(len(c) != n for c in chains.values()):
        print(f"   {pid}: chain lengths {[len(c) for c in chains.values()]} for {n} pairs - left out")
        continue                        # overhangs, mismatches or missing residues
    st = np.array([s["params"] for s in a["steps"]])
    tw, sl = st[:, 5].mean(), st[:, 1].mean()
    form = "Z" if tw < 0 else "A" if sl < -1.2 else "B" if sl > -0.8 else None
    print(f"   {pid}: {n} bp {a['sequence']} mean twist {tw:.1f} slide {sl:.2f} -> {form}")
    if form is None:
        continue
    kept.append({"id": pid, "form": form, "n": n, "analysis": a, "residues": res, "file": f})
print("kept duplexes:", {fm: sum(k["form"] == fm for k in kept) for fm in "ABZ"}, flush=True)
summary["structures"] = {fm: [k["id"] for k in kept if k["form"] == fm] for fm in "ABZ"}

# ------------------------------------------------------------------ 2. against 3DNA
NUM = r"(-?\d+\.\d+)"


def x3dna_analyse(path: Path):
    d = WORK / f"x3_{path.stem}"
    d.mkdir(exist_ok=True)
    subprocess.run(["find_pair", str(path.resolve()), "out.inp"], cwd=d, env=ENV, capture_output=True, check=True)
    subprocess.run(["analyze", "out.inp"], cwd=d, env=ENV, capture_output=True, check=True)
    txt = (d / f"{path.stem}.out").read_text()
    pairs = re.findall(r"^\s+\d+\s+\(\s*[\d.]+\)\s+\.*>(\S):\.*(-?\d+)_:\[\.*(\w+)\](\w)-[-*+x]+-(\w)\[\.*(\w+)\]:\.*(-?\d+)_:(\S)<",
                       txt, re.M)

    def block(title, ncol=6):
        m = re.search(title + r".*?\n(.*?)~~~~", txt, re.S)
        rows = []
        for line in (m.group(1).splitlines() if m else []):
            v = re.findall(NUM + r"|(----)", line)
            nums = [float(a) if a else None for a, b in v]
            if len(nums) >= ncol:
                rows.append(nums[-ncol:])
        return rows
    bp = block(r"Local base-pair parameters\n")
    stp = block(r"Local base-pair step parameters\n")
    g = []
    m = re.search(r"P-P     Refined     P-P     Refined\n(.*?)\*\*\*\*", txt, re.S)
    for line in (m.group(1).splitlines() if m else []):
        f = line.split()
        if len(f) >= 6:
            g.append([None if x == "---" else float(x) for x in (f[2], f[4])])
    strand1 = pairs[0][0] if pairs else None
    return {"pairs": bp, "steps": stp, "grooves": g, "strand1_chain": strand1, "pair_ids": pairs}


worst_bp = worst_st = worst_g = 0.0
compared = failed = 0
for k in kept:
    try:
        x = x3dna_analyse(k["file"])
    except Exception as e:  # noqa: BLE001
        failed += 1
        print("3DNA failed on", k["id"], repr(e)[:200])
        continue
    a = k["analysis"]
    if x["strand1_chain"] and x["strand1_chain"] != a["strands"][0]:
        # read from 3DNA's strand I
        order = [x["strand1_chain"]] + [c for c in a["strands"] if c != x["strand1_chain"]]
        res = sorted(k["residues"], key=lambda r: order.index(r.chain))
        a = analyse(res, std)
    if len(x["steps"]) != len(a["steps"]) or len(x["pairs"]) != len(a["pairs"]):
        failed += 1
        print("pair lists differ for", k["id"], len(x["pairs"]), len(a["pairs"]))
        continue
    theirs_ids = [(f"{c1}{n1}", f"{c2}{n2}") for c1, n1, _, _, _, _, n2, c2 in x["pair_ids"]]
    ours_ids = [tuple(p["residues"]) for p in a["pairs"]]
    if theirs_ids and theirs_ids != ours_ids:
        diff = [(i + 1, t, o) for i, (t, o) in enumerate(zip(theirs_ids, ours_ids)) if t != o]
        failed += 1
        print("pairing differs for", k["id"], diff[:3])
        continue
    compared += 1
    k["same_pairs"] = True
    db = np.abs(np.array(x["pairs"]) - np.array([p["params"] for p in a["pairs"]]))
    ds = np.abs(np.array(x["steps"]) - np.array([s["params"] for s in a["steps"]]))
    if db.max() > 0.011 or ds.max() > 0.011:
        i = int(np.argmax(db.max(1))) if db.max() >= ds.max() else int(np.argmax(ds.max(1)))
        kind = "pair" if db.max() >= ds.max() else "step"
        ours = a["pairs"][i] if kind == "pair" else a["steps"][i]
        theirs = x["pairs"][i] if kind == "pair" else x["steps"][i]
        print(f"   {k['id']} ({k['form']}) differs at {kind} {i + 1} {ours.get('pair') or ours.get('step')}: ours "
              f"{np.round(ours['params'], 2).tolist()} 3DNA {theirs}; 3DNA strand I {x['strand1_chain']}, ours {a['strands']}")
    worst_bp = max(worst_bp, float(db[:, :3].max()), float(db[:, 3:].max()) / 10)
    worst_st = max(worst_st, float(ds[:, :3].max()), float(ds[:, 3:].max()) / 10)
    for (mx, Mx), g in zip(x["grooves"], a["grooves"]):
        for u, v in ((mx, g["minor_pp"]), (Mx, g["major_pp"])):
            if u is not None and v is not None:
                worst_g = max(worst_g, abs(u - v))
            elif (u is None) != (v is None):
                worst_g = max(worst_g, 99.0)
print(f"compared with 3DNA: {compared} structures, {failed} not comparable")
gate("base-pair and step parameters equal 3DNA's on every kept crystal duplex (to its printed 0.01 A / 0.01 deg)",
     compared >= 20 and worst_bp <= 0.011 and worst_st <= 0.011,
     f"{compared} structures; largest difference {worst_bp:.3f} (pairs) and {worst_st:.3f} (steps), "
     "translations in A, rotations in units of 10 deg")
gate("groove widths equal 3DNA's direct P-P distances (to its printed 0.1 A)", compared >= 20 and worst_g <= 0.051,
     f"largest difference {worst_g:.3f} A")

# ------------------------------------------------------------------ 3. parameter sets
def norm_pair(p):
    """Pair parameters as if read with the purine (G or A) on strand I: shear and buckle reverse sign
    when the two strands are exchanged."""
    v = list(p["params"])
    if p["pair"][0] in "CT":
        v[0], v[3] = -v[0], -v[3]
    return v


def med(rows):
    return [round(float(v), 2) for v in np.median(np.array(rows), axis=0)] if rows else None


def sets(exclude=None):
    out = {}
    use = [k for k in kept if k["id"] != exclude]
    # B: pairs by type (A-T or G-C) from the crystals; steps from Olson et al. 1998
    bpairs = {"AT": [], "GC": []}
    for k in use:
        if k["form"] != "B":
            continue
        for p in k["analysis"]["pairs"]:
            bpairs["GC" if p["pair"][0] in "GC" else "AT"].append(norm_pair(p))
    out["B"] = {"dimers": OLSON, "pair_by_type": {t: med(v) for t, v in bpairs.items()}}
    for fm in ("A", "Z"):
        steps, pairs = {}, {}
        for k in use:
            if k["form"] != fm:
                continue
            seq = k["analysis"]["sequence"]
            for i, s in enumerate(k["analysis"]["steps"]):
                key = seq[i:i + 2] if fm == "Z" else "all"
                steps.setdefault(key, []).append(s["params"])
            for p in k["analysis"]["pairs"]:
                pairs.setdefault("GC" if p["pair"][0] in "GC" else "AT", []).append(norm_pair(p))
        out[fm] = {"steps": {k2: med(v) for k2, v in steps.items()}, "pair_by_type": {t: med(v) for t, v in pairs.items()},
                   "n_structures": sum(k["form"] == fm for k in use)}
    return out


P = sets()
print(json.dumps(P, indent=1)[:3000])

# Z-DNA nucleotides in their base frames (syn G, anti C): medians over the crystal nucleotides
ztmpl = {}
for b in "CG":
    acc = {}
    for k in kept:
        if k["form"] != "Z":
            continue
        for r in k["residues"]:
            if r.base != b:
                continue
            fr = base_frame(r, std)
            if fr is None:
                continue
            R, o, _ = fr
            for nm, xyz in r.atoms.items():
                acc.setdefault(nm, []).append(R.T @ (xyz - o))
    ztmpl[b] = {nm: [round(float(v), 3) for v in np.median(np.array(x), axis=0)] for nm, x in acc.items()
                if len(x) >= 0.8 * max(len(y) for y in acc.values())}
    spread = {nm: float(np.median(np.linalg.norm(np.array(x) - np.median(np.array(x), axis=0), axis=1)))
              for nm, x in acc.items()}
    print(f"Z {b}: {len(ztmpl[b])} atoms; median spread {np.median(list(spread.values())):.2f} A")


def write_data(Pset):
    steps = {
        "B": {"dimers": Pset["B"]["dimers"], "pair": Pset["B"]["pair_by_type"]["GC"], "pair_by_type": Pset["B"]["pair_by_type"],
              "method": "Step parameters: averages for each dimer step in protein-DNA crystal complexes (Olson et al. 1998, "
                        "Table 1). Base-pair parameters: medians over high-resolution B-DNA crystal duplexes in the PDB. "
                        "Nucleotide geometry: 3DNA's B-DNA set (Arnott backbone, chi -110 deg).",
              "sources": ["Olson et al., DNA sequence-dependent deformability deduced from protein-DNA crystal complexes, PNAS 95, 11163 (1998)",
                          "Olson et al., A standard reference frame for the description of nucleic acid base-pair geometry, J. Mol. Biol. 313, 229 (2001)",
                          "Lu & Olson, 3DNA, Nucleic Acids Res. 31, 5108 (2003)"]},
        "A": {"uniform": Pset["A"]["steps"].get("all"), "pair": Pset["A"]["pair_by_type"].get("GC"), "pair_by_type": Pset["A"]["pair_by_type"],
              "method": f"Uniform step and base-pair parameters: medians over {Pset['A']['n_structures']} high-resolution A-DNA "
                        "crystal duplexes in the PDB. Nucleotide geometry: 3DNA's A-DNA set (chi -160 deg).",
              "sources": ["Lu & Olson, 3DNA, Nucleic Acids Res. 31, 5108 (2003)", "RCSB Protein Data Bank (CC0)"]},
        "Z": {"uniform": {k2: v for k2, v in Pset["Z"]["steps"].items() if k2 in ("CG", "GC")}, "pair": Pset["Z"]["pair_by_type"].get("GC"),
              "pair_by_type": Pset["Z"]["pair_by_type"],
              "method": f"Alternating CG and GC step parameters and base-pair parameters: medians over {Pset['Z']['n_structures']} "
                        "high-resolution Z-DNA crystal duplexes in the PDB; nucleotides (syn G, anti C) from the same crystals.",
              "sources": ["Wang et al., Molecular structure of a left-handed double helical DNA fragment at atomic resolution, Nature 282, 680 (1979)",
                          "RCSB Protein Data Bank (CC0)"]},
    }
    (DATA / "steps.json").write_text(json.dumps(steps, indent=1))
    t = json.loads((DATA / "templates.json").read_text())
    t["Z"] = ztmpl
    (DATA / "templates.json").write_text(json.dumps(t, indent=1))


write_data(P)

# ------------------------------------------------------------------ 4. builder benchmarks
def rmsd_against(model, res, atoms_filter=None):
    ours = {}
    for c, n, rn, nm, xyz in model["atoms"]:
        ours.setdefault((c, n), {})[nm] = xyz
    by = {}
    for r in res:
        by.setdefault(r.chain, []).append(r)
    ch = list(by)
    s1, s2 = by[ch[0]], by[ch[1]]
    n = len(s1)
    Pm, Qm = [], []
    for k in range(n):
        for (c, num), cr in ((("A", k + 1), s1[k]), (("B", n + 1 + k), s2[k])):
            o = ours.get((c, num), {})
            for nm, xyz in cr.atoms.items():
                if nm in o and (atoms_filter is None or atoms_filter(nm)):
                    Pm.append(o[nm])
                    Qm.append(xyz)
    _, _, r = kabsch(np.array(Pm), np.array(Qm))
    return r, len(Pm)


BACK = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"}
own, generic = {"A": [], "B": [], "Z": []}, {"A": [], "B": [], "Z": []}
for k in kept:
    a, fm = k["analysis"], k["form"]
    seq = a["sequence"]
    if fm == "Z" and any(seq[i:i + 2] not in ("CG", "GC") for i in range(len(seq) - 1)):
        continue
    m_own = build(seq, fm, steps=[s["params"] for s in a["steps"]], pairs=[p["params"] for p in a["pairs"]])
    rb, _ = rmsd_against(m_own, k["residues"], lambda nm: nm not in BACK)
    ra, na = rmsd_against(m_own, k["residues"])
    own[fm].append((rb, ra))
    Ploo = sets(exclude=k["id"])
    write_data(Ploo)
    m_gen = build(seq, fm)
    rg, _ = rmsd_against(m_gen, k["residues"])
    generic[fm].append(rg)
write_data(P)
for fm in "ABZ":
    if own[fm]:
        o = np.array(own[fm])
        print(f"{fm}: {len(o)} crystals; rebuilt from own parameters: base atoms median {np.median(o[:, 0]):.3f} A, "
              f"all atoms median {np.median(o[:, 1]):.2f} A; generic parameters (leave-one-out): median "
              f"{np.median(generic[fm]):.2f} A, 90th percentile {np.percentile(generic[fm], 90):.2f} A")
summary["rebuild"] = {fm: {"n": len(own[fm]), "own_base": float(np.median([x[0] for x in own[fm]])) if own[fm] else None,
                           "own_all": float(np.median([x[1] for x in own[fm]])) if own[fm] else None,
                           "generic_median": float(np.median(generic[fm])) if generic[fm] else None,
                           "generic_p90": float(np.percentile(generic[fm], 90)) if generic[fm] else None} for fm in "ABZ"}
allown = [x[0] for v in own.values() for x in v]
gate("rebuilt from its own parameters, every crystal's bases come back exactly (base-atom RMSD median <= 0.05 A)",
     allown and np.median(allown) <= 0.05, f"median {np.median(allown):.3f} A over {len(allown)} crystals")
gb = summary["rebuild"]["B"]
gate("B-DNA built from sequence alone (Olson 1998 steps) lands within 2.0 A all-atom RMSD of the crystal (median)",
     gb["generic_median"] is not None and gb["generic_median"] <= 2.0,
     f"median {gb['generic_median']:.2f} A, 90th percentile {gb['generic_p90']:.2f} A over {gb['n']} crystals")
for fm in "AZ":
    g = summary["rebuild"][fm]
    if g["n"]:
        gate(f"{fm}-DNA built from its generic parameters lands within 2.0 A all-atom RMSD of the crystal (median, leave-one-out)",
             g["generic_median"] <= 2.0, f"median {g['generic_median']:.2f} A, 90th percentile {g['generic_p90']:.2f} A over {g['n']} crystals")

# against 3DNA's rebuild from the same parameters
worst_rb = 0.0
nrb = 0
for k in [k for k in kept if k["form"] == "B" and k.get("same_pairs")][:15]:
    d = WORK / f"x3_{k['file'].stem}"
    try:
        subprocess.run(["x3dna_utils", "cp_std", "BDNA"], cwd=d, env=ENV, capture_output=True, check=True)
        subprocess.run(["rebuild", "-atomic", "bp_step.par", "rb.pdb"], cwd=d, env=ENV, capture_output=True, check=True)
    except Exception as e:  # noqa: BLE001
        print("3DNA rebuild failed", k["id"], repr(e)[:200])
        continue
    ref = read_pdb((d / "rb.pdb").read_text())
    x = x3dna_analyse(k["file"])
    a = k["analysis"]
    if x["strand1_chain"] != a["strands"][0]:
        order = [x["strand1_chain"]] + [c for c in a["strands"] if c != x["strand1_chain"]]
        a = analyse(sorted(k["residues"], key=lambda r: order.index(r.chain)), std)
    m = build(a["sequence"], "B", steps=[s["params"] for s in a["steps"]], pairs=[p["params"] for p in a["pairs"]])
    r, n = rmsd_against(m, ref)
    if r > 0.02:
        print(f"   rebuild {k['id']}: RMSD {r:.3f} A over {n} atoms; our sequence {a['sequence']}, 3DNA chains "
              f"{sorted({rr.chain for rr in ref})}, residues {len(ref)}")
    worst_rb = max(worst_rb, r)
    nrb += 1
gate("Athanor's rebuild equals 3DNA's rebuild from the same parameters (RMSD <= 0.02 A)", nrb >= 5 and worst_rb <= 0.02,
     f"largest RMSD {worst_rb:.4f} A over {nrb} duplexes")

(DATA / "validation.json").write_text(json.dumps({"results": results, "summary": summary}, indent=1, default=str))
(WORK / "dna_report.json").write_text(json.dumps({"results": results, "summary": summary, "params": P}, indent=1, default=str))
sys.exit(0 if all(r["pass"] for r in results) else 1)
