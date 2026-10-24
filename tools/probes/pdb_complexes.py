"""Probe: chain layout of candidate protein-DNA complexes and nucleosome structures."""
import json, urllib.request, collections
def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "athanor-ci"}), timeout=120) as r:
        return r.read()
NUC = {"DA", "DC", "DG", "DT", "DU", "A", "C", "G", "U"}
for pid in ("1KX5", "1ZBB", "6HKT", "1YTB", "1CGP", "1IHF", "1AAY", "1TSR", "1RVA", "1EFA", "1CDW", "3LZ0", "1A1T", "1HCQ", "3DNB"):
    try:
        txt = get(f"https://files.rcsb.org/download/{pid}.pdb").decode()
    except Exception as e:
        print(pid, "no PDB-format file:", e); continue
    chains = collections.OrderedDict()
    for l in txt.splitlines():
        if l.startswith("ENDMDL"): break
        if l.startswith(("ATOM", "HETATM")):
            c, rn, num = l[21], l[17:20].strip(), l[22:27]
            d = chains.setdefault(c, {"res": collections.OrderedDict(), "nuc": 0, "atoms": 0})
            if num not in d["res"]:
                d["res"][num] = rn
                d["nuc"] += rn in NUC
            d["atoms"] += 1
    comp = [l[10:].strip() for l in txt.splitlines() if l.startswith("COMPND") and ("MOLECULE" in l or "CHAIN" in l)]
    print(f"== {pid}: {' | '.join(comp)[:400]}")
    for c, d in chains.items():
        seq = "".join(r[-1] for r in d["res"].values() if r in NUC)
        print(f"   chain {c}: {len(d['res'])} residues, {d['nuc']} nucleotides, {d['atoms']} atoms  {seq[:80]}")
