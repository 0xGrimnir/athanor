"""Probe: sources for the remaining DNA Structure features."""
import json, os, re, subprocess, urllib.request, glob, zipfile, io

def get(url, n=None):
    req = urllib.request.Request(url, headers={"User-Agent": "athanor-ci"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read(n) if n else r.read()

# psi4 S22 database, dftd3
import psi4
base = os.path.dirname(psi4.__file__)
cands = glob.glob(os.path.join(base, "..", "..", "..", "share", "psi4", "databases", "S22*.py")) + glob.glob(os.path.join(os.environ.get("CONDA_PREFIX", ""), "share", "psi4", "databases", "S22*.py"))
print("psi4", psi4.__version__, sorted(set(os.path.basename(c) for c in cands)))
if cands:
    s = open(sorted(cands)[0]).read()
    print(s[:3000])
    m = re.search(r"BIND_S22A\s*=\s*\{.*?\}", s, re.S) or re.search(r"BIND_S22B\s*=\s*\{.*?\}", s, re.S)
    print(m.group(0)[:1500] if m else "no BIND table found")
    print([x for x in re.findall(r"^(BIND_\w+)\s*=", s, re.M)])
for exe in ("dftd3", "s-dftd3"):
    print(exe, subprocess.run(["which", exe], capture_output=True, text=True).stdout.strip())
try:
    import dftd3; print("dftd3 python", dftd3.__version__)
except Exception as e:
    print("dftd3 python:", e)
# Zenodo Deep DNAshape processed data
try:
    meta = json.loads(get("https://zenodo.org/api/records/10403307"))
    print("zenodo", meta["metadata"]["title"], meta["metadata"].get("license"), [(f["key"], f["size"]) for f in meta["files"]])
    z = zipfile.ZipFile(io.BytesIO(get(meta["files"][0]["links"]["self"])))
    for i in z.infolist()[:80]:
        print("  ", i.filename, i.file_size)
    for i in z.infolist():
        if not i.is_dir() and i.file_size < 50_000_000:
            with z.open(i) as fh:
                print("HEAD", i.filename, fh.read(300))
            break
except Exception as e:
    print("zenodo failed", e)
# GEO GSE117876 (Dip-C)
try:
    html = get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117876/suppl/").decode()
    print("GEO suppl:", re.findall(r'href="([^"]+)"', html)[:40])
except Exception as e:
    print("GEO failed", e)
# PDB entries
for pid in ("1KX5", "1ZBB", "6HKT", "1YTB", "1CGP", "1IHF", "1AAY", "1RVA", "1EFA", "1TSR", "3LZ0"):
    try:
        d = json.loads(get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}"))
        print(pid, d["struct"]["title"][:90], "|", d["rcsb_entry_info"].get("resolution_combined"), d["rcsb_entry_info"].get("deposited_atom_count"),
              d["exptl"][0]["method"], "| pdb format:", "yes" if d["rcsb_entry_info"].get("pdbx_structure_determination_methodology") else "?")
    except Exception as e:
        print(pid, "failed", e)
