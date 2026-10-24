"""Probe: data for the multiscale views (Dip-C 3D genomes, PDB complexes)."""
import gzip, io, json, re, tarfile, urllib.request

def get(url, n=None):
    req = urllib.request.Request(url, headers={"User-Agent": "athanor-ci"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read(n) if n else r.read()

for gse in ("GSE117876", "GSE117874", "GSE117109"):
    try:
        fl = get(f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/{gse}/suppl/filelist.txt").decode()
        lines = fl.splitlines()
        print(gse, len(lines), "files")
        for l in lines[:12]:
            print("  ", l)
        pat = [l for l in lines if "3dg" in l.lower()]
        print("  3dg files:", len(pat))
        for l in pat[:25]:
            print("    ", l)
    except Exception as e:
        print(gse, "failed", e)
# GEO sample pages for GM12878 Dip-C (GSE117109 series matrix)
try:
    t = get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117109/matrix/GSE117109_series_matrix.txt.gz")
    txt = gzip.decompress(t).decode(errors="replace")
    for l in txt.splitlines():
        if l.startswith(("!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1", "!Series_title", "!Sample_supplementary_file_1")):
            print(l[:600])
except Exception as e:
    print("matrix failed", e)
