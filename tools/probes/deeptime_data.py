"""Deep Time benchmark data: can CI fetch just the pieces it needs?

Checks range requests against the AADR's .geno (so one individual's genotypes
can be read without the 5 GB file), the ENA BAMs of a few 1240K-capture
individuals, and remote region reads from the published archaic BAMs and VCFs.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import urllib.request

BASE = "https://dataverse.harvard.edu/api/access/datafile"
AADR = {"geno": 13994829, "ind": 13994513, "snp": 13994514, "anno": 13994515}
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
WORK = "/tmp/dtprobe"
os.makedirs(WORK, exist_ok=True)


def head(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "athanor-probe"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, dict(r.headers)
    except Exception as e:  # noqa: BLE001
        return None, {"error": str(e)}


def ranged(url, start, length):
    req = urllib.request.Request(url, headers={"User-Agent": "athanor-probe",
                                               "Range": f"bytes={start}-{start + length - 1}"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.status, r.read()


print("#### AADR files")
for key, fid in AADR.items():
    url = f"{BASE}/{fid}"
    st, h = head(url)
    print(f"{key:5s} {st} size={h.get('Content-Length')} ranges={h.get('Accept-Ranges')} "
          f"type={h.get('Content-Type')} {h.get('error', '')}")

print("\n#### range request on the .geno (first 64 bytes, then an interior slice)")
try:
    st, data = ranged(f"{BASE}/{AADR['geno']}", 0, 64)
    print("status", st, "bytes", len(data), "header", data.split(b"\x00")[0][:60])
    st2, d2 = ranged(f"{BASE}/{AADR['geno']}", 48 + 5 * 308_254, 320)
    print("interior status", st2, "bytes", len(d2), "first bytes", d2[:16].hex())
except Exception as e:  # noqa: BLE001
    print("range request failed:", e)

print("\n#### small 1240K BAMs on ENA (PRJEB11450)")
fields = "run_accession,sample_alias,submitted_ftp,submitted_bytes,submitted_md5"
rows = urllib.request.urlopen(
    urllib.request.Request(f"{ENA}?accession=PRJEB11450&result=read_run&fields={fields}&format=tsv",
                           headers={"User-Agent": "athanor-probe"}), timeout=120).read().decode().splitlines()
hdr = rows[0].split("\t")
recs = []
for line in rows[1:]:
    f = line.split("\t")
    d = dict(zip(hdr, f))
    try:
        size = min(int(x) for x in d["submitted_bytes"].split(";") if x)
    except ValueError:
        continue
    recs.append((size, d["sample_alias"], d["submitted_ftp"], d.get("submitted_md5", "")))
recs.sort()
print(f"{len(recs)} runs; the ten smallest:")
for size, alias, ftp, md5 in recs[:10]:
    print(f"  {alias:8s} {size / 1e6:8.1f} MB  {ftp}  {md5[:40]}")

print("\n#### downloading the smallest and reading it")
size, alias, ftp, _ = recs[0]
url = "https://" + ftp.split(";")[0]
out = f"{WORK}/{alias}.bam"
try:
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "athanor-probe"}),
                                timeout=1200) as r, open(out, "wb") as fh:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            fh.write(b)
    print("downloaded", out, os.path.getsize(out), "bytes")
    for cmd in (["samtools", "view", "-H", out], ["samtools", "index", out],
                ["samtools", "idxstats", out], ["samtools", "flagstat", out]):
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        o = r.stdout.splitlines()
        print(" ".join(cmd[:2]), "->", r.returncode, "|",
              " ; ".join(o[:6] if cmd[1] != "idxstats" else [x for x in o if not x.startswith("*")][:4]),
              r.stderr[-200:])
    r = subprocess.run(["samtools", "view", out, "-c"], capture_output=True, text=True)
    print("reads:", r.stdout.strip())
except Exception as e:  # noqa: BLE001
    print("download failed:", e)

print("\n#### AADR annotation columns for these individuals")
try:
    anno = urllib.request.urlopen(urllib.request.Request(f"{BASE}/{AADR['anno']}",
                                                         headers={"User-Agent": "athanor-probe"}),
                                  timeout=600).read().decode("utf-8", "replace").splitlines()
    hdr = anno[0].split("\t")
    want = [i for i, h in enumerate(hdr) if re.search(r"genetic id|molecular sex|damage|angsd|1240k snpset|"
                                                      r"y haplogroup in terminal|mtdna haplogroup|mean coverage|"
                                                      r"assessment", h, re.I)]
    print("columns:", [(i, hdr[i][:60]) for i in want])
    names = {a for _, a, _, _ in recs[:10]}
    shown = 0
    for line in anno[1:]:
        f = line.split("\t")
        if f and f[0] in names or (f and f[0].split(".")[0] in names):
            print("  " + " | ".join(f[i][:28] if i < len(f) else "" for i in want))
            shown += 1
            if shown > 12:
                break
    print("rows shown:", shown, "of", len(anno) - 1, "annotation rows")
except Exception as e:  # noqa: BLE001
    print("anno failed:", e)

print("\n#### archaic genomes: remote region reads")
ALT_BAM = ("http://cdna.eva.mpg.de/neandertal/Vindija/bam/Pruefer_etal_2017/Altai_realigned/",
           "http://cdna.eva.mpg.de/neandertal/Vindija/bam/Pruefer_etal_2017/Vindija33.19/")
for base in ALT_BAM:
    try:
        html = urllib.request.urlopen(urllib.request.Request(base, headers={"User-Agent": "athanor-probe"}),
                                      timeout=120).read().decode("utf-8", "replace")
        names = re.findall(r'<a href="([^"?/][^"]*)">', html)
        print(base, "->", names[:8], "..." if len(names) > 8 else "")
    except Exception as e:  # noqa: BLE001
        print(base, "listing failed:", e)
for url in ("http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz",
            "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr21_mq25_mapab100.vcf.gz"):
    r = subprocess.run(["bcftools", "view", "-H", "-r", "21:30000000-30002000", url], capture_output=True,
                       text=True, timeout=1800)
    lines = r.stdout.splitlines()
    print(url.split("/")[-2], "->", r.returncode, len(lines), "records")
    for line in lines[:4]:
        print("   ", line[:200])
    if r.returncode:
        print("    err:", r.stderr[-300:])
print("\n#### filter-bed README")
for url in ("http://cdna.eva.mpg.de/neandertal/Vindija/FilterBed/README",
            "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/README",
            "http://cdna.eva.mpg.de/neandertal/Vindija/bam/Pruefer_etal_2017/README"):
    try:
        print(f"-- {url}\n" + urllib.request.urlopen(urllib.request.Request(
            url, headers={"User-Agent": "athanor-probe"}), timeout=60).read().decode("utf-8", "replace"))
    except Exception as e:  # noqa: BLE001
        print(url, "failed:", e)
