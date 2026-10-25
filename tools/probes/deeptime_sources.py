"""Deep Time step 1: which ancient-genome data can CI reach, how big is it, and can samtools read it remotely?

Lists the MPI EVA archaic-genome directories (BAMs, VCFs, filter beds), ENA file reports for candidate
ancient-human studies, tries a remote region read with samtools, and checks bioconda for psmc and angsd.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import urllib.request

UA = {"User-Agent": "athanor-probe"}


def get(url, n=None):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read(n) if n else r.read()


def listing(url):
    """(name, size) pairs from an Apache-style index."""
    try:
        html = get(url).decode("utf-8", "replace")
    except Exception as e:  # noqa: BLE001
        print(f"  !! {url}: {e}")
        return []
    out = []
    for m in re.finditer(r'<a href="([^"?/][^"]*)">[^<]*</a>\s*</td>\s*<td[^>]*>([^<]*)</td>\s*<td[^>]*>\s*([^<]*)</td>', html):
        out.append((m.group(1), m.group(3).strip()))
    if not out:   # pre-formatted listing
        for m in re.finditer(r'<a href="([^"?/][^"]*)">[^<]*</a>\s+(\S+ \S+)\s+(\S+)', html):
            out.append((m.group(1), m.group(3)))
    return out


def walk(url, depth=2, show=40):
    items = listing(url)
    print(f"\n== {url} ({len(items)} entries)")
    for name, size in items[:show]:
        print(f"  {name:70s} {size}")
    if len(items) > show:
        print(f"  ... {len(items) - show} more; last: {items[-1]}")
    if depth > 1:
        for name, _ in items:
            if name.endswith("/"):
                walk(url + name, depth - 1, show)


print("#### MPI EVA archaic genomes")
for d in ("neandertal/altai/", "neandertal/Vindija/", "neandertal/Chagyrskaya/", "denisova/", "ust_ishim/", "ust-ishim/",
          "neandertal/Mezmaiskaya1/", "neandertal/Denisova17/"):
    walk("http://cdna.eva.mpg.de/" + d, depth=3, show=30)
walk("http://cdna.eva.mpg.de/", depth=1, show=80)

print("\n#### ENA file reports")
F_RUN = "run_accession,sample_alias,library_strategy,fastq_bytes,submitted_ftp,submitted_bytes,read_count,base_count"
F_ANA = "analysis_accession,sample_alias,analysis_type,submitted_ftp,submitted_bytes"
for acc in ("PRJEB6622", "PRJEB6272", "PRJEB11450", "PRJEB7618", "PRJEB9783", "PRJEB21037", "PRJEB29700"):
    for result, fields in (("analysis", F_ANA), ("read_run", F_RUN)):
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result={result}&fields={fields}&format=tsv"
        try:
            rows = get(url).decode().splitlines()
        except Exception as e:  # noqa: BLE001
            print(f"{acc} {result}: {e}")
            continue
        print(f"\n{acc} {result}: {len(rows) - 1} rows")
        for r in rows[:8]:
            print("  " + r[:400])
        try:
            info = json.loads(get(f"https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=study_accession={acc}"
                                  f"&fields=study_title&format=json"))
            print("  title:", info[0]["study_title"] if info else "?")
        except Exception as e:  # noqa: BLE001
            print("  title:", e)

print("\n#### samtools remote reads")
cands = []
for base in ("http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/bam/",
             "http://cdna.eva.mpg.de/neandertal/Vindija/bam/", "http://cdna.eva.mpg.de/neandertal/Chagyrskaya/BAM/"):
    for name, size in listing(base):
        if name.endswith(".bam") and re.search(r"(^|[._])(21|chr21)[._]", name):
            cands.append(base + name)
print("chr21 BAMs:", cands)
for url in cands[:3]:
    for cmd in (["samtools", "view", "-H", url], ["samtools", "view", "-c", url, "21:30000000-30100000"],
                ["samtools", "view", url, "21:30000000-30000100"]):
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        print(" ".join(cmd[:3]), cmd[-1] if len(cmd) > 4 else "", "->", r.returncode)
        print("   out:", "\n   ".join(r.stdout.splitlines()[:6])[:1500])
        print("   err:", r.stderr[-400:])
for cmd in (["samtools", "--version"], ["bcftools", "--version"]):
    print(subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()[0])

print("\n#### bioconda")
mamba = os.environ.get("MAMBA_EXE", "micromamba")
for pkg in ("psmc", "angsd", "bwa-mem2", "minimap2", "msmc2", "smcpp", "pyrho", "admixtools", "yleaf", "haplogrep3",
            "mapdamage2", "pmdtools", "htsbox", "sequencetools", "eigensoft", "plink2", "atlas", "snpad", "contammix",
            "schmutzi", "hapcon"):
    r = subprocess.run([mamba, "search", "-c", "conda-forge", "-c", "bioconda", pkg, "--json"], capture_output=True, text=True,
                       timeout=300)
    try:
        res = json.loads(r.stdout)["result"]["pkgs"]
        vers = sorted({p["version"] for p in res if p["name"] == pkg})
        print(f"  {pkg:16s} {vers[-5:] if vers else 'NOT FOUND'}")
    except Exception as e:  # noqa: BLE001
        print(f"  {pkg:16s} ? {e} {r.stderr[-200:]}")
