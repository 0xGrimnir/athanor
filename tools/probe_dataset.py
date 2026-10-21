"""List a public dataset's files (ENA) and deposited sequences (NCBI) from CI.

Usage: python tools/probe_dataset.py PRJEB103048 PZ244049 PZ244060
"""
import json
import sys
import time
import urllib.request

proj, first, last = sys.argv[1], sys.argv[2], sys.argv[3]
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={a}&result={r}&fields={f}&format=tsv"
for result, fields in (("analysis", "analysis_accession,sample_accession,sample_alias,analysis_type,analysis_title,submitted_ftp,submitted_bytes,submitted_md5"),
                       ("read_run", "run_accession,sample_accession,sample_alias,library_strategy,library_source,library_selection,instrument_model,read_count,base_count,fastq_ftp,fastq_bytes,submitted_ftp,submitted_bytes")):
    try:
        txt = urllib.request.urlopen(ENA.format(a=proj, r=result, f=fields), timeout=60).read().decode()
    except Exception as e:  # noqa: BLE001
        txt = f"ERROR {e}"
    print(f"== ENA {result} ({len(txt.splitlines()) - 1} rows)")
    print(txt[:20000])

prefix, n0, n1 = first[:2], int(first[2:]), int(last[2:])
ids = ",".join(f"{prefix}{i}" for i in range(n0, n1 + 1))
url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={ids}&retmode=json"
time.sleep(0.5)
d = json.loads(urllib.request.urlopen(url, timeout=60).read())
print("== GenBank")
for uid in d["result"].get("uids", []):
    r = d["result"][uid]
    print(r.get("accessionversion"), r.get("slen"), r.get("title"))
