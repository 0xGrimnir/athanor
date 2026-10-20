"""Fetch the reference mitochondrial genomes from NCBI (E-utilities) and record
exactly which accession versions were used, with SHA-256 of each file."""
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

REFS = {
    # accession: (organism, role)
    "NC_007596": ("Mammuthus primigenius (woolly mammoth)", "truth source for the ancient-DNA demo"),
    "NC_000934": ("Loxodonta africana (African savanna elephant)", "mapping reference"),
    "NC_005129": ("Elephas maximus (Asian elephant)", "mapping reference"),
    "NC_012920": ("Homo sapiens (revised Cambridge Reference Sequence)", "contaminant source / human mtDNA"),
}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype={rt}&retmode=text"

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
manifest = {"source": "NCBI nuccore via E-utilities efetch", "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "genomes": []}
for acc, (organism, role) in REFS.items():
    files = {}
    for rt, ext in (("fasta", "fa"), ("gb", "gb")):
        for attempt in range(5):
            try:
                data = urllib.request.urlopen(EUTILS.format(acc=acc, rt=rt), timeout=60).read()
                break
            except Exception:
                time.sleep(3 * (attempt + 1))
        else:
            raise SystemExit(f"could not fetch {acc} ({rt})")
        (out / f"{acc}.{ext}").write_bytes(data)
        files[ext] = hashlib.sha256(data).hexdigest()
        time.sleep(0.4)  # NCBI: at most 3 requests per second without a key
    header = (out / f"{acc}.fa").read_text().splitlines()[0]
    version = header[1:].split()[0]
    seq = "".join((out / f"{acc}.fa").read_text().splitlines()[1:])
    manifest["genomes"].append({"accession": acc, "version": version, "organism": organism, "role": role,
                                "length": len(seq), "header": header[1:], "sha256": files})
(out / "manifest.json").write_text(json.dumps(manifest, indent=1))
print(json.dumps(manifest, indent=1))
