"""Download the published data for the reproduction test (runs on the CI host).

Dataset: woolly mammoth ivory from Hohle Fels (Scientific Reports, 2026;
https://www.nature.com/articles/s41598-026-46761-x). Reads: ENA project
PRJEB103048 (Illumina HiSeq 4000, 1x75, partial-UDG double-stranded
libraries). Published mitogenomes: GenBank PZ244049-PZ244060.

For each sample whose GenBank record is a complete mitogenome, every
deposited BAM of reads realigned to the mitochondrial genome
("MT_realigned", shotgun and capture runs) is downloaded and checked
against ENA's MD5, together with the GenBank sequence.

Usage: python fetch_reproduction.py OUT_DIR
"""
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

PROJECT = "PRJEB103048"
GENBANK = [f"PZ2440{i}" for i in range(49, 61)]
ENA = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={a}&result=read_run"
       "&fields=run_accession,sample_alias,library_strategy,read_count,submitted_ftp,submitted_md5,submitted_bytes"
       "&format=tsv")
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ids}&rettype=fasta&retmode=text"


def get(url: str, tries: int = 4) -> bytes:
    for k in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=120).read()
        except Exception:  # noqa: BLE001
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))


def main(out: Path):
    out.mkdir(parents=True, exist_ok=True)
    fasta = get(EFETCH.format(ids=",".join(GENBANK))).decode()
    records, cur = {}, None
    for line in fasta.splitlines():
        if line.startswith(">"):
            cur = line[1:]
            records[cur] = []
        elif cur:
            records[cur].append(line.strip())
    genomes = {}
    for header, seq in records.items():
        acc = header.split()[0]
        if "complete genome" not in header or " isolate " not in header:
            continue
        isolate = header.split(" isolate ")[1].split()[0]
        genomes[isolate] = {"accession": acc, "title": header, "sequence": "".join(seq)}
    runs = [dict(zip(lines[0].split("\t"), l.split("\t")))
            for lines in [get(ENA.format(a=PROJECT)).decode().splitlines()] for l in lines[1:] if l.strip()]
    samples = []
    for isolate, g in sorted(genomes.items()):
        files = []
        for r in runs:
            if r["sample_alias"] != isolate:
                continue
            for url, md5, size in zip(r["submitted_ftp"].split(";"), r["submitted_md5"].split(";"),
                                      r["submitted_bytes"].split(";")):
                if not url.endswith(".bam") or "MT_realigned" not in url:
                    continue
                dest = out / isolate / f"{r['run_accession']}.bam"
                dest.parent.mkdir(parents=True, exist_ok=True)
                data = get("https://" + url)
                got = hashlib.md5(data).hexdigest()
                if got != md5:
                    raise SystemExit(f"MD5 mismatch for {url}: {got} != {md5}")
                dest.write_bytes(data)
                files.append({"run": r["run_accession"], "strategy": r["library_strategy"],
                              "reads": int(r["read_count"] or 0), "url": "https://" + url, "md5": md5,
                              "bytes": int(size), "path": str(dest.relative_to(out))})
        if files:
            (out / isolate / "deposited.fa").write_text(f">{g['title']}\n{g['sequence']}\n")
            samples.append({"isolate": isolate, "genbank": g["accession"], "title": g["title"],
                            "length": len(g["sequence"]), "files": files})
            print(isolate, g["accession"], len(g["sequence"]), [(f["run"], f["strategy"], f["reads"]) for f in files])
    (out / "manifest.json").write_text(json.dumps({"project": PROJECT, "samples": samples}, indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
