"""Reference sequence around each Personal Genome trait variant, so the app can show a variant in a helix.

For every variant in tools/pgdata/traits_spec.py: its GRCh37 position and alleles (Ensembl GRCh37 REST,
variation endpoint) and the forward-strand reference sequence 10 bases either side (sequence/region
endpoint). The base at the variant's position must equal the reference allele Ensembl reports, or the
script fails. GRCh37 is frozen, so the output is stable and CI checks that the committed file regenerates.

Usage: python flanks.py OUT_JSON
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pgdata"))
from traits_spec import VARIANTS  # noqa: E402

FLANK = 10
API = "https://grch37.rest.ensembl.org"


def get(url, tries=4):
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json", "User-Agent": "athanor-ci"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(3 * (k + 1))
    raise RuntimeError(f"{url}: {err!r}")


out, fails = {}, []
for rs in sorted(VARIANTS):
    d = get(f"{API}/variation/human/{rs}")
    maps = [m for m in d.get("mappings", []) if re.fullmatch(r"\d+|X|Y", str(m.get("seq_region_name")))]
    if len(maps) != 1:
        fails.append(f"{rs}: {len(maps)} primary-assembly mappings")
        continue
    m = maps[0]
    alleles = m["allele_string"].split("/")
    if m["start"] != m["end"] or any(len(a) != 1 for a in alleles):
        print(f"{rs}: not a single-base variant ({m['allele_string']}) - left out")
        continue
    chrom, pos = m["seq_region_name"], int(m["start"])
    s = get(f"{API}/sequence/region/human/{chrom}:{pos - FLANK}..{pos + FLANK}:1")["seq"].upper()
    if len(s) != 2 * FLANK + 1 or s[FLANK] != alleles[0]:
        fails.append(f"{rs}: reference base {s[FLANK:FLANK + 1]} at {chrom}:{pos} is not Ensembl's reference allele {alleles[0]}")
        continue
    out[rs] = {"chrom": chrom, "pos37": pos, "ref": alleles[0], "alts": alleles[1:],
               "left": s[:FLANK], "right": s[FLANK + 1:]}
    print(f"{rs:12s} {chrom}:{pos} {s[:FLANK]}[{alleles[0]}/{','.join(alleles[1:])}]{s[FLANK + 1:]}", flush=True)
    time.sleep(0.15)

Path(sys.argv[1]).write_text(json.dumps({
    "source": "Ensembl GRCh37 REST (variation and sequence/region endpoints), forward strand",
    "flank": FLANK, "variants": out}, indent=1, sort_keys=True) + "\n")
print(f"{len(out)} variants; {len(fails)} failure(s)")
for f in fails:
    print("FAIL:", f)
sys.exit(1 if fails else 0)
