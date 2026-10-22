"""Decode the AADR mitogenome repository to FASTA.

The repository stores one record per line: an identifier, a tab, then one
character per position of the mitochondrial reference. 'Q' means "the same as
the reference here"; any other character is the base itself; a record shorter
than the reference is padded. This is the scheme in the AADR's own
mtdna_uncompress script (written for Python 2, so it is reimplemented here
rather than run), and the reference is read out of that script so the exact
sequence the authors used is the one applied.

Usage: python mt_decode.py MT_REPO_FA_GZ UNCOMPRESS_PY OUT_FA
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

repo, script, out = (Path(a) for a in sys.argv[1:4])

src = script.read_text(errors="replace")
m = re.search(r'rsrs\s*=\s*"""(.*?)"""', src, re.S)
if not m:
    raise SystemExit("could not find the reference sequence in the uncompress script")
ref = "".join(m.group(1).split())
if len(ref) < 16000:
    raise SystemExit(f"the reference sequence read from the script is only {len(ref)} long")

n = 0
with gzip.open(repo, "rt") as fh, open(out, "w") as fo:
    for line in fh:
        line = line.strip()
        if not line or "\t" not in line:
            continue
        ident, seq = line.split("\t", 1)
        ident = ident.lstrip(">")
        seq = seq + "-" * (len(ref) - len(seq))
        decoded = "".join(r if c == "Q" else c for c, r in zip(seq, ref))
        fo.write(f">{ident}\n" + "\n".join(decoded[i:i + 60] for i in range(0, len(decoded), 60)) + "\n")
        n += 1
print(f"decoded {n} mitogenomes of {len(ref)} positions each")
