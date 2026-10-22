"""Peek at the AADR MT repository FASTA (header format, sequence length)."""
import gzip, io, urllib.request, zlib
req = urllib.request.Request("https://dataverse.harvard.edu/api/access/datafile/13994529", headers={"User-Agent": "athanor-ci"})
data = urllib.request.urlopen(req, timeout=300).read()
txt = gzip.decompress(data).decode(errors="replace")
lines = txt.splitlines()
heads = [i for i, l in enumerate(lines) if l.startswith(">")]
print("records", len(heads))
for h in heads[:5]:
    seq = []
    j = h + 1
    while j < len(lines) and not lines[j].startswith(">"):
        seq.append(lines[j]); j += 1
    s = "".join(seq)
    print(repr(lines[h][:200]), len(s), s[:60], {c: s.count(c) for c in "ACGTN-"})
import collections
L = collections.Counter()
for k, h in enumerate(heads[:3000]):
    e = heads[k + 1] if k + 1 < len(heads) else len(lines)
    L[sum(len(x) for x in lines[h + 1:e])] += 1
print("length distribution (first 3000)", L.most_common(8))
