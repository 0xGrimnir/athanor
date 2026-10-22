"""The AADR's mtDNA uncompress script, and candidate labels for the
Western hunter-gatherer source."""
import csv, io, re, urllib.request

def get(fid):
    return urllib.request.urlopen(urllib.request.Request(
        "https://dataverse.harvard.edu/api/access/datafile/" + str(fid),
        headers={"User-Agent": "athanor-ci"}), timeout=300).read()

print("== mtdna_uncompress script")
print(get(13994857).decode(errors="replace")[:6000])

print("\n== hunter-gatherer candidate groups (western Europe, pre-farming)")
anno = get(13994515).decode("utf-8", errors="replace")
rows = list(csv.reader(io.StringIO(anno), delimiter="\t"))
hdr = rows[0]
gi = next(i for i, h in enumerate(hdr) if h.lower().startswith("group id"))
di = next(i for i, h in enumerate(hdr) if h.lower().startswith("date mean in bp"))
si = next(i for i, h in enumerate(hdr) if "1240k snpset" in h.lower())
ai = next(i for i, h in enumerate(hdr) if h.lower() == "assessment")
import collections
c = collections.Counter()
best = {}
for r in rows[1:]:
    if len(r) <= max(gi, di, si, ai):
        continue
    g = r[gi]
    try:
        bp = float(r[di]); snps = float(r[si] or 0)
    except ValueError:
        continue
    if not r[ai].upper().startswith("PASS") or snps < 30000:
        continue
    if re.search(r"Mesolithic|_HG|HG_|Epipaleolithic|LateUP|UP$|Magdalenian|Azilian", g) and bp > 6000:
        c[g] += 1
        best[g] = max(best.get(g, 0), snps)
for g, n in sorted(c.items(), key=lambda x: -x[1])[:60]:
    print(f"{n:4d}  {g}   (best {int(best[g])} SNPs)")
