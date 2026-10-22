"""Inspect the AADR v66.p1 1240K annotation and .ind files: columns, present-day
groups, and candidate labels for the qpAdm/f4 populations."""
import collections, csv, io, re, urllib.request

BASE = "https://dataverse.harvard.edu/api/access/datafile/"
def get(fid):
    return urllib.request.urlopen(urllib.request.Request(BASE + str(fid), headers={"User-Agent": "athanor-ci"}), timeout=600).read()

anno = get(13994515).decode("utf-8", errors="replace")
rows = list(csv.reader(io.StringIO(anno), delimiter="\t"))
hdr = rows[0]
print("columns:")
for i, h in enumerate(hdr):
    print(i, repr(h[:140]), "| e.g.", repr(rows[1][i][:60]), "|", repr(rows[5000][i][:60]))
print("n rows", len(rows) - 1)
ind = get(13994513).decode()
print("ind head", ind.splitlines()[:3], "n", len(ind.splitlines()))
snp = urllib.request.urlopen(urllib.request.Request(BASE + "13994514", headers={"User-Agent": "athanor-ci", "Range": "bytes=0-2000"}), timeout=120).read().decode()
print("snp head", snp.splitlines()[:5])
gidx = next(i for i, h in enumerate(hdr) if h.lower().startswith("group id"))
didx = next(i for i, h in enumerate(hdr) if "date mean in bp" in h.lower())
groups = collections.Counter(r[gidx] for r in rows[1:] if len(r) > gidx)
modern = collections.Counter(r[gidx] for r in rows[1:] if len(r) > didx and r[didx].strip() in ("0", "0.0"))
print("groups", len(groups), "present-day groups", len(modern), "present-day individuals", sum(modern.values()))
print("present-day groups (top 300):", modern.most_common(300))
pats = ["Mbuti", "Yoruba", "Chimp", "Altai", "Vindija", "Chagyrskaya", "Denisova", "Ust", "Kostenki", "MA1|Malta", "Han", "Papuan",
        "Onge", "Chukchi", "Karitiana", "Loschbour", "Luxembourg", "WHG", "Turkey.*N|Barcin|Anatolia", "Yamnaya", "Ganj|Iran.*N",
        "Kotias|CHG|Satsurblia", "Natufian", "French", "English|GBR|CEU", "Sardinian", "Mota", "Mixe|Pima|Surui", "Sintashta"]
for p in pats:
    m = [(g, n) for g, n in groups.items() if re.search(p, g)]
    print(f"[{p}]", sorted(m, key=lambda x: -x[1])[:25])
