"""TfReg (Weber, UFMG) source package via the openSUSE build service mirror: which experimental
melting-temperature files it carries (probe)."""
import io
import re
import tarfile
import urllib.request

ROOT = "https://download.opensuse.org/repositories/home:/drgweber/"


def page(u):
    return urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "athanor-ci"}), timeout=90).read()


found = []
top = [l for l in re.findall(r'href="\./([^"]+/)"', page(ROOT).decode()) if l]
print("dists:", top)
for d in top:
    for sub in ("", "src/", "all/", "x86_64/"):
        u = ROOT + d + sub
        try:
            html = page(u).decode("utf-8", "replace")
        except Exception:  # noqa: BLE001
            continue
        for l in re.findall(r'href="([^"]+)"', html):
            if "tfreg" in l.lower() and not l.endswith("/"):
                found.append(u + l.lstrip("./"))
found = sorted(set(found))
print("\n".join(found[:80]))
src = next((f for f in found if f.endswith((".tar.gz", ".tar.xz", ".tar.bz2")) and "orig" in f), None) or \
    next((f for f in found if f.endswith((".tar.gz", ".tar.xz", ".tar.bz2"))), None)
print("source:", src)
if src:
    data = page(src)
    with tarfile.open(fileobj=io.BytesIO(data)) as t:
        names = t.getnames()
        print(len(names), "files")
        dats = [n for n in names if n.endswith(".dat") or "COPYING" in n or "LICENSE" in n.upper()]
        print("\n".join(dats))
        for n in dats:
            low = n.lower()
            if any(k in low for k in ("owczarzy", "santalucia", "allawi")):
                f = t.extractfile(n)
                txt = f.read().decode("utf-8", "replace").splitlines() if f else []
                print(f"--- {n}: {len(txt)} lines")
                print("\n".join(txt[:10]))
