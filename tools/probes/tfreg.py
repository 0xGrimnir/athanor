"""TfReg (Weber, UFMG) source package: which experimental melting-temperature files it carries (probe)."""
import io
import re
import tarfile
import urllib.request

base = "https://bioinf.fisica.ufmg.br/software/tfreg-8.0/"
html = urllib.request.urlopen(urllib.request.Request(base, headers={"User-Agent": "athanor-ci"}), timeout=60).read().decode("utf-8", "replace")
links = re.findall(r'href="([^"]+)"', html)
print("links:", links)
tar = next((l for l in links if l.endswith((".tar.gz", ".tgz", ".tar.bz2", ".tar.xz"))), None)
if tar:
    url = tar if tar.startswith("http") else base + tar
    data = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "athanor-ci"}), timeout=180).read()
    print(url, len(data))
    with tarfile.open(fileobj=io.BytesIO(data)) as t:
        names = t.getnames()
        print(len(names), "files")
        dats = [n for n in names if n.endswith(".dat") or "COPYING" in n or "LICENSE" in n.upper()]
        print("\n".join(dats))
        for n in dats:
            if "owczarzy" in n.lower() or "santalucia" in n.lower() or "dna" in n.lower():
                f = t.extractfile(n)
                if f:
                    txt = f.read().decode("utf-8", "replace").splitlines()
                    print(f"--- {n}: {len(txt)} lines")
                    print("\n".join(txt[:12]))
