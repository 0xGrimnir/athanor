"""Sources for DNA Structure (probe): 3DNA availability, PDB entries, NAKB annotations."""
import urllib.request


def get(url, n=3000):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "athanor-ci"})
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read()
        print(f"\n### {url}  ({len(body)} bytes, {r.headers.get('Content-Type')})")
        print(body[:n].decode("utf-8", "replace"))
    except Exception as e:  # noqa: BLE001
        print(f"\n### {url}\nERROR {e!r}")


get("https://aur.archlinux.org/cgit/aur.git/plain/PKGBUILD?h=x3dna", 4000)
get("https://x3dna.org/downloads", 1500)
get("http://forum.x3dna.org/downloads/3dna/", 1500)
get("https://files.rcsb.org/download/1BNA.cif", 800)
get("https://files.rcsb.org/download/355D.cif", 400)
get("https://data.rcsb.org/rest/v1/core/entry/355D", 1500)
get("https://www.nakb.org/", 1500)
get("https://www.nakb.org/api/", 1500)
get("https://www.nakb.org/dssr/355D.json", 1500)
get("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/355d", 800)
get("http://protein.bio.puc.cl/melting-temperatures.html", 4000)
get("http://protein.bio.puc.cl/", 1500)
