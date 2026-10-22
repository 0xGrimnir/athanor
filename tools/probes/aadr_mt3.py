"""The decoding logic of the AADR mtDNA uncompress script (the part after the
embedded reference sequence), plus a sample record."""
import re, urllib.request

def get(fid, rng=None):
    h = {"User-Agent": "athanor-ci"}
    if rng:
        h["Range"] = rng
    return urllib.request.urlopen(urllib.request.Request(
        "https://dataverse.harvard.edu/api/access/datafile/" + str(fid), headers=h), timeout=300).read()

src = get(13994857).decode(errors="replace")
print("script length", len(src))
# everything that is not the embedded reference block
body = re.sub(r'"""[^"]*"""', '<<DOCSTRING>>', src, flags=re.S)
body = re.sub(r'\n[ACGTNacgtn]{40,}', '', body)
print(body[-5000:])
