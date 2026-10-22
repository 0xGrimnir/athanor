"""First bytes of the AADR 1240K .geno (format check) and redirect/range support."""
import urllib.request
req = urllib.request.Request("https://dataverse.harvard.edu/api/access/datafile/13994829",
                             headers={"User-Agent": "athanor-ci", "Range": "bytes=0-199"})
r = urllib.request.urlopen(req, timeout=120)
print("status", r.status, "url", r.geturl()[:120])
print("headers", dict(r.headers))
b = r.read()
print(len(b), b[:120])
