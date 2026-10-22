"""List the AADR Dataverse dataset (versions, files, licence) and the Reich lab
curated-release index, to choose the reference files for Personal Genome."""
import json
import urllib.request

def get(url, n=None):
    try:
        r = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "athanor-ci"}), timeout=120)
        data = r.read(n) if n else r.read()
        return data
    except Exception as e:  # noqa: BLE001
        return f"ERROR {e}".encode()

ds = get("https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/FFIDCW")
try:
    d = json.loads(ds)["data"]
    v = d["latestVersion"]
    print("version", v.get("versionNumber"), v.get("versionMinorNumber"), v.get("releaseTime"))
    print("license", json.dumps(v.get("license")), "terms", (v.get("termsOfUse") or "")[:2000])
    for f in v["files"]:
        df = f["dataFile"]
        print(df["id"], df.get("filename"), df.get("filesize"), df.get("md5") or df.get("checksum"), f.get("directoryLabel"))
except Exception as e:  # noqa: BLE001
    print("dataverse parse error", e, ds[:2000])
idx = get("https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/")
print("== reichdata index\n", idx[:6000].decode(errors="replace"))
