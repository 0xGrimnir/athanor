"""All GWAS Catalog associations (trait, effect allele, beta, OR, author) for a few variants."""
import json, urllib.request
def get(u):
    with urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "athanor-ci", "Accept": "application/json"}), timeout=120) as r:
        return json.loads(r.read())
for rs in ["rs174537", "rs10427255", "rs2472297", "rs12203592", "rs12821256", "rs11803731", "rs4410790", "rs2282679", "rs10741657"]:
    page, rows = 0, []
    while True:
        d = get(f"https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id={rs}&size=100&page={page}")
        rows += (d.get("_embedded") or {}).get("associations", [])
        if page + 1 >= (d.get("page") or {}).get("totalPages", 1):
            break
        page += 1
    print(f"== {rs}: {len(rows)}")
    for a in rows:
        tr = (a.get("reported_trait") or [""])[0]
        efo = ",".join(e.get("efo_trait", "") for e in a.get("efo_traits", []))
        print(f"   {a.get('snp_effect_allele')} beta={a.get('beta')} OR={a.get('or_per_copy_num')} p={a.get('p_value')} | {tr[:70]} | {efo[:60]} | {a.get('first_author')} {a.get('pubmed_id')}")
