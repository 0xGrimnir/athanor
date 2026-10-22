"""Compact listing for the trait-report variants: GWAS Catalog v2 associations
(all fields of the first few), VEP amino acids per allele, 1000 Genomes
superpopulation frequencies, and PubMed titles for candidate references."""
import json
import urllib.request


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "athanor-ci", "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode("utf-8", errors="replace"))
    except Exception as e:  # noqa: BLE001
        return {"error": repr(e)}


GWAS = ["rs2282679", "rs10741657", "rs602662", "rs601338", "rs174537", "rs11803731", "rs12203592", "rs10427255",
        "rs72921001", "rs4481887", "rs12821256", "rs762551", "rs4410790", "rs2472297"]
for rs in GWAS:
    d = get(f"https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id={rs}&size=50")
    assoc = (d.get("_embedded") or {}).get("associations", [])
    print(f"== {rs}: {len(assoc)} associations {d.get('error', '')}")
    for a in assoc[:40]:
        a = {k: v for k, v in a.items() if k not in ("_links", "bg_efo_traits", "multi_snp_haplotype", "snp_interaction",
                                                     "pvalue_description", "pvalue_mantissa", "pvalue_exponent")}
        print("  ", json.dumps(a)[:700])

VEP = ["rs713598", "rs1726866", "rs10246939", "rs17822931", "rs3827760", "rs1805007", "rs1805008", "rs1805009",
       "rs1805006", "rs11547464", "rs1815739", "rs1800407", "rs1229984", "rs16891982", "rs1426654", "rs11803731"]
for rs in VEP:
    d = get(f"https://grch37.rest.ensembl.org/vep/human/id/{rs}")
    out = set()
    if isinstance(d, list):
        for v in d:
            for t in v.get("transcript_consequences", []):
                if t.get("amino_acids") and t.get("biotype") == "protein_coding":
                    out.add((t.get("gene_symbol"), t.get("variant_allele"), t.get("amino_acids"), t.get("protein_start"),
                             t.get("strand")))
    print(f"VEP {rs}: {sorted(out)[:8]} {d.get('error', '') if isinstance(d, dict) else ''}")

FREQ = ["rs12913832", "rs4988235", "rs41380347", "rs145946881", "rs41525747", "rs869051967", "rs2187668", "rs7454108",
        "rs1426654", "rs16891982", "rs762551", "rs671", "rs1229984", "rs17822931", "rs3827760", "rs1815739"]
for rs in FREQ:
    d = get(f"https://grch37.rest.ensembl.org/variation/human/{rs}?pops=1")
    if "error" in d:
        print(f"FREQ {rs}: {d['error']}")
        continue
    m = [(x["seq_region_name"], x["start"], x["allele_string"]) for x in d.get("mappings", [])]
    f = {}
    for p in d.get("populations", []):
        pop = p.get("population", "")
        if pop.startswith("1000GENOMES:phase_3:") and pop.split(":")[-1] in ("AFR", "AMR", "EAS", "EUR", "SAS"):
            f.setdefault(pop.split(":")[-1], {})[p["allele"]] = round(p["frequency"], 3)
    print(f"FREQ {rs}: {m} {f} minor={d.get('minor_allele')} maf={d.get('MAF')}")

PMIDS = {"11788828": "Enattah lactase", "17159977": "Tishkoff lactase Africa", "19048317": "Ingram lactase",
         "19320537": "Brooks flushing", "17718394": "Edenberg alcohol", "10233211": "Sachse CYP1A2",
         "12595690": "Kim PTC", "15723789": "Bufe TAS2R38", "16444273": "Yoshiura earwax", "19710691": "Martin axillary",
         "10192379": "North ACTN3", "12879365": "Yang ACTN3 athletes", "18065779": "Fujimoto EDAR",
         "23415220": "Kamberov EDAR", "19896111": "Medland hair", "7581459": "Valverde MC1R",
         "14986829": "Duffy MC1R", "18252222": "Sturm HERC2", "18172690": "Eiberg HERC2", "16357253": "Lamason SLC24A5",
         "18483556": "Han pigmentation", "17952075": "Sulem pigmentation", "20585627": "Eriksson 23andMe",
         "20541252": "Wang vitamin D", "18776911": "Hazra B12", "19148276": "Tanaka fatty acids",
         "18509540": "Monsuur celiac", "24267888": "Praetorius IRF4", "21490707": "Cornelis coffee"}
d = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(PMIDS))
for k, exp in PMIDS.items():
    r = (d.get("result") or {}).get(k, {})
    print(f"PMID {k} [{exp}]: {r.get('sortfirstauthor')} | {r.get('source')} {r.get('pubdate')} | {r.get('title')} | "
          f"{[x['value'] for x in r.get('articleids', []) if x.get('idtype') == 'doi']}")
