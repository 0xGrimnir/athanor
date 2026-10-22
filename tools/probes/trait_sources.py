"""Show what the public sources for the trait reports return: Ensembl GRCh37
variation (positions, alleles, 1000 Genomes frequencies), Ensembl VEP (amino
acid per allele), the GWAS Catalog (old and v2 REST APIs) and NCBI
E-utilities, for a few example variants."""
import json
import urllib.request


def get(url, accept="application/json"):
    req = urllib.request.Request(url, headers={"User-Agent": "athanor-ci", "Accept": accept,
                                               "Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.status, r.read().decode("utf-8", errors="replace")
    except Exception as e:  # noqa: BLE001
        return None, repr(e)


def show(label, url, n=3000):
    st, body = get(url)
    print(f"===== {label}: {url} -> {st}")
    print(body[:n])
    print()


show("ensembl37 variation", "https://grch37.rest.ensembl.org/variation/human/rs671?pops=1", 2500)
st, body = get("https://grch37.rest.ensembl.org/variation/human/rs4988235?pops=1")
try:
    d = json.loads(body)
    print("rs4988235 mappings", d.get("mappings"))
    print("populations (1000GENOMES:phase_3 only):",
          [p for p in d.get("populations", []) if p.get("population", "").startswith("1000GENOMES:phase_3:")][:20])
except Exception as e:  # noqa: BLE001
    print("parse failed", e)
show("ensembl37 vep", "https://grch37.rest.ensembl.org/vep/human/id/rs671", 2500)
show("gwas old snp associations", "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs4481887/associations?projection=associationBySnp", 4000)
show("gwas v2 associations", "https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id=rs4481887", 4000)
show("gwas v2 snp", "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms/rs4481887", 2000)
show("eutils esummary", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=16444273,20585627&retmode=json", 2500)
