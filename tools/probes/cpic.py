"""Shapes of the CPIC API responses and the other pharmacogenomic sources (probe)."""
import json
import urllib.request

API = "https://api.cpicpgx.org/v1/"


def get(url, n=3000):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=60) as r:
            body = r.read()
        print(f"\n### {url}  ({len(body)} bytes)")
        try:
            d = json.loads(body)
            if isinstance(d, list):
                print(f"rows: {len(d)}")
                print(json.dumps(d[:2], indent=1)[:n])
            else:
                print(json.dumps(d, indent=1)[:n])
            return d
        except ValueError:
            print(body[:n].decode("utf-8", "replace"))
    except Exception as e:  # noqa: BLE001 - a probe reports everything
        print(f"\n### {url}\nERROR {e!r}")


get(API + "gene?select=symbol,chr,lookupmethod,genesequenceid,chromosequenceid,mrnasequenceid,"
    "allelesLastModified,functionmethod&lookupmethod=not.is.null", 6000)
pairs = get(API + "pair?select=*&cpiclevel=eq.A&removed=eq.false", 2500)
if pairs:
    print("level A genes:", sorted({p.get("genesymbol") for p in pairs}))
get(API + "pair_view?limit=2")
get(API + "allele?genesymbol=eq.CYP2C19&limit=3")
get(API + "allele?genesymbol=eq.DPYD&limit=3")
get(API + "allele_definition?genesymbol=eq.CYP2C19&name=eq.*2&select=*,allele_location_value(*,sequence_location(*))", 5000)
get(API + "sequence_location?geneSymbol=eq.CYP2C19&limit=3")
get(API + "diplotype?genesymbol=eq.CYP2C9&diplotype=eq.*1/*3", 3000)
get(API + "gene_result?genesymbol=eq.CYP2C19&limit=2")
get(API + "gene_result_lookup?limit=2")
get(API + "drug?name=eq.clopidogrel")
get(API + "recommendation?drugid=eq.RxNorm:32968&limit=2", 4000)
get(API + "guideline?select=*,publication(title,pmid,year)&limit=2", 3000)
get(API + "test_alert?limit=1")
get(API + "allele_frequency?limit=2")
get(API + "change_log?limit=2")
# 1000 Genomes 30x phased panel and PharmCAT releases
get("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
    "20220422_3202_phased_SNV_INDEL_SV/", 4000)
get("https://api.github.com/repos/PharmGKB/PharmCAT/releases?per_page=1", 3000)
get("https://api.github.com/repos/cpicpgx/cpic-data/releases?per_page=1", 3000)
get("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/", 3000)
