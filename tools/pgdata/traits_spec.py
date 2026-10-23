"""The variants behind the trait reports, and how each effect allele is pinned down.

Every effect allele is written on the GRCh37 + strand and checked in CI
against an independent source before the pack is built:
  aa      the amino acid the allele codes (Ensembl VEP), for coding variants
  expect  the allele's frequency range in 1000 Genomes superpopulations
          (Ensembl), which catches a flipped strand or a swapped allele
  gwas    keywords for the GWAS Catalog association whose effect allele and
          direction are recorded with the report (and must agree)
  direction  for an effect allele given here, whether the GWAS Catalog must
          record it as raising ("increase") or lowering ("decrease") the trait
  orient  for an effect allele left to the catalog: take the allele it
          records as lowering ("lower") the trait
  alt     the allele must be the non-reference allele in Ensembl
A check that fails stops the build.
"""

REFS = {
    "enattah2002": "11788828", "tishkoff2007": "17159977", "brooks2009": "19320537", "edenberg2007": "17718394",
    "sachse1999": "10233211", "cornelis2011": "21490707", "kim2003": "12595690", "yoshiura2006": "16444273",
    "north1999": "10192379", "yang2003": "12879365", "fujimoto2008": "18065779", "kamberov2013": "23415220",
    "medland2009": "19896111", "valverde1995": "7581459", "sturm2008": "18252222", "eiberg2008": "18172690",
    "lamason2005": "16357253", "han2008": "18483556", "sulem2007": "17952075", "eriksson2010": "20585627",
    "wang2010": "20541252", "hazra2008": "18776911", "tanaka2009": "19148276", "monsuur2008": "18509540",
    "praetorius2013": "24267888",
}

VARIANTS = {
    # appearance
    "rs12913832": {"gene": "HERC2", "effect": "G", "label": "HERC2 rs12913832 (blue-eye allele)",
                   "expect": {"EUR": (0.5, 0.85), "AFR": (0, 0.05), "EAS": (0, 0.05)}},
    "rs1800407": {"gene": "OCA2", "effect": "T", "label": "OCA2 R419Q", "aa": "Q"},
    "rs1805007": {"gene": "MC1R", "effect": "T", "label": "MC1R R151C", "aa": "C"},
    "rs1805008": {"gene": "MC1R", "effect": "T", "label": "MC1R R160W", "aa": "W"},
    "rs1805009": {"gene": "MC1R", "effect": "C", "label": "MC1R D294H", "aa": "H"},
    "rs1805006": {"gene": "MC1R", "effect": "A", "label": "MC1R D84E", "aa": "E"},
    "rs11547464": {"gene": "MC1R", "effect": "A", "label": "MC1R R142H", "aa": "H"},
    "rs12203592": {"gene": "IRF4", "effect": "T", "label": "IRF4 rs12203592", "gwas": ["freckl"], "direction": "increase"},
    "rs12821256": {"gene": "KITLG", "effect": "C", "label": "KITLG rs12821256 (blond-hair allele)", "gwas": ["blond"], "direction": "increase"},
    "rs1426654": {"gene": "SLC24A5", "effect": "A", "label": "SLC24A5 A111T", "aa": "T",
                  "expect": {"EUR": (0.95, 1.0), "EAS": (0, 0.05), "AFR": (0, 0.2)}},
    "rs16891982": {"gene": "SLC45A2", "effect": "G", "label": "SLC45A2 L374F", "aa": "F",
                   "expect": {"EUR": (0.85, 1.0), "AFR": (0, 0.1), "EAS": (0, 0.05)}},
    "rs17822931": {"gene": "ABCC11", "effect": "T", "label": "ABCC11 G180R (dry-earwax allele)", "aa": "R",
                   "expect": {"EAS": (0.6, 0.95), "AFR": (0, 0.05)}},
    "rs3827760": {"gene": "EDAR", "effect": "G", "label": "EDAR V370A", "aa": "A",
                  "expect": {"EAS": (0.7, 0.95), "EUR": (0, 0.05), "AFR": (0, 0.05)}},
    "rs11803731": {"gene": "TCHH", "effect": "T", "label": "TCHH L790M", "aa": "M",
                   "gwas": ["hair shape", "hair curvature", "hair morphology"], "direction": "decrease"},
    # food and drink
    "rs4988235": {"gene": "MCM6/LCT", "effect": "A", "label": "LCT -13910 C>T (European lactase persistence)",
                  "expect": {"EUR": (0.4, 0.8), "EAS": (0, 0.02)}},
    "rs41380347": {"gene": "MCM6/LCT", "effect": "C", "label": "LCT -13915 T>G (Arabian lactase persistence)", "alt": True},
    "rs145946881": {"gene": "MCM6/LCT", "effect": "G", "label": "LCT -14010 G>C (East African lactase persistence)", "alt": True},
    "rs41525747": {"gene": "MCM6/LCT", "effect": "C", "label": "LCT -13907 C>G (East African lactase persistence)", "alt": True},
    "rs869051967": {"gene": "MCM6/LCT", "effect": "C", "label": "LCT -14009 T>G (East African lactase persistence)", "alt": True},
    "rs671": {"gene": "ALDH2", "effect": "A", "label": "ALDH2 E504K (ALDH2*2)", "aa": "K",
              "expect": {"EAS": (0.1, 0.35), "EUR": (0, 0.01)}},
    "rs1229984": {"gene": "ADH1B", "effect": "T", "label": "ADH1B R48H (ADH1B*2)", "aa": "H",
                  "expect": {"EAS": (0.5, 0.9), "EUR": (0, 0.1)}},
    "rs762551": {"gene": "CYP1A2", "effect": "A", "label": "CYP1A2*1F (-163C>A)",
                 "expect": {"EUR": (0.55, 0.8)}},
    "rs2472297": {"gene": "CYP1A1/CYP1A2", "effect": "T", "label": "rs2472297", "gwas": ["caffeine level", "caffeine measurement", "caffeine metabolite"], "direction": "decrease"},
    "rs4410790": {"gene": "AHR", "effect": "C", "label": "rs4410790", "gwas": ["coffee"], "direction": "increase"},
    "rs713598": {"gene": "TAS2R38", "effect": "G", "label": "TAS2R38 A49P", "aa": "P"},
    "rs1726866": {"gene": "TAS2R38", "effect": None, "label": "TAS2R38 V262A", "aa": "A"},
    "rs10246939": {"gene": "TAS2R38", "effect": "C", "label": "TAS2R38 I296V", "aa": "V"},
    "rs2282679": {"gene": "GC", "effect": None, "label": "GC rs2282679", "gwas": ["vitamin d level", "25-hydroxyvitamin d"],
                  "orient": "lower"},
    "rs10741657": {"gene": "CYP2R1", "effect": None, "label": "CYP2R1 rs10741657",
                   "gwas": ["vitamin d level", "25-hydroxyvitamin d"], "orient": "lower"},
    "rs601338": {"gene": "FUT2", "effect": "A", "label": "FUT2 W143X (non-secretor)", "aa": "*",
                 "expect": {"EUR": (0.35, 0.6)}},
    "rs174537": {"gene": "FADS1", "effect": None, "label": "FADS1 rs174537",
                 "gwas": ["arachidonic acid", "arachidonoyl", "eicosapentaeno", "omega-6 polyunsaturated"], "orient": "lower"},
    "rs2187668": {"gene": "HLA-DQ", "effect": "T", "label": "HLA-DQ2.5 tag (rs2187668)",
                  "expect": {"EUR": (0.05, 0.2)}},
    "rs7454108": {"gene": "HLA-DQ", "effect": "C", "label": "HLA-DQ8 tag (rs7454108)", "alt": True},
    # fitness
    "rs1815739": {"gene": "ACTN3", "effect": "T", "label": "ACTN3 R577X", "expect": {"EUR": (0.35, 0.55), "AFR": (0, 0.2)}},
}
