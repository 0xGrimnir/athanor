"""Pharmacogenomic caller: CPIC table semantics, the unphased matcher, the
variant-list genes, chip input, recommendations and the 1000 Genomes check.

The gene tables here are small made-up examples shaped like CPIC's."""
import numpy as np

from athanor_core.personal.pgx import (NOVEL, UNTESTED, Gene, Panel, PgxGenotyper, pgx_cards,
                                       recommendations)
from athanor_core.personal.rawfile import CODE, RawGenotypes


def pos(rsid, p, alleles, kind="snv", chip=None, **kw):
    return {"rsid": rsid, "pos38": p + 1000, "pos37": p, "hgvs": f"g.{p}", "name": None, "kind": kind,
            "alleles": alleles, "genome": 0, "chip": chip, **kw}


def lookups(pairs, lookup="PHENOTYPE"):
    return [{"key": k, "result": r, "score": s, "description": None,
             "priority": "Normal/Routine/Low Risk" if r.startswith("Normal") else "Abnormal/Priority/High Risk"}
            for k, r, s in pairs]


# a TPMT-like gene: *3A = *3B + *3C in cis; *9 defined with an IUPAC-style set
TPMT = {"chrom": 6, "chrom_name": "chr6", "lookup": "PHENOTYPE",
        "positions": [pos("rs1", 100, ["C", "T"]), pos("rs2", 200, ["T", "C"]), pos("rs3", 300, ["A", "G", "T"])],
        "alleles": [
            {"name": "*1", "ref": True, "def": [[0], [0], [0]], "function": "Normal function", "freq": {"European": 0.9}},
            {"name": "*3A", "def": [[1], [1], None], "function": "No function", "freq": {"European": 0.04}},
            {"name": "*3B", "def": [[1], None, None], "function": "No function", "freq": {"European": 0.001}},
            {"name": "*3C", "def": [None, [1], None], "function": "No function", "freq": {"European": 0.004}},
            {"name": "*9", "def": [None, None, [1, 2]], "function": "Decreased function", "freq": {"European": 0.02}},
        ],
        "results": lookups([({"Normal function": 2}, "Normal Metabolizer", None),
                            ({"Normal function": 1, "No function": 1}, "Intermediate Metabolizer", None),
                            ({"Normal function": 1, "Decreased function": 1}, "Intermediate Metabolizer", None),
                            ({"No function": 2}, "Poor Metabolizer", None),
                            ({"No function": 1, "Decreased function": 1}, "Poor Metabolizer", None),
                            ({"Decreased function": 2}, "Intermediate Metabolizer", None)]),
        "phenotypes": {}}

DPYD = {"chrom": 1, "chrom_name": "chr1", "lookup": "ACTIVITY_SCORE",
        "positions": [pos("rs3918290", 1000, ["C", "T"]), pos("rs56038477", 2000, ["C", "T"]),
                      pos("rs75017182", 2500, ["G", "C"], primary=True), pos("rs67376798", 3000, ["T", "A"])],
        "alleles": [
            {"name": "Reference", "ref": True, "def": [[0]] * 4, "activity": "1.0", "function": "Normal function"},
            {"name": "c.1905+1G>A (*2A)", "def": [[1], None, None, None], "activity": "0.0", "function": "No function"},
            {"name": "c.1129-5923C>G", "def": [None, None, [1], None], "activity": "0.5", "function": "Decreased function"},
            {"name": "HapB3", "def": [None, [1], [1], None], "activity": "0.5", "function": "Decreased function"},
            {"name": "c.2846A>T", "def": [None, None, None, [1]], "activity": "0.5", "function": "Decreased function"},
        ],
        "results": lookups([({"1.0": 2}, "Normal Metabolizer", "2.0"), ({"1.0": 1, "0.5": 1}, "Intermediate Metabolizer", "1.5"),
                            ({"1.0": 1, "0.0": 1}, "Intermediate Metabolizer", "1.0"), ({"0.5": 2}, "Intermediate Metabolizer", "1.0"),
                            ({"0.5": 1, "0.0": 1}, "Poor Metabolizer", "0.5"), ({"0.0": 2}, "Poor Metabolizer", "0.0")]),
        "phenotypes": {}}


def G(gene, calls, default=0):
    out = np.full((len(gene.pos), 2), default, dtype=np.uint8)
    for k, p in enumerate(gene.pos):
        if p["rsid"] in calls:
            v = calls[p["rsid"]]
            out[k] = UNTESTED if v is None else [p["alleles"].index(x) for x in v]
    return out


def test_reference_het_hom_and_longest_allele():
    g = Gene("TPMT", TPMT)
    assert g.call(G(g, {}))["diplotypes"] == ["*1/*1"]
    c = g.call(G(g, {"rs1": "CT"}))
    assert c["diplotypes"] == ["*1/*3B"] and c["phenotype"] == "Intermediate Metabolizer"
    assert g.call(G(g, {"rs1": "TT", "rs2": "CC"}))["diplotypes"] == ["*3A/*3A"]
    # unphased double heterozygote: *1/*3A scores higher than *3B/*3C (as PharmCAT), and the
    # lower-scoring alternative with a different phenotype is reported
    c = g.call(G(g, {"rs1": "CT", "rs2": "TC"}), alternatives=True)
    assert c["diplotypes"] == ["*1/*3A"] and c["phenotype"] == "Intermediate Metabolizer"
    assert {"diplotype": "*3B/*3C", "phenotype": "Poor Metabolizer"} in c["alternatives"]


def test_sets_in_definitions_and_missing_positions():
    g = Gene("TPMT", TPMT)
    assert g.call(G(g, {"rs3": "AG"}))["diplotypes"] == ["*1/*9"]
    assert g.call(G(g, {"rs3": "AT"}))["diplotypes"] == ["*1/*9"]     # either listed base defines *9
    # rs1 not read: *3A and *3C both fit a C/T at rs2 with the same score
    c = g.call(G(g, {"rs1": None, "rs2": "TC"}))
    assert set(c["diplotypes"]) == {"*1/*3A", "*1/*3C"} and c["phenotype"] == "Intermediate Metabolizer"
    assert g.call(G(g, {}, default=UNTESTED))["status"] == "not tested"


def test_unlisted_base():
    g = Gene("CYPX", TPMT)
    geno = G(g, {})
    geno[0] = [0, NOVEL]
    assert g.call(geno)["status"].startswith("no call")
    t = Gene("TPMT", TPMT)                    # a toxicity gene: read as the reference
    c = t.call(geno)
    assert c["diplotypes"] == ["*1/*1"] and c["notes"]


def test_haploid():
    g = Gene("TPMT", TPMT)
    c = g.call(G(g, {"rs1": "TT", "rs2": "CC"}), haploid=True)
    assert c["diplotypes"] == ["*3A"]


def test_dpyd_two_lowest_and_hapb3():
    g = Gene("DPYD", DPYD)
    c = g.call(G(g, {}))
    assert c["diplotypes"] == ["Reference/Reference"] and c["score"] == "2.0"
    c = g.call(G(g, {"rs56038477": "CT", "rs75017182": "GC"}))
    assert [v["name"] for v in c["variants"]] == ["HapB3"] and c["score"] == "1.5"
    c = g.call(G(g, {"rs3918290": "CT", "rs56038477": "CT", "rs75017182": "GC", "rs67376798": "TA"}))
    assert c["phenotype"] == "Poor Metabolizer" and c["score"] == "0.5"      # *2A (0) + one 0.5 variant
    c = g.call(G(g, {"rs67376798": "AA"}))
    assert c["score"] == "1.0"                                                 # homozygous counts twice
    # intronic SNP missing: HapB3 from the exonic SNP, with a note
    c = g.call(G(g, {"rs56038477": "CT", "rs75017182": None}))
    assert c["score"] == "1.5" and any("rs56038477" in n for n in c["notes"])
    # the intronic SNP decides when present
    c = g.call(G(g, {"rs56038477": "CT", "rs75017182": "GG"}))
    assert c["score"] == "2.0"


def raw(rows, vendor="23andMe"):
    rs = [int(r[0][2:]) for r in rows]
    ch = [r[1] for r in rows]
    ps = [r[2] for r in rows]
    a1 = [CODE.get(r[3][0], 255) for r in rows]
    a2 = [CODE.get(r[3][-1], 255) for r in rows]
    return RawGenotypes(vendor=vendor, file_name="t.txt", sha256="", rsid=np.array(rs), chrom=np.array(ch, dtype=np.int8),
                        pos=np.array(ps), a1=np.array(a1, dtype=np.uint8), a2=np.array(a2, dtype=np.uint8))


def test_chip_input_strand_and_indels():
    d = {**TPMT, "positions": TPMT["positions"] + [pos("rs9", 400, ["T", "delT"], kind="indel", chip={"D": 1, "I": 0})]}
    d["alleles"] = [dict(a, **{"def": a["def"] + ([[0]] if a.get("ref") else [None])}) for a in TPMT["alleles"]]
    d["alleles"].append({"name": "*X", "def": [None, None, None, [1]], "function": "No function"})
    g = Gene("TPMT", d)
    r = raw([("rs1", 6, 100, "CT"), ("rs2", 6, 200, "TT"), ("rs9", 6, 400, "DI")])
    gt = PgxGenotyper(r, r.pos, "the reference's + strand")
    geno = gt.genotype(g, 6)
    assert geno[0].tolist() == [0, 1] and geno[3].tolist() == [1, 0] and geno[2].tolist() == [UNTESTED, UNTESTED]
    assert g.call(geno)["diplotypes"] == ["*3B/*X"]
    # the same calls on the opposite strand
    r2 = raw([("rs1", 6, 100, "GA"), ("rs2", 6, 200, "AA")])
    g2 = PgxGenotyper(r2, r2.pos, "the opposite strand").genotype(g, 6)
    assert g2[0].tolist() == [0, 1] and g2[1].tolist() == [0, 0]
    # a chip row for an indel at the same coordinate as a SNV does not hide the SNV
    r4 = raw([("rs77", 6, 100, "DI"), ("rs1", 6, 100, "TT")])
    assert PgxGenotyper(r4, r4.pos, "the reference's + strand").genotype(g, 6)[0].tolist() == [1, 1]
    # found by rsID when the position differs (an indel placed elsewhere by the chip)
    r3 = raw([("rs9", 6, 399, "DD")])
    assert PgxGenotyper(r3, r3.pos, "the reference's + strand").genotype(g, 6)[3].tolist() == [1, 1]


def test_recommendations_use_no_result_for_unknown_genes():
    pgx = {"drugs": [{"name": "drugA", "genes": ["TPMT", "CYP2D6"], "guideline": {},
                      "recs": [{"key": {"TPMT": "Poor Metabolizer", "CYP2D6": "No Result"}, "recommendation": "reduce"},
                               {"key": {"TPMT": "Poor Metabolizer", "CYP2D6": "1.0"}, "recommendation": "reduce more"},
                               {"key": {"TPMT": "Normal Metabolizer", "CYP2D6": "No Result"}, "recommendation": "standard"}]},
                     {"name": "drugB", "genes": ["DPYD"], "guideline": {},
                      "recs": [{"key": {"DPYD": "1.5"}, "recommendation": "half"}, {"key": {"DPYD": "2.0"}, "recommendation": "full"}]}]}
    out = recommendations(pgx, {"TPMT": "Poor Metabolizer", "DPYD": "1.50"}, {})
    a, b = out
    assert [r["recommendation"] for r in a["rows"]] == ["reduce"] and a["unknown"] == ["CYP2D6"]
    assert [r["recommendation"] for r in b["rows"]] == ["half"]


def panel_for(gene_def, n=400, seed=1):
    """Haplotypes drawn from the definitions: *1 80%, *3A 10%, *3C 6%, *9 4%."""
    rng = np.random.default_rng(seed)
    names = ["*1", "*3A", "*3C", "*9"]
    pick = rng.choice(4, size=2 * n, p=[0.8, 0.1, 0.06, 0.04])
    al = {a["name"]: a for a in gene_def["alleles"]}
    H = np.zeros((2 * n, len(gene_def["positions"])), dtype=np.uint8)
    for h, i in enumerate(pick):
        for k, s in enumerate(al[names[i]]["def"]):
            H[h, k] = 0 if s is None else s[0]
    return {"samples": np.array([f"S{i}" for i in range(n)]), "superpop": np.array(["EUR", "AFR"] * (n // 2)),
            "male": np.zeros(n, dtype=bool), "related": np.zeros(n, dtype=bool), "hap_TPMT": H}


def test_panel_agreement_drops_when_positions_are_missing():
    g = Gene("TPMT", TPMT)
    panel = Panel(panel_for(TPMT), {"TPMT": TPMT})
    allpos = np.ones(3, dtype=bool)
    c = g.call(G(g, {}))
    full = panel.agreement(g, 6, allpos, c)
    assert full["phenotype_agrees"] == 1.0 and full["positions_used"] == 3
    # without rs2 a *3C carrier looks like *1: the reference call is right less often
    part = np.array([True, False, True])
    geno = G(g, {"rs2": None})
    red = panel.agreement(g, 6, part, g.call(geno))
    assert red["same_call"] > 0 and red["phenotype_agrees"] < 1.0
    assert abs(red["phenotype_agrees"] - 0.8 ** 2 / (1 - 0.1 - 0.04) ** 2 * 1) < 0.15


def test_cards_end_to_end():
    pgx = {"genes": {"TPMT": TPMT, "DPYD": DPYD}, "references": {},
           "drugs": [{"name": "mercaptopurine", "genes": ["TPMT"], "guideline": {"name": "thiopurines"},
                      "recs": [{"key": {"TPMT": "Intermediate Metabolizer"}, "recommendation": "start at 30-80%"}]}]}
    r = raw([("rs1", 6, 100, "CT"), ("rs2", 6, 200, "TC"), ("rs3", 6, 300, "AA"), ("rs3918290", 1, 1000, "CC")])

    def make(cid, title, group, grade, summary, result, limits, method, sources, ok=True):
        return {"id": cid, "title": title, "grade": grade, "summary": summary, "result": result, "ok": ok}
    cards = pgx_cards(pgx, Panel(panel_for(TPMT), pgx["genes"]), r, r.pos, "the reference's + strand", "XX", make)
    ids = [c["id"] for c in cards]
    assert ids[0] == "medicines" and "pgx_tpmt" in ids and "pgx_dpyd" in ids
    tp = next(c for c in cards if c["id"] == "pgx_tpmt")
    assert tp["result"]["call"]["diplotypes"] == ["*1/*3A"] and "1000 Genomes" in tp["summary"]
    s = cards[0]
    assert s["result"]["flagged"] == ["TPMT"]
    assert s["result"]["drugs"][0]["rows"][0]["recommendation"] == "start at 30-80%"


def test_frequency_choice_for_unphased_nat2_like_genes():
    # a double heterozygote fits *1/*5 (*5 carries both variants) and *4/*16 (one each)
    d = {"chrom": 8, "chrom_name": "chr8", "lookup": "PHENOTYPE",
         "positions": [pos("rs1", 10, ["T", "C"]), pos("rs2", 20, ["G", "A"])],
         "alleles": [{"name": "*1", "ref": True, "def": [[0], [0]], "function": "Normal function"},
                     {"name": "*4", "def": [None, [1]], "function": "Normal function"},
                     {"name": "*5", "def": [[1], [1]], "function": "Decreased function"},
                     {"name": "*16", "def": [[1], None], "function": "Decreased function"}],
         "results": lookups([({"Normal function": 2}, "Rapid Metabolizer", None),
                             ({"Normal function": 1, "Decreased function": 1}, "Intermediate Metabolizer", None),
                             ({"Decreased function": 2}, "Poor Metabolizer", None)]),
         "phenotypes": {}, "hap_freq": {"all": {"*1": 0.05, "*4": 0.5, "*5": 0.4, "*16": 0.05}}}
    g = Gene("NAT2", d)
    c = g.call(G(g, {"rs1": "TC", "rs2": "GA"}))
    # by score *1/*5 would win; by frequency *4/*16 (2 x 0.5 x 0.05 = 0.05) beats *1/*5 (2 x 0.05 x 0.4 = 0.04)
    assert c["chosen_by"] == "haplotype frequency" and c["diplotypes"] == ["*4/*16"]
    assert abs(c["diplotype_probability"] - 0.025 / 0.045) < 1e-9
    assert set(c["phenotype_probabilities"]) == {"Intermediate Metabolizer"}


def test_unread_common_allele_makes_the_call_unknown():
    d = {**TPMT, "alleles": [dict(a, freq={"European": 0.06}) if a["name"] == "*3C" else a for a in TPMT["alleles"]]}
    pgx = {"genes": {"TPMT": d}, "references": {},
           "drugs": [{"name": "mercaptopurine", "genes": ["TPMT"], "guideline": {},
                      "recs": [{"key": {"TPMT": "Normal Metabolizer"}, "recommendation": "standard"}]}]}
    r = raw([("rs1", 6, 100, "CC"), ("rs3", 6, 300, "AA")])           # rs2 (*3C) not on the chip

    def make(cid, title, group, grade, summary, result, limits, method, sources, ok=True):
        return {"id": cid, "grade": grade, "summary": summary, "result": result}
    cards = pgx_cards(pgx, None, r, r.pos, "the reference's + strand", "XX", make)
    tp = next(c for c in cards if c["id"] == "pgx_tpmt")
    assert tp["result"]["call"]["status"] == "unreliable" and tp["grade"] == "C"
    assert "*3C" in tp["summary"] and "unknown" in tp["summary"]
    # an unknown gene selects no recommendation row
    assert cards[0]["result"]["drugs"][0]["depends"] and not cards[0]["result"]["drugs"][0]["rows"]
