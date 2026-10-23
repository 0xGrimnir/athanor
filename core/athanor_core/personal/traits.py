"""Trait reports: appearance, food and drink, fitness.

Each report reads a handful of well-studied variants. Effect alleles,
positions and population frequencies come from the reference pack's
traits.json, which CI builds from Ensembl (GRCh37, VEP, 1000 Genomes) and the
GWAS Catalog, checking every allele against an independent source
(tools/pgdata/build_traits.py). The wording of each report states only what
the cited studies support, and every report says what it cannot tell you.
"""
from __future__ import annotations

import numpy as np

from .rawfile import BASES, MISSING, RawGenotypes

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


class Genotyper:
    """Your genotype at any GRCh37 position, on the + strand."""

    def __init__(self, g: RawGenotypes, pos37: np.ndarray, strand: str):
        self.g = g
        self.flip = {"the reference's + strand": False, "the opposite strand": True}.get(strand)
        called = g.a1 != MISSING
        self.by_pos = {}
        self.by_rs = {}
        for i in np.flatnonzero(called & (pos37 > 0)):
            self.by_pos.setdefault((int(g.chrom[i]), int(pos37[i])), int(i))
        for i in np.flatnonzero(called & (g.rsid > 0)):
            self.by_rs.setdefault(int(g.rsid[i]), int(i))

    def call(self, rec: dict) -> tuple[int | None, str | None, str]:
        """(copies of the effect allele, genotype on the + strand, note)."""
        i = self.by_pos.get((rec["chrom"], rec["pos37"]))
        if i is None:
            i = self.by_rs.get(int(rec["rsid"][2:]))
        if i is None:
            return None, None, "not on your chip"
        b1, b2 = self.g.a1[i], self.g.a2[i]
        if b1 > 3 or b2 > 3:
            return None, None, "not a base call"
        x1, x2 = BASES[b1], BASES[b2]
        allowed = {rec["ref"], *rec["alts"]}
        ambiguous = COMP.get(rec["effect"]) == rec["other"]
        if self.flip is True:
            x1, x2 = COMP[x1], COMP[x2]
        elif self.flip is None:
            if ambiguous:
                return None, None, "strand of the file undetermined at an A/T or C/G variant"
            if not {x1, x2} <= allowed and {COMP[x1], COMP[x2]} <= allowed:
                x1, x2 = COMP[x1], COMP[x2]
        if not {x1, x2} <= allowed:
            return None, f"{x1}{x2}", "alleles do not match the reference"
        return int(x1 == rec["effect"]) + int(x2 == rec["effect"]), "".join(sorted((x1, x2))), ""


def _ref(traits: dict, keys: list[str]) -> list[dict]:
    out = []
    for k in keys:
        r = traits["references"].get(k)
        if r:
            out.append({"title": f"{r['first_author']} et al., {r['title']}, {r['journal']} ({r['year']})",
                        "url": r["url"]})
    return out


def _gwas_refs(recs: list[dict]) -> list[dict]:
    out, seen = [], set()
    for rec in recs:
        for h in (rec.get("gwas") or [])[:1]:
            if h.get("pubmed_id") and h["pubmed_id"] not in seen:
                seen.add(h["pubmed_id"])
                out.append({"title": f"{h['first_author']} et al., GWAS Catalog {h['accession']}: {h['trait']}",
                            "url": f"https://www.ebi.ac.uk/gwas/studies/{h['accession']}"})
    return out


def _row(rec, copies, gt, note):
    return {"rsid": rec["rsid"], "gene": rec["gene"], "variant": rec["label"], "genotype": gt,
            "effect_allele": rec["effect"], "copies": copies, "note": note, "freq": rec.get("freq"),
            "gwas": (rec.get("gwas") or [None])[0]}


def _plural(n, one, many):
    return one if n == 1 else many


def trait_cards(traits: dict, gt: Genotyper, make_card) -> list[dict]:
    """make_card(cid, title, group, grade, summary, result, limits, method, sources, ok) -> card"""
    V = traits["variants"]
    cards = []

    def look(*rsids):
        rows, cnt = [], {}
        for rs in rsids:
            rec = V.get(rs)
            if not rec:
                continue
            c, g_, note = gt.call(rec)
            rows.append(_row(rec, c, g_, note))
            cnt[rs] = c
        return rows, cnt

    def add(cid, title, group, grade, summary, rows, limits, method, refs, extra=None, ok=True):
        res = {"variants": rows, **(extra or {})}
        cards.append(make_card(cid, title, group, grade, summary, res, limits, method, refs, ok))

    std_method = ("Your genotype at each listed variant, read on the GRCh37 + strand; effect alleles and 1000 Genomes "
                  "frequencies from Ensembl, checked in CI against the amino acid each allele codes or its known "
                  "population frequencies.")

    # ---------------------------------------------------------------- appearance
    rows, c = look("rs12913832", "rs1800407")
    n = c.get("rs12913832")
    if n is None:
        s, ok = "The key eye-colour variant (HERC2 rs12913832) is not on your chip.", False
    else:
        ok = True
        s = {2: "Two copies of the blue-eye variant in HERC2: most people with this genotype have blue or other light eyes (grey, light green); brown is uncommon.",
             1: "One copy of the blue-eye variant in HERC2: most people with this genotype have brown or hazel eyes; green is possible and blue uncommon.",
             0: "No copies of the blue-eye variant in HERC2: brown eyes are by far the most likely."}[n]
        if n >= 1 and (c.get("rs1800407") or 0) >= 1:
            s += " You also carry OCA2 R419Q, which shifts the odds towards green or hazel."
    add("eye_colour", "Eye colour", "Appearance", "A", s, rows,
        ["Blue versus brown is predicted well; green, hazel and other in-between colours much less so, because many other genes adjust the shade.",
         "The studies are mostly of people with European ancestry; elsewhere brown eyes are nearly universal whatever the genotype."],
        "A single variant in an enhancer within HERC2 controls how much OCA2 is expressed in the iris (Sturm et al. 2008; Eiberg et al. 2008); OCA2 R419Q modifies the shade. " + std_method,
        _ref(traits, ["sturm2008", "eiberg2008", "sulem2007"]), ok=ok)

    r_ids = ["rs1805006", "rs11547464", "rs1805007", "rs1805008", "rs1805009"]
    rows, c = look(*r_ids)
    typed = [r for r in r_ids if c.get(r) is not None]
    R = sum(c[r] for r in typed)
    if not typed:
        s, ok = "None of the main MC1R red-hair variants is on your chip.", False
    else:
        ok = True
        s = ("Two red-hair variants in MC1R: most people with two such variants have red or reddish hair and fair skin that burns easily." if R >= 2 else
             "One red-hair variant in MC1R (a carrier): red hair is uncommon with one copy, but freckling and sunburn are more likely, and a child can have red hair if the other parent is also a carrier." if R == 1 else
             f"None of the {len(typed)} main MC1R red-hair variants on your chip.")
    add("red_hair", "Red hair (MC1R)", "Appearance", "A", s, rows,
        ["MC1R has many rarer variants that chips do not test, so 'none found' does not rule out red hair or carrier status.",
         "Chip data do not show whether two variants sit on the same copy of the gene; in practice two red-hair variants are almost always on different copies.",
         f"Variants on your chip: {len(typed)} of {len(r_ids)}."],
        "Count of the MC1R variants with strong effects on red hair (D84E, R142H, R151C, R160W, D294H), which reduce the receptor's signalling so pigment cells make red-yellow pheomelanin instead of dark eumelanin (Valverde et al. 1995). " + std_method,
        _ref(traits, ["valverde1995", "sulem2007"]), {"red_hair_variants": R}, ok=ok)

    rows, c = look("rs12203592")
    t = c.get("rs12203592")
    parts = []
    if t is not None:
        parts.append(f"{t} {_plural(t, 'copy', 'copies')} of the IRF4 allele linked to more freckling and sun sensitivity")
    if typed:
        parts.append(f"{R} MC1R red-hair {_plural(R, 'variant', 'variants')}")
    level = (t or 0) + R
    s = ("; ".join(parts) + ". " + ("More freckle-prone than most people." if level >= 2 else
                                    "Somewhat more freckle-prone than average." if level == 1 else
                                    "No increase in freckling expected from these variants.")) if parts else "The freckling variants are not on your chip."
    add("freckles", "Freckles", "Appearance", "B", s, rows + [r for r in look(*r_ids)[0]],
        ["Freckling depends heavily on sun exposure, and on other genes not counted here."],
        "IRF4 rs12203592 lowers IRF4 expression in pigment cells, affecting tanning and freckling (Praetorius et al. 2013); MC1R red-hair variants add to freckling. The effect allele's association is checked against the GWAS Catalog. " + std_method,
        _ref(traits, ["praetorius2013", "han2008"]) + _gwas_refs([V.get("rs12203592", {})]), ok=bool(parts))

    rows, c = look("rs12821256")
    k = c.get("rs12821256")
    rec = V.get("rs12821256", {})
    orv = ((rec.get("gwas") or [{}])[0]).get("or")
    s = ("Not on your chip." if k is None else
         f"{k} {_plural(k, 'copy', 'copies')} of the KITLG variant that raises the odds of blond rather than brown hair in Europeans"
         + (f" (odds ratio about {float(orv):.1f} per copy)." if orv not in (None, '-') else "."))
    add("blond_hair", "Blond hair variant (KITLG)", "Appearance", "C", s, rows,
        ["One of many variants that affect hair lightness; on its own it does not predict hair colour.",
         "Studied in people of European ancestry; childhood blond hair often darkens with age."],
        "A regulatory variant upstream of KITLG associated with blond hair (Sulem et al. 2007); effect size from the GWAS Catalog. " + std_method,
        _ref(traits, ["sulem2007"]), ok=k is not None)

    rows, c = look("rs1426654", "rs16891982")
    a, b = c.get("rs1426654"), c.get("rs16891982")
    if a is None and b is None:
        s, ok = "Neither variant is on your chip.", False
    else:
        ok = True
        bits = []
        if a is not None:
            bits.append(f"SLC24A5: {a} of 2 copies of the light-skin allele")
        if b is not None:
            bits.append(f"SLC45A2: {b} of 2")
        s = "; ".join(bits) + ". These two variants have the largest known single effects on the lighter skin common in Europe; the SLC24A5 allele is also common in the Middle East, North Africa and South Asia."
    add("skin_pigmentation", "Skin pigmentation variants", "Appearance", "B", s, rows,
        ["Skin colour is shaped by many genes, and these two do not predict a shade.",
         "The lighter skin common in East Asia arose mostly through different variants, which these do not measure."],
        "SLC24A5 A111T (Lamason et al. 2005) and SLC45A2 L374F (Sulem et al. 2007), both of which reduce melanin production. " + std_method,
        _ref(traits, ["lamason2005", "sulem2007"]), ok=ok)

    rows, c = look("rs17822931")
    k = c.get("rs17822931")
    s = ("Not on your chip." if k is None else
         "Two copies of the ABCC11 dry-earwax allele: dry, flaky earwax." if k == 2 else
         "Wet earwax (the dry type needs two copies of the ABCC11 allele" + ("; you carry one)." if k == 1 else ")."))
    add("earwax", "Earwax type", "Appearance", "A", s, rows,
        ["A rare variant elsewhere in ABCC11 can also give dry earwax; chips usually do not test it."],
        "ABCC11 538G>A (G180R): the Arg form does not reach the cell membrane of the glands, giving dry earwax when both copies carry it (Yoshiura et al. 2006). " + std_method,
        _ref(traits, ["yoshiura2006"]), ok=k is not None)

    rows, c = look("rs3827760")
    k = c.get("rs3827760")
    s = ("Not on your chip." if k is None else
         f"{k} {_plural(k, 'copy', 'copies')} of EDAR 370A, the variant common in East Asia and the Americas that is linked to thicker, straighter hair."
         if k else "No copies of EDAR 370A; this variant is rare outside East Asia and the Americas.")
    add("hair_thickness", "Hair thickness (EDAR)", "Appearance", "B", s, rows,
        ["Hair thickness also depends on other genes; this variant explains part of the difference between East Asian and other hair."],
        "EDAR V370A is associated with thicker hair shafts (Fujimoto et al. 2008); mice carrying it grow thicker hair and more sweat glands (Kamberov et al. 2013). " + std_method,
        _ref(traits, ["fujimoto2008", "kamberov2013"]), ok=k is not None)

    rows, c = look("rs11803731")
    k = c.get("rs11803731")
    s = ("Not on your chip." if k is None else
         f"{k} {_plural(k, 'copy', 'copies')} of the TCHH variant associated with straighter hair in people of European ancestry.")
    add("hair_curl", "Straight or curly hair (TCHH)", "Appearance", "C", s, rows,
        ["A small effect: many genes shape hair form, and this variant alone does not predict it."],
        "TCHH L790M in trichohyalin, a protein of the hair follicle, associated with straight hair in Europeans (Medland et al. 2009); direction checked against the GWAS Catalog. " + std_method,
        _ref(traits, ["medland2009"]) + _gwas_refs([V.get("rs11803731", {})]), ok=k is not None)

    # ---------------------------------------------------------------- food and drink
    lct = ["rs4988235", "rs41380347", "rs145946881", "rs41525747", "rs869051967"]
    rows, c = look(*lct)
    typed = [r for r in lct if c.get(r) is not None]
    carriers = [V[r]["label"] for r in typed if c[r] >= 1]
    if not typed:
        s, ok = "None of the lactase-persistence variants is on your chip.", False
    elif carriers:
        ok = True
        s = f"Likely lactase persistent: you carry {', '.join(carriers)}. One copy is enough to keep producing lactase into adulthood, so milk is usually digested well."
    else:
        ok = True
        s = "No lactase-persistence variant found. Most adults without one gradually lose the ability to digest lactose; symptoms vary, and many tolerate small amounts, yogurt and hard cheese."
        if not any(r in typed for r in lct[1:]):
            s += " The African and Arabian variants are not on your chip, so if your ancestry is from East Africa or Arabia this result is less certain."
    add("lactose", "Lactose tolerance", "Food and drink", "A", s, rows,
        ["Genotype predicts whether lactase production continues, not how much milk you can comfortably drink; gut bacteria and habit matter too.",
         f"Variants on your chip: {len(typed)} of {len(lct)}."],
        "Variants in an enhancer of the lactase gene (LCT) that keep it active after weaning: the European -13910T (Enattah et al. 2002) and the African and Arabian variants (Tishkoff et al. 2007). Lactase persistence is dominant. " + std_method,
        _ref(traits, ["enattah2002", "tishkoff2007"]), ok=ok)

    rows, c = look("rs671", "rs1229984")
    a, b = c.get("rs671"), c.get("rs1229984")
    bits = []
    if a is None:
        bits.append("ALDH2 E504K is not on your chip")
    elif a == 2:
        bits.append("Two copies of ALDH2*2: aldehyde dehydrogenase is nearly inactive, so even a little alcohol causes strong flushing, nausea and a racing heart; most people with this genotype drink little or none")
    elif a == 1:
        bits.append("One copy of ALDH2*2: the enzyme that clears acetaldehyde works at a fraction of normal, causing the 'Asian flush' with alcohol. Drinkers with this variant have a substantially raised risk of oesophageal cancer from alcohol")
    else:
        bits.append("No copies of ALDH2*2: acetaldehyde is cleared normally")
    if b is not None:
        bits.append(f"{b} {_plural(b, 'copy', 'copies')} of ADH1B*2, which converts alcohol to acetaldehyde faster" +
                    (" (linked to a lower risk of alcohol dependence, and can add to flushing)" if b else ""))
    add("alcohol_flush", "Alcohol flush", "Food and drink", "A", ". ".join(bits) + ".", rows,
        ["Flushing has other causes (some medicines, rosacea); a normal genotype does not make alcohol safe."],
        "ALDH2 E504K (ALDH2*2) and ADH1B R48H (ADH1B*2), the main inherited differences in alcohol metabolism (Edenberg 2007); cancer risk from Brooks et al. 2009. " + std_method,
        _ref(traits, ["edenberg2007", "brooks2009"]), ok=a is not None)

    rows, c = look("rs762551", "rs2472297", "rs4410790")
    a, b, d = c.get("rs762551"), c.get("rs2472297"), c.get("rs4410790")
    bits = []
    if a is not None:
        bits.append({2: "CYP1A2*1F/*1F: the 'highly inducible' form, which clears caffeine fastest in smokers and heavy coffee drinkers",
                     1: "One copy of CYP1A2*1F: intermediate inducibility of the main caffeine-clearing enzyme",
                     0: "No copies of CYP1A2*1F: the less inducible form, so caffeine clearance rises less with smoking or heavy coffee use"}[a])
    if b is not None:
        bits.append(f"{b} {_plural(b, 'copy', 'copies')} of the CYP1A1/CYP1A2 allele linked to lower blood caffeine (faster clearance)")
    if d is not None:
        bits.append(f"{d} {_plural(d, 'copy', 'copies')} of the AHR allele linked to drinking slightly more coffee")
    add("caffeine", "Caffeine metabolism", "Food and drink", "C", (". ".join(bits) + ".") if bits else "Not on your chip.", rows,
        ["Each variant shifts caffeine clearance or intake only a little; smoking, pregnancy, oral contraceptives and some medicines change it far more."],
        "CYP1A2*1F raises how much the enzyme is induced (Sachse et al. 1999); AHR and CYP1A1/CYP1A2 variants from genome-wide studies of coffee intake and caffeine levels (Cornelis et al. 2011), directions checked against the GWAS Catalog. " + std_method,
        _ref(traits, ["sachse1999", "cornelis2011"]) + _gwas_refs([V.get("rs2472297", {}), V.get("rs4410790", {})]),
        ok=bool(bits))

    rows, c = look("rs713598", "rs1726866", "rs10246939")
    p, a, v = c.get("rs713598"), c.get("rs1726866"), c.get("rs10246939")
    if p is None:
        s, ok, hap = "The main bitter-taste variant (TAS2R38 A49P) is not on your chip.", False, None
    else:
        ok = True
        consistent = a is not None and v is not None and p == a == v
        hap = p
        s = {2: "Two copies of the taster form of TAS2R38 (PAV/PAV): the bitter compounds PTC and PROP taste strongly bitter; some people with this genotype also find broccoli, kale or grapefruit notably bitter.",
             1: "One taster copy of TAS2R38 (PAV/AVI): PTC and PROP taste bitter, usually less intensely than with two copies.",
             0: "No taster copies of TAS2R38 (AVI/AVI): a 'non-taster', for whom PTC and PROP taste weak or not bitter at all."}[p]
        if not consistent:
            s += " Your three TAS2R38 variants do not form the two common haplotypes, so you may carry a rarer form; the call is based on A49P, the most important of the three."
    add("bitter_taste", "Bitter taste (TAS2R38)", "Food and drink", "A", s, rows,
        ["TAS2R38 explains most of the difference in tasting PTC and PROP, but there are about 25 bitter receptors; liking for bitter foods depends on much more."],
        "The bitter-taste receptor TAS2R38: haplotypes of three amino-acid changes (A49P, V262A, I296V) determine sensitivity to PTC (Kim et al. 2003). " + std_method,
        _ref(traits, ["kim2003"]), {"taster_copies": hap}, ok=ok)

    rows, c = look("rs2282679", "rs10741657")
    typed = [r for r in ("rs2282679", "rs10741657") if c.get(r) is not None]
    low = sum(c[r] for r in typed)
    s = (f"{low} of {2 * len(typed)} alleles at {', '.join(V[r]['gene'] for r in typed)} are the ones linked to lower blood vitamin D. "
         "Each lowers the average level slightly; sunlight, diet and supplements matter far more." if typed else "Not on your chip.")
    add("vitamin_d", "Vitamin D levels", "Food and drink", "B", s, rows,
        ["These variants shift average blood levels a little; they do not say whether you are deficient. Only a blood test can."],
        "Variants in GC (vitamin D-binding protein) and CYP2R1 (25-hydroxylase) from genome-wide association studies of blood 25-hydroxyvitamin D (Wang et al. 2010); the lowering allele is taken from the GWAS Catalog. " + std_method,
        _ref(traits, ["wang2010"]) + _gwas_refs([V.get("rs2282679", {}), V.get("rs10741657", {})]), ok=bool(typed))

    rows, c = look("rs601338")
    k = c.get("rs601338")
    s = ("Not on your chip." if k is None else
         "Two copies of FUT2 W143X: a 'non-secretor', who does not release blood-group antigens into saliva and gut mucus. Non-secretors have higher vitamin B12 levels on average." if k == 2 else
         "A 'secretor' (non-secretor status needs two copies of FUT2 W143X" + ("; you carry one)." if k == 1 else ").") +
         " Secretors have slightly lower vitamin B12 levels on average than non-secretors.")
    add("vitamin_b12", "Vitamin B12 and secretor status", "Food and drink", "B", s, rows,
        ["Average differences in blood B12, not deficiency: diet (B12 comes only from animal foods or supplements), absorption and age matter more."],
        "FUT2 controls secretion of ABO antigens into body fluids; its common stop variant W143X is associated with plasma vitamin B12 (Hazra et al. 2008). " + std_method,
        _ref(traits, ["hazra2008"]), ok=k is not None)

    rows, c = look("rs174537")
    k = c.get("rs174537")
    rec = V.get("rs174537", {})
    s = ("Not on your chip." if k is None else
         f"{k} {_plural(k, 'copy', 'copies')} of the FADS1 allele linked to lower blood levels of long-chain omega-6 and omega-3 fats such as arachidonic acid: the body makes these from plant fats less efficiently with each copy.")
    add("omega3", "Omega-3 and omega-6 conversion (FADS1)", "Food and drink", "B", s, rows,
        ["An effect on blood fatty-acid levels, not a dietary recommendation; eating oily fish supplies long-chain omega-3 directly."],
        "FADS1 encodes the delta-5 desaturase that makes arachidonic acid and EPA; rs174537 is the strongest genetic influence on their blood levels (Tanaka et al. 2009). The lowering allele is taken from the GWAS Catalog. " + std_method,
        _ref(traits, ["tanaka2009"]) + _gwas_refs([rec]), ok=k is not None)

    rows, c = look("rs2187668", "rs7454108")
    d2, d8 = c.get("rs2187668"), c.get("rs7454108")
    if d2 is None and d8 is None:
        s, ok = "The HLA tag variants are not on your chip.", False
    else:
        ok = True
        found = [n for n, x in (("HLA-DQ2.5", d2), ("HLA-DQ8", d8)) if x]
        if found:
            s = (f"You carry {' and '.join(found)}, the HLA types found in most people with celiac disease. They are also common in people "
                 "who never develop it, so this alone means little; it matters only alongside symptoms and blood tests.")
        elif d2 is not None and d8 is not None:
            s = ("Neither HLA-DQ2.5 nor HLA-DQ8 found. These are carried by the large majority of people with celiac disease, so "
                 "celiac disease is unlikely; ruling it out is what this test is most useful for.")
        else:
            s = "Only one of the two tag variants is on your chip, and it is negative; the other HLA type could not be checked."
    add("celiac_hla", "Celiac disease HLA types", "Food and drink", "A", s, rows,
        ["Tag SNPs stand in for HLA types and were validated mainly in Europeans; elsewhere they are less reliable.",
         "HLA-DQ2.2 and DQ7, which account for a small share of celiac disease, are not assessed.",
         "This is not a diagnosis; anyone with symptoms should be tested before cutting out gluten."],
        "Single SNPs that tag the HLA-DQ2.5 and HLA-DQ8 haplotypes (Monsuur et al. 2008). " + std_method,
        _ref(traits, ["monsuur2008"]), ok=ok)

    # ---------------------------------------------------------------- fitness
    rows, c = look("rs1815739")
    k = c.get("rs1815739")
    eur = (V.get("rs1815739", {}).get("freq") or {}).get("EUR")
    s = ("Not on your chip." if k is None else
         {2: "XX: no alpha-actinin-3 in your fast-twitch muscle fibres" + (f" (about {100 * eur * eur:.0f}% of Europeans share this)" if eur else "") + ". This genotype is rare among elite sprint and power athletes; for everyday fitness the effect is small.",
          1: "RX: alpha-actinin-3 is produced normally in your fast-twitch fibres (one working copy is enough).",
          0: "RR: alpha-actinin-3 is produced normally in your fast-twitch fibres; RR is slightly more common among elite sprint and power athletes."}[k])
    add("actn3", "Sprint and power muscle variant (ACTN3)", "Fitness", "B", s, rows,
        ["A real but small effect: training, other genes and many other factors decide athletic performance.",
         "People with XX are healthy; the protein's role is covered by alpha-actinin-2."],
        "ACTN3 R577X, a common stop variant (North et al. 1999), associated with elite athletic performance (Yang et al. 2003). " + std_method,
        _ref(traits, ["north1999", "yang2003"]), ok=k is not None)
    return cards
