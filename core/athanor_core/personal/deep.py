"""Personal Genome: file check, deep ancestry, modern ancestry and trait reports.

Every report is a card with an evidence grade (A established, B moderate,
C weak or for fun), the result, what it cannot tell you, the method, sources,
and how it was validated (from the reference pack's CI record).
"""
from __future__ import annotations

import math
import time

import numpy as np

from . import fstats
from .haplo import Tree
from .modern_cards import modern_cards
from .pgx import Panel, pgx_cards
from .traits import Genotyper, trait_cards
from .pca import nearest, normalise, project
from .rawfile import RawGenotypes, file_check, parse
from .refdata import RefData, align, chromosome_calls

Z95 = 1.959964
MIN_SUPPORT = 5        # matching markers along the called lineage needed to report a haplogroup

SRC = {
    "aadr": {"title": "Mallick et al., The Allen Ancient DNA Resource (AADR), Sci. Data 11, 182 (2024)",
             "url": "https://doi.org/10.1038/s41597-024-03031-7"},
    "patterson2012": {"title": "Patterson et al., Ancient admixture in human history, Genetics 192, 1065 (2012)",
                      "url": "https://doi.org/10.1534/genetics.112.145037"},
    "haak2015": {"title": "Haak et al., Massive migration from the steppe was a source for Indo-European languages in Europe, Nature 522, 207 (2015)",
                 "url": "https://doi.org/10.1038/nature14317"},
    "lazaridis2014": {"title": "Lazaridis et al., Ancient human genomes suggest three ancestral populations for present-day Europeans, Nature 513, 409 (2014)",
                      "url": "https://doi.org/10.1038/nature13673"},
    "harney2021": {"title": "Harney et al., Assessing the performance of qpAdm, Genetics 217, iyaa045 (2021)",
                   "url": "https://doi.org/10.1093/genetics/iyaa045"},
    "narasimhan2019": {"title": "Narasimhan et al., The formation of human populations in South and Central Asia, Science 365, eaat7487 (2019)",
                       "url": "https://doi.org/10.1126/science.aat7487"},
    "petr2019": {"title": "Petr et al., Limits of long-term selection against Neandertal introgression, PNAS 116, 1639 (2019)",
                 "url": "https://doi.org/10.1073/pnas.1814338116"},
    "prufer2017": {"title": "Prüfer et al., A high-coverage Neandertal genome from Vindija Cave in Croatia, Science 358, 655 (2017)",
                   "url": "https://doi.org/10.1126/science.aao1887"},
    "meyer2012": {"title": "Meyer et al., A high-coverage genome sequence from an archaic Denisovan individual, Science 338, 222 (2012)",
                  "url": "https://doi.org/10.1126/science.1224344"},
    "patterson2006": {"title": "Patterson, Price & Reich, Population structure and eigenanalysis, PLoS Genet. 2, e190 (2006)",
                      "url": "https://doi.org/10.1371/journal.pgen.0020190"},
    "phylotree": {"title": "van Oven & Kayser, Updated comprehensive phylogenetic tree of global human mitochondrial DNA variation, Hum. Mutat. 30, E386 (2009); PhyloTree Build 17",
                  "url": "https://doi.org/10.1002/humu.20921"},
    "haplogrep3": {"title": "Schönherr et al., Haplogrep 3, Nucleic Acids Res. 51, W263 (2023)",
                   "url": "https://doi.org/10.1093/nar/gkad284"},
    "yleaf": {"title": "Ralf et al., Yleaf: software for human Y-chromosomal haplogroup inference, Mol. Biol. Evol. 35, 1291 (2018)",
              "url": "https://doi.org/10.1093/molbev/msy032"},
    "yfull": {"title": "YFull YTree", "url": "https://www.yfull.com/tree/"},
    "1000g": {"title": "The 1000 Genomes Project Consortium, A global reference for human genetic variation, Nature 526, 68 (2015)",
              "url": "https://doi.org/10.1038/nature15393"},
    "hgdp": {"title": "Bergström et al., Insights into human genetic variation and population history from 929 diverse genomes, Science 367, eaay5012 (2020)",
             "url": "https://doi.org/10.1126/science.aay5012"},
    "alexander2009": {"title": "Alexander, Novembre & Lange, Fast model-based estimation of ancestry in unrelated individuals, Genome Res. 19, 1655 (2009)",
                      "url": "https://doi.org/10.1101/gr.094052.109"},
    "maples2013": {"title": "Maples et al., RFMix: a discriminative modeling approach for rapid and robust local-ancestry inference, Am. J. Hum. Genet. 93, 278 (2013)",
                   "url": "https://doi.org/10.1016/j.ajhg.2013.06.020"},
    "sankararaman2008": {"title": "Sankararaman et al., Estimating local ancestry in admixed populations, Am. J. Hum. Genet. 82, 290 (2008)",
                         "url": "https://doi.org/10.1016/j.ajhg.2007.09.022"},
    "reich2012": {"title": "Reich et al., Reconstructing Native American population history, Nature 488, 370 (2012)",
                  "url": "https://doi.org/10.1038/nature11258"},
    "purcell2007": {"title": "Purcell et al., PLINK: a tool set for whole-genome association and population-based linkage analyses, Am. J. Hum. Genet. 81, 559 (2007)",
                    "url": "https://doi.org/10.1086/519795"},
    "ceballos2018": {"title": "Ceballos et al., Runs of homozygosity: windows into population history and trait architecture, Nat. Rev. Genet. 19, 220 (2018)",
                     "url": "https://doi.org/10.1038/nrg.2017.109"},
}


def _card(cid, title, group, grade, summary, result, limits, method, sources, ref: RefData, ok=True):
    val = (ref.manifest.get("validation") or {}).get(cid)
    return {"id": cid, "title": title, "group": group, "grade": grade, "ok": ok, "summary": summary,
            "result": result, "limits": limits, "method": method,
            "sources": [s if isinstance(s, dict) else SRC[s] for s in sources if isinstance(s, dict) or s in SRC],
            "validation": val}


def _freq_matrix(ref: RefData, names: list[str], user: np.ndarray, snp_mask: np.ndarray):
    """Frequencies (populations + user as last row) on the SNPs in snp_mask."""
    fr = ref.freqs
    idx = [fr["names"].index(n) for n in names]
    alt = fr["alt"][idx][:, snp_mask].astype(float)
    tot = fr["total"][idx][:, snp_mask].astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        F = np.where(tot > 0, alt / tot, np.nan)
    u = user[snp_mask].astype(float)
    u = np.where(u < 0, np.nan, u / 2)
    return np.vstack([F, u])


# ------------------------------------------------------------------ cards
def ancient_ancestry(ref: RefData, dos: np.ndarray, auto: np.ndarray) -> dict:
    models = ref.manifest["models"]
    block = ref.snps["block"][auto]
    results = []
    for mdl in models:
        names = mdl["sources"] + mdl["right"]
        F = _freq_matrix(ref, names, dos, auto)
        t = len(names)
        try:
            r = fstats.qpadm(F, block, t, list(range(len(mdl["sources"]))),
                             list(range(len(mdl["sources"]), len(names))))
        except Exception as e:  # noqa: BLE001
            results.append({"model": mdl["id"], "title": mdl["title"], "error": str(e)})
            continue
        w = np.array(r["weights"])
        se = np.array(r["se"])
        # ADMIXTOOLS' feasibility: proportions between 0 and 1 (a small tolerance for noise at 0)
        feasible = bool(np.all(w > -0.02) and np.all(w < 1.02))
        fit = "fits" if (r["p"] >= 0.05 and feasible) else ("marginal" if (r["p"] >= 0.01 and feasible) else "does not fit")
        imprecise = bool(np.max(se) > 0.15)
        results.append({"model": mdl["id"], "title": mdl["title"], "description": mdl.get("description", ""),
                        "components": [{"source": s, "label": mdl["labels"][s], "proportion": float(wi), "se": float(si),
                                        "ci95": [float(wi - Z95 * si), float(wi + Z95 * si)]}
                                       for s, wi, si in zip(mdl["sources"], w, se)],
                        "p": r["p"], "chi2": r["chi2"], "df": r["df"], "snps": r["snps"], "blocks": r["blocks"],
                        "fit": fit, "feasible": feasible, "imprecise": imprecise, "right": mdl["right"]})
    fitting = [r for r in results if r.get("fit") in ("fits", "marginal")]
    best = max(fitting, key=lambda r: r["p"]) if fitting else None
    if best:
        parts = ", ".join(f"{c['label']} {100 * c['proportion']:.0f}%" for c in best["components"])
        summary = f"{best['title']}: {parts} (model p = {best['p']:.2f})"
    else:
        summary = "None of the ancient-source models fits your genome; see the PCA and closest ancient individuals instead."
    return _card("ancient_ancestry", "Ancient ancestry: hunter-gatherers, farmers, herders", "Deep ancestry", "A",
                 summary, {"models": results, "best": best["model"] if best else None},
                 ["These proportions describe which ancient groups' DNA best explains yours, not a family tree: each source stands for a population, not individual ancestors.",
                  "A model can fit and still be one of several that fit; a model that does not fit means your ancestry includes sources it leaves out.",
                  "Genotyping chips cover a fraction of the reference SNPs, so the intervals are wider than in published population studies."],
                 "qpAdm (Haak et al. 2015): f4-statistics between your genome, the ancient sources and a fixed set of outgroup populations, fitted by generalised least squares with a 5 cM block jackknife. The model's p-value says whether the sources explain your genome's relationships to all outgroups.",
                 ["haak2015", "harney2021", "patterson2012", "lazaridis2014", "narasimhan2019", "aadr"], ref)


def pca_cards(ref: RefData, dos: np.ndarray) -> tuple[dict, dict]:
    anc = ref.ancients
    out = {}
    for name in ref.manifest["pca"]:
        P = ref.pca(name["id"])
        idx = P["snp_idx"]
        z = normalise(dos[idx], P["p"], P["scale"])
        try:
            x, n = project(z, P["loadings"])
        except ValueError as e:
            out[name["id"]] = {"error": str(e)}
            continue
        k = len(x)
        ev = P["eigvals"][:k]
        coords = np.array([a["pc"].get(name["id"]) or [np.nan] * k for a in anc], dtype=float)[:, :k]
        valid = ~np.isnan(coords).any(1)
        # Euclidean distance on the component scores (each axis keeps its own spread)
        ii, dd = nearest(x, coords[valid]) if valid.any() else ([], [])
        vidx = np.flatnonzero(valid)
        near = [{**{k2: anc[vidx[i]][k2] for k2 in ("id", "group", "date_bp", "date_sd", "locality", "country", "lat", "lon",
                                                    "mt", "y", "publication", "doi")},
                 "distance": float(d), "pc": coords[vidx[i]].tolist()} for i, d in zip(ii, dd)]
        out[name["id"]] = {"title": name["title"], "coords": x.tolist(), "snps": n,
                           "eigvals": ev.tolist(),
                           "reference": {"coords": P["ref_coords"][:, :k].tolist(), "group": P["ref_group"].tolist(),
                                         "region": P["ref_region"].tolist()},
                           "nearest_ancients": near}
    world = out.get("world", {})
    west = out.get("west", {})
    pca_card = _card("pca", "Your place among present-day and ancient people (PCA)", "Deep ancestry", "A",
                     "Your genome placed on principal-component maps built from present-day reference populations, with ancient individuals projected onto the same axes.",
                     {"maps": out},
                     ["Distances on a PCA summarise only the main axes of variation; two people close together can still differ in ancestry the axes do not show.",
                      "Ancient genomes are projected from sparse data, so their positions are less certain than present-day ones."],
                     "Principal components of present-day reference genotypes, normalised per SNP as in EIGENSOFT; your genome and the ancient genomes are placed by least squares on the SNPs each has (as smartpca's lsqproject).",
                     ["patterson2006", "lazaridis2014", "aadr", "1000g"], ref)
    near = (west.get("nearest_ancients") or world.get("nearest_ancients") or [])
    near_card = _card("closest_ancients", "Your closest ancient individuals", "Deep ancestry", "B",
                      (f"Closest on the {('West Eurasian' if west.get('nearest_ancients') else 'worldwide')} map: "
                       + ", ".join(f"{n['group']} ({_fmt_bp(n['date_bp'])})" for n in near[:3])) if near else "No ancient individuals could be compared.",
                      {"west": west.get("nearest_ancients", []), "world": world.get("nearest_ancients", [])},
                      ["Closest on the PCA, not your direct ancestors: thousands of ancient people from the same population would sit in the same place.",
                       "Only individuals with published genomes appear; most of the past is unsampled."],
                      "Euclidean distance between component scores on the PCA (each axis keeps its natural spread, so the main axes count most).",
                      ["aadr", "lazaridis2014"], ref)
    return pca_card, near_card


def _fmt_bp(bp) -> str:
    if bp is None:
        return "undated"
    y = 1950 - int(bp)
    return f"{-y:,} BCE" if y <= 0 else f"{y:,} CE"


def archaic_cards(ref: RefData, dos: np.ndarray, auto: np.ndarray) -> list[dict]:
    roles = ref.manifest["archaic_roles"]
    block = ref.snps["block"][auto]
    cards = []
    stats_ref = ref.reference_stats
    # Neanderthal: alpha = f4(Altai, Chimp; X, Africa) / f4(Altai, Chimp; Vindija, Africa)
    names = [roles["altai"], roles["chimp"], roles["vindija"], roles["africa"]]
    F = _freq_matrix(ref, names, dos, auto)
    try:
        r = fstats.f4_ratio(F, block, (0, 1, 4, 3), (0, 1, 2, 3))
        est, se = r["alpha"], r["se"]
        nean = {"estimate": est, "se": se, "ci95": [est - Z95 * se, est + Z95 * se], "snps": r["snps"],
                "reference_groups": stats_ref.get("neanderthal", [])}
        summary = f"{100 * est:.2f}% Neanderthal (95% interval {100 * (est - Z95 * se):.2f}–{100 * (est + Z95 * se):.2f}%)"
        ok = True
    except Exception as e:  # noqa: BLE001
        nean, summary, ok = {"error": str(e)}, "Not enough overlapping SNPs for an estimate.", False
    cards.append(_card("neanderthal", "Neanderthal ancestry", "Deep ancestry", "A", summary, nean,
                       ["Almost everyone with ancestry from outside Africa carries about 1.5 to 2.5% Neanderthal DNA; differences between people are small, so compare the interval with the reference groups rather than reading the number as a rank.",
                        "The estimate uses chip SNPs only; whole-genome data narrows the interval."],
                       "f4-ratio (Petr et al. 2019): the excess of alleles you share with the Altai Neanderthal, relative to Africans, divided by the excess shared by the Vindija Neanderthal; a 5 cM block jackknife gives the interval.",
                       ["petr2019", "prufer2017", "patterson2012", "aadr"], ref, ok))
    # Denisovan affinity: D(Africa, X; Denisova, Neanderthal)
    names = [roles["africa"], roles["denisova"], roles["vindija"]]
    F = _freq_matrix(ref, names, dos, auto)
    try:
        d = fstats.d_stat(F, block, 0, 3, 1, 2)
        den = {**d, "reference_groups": stats_ref.get("denisovan", [])}
        regs = (stats_ref.get("denisovan") or {}).get("regions", [])
        fits = [r["region"] for r in regs if r.get("sd") and r["n"] >= 3 and abs(d["D"] - r["mean"]) <= 2 * r["sd"]]
        summary = (f"D = {d['D']:+.4f} (Z {d['z']:+.1f}): " +
                   (f"within the usual range of {', '.join(fits)}." if fits else "outside the usual range of every reference region.")
                   + " Lower values than your region's would point to Denisovan ancestry.")
        ok = True
    except Exception as e:  # noqa: BLE001
        den, summary, ok = {"error": str(e)}, "Not enough overlapping SNPs.", False
    cards.append(_card("denisovan", "Denisovan affinity", "Deep ancestry", "B", summary, den,
                       ["This is an affinity statistic, not a percentage. Neanderthal ancestry pushes it up and Denisovan ancestry pulls it down, so on its own it mixes the two; its meaning comes from where you sit relative to the reference groups. Oceanian groups, with the most Denisovan ancestry, have the lowest values among non-Africans.",
                        "Detectable Denisovan ancestry is concentrated in Oceania and, at lower levels, in Asia and the Americas."],
                       "D-statistic D(Africans, you; Denisovan, Neanderthal) with a 5 cM block jackknife: positive when you share more derived alleles with the Neanderthal genome than with the Denisovan genome, relative to Africans.",
                       ["meyer2012", "patterson2012", "aadr"], ref, ok))
    # Neanderthal-derived alleles you carry
    A = ref.archaic
    sidx, allele = A["snp_idx"], A["allele"]
    u = dos[sidx].astype(np.int64)
    typed = u >= 0
    cop = np.where(typed, np.where(allele == 1, u, 2 - u), 0)
    groups = sorted(k for k in A if k.startswith("freq_"))
    expected = {k[5:]: float((2 * A[k][typed]).sum()) for k in groups}
    s = ref.snps
    rows = [{"rsid": int(s["rsid"][sidx[j]]), "chrom": int(s["chrom"][sidx[j]]), "pos": int(s["pos"][sidx[j]]),
             "copies": int(cop[j]), **{k[5:]: float(A[k][j]) for k in groups}}
            for j in np.flatnonzero(cop > 0)[:400]]
    copies = cop[typed]
    total = int(copies.sum())
    cards.append(_card("archaic_alleles", "Neanderthal-derived variants you carry", "Deep ancestry", "B",
                       f"{total:,} copies of Neanderthal-matching alleles at {int(typed.sum()):,} typed sites "
                       + ("(" + ", ".join(f"{g} average {v:,.0f}" for g, v in expected.items()) + ")" if expected else ""),
                       {"typed_sites": int(typed.sum()), "copies": total, "sites_carried": int((copies > 0).sum()),
                        "expected_copies_by_group": expected, "variants": rows,
                        "definition": ref.manifest.get("archaic_definition")},
                       ["An allele that matches the Neanderthal genomes and is nearly absent in Africa is very likely Neanderthal-derived, but a few will have other histories.",
                        "Effects on traits are known for only a handful of these variants and are mostly small."],
                       "Sites where the Altai, Vindija and Chagyrskaya Neanderthals carry an allele that has under 1% frequency in sub-Saharan African reference genomes; group averages are the expected copies at the same sites from 1000 Genomes frequencies.",
                       ["prufer2017", "1000g", "aadr"], ref))
    return cards


def haplogroup_cards(ref: RefData, g: RawGenotypes, pos37: np.ndarray, sex: str) -> list[dict]:
    anc = ref.ancients
    cards = []
    for which, chrom, title, srcs in (("mt", 26, "Maternal line: mtDNA haplogroup", ["phylotree", "haplogrep3", "aadr"]),
                                      ("y", 24, "Paternal line: Y haplogroup", ["yleaf", "yfull", "aadr"])):
        if which == "y" and sex == "XX":
            cards.append(_card("y_haplogroup", title, "Deep ancestry", "A",
                               "No Y chromosome in this file (genetic sex XX). A father's or brother's file would show your paternal line.",
                               {"call": None}, [], "", srcs, ref))
            continue
        T = Tree(ref.tree(which))
        obs, orient = T.orient(chromosome_calls(g, chrom, pos37))
        res = T.call(obs)
        res["strand_check"] = orient
        # A call resting on a handful of markers is not reported as a haplogroup
        weak = res.get("call") is not None and res.get("node") is not None and res.get("path_derived", 0) < MIN_SUPPORT
        if weak:
            res = {**res, "weak_call": res["call"], "call": None,
                   "reason": (f"Not enough evidence: the best-supported branch, {res['call']}, rests on only "
                              f"{res['path_derived']} matching marker{'s' if res['path_derived'] != 1 else ''} "
                              f"of {res['typed_mutations']} typed. Of the {orient['on_tree']} markers in your file at "
                              f"positions on this tree, {orient['direct']} matched the tree's alleles directly, "
                              f"{orient['complemented']} after switching strand, and {orient['unmatched']} matched neither.")}
        node = res.get("node") if not weak else None
        carriers, timeline = [], []
        if node is not None:
            key = f"{which}_node"
            same, related = [], []
            for a in anc:
                an = a.get(key)
                if an is None or an < 0:
                    continue
                if T.is_ancestor(node, an):
                    same.append(a)
            same.sort(key=lambda a: -(a.get("date_bp") or 0))
            carriers = [{k: a.get(k) for k in ("id", "group", "date_bp", "date_sd", "locality", "country", "lat", "lon",
                                               "mt", "y", "publication", "doi")} for a in same[:200]]
            # oldest sampled carrier of each branch on the line
            path = T.path(node)
            best = {i: None for i in path}
            for a in anc:
                an = a.get(key)
                if an is None or an < 0 or a.get("date_bp") is None:
                    continue
                j = an
                while j >= 0:
                    if j in best and (best[j] is None or a["date_bp"] > best[j]["date_bp"]):
                        best[j] = a
                    j = int(T.parent[j])
            timeline = [{"branch": T.nodes[i]["name"],
                         "oldest": ({k: best[i].get(k) for k in ("id", "group", "date_bp", "date_sd", "locality", "country",
                                                                "lat", "lon", "publication")} if best[i] else None)}
                        for i in path]
        cid = f"{which}_haplogroup"
        if res.get("call") is None:
            summary = res.get("reason", "No call")
            cards.append(_card(cid, title, "Deep ancestry", "A", summary,
                               {**res, "tree": T.name, "carriers": [], "timeline": [], "markers_in_file": len(obs)},
                               ["Chip files type a limited set of these positions; a whole-genome or dedicated mtDNA/Y test resolves the lineage."],
                               "Your genotypes at the positions defining each branch of the tree; a lineage is reported only when at least "
                               f"{MIN_SUPPORT} markers along it match.", srcs, ref, ok=False))
            continue
        else:
            summary = (f"{res['call']} ({res['path_derived']} supporting and {res['path_ancestral']} conflicting markers "
                       f"of {res['typed_mutations']} typed); {len(carriers)} ancient genomes carry it or a branch below it")
        cards.append(_card(cid, title, "Deep ancestry", "A", summary,
                           {**res, "tree": T.name, "carriers": carriers, "timeline": timeline,
                            "markers_in_file": len(obs)},
                           ["One line among thousands of ancestors: the haplogroup follows only your mother's mother's line (mtDNA) or father's father's line (Y).",
                            "Chips type a limited set of markers, so the call can stop above the most detailed branch; whole-genome data resolves further.",
                            "Ancient carriers share a common ancestor with you on this line; they are not necessarily your ancestors."],
                           "Your genotypes at the positions defining each branch of the tree; the call is the branch whose path from the root is best supported. Ancient carriers are AADR individuals whose published haplogroup is this branch or below it; the timeline gives, for each branch on your line, the oldest dated genome carrying it, a minimum age for that branch.",
                           srcs, ref))
    return cards


def analyse(path: str, ref: RefData | None = None, on_step=None, include=("deep", "modern", "traits", "medicines")) -> dict:
    ref = ref or RefData()
    t0 = time.time()

    def step(msg, f):
        if on_step:
            on_step({"task": "personal", "stage": msg, "fraction": f, "elapsed_seconds": round(time.time() - t0)})

    step("reading the file", 0.02)
    g = parse(path)
    fc = file_check(g)
    step("matching markers to the reference", 0.1)
    al = align(g, ref)
    dos = al["dosage"]
    auto = (ref.snps["chrom"] >= 1) & (ref.snps["chrom"] <= 22) & (dos >= 0)
    fc["reference_match"] = al["stats"]
    cards = [_card("file_check", "Your file", "Your file", "A",
                   f"{fc['vendor']} file, {fc['markers']:,} markers, call rate {100 * fc['call_rate']:.1f}%, "
                   f"{al['stats']['build']}, genetic sex {fc['genetic_sex']}; {int(auto.sum()):,} autosomal SNPs usable for ancestry",
                   fc, ["Genetic sex here means the X and Y chromosomes found in the file."],
                   "Vendor and build from the file header and from agreement of marker positions with the GRCh37 reference; strand checked against the reference alleles.",
                   [], ref)]
    if auto.sum() < 20000:
        raise ValueError(f"only {int(auto.sum())} autosomal SNPs match the reference; the file may be incomplete or not a genotype file")
    step("ancient ancestry models (qpAdm)", 0.15)
    cards.append(ancient_ancestry(ref, dos, auto))
    step("principal components and closest ancient individuals", 0.25)
    cards.extend(pca_cards(ref, dos))
    step("Neanderthal and Denisovan ancestry", 0.33)
    cards.extend(archaic_cards(ref, dos, auto))
    step("haplogroups", 0.4)
    cards.extend(haplogroup_cards(ref, g, al["positions37"], fc["genetic_sex"]))

    def make(cid, title, group, grade, summary, result, limits, method, sources, ok=True):
        return _card(cid, title, group, grade, summary, result, limits, method, sources, ref, ok)

    if "modern" in include:
        cards.extend(modern_cards(ref, g, al, fc, make, step))
    if "traits" in include and ref.traits:
        step("appearance, food and drink, fitness", 0.9)
        cards.extend(trait_cards(ref.traits, Genotyper(g, al["positions37"], al["stats"]["strand"]), make))
    if "medicines" in include and ref.pgx:
        step("medicines", 0.9)
        panel = Panel(ref.pgx_panel, ref.pgx["genes"]) if ref.pgx_panel is not None else None
        cards.extend(pgx_cards(ref.pgx, panel, g, al["positions37"], al["stats"]["strand"], fc["genetic_sex"],
                               make, step))
    step("done", 1.0)
    return {"version": 2, "file": {k: fc[k] for k in ("vendor", "file", "sha256", "markers")},
            "reference_pack": {k: ref.manifest.get(k) for k in ("name", "built", "sources")},
            "seconds": round(time.time() - t0, 1), "cards": cards}
