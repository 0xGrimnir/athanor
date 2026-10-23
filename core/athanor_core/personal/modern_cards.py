"""Modern-ancestry report cards: continental and regional composition,
chromosome painting, Indigenous American ancestry, runs of homozygosity and
X-chromosome ancestry. Methods are in modern.py; the reference panels,
detection floors and painting calibration come from the reference pack."""
from __future__ import annotations

import numpy as np

from .modern import composition, decode_freq, paint, roh, roh_summary, segments
from .rawfile import BASES, MISSING, RawGenotypes

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
RELATIONSHIPS = [("Second cousins", 1 / 64), ("First cousins", 1 / 16),
                 ("Half-siblings, uncle and niece, or double first cousins", 1 / 8),
                 ("Siblings, or parent and child", 1 / 4)]


def _pct(x):
    return f"{100 * x:.0f}%" if x >= 0.995 or x >= 0.1 else f"{100 * x:.1f}%"


def _ci(b):
    return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


def _cont_freqs(P: np.ndarray, cont: np.ndarray, n: np.ndarray, conts: list[int]) -> np.ndarray:
    out = np.empty((len(conts), P.shape[1]))
    for k, c in enumerate(conts):
        rows = np.flatnonzero(cont == c)
        w = n[rows].astype(float)
        out[k] = (w[:, None] * P[rows]).sum(0) / w.sum()
    return out


def modern_cards(ref, g: RawGenotypes, al: dict, fc: dict, make_card, step=lambda *a: None) -> list[dict]:
    M = ref.modern
    if M is None:
        return []
    man = ref.manifest["modern"]
    rids = [str(x) for x in M["region_ids"]]
    rmeta = {r["id"]: r for r in man["regions"]}
    cmeta = man["continents"]
    cont_ids = [str(x) for x in M["continent_ids"]]
    ctitle = {c["id"]: c["title"] for c in cmeta}
    cfloor = {c["id"]: c["detection_floor"] for c in cmeta}
    rcont = M["continent"].astype(int)
    n_ind = M["n_individuals"]
    s = ref.snps
    dos = al["dosage"]
    idx = np.flatnonzero(dos >= 0)
    P = decode_freq(M["P"][:, idx])
    cards = []

    # ------------------------------------------------------------ composition
    step("regional ancestry (supervised admixture)", 0.5)
    res = composition(dos[idx], P, s["block"][idx], reps=100)
    q, boot = res["q"], res["boot"]
    regions = []
    for k, r in enumerate(rids):
        fl = rmeta[r]["detection_floor"]
        regions.append({"id": r, "title": rmeta[r]["title"], "continent": cont_ids[rcont[k]],
                        "proportion": float(q[k]), "ci95": _ci(boot[:, k]), "floor": fl,
                        "detected": bool(q[k] >= max(fl, 0.005)),
                        "reference_groups": rmeta[r]["groups"], "group_names": rmeta[r].get("group_names", {}),
                        "reference_individuals": rmeta[r]["individuals"]})
    conts = []
    for c, cid in enumerate(cont_ids):
        cols = np.flatnonzero(rcont == c)
        v = float(q[cols].sum())
        conts.append({"id": cid, "title": ctitle[cid], "proportion": v, "ci95": _ci(boot[:, cols].sum(1)),
                      "floor": cfloor[cid], "detected": bool(v >= max(cfloor[cid], 0.005))})
    conts.sort(key=lambda x: -x["proportion"])
    regions.sort(key=lambda x: -x["proportion"])
    low = res["snps"] < man["loo"]["snps"]
    undetected_c = sum(c["proportion"] for c in conts if not c["detected"])
    det_c = [c for c in conts if c["detected"]]
    summary_c = ", ".join(f"{c['title']} {_pct(c['proportion'])} ({_pct(c['ci95'][0])}–{_pct(c['ci95'][1])})" for c in det_c)
    if undetected_c > 0.005:
        summary_c += f"; {_pct(undetected_c)} spread across other continents below the detection limits"
    common_limits = [
        "Proportions describe which present-day reference populations your DNA most resembles, not where your ancestors lived at any particular time.",
        "Each region is defined by its reference genomes (listed below); places without reference genomes cannot appear: " + "; ".join(man.get("not_covered", [])) + ".",
        f"Results for {res['snps']:,} SNPs shared with the reference." + (
            f" That is fewer than the {man['loo']['snps']:,} the validation tested, so the intervals may understate the uncertainty." if low else "")]
    result = {"continents": conts, "regions": regions, "snps": res["snps"], "iterations": res["iterations"]}
    cards.append(make_card("continental_ancestry", "Continental ancestry", "Modern ancestry", "A", summary_c, result,
                           common_limits,
                           "Supervised admixture (as in ADMIXTURE's projection mode): the mixture of reference populations' allele frequencies that best explains your genotypes, by maximum likelihood. Each panel's frequencies are shrunk towards the pooled frequency of all panels, so that panels of different sizes compete fairly. 95% intervals from 100 bootstrap resamples of ~5 cM blocks. A proportion below a region's detection floor (the 99th percentile assigned to that region in left-out people from elsewhere) is not reported as found.",
                           ["alexander2009", "1000g", "hgdp", "aadr"]))
    det_r = [r for r in regions if r["detected"]]
    undetected_r = sum(r["proportion"] for r in regions if not r["detected"])
    summary_r = ", ".join(f"{r['title']} {_pct(r['proportion'])}" for r in det_r)
    if undetected_r > 0.005:
        summary_r += f"; {_pct(undetected_r)} below the detection limits"
    cards.append(make_card("regional_ancestry", "Regional ancestry", "Modern ancestry", "B", summary_r, result,
                           common_limits + [
                               "Neighbouring regions share much of their history, so splits between them are less certain than continental totals; regions that the validation could not tell apart have been merged.",
                               "The regional reference sets are small (tens of people each), so a region stands for its reference groups, not a whole country."],
                           "The same model at regional level: each region's own reference panel. Regions that people left out of the panels could not be assigned back to (on average under 60%, on 80,000 SNPs) were merged with their nearest neighbour.",
                           ["alexander2009", "1000g", "hgdp", "aadr"]))

    # ------------------------------------------------------------ painting
    step("chromosome painting", 0.62)
    qc = {c["id"]: c["proportion"] for c in conts}
    used = [c for c in cont_ids if qc[c] >= max(cfloor[c], 0.01)]
    if not used:
        used = [conts[0]["id"]]
    ci = [cont_ids.index(c) for c in used]
    chrom = s["chrom"][idx].astype(np.int64)
    pos = s["pos"][idx].astype(np.int64)
    if len(used) == 1:
        a = np.zeros(len(idx), np.int8); b = a.copy(); post = np.ones(len(idx), np.float32)
    else:
        Pc = _cont_freqs(P, rcont, n_ind, ci)
        qn = np.array([qc[c] for c in used]); qn /= qn.sum()
        pp = man["painting"]
        pr = paint(dos[idx], Pc, chrom, s["cm"][idx].astype(np.float64), qn,
                   generations=pp["generations"], temper=pp["temper"])
        a, b, post = pr["a"], pr["b"], pr["posterior"]
    segs = segments(chrom, pos, a, b, post)
    for sg in segs:
        sg["a"] = used[sg["a"]]; sg["b"] = used[sg["b"]]
    share = {c: float(((a == k).sum() + (b == k).sum()) / (2 * len(a))) for k, c in enumerate(used)}
    ends = {int(c): int(pos[chrom == c].max()) for c in np.unique(chrom)}
    paint_sum = ", ".join(f"{ctitle[c]} {_pct(v)}" for c, v in sorted(share.items(), key=lambda x: -x[1]))
    pr_summary = ("All of your genome is assigned to one continent: " + ctitle[used[0]] + "." if len(used) == 1 else
                  f"{len(segs)} segments across your 22 autosomes; shares along the painting: {paint_sum}.")
    cards.append(make_card("chromosome_painting", "Chromosome painting", "Modern ancestry", "B", pr_summary,
                           {"continents": [{"id": c, "title": ctitle[c]} for c in used], "segments": segs,
                            "chromosome_ends": ends, "mean_posterior": float(np.mean(post)) if len(post) else None},
                           ["Chip data are unphased: each segment gives the pair of ancestries of your two copies at that point, but not which copy is which, so the two colours in a segment are not 'mother' and 'father'.",
                            "Painting is at continental level only; regional differences are too small to follow along a chromosome with chip data.",
                            "Short segments (under 20 SNPs) are absorbed into their neighbours."],
                           f"Hidden Markov model over the pair of continental ancestries along each chromosome, with switches at recombination ({man['painting']['generations']:.0f} generations since mixing assumed) and SNP evidence down-weighted (factor {man['painting']['temper']}) to allow for linkage; the factor was fitted on simulated genomes built from left-out reference individuals.",
                           ["maples2013", "sankararaman2008", "1000g", "hgdp"]))

    # ------------------------------------------------------------ Indigenous American
    am = next((c for c in conts if c["id"] == "americas"), None)
    if am is not None:
        am_segs = [sg for sg in segs if "americas" in (sg["a"], sg["b"])]
        mb = sum((sg["end"] - sg["start"]) / 1e6 * ((sg["a"] == "americas") + (sg["b"] == "americas")) for sg in am_segs)
        if am["detected"]:
            s_am = (f"About {_pct(am['proportion'])} of your DNA (95% interval {_pct(am['ci95'][0])}–{_pct(am['ci95'][1])}) "
                    f"matches the Indigenous American reference populations, in {len(am_segs)} painted segments covering about {mb:.0f} Mb of chromosome copies.")
        else:
            s_am = (f"No Indigenous American ancestry detected: the estimate, {_pct(am['proportion'])}, is below the detection limit of "
                    f"{_pct(am['floor'])} (the most assigned to this group in people with no Indigenous American ancestry).")
        cards.append(make_card("indigenous_american", "Indigenous American ancestry", "Modern ancestry", "B", s_am,
                               {**am, "segments": am_segs, "segment_mb": mb,
                                "reference_groups": next((r["reference_groups"] for r in regions if r["id"] == "americas"), [])},
                               ["Genetic ancestry is not tribal membership. Each Nation decides its own citizenship, usually by documented descent from its rolls; no DNA test can establish it.",
                                "The reference genomes come mostly from Central and South America (Karitiana, Surui, Maya, Pima and others); Indigenous North American ancestry is detected as Indigenous American but cannot be placed regionally.",
                                "A single ancestor six or more generations back may have passed on little or none of their DNA, so a family story can be true even when nothing is detected.",
                                "East Asian and Siberian ancestry share deep history with Indigenous Americans; the detection limit was set so that people of East Asian ancestry are not reported as Indigenous American."],
                               "The Indigenous American component of the continental model, with its detection limit from left-out reference individuals, and the segments assigned to it by chromosome painting.",
                               ["alexander2009", "reich2012", "hgdp"]))

    # ------------------------------------------------------------ runs of homozygosity
    step("runs of homozygosity", 0.7)
    pos37 = al["positions37"]
    m = (g.chrom >= 1) & (g.chrom <= 22) & (pos37 > 0)
    order = np.lexsort((pos37[m], g.chrom[m]))
    ch = g.chrom[m][order].astype(np.int64)
    ps = pos37[m][order].astype(np.int64)
    a1, a2 = g.a1[m][order], g.a2[m][order]
    dup = np.r_[False, (ch[1:] == ch[:-1]) & (ps[1:] == ps[:-1])]
    ch, ps, a1, a2 = ch[~dup], ps[~dup], a1[~dup], a2[~dup]
    miss = a1 == MISSING
    het = (~miss) & (a1 != a2)
    rsegs = roh(ch, ps, het, miss)
    rs = roh_summary(rsegs)
    rr = ref.roh_reference.get("groups", {})
    top_titles = [r["title"] for r in det_r[:2]]
    near = {gname: v for gname, v in rr.items() if v.get("region") in top_titles}
    long_mb = rs["by_length_mb"]["10+"]
    if long_mb == 0:
        s_roh = f"{rs['total_mb']:.0f} Mb in {rs['segments']} runs, none 10 Mb or longer: no sign that your parents were related."
    elif long_mb < 50:
        s_roh = (f"{rs['total_mb']:.0f} Mb in {rs['segments']} runs, {long_mb:.0f} Mb of it in runs 10 Mb or longer: a modest amount, "
                 "seen in many people and consistent with a shared ancestor several generations back.")
    else:
        s_roh = (f"{rs['total_mb']:.0f} Mb in {rs['segments']} runs, {long_mb:.0f} Mb in runs 10 Mb or longer ({100 * long_mb / 2881:.1f}% of the autosomes). "
                 "This is in the range seen when parents are related (for example as second cousins or closer) or in small, isolated communities. On its own it says nothing about health.")
    cards.append(make_card("roh", "Runs of homozygosity", "Modern ancestry", "A", s_roh,
                           {**rs, "segments_list": rsegs, "markers": int(len(ps)),
                            "reference_groups": near or rr, "relationships": [{"parents": n, "expected_froh": f} for n, f in RELATIONSHIPS]},
                           ["Long runs come from a recent shared ancestor of your two parents; shorter runs (1-5 Mb) mostly reflect the size of ancestral populations, and are larger in historically small or isolated groups.",
                            "Expected values for related parents are averages; actual amounts vary widely around them.",
                            f"Measured on {len(ps):,} markers in your file; runs are harder to detect at lower densities."],
                           "PLINK 1.9's scanning-window method (--homozyg: windows of 50 SNPs with at most 1 heterozygous and 5 missing calls; runs of at least 100 SNPs and 1,000 kb, no gap over 1,000 kb, at least one SNP per 50 kb), checked against PLINK in CI. Reference groups were measured the same way on half of the 1240K SNPs, about the density of a consumer chip.",
                           ["purcell2007", "ceballos2018", "aadr"]))

    # ------------------------------------------------------------ X chromosome
    step("X-chromosome ancestry", 0.76)
    xs_pos = M["x_pos"].astype(np.int64)
    x_card = None
    if len(xs_pos):
        male = fc.get("genetic_sex") == "XY"
        flip = {"the reference's + strand": False, "the opposite strand": True}.get(al["stats"].get("strand"))
        xm = (g.chrom == 23) & (pos37 > 0) & (g.a1 != MISSING)
        look = {int(p_): i for i, p_ in zip(np.flatnonzero(xm), pos37[xm])}
        gx = np.full(len(xs_pos), -1.0)
        nx = np.zeros(len(xs_pos))
        refb, altb = M["x_ref"], M["x_alt"]
        for j, p_ in enumerate(xs_pos):
            i = look.get(int(p_))
            if i is None or g.a1[i] > 3 or g.a2[i] > 3:
                continue
            r_, a_ = BASES[refb[j]], BASES[altb[j]]
            x1, x2 = BASES[g.a1[i]], BASES[g.a2[i]]
            amb = COMP[r_] == a_
            if flip is True:
                x1, x2 = COMP[x1], COMP[x2]
            elif flip is None and (amb or not {x1, x2} <= {r_, a_}):
                continue
            if not {x1, x2} <= {r_, a_}:
                continue
            k = (x1 == a_) + (x2 == a_)
            if male:
                if x1 != x2:
                    continue
                gx[j], nx[j] = k // 2, 1
            else:
                gx[j], nx[j] = k, 2
        ok = gx >= 0
        if ok.sum() >= 1500:
            conts_used = [c for c in cont_ids if qc[c] >= max(cfloor[c], 0.01)] or [conts[0]["id"]]
            cidx = [cont_ids.index(c) for c in cont_ids]
            Px = decode_freq(M["Px"][:, ok])
            Pxc = _cont_freqs(Px, rcont, n_ind, cidx)
            xr = composition(gx[ok].astype(np.int8), Pxc, M["x_block"][ok], n_alleles=nx[ok], reps=100)
            xc = [{"id": c, "title": ctitle[c], "proportion": float(xr["q"][k]), "ci95": _ci(xr["boot"][:, k]),
                   "autosomal": qc[c]} for k, c in enumerate(cont_ids)]
            xc.sort(key=lambda x: -x["proportion"])
            shown = [x for x in xc if x["proportion"] >= 0.02 or x["id"] in conts_used]
            s_x = ("X chromosome: " + ", ".join(f"{x['title']} {_pct(x['proportion'])}" for x in shown)
                   + ". Autosomes: " + ", ".join(f"{ctitle[c]} {_pct(qc[c])}" for c in conts_used) + ".")
            x_card = make_card("x_ancestry", "X-chromosome ancestry", "Modern ancestry", "B", s_x,
                               {"continents": xc, "snps": int(ok.sum()), "male": male},
                               [("You have one X, inherited from your mother, so this is her side of the family only."
                                 if male else
                                 "You have two X chromosomes, one from each parent; your father's came from his mother, so this side of the family is followed mostly through women."),
                                "Women carry two X chromosomes and men one, so ancestry on the X over-represents female ancestors; a difference from your autosomes can reflect who partnered with whom in the past.",
                                f"Measured on {int(ok.sum()):,} X-chromosome SNPs (outside the pseudo-autosomal regions), with wider intervals than the autosomes."],
                               "The continental admixture model applied to the X chromosome alone, counting one copy in men; 95% intervals from a block bootstrap.",
                               ["alexander2009", "1000g", "hgdp"])
    if x_card is None:
        x_card = make_card("x_ancestry", "X-chromosome ancestry", "Modern ancestry", "B",
                           "Not enough X-chromosome markers in your file match the reference panel.", {"snps": 0}, [], "",
                           [], ok=False)
    cards.append(x_card)
    return cards
