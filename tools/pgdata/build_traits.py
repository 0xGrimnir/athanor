"""Resolve the trait-report variants against public sources and write
traits.json into the reference pack.

For every variant in traits_spec.VARIANTS: GRCh37 position and alleles and
1000 Genomes superpopulation frequencies (Ensembl GRCh37 REST), the amino
acid each allele codes (Ensembl VEP), and GWAS Catalog associations matching
the spec's keywords. Each check in the spec must pass. References are looked
up in PubMed so titles, journals and DOIs come from the record, not from
hand-typed citations.

Usage: python build_traits.py PACK_DIR
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from traits_spec import REFS, VARIANTS  # noqa: E402

OUT = Path(sys.argv[1])
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
SUPER = ("AFR", "AMR", "EAS", "EUR", "SAS")
fails = []


def get(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "athanor-ci", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode("utf-8", errors="replace"))
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(3 * (k + 1))
    raise RuntimeError(f"{url}: {err!r}")


def fail(msg):
    print("FAIL:", msg, flush=True)
    fails.append(msg)


def ensembl(rs):
    d = get(f"https://grch37.rest.ensembl.org/variation/human/{rs}?pops=1")
    maps = [m for m in d.get("mappings", []) if re.fullmatch(r"\d+|X", str(m.get("seq_region_name")))]
    if len(maps) != 1:
        raise RuntimeError(f"{rs}: {len(maps)} primary-assembly mappings")
    m = maps[0]
    alleles = m["allele_string"].split("/")
    freqs = {}
    for p in d.get("populations", []):
        pop = p.get("population", "")
        if pop.startswith("1000GENOMES:phase_3:") and pop.split(":")[-1] in SUPER:
            freqs.setdefault(pop.split(":")[-1], {})[p["allele"]] = float(p["frequency"])
    chrom = 23 if m["seq_region_name"] == "X" else int(m["seq_region_name"])
    return {"chrom": chrom, "pos": int(m["start"]), "ref": alleles[0], "alts": alleles[1:], "freqs": freqs,
            "minor": d.get("minor_allele"), "maf": d.get("MAF")}


def vep_aa(rs, ref):
    """allele -> amino acid on protein-coding transcripts (the reference allele
    codes the first amino acid of each change)."""
    d = get(f"https://grch37.rest.ensembl.org/vep/human/id/{rs}")
    out = {}
    for v in d if isinstance(d, list) else []:
        for t in v.get("transcript_consequences", []):
            aa = t.get("amino_acids")
            if not aa or t.get("biotype") != "protein_coding":
                continue
            parts = aa.split("/")
            if len(parts) == 2:
                out.setdefault(ref, set()).add(parts[0])
                out.setdefault(t["variant_allele"], set()).add(parts[1])
            elif len(parts) == 1:
                out.setdefault(t["variant_allele"], set()).add(parts[0])
    return {k: sorted(v) for k, v in out.items()}


def gwas(rs, keywords, alleles):
    page, rows = 0, []
    while True:
        d = get(f"https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id={rs}&size=100&page={page}")
        rows += (d.get("_embedded") or {}).get("associations", [])
        if page + 1 >= (d.get("page") or {}).get("totalPages", 1):
            break
        page += 1
    hits = []
    for a in rows:
        text = " ".join((a.get("reported_trait") or []) + [e.get("efo_trait", "") for e in a.get("efo_traits", [])]).lower()
        if not any(k in text for k in keywords):
            continue
        ea = None
        for x in a.get("snp_allele") or []:
            if x.get("rs_id") == rs:
                ea = (x.get("effect_allele") or "").upper()
        if not ea or ea == "?":
            continue
        if ea not in alleles and COMP.get(ea) in alleles:
            ea = COMP[ea]                                   # reported on the other strand
        if ea not in alleles:
            continue
        beta = str(a.get("beta") or "")
        orv = a.get("or_per_copy_num")
        if "increase" in beta:
            direction = "increase"
        elif "decrease" in beta:
            direction = "decrease"
        elif orv not in (None, "-") and float(orv) != 1:
            direction = "increase" if float(orv) > 1 else "decrease"
        else:
            continue
        hits.append({"effect_allele": ea, "direction": direction, "beta": beta if beta != "-" else None,
                     "or": orv, "p": a.get("p_value"), "trait": (a.get("reported_trait") or [""])[0],
                     "efo": [e.get("efo_trait") for e in a.get("efo_traits", [])], "pubmed_id": a.get("pubmed_id"),
                     "first_author": a.get("first_author"), "accession": a.get("accession_id")})
    hits.sort(key=lambda h: h["p"] if h["p"] is not None else 1)
    return hits


def other_allele(e, info):
    cands = [a for a in [info["ref"]] + info["alts"] if a != e and len(a) == 1]
    if len(cands) == 1:
        return cands[0]
    # several alternative alleles: the one common in 1000 Genomes
    best = max(cands, key=lambda a: max((info["freqs"].get(s, {}).get(a, 0) for s in SUPER), default=0))
    return best


out = {}
for rs, spec in VARIANTS.items():
    try:
        info = ensembl(rs)
    except Exception as e:  # noqa: BLE001
        fail(f"{rs}: Ensembl lookup failed: {e}")
        continue
    alleles = [info["ref"]] + info["alts"]
    rec = {"rsid": rs, "gene": spec["gene"], "label": spec["label"], "chrom": info["chrom"], "pos37": info["pos"],
           "ref": info["ref"], "alts": info["alts"]}
    effect = spec.get("effect")
    if "aa" in spec:
        aa = vep_aa(rs, info["ref"])
        rec["amino_acids"] = aa
        coding = [a for a, v in aa.items() if spec["aa"] in v]
        if effect is None:
            if len(coding) != 1:
                fail(f"{rs}: {len(coding)} alleles code {spec['aa']} ({aa})")
                continue
            effect = coding[0]
        elif effect not in coding and aa:
            fail(f"{rs}: allele {effect} does not code {spec['aa']} (VEP: {aa})")
        elif not aa:
            print(f"note: {rs}: VEP gave no amino-acid changes; relying on the other checks")
    if spec.get("gwas"):
        hits = gwas(rs, spec["gwas"], alleles)
        rec["gwas"] = hits[:8]
        if not hits:
            fail(f"{rs}: no GWAS Catalog association with a direction for {spec['gwas']}")
            continue
        best = hits[0]
        if effect is None and spec.get("orient") == "lower":
            effect = best["effect_allele"] if best["direction"] == "decrease" else other_allele(best["effect_allele"], info)
            rec["direction"] = "lower"
        elif effect is not None and spec.get("direction"):
            same = best["effect_allele"] == effect
            implied = best["direction"] if same else ("decrease" if best["direction"] == "increase" else "increase")
            if implied != spec["direction"]:
                fail(f"{rs}: the GWAS Catalog's strongest association ({best['trait']}, {best['first_author']}) "
                     f"implies allele {effect} {implied}s the trait; the spec says {spec['direction']}")
            rec["direction"] = spec["direction"]
        agree = sum(1 for h in hits if (h["effect_allele"] == best["effect_allele"]) == (h["direction"] == best["direction"]))
        rec["gwas_agreement"] = f"{agree} of {len(hits)} matching associations agree in direction"
    if effect is None:
        fail(f"{rs}: no effect allele")
        continue
    if effect not in alleles:
        fail(f"{rs}: effect allele {effect} not among Ensembl alleles {alleles}")
        continue
    if spec.get("alt") and effect not in info["alts"]:
        fail(f"{rs}: {effect} is not an alternative allele ({alleles})")
    rec["effect"] = effect
    rec["other"] = other_allele(effect, info)
    rec["freq"] = {s: info["freqs"].get(s, {}).get(effect, 0.0 if info["freqs"].get(s) else None) for s in SUPER}
    for s, (lo, hi) in (spec.get("expect") or {}).items():
        f = rec["freq"].get(s)
        if f is None or not (lo <= f <= hi):
            fail(f"{rs}: allele {effect} frequency in {s} is {f}, expected {lo}-{hi}")
    out[rs] = rec
    print(f"{rs:12s} {spec['gene']:14s} effect {effect} other {rec['other']} chr{info['chrom']}:{info['pos']} "
          f"freq {rec['freq']} {rec.get('gwas_agreement', '')}", flush=True)
    time.sleep(0.2)

ids = ",".join(REFS.values())
d = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id={ids}")
refs = {}
for key, pmid in REFS.items():
    r = (d.get("result") or {}).get(pmid)
    if not r or "title" not in r:
        fail(f"PubMed {pmid} ({key}) not found")
        continue
    doi = next((x["value"] for x in r.get("articleids", []) if x.get("idtype") == "doi"), None)
    authors = r.get("authors") or []
    first = authors[0]["name"] if authors else r.get("sortfirstauthor", "")
    refs[key] = {"pmid": pmid, "title": r["title"].rstrip("."), "journal": r.get("source"),
                 "year": (r.get("pubdate") or "")[:4], "first_author": first, "doi": doi,
                 "url": f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "traits.json").write_text(json.dumps({"built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                             "sources": ["Ensembl GRCh37 REST (variation, VEP)", "1000 Genomes phase 3",
                                                         "GWAS Catalog REST API v2", "NCBI PubMed"],
                                             "variants": out, "references": refs}, indent=1))
print(f"{len(out)} variants, {len(refs)} references; {len(fails)} failure(s)")
sys.exit(1 if fails else 0)
