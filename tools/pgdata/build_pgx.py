"""Pharmacogenomics data for the reference pack: CPIC tables and the 1000
Genomes high-coverage genotypes at every CPIC-defined position.

Sources (all fetched here, in CI):
  CPIC database API (https://api.cpicpgx.org, CC0): level-A gene-drug pairs,
      allele definitions, allele functions and frequencies, the function ->
      phenotype tables, prescribing recommendations, guideline publications.
  Ensembl REST: the GRCh38 reference sequence around each defined position
      (to place indels and check every CPIC reference allele against the
      genome), and GRCh37 positions of each rsID (to check the liftover).
  UCSC hg38ToHg19 chain: GRCh37 positions, which consumer chips report.
  1000 Genomes high-coverage phased panel (Byrska-Bishop et al. 2022):
      genotypes of 3,202 people at the defined positions, read with bcftools
      straight from the EBI FTP site.
  NCBI PubMed: titles of the guideline publications and of the papers cited.

Writes PACK/pgx.json and PACK/pgx_1kg.npz, and VALDATA/pgx_1kg.vcf (every
defined position, unphased, in PharmCAT's own representation when its
positions file is given) for the comparison with PharmCAT.

Usage: python build_pgx.py PACK VALDATA CHAIN.gz BCFTOOLS [PHARMCAT_POSITIONS.vcf]
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.personal.liftover import Liftover  # noqa: E402

PACK, VAL, CHAIN, BCFTOOLS = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
PC_POS = sys.argv[5] if len(sys.argv) > 5 else None
API = "https://api.cpicpgx.org/v1/"
KG = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
      "20220422_3202_phased_SNV_INDEL_SV/")
KG_PED = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
          "20130606_g1k_3202_samples_ped_population.txt")
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT",
         "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
CHROM = {**{f"chr{i}": i for i in range(1, 23)}, "chrX": 23, "chrY": 24, "chrM": 26, "chrMT": 26}
REFS = {  # key: (PMID, words the PubMed title must contain)
    "relling2011": ("21270786", ["Clinical Pharmacogenetics Implementation Consortium"]),
    "pharmcat": ("31306493", ["PharmCAT"]),
    "weedon2021": ("33589468", ["SNP chips", "rare"]),
    "1000g30x": ("36055201", ["1000 Genomes", "602 trios"]),
}
NOT_DEFINED = {"HLA-A": "HLA alleles need HLA typing; CPIC has no SNP definitions for them",
               "HLA-B": "HLA alleles need HLA typing; CPIC has no SNP definitions for them"}
PRIMARY = {"DPYD": {"rs75017182"}}    # HapB3: the intronic SNP decides (CPIC/PharmCAT)
fails = []


def fail(msg):
    print("FAIL:", msg, flush=True)
    fails.append(msg)


def get(url, tries=5, raw=False, headers=None):
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "athanor-ci", "Accept": "application/json",
                                                       **(headers or {})})
            with urllib.request.urlopen(req, timeout=180) as r:
                body = r.read()
            return body if raw else json.loads(body.decode("utf-8", errors="replace"))
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(4 * (k + 1))
    raise RuntimeError(f"{url}: {err!r}")


def api(path):
    """All rows of a CPIC API query (PostgREST pages with limit/offset)."""
    out, off, page = [], 0, 5000
    sep = "&" if "?" in path else "?"
    while True:
        rows = get(f"{API}{path}{sep}limit={page}&offset={off}")
        out += rows
        if len(rows) < page:
            return out
        off += page


def inlist(xs):
    return "in.(" + ",".join(urllib.parse.quote(f'"{x}"', safe='"') for x in xs) + ")"


# ------------------------------------------------------------------ CPIC
t0 = time.time()
pairs = api("pair_view?cpiclevel=eq.A")
genes = sorted({p["genesymbol"] for p in pairs})
print(f"CPIC level A: {len(pairs)} pairs, genes {genes}", flush=True)
gene_rows = {g["symbol"]: g for g in api(f"gene?symbol={inlist(genes)}")}
skipped = {}
gdata = {}
for sym in genes:
    defs = api(f"allele_definition?genesymbol=eq.{urllib.parse.quote(sym)}&select=*,allele_location_value(*,sequence_location(*))")
    if not defs:
        skipped[sym] = NOT_DEFINED.get(sym, "no allele definitions in CPIC")
        print(f"{sym}: skipped ({skipped[sym]})")
        continue
    alleles = api(f"allele?genesymbol=eq.{urllib.parse.quote(sym)}")
    gdata[sym] = {"defs": defs, "alleles": alleles}
    print(f"{sym}: {len(defs)} definitions, {len(alleles)} allele rows", flush=True)

results = [r for r in api(f"gene_result?genesymbol={inlist(list(gdata))}")]
lookups = api("gene_result_lookup")
by_ph = {}
for l in lookups:
    by_ph.setdefault(l["phenotypeid"], []).append(l)
guideline_ids = sorted({p["guidelinename"] for p in pairs})
guides = api("guideline?select=*,publication(title,pmid,year)")
gid = {g["name"]: g for g in guides}
recs = api(f"recommendation?guidelineid={inlist(sorted({gid[n]['id'] for n in guideline_ids if n in gid}))}")
drug_ids = sorted({p["drugid"] for p in pairs})
drugs = {d["drugid"]: d for d in api(f"drug?drugid={inlist(drug_ids)}")}
print(f"{len(results)} gene results, {len(lookups)} lookups, {len(recs)} recommendations, {len(drugs)} drugs "
      f"({time.time() - t0:.0f} s)", flush=True)
VAL.mkdir(parents=True, exist_ok=True)
(VAL / "pgx_cpic_raw.json").write_text(json.dumps({"pairs": pairs, "genes": gene_rows, "data": gdata,
                                                   "results": results, "lookups": lookups, "recs": recs,
                                                   "drugs": drugs, "guidelines": guides}))

# ------------------------------------------------------------------ genome
lift = Liftover(CHAIN)


def ensembl_seq(chrom: str, start: int, end: int) -> str:
    c = chrom.replace("chr", "").replace("M", "MT") if chrom != "chrM" else "MT"
    s = get(f"https://rest.ensembl.org/sequence/region/human/{c}:{start}..{end}:1?content-type=text/plain",
            raw=True, headers={"Accept": "text/plain"})
    s = s.decode().strip().upper()
    if len(s) != end - start + 1:
        raise RuntimeError(f"sequence {chrom}:{start}-{end} has length {len(s)}")
    return s


def ensembl37(rsids: list[str]) -> dict:
    out = {}
    for i in range(0, len(rsids), 150):
        body = json.dumps({"ids": rsids[i:i + 150]}).encode()
        for k in range(5):
            try:
                req = urllib.request.Request("https://grch37.rest.ensembl.org/variation/homo_sapiens", data=body,
                                             headers={"Content-Type": "application/json", "Accept": "application/json",
                                                      "User-Agent": "athanor-ci"})
                with urllib.request.urlopen(req, timeout=180) as r:
                    out.update(json.loads(r.read()))
                break
            except Exception as e:  # noqa: BLE001
                if k == 4:
                    fail(f"Ensembl GRCh37 lookup failed: {e!r}")
                time.sleep(5 * (k + 1))
    return out


REP = re.compile(r"^([ACGT]+)\((\d+)\)$")


def expand(s: str, ref_is_snv: bool) -> list[str]:
    """CPIC allele cell -> concrete allele strings."""
    s = s.strip()
    if " or " in s:
        return [x for part in s.split(" or ") for x in expand(part, ref_is_snv)]
    if ref_is_snv and len(s) == 1 and s in IUPAC:
        return list(IUPAC[s])
    return [s]


def allele_seq(a: str, refcell: str, window: str, off: int):
    """Local sequence after replacing the CPIC event at window offset ``off``
    (0-based position of the CPIC position in the window) with allele ``a``.
    Returns (sequence, span) or (None, why)."""
    m_ref = REP.match(refcell)
    if refcell == "del":
        start, span = off + 1, 0                         # insertion after the position
    elif m_ref:
        unit = m_ref.group(1)
        n = 0
        while window[off + n * len(unit): off + (n + 1) * len(unit)] == unit:
            n += 1
        start, span = off, n * len(unit)
    elif refcell.startswith("del") and re.fullmatch(r"del[ACGT]+", refcell):
        start, span = off, len(refcell) - 3
    elif re.fullmatch(r"[ACGT]+", refcell):
        start, span = off, len(refcell)
    else:
        return None, f"reference cell {refcell!r} not understood"
    m = REP.match(a)
    if m:
        rep = m.group(1) * int(m.group(2))
    elif a == "del" or re.fullmatch(r"del[ACGT]*", a):
        rep = ""
    elif re.fullmatch(r"ins[ACGT]+", a):
        rep = a[3:] if span else a[3:]
    elif re.fullmatch(r"[ACGT]+", a):
        rep = a
    else:
        return None, f"allele {a!r} not understood"
    return window[:start] + rep + window[start + span:], span


# ------------------------------------------------------------------ per gene
out_genes = {}
positions_all = []          # (gene, k, chrom, pos38, window_start, window, seqs per allele)
for sym, gd in gdata.items():
    grow = gene_rows.get(sym, {})
    chrom_name = grow.get("chr") or ""
    chrom = CHROM.get(chrom_name)
    if chrom is None:
        fail(f"{sym}: chromosome {chrom_name!r}")
        continue
    defs = gd["defs"]
    refdef = [d for d in defs if d.get("matchesreferencesequence")]
    # CPIC's reference row is the first row of the definition table; the API flags the row matching the
    # reference sequence. Where none is flagged, use the allele named *1 or 'Reference'.
    named = {d["name"]: d for d in defs}
    ref_row = next((named[n] for n in ("*1", "Reference", "B (reference)") if n in named), None) or (refdef[0] if refdef else None)
    if ref_row is None:
        fail(f"{sym}: no reference allele row")
        continue
    locs = {}
    for d in defs:
        for v in d["allele_location_value"]:
            locs[v["locationid"]] = v["sequence_location"]
    loc_ids = sorted(locs, key=lambda i: (locs[i]["position"], i))
    refcell = {}
    for v in ref_row["allele_location_value"]:
        refcell[v["locationid"]] = v["variantallele"]
    for lid in loc_ids:
        if lid not in refcell:
            m = re.match(r"g\.\d+([ACGT])>", locs[lid].get("chromosomelocation") or "")
            if m:
                refcell[lid] = m.group(1)
            else:
                fail(f"{sym}: no reference cell for {locs[lid].get('chromosomelocation')}")
                refcell[lid] = None
    p38 = np.array([locs[l]["position"] for l in loc_ids], dtype=np.int64)
    lo, hi = int(p38.min()) - 120, int(p38.max()) + 120
    try:
        region = ensembl_seq(chrom_name, lo, hi)
    except Exception as e:  # noqa: BLE001
        fail(f"{sym}: reference sequence: {e!r}")
        continue
    if chrom == 26:
        p37 = p38.copy()                           # rCRS in both builds
    else:
        p37 = lift.map(chrom_name, p38)
    positions = []
    cell_index = []
    for k, lid in enumerate(loc_ids):
        L = locs[lid]
        rc = refcell[lid] or ""
        one = len(rc) == 1 and rc in IUPAC
        rset = expand(rc, one)
        cells = set()
        for d in defs:
            for v in d["allele_location_value"]:
                if v["locationid"] == lid:
                    cells.update(expand(v["variantallele"], one))
        p = int(L["position"])
        first = region[p - lo] if (one and len(rset) > 1 and region[p - lo] in rset) else rset[0]
        concrete = [first] + [c for c in rset if c != first] + sorted(c for c in cells if c not in rset)
        snv = all(len(c) == 1 and c in "ACGT" for c in concrete)
        off = p - lo
        if snv:
            gbase = region[off]
            if gbase not in concrete:
                concrete.append(gbase)
            genome = concrete.index(gbase)
            seqs = None
            kind = "snv"
        else:
            kind = "indel"
            w0 = max(0, off - 60)
            window = region[w0: off + 120]
            seqs = []
            for a in concrete:
                s, info = allele_seq(a, rc, window, off - w0)
                if s is None:
                    fail(f"{sym} {L.get('chromosomelocation')}: {info}")
                seqs.append(s)
            gwin = window
            genome = next((j for j, s in enumerate(seqs) if s == gwin), None)
            if genome is None:
                fail(f"{sym} {L.get('chromosomelocation')}: no CPIC allele matches the GRCh38 sequence "
                     f"({rc!r} at {chrom_name}:{p})")
        if snv and region[off] not in rset:
            print(f"  {sym} {L.get('dbsnpid') or L.get('chromosomelocation')}: CPIC reference row has {rc}, "
                  f"GRCh38 has {region[off]}")
        chip = None
        if kind == "indel" and seqs and all(s is not None for s in seqs) and len(concrete) == 2:
            short, long_ = sorted(range(2), key=lambda j: len(seqs[j]))
            if len(seqs[short]) != len(seqs[long_]):
                chip = {"D": short, "I": long_}
        rec = {"rsid": L.get("dbsnpid"), "pos38": p, "pos37": int(p37[k]) if p37[k] > 0 else None,
               "hgvs": L.get("chromosomelocation"), "name": L.get("name"), "protein": L.get("proteinlocation"),
               "kind": kind, "alleles": concrete, "refset": [concrete.index(c) for c in rset], "genome": genome,
               "chip": chip}
        if L.get("dbsnpid") in PRIMARY.get(sym, set()):
            rec["primary"] = True
        positions.append(rec)
        positions_all.append((sym, k, chrom_name, p, seqs, window if kind == "indel" else None,
                              (w0 + lo) if kind == "indel" else None))
        cell_index.append({c: concrete.index(c) for c in concrete})
    # allele rows -> function, activity, frequency
    arows = {a["name"]: a for a in gd["alleles"]}
    alleles = []
    for d in sorted(defs, key=lambda d: (d is not ref_row, d["name"])):
        vals = {v["locationid"]: v["variantallele"] for v in d["allele_location_value"]}
        df = []
        for k, lid in enumerate(loc_ids):
            if d is ref_row:
                df.append(positions[k]["refset"])
            elif lid in vals:
                rc = refcell[lid] or ""
                df.append(sorted({cell_index[k][c] for c in expand(vals[lid], len(rc) == 1 and rc in IUPAC)}))
            else:
                df.append(None)
        a = arows.get(d["name"], {})
        alleles.append({"name": d["name"], "ref": d is ref_row, "structural": bool(d.get("structuralvariation")),
                        "def": df, "function": a.get("clinicalfunctionalstatus"), "activity": a.get("activityvalue"),
                        "freq": a.get("frequency") or {}, "strength": a.get("strength")})
    gres = [r for r in results if r["genesymbol"] == sym]
    rlist, phen = [], {}
    for r in gres:
        phen.setdefault(r["result"], {"ehrpriority": r.get("ehrpriority"), "consultation": r.get("consultationtext"),
                                      "freq": r.get("frequency") or {}})
        for l in by_ph.get(r["id"], []):
            rlist.append({"key": l["lookupkey"], "result": r["result"], "score": l.get("totalactivityscore"),
                          "description": l.get("description"), "priority": r.get("ehrpriority")})
    out_genes[sym] = {"chrom": chrom, "chrom_name": chrom_name, "lookup": grow.get("lookupmethod"),
                      "positions": positions, "alleles": alleles, "results": rlist, "phenotypes": phen}
    n_ind = sum(p["kind"] == "indel" for p in positions)
    print(f"{sym:8s} {chrom_name:5s} {len(positions):4d} positions ({n_ind} indels), {len(alleles)} alleles, "
          f"{len(rlist)} result rows, lookup {grow.get('lookupmethod')}", flush=True)
    if not rlist:
        fail(f"{sym}: no phenotype lookup rows")

# GRCh37 check against Ensembl for every SNV with an rsID
rs = sorted({p["rsid"] for g in out_genes.values() for p in g["positions"] if p["kind"] == "snv" and p["rsid"]
             and g["chrom"] != 26})
e37 = ensembl37(rs)
checked = bad = 0
for sym, g in out_genes.items():
    for p in g["positions"]:
        if p["kind"] != "snv" or not p["rsid"] or g["chrom"] == 26:
            continue
        rec = e37.get(p["rsid"])
        if not rec:
            continue
        maps = [m for m in rec.get("mappings", []) if str(m.get("seq_region_name")) == g["chrom_name"].replace("chr", "")]
        if not maps:
            continue
        checked += 1
        if p["pos37"] != int(maps[0]["start"]):
            bad += 1
            print(f"  {sym} {p['rsid']}: liftover {p['pos37']}, Ensembl GRCh37 {maps[0]['start']}")
print(f"GRCh37 positions: {checked} SNVs checked against Ensembl, {bad} differ", flush=True)
if checked < 0.8 * len(rs) or bad:
    fail(f"GRCh37 liftover: {bad} of {checked} SNV positions differ from Ensembl ({len(rs)} rsIDs)")

# ------------------------------------------------------------------ 1000 Genomes
PA = {(e[0], e[1]): e for e in positions_all}
ped = get(KG_PED, raw=True).decode().splitlines()
hdr = ped[0].split()
rows = [dict(zip(hdr, l.split())) for l in ped[1:] if l.strip()]
info = {r["SampleID"]: r for r in rows}
children = {r["SampleID"] for r in rows if r.get("FatherID", "0") in info and r.get("MotherID", "0") in info}


def kg_file(chrom_name):
    if chrom_name == "chrX":
        return KG + "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    return KG + f"1kGP_high_coverage_Illumina.{chrom_name}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


samples = None
haps = {}
for sym, g in out_genes.items():
    if g["chrom"] == 26:
        print(f"{sym}: the 1000 Genomes panel has no mitochondrial genotypes")
        continue
    P = g["positions"]
    lo = min(p["pos38"] for p in P) - 150
    hi = max(p["pos38"] for p in P) + 150
    url = kg_file(g["chrom_name"])
    cmd = [BCFTOOLS, "view", "-H", "-r", f"{g['chrom_name']}:{lo}-{hi}", url]
    for k in range(4):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0:
            break
        time.sleep(10 * (k + 1))
    if r.returncode != 0:
        fail(f"{sym}: bcftools {r.stderr[-300:]}")
        continue
    if samples is None:
        hs = subprocess.run([BCFTOOLS, "query", "-l", url], capture_output=True, text=True).stdout.split()
        samples = hs
    n = len(samples)
    H = np.full((2 * n, len(P)), 255, dtype=np.uint8)
    for k, p in enumerate(P):
        H[:, k] = p["genome"] if p["genome"] is not None else 255
    by_pos = {}
    for k, p in enumerate(P):
        by_pos.setdefault(p["pos38"], []).append(k)
    indel_k = [k for k, p in enumerate(P) if p["kind"] == "indel"]
    lines = r.stdout.splitlines()
    used = 0
    for line in lines:
        f = line.split("\t", 9)
        pos, ref, alts = int(f[1]), f[3], f[4].split(",")
        if any(a.startswith("<") for a in alts):
            continue                                       # structural variant records
        targets = []
        if len(ref) == 1 and all(len(a) == 1 for a in alts):
            for k in by_pos.get(pos, []):
                if P[k]["kind"] == "snv":
                    targets.append((k, {j + 1: (P[k]["alleles"].index(a) if a in P[k]["alleles"] else 15)
                                        for j, a in enumerate(alts)}))
        else:
            for k in indel_k:
                p = P[k]
                if abs(pos - p["pos38"]) > 80:
                    continue
                # apply the record alone to the position's window and compare sequences
                _, _, _, _, seqs, window, w0 = PA[(sym, k)]
                o = pos - w0
                if o < 0 or window[o:o + len(ref)] != ref or not seqs:
                    continue
                m = {}
                for j, a in enumerate(alts):
                    s = window[:o] + a + window[o + len(ref):]
                    if s in seqs:
                        m[j + 1] = seqs.index(s)
                if m:
                    targets.append((k, m))
        if not targets:
            continue
        gts = f[9].split("\t")
        used += 1
        for k, m in targets:
            for si, gt in enumerate(gts):
                a = re.split(r"[|/]", gt.split(":")[0])
                for h, x in enumerate(a[:2]):
                    if x not in (".", "0") and int(x) in m:
                        H[2 * si + h, k] = m[int(x)]
                if len(a) == 1:                            # haploid (male X)
                    H[2 * si + 1, k] = H[2 * si, k]
    haps[sym] = H
    carriers = int((H != np.array([p["genome"] if p["genome"] is not None else 255 for p in P])).any(1).sum())
    print(f"{sym}: {len(lines)} 1000 Genomes records in the region, {used} at CPIC positions; "
          f"{carriers} of {2 * n} haplotypes differ from GRCh38 somewhere", flush=True)

pop = np.array([info.get(s, {}).get("Superpopulation", "") for s in samples or []])
male = np.array([info.get(s, {}).get("Sex", "") in ("1", "male") for s in samples or []])
related = np.array([s in children for s in samples or []])
np.savez_compressed(PACK / "pgx_1kg.npz", samples=np.array(samples or []), superpop=pop, male=male, related=related,
                    **{f"hap_{s}": h for s, h in haps.items()})
print(f"1000 Genomes: {len(samples or [])} samples, {int((~related).sum())} unrelated", flush=True)

# VCF for PharmCAT (unphased; PharmCAT's representation)
if PC_POS and samples:
    pcl = [l.rstrip("\n").split("\t") for l in open(PC_POS) if not l.startswith("#")]
    vout = open(VAL / "pgx_1kg.vcf", "w")
    vout.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    for c in sorted({g["chrom_name"] for g in out_genes.values()}, key=lambda c: CHROM.get(c, 99)):
        vout.write(f"##contig=<ID={c},assembly=GRCh38.p13,species=\"Homo sapiens\">\n")
    vout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    mapped = unmapped = 0
    written = {}
    for f in pcl:
        chrom_name, pos, ref, alts = f[0], int(f[1]), f[3], f[4].split(",")
        hit = None
        for sym, g in out_genes.items():
            if g["chrom_name"] != chrom_name or sym not in haps:
                continue
            for k, p in enumerate(g["positions"]):
                if p["kind"] == "snv" and p["pos38"] == pos and len(ref) == 1:
                    amap = {j: ([ref] + alts).index(a) if a in [ref] + alts else None for j, a in enumerate(p["alleles"])}
                    hit = (sym, k, amap)
                elif p["kind"] == "indel" and abs(p["pos38"] - pos) <= 80 and len(ref) + max(map(len, alts)) > 2:
                    _, _, _, _, seqs, window, w0 = PA[(sym, k)]
                    o = pos - w0
                    if not seqs or o < 0 or window[o:o + len(ref)] != ref:
                        continue
                    vs = [window[:o] + a + window[o + len(ref):] for a in [ref] + alts]
                    amap = {j: (vs.index(s) if s in vs else None) for j, s in enumerate(seqs)}
                    if sum(v is not None for v in amap.values()) >= 2:
                        hit = (sym, k, amap)
                if hit:
                    break
            if hit:
                break
        if not hit:
            unmapped += 1
            continue
        mapped += 1
        sym, k, amap = hit
        written.setdefault(sym, []).append(k)
        H = haps[sym][:, k]
        gts = []
        for si in range(len(samples)):
            a, b = amap.get(int(H[2 * si])), amap.get(int(H[2 * si + 1]))
            if a is None or b is None:
                gts.append("./.")
            elif g["chrom"] == 23 and male[si]:
                gts.append(f"{a}")
            else:
                gts.append(f"{min(a, b)}/{max(a, b)}")
        vout.write(f"{chrom_name}\t{pos}\t{f[2]}\t{ref}\t{','.join(alts)}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
    vout.close()
    (VAL / "pgx_vcf_positions.json").write_text(json.dumps(written))
    print(f"PharmCAT positions: {mapped} written with 1000 Genomes genotypes, {unmapped} not in the CPIC genes used here")

# ------------------------------------------------------------------ drugs
drug_out = []
for did in drug_ids:
    d = drugs.get(did)
    if not d:
        continue
    prs = [p for p in pairs if p["drugid"] == did]
    rr = [r for r in recs if r["drugid"] == did]
    gene_set = sorted({g for r in rr for g in (r.get("lookupkey") or {})} | {p["genesymbol"] for p in prs})
    gname = prs[0]["guidelinename"]
    gl = gid.get(gname, {})
    drug_out.append({"name": d["name"], "drugid": did, "genes": gene_set, "level_a": sorted({p["genesymbol"] for p in prs}),
                     "guideline": {"name": gname, "url": prs[0].get("guidelineurl"),
                                   "publications": [{"pmid": x.get("pmid"), "title": x.get("title"), "year": x.get("year")}
                                                    for x in gl.get("publication") or []]},
                     "recs": [{"key": r.get("lookupkey") or {}, "population": (r.get("population") or "").strip(),
                               "implications": r.get("implications"), "recommendation": r.get("drugrecommendation"),
                               "classification": r.get("classification"), "comments": r.get("comments")} for r in rr]})
print(f"{len(drug_out)} drugs with level-A pairs", flush=True)

# ------------------------------------------------------------------ references
ids = ",".join(v[0] for v in REFS.values())
d = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id={ids}")
refs = {}
for key, (pmid, words) in REFS.items():
    r = (d.get("result") or {}).get(pmid)
    if not r or "title" not in r:
        fail(f"PubMed {pmid} ({key}) not found")
        continue
    if not all(w.lower() in r["title"].lower() for w in words):
        fail(f"PubMed {pmid} ({key}) title {r['title']!r} lacks {words}")
    doi = next((x["value"] for x in r.get("articleids", []) if x.get("idtype") == "doi"), None)
    authors = r.get("authors") or []
    refs[key] = {"pmid": pmid, "title": r["title"].rstrip("."), "journal": r.get("source"),
                 "year": (r.get("pubdate") or "")[:4], "first_author": authors[0]["name"] if authors else "",
                 "doi": doi, "url": f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"}

cpic_release = None
try:
    cpic_release = get("https://api.github.com/repos/cpicpgx/cpic-data/releases?per_page=1")[0]["tag_name"]
except Exception:  # noqa: BLE001
    pass
PACK.mkdir(parents=True, exist_ok=True)
(PACK / "pgx.json").write_text(json.dumps({
    "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "source": {"name": "CPIC database API", "url": "https://api.cpicpgx.org", "accessed": time.strftime("%Y-%m-%d"),
               "release": cpic_release, "licence": "CC0 1.0",
               "notice": "CPIC guidelines and data are subject to updates; see cpicpgx.org for the current content."},
    "genes": out_genes, "skipped": skipped, "drugs": drug_out, "references": refs}, indent=1))
print(f"pgx.json: {len(out_genes)} genes, {len(drug_out)} drugs; {len(fails)} failure(s) ({time.time() - t0:.0f} s)")
sys.exit(1 if fails else 0)
