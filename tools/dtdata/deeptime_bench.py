"""Deep Time step 1 against published ancient genomes.

Two comparisons, both against numbers other people published:

1. 1240K-capture individuals (Mathieson et al., Nature 528, 499, 2015; reads
   ENA PRJEB11450). Athanor reads each BAM and calls pseudo-haploid genotypes
   at the 1240K sites; these are compared with the genotypes the Allen Ancient
   DNA Resource distributes for the same individuals, and the genetic sex with
   the sex published in its annotations. The AADR's genotype file is read with
   HTTP range requests, one individual's row at a time, so the 5 GB file is
   never downloaded.

2. A high-coverage archaic genome (Pruefer et al., Science 358, 655, 2017).
   A region of the published alignment is read remotely; Athanor calls diploid
   genotypes and compares them with the study's own call set for the same
   region, inside the study's own filter mask, and compares heterozygosity
   computed from both.

Usage: python deeptime_bench.py WORK_DIR [--individuals N] [--region CHR:START-END]
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))
from athanor_core.deeptime import calls as C           # noqa: E402
from athanor_core.deeptime.bam import Genome           # noqa: E402
from athanor_core.deeptime.sites import damage_profile, pileup   # noqa: E402

UA = {"User-Agent": "athanor-ci"}
DATAVERSE = "https://dataverse.harvard.edu/api/access/datafile"
AADR_FILES = {"geno": 13994829, "ind": 13994513, "snp": 13994514, "anno": 13994515}
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
EVA = "http://cdna.eva.mpg.de/neandertal/Vindija"
BASES = "ACGT"
CODE = {b: i for i, b in enumerate(BASES)}

WORK = Path(sys.argv[1])
WORK.mkdir(parents=True, exist_ok=True)
args = sys.argv[2:]


def opt(name, default):
    return args[args.index(name) + 1] if name in args else default


N_IND = int(opt("--individuals", "4"))
REGION = opt("--region", "21:30000000-32000000")
ARCHAIC = opt("--archaic", "Altai")
results, report = [], {}


def gate(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def info(msg):
    print(f"[INFO] {msg}", flush=True)


def fetch(url: str, dest: Path, expect_md5: str | None = None, tries: int = 3) -> Path:
    if dest.exists() and dest.stat().st_size:
        return dest
    for k in range(tries):
        try:
            t0 = time.time()
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=1800) as r, \
                    open(dest.with_suffix(dest.suffix + ".part"), "wb") as fh:
                h = hashlib.md5()
                n = 0
                while True:
                    b = r.read(1 << 22)
                    if not b:
                        break
                    fh.write(b)
                    h.update(b)
                    n += len(b)
            if expect_md5 and h.hexdigest() != expect_md5:
                raise ValueError(f"md5 {h.hexdigest()} != {expect_md5}")
            dest.with_suffix(dest.suffix + ".part").rename(dest)
            info(f"fetched {dest.name} {n / 1e6:.1f} MB in {time.time() - t0:.0f} s")
            return dest
        except Exception as e:  # noqa: BLE001
            info(f"fetch {url} failed ({e}); attempt {k + 1} of {tries}")
            time.sleep(5)
    raise RuntimeError(f"could not fetch {url}")


class AADRRemote:
    """The AADR's 1240K release, with the genotype file read by range request."""

    def __init__(self, work: Path):
        self.snp_path = fetch(f"{DATAVERSE}/{AADR_FILES['snp']}", work / "aadr.snp")
        self.ind_path = fetch(f"{DATAVERSE}/{AADR_FILES['ind']}", work / "aadr.ind")
        self.anno_path = fetch(f"{DATAVERSE}/{AADR_FILES['anno']}", work / "aadr.anno")
        chrom, pos, a1, a2 = [], [], [], []
        for line in self.snp_path.read_text().splitlines():
            f = line.split()
            chrom.append(int(f[1])); pos.append(int(f[3]))
            a1.append(CODE.get(f[4], -1)); a2.append(CODE.get(f[5], -1))
        self.chrom = np.array(chrom, dtype=np.int16)
        self.pos = np.array(pos, dtype=np.int64)
        self.a1 = np.array(a1, dtype=np.int8)
        self.a2 = np.array(a2, dtype=np.int8)
        self.nsnp = len(pos)
        self.rlen = (self.nsnp + 3) // 4
        self.ids = [line.split()[0] for line in self.ind_path.read_text().splitlines()]
        rows = list(csv.reader(open(self.anno_path, encoding="utf-8", errors="replace"), delimiter="\t"))
        self.hdr = rows[0]
        self.anno = {r[0]: r for r in rows[1:] if r}
        self.geno_url = f"{DATAVERSE}/{AADR_FILES['geno']}"
        head = self._range(0, 48).split(b"\x00")[0].decode().split()
        assert head[0] == "TGENO" and int(head[2]) == self.nsnp, head
        self.nind = int(head[1])
        info(f"AADR: {self.nsnp:,} sites, {self.nind:,} individuals, {len(self.anno):,} annotation rows")

    def _range(self, start: int, length: int) -> bytes:
        req = urllib.request.Request(self.geno_url, headers={**UA, "Range": f"bytes={start}-{start + length - 1}"})
        with urllib.request.urlopen(req, timeout=1800) as r:
            if r.status != 206:
                raise RuntimeError(f"the server ignored the range request (status {r.status})")
            return r.read()

    def column(self, name_pattern: str) -> int | None:
        return next((i for i, h in enumerate(self.hdr) if re.search(name_pattern, h, re.I)), None)

    def field(self, gid: str, pattern: str):
        r = self.anno.get(gid)
        i = self.column(pattern)
        if not r or i is None:
            return None
        v = r[i].strip() if i < len(r) else ""
        return None if v in ("", "..", "n/a", "N/A") else v

    def genetic_id(self, alias: str) -> str | None:
        """The AADR's own identifier for a sample: the alias, or the alias with
        one of its data-type suffixes (.AG shotgun, .SG, .DG, .TW)."""
        if alias in self.anno and alias in self.ids:
            return alias
        cands = [g for g in self.ids if g.split(".")[0] == alias and g in self.anno]
        if not cands:
            return None
        return sorted(cands, key=lambda g: (len(g.split(".")), g))[0]

    def dosage(self, gid: str) -> np.ndarray:
        """Alternative-allele (allele 2) dosage for one individual, -1 missing."""
        i = self.ids.index(gid)
        raw = np.frombuffer(self._range(48 + i * self.rlen, self.rlen), dtype=np.uint8)
        g = np.empty(self.rlen * 4, dtype=np.int8)
        for k in range(4):
            g[k::4] = (raw >> (6 - 2 * k)) & 3
        g = g[: self.nsnp]
        return np.where(g == 3, -1, 2 - g).astype(np.int8)


def ena_runs(project: str) -> list[dict]:
    fields = "run_accession,sample_alias,submitted_ftp,submitted_bytes,submitted_md5"
    url = f"{ENA}?accession={project}&result=read_run&fields={fields}&format=tsv"
    rows = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=600).read().decode().splitlines()
    hdr = rows[0].split("\t")
    out = []
    for line in rows[1:]:
        d = dict(zip(hdr, line.split("\t")))
        ftp = [x for x in d["submitted_ftp"].split(";") if x.endswith(".bam") and "1240k" in x.lower()]
        if not ftp:
            continue
        sizes = d["submitted_bytes"].split(";")
        md5s = d.get("submitted_md5", "").split(";")
        i = d["submitted_ftp"].split(";").index(ftp[0])
        out.append({"run": d["run_accession"], "alias": d["sample_alias"], "url": "https://" + ftp[0],
                    "bytes": int(sizes[i]) if i < len(sizes) and sizes[i] else 0,
                    "md5": md5s[i] if i < len(md5s) else ""})
    out.sort(key=lambda r: r["bytes"])
    return out


def run(cmd: list, timeout: int = 7200, stdout=None) -> tuple[int, str, str]:
    """Binary-safe: these servers sometimes put raw bytes in stderr."""
    r = subprocess.run([str(c) for c in cmd], stdout=stdout or subprocess.PIPE, stderr=subprocess.PIPE,
                       timeout=timeout)
    out = "" if stdout else r.stdout.decode("utf-8", "replace")
    return r.returncode, out, r.stderr.decode("utf-8", "replace")


def samtools(*cmd):
    rc, out, err = run(["samtools", *cmd])
    if rc:
        raise RuntimeError(f"samtools {cmd[0]} failed: {err[-400:]}")
    return out


def genotype_bam(path: Path, aadr: AADRRemote, mode: str = "pseudo_haploid", trim: int = 0,
                 min_qual: int = 20, min_mapq: int = 30, min_depth: int = 8) -> dict:
    """Calls at every 1240K site on the autosomes and the X."""
    dos = np.full(aadr.nsnp, -1, dtype=np.int8)
    dmg = {}
    covered = 0
    with Genome(str(path), min_mapq=min_mapq) as g:
        build = g.build
        sex = g.sex()
        for c in range(1, 24):
            idx = np.nonzero(aadr.chrom == c)[0]
            if not len(idx):
                continue
            idx = idx[np.argsort(aadr.pos[idx])]
            name = "X" if c == 23 else str(c)
            if g.name_of(name) is None:
                continue
            pu = pileup(g, name, aadr.pos[idx])
            if mode == "diploid":
                r = C.panel_genotypes(pu, aadr.a1[idx], aadr.a2[idx], min_qual=min_qual, trim=trim,
                                      min_depth=min_depth)
            else:
                r = C.pseudo_haploid(pu, aadr.a1[idx], aadr.a2[idx], seed=100 + c, min_qual=min_qual, trim=trim)
            dos[idx] = r["dosage"]
            covered += int(r.get("sites_covered", 0))
            if c <= 22:
                d = damage_profile(pu, aadr.a1[idx], aadr.a2[idx], min_qual=min_qual)
                for key in ("ct5", "ga3"):
                    a = dmg.setdefault(key, [np.zeros(25), np.zeros(25)])
                    tot = np.array(d[key]["observations"], dtype=float)
                    freq = np.array([0.0 if x is None else x for x in d[key]["frequency"]], dtype=float)
                    a[0] += freq * tot
                    a[1] += tot
    damage = {}
    for key, (hit, tot) in dmg.items():
        with np.errstate(invalid="ignore", divide="ignore"):
            f = np.where(tot > 0, hit / np.maximum(tot, 1), np.nan)
        bg = float(hit[15:].sum() / tot[15:].sum()) if tot[15:].sum() else float("nan")
        damage[key] = {"terminal": float(f[0]), "background": bg, "excess": float(f[0] - bg),
                       "observations_at_first_base": int(tot[0]),
                       "profile": [round(float(x), 4) for x in f[:10]]}
    return {"dosage": dos, "sites_covered": covered, "sites_called": int((dos >= 0).sum()),
            "build": build, "sex": sex, "damage": damage}


def compare_dosage(mine: np.ndarray, theirs: np.ndarray) -> dict:
    both = (mine >= 0) & (theirs >= 0)
    if not both.any():
        return {"shared": 0}
    m, t = mine[both], theirs[both]
    agree = m == t
    # pseudo-haploid calls are homozygous; a published heterozygote (from a
    # diploid caller) can only half-agree, so count those separately
    het = t == 1
    return {"shared": int(both.sum()), "agree": float(agree.mean()),
            "agree_excluding_published_heterozygotes": float(agree[~het].mean()) if (~het).any() else None,
            "published_heterozygous_share": float(het.mean()),
            "mine_called": int((mine >= 0).sum()), "theirs_called": int((theirs >= 0).sum()),
            "mine_only": int(((mine >= 0) & (theirs < 0)).sum()),
            "theirs_only": int(((mine < 0) & (theirs >= 0)).sum())}


def step_1240k(aadr: AADRRemote) -> None:
    runs = ena_runs("PRJEB11450")
    info(f"{len(runs)} BAMs in PRJEB11450; smallest {runs[0]['bytes'] / 1e6:.0f} MB")
    done, rows = 0, []
    for r in runs:
        if done >= N_IND:
            break
        gid = aadr.genetic_id(r["alias"])
        if gid is None:
            info(f"{r['alias']}: no AADR annotation row; skipped")
            continue
        bam = fetch(r["url"], WORK / f"{r['alias']}.bam", r["md5"] or None)
        if not (bam.with_suffix(".bam.bai")).exists():
            samtools("index", bam)
        got = genotype_bam(bam, aadr)
        theirs = aadr.dosage(gid)
        cmp = compare_dosage(got["dosage"], theirs)
        row = {"individual": gid, "run": r["run"], "bytes": r["bytes"], "build": got["build"],
               "published": {"sex": aadr.field(gid, r"^molecular sex"),
                             "snps_hit_1240k": aadr.field(gid, r"1240k snpset"),
                             "coverage_1240k": aadr.field(gid, r"^mean coverage on 1"),
                             "damage": aadr.field(gid, r"^damage rate in first"),
                             "y_haplogroup": aadr.field(gid, r"^y haplogroup in terminal"),
                             "mt_haplogroup": aadr.field(gid, r"^mtdna haplogroup")},
               "ours": {"sites_called": got["sites_called"], "sites_covered": got["sites_covered"],
                        "sex": got["sex"]["call"], "r_y": got["sex"].get("r_y"),
                        "damage": got["damage"]},
               "genotypes": cmp}
        rows.append(row)
        print(json.dumps(row, default=float)[:1200], flush=True)
        done += 1
        os.remove(bam)                     # the runner's disk is small
        for extra in (bam.with_suffix(".bam.bai"),):
            if extra.exists():
                os.remove(extra)
    report["individuals"] = rows
    if not rows:
        gate("1240K individuals: genotypes agree with the AADR's published calls", False, "no individuals ran")
        return
    worst = min(r["genotypes"].get("agree_excluding_published_heterozygotes") or 0 for r in rows)
    gate("Pseudo-haploid calls from the reads agree with the AADR's published genotypes for the same individuals "
         "(sites called by both, excluding the published heterozygotes a single allele cannot match): at least 97%",
         worst >= 0.97,
         "; ".join(f"{r['individual']} {100 * (r['genotypes'].get('agree_excluding_published_heterozygotes') or 0):.2f}% "
                   f"of {r['genotypes']['shared']:,}" for r in rows))
    sexes = [(r["individual"], r["ours"]["sex"], (r["published"]["sex"] or "?")[:1]) for r in rows]
    ok = all(mine == {"F": "XX", "M": "XY"}.get(theirs, "?") for _, mine, theirs in sexes)
    gate("Genetic sex matches the sex published for each individual", ok,
         "; ".join(f"{i}: {mine} vs {t}" for i, mine, t in sexes))
    dmg = []
    for r in rows:
        ours = [r["ours"]["damage"].get(k, {}).get("excess") for k in ("ct5", "ga3")]
        mine = None if any(x is None for x in ours) else float(np.mean(ours))
        try:
            theirs = float(r["published"]["damage"])
        except (TypeError, ValueError):
            theirs = None
        dmg.append((r["individual"], mine, theirs))
        r["ours"]["terminal_deamination"] = mine
    have = [(i, m, t) for i, m, t in dmg if m is not None and t is not None]
    ok = bool(have) and all(abs(m - t) <= max(0.02, 0.4 * t) for _, m, t in have)
    gate("Terminal deamination matches the damage rate published for each individual (within 40%, or 0.02)",
         ok, "; ".join(f"{i}: {m:.3f} vs {t:.3f}" for i, m, t in have) or "no published damage rates to compare")
    measurable = [(i, m) for i, m, _ in dmg if m is not None]
    gate("Terminal deamination is measurable on every individual", len(measurable) == len(rows),
         "; ".join(f"{i}: {m:.3f}" for i, m in measurable))


def published_calls(url: str, region: str, min_gq: int = 30) -> dict:
    """Genotypes from a published VCF: position -> (filter, bases). Only calls
    the study itself reports confidently (its own genotype quality)."""
    out = {}
    cmd = ["bcftools", "query", "-r", region, "-f", "%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%GQ]\n", url]
    rc, text, err = run(cmd)
    if rc:
        raise RuntimeError(f"bcftools failed: {err[-400:]}")
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        pos, ref, alt, filt, gt, gq = int(f[0]), f[1], f[2], f[3], f[4], f[5]
        try:
            if float(gq) < min_gq:
                continue
        except ValueError:
            continue
        if len(ref) != 1 or ref not in BASES:
            continue
        alleles = [ref] + [a for a in alt.split(",") if a != "."]
        if any(len(a) != 1 or a not in BASES for a in alleles[1:]):
            continue
        g = re.split(r"[/|]", gt)
        if len(g) != 2 or any(x == "." for x in g):
            continue
        try:
            bases = tuple(sorted(alleles[int(x)] for x in g))
        except (ValueError, IndexError):
            continue
        out[pos] = (filt, bases)
    return out


def step_archaic() -> None:
    chrom, span = REGION.split(":")
    start, end = (int(x) for x in span.split("-"))
    folder = {"Altai": "Altai_realigned", "Vindija": "Vindija33.19"}[ARCHAIC]
    vcf = f"{EVA}/VCF/{'Altai' if ARCHAIC == 'Altai' else 'Vindija33.19'}/chr{chrom}_mq25_mapab100.vcf.gz"
    mask = f"{EVA}/FilterBed/{'Altai' if ARCHAIC == 'Altai' else 'Vindija33.19'}/chr{chrom}_mask.bed.gz"
    listing = urllib.request.urlopen(urllib.request.Request(f"{EVA}/bam/Pruefer_etal_2017/{folder}/",
                                                            headers=UA), timeout=300).read().decode("utf-8", "replace")
    names = re.findall(r'<a href="([^"?/][^"]*\.bam)">', listing)
    want = next((n for n in names if re.search(rf"(^|\D){chrom}(\D|$)", n)), None)
    if want is None:
        gate(f"{ARCHAIC}: a chromosome-{chrom} alignment is published", False, f"found {names[:6]}")
        return
    url = f"{EVA}/bam/Pruefer_etal_2017/{folder}/{want}"
    info(f"streaming {REGION} from {want}")
    local = WORK / f"{ARCHAIC}_{chrom}_{start}_{end}.bam"
    if not local.exists():
        with open(local, "wb") as fh:
            rc, _, err = run(["samtools", "view", "-b", url, f"{chrom}:{start}-{end}"], stdout=fh)
        if rc:
            gate(f"{ARCHAIC}: the published alignment can be read remotely", False, err[-300:])
            return
        samtools("index", local)
    info(f"{local.name}: {local.stat().st_size / 1e6:.1f} MB")
    theirs = published_calls(vcf, f"{chrom}:{start}-{end}")
    info(f"{len(theirs):,} published genotypes in the region")
    pos = np.array(sorted(theirs), dtype=np.int64)
    maskf = fetch(mask, WORK / f"{ARCHAIC}_chr{chrom}_mask.bed.gz")
    callable_mask = C.bed_mask(str(maskf), chrom, start - 1, end)
    info(f"mask: {callable_mask.mean():.1%} of the region")
    with Genome(str(local), min_mapq=25) as g:
        pu = pileup(g, chrom, pos, max_depth=200)
        counts = pu.filtered(min_qual=25, trim=0).counts()
        dip = C.diploid(counts, min_depth=8, min_quality=30)
        region = C.region_counts(g, chrom, start - 1, end, min_qual=25)
    in_mask = callable_mask[pos - start]
    pass_only = np.array([theirs[int(p)][0] in ("PASS", ".", "") for p in pos])
    theirs_bases = [theirs[int(p)][1] for p in pos]
    ours = [tuple(sorted((BASES[int(a)], BASES[int(b)]))) if ok else None
            for a, b, ok in zip(dip["a1"], dip["a2"], dip["called"])]
    both = np.array([o is not None for o in ours]) & in_mask & pass_only
    agree = np.array([o == t for o, t, k in zip(ours, theirs_bases, both)])
    n = int(both.sum())
    share = float(agree[both].mean()) if n else 0.0
    het_theirs = np.array([t[0] != t[1] for t in theirs_bases])
    het_ours = np.array([o is not None and o[0] != o[1] for o in ours])
    report["archaic"] = {
        "individual": ARCHAIC, "region": REGION, "alignment": want, "vcf": vcf, "mask": mask,
        "bam_bytes": local.stat().st_size, "published_sites": len(theirs),
        "mask_fraction": float(callable_mask.mean()), "sites_compared": n,
        "genotypes_agree": share,
        "sites_we_call_but_they_filter": int((np.array([o is not None for o in ours]) & ~pass_only).sum()),
        "heterozygous_ours": int((het_ours & both).sum()), "heterozygous_theirs": int((het_theirs & both).sum()),
        "heterozygosity_region": C.heterozygosity(region, callable_mask=callable_mask, min_depth=8),
    }
    hp = int((het_theirs & both).sum()) / max(1, n)
    ho = int((het_ours & both).sum()) / max(1, n)
    print(json.dumps(report["archaic"], default=float)[:1500], flush=True)
    gate(f"{ARCHAIC} Neanderthal, {REGION}: diploid genotypes agree with the study's own call set at the sites both "
         "call, inside the study's filter mask (at least 99%)", n > 10_000 and share >= 0.99,
         f"{100 * share:.3f}% of {n:,} sites")
    gate(f"{ARCHAIC} Neanderthal, {REGION}: heterozygous sites per callable site within 25% of the published call "
         "set's", n > 10_000 and hp > 0 and abs(ho - hp) <= 0.25 * hp,
         f"ours {ho:.2e}, theirs {hp:.2e} over {n:,} shared sites")


def main():
    t0 = time.time()
    aadr = AADRRemote(WORK)
    report["aadr"] = {"sites": aadr.nsnp, "individuals": aadr.nind}
    step_1240k(aadr)
    step_archaic()
    report["results"] = results
    report["seconds"] = round(time.time() - t0)
    (WORK / "deeptime_report.json").write_text(json.dumps(report, indent=1, default=float))
    print(f"total {report['seconds']} s", flush=True)
    sys.exit(0 if all(r["pass"] for r in results) else 1)


if __name__ == "__main__":
    main()
