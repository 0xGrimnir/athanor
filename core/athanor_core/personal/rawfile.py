"""Consumer raw DNA files: AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA,
Living DNA, and single-sample VCF.

Formats (as the vendors document them):
  23andMe     tab-separated, '#' comments: rsid, chromosome, position, genotype
              ("AG", "--" no call, one letter on haploid X/Y/MT in males;
              chromosomes 1-22, X, Y, MT); GRCh37, + strand
  AncestryDNA tab-separated, '#' comments, header line: rsid, chromosome,
              position, allele1, allele2 ("0" no call; chromosomes 1-22,
              23 = X, 24 = Y, 25 = pseudoautosomal XY, 26 = MT); GRCh37, + strand
  MyHeritage  comma-separated, quoted, '#' comments, header RSID,CHROMOSOME,POSITION,RESULT
  FTDNA       the same CSV layout as MyHeritage (Family Finder), no comment block
  Living DNA  like 23andMe, with a Living DNA comment header
  VCF         one sample; GT field

Archives (.zip as downloaded from the vendors) and .gz files are opened
directly. Everything here runs on the user's computer.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

BASES = "ACGTID"            # I/D: 23andMe insertion/deletion calls
CODE = {b: i for i, b in enumerate(BASES)}
MISSING = 255
CHROM = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26,
         "23": 23, "24": 24, "25": 25, "26": 26, "0": 0}
CHROM_NAME = {**{i: str(i) for i in range(1, 23)}, 23: "X", 24: "Y", 25: "XY", 26: "MT", 0: "0"}


@dataclass
class RawGenotypes:
    vendor: str
    file_name: str
    sha256: str
    rsid: np.ndarray          # int64; -1 when the marker has no rs number
    chrom: np.ndarray         # int8: 1-22, 23 X, 24 Y, 25 XY, 26 MT
    pos: np.ndarray           # int64, as given in the file
    a1: np.ndarray            # uint8 base code, 255 = no call
    a2: np.ndarray            # uint8; equals a1 for haploid calls
    header: list = field(default_factory=list)
    build_declared: str | None = None

    def __len__(self) -> int:
        return len(self.pos)


def _open_text(path: Path) -> tuple[str, str]:
    """Return (text, inner name) for plain, .gz and .zip files."""
    data = path.read_bytes()
    name = path.name
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
        name = re.sub(r"\.gz$", "", name)
    elif data[:4] == b"PK\x03\x04":
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            members = [m for m in z.namelist() if not m.endswith("/") and not m.startswith("__MACOSX")]
            members.sort(key=lambda m: (not re.search(r"\.(txt|csv|tsv|vcf)$", m, re.I), -z.getinfo(m).file_size))
            if not members:
                raise ValueError("the archive is empty")
            name = members[0]
            data = z.read(name)
            if data[:2] == b"\x1f\x8b":
                data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace"), name


def _vendor(header_text: str, first_data: str, is_vcf: bool) -> str:
    h = header_text.lower()
    if is_vcf:
        return "VCF"
    if "23andme" in h:
        return "23andMe"
    if "ancestrydna" in h or "ancestry.com" in h:
        return "AncestryDNA"
    if "myheritage" in h:
        return "MyHeritage"
    if "living dna" in h or "livingdna" in h:
        return "Living DNA"
    if first_data.upper().startswith('RSID,') or first_data.upper().startswith('"RSID"'):
        return "FamilyTreeDNA"
    return "unknown"


def _code(b: str) -> int:
    return CODE.get(b, MISSING)


def parse(path: str | Path) -> RawGenotypes:
    path = Path(path)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    text, inner = _open_text(path)
    lines = text.splitlines()
    header = [l for l in lines[:400] if l.startswith("#")]
    body_start = next((i for i, l in enumerate(lines) if l.strip() and not l.startswith("#")), len(lines))
    is_vcf = any(l.startswith("##fileformat=VCF") for l in header)
    vendor = _vendor("\n".join(header), lines[body_start] if body_start < len(lines) else "", is_vcf)
    build = None
    for l in header:
        m = re.search(r"(?:build|grch|hg)\s*[-:]?\s*(36|37|38|19|18)\b", l, re.I)
        if m:
            build = {"19": "37", "18": "36"}.get(m.group(1), m.group(1))
            break
    if is_vcf:
        return _parse_vcf(lines, vendor, path.name, sha, header)

    rs, ch, ps, g1, g2 = [], [], [], [], []
    csv_mode = "," in lines[body_start] if body_start < len(lines) else False
    for l in lines[body_start:]:
        if not l or l.startswith("#"):
            continue
        parts = [p.strip().strip('"') for p in (l.split(",") if csv_mode else l.split("\t"))]
        if len(parts) < 4:
            parts = l.split()
            if len(parts) < 4:
                continue
        if parts[0].lower() in ("rsid", "snp", "name"):
            continue              # column header
        c = CHROM.get(parts[1].upper().replace("CHR", ""))
        if c is None:
            continue
        try:
            p = int(parts[2])
        except ValueError:
            continue
        if len(parts) >= 5 and vendor == "AncestryDNA":
            x, y = parts[3], parts[4]
        else:
            g = parts[3].upper()
            if len(g) == 2:
                x, y = g[0], g[1]
            elif len(g) == 1:
                x = y = g
            else:
                x = y = "-"
        rid = parts[0].lower()
        rs.append(int(rid[2:]) if rid.startswith("rs") and rid[2:].isdigit() else -1)
        ch.append(c); ps.append(p); g1.append(_code(x.upper())); g2.append(_code(y.upper()))
    if not ps:
        raise ValueError("no genotype lines recognised; expected a raw data file from a DNA testing company or a VCF")
    return RawGenotypes(vendor=vendor, file_name=inner, sha256=sha, rsid=np.array(rs, dtype=np.int64),
                        chrom=np.array(ch, dtype=np.int8), pos=np.array(ps, dtype=np.int64),
                        a1=np.array(g1, dtype=np.uint8), a2=np.array(g2, dtype=np.uint8),
                        header=header[:40], build_declared=build)


def _parse_vcf(lines, vendor, name, sha, header) -> RawGenotypes:
    rs, ch, ps, g1, g2 = [], [], [], [], []
    samples = 0
    for l in lines:
        if l.startswith("##"):
            continue
        if l.startswith("#CHROM"):
            samples = len(l.split("\t")) - 9
            if samples != 1:
                raise ValueError(f"the VCF has {samples} samples; Personal Genome reads one person's VCF")
            continue
        f = l.split("\t")
        if len(f) < 10:
            continue
        c = CHROM.get(f[0].upper().replace("CHR", ""))
        if c is None:
            continue
        ref, alts = f[3].upper(), f[4].upper().split(",")
        if len(ref) != 1 or any(len(a) != 1 for a in alts if a not in (".", "*")):
            continue      # SNVs only
        keys = f[8].split(":")
        if "GT" not in keys:
            continue
        gt = f[9].split(":")[keys.index("GT")]
        al = [ref] + alts
        idx = re.split(r"[/|]", gt)
        try:
            calls = [al[int(i)] for i in idx if i != "."]
        except (ValueError, IndexError):
            calls = []
        if not calls:
            x = y = "-"
        elif len(calls) == 1:
            x = y = calls[0]
        else:
            x, y = calls[0], calls[1]
        rid = f[2].lower()
        rs.append(int(rid[2:]) if rid.startswith("rs") and rid[2:].isdigit() else -1)
        ch.append(c); ps.append(int(f[1])); g1.append(_code(x)); g2.append(_code(y))
    build = None
    for l in header:
        m = re.search(r"(GRCh3[78]|hg19|hg38|b37)", l, re.I)
        if m:
            build = {"grch37": "37", "hg19": "37", "b37": "37", "grch38": "38", "hg38": "38"}[m.group(1).lower()]
            break
    return RawGenotypes(vendor=vendor, file_name=name, sha256=sha, rsid=np.array(rs, dtype=np.int64),
                        chrom=np.array(ch, dtype=np.int8), pos=np.array(ps, dtype=np.int64),
                        a1=np.array(g1, dtype=np.uint8), a2=np.array(g2, dtype=np.uint8),
                        header=[h for h in header if not h.startswith("##contig")][:40], build_declared=build)


def file_check(g: RawGenotypes) -> dict:
    """Markers, call rate per chromosome group, and genetic sex from the X and
    Y chromosomes (males carry one X: almost no heterozygous X calls outside
    the pseudoautosomal regions, and they have Y calls)."""
    called = g.a1 != MISSING
    auto = (g.chrom >= 1) & (g.chrom <= 22)
    x = g.chrom == 23
    y = g.chrom == 24
    mt = g.chrom == 26
    xc = x & called
    het_x = float((g.a1[xc] != g.a2[xc]).mean()) if xc.any() else None
    y_rate = float(called[y].mean()) if y.any() else None
    if het_x is None:
        sex, why = "unknown", "the file has no X-chromosome markers"
    elif het_x < 0.02 and (y_rate is None or y_rate > 0.5):
        sex, why = "XY", f"{100 * het_x:.2f}% of X markers heterozygous; Y markers called"
    elif het_x > 0.05 and (y_rate is None or y_rate < 0.2):
        sex, why = "XX", f"{100 * het_x:.1f}% of X markers heterozygous; Y markers not called"
    else:
        sex, why = "unclear", f"X heterozygosity {100 * het_x:.1f}%, Y call rate {100 * (y_rate or 0):.0f}%"
    het_auto = float((g.a1[auto & called] != g.a2[auto & called]).mean()) if (auto & called).any() else None
    return {"vendor": g.vendor, "file": g.file_name, "sha256": g.sha256, "markers": int(len(g)),
            "call_rate": float(called.mean()) if len(g) else 0.0,
            "autosomal_markers": int(auto.sum()), "autosomal_call_rate": float(called[auto].mean()) if auto.any() else None,
            "x_markers": int(x.sum()), "y_markers": int(y.sum()), "mt_markers": int(mt.sum()),
            "y_call_rate": y_rate, "x_heterozygosity": het_x, "autosomal_heterozygosity": het_auto,
            "genetic_sex": sex, "genetic_sex_basis": why, "build_declared": g.build_declared,
            "header": [h[:200] for h in g.header[:12]]}
