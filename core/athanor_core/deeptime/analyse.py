"""An ancient genome from aligned reads: authentication, sex, genotypes, and
the deep-ancestry analyses the Personal Genome tier already does.

The genotypes are pseudo-haploid by default, as in the published panels of
ancient individuals, so the ancestry models, principal components, archaic
ancestry and haplogroup calls are the same methods run on the same reference
data as for a present-day genotype file.
"""
from __future__ import annotations

import hashlib
import os
import time

import numpy as np

from ..personal import deep as pdeep
from ..personal.rawfile import RawGenotypes
from ..personal.refdata import RefData
from . import calls as C
from . import damage as D
from .bam import Genome, length_summary, sex_from_panel
from .sites import pileup

SOURCES = {
    "briggs2007": {"title": "Briggs et al., Patterns of damage in genomic DNA sequences from a Neandertal, "
                            "PNAS 104, 14616 (2007)", "url": "https://doi.org/10.1073/pnas.0704665104"},
    "skoglund2013": {"title": "Skoglund et al., Accurate sex identification of ancient human remains using DNA "
                              "shotgun sequencing, J. Archaeol. Sci. 40, 4477 (2013)",
                     "url": "https://doi.org/10.1016/j.jas.2013.07.004"},
    "mittnik2016": {"title": "Mittnik et al., A molecular approach to the sexing of the triple burial at the Upper "
                             "Paleolithic site of Dolni Vestonice, PLoS ONE 11, e0163019 (2016)",
                    "url": "https://doi.org/10.1371/journal.pone.0163019"},
    "skoglund2014": {"title": "Skoglund et al., Genomic diversity and admixture differs for Stone-Age Scandinavian "
                              "foragers and farmers, Science 344, 747 (2014)",
                     "url": "https://doi.org/10.1126/science.1253448"},
    "li2011": {"title": "Li, A statistical framework for SNP calling, mutation discovery, association mapping and "
                        "population genetical parameter estimation from sequencing data, Bioinformatics 27, 2987 "
                        "(2011)", "url": "https://doi.org/10.1093/bioinformatics/btr509"},
    "prufer2018": {"title": "Pruefer, snpAD: an ancient DNA genotype caller, Bioinformatics 34, 4165 (2018)",
                   "url": "https://doi.org/10.1093/bioinformatics/bty507"},
    "prufer2017": {"title": "Pruefer et al., A high-coverage Neandertal genome from Vindija Cave in Croatia, "
                            "Science 358, 655 (2017)", "url": "https://doi.org/10.1126/science.aao1887"},
}

Y_CHROM, MT_CHROM = 24, 26


def _sha256(path: str, limit: int = 1 << 28) -> dict:
    """Hash of the file, or of its first 256 MB when it is larger (a whole
    ancient genome can be hundreds of gigabytes)."""
    h, n = hashlib.sha256(), 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
            n += len(chunk)
            if n >= limit:
                return {"sha256_first_256mb": h.hexdigest(), "bytes_hashed": n,
                        "bytes_total": os.path.getsize(path)}
    return {"sha256": h.hexdigest(), "bytes_total": n}


def _panel(ref: RefData) -> dict:
    s = ref.snps
    return {"chrom": s["chrom"], "pos": s["pos"], "ref": s["ref"], "alt": s["alt"]}


def genotype_panel(g: Genome, ref: RefData, mode: str = "pseudo_haploid", seed: int = 1, min_qual: int = 20,
                   trim: int = 0, transversions_only: bool = False, min_depth: int = 8,
                   on_step=None) -> dict:
    """Calls at every panel site, chromosome by chromosome."""
    p = _panel(ref)
    dos = np.full(len(p["pos"]), -1, dtype=np.int8)
    depth_sum = depth_n = 0
    acc = None
    stats = {"sites": int(len(p["pos"])), "sites_covered": 0, "sites_called": 0, "chromosomes": [], "method": ""}
    depth_by_chrom, sites_by_chrom = {}, {}
    chroms = [c for c in range(1, 24) if (p["chrom"] == c).any()]
    for k, c in enumerate(chroms):
        idx = np.nonzero(p["chrom"] == c)[0]
        idx = idx[np.argsort(p["pos"][idx])]
        name = "X" if c == 23 else str(c)
        if g.name_of(name) is None:
            continue
        if on_step:
            on_step({"task": "deeptime", "stage": f"reading chromosome {name}",
                     "fraction": 0.1 + 0.5 * k / max(1, len(chroms))})
        pu = pileup(g, name, p["pos"][idx].astype(np.int64))
        if mode == "diploid":
            r = C.panel_genotypes(pu, p["ref"][idx], p["alt"][idx], min_qual=min_qual, trim=trim,
                                  min_depth=min_depth)
        else:
            r = C.pseudo_haploid(pu, p["ref"][idx], p["alt"][idx], seed=seed + c, min_qual=min_qual, trim=trim,
                                 transversions_only=transversions_only)
        dos[idx] = r["dosage"]
        covered = int(r.get("sites_covered", 0))
        stats["sites_covered"] += covered
        stats["sites_called"] += int((r["dosage"] >= 0).sum())
        stats["method"] = r["method"]
        if r.get("mean_depth_at_covered"):
            depth_sum += r["mean_depth_at_covered"] * covered
            depth_n += covered
        observations = int(len(pu.filtered(min_qual=min_qual, trim=trim).site))
        depth_by_chrom[name] = observations
        sites_by_chrom[name] = int(len(idx))
        stats["chromosomes"].append({"chromosome": name, "sites": int(len(idx)), "covered": covered,
                                     "called": int((r["dosage"] >= 0).sum()), "observations": observations})
        if c <= 22:
            acc = D.accumulate_panel(pu, p["ref"][idx], p["alt"][idx], acc, min_qual=min_qual)
    stats["mean_depth_at_covered_sites"] = (depth_sum / depth_n) if depth_n else 0.0
    return {"dosage": dos, "stats": stats, "depth_by_chrom": depth_by_chrom, "sites_by_chrom": sites_by_chrom,
            "damage": D.finish_panel(acc) if acc else {"terminal_deamination": None, "method": D.PANEL_METHOD}}


def _raw_from_haploid(ycalls: dict, mtcalls: dict) -> RawGenotypes:
    code = {"A": 0, "C": 1, "G": 2, "T": 3}
    pos, chrom, a = [], [], []
    for p, b in sorted(ycalls.items()):
        pos.append(p); chrom.append(Y_CHROM); a.append(code[b])
    for p, b in sorted(mtcalls.items()):
        pos.append(p); chrom.append(MT_CHROM); a.append(code[b])
    arr = np.array(a, dtype=np.uint8)
    return RawGenotypes(vendor="aligned reads", file_name="", sha256="", rsid=np.full(len(pos), -1, dtype=np.int64),
                        chrom=np.array(chrom, dtype=np.int8), pos=np.array(pos, dtype=np.int64),
                        a1=arr, a2=arr.copy())


def analyse(path: str, ref: RefData | None = None, mode: str = "pseudo_haploid", seed: int = 1,
            min_mapq: int = 25, min_qual: int = 20, trim: int = 0, transversions_only: bool = False,
            min_depth: int = 8, max_reads_for_lengths: int = 500_000, reference_fasta: str | None = None,
            damage_chromosomes: list[str] | None = None, damage_max_reads: int = 300_000, on_step=None) -> dict:
    ref = ref or RefData()
    t0 = time.time()

    def step(msg, f):
        if on_step:
            on_step({"task": "deeptime", "stage": msg, "fraction": f, "elapsed_seconds": round(time.time() - t0)})

    step("opening the alignment", 0.02)
    with Genome(path, min_mapq=min_mapq) as g:
        if not g.indexed:
            raise ValueError("the BAM has no index (.bai); index it with samtools index")
        if g.build != "GRCh37":
            raise ValueError(f"the reads are mapped to {g.build}; this step needs GRCh37 (hs37d5, hg19 or "
                             "1000 Genomes naming), the build the reference panel is on")
        counts = g.counts()
        sex = g.sex()
        step("molecule lengths", 0.06)
        lp = g.lengths_and_damage_pass(max_reads=max_reads_for_lengths)
        lens = length_summary(lp["length_histogram"])
        step("reading the panel's positions", 0.1)
        gt = genotype_panel(g, ref, mode=mode, seed=seed, min_qual=min_qual, trim=trim,
                            transversions_only=transversions_only, min_depth=min_depth, on_step=on_step)
        panel_sex = sex_from_panel(gt["depth_by_chrom"], gt["sites_by_chrom"])
        sex = {**panel_sex, "whole_chromosome_counts": sex} if panel_sex.get("call") else \
              {**panel_sex, "whole_chromosome_counts": sex,
               "call": sex.get("call"), "fell_back": sex.get("call") is not None,
               "method": panel_sex.get("method", "") + (
                   " With too little coverage at the panel's X positions, the call falls back to the read counts "
                   "over whole chromosomes (R_y), which hold for shotgun sequencing." if sex.get("call") else "")}
        if reference_fasta:
            step("deamination against the reference", 0.58)
            try:
                gt["damage"] = D.from_reference(g, reference_fasta, chroms=damage_chromosomes,
                                                max_reads=damage_max_reads, min_qual=min_qual, on_step=on_step)
                gt["damage"]["source"] = "reference"
            except Exception as e:                                        # noqa: BLE001
                gt["damage"]["reference_error"] = str(e)
        gt["damage"].setdefault("source", "panel")
        step("maternal and paternal lines", 0.62)
        ytree = ref.tree("y")
        mt_ok = g.mt_usable
        ycalls, ystats = ({}, {"positions": 0, "called": 0, "note": "no Y chromosome in the alignment"})
        if g.name_of("Y") is not None and sex.get("call") != "XX":
            from ..personal.haplo import Tree
            ycalls, ystats = C.haploid_calls(g, "Y", np.array(sorted(Tree(ytree).positions), dtype=np.int64),
                                             min_qual=min_qual, trim=trim)
        mtcalls, mtstats = ({}, {"positions": 0, "called": 0,
                                 "note": None if mt_ok else "the alignment's mitochondrial reference is not the rCRS, "
                                                             "so its coordinates do not match the mtDNA tree"})
        if mt_ok and g.name_of("MT") is not None:
            from ..personal.haplo import Tree
            mtcalls, mtstats = C.haploid_calls(g, "MT", np.array(sorted(Tree(ref.tree("mt")).positions),
                                                                 dtype=np.int64), min_qual=min_qual, trim=trim)
        file_info = {"file": os.path.basename(path), **_sha256(path), "build": g.build,
                     "build_detail": g.build_detail, "read_groups": g.read_groups(),
                     "programs": g.programs(), "reads_mapped": int(sum(counts.values())),
                     "reads_by_chromosome": counts}

    dos = gt["dosage"]
    s = ref.snps
    auto = (s["chrom"] >= 1) & (s["chrom"] <= 22) & (dos >= 0)
    dmg = gt["damage"].get("terminal_deamination")
    cards = [pdeep._card("deeptime_file", "The alignment", "The sample", "A",
                         f"{file_info['reads_mapped']:,} mapped reads on {g.build}; "
                         f"{gt['stats']['sites_called']:,} of {gt['stats']['sites']:,} panel sites called "
                         f"({int(auto.sum()):,} autosomal, usable for ancestry); median molecule length "
                         f"{lens.get('median', 0):.0f} bp",
                         {**file_info, "panel": gt["stats"], "lengths": lens, "length_pass": lp},
                         ["A BAM says which reference it was mapped to, not how it was mapped; read groups and "
                          "programs are reported as the file states them.",
                          "Molecule lengths are read lengths, which are the molecules' lengths only when adapters "
                          "were trimmed and overlapping pairs merged, as ancient-DNA pipelines do."],
                         "Chromosome names and lengths against GRCh37 and GRCh38; read counts from the BAM index; "
                         "lengths from a pass over the autosomes.", [SOURCES["skoglund2014"], "aadr"], ref),
             pdeep._card("deeptime_authenticity", "Is this DNA ancient?", "The sample", "A",
                         (f"terminal deamination {100 * dmg:.1f}% outside CpG sites "
                          f"({100 * (gt['damage'].get('terminal_deamination_including_cpg') or dmg):.1f}% with them): "
                          f"C→T at the 5′ end and G→A at the 3′ end above background; median molecule "
                          f"{lens.get('median', 0):.0f} bp"
                          if dmg is not None else "not enough covered sites to measure deamination"),
                         {"damage": gt["damage"], "lengths": lens,
                          "reverse_fraction": lp["reverse_fraction"]},
                         ([] if gt["damage"].get("source") == "reference" else
                          ["Without the reference sequence the rate can only be measured at the panel's own "
                           "positions, where there are few observations and CpG sites, which deaminate several "
                           "times faster, are over-represented; the non-CpG figure is the comparable one. Choose "
                           "the reference FASTA the reads were mapped to for the full measurement."]) +
                         ["A library treated with uracil-DNA glycosylase shows little or no damage even when the "
                          "DNA is genuinely ancient.",
                          "Damage and short molecules show that DNA is old and degraded. They cannot show whose it "
                          "is: present-day contamination is measured separately."],
                         gt["damage"]["method"], [SOURCES["briggs2007"], SOURCES["skoglund2014"]], ref,
                         ok=dmg is not None),
             pdeep._card("deeptime_sex", "Genetic sex", "The sample", "A",
                         ((f"{sex['call']} (R_x = {sex['r_x']:.2f}, interval {sex['r_x_ci'][0]:.2f}–"
                           f"{sex['r_x_ci'][1]:.2f}, at the panel's positions)" if sex.get("r_x") is not None
                           else f"{sex['call']} (from read counts over whole chromosomes)") if sex.get("call")
                          else f"not called: {sex.get('reason')}"),
                         sex, ["Genetic sex here means which sex chromosomes the reads come from.",
                               "Very low coverage, or contamination by DNA of the other sex, widens the interval "
                               "until neither range fits."],
                         sex.get("method", ""), [SOURCES["skoglund2013"], SOURCES["mittnik2016"]], ref,
                         ok=sex.get("call") is not None)]
    if int(auto.sum()) < 10000:
        return {"version": 1, "seconds": round(time.time() - t0, 1), "cards": cards,
                "reference_pack": {k: ref.manifest.get(k) for k in ("name", "built")},
                "note": f"only {int(auto.sum())} autosomal panel sites were called; too few for the ancestry "
                        "analyses, which need at least 10,000"}
    step("ancient ancestry models (qpAdm)", 0.68)
    cards.append(pdeep.ancient_ancestry(ref, dos, auto))
    step("principal components and closest ancient individuals", 0.78)
    cards.extend(pdeep.pca_cards(ref, dos))
    step("Neanderthal and Denisovan ancestry", 0.86)
    cards.extend(pdeep.archaic_cards(ref, dos, auto))
    step("haplogroups", 0.92)
    raw = _raw_from_haploid(ycalls, mtcalls)
    cards.extend(pdeep.haplogroup_cards(ref, raw, raw.pos, sex.get("call") or "unknown"))
    for c in cards:
        if c["id"] in ("y_haplogroup", "mt_haplogroup"):
            c["result"] = {**c["result"], "reads": ystats if c["id"] == "y_haplogroup" else mtstats}
    step("done", 1.0)
    return {"version": 1, "file": file_info, "seconds": round(time.time() - t0, 1),
            "reference_pack": {k: ref.manifest.get(k) for k in ("name", "built", "sources")},
            "genotypes": {"mode": mode, "seed": seed, **gt["stats"]}, "cards": cards}
