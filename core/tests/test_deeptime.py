"""Deep Time: an ancient genome simulated from a known reference, with known
genotypes and known deamination, read back through the pipeline."""
import numpy as np
import pytest

pysam = pytest.importorskip("pysam")

from athanor_core.deeptime import calls as C                      # noqa: E402
from athanor_core.deeptime import damage as D                     # noqa: E402
from athanor_core.deeptime.bam import (GRCH37, Genome, length_summary, plain,    # noqa: E402
                                       sex_from_panel, wilson)
from athanor_core.deeptime.sites import pileup    # noqa: E402

BASES = "ACGT"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
CHROM_LEN = 1_000_000
DAMAGE0, LAM = 0.30, 0.3


def _reference(rng, n=CHROM_LEN):
    return "".join(rng.choice(list(BASES), n))


def _sites(rng, ref_seq, n=8000, spacing=120):
    """Panel sites every ``spacing`` bases: reference allele from the sequence,
    an alternative allele drawn from the other three bases."""
    pos = np.arange(1, n + 1) * spacing
    ref = np.array([BASES.index(ref_seq[p - 1]) for p in pos], dtype=np.int8)
    alt = np.array([(r + 1 + rng.integers(0, 3)) % 4 for r in ref], dtype=np.int8)
    return pos.astype(np.int64), ref, alt


def _simulate(tmp_path, rng, ref_seq, pos, ref, alt, depth=16, read_len=60, damage=DAMAGE0, with_y=0,
              cpg_only=False):
    """Reads copied from the reference, carrying each site's genotype, damaged
    as ancient DNA is: C->T towards the molecule's 5' end and G->A towards its
    3' end, on whichever strand the molecule was sequenced from."""
    truth = rng.integers(0, 3, len(pos))          # copies of the alternative allele
    genome = [list(ref_seq), list(ref_seq)]       # two chromosome copies
    for p, r, a, g in zip(pos, ref, alt, truth):
        genome[0][p - 1] = BASES[a] if g >= 1 else BASES[r]
        genome[1][p - 1] = BASES[a] if g == 2 else BASES[r]
    reads = []
    for p in pos:
        for _ in range(int(rng.poisson(depth))):
            hap = genome[int(rng.integers(0, 2))]
            off = int(rng.integers(0, read_len))
            start = max(0, min(CHROM_LEN - read_len, p - 1 - off))
            seq = list(hap[start:start + read_len])
            rev = bool(rng.integers(0, 2))
            original = list(seq)
            for j in range(len(seq)):
                d5 = (len(seq) - 1 - j) if rev else j
                d3 = j if rev else (len(seq) - 1 - j)
                frm, to = ("G", "A") if rev else ("C", "T")
                frm3, to3 = ("C", "T") if rev else ("G", "A")
                # a cytosine on the molecule is followed by a guanine; stored
                # reverse-complemented that is a G preceded by a C
                cpg5 = (j > 0 and original[j - 1] == "C") if rev else (j + 1 < len(seq) and original[j + 1] == "G")
                cpg3 = (j + 1 < len(seq) and original[j + 1] == "G") if rev else (j > 0 and original[j - 1] == "C")
                if seq[j] == frm and (cpg5 or not cpg_only) and rng.random() < damage * (1 - LAM) ** d5:
                    seq[j] = to
                elif seq[j] == frm3 and (cpg3 or not cpg_only) and rng.random() < damage * (1 - LAM) ** d3:
                    seq[j] = to3
            reads.append(("1", start, "".join(seq), rev))
    for _ in range(with_y):
        reads.append(("Y", int(rng.integers(0, 50_000)), _reference(rng, read_len), False))
    path = str(tmp_path / "ancient.bam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "1", "LN": CHROM_LEN}, {"SN": "X", "LN": 60_000}, {"SN": "Y", "LN": 50_000}],
              "RG": [{"ID": "L1", "SM": "ANC1", "LB": "lib1", "PL": "ILLUMINA"}]}
    order = {"1": 0, "X": 1, "Y": 2}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for n, (chrom, start, seq, rev) in enumerate(sorted(reads, key=lambda r: (order[r[0]], r[1]))):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence = f"r{n}", seq
            a.reference_id, a.reference_start = order[chrom], start
            a.flag = 16 if rev else 0
            a.mapping_quality = 37
            a.cigarstring = f"{len(seq)}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.set_tag("RG", "L1")
            out.write(a)
    pysam.index(path)
    return path, truth


@pytest.fixture(scope="module")
def ancient(tmp_path_factory):
    rng = np.random.default_rng(11)
    ref_seq = _reference(rng)
    pos, ref, alt = _sites(rng, ref_seq)
    d = tmp_path_factory.mktemp("deeptime")
    path, truth = _simulate(d, rng, ref_seq, pos, ref, alt, with_y=3000)
    fasta = d / "ref.fa"
    fasta.write_text(">1\n" + "\n".join(ref_seq[i:i + 60] for i in range(0, len(ref_seq), 60)) + "\n")
    pysam.faidx(str(fasta))
    return {"path": path, "pos": pos, "ref": ref, "alt": alt, "truth": truth, "sequence": ref_seq,
            "fasta": str(fasta)}


def test_damage_against_the_reference_sequence(ancient):
    """With the reference at hand every aligned base counts, so the profile is
    the whole curve, split by dinucleotide context."""
    with Genome(ancient["path"]) as g:
        d = D.from_reference(g, ancient["fasta"], chroms=["1"])
    model = DAMAGE0 * (1 - LAM) ** np.arange(10)
    for key in ("ct5", "ga3"):
        prof = d[key]
        assert prof["observations"][0] > 20_000           # far more than the panel's positions give
        assert prof["background"] < 0.01
        assert abs(prof["terminal_excess"] - DAMAGE0) < 0.02
        obs = np.array(prof["frequency"][:10], dtype=float) - prof["background"]
        assert np.abs(obs - model).mean() < 0.01
        # the simulation damages every cytosine alike, so the two contexts must agree
        for context in ("cpg", "not_cpg"):
            assert abs(prof[context]["terminal_excess"] - DAMAGE0) < 0.05, (key, context)
        assert prof["cpg"]["observations"][0] + prof["not_cpg"]["observations"][0] == prof["observations"][0]
    assert abs(d["terminal_deamination"] - DAMAGE0) < 0.02
    assert d["reads_used"] > 1000 and d["chromosomes"] == ["1"]


def test_names_and_intervals():
    assert plain("chr1") == "1" and plain("chrM") == "MT" and plain("MT") == "MT" and plain("chrX") == "X"
    lo, hi = wilson(0, 1000)
    assert lo == 0 and hi < 0.01
    lo, hi = wilson(500, 1000)
    assert lo < 0.5 < hi


def test_build_is_read_off_the_chromosome_lengths(tmp_path):
    path = str(tmp_path / "b.bam")
    sq = [{"SN": c, "LN": L} for c, L in GRCH37.items()]
    with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}) as out:
        pass
    pysam.index(path)
    with Genome(path) as g:
        assert g.build == "GRCh37"
        assert g.build_detail["naming"] == "plain"
        assert g.mt_usable is False            # no MT in this header


def test_damage_is_measured_at_the_simulated_rate(ancient):
    """Without a reference sequence the panel's own positions are the only
    places a reference base is known; the rate must still come back right."""
    with Genome(ancient["path"]) as g:
        p = pileup(g, "1", ancient["pos"])
        d = D.from_panel(p, ancient["ref"], ancient["alt"])
    model = DAMAGE0 * (1 - LAM) ** np.arange(10)
    for key in ("ct5", "ga3"):
        prof = d[key]
        assert prof["observations"][0] > 100
        assert prof["background"] < 0.02
        assert abs(prof["terminal_excess"] - DAMAGE0) < 0.08
        obs = np.array(prof["frequency"][:10], dtype=float) - prof["background"]
        assert np.abs(obs - model).mean() < 0.03          # the whole profile, not just the first base
        # the profile falls off towards the middle of the molecule
        assert prof["frequency"][0] > prof["frequency"][3] > prof["frequency"][10]
    assert abs(d["terminal_deamination"] - DAMAGE0) < 0.08


def test_only_cpg_cytosines_damaged(tmp_path_factory):
    """Damage only the cytosines that are followed by a guanine, on both
    strands. The CpG bucket must carry all of it and the other none, which a
    swapped neighbour or a reference-versus-molecule mix-up could not do."""
    rng = np.random.default_rng(23)
    ref_seq = _reference(rng)
    pos, ref, alt = _sites(rng, ref_seq, n=2500, spacing=120)
    d = tmp_path_factory.mktemp("cpg")
    path, _ = _simulate(d, rng, ref_seq, pos, ref, alt, depth=20, damage=0.8, cpg_only=True)
    fasta = d / "ref.fa"
    fasta.write_text(">1\n" + "\n".join(ref_seq[i:i + 60] for i in range(0, len(ref_seq), 60)) + "\n")
    pysam.faidx(str(fasta))
    with Genome(path) as g:
        d_ref = D.from_reference(g, str(fasta), chroms=["1"])
    for key in ("ct5", "ga3"):
        prof = d_ref[key]
        assert prof["cpg"]["observations"][0] > 200 and prof["not_cpg"]["observations"][0] > 2000
        assert prof["cpg"]["terminal_excess"] > 0.5, (key, prof["cpg"]["terminal_excess"])
        assert abs(prof["not_cpg"]["terminal_excess"]) < 0.02, (key, prof["not_cpg"]["terminal_excess"])


def test_pseudo_haploid_calls_are_unbiased(ancient):
    with Genome(ancient["path"]) as g:
        p = pileup(g, "1", ancient["pos"])
    r = C.pseudo_haploid(p, ancient["ref"], ancient["alt"], seed=5, trim=5)
    dos, truth = r["dosage"], ancient["truth"]
    called = dos >= 0
    assert called.mean() > 0.9
    hom = called & (truth != 1)
    assert (dos[hom] == truth[hom]).mean() > 0.98
    het = called & (truth == 1)
    share = (dos[het] == 2).mean()
    assert 0.42 < share < 0.58                      # one allele at random from a heterozygote
    same = C.pseudo_haploid(p, ancient["ref"], ancient["alt"], seed=5, trim=5)
    assert np.array_equal(dos, same["dosage"])      # seeded: the same calls again


def test_diploid_genotypes_match_the_truth_when_ends_are_trimmed(ancient):
    with Genome(ancient["path"]) as g:
        p = pileup(g, "1", ancient["pos"])
    raw = C.panel_genotypes(p, ancient["ref"], ancient["alt"], min_depth=6)
    trimmed = C.panel_genotypes(p, ancient["ref"], ancient["alt"], min_depth=6, trim=6)
    truth = ancient["truth"]
    for r, floor in ((raw, 0.90), (trimmed, 0.98)):
        ok = r["dosage"] >= 0
        assert ok.mean() > 0.8
        assert (r["dosage"][ok] == truth[ok]).mean() > floor
    rawok = raw["dosage"] >= 0
    trok = trimmed["dosage"] >= 0
    assert ((trimmed["dosage"][trok] == truth[trok]).mean()
            > (raw["dosage"][rawok] == truth[rawok]).mean())       # trimming removes damage errors


def test_transversions_only_drops_the_sites_damage_can_reach(ancient):
    with Genome(ancient["path"]) as g:
        p = pileup(g, "1", ancient["pos"])
    all_sites = C.pseudo_haploid(p, ancient["ref"], ancient["alt"], seed=1)
    tv = C.pseudo_haploid(p, ancient["ref"], ancient["alt"], seed=1, transversions_only=True)
    assert 0.5 < tv["sites_called"] / all_sites["sites_called"] < 0.8
    truth = ancient["truth"]
    ok = tv["dosage"] >= 0
    hom = ok & (truth != 1)
    assert (tv["dosage"][hom] == truth[hom]).mean() > 0.99


def test_heterozygosity_recovers_the_simulated_rate(ancient):
    """Only the panel's heterozygous sites differ between the two chromosome
    copies, so the heterozygosity of the region is known. Deamination mimics
    heterozygous transitions, so the transversion rate is the one that must
    come back right, and the total must come back too high."""
    pos, truth = ancient["pos"], ancient["truth"]
    start, end = 1, int(pos[-1])
    with Genome(ancient["path"]) as g:
        r = C.heterozygosity_over_region(g, "1", start - 1, end, min_depth=8, min_quality=30)
        counts = C.region_counts(g, "1", start - 1, end)          # an independent count of the depths
    floor = r["depth_floor"]
    assert floor >= 10                                            # a homozygote must be callable too
    deep_positions = int((counts.sum(1) >= floor).sum())
    # the denominator is the positions deep enough, less the few the caller cannot settle
    # (damage leaves some sites genuinely ambiguous); it can never exceed them
    assert r["callable_sites"] <= deep_positions
    assert r["callable_sites"] > 0.9 * deep_positions
    deep = counts.sum(1)[pos - start] >= floor
    het = deep & (truth == 1)
    transition = np.array([{(r_, a_)} <= {(0, 2), (2, 0), (1, 3), (3, 1)}
                           for r_, a_ in zip(ancient["ref"], ancient["alt"])])
    expected_tv = int((het & ~transition).sum())
    assert abs(r["transversions"] - expected_tv) < 0.1 * expected_tv
    assert r["transversion_heterozygosity"] == pytest.approx(expected_tv / r["callable_sites"], rel=0.1)
    lo, hi = r["transversion_heterozygosity_ci"]
    assert lo < expected_tv / r["callable_sites"] < hi
    # damage turns homozygous C/G sites into apparent C/T and G/A heterozygotes
    assert r["transitions"] > 1.15 * int((het & transition).sum())
    assert r["heterozygous_sites"] > int(het.sum())


def test_heterozygosity_denominator_does_not_depend_on_the_genotype():
    """A quality floor is also a depth floor, and a stricter one for
    homozygotes. If the two are left inconsistent, homozygous sites fall out of
    the denominator while heterozygous ones stay in the numerator, and every
    rate computed from them comes out too high."""
    rng = np.random.default_rng(2)
    n = 60_000
    depth = rng.poisson(20, n)
    is_het = rng.random(n) < 2e-3
    counts = np.zeros((n, 4), dtype=np.int64)
    for i in range(n):
        if depth[i] == 0:
            continue
        if is_het[i]:
            k = rng.binomial(depth[i], 0.5)
            counts[i, 0], counts[i, 1] = k, depth[i] - k
        else:
            counts[i, 0] = depth[i]
    r = C.heterozygosity(counts, min_depth=8, min_quality=30)
    truth = is_het[depth >= r["depth_floor"]].mean()
    assert 0.9 < r["heterozygosity"] / truth < 1.1
    assert r["depth_floor"] == C.minimum_depth_for_a_homozygote(0.005, 30.0)
    # the floor rises with the quality asked for, and the rate stays right
    strict = C.heterozygosity(counts, min_depth=8, min_quality=60)
    assert strict["depth_floor"] > r["depth_floor"]
    assert 0.9 < strict["heterozygosity"] / is_het[depth >= strict["depth_floor"]].mean() < 1.1


def test_base_qualities_count_for_what_they_say():
    """A Q10 base is ten times as likely to be wrong as a Q20 one, so a pile of
    bad bases must not call a heterozygote the way good ones would."""
    from athanor_core.deeptime.sites import Pileup
    n_sites = 2

    def pile(quals):
        site = np.zeros(len(quals), dtype=np.int32)
        base = np.array([0] * (len(quals) // 2) + [1] * (len(quals) - len(quals) // 2), dtype=np.int8)
        z = np.zeros(len(quals), dtype=np.int16)
        return Pileup(site, base, np.array(quals, dtype=np.int16), z, z, np.zeros(len(quals), bool),
                      n_sites, len(quals), len(quals))

    good = C.diploid(pileup=pile([37] * 12), n_sites=n_sites, min_depth=8, min_quality=30)
    bad = C.diploid(pileup=pile([3] * 12), n_sites=n_sites, min_depth=8, min_quality=30)
    assert good["called"][0] and (good["a1"][0], good["a2"][0]) == (0, 1)
    assert good["quality"][0] > bad["quality"][0]
    assert not bad["called"][0]


def _sex_bam(path, rng, counts: dict, duplicates: int = 0, secondary: int = 0, mapq: int = 37):
    """A BAM with a GRCh37-shaped header and a chosen number of reads on each
    named chromosome, plus flagged reads that must not be counted."""
    sq = [{"SN": c, "LN": L} for c, L in GRCH37.items()]
    reads = []
    for chrom, n in counts.items():
        for _ in range(n):
            reads.append((chrom, int(rng.integers(1_000_000, 2_000_000)), 0, mapq))
    for chrom, n, flag in (("Y", duplicates, 1024), ("Y", secondary, 256)):
        for _ in range(n):
            reads.append((chrom, int(rng.integers(1_000_000, 2_000_000)), flag, mapq))
    order = {c["SN"]: i for i, c in enumerate(sq)}
    with pysam.AlignmentFile(str(path), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}) as out:
        for n, (chrom, start, flag, mq) in enumerate(sorted(reads, key=lambda r: (order[r[0]], r[1]))):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence = f"s{n}", "ACGT" * 15
            a.reference_id, a.reference_start, a.flag, a.mapping_quality = order[chrom], start, flag, mq
            a.cigarstring = "60M"
            a.query_qualities = pysam.qualitystring_to_array("I" * 60)
            out.write(a)
    pysam.index(str(path))
    return str(path)


@pytest.mark.parametrize("nx, ny, expect", [(20_000, 20, "XX"),      # R_y 0.001
                                            (10_000, 1_600, "XY"),   # R_y 0.14
                                            (10_000, 450, None)])    # R_y 0.043: between the thresholds
def test_sex_from_whole_chromosome_counts_needs_the_same_filters(tmp_path, nx, ny, expect):
    rng = np.random.default_rng(4)
    bam = _sex_bam(tmp_path / f"s{nx}_{ny}.bam", rng, {"1": 50_000, "2": 50_000, "X": nx, "Y": ny},
                   duplicates=5_000, secondary=5_000)
    with Genome(bam) as g:
        unfiltered = g.sex()
        filtered = g.sex(filtered=True)
    # the index counts the duplicates and secondary alignments, so no call is made from them
    assert unfiltered["call"] is None and "index" in unfiltered["reason"]
    assert unfiltered["reads_Y"] > filtered["reads_Y"]
    assert filtered["reads_Y"] == ny                  # flagged reads dropped
    assert filtered["call"] == expect, (filtered["r_y"], filtered["r_y_ci"])
    if expect is None:
        assert "spans the thresholds" in filtered["reason"]


def test_the_r_x_interval_carries_the_counting_noise():
    """The X's own count sits in every per-autosome ratio, so the spread
    between those ratios cannot measure how well it is known. At low coverage
    that noise is the whole uncertainty."""
    rng = np.random.default_rng(5)
    for per_site, truth in ((2.0, 0.5), (2.0, 1.0), (10.0, 0.5)):
        sites = {str(c): 1000 for c in range(1, 23)} | {"X": 1000}
        misses = 0
        for _ in range(200):
            depth = {c: int(rng.poisson(per_site * n)) for c, n in sites.items() if c != "X"}
            depth["X"] = int(rng.poisson(per_site * truth * sites["X"]))
            lo, hi = sex_from_panel(depth, sites)["r_x_ci"]
            misses += not (lo <= truth <= hi)
        assert misses / 200 < 0.12, (per_site, truth, misses)       # a 95% interval, not a decorative one
    # and the interval widens as coverage falls
    thin = sex_from_panel({str(c): 1000 for c in range(1, 23)} | {"X": 500}, {str(c): 1000 for c in range(1, 23)} | {"X": 1000})
    thick = sex_from_panel({str(c): 100_000 for c in range(1, 23)} | {"X": 50_000},
                           {str(c): 1000 for c in range(1, 23)} | {"X": 1000})
    assert (thin["r_x_ci"][1] - thin["r_x_ci"][0]) > 5 * (thick["r_x_ci"][1] - thick["r_x_ci"][0])


def test_molecule_lengths_and_strand_balance(ancient):
    with Genome(ancient["path"]) as g:
        lp = g.lengths_and_damage_pass(chroms=["1"], max_reads=20_000)
    summary = length_summary(lp["length_histogram"])
    assert summary["median"] == 60 and summary["mean"] == pytest.approx(60, abs=0.01)
    assert 0.45 < lp["reverse_fraction"] < 0.55


def test_bed_mask_marks_only_its_own_intervals(tmp_path):
    bed = tmp_path / "m.bed"
    bed.write_text("1\t100\t200\n1\t300\t310\n2\t0\t1000\n")
    m = C.bed_mask(str(bed), "1", 150, 320)
    assert m[:50].all()                 # 150-200
    assert not m[50:150].any()          # 200-300
    assert m[150:160].all()             # 300-310
    assert not m[160:].any()


def test_bed_mask_reads_the_files_people_actually_have(tmp_path):
    """Gzipped, chr-prefixed, with a track line, and with touching and
    overlapping intervals; and a contig whose name starts with t or b."""
    import gzip as _gzip
    plain_bed = tmp_path / "p.bed"
    plain_bed.write_text("track name=\"mask\"\nbrowser hide all\n# a comment\n"
                         "chr1\t100\t200\nchr1\t200\t250\nchr1\t240\t260\ntig1\t0\t1000\n")
    gz = tmp_path / "p.bed.gz"
    gz.write_bytes(_gzip.compress(plain_bed.read_bytes()))
    for path in (plain_bed, gz):
        m = C.bed_mask(str(path), "1", 90, 300)
        assert not m[:10].any()
        assert m[10:170].all()          # 100-260, the three intervals joined
        assert not m[170:].any()
    assert C.bed_mask(str(plain_bed), "tig1", 0, 10).all()      # not mistaken for a track line


# --- the driver, end to end -------------------------------------------------

def _grch37_bam(tmp_path, rng, panel, depth=12, sex="XY"):
    """A BAM with the real GRCh37 chromosome lengths, carrying reads at the
    panel's own positions, so analyse() can run as it does on real data."""
    sq = [{"SN": c, "LN": L} for c, L in GRCH37.items()] + [{"SN": "MT", "LN": 16569}]
    order = {c["SN"]: i for i, c in enumerate(sq)}
    reads, truth = [], {}
    for chrom, positions, ref, alt in panel:
        for p, r_, a_ in zip(positions, ref, alt):
            g = int(rng.integers(0, 3))
            truth[(chrom, int(p))] = g
            alleles = [a_] * g + [r_] * (2 - g)
            d = int(rng.poisson(depth if chrom != "X" or sex == "XX" else depth / 2))
            for _ in range(d):
                al = alleles[int(rng.integers(0, 2))]
                off = int(rng.integers(0, 50))
                seq = list(BASES[int(x)] for x in rng.integers(0, 4, 60))
                seq[off] = BASES[al]
                reads.append((chrom, max(0, p - 1 - off), "".join(seq), bool(rng.integers(0, 2))))
    if sex == "XY":
        for _ in range(4000):
            reads.append(("Y", int(rng.integers(2_700_000, 28_000_000)),
                          "".join(BASES[int(x)] for x in rng.integers(0, 4, 60)), False))
    path = str(tmp_path / f"genome_{sex}.bam")
    with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq,
                                                 "RG": [{"ID": "L1", "SM": "A1", "LB": "lib"}]}) as out:
        for n, (chrom, start, seq, rev) in enumerate(sorted(reads, key=lambda r: (order[r[0]], r[1]))):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence = f"r{n}", seq
            a.reference_id, a.reference_start = order[chrom], start
            a.flag = 16 if rev else 0
            a.mapping_quality = 37
            a.cigarstring = "60M"
            a.query_qualities = pysam.qualitystring_to_array("I" * 60)
            a.set_tag("RG", "L1")
            out.write(a)
    pysam.index(path)
    return path, truth


class _Panel:
    """The smallest stand-in for the Personal Genome reference pack: enough for
    analyse() to call genotypes and report, without the ancestry analyses."""

    def __init__(self, snps):
        self.snps = snps
        self.manifest = {"name": "test-pack", "built": "2026-10-25", "validation": {}}

    def tree(self, which):
        return {"tree": which, "nodes": [{"name": "root", "parent": -1, "muts": []}]}


def _small_panel(rng, per_chrom=400):
    chrom, pos, ref, alt = [], [], [], []
    for c in list(range(1, 23)) + [23]:
        p = np.sort(rng.choice(np.arange(1_000_000, 5_000_000), per_chrom, replace=False))
        r = rng.integers(0, 4, per_chrom)
        a = (r + 1 + rng.integers(0, 3, per_chrom)) % 4
        chrom += [c] * per_chrom; pos += list(p); ref += list(r); alt += list(a)
    return {"chrom": np.array(chrom, dtype=np.int8), "pos": np.array(pos, dtype=np.int64),
            "ref": np.array(ref, dtype=np.uint8), "alt": np.array(alt, dtype=np.uint8),
            "rsid": np.full(len(pos), -1, dtype=np.int64)}


@pytest.mark.parametrize("sex", ["XY", "XX"])
def test_analyse_runs_end_to_end_and_its_cards_say_what_was_done(tmp_path, sex):
    from athanor_core.deeptime.analyse import analyse
    rng = np.random.default_rng(13)
    snps = _small_panel(rng)
    panel = []
    for c in sorted(set(snps["chrom"])):
        m = snps["chrom"] == c
        panel.append(("X" if c == 23 else str(c), snps["pos"][m], snps["ref"][m], snps["alt"][m]))
    path, _ = _grch37_bam(tmp_path, rng, panel, sex=sex)
    res = analyse(path, ref=_Panel(snps), mode="diploid", min_depth=6)
    ids = {c["id"]: c for c in res["cards"]}
    assert {"deeptime_file", "deeptime_authenticity", "deeptime_sex"} <= set(ids)
    f = ids["deeptime_file"]["result"]
    assert f["build"] == "GRCh37"
    assert f["panel"]["sites_covered"] > 0                     # diploid mode used to report zero
    assert f["panel"]["mean_depth_at_covered_sites"] > 1
    assert sum(x["observations"] for x in f["panel"]["chromosomes"]) > 0
    s = ids["deeptime_sex"]["result"]
    assert s["call"] == sex, (s["r_x"], s.get("reason"))
    assert not s["fell_back"]                                  # the panel's own positions settled it
    assert f"R_x = {s['r_x']:.2f}" in ids["deeptime_sex"]["summary"]
    # too few autosomal calls for the ancestry models: it says so rather than running them
    assert "note" in res and "10,000" in res["note"]
    assert all("contamination of a nuclear genome is not measured" in " ".join(ids["deeptime_authenticity"]["limits"])
               for _ in [0])


def test_analyse_refuses_a_genome_on_the_wrong_build(tmp_path):
    """The panel's positions are GRCh37 coordinates, so reads mapped to GRCh38
    would be read at the wrong places. The build comes off the chromosome
    lengths, and the module refuses rather than reporting nonsense."""
    from athanor_core.deeptime.analyse import analyse
    from athanor_core.deeptime.bam import GRCH38
    rng = np.random.default_rng(14)
    snps = _small_panel(rng, per_chrom=20)
    path = str(tmp_path / "grch38.bam")
    sq = [{"SN": c, "LN": L} for c, L in GRCH38.items()]
    with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}) as out:
        pass
    pysam.index(path)
    with Genome(path) as g:
        assert g.build == "GRCh38"
    with pytest.raises(ValueError, match="GRCh37"):
        analyse(path, ref=_Panel(snps))


def _bam_from_reads(path, reads, chrom_len=CHROM_LEN, with_qualities=True):
    sq = [{"SN": "1", "LN": chrom_len}]
    with pysam.AlignmentFile(str(path), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}) as out:
        for n, (start, seq, rev, cigar, name, flag) in enumerate(sorted(reads, key=lambda r: r[0])):
            a = pysam.AlignedSegment()
            a.query_name = name or f"r{n}"
            a.query_sequence = seq
            a.reference_id, a.reference_start = 0, start
            a.flag = (16 if rev else 0) | flag
            a.mapping_quality = 37
            a.cigarstring = cigar
            if with_qualities:
                a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            out.write(a)
    pysam.index(str(path))
    return str(path)


def test_hard_clipped_reads_keep_their_place_in_the_molecule(tmp_path):
    """A hard-clipped base is gone from the record but was in the molecule, so
    the bases that remain are not at the molecule's end."""
    seq = "ACGTACGTAC" * 2
    bam = _bam_from_reads(tmp_path / "clip.bam", [
        (1000, seq, False, "20M", "plain", 0),
        (2000, seq, False, "5H20M", "hardclip", 0),
        (3000, seq, True, "20M5H", "hardclip_rev", 0),
    ])
    with Genome(bam) as g:
        p = pileup(g, "1", np.array([1001, 2001, 3020], dtype=np.int64))
    by_site = {int(s): (int(a), int(b)) for s, a, b in zip(p.site, p.d5, p.d3)}
    assert by_site[0] == (0, 19)              # no clipping: the first base is the molecule's first
    assert by_site[1] == (5, 19)              # five bases were cut from the 5' side
    assert by_site[2][0] == 5                 # reverse read, clip on the molecule's 5' side


def test_reads_without_base_qualities_are_not_taken_for_good_ones(tmp_path):
    seq = "ACGT" * 5
    bam = _bam_from_reads(tmp_path / "noqual.bam", [(1000, seq, False, "20M", "a", 0)], with_qualities=False)
    with Genome(bam) as g:
        p = pileup(g, "1", np.array([1001], dtype=np.int64))
    assert p.reads_without_base_qualities == 1
    assert len(p.filtered(min_qual=20).site) == 0        # not silently given Q30
    assert len(p.filtered(min_qual=0).site) == 1         # available if asked for explicitly


def test_overlapping_mates_count_as_one_molecule(tmp_path):
    """Both reads of a pair cover the same molecule; counting both doubles the
    depth and manufactures confidence."""
    seq = "ACGT" * 5
    reads = [(1000, seq, False, "20M", "pair1", 1), (1000, seq, True, "20M", "pair1", 1),
             (1000, seq, False, "20M", "other", 0)]
    bam = _bam_from_reads(tmp_path / "mates.bam", reads)
    with Genome(bam) as g:
        both = pileup(g, "1", np.array([1001], dtype=np.int64), merge_mates=False)
        merged = pileup(g, "1", np.array([1001], dtype=np.int64))
    assert int(both.counts().sum()) == 3
    assert int(merged.counts().sum()) == 2               # the pair counts once, the other read stands


def test_the_depth_cap_takes_a_fair_sample(tmp_path):
    """Keeping the first reads in coordinate order would keep the ones whose
    bases sit near their 3' ends."""
    rng = np.random.default_rng(31)
    site = 2000
    reads = []
    for k in range(200):
        off = k % 20                                      # where the site sits in the read
        start = site - 1 - off
        reads.append((start, "ACGT" * 5, False, "20M", f"r{k}", 0))
    bam = _bam_from_reads(tmp_path / "deep.bam", reads)
    with Genome(bam) as g:
        capped = pileup(g, "1", np.array([site], dtype=np.int64), max_depth=40, seed=7)
        full = pileup(g, "1", np.array([site], dtype=np.int64), max_depth=0)
    assert len(capped.site) == 40 and len(full.site) == 200
    # the kept reads span the molecule, rather than clustering at one end
    assert capped.d5.min() < 3 and capped.d5.max() > 16
    assert abs(float(capped.d5.mean()) - float(full.d5.mean())) < 3
