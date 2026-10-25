"""Deep Time: an ancient genome simulated from a known reference, with known
genotypes and known deamination, read back through the pipeline."""
import numpy as np
import pytest

pysam = pytest.importorskip("pysam")

from athanor_core.deeptime import calls as C                      # noqa: E402
from athanor_core.deeptime.bam import GRCH37, Genome, length_summary, plain, wilson   # noqa: E402
from athanor_core.deeptime.sites import damage_profile, pileup    # noqa: E402

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


def _simulate(tmp_path, rng, ref_seq, pos, ref, alt, depth=16, read_len=60, damage=DAMAGE0, with_y=0):
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
            for j in range(len(seq)):
                d5 = (len(seq) - 1 - j) if rev else j
                d3 = j if rev else (len(seq) - 1 - j)
                frm, to = ("G", "A") if rev else ("C", "T")
                frm3, to3 = ("C", "T") if rev else ("G", "A")
                if seq[j] == frm and rng.random() < damage * (1 - LAM) ** d5:
                    seq[j] = to
                elif seq[j] == frm3 and rng.random() < damage * (1 - LAM) ** d3:
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
    return {"path": path, "pos": pos, "ref": ref, "alt": alt, "truth": truth, "sequence": ref_seq}


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
    with Genome(ancient["path"]) as g:
        p = pileup(g, "1", ancient["pos"])
        d = damage_profile(p, ancient["ref"], ancient["alt"])
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
        counts = C.region_counts(g, "1", start - 1, end)
    r = C.heterozygosity(counts, min_depth=8, min_quality=30)
    assert r["callable_sites"] > 100_000
    # the heterozygous sites that are deep enough to be called
    deep = counts.sum(1)[pos - start] >= 8
    het = deep & (truth == 1)
    transition = np.array([{(r_, a_)} <= {(0, 2), (2, 0), (1, 3), (3, 1)}
                           for r_, a_ in zip(ancient["ref"], ancient["alt"])])
    expected_tv = int((het & ~transition).sum())
    assert abs(r["transversions"] - expected_tv) < 0.1 * expected_tv
    assert r["transversion_heterozygosity"] == pytest.approx(expected_tv / r["callable_sites"], rel=0.1)
    # damage turns homozygous C/G sites into apparent C/T and G/A heterozygotes
    assert r["transitions"] > 1.3 * int((het & transition).sum())
    assert r["heterozygous_sites"] > int(het.sum())


def test_sex_from_the_sex_chromosome_read_counts(ancient):
    with Genome(ancient["path"]) as g:
        s = g.sex()
        assert s["reads_Y"] == 3000
        assert s["call"] == "XY"                   # no X reads at all: R_y = 1
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
