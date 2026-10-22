"""Personal Genome: file parsing, build conversion, f-statistics, haplogroups, PCA."""
import numpy as np
import pytest

from athanor_core.personal import fstats
from athanor_core.personal.haplo import Tree, build_mt_tree
from athanor_core.personal.pca import normalise, project
from athanor_core.personal.rawfile import MISSING, file_check, parse
from athanor_core.personal.writers import write


def _toy(tmp_path, fmt, n=400, seed=0):
    rng = np.random.default_rng(seed)
    rsid = np.arange(1, n + 1)
    chrom = np.repeat(np.arange(1, 5), n // 4)
    pos = np.sort(rng.choice(np.arange(10000, 900000), n, replace=False))
    ref = rng.integers(0, 4, n)
    alt = (ref + rng.integers(1, 4, n)) % 4
    dos = rng.integers(-1, 3, n)
    p = tmp_path / f"f_{fmt}.txt"
    write(p, fmt, rsid, chrom, pos, ref, alt, dos)
    return p, rsid, chrom, pos, ref, alt, dos


@pytest.mark.parametrize("fmt,vendor", [("23andme", "23andMe"), ("ancestry", "AncestryDNA"),
                                        ("myheritage", "MyHeritage"), ("ftdna", "FamilyTreeDNA"),
                                        ("livingdna", "Living DNA"), ("vcf", "VCF")])
def test_every_vendor_format_round_trips(tmp_path, fmt, vendor):
    p, rsid, chrom, pos, ref, alt, dos = _toy(tmp_path, fmt)
    g = parse(p)
    assert g.vendor == vendor
    assert len(g) == len(pos)
    assert list(g.rsid) == list(rsid)
    assert list(g.pos) == list(pos)
    called = dos >= 0
    # genotypes come back as the same unordered pair of bases
    got = np.sort(np.vstack([g.a1, g.a2]), axis=0)
    want_lo = np.where(dos == 2, alt, ref)
    want_hi = np.where(dos == 0, ref, alt)
    want = np.sort(np.vstack([want_lo, want_hi]), axis=0)
    assert np.array_equal(got[:, called], want[:, called])
    assert np.all(g.a1[~called] == MISSING)


def test_gzip_and_zip_are_opened(tmp_path):
    import gzip
    import zipfile
    p, *_ = _toy(tmp_path, "23andme")
    raw = p.read_bytes()
    gz = tmp_path / "f.txt.gz"
    gz.write_bytes(gzip.compress(raw))
    zp = tmp_path / "f.zip"
    with zipfile.ZipFile(zp, "w") as z:
        z.writestr("AncestryDNA.txt", raw)
    assert len(parse(gz)) == len(parse(p))
    assert len(parse(zp)) == len(parse(p))


def test_file_check_reads_sex_from_x_and_y(tmp_path):
    rng = np.random.default_rng(2)
    n = 600
    rsid = np.arange(1, n + 1)
    chrom = np.array([1] * 400 + [23] * 150 + [24] * 50)
    pos = np.arange(1, n + 1) * 100
    ref = rng.integers(0, 4, n)
    alt = (ref + 1) % 4
    male = np.where(chrom == 23, rng.choice([0, 2], n), rng.integers(0, 3, n))
    p = tmp_path / "m.txt"
    write(p, "23andme", rsid, chrom, pos, ref, alt, male)
    fc = file_check(parse(p))
    assert fc["genetic_sex"] == "XY"
    female = np.where(chrom == 24, -1, rng.integers(0, 3, n))
    female = np.where(chrom == 23, rng.integers(0, 3, n), female)
    p2 = tmp_path / "f.txt"
    write(p2, "23andme", rsid, chrom, pos, ref, alt, female)
    assert file_check(parse(p2))["genetic_sex"] == "XX"


def test_f4_is_zero_without_admixture_and_detects_it_with():
    rng = np.random.default_rng(5)
    S = 40000
    anc = rng.beta(1, 1, S)
    def drift(p, t=0.05):
        return np.clip(p + rng.normal(0, t * np.sqrt(p * (1 - p) + 1e-9), S), 1e-4, 1 - 1e-4)
    a, b = drift(anc), drift(anc)
    c, d = drift(anc), drift(anc)
    block = np.repeat(np.arange(40), S // 40)
    F = np.vstack([a, b, c, d])
    r = fstats.f4(F, block, 0, 1, 2, 3)
    assert abs(r["z"]) < 4, r
    # make A share drift with C: f4(A, B; C, D) must become positive
    shared = rng.normal(0, 0.02, S)
    F2 = np.vstack([np.clip(a + shared, 1e-4, 1 - 1e-4), b, np.clip(c + shared, 1e-4, 1 - 1e-4), d])
    r2 = fstats.f4(F2, block, 0, 1, 2, 3)
    assert r2["f4"] > 0 and r2["z"] > 6, r2


def test_qpadm_recovers_known_mixture():
    """Populations are sampled as genotypes, not exact frequencies: a real data
    set always carries sampling noise, and the jackknife needs it to estimate
    the covariance the fit is weighted by."""
    rng = np.random.default_rng(11)
    S, N = 60000, 40
    block = np.repeat(np.arange(40), S // 40)
    root = rng.beta(1, 1, S)

    def drift(p, t):
        return np.clip(p + rng.normal(0, t * np.sqrt(p * (1 - p) + 1e-9), S), 1e-3, 1 - 1e-3)

    def sample(p):
        return rng.binomial(2, p, (N, S)).mean(axis=0) / 2

    base = [drift(root, 0.25) for _ in range(3)]
    sources = [drift(b, 0.15) for b in base]
    rights = [drift(root, 0.3)] + [drift(b, 0.3) for b in base] + [drift(root, 0.35)]
    w = np.array([0.5, 0.3, 0.2])
    target = np.clip(sum(wi * s for wi, s in zip(w, sources)), 1e-3, 1 - 1e-3)
    F = np.vstack([sample(x) for x in sources + rights + [target]])
    r = fstats.qpadm(F, block, len(F) - 1, [0, 1, 2], [3, 4, 5, 6, 7])
    assert np.allclose(r["weights"], w, atol=0.05), (r["weights"], r["se"])
    assert r["p"] > 0.01, r
    assert sum(r["weights"]) == pytest.approx(1.0)


def test_f4_ratio_matches_a_known_proportion():
    rng = np.random.default_rng(13)
    S = 80000
    block = np.repeat(np.arange(40), S // 40)

    def drift(p, t):
        return np.clip(p + rng.normal(0, t * np.sqrt(p * (1 - p) + 1e-9), S), 1e-3, 1 - 1e-3)

    root = rng.beta(1, 1, S)
    archaic_anc = drift(root, 0.4)
    a_hi = drift(archaic_anc, 0.1)        # the "Altai" genome
    a_lo = drift(archaic_anc, 0.1)        # the "Vindija" genome
    human = drift(root, 0.3)
    africa = drift(human, 0.2)
    out = drift(root, 0.5)                 # the "Chimp" outgroup
    alpha = 0.3
    mixed = np.clip(alpha * a_lo + (1 - alpha) * africa, 1e-3, 1 - 1e-3)
    F = np.vstack([a_hi, out, a_lo, africa, mixed])
    r = fstats.f4_ratio(F, block, (0, 1, 4, 3), (0, 1, 2, 3))
    assert abs(r["alpha"] - alpha) < 0.05, r


def test_haplogroup_tree_calls_the_defining_branch():
    xml = """<?xml version="1.0"?><phylotree>
    <haplogroup name="R"><details></details>
      <haplogroup name="R1"><details><poly>100A</poly></details>
        <haplogroup name="R1a"><details><poly>200C</poly></details></haplogroup>
        <haplogroup name="R1b"><details><poly>300T</poly></details></haplogroup>
      </haplogroup>
      <haplogroup name="R2"><details><poly>400T</poly></details></haplogroup>
    </haplogroup></phylotree>"""
    root_seq = "G" * 500          # every mutation above is read against this root
    tree = Tree(build_mt_tree(xml, root_seq, "100A\t3.0\n"))
    assert set(tree.positions) == {100, 200, 300, 400}
    call = tree.call({100: "A", 200: "C", 300: "G", 400: "G"})
    assert call["call"] == "R1a"
    assert [p["name"] for p in call["path"]] == ["R", "R1", "R1a"]
    # one stray derived call off the lineage must not outweigh a consistent path
    stray = tree.call({100: "A", 200: "C", 400: "T"})
    assert stray["call"] == "R1a"
    # ancestral everywhere: stays at the root
    assert tree.call({100: "G", 200: "G", 300: "G", 400: "G"})["call"] == "R"
    assert tree.is_ancestor(tree.index["R1"], tree.index["R1a"])
    assert not tree.is_ancestor(tree.index["R1a"], tree.index["R1"])


def test_pca_projection_reproduces_reference_coordinates():
    rng = np.random.default_rng(17)
    S, N, K = 3000, 40, 4
    p = rng.uniform(0.1, 0.9, S)
    g = rng.binomial(2, p[None, :], (N, S)).astype(float)
    scale = np.sqrt(p * (1 - p))
    Z = (g - 2 * p) / scale
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    coords = (U * s)[:, :K]
    loadings = Vt[:K].T
    z = normalise(g[0], p, scale)
    x, n = project(z, loadings, min_snps=100)
    assert n == S
    assert np.allclose(x, coords[0], atol=1e-6)
    # with a quarter of the SNPs missing the placement is close, not exact
    z2 = z.copy()
    z2[rng.random(S) < 0.25] = np.nan
    x2, n2 = project(z2, loadings, min_snps=100)
    assert n2 < S
    assert np.linalg.norm(x2 - coords[0]) < 0.25 * np.linalg.norm(coords[0])


def test_liftover_maps_through_a_chain(tmp_path):
    from athanor_core.personal.liftover import Liftover
    chain = tmp_path / "c.chain"
    # one chain, two blocks: target 1000..1099 -> query 2000..2099, then a gap
    chain.write_text("chain 100 chr1 1000000 + 1000 1200 chr1 1000000 + 2000 2200 1\n100\t50\t50\n100\n\n")
    lo = Liftover(chain)
    got = lo.map("chr1", np.array([1001, 1100, 1101, 1151, 1251]))
    assert got[0] == 2001            # inside the first block
    assert got[1] == 2100            # last base of the first block
    assert got[2] == -1              # inside the gap
    assert got[3] == 2151            # inside the second block
    assert got[4] == -1              # past the last block
