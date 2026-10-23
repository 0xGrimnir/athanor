"""Modern ancestry: regional composition, chromosome painting, runs of
homozygosity and X-chromosome ancestry.

Composition is supervised admixture with fixed reference frequencies, the
model of ADMIXTURE's projection mode (Alexander, Novembre & Lange, Genome
Res. 19, 1655, 2009) and of the "frappe" EM (Tang et al., Genet. Epidemiol.
28, 289, 2005): each allele in your genome is drawn from one of K reference
populations, chosen independently with probabilities q, so at SNP m

    g_m ~ Binomial(n_m, f_m),   f_m = sum_k q_k P_km,

with n_m = 2 (1 on a male X). q is the maximum-likelihood estimate, found by
EM accelerated with SQUAREM (Varadhan & Roland, Scand. J. Stat. 35, 335,
2008); intervals come from a bootstrap over ~5 cM blocks, so linked SNPs
are resampled together. Panel frequencies are shrunk towards the pooled
frequency (panel_freqs), so that panels of different sizes compete fairly.

Chromosome painting is a hidden Markov model along each chromosome whose
hidden state is the pair of continental ancestries of your two copies at
that point (unordered: chip data are unphased, so which copy is which is
unknown). Each copy's ancestry switches at recombination events at rate
T generations per Morgan, to a continent drawn from your composition, as in
the admixture models behind LAMP and RFMix (Sankararaman et al., AJHG 82,
290, 2008; Maples et al., AJHG 93, 278, 2013). Emissions are the genotype
probabilities under the two continents' allele frequencies; because nearby
SNPs are correlated, log-emissions are down-weighted by a factor fitted on
simulated genomes with known segments (the reference pack records it).

Runs of homozygosity follow PLINK 1.9's --homozyg scanning-window algorithm
(Purcell et al., AJHG 81, 559, 2007; Chang et al., GigaScience 4, 7, 2015)
with its parameters made explicit; CI checks the segments against PLINK.
"""
from __future__ import annotations

from contextlib import nullcontext

import numpy as np

EPS = 1e-6


def _one_thread():
    """These estimators do many small matrix-vector products, which a
    multi-threaded BLAS runs many times slower than a single thread."""
    try:
        from threadpoolctl import threadpool_limits
        return threadpool_limits(1)
    except Exception:  # noqa: BLE001
        return nullcontext()


# ------------------------------------------------------------------ frequency encoding
def encode_freq(p: np.ndarray) -> np.ndarray:
    """uint8 arcsine-square-root encoding: fine steps near 0 and 1 (where a
    rare allele's frequency matters most), at most ~0.003 error at 0.5."""
    t = np.arcsin(np.sqrt(np.clip(p, 0.0, 1.0))) * (2 / np.pi)
    return np.round(t * 254).astype(np.uint8)


def decode_freq(v: np.ndarray) -> np.ndarray:
    p = np.sin(v.astype(np.float64) * (np.pi / 508.0)) ** 2
    return np.clip(p, 1e-5, 1 - 1e-5)


# ------------------------------------------------------------------ supervised admixture
def _loglik(g, n, w, f):
    return float(np.sum(w * (g * np.log(f) + (n - g) * np.log1p(-f))))


def _em_step(q, P, g, n, w, wn):
    f = np.clip(q @ P, EPS, 1 - EPS)
    a = w * g / f
    b = w * (n - g) / (1 - f)
    up = q * (P @ a + (b.sum() - P @ b))
    return up / wn


def admix(g: np.ndarray, P: np.ndarray, n: np.ndarray | None = None, w: np.ndarray | None = None,
          q0: np.ndarray | None = None, tol: float = 1e-6, max_iter: int = 5000) -> dict:
    """Maximum-likelihood ancestry proportions q (K,) for alt-allele counts g
    (M,) given reference frequencies P (K, M). Missing SNPs must already be
    removed. n: alleles per SNP (2, or 1 on a male X); w: SNP weights (the
    bootstrap's block multiplicities)."""
    K = P.shape[0]
    g = g.astype(np.float64)
    n = np.full(g.shape, 2.0) if n is None else n.astype(np.float64)
    w = np.ones(g.shape) if w is None else w.astype(np.float64)
    wn = float(np.sum(w * n))
    q = np.full(K, 1.0 / K) if q0 is None else np.clip(q0, 1e-6, None) / np.clip(q0, 1e-6, None).sum()
    ll = _loglik(g, n, w, np.clip(q @ P, EPS, 1 - EPS))
    it = 0
    while it < max_iter:
        q1 = _em_step(q, P, g, n, w, wn)
        q2 = _em_step(q1, P, g, n, w, wn)
        r = q1 - q
        v = q2 - q1 - r
        nv = np.sqrt(v @ v)
        it += 2
        if nv < 1e-14:
            qn = q2
        else:
            # SQUAREM step length; shortened (towards the plain double EM step,
            # alpha = -1) until every proportion stays positive, so no
            # component is pushed to zero where the multiplicative EM update
            # could not bring it back
            alpha = min(-np.sqrt(r @ r) / nv, -1.0)
            while True:
                qn = q - 2 * alpha * r + alpha * alpha * v
                if np.all(qn > 0) or alpha >= -1.0:
                    break
                alpha = (alpha - 1.0) / 2
            qn = np.clip(qn, 1e-300, None)
            qn /= qn.sum()
            qn = _em_step(qn, P, g, n, w, wn)              # stabilising step
            it += 1
        lln = _loglik(g, n, w, np.clip(qn @ P, EPS, 1 - EPS))
        if lln < ll - 1e-9:                               # extrapolation overshot: keep the plain EM step
            qn = q2
            lln = _loglik(g, n, w, np.clip(qn @ P, EPS, 1 - EPS))
        done = float(np.max(np.abs(qn - q))) < tol
        q, ll = qn, lln
        if done:
            break
    return {"q": q, "loglik": ll, "iterations": it}


def bootstrap(g, P, block, q_hat, n=None, reps: int = 100, seed: int = 1) -> np.ndarray:
    """(reps, K) proportions from resampling ~5 cM blocks with replacement."""
    rng = np.random.default_rng(seed)
    ub, inv = np.unique(block, return_inverse=True)
    out = np.empty((reps, P.shape[0]))
    for r in range(reps):
        cnt = np.bincount(rng.integers(0, len(ub), len(ub)), minlength=len(ub)).astype(np.float64)
        w = cnt[inv]
        keep = w > 0
        out[r] = admix(g[keep], P[:, keep], None if n is None else n[keep], w[keep], q0=q_hat, tol=1e-4)["q"]
    return out


def panel_freqs(alt: np.ndarray, tot: np.ndarray, kappa: float = 2.0) -> np.ndarray:
    """Reference-panel allele frequencies, shrunk towards the pooled frequency
    of all panels: p_k = (a_k + kappa * p_all) / (t_k + kappa).

    A fixed pseudo-count such as (a + 0.5) / (t + 1) gives each panel a floor
    frequency that depends on its size, and at the many SNPs where nobody in
    a region carries an allele that floor alone decides between panels: a
    person is then pulled towards the largest panels, whatever their
    ancestry. Shrinking towards the pooled frequency (a Beta prior centred on
    it, worth kappa alleles) gives every panel the same tiny value where all
    lack the allele, so panel size no longer matters there.
    """
    alt = alt.astype(np.float64)
    tot = tot.astype(np.float64)
    pall = (alt.sum(0) + 0.5) / (tot.sum(0) + 1.0)
    return (alt + kappa * pall) / (tot + kappa)


def composition(dos: np.ndarray, P: np.ndarray, block: np.ndarray, n_alleles: np.ndarray | None = None,
                reps: int = 100, seed: int = 1) -> dict:
    """dos: alt-allele dosage per reference SNP (-1 missing)."""
    with _one_thread():
        return _composition(dos, P, block, n_alleles, reps, seed)


def _composition(dos, P, block, n_alleles, reps, seed):
    m = dos >= 0
    g = dos[m].astype(np.float64)
    n = None if n_alleles is None else n_alleles[m]
    Pm = np.ascontiguousarray(P[:, m])
    fit = admix(g, Pm, n)
    # The bootstrap refits only the components the point estimate uses
    # (above 0.1%): the others are reported as not detected, and dropping them
    # makes each refit several times faster.
    act = np.flatnonzero(fit["q"] > 1e-3)
    bs = np.zeros((reps, P.shape[0]))
    if reps:
        bs[:, act] = bootstrap(g, Pm[act], block[m], fit["q"][act], n, reps=reps, seed=seed)
    return {"q": fit["q"], "loglik": fit["loglik"], "iterations": fit["iterations"], "snps": int(m.sum()),
            "boot": bs}


def interval(boot: np.ndarray, cols: list[int] | np.ndarray | None = None) -> tuple[float, float]:
    x = boot if cols is None else boot[:, cols].sum(1)
    if x.ndim > 1:
        x = x.sum(1)
    return float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))


# ------------------------------------------------------------------ chromosome painting
def _emission(g: np.ndarray, Pc: np.ndarray, err: float = 0.01) -> np.ndarray:
    """(M, K, K) probabilities of genotype g (0/1/2, -1 missing) when one copy
    comes from continent a and the other from b."""
    pa = Pc.T[:, :, None]          # (M, K, 1)
    pb = Pc.T[:, None, :]          # (M, 1, K)
    p0 = (1 - pa) * (1 - pb)
    p2 = pa * pb
    p1 = 1 - p0 - p2
    e = np.where(g[:, None, None] == 0, p0, np.where(g[:, None, None] == 2, p2, p1))
    e = (1 - err) * e + err / 3
    e = np.where(g[:, None, None] < 0, 1.0, e)
    return e


def paint(dos: np.ndarray, Pc: np.ndarray, chrom: np.ndarray, cm: np.ndarray, q: np.ndarray,
          generations: float = 10.0, temper: float = 0.5, err: float = 0.01) -> dict:
    """Posterior probabilities of the (unordered) pair of continental
    ancestries at each SNP, per chromosome.

    dos: genotypes at the SNPs used (0/1/2, -1 missing); Pc: (K, M) continental
    frequencies; q: (K,) your proportions, the prior for each copy's ancestry.
    Returns per-SNP MAP pair (a, b) with a <= b, and its posterior.
    """
    with _one_thread():
        return _paint(dos, Pc, chrom, cm, q, generations, temper, err)


def _paint(dos, Pc, chrom, cm, q, generations, temper, err):
    K = Pc.shape[0]
    q = np.clip(q, 1e-4, None)
    q = q / q.sum()
    best_a = np.zeros(len(dos), np.int8)
    best_b = np.zeros(len(dos), np.int8)
    post_best = np.zeros(len(dos), np.float32)
    marg = np.zeros((len(dos), K), np.float32)       # expected copies of each continent (0..2)
    prior = np.outer(q, q)
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        if len(idx) < 2:
            continue
        E = _emission(dos[idx], Pc[:, idx], err) ** temper            # (L, K, K)
        d = np.diff(cm[idx]) / 100.0                                    # Morgans
        stay = np.exp(-generations * np.clip(d, 0, None))
        L = len(idx)
        fw = np.empty((L, K, K))
        a = prior * E[0]
        s = a.sum()
        fw[0] = a / s
        for i in range(1, L):
            e_ = stay[i - 1]
            j = 1 - e_
            A = fw[i - 1]
            r = A.sum(1)
            cc = A.sum(0)
            nxt = (e_ * e_) * A + (e_ * j) * (np.outer(r, q) + np.outer(q, cc)) + (j * j) * np.outer(q, q)
            a = nxt * E[i]
            s = a.sum()
            fw[i] = a / s
        bw = np.ones((K, K))
        post = np.empty((L, K, K))
        post[L - 1] = fw[L - 1]
        for i in range(L - 2, -1, -1):
            e_ = stay[i]
            j = 1 - e_
            B = bw * E[i + 1]
            # transpose of the transition applied to B
            r = B @ q                      # sum over b' of B[a', b'] q[b']
            cc = q @ B                     # sum over a' of q[a'] B[a', b']
            tot = q @ B @ q
            nb = (e_ * e_) * B + (e_ * j) * (r[:, None] + cc[None, :]) + (j * j) * tot
            nb /= nb.sum()
            bw = nb
            p = fw[i] * bw
            post[i] = p / p.sum()
        sym = post + np.transpose(post, (0, 2, 1))
        iu = np.triu_indices(K)
        flat = sym[:, iu[0], iu[1]].copy()
        flat[:, iu[0] == iu[1]] /= 2                      # diagonal counted once
        k = flat.argmax(1)
        best_a[idx] = iu[0][k]
        best_b[idx] = iu[1][k]
        post_best[idx] = flat[np.arange(L), k]
        marg[idx] = post.sum(2) + post.sum(1)
    return {"a": best_a, "b": best_b, "posterior": post_best, "copies": marg}


def segments(chrom, pos, a, b, post, min_snps: int = 20) -> list[dict]:
    """Collapse per-SNP calls into runs of the same pair, dropping runs shorter
    than ``min_snps`` (absorbed into their neighbours)."""
    out = []
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        key = a[idx].astype(int) * 64 + b[idx]
        # smooth: runs shorter than min_snps take the previous run's state
        starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
        ends = np.r_[starts[1:], len(idx)]
        runs = [[int(key[s]), int(s), int(e)] for s, e in zip(starts, ends)]
        merged = []
        for r in runs:
            if merged and (r[2] - r[1] < min_snps or merged[-1][0] == r[0]):
                merged[-1][2] = r[2]
            else:
                merged.append(r)
        for k, s, e in merged:
            sl = idx[s:e]
            out.append({"chrom": int(c), "start": int(pos[sl[0]]), "end": int(pos[sl[-1]]),
                        "a": k // 64, "b": k % 64, "snps": int(e - s),
                        "posterior": float(np.mean(post[sl]))})
    return out


# ------------------------------------------------------------------ runs of homozygosity
ROH_DEFAULTS = {"window_snp": 50, "window_het": 1, "window_missing": 5, "window_threshold": 0.05,
                "min_snp": 100, "min_kb": 1000, "max_density_kb": 50, "max_gap_kb": 1000, "max_het": 999}


def roh(chrom: np.ndarray, pos: np.ndarray, het: np.ndarray, missing: np.ndarray, **kw) -> list[dict]:
    """PLINK 1.9 --homozyg. chrom/pos sorted; het/missing boolean per SNP.

    Each run of SNPs with no gap above max_gap_kb is scanned with windows of
    window_snp consecutive SNPs; a window is homozygous if it has at most
    window_het heterozygous and window_missing missing calls. A SNP is a
    candidate when at least window_threshold of the windows covering it are
    homozygous; consecutive candidates form a segment, kept if it has at least
    min_snp SNPs, spans at least min_kb, averages no more than max_density_kb
    per SNP and has at most max_het heterozygous calls.
    """
    p = {**ROH_DEFAULTS, **kw}
    W = p["window_snp"]
    out = []
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        cp = pos[idx]
        # split at large gaps
        brk = np.flatnonzero(np.diff(cp) > p["max_gap_kb"] * 1000) + 1
        for part in np.split(np.arange(len(idx)), brk):
            L = len(part)
            if L < W:
                continue
            h = het[idx[part]].astype(np.int32)
            ms = missing[idx[part]].astype(np.int32)
            ch = np.r_[0, np.cumsum(h)]
            cmi = np.r_[0, np.cumsum(ms)]
            nw = L - W + 1
            wh = ch[W:W + nw] - ch[:nw]
            wm = cmi[W:W + nw] - cmi[:nw]
            hom_w = ((wh <= p["window_het"]) & (wm <= p["window_missing"])).astype(np.int32)
            # windows covering SNP i: starts max(0, i-W+1) .. min(i, nw-1)
            cw = np.r_[0, np.cumsum(hom_w)]
            i = np.arange(L)
            lo = np.maximum(0, i - W + 1)
            hi = np.minimum(i, nw - 1)
            cover = hi - lo + 1
            frac = (cw[hi + 1] - cw[lo]) / cover
            cand = frac >= p["window_threshold"] - 1e-12
            # PLINK only counts called SNPs at the segment's ends
            st = np.flatnonzero(np.diff(np.r_[0, cand.astype(np.int8), 0]) == 1)
            en = np.flatnonzero(np.diff(np.r_[0, cand.astype(np.int8), 0]) == -1)
            for s, e in zip(st, en):
                seg = part[s:e]
                called = ~missing[idx[seg]]
                if not called.any():
                    continue
                cs = seg[called]
                n = len(cs)
                start, end = int(cp[cs[0]]), int(cp[cs[-1]])
                kb = (end - start + 1) / 1000
                nhet = int(het[idx[seg]].sum())
                if n >= p["min_snp"] and kb >= p["min_kb"] and kb / n <= p["max_density_kb"] and nhet <= p["max_het"]:
                    out.append({"chrom": int(c), "start": start, "end": end, "kb": kb, "snps": n, "het": nhet})
    return out


def roh_summary(segs: list[dict], autosome_mb: float = 2881.0) -> dict:
    kb = np.array([s["kb"] for s in segs]) if segs else np.zeros(0)
    def tot(lo, hi=np.inf):
        return float(kb[(kb >= lo * 1000) & (kb < hi * 1000)].sum() / 1000)
    return {"segments": len(segs), "total_mb": float(kb.sum() / 1000),
            "froh": float(kb.sum() / 1000 / autosome_mb),
            "by_length_mb": {"1-2": tot(1, 2), "2-5": tot(2, 5), "5-10": tot(5, 10), "10+": tot(10)},
            "longest_mb": float(kb.max() / 1000) if len(kb) else 0.0,
            "long_segments": int((kb >= 10000).sum())}
