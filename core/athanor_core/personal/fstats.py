"""f-statistics with a weighted block jackknife: f4, f4-ratios, D-statistics, qpAdm.

Definitions follow Patterson et al., "Ancient admixture in human history",
Genetics 192, 1065 (2012):
  f4(A, B; C, D) = mean over SNPs of (a - b)(c - d), with a, b, c, d the allele
  frequencies; D(A, B; C, D) = sum (a - b)(c - d) / sum (a + b - 2ab)(c + d - 2cd).
Standard errors use the weighted block jackknife (Busing, Meijer & van der
Leeden, Stat. Comput. 9, 3, 1999) over blocks of about 5 cM, as in ADMIXTOOLS.

qpAdm (Haak et al., Nature 522, 207, 2015; Harney et al., Genetics 217,
iyaa045, 2021): a target T is modelled as a mixture of sources S_1..S_n.
If T = sum_i w_i S_i with sum_i w_i = 1, then for every pair of "right"
(outgroup) populations R_0, R_j, sum_i w_i f4(T, S_i; R_0, R_j) = 0. The
weights are fitted by generalised least squares on these f4 statistics,
using their jackknife covariance; the fit's chi-squared with
(m - n) degrees of freedom (m right populations) gives the p-value for the
model. Weights' standard errors come from the jackknife.

Frequencies are arrays (populations x SNPs) with NaN where missing; an
individual enters as its genotype dosage / 2.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def blocks_from_cm(chrom: np.ndarray, cm: np.ndarray, size_cm: float = 5.0) -> np.ndarray:
    """Block index per SNP: consecutive windows of ``size_cm`` on each chromosome."""
    key = chrom.astype(np.int64) * 10_000 + np.floor(cm / size_cm).astype(np.int64)
    _, inv = np.unique(key, return_inverse=True)
    return inv.astype(np.int32)


class Jackknife:
    """Per-block sums of numerators (and denominators) for a set of statistics."""

    def __init__(self, block: np.ndarray, nblocks: int | None = None):
        self.block = block
        self.nb = int(nblocks if nblocks is not None else block.max() + 1)

    def sums(self, values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        """values: (k, S) with NaN allowed where mask is False. Returns (k, nb)."""
        v = np.where(mask, values, 0.0)
        out = np.empty((v.shape[0], self.nb))
        for i in range(v.shape[0]):
            out[i] = np.bincount(self.block, weights=v[i], minlength=self.nb)
        return out

    def counts(self, mask: np.ndarray) -> np.ndarray:
        return np.bincount(self.block, weights=mask.astype(float), minlength=self.nb)


def jackknife_estimate(num: np.ndarray, den: np.ndarray, func=None):
    """Weighted block jackknife of ratio statistics.

    num: (k, nb) block sums of numerators; den: (nb,) or (k, nb) block sums of
    denominators (SNP counts for means). ``func`` maps the vector of k ratios
    to the statistic(s) of interest (default: the ratios themselves).
    Returns (estimate, covariance, n_blocks_used)."""
    den = np.broadcast_to(den, num.shape)
    keep = den[0] > 0
    num, den = num[:, keep], den[:, keep]
    g = num.shape[1]
    if g < 2:
        raise ValueError("too few blocks with data for a jackknife")
    f = func or (lambda x: x)
    tot_n, tot_d = num.sum(1), den.sum(1)
    theta = np.atleast_1d(f(tot_n / tot_d))
    w = den[0] / den[0].sum()                         # block weight by SNP count
    loo = np.array([np.atleast_1d(f((tot_n - num[:, j]) / (tot_d - den[:, j]))) for j in range(g)])
    h = 1.0 / w
    tj = g * theta - ((1 - w)[:, None] * loo).sum(0)
    pseudo = h[:, None] * theta[None, :] - (h - 1)[:, None] * loo
    dev = pseudo - tj[None, :]
    cov = (dev[:, :, None] * dev[:, None, :] / (h - 1)[:, None, None]).sum(0) / g
    return theta, cov, g


def _f4_values(F: np.ndarray, quads: list[tuple[int, int, int, int]]) -> tuple[np.ndarray, np.ndarray]:
    vals = np.empty((len(quads), F.shape[1]))
    ok = np.ones(F.shape[1], bool)
    for k, (a, b, c, d) in enumerate(quads):
        vals[k] = (F[a] - F[b]) * (F[c] - F[d])
    used = sorted({i for q in quads for i in q})
    for i in used:
        ok &= ~np.isnan(F[i])
    return vals, ok


def f4(F: np.ndarray, block: np.ndarray, a: int, b: int, c: int, d: int) -> dict:
    jk = Jackknife(block)
    vals, ok = _f4_values(F, [(a, b, c, d)])
    est, cov, g = jackknife_estimate(jk.sums(vals, ok), jk.counts(ok))
    se = float(np.sqrt(cov[0, 0]))
    return {"f4": float(est[0]), "se": se, "z": float(est[0] / se) if se > 0 else 0.0, "snps": int(ok.sum()), "blocks": g}


def f4_ratio(F: np.ndarray, block: np.ndarray, num: tuple, den: tuple) -> dict:
    """alpha = f4(num) / f4(den), on the SNPs where all populations are present."""
    jk = Jackknife(block)
    vals, ok = _f4_values(F, [num, den])
    sums = jk.sums(vals, ok)
    cnt = jk.counts(ok)
    est, cov, g = jackknife_estimate(sums, cnt, func=lambda x: np.array([x[0] / x[1]]))
    se = float(np.sqrt(cov[0, 0]))
    return {"alpha": float(est[0]), "se": se, "snps": int(ok.sum()), "blocks": g}


def d_stat(F: np.ndarray, block: np.ndarray, a: int, b: int, c: int, d: int) -> dict:
    jk = Jackknife(block)
    A, B, C, D = F[a], F[b], F[c], F[d]
    ok = ~(np.isnan(A) | np.isnan(B) | np.isnan(C) | np.isnan(D))
    num = (A - B) * (C - D)
    den = (A + B - 2 * A * B) * (C + D - 2 * C * D)
    sums = jk.sums(np.vstack([num, den]), ok)
    cnt = jk.counts(ok)
    est, cov, g = jackknife_estimate(sums, cnt, func=lambda x: np.array([x[0] / x[1]]))
    se = float(np.sqrt(cov[0, 0]))
    return {"D": float(est[0]), "se": se, "z": float(est[0] / se) if se > 0 else 0.0, "snps": int(ok.sum()), "blocks": g}


def _gls_weights(X: np.ndarray, Sinv: np.ndarray) -> np.ndarray:
    """min over w (sum w = 1) of r' Sinv r with r = X' w; X is (n, m-1)."""
    A = X @ Sinv @ X.T
    one = np.ones(X.shape[0])
    try:
        z = np.linalg.solve(A, one)
    except np.linalg.LinAlgError:
        z = np.linalg.lstsq(A, one, rcond=None)[0]
    return z / z.sum()


def qpadm(F: np.ndarray, block: np.ndarray, target: int, sources: list[int], right: list[int],
          iterations: int = 3) -> dict:
    """Fit target = sum w_i sources_i against right populations right[0..m-1]."""
    n, m = len(sources), len(right)
    if m <= n:
        raise ValueError("qpAdm needs more right populations than sources")
    quads = [(target, s, right[0], r) for s in sources for r in right[1:]]
    vals, ok = _f4_values(F, quads)
    jk = Jackknife(block)
    sums = jk.sums(vals, ok)
    cnt = jk.counts(ok)
    keep = cnt > 0
    sums, cnt = sums[:, keep], cnt[keep]
    g = len(cnt)
    tot = sums.sum(1) / cnt.sum()
    X = tot.reshape(n, m - 1)
    loo = (sums.sum(1)[:, None] - sums) / (cnt.sum() - cnt)[None, :]     # (n(m-1), g)
    Xloo = loo.T.reshape(g, n, m - 1)
    wts = cnt / cnt.sum()
    h = 1.0 / wts

    def jk_cov(vals_full: np.ndarray, vals_loo: np.ndarray) -> np.ndarray:
        tj = g * vals_full - ((1 - wts)[:, None] * vals_loo).sum(0)
        pseudo = h[:, None] * vals_full[None, :] - (h - 1)[:, None] * vals_loo
        dev = pseudo - tj[None, :]
        return (dev[:, :, None] * dev[:, None, :] / (h - 1)[:, None, None]).sum(0) / g

    # covariance of the full f4 matrix entries
    cov_X = jk_cov(X.ravel(), loo.T)
    Sinv = np.eye(m - 1)
    w = _gls_weights(X, Sinv)
    for _ in range(iterations):
        # covariance of r = X' w for the current weights (w held fixed)
        Wm = np.kron(w[:, None], np.eye(m - 1)).T            # (m-1, n(m-1)): r = Wm vec(X)
        cov_r = Wm @ cov_X @ Wm.T
        Sinv = np.linalg.pinv(cov_r)
        w = _gls_weights(X, Sinv)
    r = X.T @ w
    Wm = np.kron(w[:, None], np.eye(m - 1)).T
    cov_r = Wm @ cov_X @ Wm.T
    chi2 = float(r @ np.linalg.pinv(cov_r) @ r)
    df = m - n
    p = float(stats.chi2.sf(chi2, df))
    w_loo = np.array([_gls_weights(Xloo[j], Sinv) for j in range(g)])
    cov_w = jk_cov(w, w_loo)
    se = np.sqrt(np.clip(np.diag(cov_w), 0, None))
    return {"weights": w.tolist(), "se": se.tolist(), "chi2": chi2, "df": df, "p": p,
            "snps": int(ok.sum()), "blocks": g}
