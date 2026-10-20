"""Damage-aware haploid consensus with per-base posterior probabilities.

Model (per reference site, true base t in {A, C, G, T}):

  P(t | reads) ∝ P(t) · Π_reads P(observed base | t)

  P(observed | t) = Σ_d P(observed | d) · P(d | t)
    d: the template base after post-mortem deamination, in read orientation.
    Deamination turns C into U (read as T). In double-stranded libraries it
    appears as C→T near read 5' ends and as G→A near 3' ends (the
    complementary strand); Briggs et al., PNAS 104, 14616 (2007).
      P(d=T | t=C) = δ5(i),  i = distance from the 5' end
      P(d=A | t=G) = δ3(j),  j = distance from the 3' end
    Sequencing error: P(observed | d) = 1 - ε if observed = d, else ε / 3,
    with ε = 10^(-Q/10) from the base quality Q.

  Prior: P(t = reference base) = 1 - θ, other bases θ / 3 each. θ (the
  divergence between the sample and the reference) is estimated from the
  data by expectation maximisation.

Bases are coded 0..3 = A, C, G, T in reference orientation; reads mapped to
the reverse strand are complemented into read orientation for the damage
terms, which are defined on the sequenced molecule.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

BASES = "ACGT"
CODE = {b: i for i, b in enumerate(BASES)}
COMP = np.array([3, 2, 1, 0])  # A<->T, C<->G


@dataclass
class Observations:
    """One row per aligned base (from all reads)."""
    site: np.ndarray       # reference position (0-based)
    base: np.ndarray       # observed base code, reference orientation
    qual: np.ndarray       # Phred base quality
    dist5: np.ndarray      # distance from the read's 5' end (0-based)
    dist3: np.ndarray      # distance from the read's 3' end (0-based)
    reverse: np.ndarray    # read mapped to the reverse strand

    def __len__(self) -> int:
        return len(self.site)


@dataclass
class Consensus:
    posterior: np.ndarray          # (L, 4) posterior over A, C, G, T
    depth: np.ndarray              # (L,) reads covering each site
    call: str                      # MAP base, or N where depth < min_depth
    theta: float                   # estimated sample-reference divergence
    notes: list = field(default_factory=list)

    @property
    def confidence(self) -> np.ndarray:
        return self.posterior.max(axis=1)


def log_likelihoods(obs: Observations, delta5: np.ndarray, delta3: np.ndarray) -> np.ndarray:
    """(n_obs, 4) log P(observed base | true base t) for t = A, C, G, T (reference orientation)."""
    eps = np.clip(10.0 ** (-obs.qual / 10.0), 1e-6, 0.75)
    # read-orientation observed base
    o = np.where(obs.reverse, COMP[obs.base], obs.base)
    d5 = delta5[np.minimum(obs.dist5, len(delta5) - 1)]
    d3 = delta3[np.minimum(obs.dist3, len(delta3) - 1)]
    n = len(obs)
    like = np.empty((n, 4))
    for t_ref in range(4):
        t = np.where(obs.reverse, COMP[t_ref], t_ref)        # true base in read orientation
        # P(d | t): mostly d = t; C->T with prob d5, G->A with prob d3
        p_same = np.where(t == 1, 1 - d5, np.where(t == 2, 1 - d3, 1.0))
        p_damaged = 1 - p_same
        damaged_to = np.where(t == 1, 3, np.where(t == 2, 0, t))
        p_obs_given = lambda d: np.where(o == d, 1 - eps, eps / 3)  # noqa: E731
        like[:, t_ref] = p_same * p_obs_given(t) + p_damaged * p_obs_given(damaged_to)
    return np.log(like)


def call_consensus(reference: str, obs: Observations, delta5: np.ndarray, delta3: np.ndarray,
                   min_depth: int = 1, theta0: float = 0.01, em_iterations: int = 20) -> Consensus:
    L = len(reference)
    ref = np.array([CODE.get(b, -1) for b in reference.upper()])
    ll = log_likelihoods(obs, delta5, delta3)
    site_ll = np.zeros((L, 4))
    np.add.at(site_ll, obs.site, ll)
    depth = np.bincount(obs.site, minlength=L)
    covered = depth >= max(min_depth, 1)
    theta = theta0
    post = np.full((L, 4), 0.25)
    for _ in range(em_iterations):
        prior = np.full((L, 4), theta / 3)
        known = ref >= 0
        prior[known, :] = theta / 3
        prior[known, ref[known]] = 1 - theta
        prior[~known, :] = 0.25
        lp = site_ll + np.log(prior)
        lp -= lp.max(axis=1, keepdims=True)
        p = np.exp(lp)
        p /= p.sum(axis=1, keepdims=True)
        post = np.where(covered[:, None], p, prior)
        sel = covered & known
        new_theta = float(np.clip(1 - post[sel, ref[sel]].mean(), 1e-5, 0.5)) if sel.any() else theta
        if abs(new_theta - theta) < 1e-7:
            theta = new_theta
            break
        theta = new_theta
    call = "".join(BASES[i] if covered[k] else "N" for k, i in enumerate(post.argmax(axis=1)))
    return Consensus(posterior=post, depth=depth, call=call, theta=theta)


def calibration(confidence: np.ndarray, correct: np.ndarray, bins: int = 10) -> dict:
    """Reliability of the posteriors: in each confidence bin, the observed
    fraction of correct calls; expected calibration error (ECE) weighted by
    bin size. Bins are equal-width on the error scale log10(1 - p)."""
    err_pred = np.clip(1 - confidence, 1e-12, 1)
    edges = np.logspace(-12, 0, bins + 1)
    idx = np.clip(np.digitize(err_pred, edges) - 1, 0, bins - 1)
    rows, ece, n = [], 0.0, len(confidence)
    for b in range(bins):
        m = idx == b
        if not m.any():
            continue
        pred, obs = float(err_pred[m].mean()), float(1 - correct[m].mean())
        rows.append({"bin_upper_error": float(edges[b + 1]), "n": int(m.sum()),
                     "predicted_error": pred, "observed_error": obs})
        ece += m.sum() / n * abs(pred - obs)
    return {"bins": rows, "expected_calibration_error": float(ece), "n_sites": int(n),
            "observed_error_rate": float(1 - correct.mean()) if n else None,
            "predicted_error_rate": float(err_pred.mean()) if n else None}
