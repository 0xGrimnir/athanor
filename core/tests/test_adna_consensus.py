"""Statistical checks of the damage-aware consensus on synthetic data with a
known truth (independent of the external aDNA tools)."""
import numpy as np

from athanor_core.adna.consensus import BASES, COMP, Observations, calibration, call_consensus


def simulate(seed=7, L=6000, theta=0.03, coverage=6.0, delta_max=0.35, decay=0.55):
    rng = np.random.default_rng(seed)
    ref = rng.integers(0, 4, L)
    truth = ref.copy()
    sub = rng.random(L) < theta
    truth[sub] = (truth[sub] + rng.integers(1, 4, sub.sum())) % 4
    delta = delta_max * decay ** np.arange(30)
    rows = {k: [] for k in ("site", "base", "qual", "dist5", "dist3", "reverse")}
    n_frag = int(coverage * L / 55)
    for _ in range(n_frag):
        flen = int(np.clip(rng.lognormal(np.log(50), 0.35), 25, 120))
        start = rng.integers(0, L - flen)
        rev = rng.random() < 0.5
        mol = truth[start:start + flen].copy()
        if rev:
            mol = COMP[mol][::-1]                       # molecule as sequenced (read orientation)
        i = np.arange(flen)
        d5, d3 = delta[np.minimum(i, 29)], delta[np.minimum(flen - 1 - i, 29)]
        ct = (mol == 1) & (rng.random(flen) < d5)
        ga = (mol == 2) & (rng.random(flen) < d3)
        mol[ct], mol[ga] = 3, 0
        q = rng.integers(15, 41, flen)
        err = rng.random(flen) < 10.0 ** (-q / 10)
        mol[err] = (mol[err] + rng.integers(1, 4, err.sum())) % 4
        obs_ref = COMP[mol][::-1] if rev else mol       # back to reference orientation
        sites = np.arange(start, start + flen)
        d5r = i[::-1] if rev else i                     # distance from 5' end, per reference position
        rows["site"].append(sites); rows["base"].append(obs_ref)
        rows["qual"].append(q[::-1] if rev else q)
        rows["dist5"].append(d5r); rows["dist3"].append(flen - 1 - d5r)
        rows["reverse"].append(np.full(flen, rev))
    obs = Observations(**{k: np.concatenate(v) for k, v in rows.items()})
    return "".join(BASES[b] for b in ref), truth, obs, delta


def test_high_coverage_consensus_is_accurate_and_recovers_divergence():
    ref, truth, obs, delta = simulate()
    c = call_consensus(ref, obs, delta, delta, min_depth=1)
    cov = c.depth > 0
    correct = np.array([BASES[t] for t in truth])[cov] == np.array(list(c.call))[cov]
    assert abs(c.theta - 0.03) < 0.006, c.theta
    assert correct.mean() > 0.995


def test_low_coverage_posteriors_are_calibrated():
    """At ~1.5x coverage many sites rest on one damaged read, so errors are
    frequent enough to test that predicted and observed error rates agree."""
    ref, truth, obs, delta = simulate(seed=3, L=60000, coverage=1.5)
    c = call_consensus(ref, obs, delta, delta, min_depth=1)
    cov = c.depth > 0
    correct = np.array([BASES[t] for t in truth])[cov] == np.array(list(c.call))[cov]
    cal = calibration(c.confidence[cov], correct)
    assert (~correct).sum() > 100, "test needs enough errors to measure calibration"
    assert cal["expected_calibration_error"] < 0.01, cal
    # predicted and observed overall error rates agree within a factor of 1.5
    ratio = cal["predicted_error_rate"] / max(cal["observed_error_rate"], 1e-9)
    assert 0.67 < ratio < 1.5, cal


def test_ignoring_damage_is_measurably_worse():
    ref, truth, obs, delta = simulate(seed=11, theta=0.02, coverage=3.0)
    aware = call_consensus(ref, obs, delta, delta)
    naive = call_consensus(ref, obs, np.zeros_like(delta), np.zeros_like(delta))
    cov = aware.depth > 0
    t = np.array([BASES[x] for x in truth])[cov]
    err_aware = (np.array(list(aware.call))[cov] != t).mean()
    err_naive = (np.array(list(naive.call))[cov] != t).mean()
    assert err_aware < err_naive
    # the naive model is overconfident: its predicted error is far below its observed error
    cal_naive = calibration(naive.confidence[cov], np.array(list(naive.call))[cov] == t)
    cal_aware = calibration(aware.confidence[cov], np.array(list(aware.call))[cov] == t)
    assert cal_aware["expected_calibration_error"] < cal_naive["expected_calibration_error"]
