"""Deamination profiles from mapDamage2 output.

mapDamage2 (Jónsson et al., Bioinformatics 29, 1682, 2013) writes the
frequency of C->T at each position from the 5' end (5pCtoT_freq.txt) and of
G->A from the 3' end (3pGtoA_freq.txt), relative to the reference. These
frequencies include genuine sample-reference differences, which do not
depend on position; that background is estimated from the innermost
positions reported and subtracted, leaving the damage-induced excess used as
δ5(i) and δ3(i) by the consensus model.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def _read_freq(path: Path) -> np.ndarray:
    vals = []
    for line in path.read_text().splitlines()[1:]:
        p = line.split()
        if len(p) >= 2:
            vals.append(float(p[1]))
    return np.array(vals)


def damage_profiles(mapdamage_dir: str | Path, background_positions: int = 5) -> dict:
    d = Path(mapdamage_dir)
    ct = _read_freq(d / "5pCtoT_freq.txt")
    ga = _read_freq(d / "3pGtoA_freq.txt")
    bg_ct = float(ct[-background_positions:].mean())
    bg_ga = float(ga[-background_positions:].mean())
    return {
        "raw_5p_CtoT": ct.tolist(), "raw_3p_GtoA": ga.tolist(),
        "background_CtoT": bg_ct, "background_GtoA": bg_ga,
        "delta5": np.clip(ct - bg_ct, 0, 1), "delta3": np.clip(ga - bg_ga, 0, 1),
    }


def bayesian_estimates(mapdamage_dir: str | Path) -> dict:
    """Posterior summaries of mapDamage2's Briggs-model parameters
    (Stats_out_MCMC_iter_summ_stat.csv), keyed by parameter name."""
    f = Path(mapdamage_dir) / "Stats_out_MCMC_iter_summ_stat.csv"
    if not f.exists():
        return {}
    import csv
    rows = list(csv.reader(f.read_text().splitlines()))
    header, out = rows[0], {}
    # file layout: first column is the statistic (Mean, Std., quantiles...), one column per parameter
    stats = {r[0].strip('"'): r[1:] for r in rows[1:] if r}
    params = [h.strip('"') for h in header[1:]]
    for i, name in enumerate(params):
        out[name] = {k: float(v[i]) for k, v in stats.items() if i < len(v) and _isfloat(v[i])}
    return out


def _isfloat(x: str) -> bool:
    try:
        float(x)
        return True
    except ValueError:
        return False


def damage_from_observations(obs, sample_codes: np.ndarray, use_site: np.ndarray | None = None,
                             max_pos: int = 25, background_positions: int = 5) -> dict:
    """Deamination profile measured against the sample's own consensus
    instead of the reference, so genuine sample-reference differences do not
    enter it. Same definitions as mapDamage2: the fraction of C (in read
    orientation) read as T at distance i from the 5' end, and of G read as A
    at distance i from the 3' end. 3' ends are counted only on whole
    molecules. ``sample_codes``: per-site base code 0..3 (or -1 = unknown);
    ``use_site``: per-site mask of trusted consensus bases."""
    from .consensus import COMP
    true_ref = sample_codes[obs.site]
    ok = true_ref >= 0
    if use_site is not None:
        ok &= use_site[obs.site]
    t = np.where(obs.reverse, COMP[np.clip(true_ref, 0, 3)], true_ref)    # read orientation
    o = np.where(obs.reverse, COMP[obs.base], obs.base)
    full = obs.full if obs.full is not None else np.ones(len(obs), bool)

    def profile(true_b, obs_b, dist, mask):
        num, den = np.zeros(max_pos), np.zeros(max_pos)
        m = ok & mask & (t == true_b) & (dist < max_pos)
        np.add.at(den, dist[m], 1)
        np.add.at(num, dist[m], (o[m] == obs_b))
        return np.divide(num, den, out=np.zeros(max_pos), where=den > 0), den

    ct, n5 = profile(1, 3, obs.dist5, np.ones(len(obs), bool))
    ga, n3 = profile(2, 0, obs.dist3, full)
    bg_ct = float(ct[-background_positions:].mean())
    bg_ga = float(ga[-background_positions:].mean())
    return {"raw_5p_CtoT": ct.tolist(), "raw_3p_GtoA": ga.tolist(), "n_5p": n5.astype(int).tolist(),
            "n_3p": n3.astype(int).tolist(), "background_CtoT": bg_ct, "background_GtoA": bg_ga,
            "delta5": np.clip(ct - bg_ct, 0, 1), "delta3": np.clip(ga - bg_ga, 0, 1),
            "measured_against": "sample consensus"}


def wilson_interval(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (float(max(0.0, c - h)), float(min(1.0, c + h)))


def terminal_signal(profile: dict, interior_from: int = 10) -> dict:
    """Evidence of post-mortem deamination: C->T at the first 5' base against
    the interior of the molecules, with 95% Wilson intervals. Genuinely
    ancient double-stranded libraries show an excess at the ends; partially
    UDG-treated libraries keep it only at the terminal one or two bases."""
    ct, n5 = np.array(profile["raw_5p_CtoT"]), np.array(profile["n_5p"])
    ga, n3 = np.array(profile["raw_3p_GtoA"]), np.array(profile["n_3p"])

    def row(f, n, i):
        k = int(round(f[i] * n[i]))
        return {"frequency": float(f[i]), "n": int(n[i]), "ci95": wilson_interval(k, int(n[i]))}

    def interior(f, n):
        k = int(round(float((f[interior_from:] * n[interior_from:]).sum())))
        nn = int(n[interior_from:].sum())
        return {"frequency": k / nn if nn else None, "n": nn, "ci95": wilson_interval(k, nn)}

    out = {"CtoT_5p_position1": row(ct, n5, 0), "CtoT_5p_position2": row(ct, n5, 1),
           "CtoT_interior": interior(ct, n5),
           "GtoA_3p_position1": row(ga, n3, 0), "GtoA_interior": interior(ga, n3)}
    lo1 = out["CtoT_5p_position1"]["ci95"][0]
    hi_int = out["CtoT_interior"]["ci95"][1]
    out["terminal_excess"] = bool(lo1 > hi_int)
    return out
