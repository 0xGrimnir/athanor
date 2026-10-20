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
