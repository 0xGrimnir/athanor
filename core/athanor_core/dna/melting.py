"""DNA duplex melting from nearest-neighbour thermodynamics.

Parameters: the unified set of SantaLucia, Proc. Natl. Acad. Sci. USA 95, 1460 (1998), Table 2
(1 M NaCl): enthalpy and entropy of each of the ten Watson-Crick nearest-neighbour stacks, two
initiation terms (by the identity of each terminal pair), and the symmetry correction for
self-complementary duplexes. Salt enters through the entropy, 0.368 cal/(K mol) x ln[Na+] per
phosphate pair (SantaLucia 1998, eq. 8; N - 1 for an N-base-pair duplex).

Two models are given:

* ``two_state`` - the duplex is either fully formed or fully dissociated; this is the model the
  parameters were fitted with and the one the melting temperature is reported from;
* ``zipper`` - every contiguous helical stretch i..j is allowed (frayed ends), so the probability that
  each base pair is closed can be shown along the helix as it melts. It reduces to the two-state model
  when only the full duplex is counted.
"""
from __future__ import annotations

import math

import numpy as np

R = 1.98720425864083    # cal/(K mol)
T0 = 273.15
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}

# (dH kcal/mol, dS cal/K/mol), SantaLucia 1998 Table 2; keys are 5'-XY-3' on the top strand
NN = {"AA": (-7.9, -22.2), "TT": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
      "CA": (-8.5, -22.7), "TG": (-8.5, -22.7), "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
      "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
      "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9)}
INIT_GC = (0.1, -2.8)      # initiation with a terminal G.C pair (per end)
INIT_AT = (2.3, 4.1)       # initiation with a terminal A.T pair (per end)
SYM = (0.0, -1.4)          # self-complementary duplexes
SALT_DS = 0.368            # cal/(K mol) per phosphate pair per ln[Na+]
SOURCE = ("SantaLucia, A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor "
          "thermodynamics, Proc. Natl. Acad. Sci. USA 95, 1460 (1998)")


def clean(seq: str) -> str:
    s = "".join(seq.split()).upper()
    bad = sorted(set(s) - set("ACGT"))
    if bad:
        raise ValueError(f"only A, C, G and T are allowed (found {''.join(bad)})")
    if len(s) < 2:
        raise ValueError("a duplex needs at least two base pairs")
    return s


def revcomp(s: str) -> str:
    return "".join(COMP[c] for c in reversed(s))


def self_complementary(s: str) -> bool:
    return s == revcomp(s)


def _end(c):
    return INIT_GC if c in "GC" else INIT_AT


def duplex_terms(s: str, na: float = 1.0) -> tuple[float, float]:
    """(dH kcal/mol, dS cal/K/mol) for the perfect duplex of s with its complement at [Na+] = na M."""
    dh = ds = 0.0
    for i in range(len(s) - 1):
        h, e = NN[s[i:i + 2]]
        dh += h
        ds += e
    for c in (s[0], s[-1]):
        h, e = _end(c)
        dh += h
        ds += e
    if self_complementary(s):
        dh += SYM[0]
        ds += SYM[1]
    ds += SALT_DS * (len(s) - 1) * math.log(na)
    return dh, ds


def tm(seq: str, ct: float = 1e-6, na: float = 1.0) -> dict:
    """Two-state melting temperature (deg C) for total strand concentration ct (M; equal amounts of the
    two strands, or ct of a self-complementary strand) at [Na+] = na (M)."""
    s = clean(seq)
    dh, ds = duplex_terms(s, na)
    sc = self_complementary(s)
    x = 1.0 if sc else 4.0
    t = dh * 1000.0 / (ds + R * math.log(ct / x)) - T0
    dg37 = dh - (37 + T0) * ds / 1000.0
    return {"sequence": s, "tm": t, "dH": dh, "dS": ds, "dG37": dg37, "self_complementary": sc,
            "ct": ct, "na": na, "source": SOURCE}


def _fraction(K: np.ndarray, ct: float, sc: bool) -> np.ndarray:
    """Fraction of strands in duplex, theta, for association constant K (1/M).
    Non-self-complementary, [A]0 = [B]0 = ct/2: K = theta / ((1-theta)^2 ct/2).
    Self-complementary, [A]0 = ct: K = theta / (2 (1-theta)^2 ct)."""
    a = K * (ct / 2 if not sc else 2 * ct)
    # a (1-t)^2 = t  ->  a t^2 - (2a+1) t + a = 0, the root in [0, 1]
    a = np.asarray(a, float)
    with np.errstate(over="ignore", invalid="ignore"):
        t = np.where(a > 1e-12, ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * a), a)
    return np.clip(t, 0.0, 1.0)


def curve(seq: str, temps_c, ct: float = 1e-6, na: float = 1.0) -> np.ndarray:
    """Two-state fraction of strands in duplex at each temperature."""
    s = clean(seq)
    dh, ds = duplex_terms(s, na)
    T = np.asarray(temps_c, float) + T0
    K = np.exp(-(dh * 1000.0 - T * ds) / (R * T))
    return _fraction(K, ct, self_complementary(s))


def zipper(seq: str, temp_c: float, ct: float = 1e-6, na: float = 1.0) -> dict:
    """Probability that each base pair is closed at temp_c, allowing any contiguous helical stretch
    i..j (ends frayed). Each stretch has the stacks it contains, initiation terms for its own end pairs,
    and the salt correction for its own phosphates."""
    s = clean(seq)
    n = len(s)
    T = temp_c + T0
    sc = self_complementary(s)
    g_stack = np.array([(NN[s[i:i + 2]][0] * 1000 - T * (NN[s[i:i + 2]][1] + SALT_DS * math.log(na)))
                        for i in range(n - 1)])
    g_end = np.array([(_end(c)[0] * 1000 - T * _end(c)[1]) for c in s])
    g_sym = (SYM[0] * 1000 - T * SYM[1]) if sc else 0.0
    cum = np.concatenate([[0.0], np.cumsum(g_stack)])
    i, j = np.triu_indices(n, k=1)                     # stretches of at least two pairs
    G = cum[j] - cum[i] + g_end[i] + g_end[j] + g_sym
    gmin = G.min()
    w = np.exp(-(G - gmin) / (R * T))
    Z = w.sum()
    K = Z * math.exp(-gmin / (R * T))
    theta = float(_fraction(np.array(K), ct, sc))
    closed = np.zeros(n)
    # sum of w over stretches covering each position k: difference array
    acc = np.zeros(n + 1)
    np.add.at(acc, i, w)
    np.add.at(acc, j + 1, -w)
    closed = np.cumsum(acc)[:n] / Z * theta
    return {"sequence": s, "temp": temp_c, "fraction_duplex": theta, "closed": closed.tolist()}
