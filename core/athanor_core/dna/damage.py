"""Ancient DNA damage in place: single-stranded overhangs and deaminated cytosines on a molecule.

Post-mortem, DNA breaks into short fragments whose ends are often single-stranded, and cytosines lose
their amino group (hydrolytic deamination), becoming uracils that a polymerase reads as thymine. Both
happen far more in single-stranded stretches, so an ancient molecule carries most of its uracils near
its ends (Briggs et al., PNAS 104, 14616, 2007). In double-stranded library sequencing this appears as
C->T at the reads' 5' ends and, from the complementary strand, G->A at their 3' ends.

The Briggs model, as published:

* each end of a molecule has a single-stranded overhang; it leaves the 5' or the 3' strand protruding
  with equal probability, and its length is geometric with parameter lambda (it reaches past position
  i with probability (1-lambda)^(i+1));
* a cytosine deaminates with probability delta_s in single-stranded DNA and delta_d in double-stranded;
* library preparation fills in 3' overhangs and keeps 5' ones, and a read copies the other strand after
  its first nick (probability nu per double-stranded position; molecules with a second nick are lost,
  so the first nick is uniform over the double-stranded stretch).

gargammel's deamSim implements this from version 1.1.5; up to 1.1.4 it drew a 5' overhang at both
ends every time and placed the nick geometrically, which doubles the terminal damage for the same
parameters (gargammel 1.1.5 release notes). Both are available here: ``model="briggs"`` (the paper)
and ``model="legacy"`` (gargammel 1.1.4, which the Ancient DNA Lab's simulations used). CI checks each
against the matching deamSim mode, position by position.

A deaminated cytosine is drawn as uracil in place: the exocyclic N4 amine becomes the O4 carbonyl
(C4=O 1.23 A) on the same ring.
"""
from __future__ import annotations

import numpy as np

from .builder import build

# Briggs et al. (2007), maximum-likelihood estimates for Neanderthal DNA: nu, lambda, delta_d, delta_s
BRIGGS_2007 = (0.024, 0.36, 0.0097, 0.68)


def p_beyond(i, lam: float):
    """Chance that an overhang (given there is one) reaches past position i (0 = the terminal base)."""
    return (1.0 - lam) ** (np.asarray(i, float) + 1.0)


def p_uracil(i_from5, i_from3, lam: float, dd: float, ds: float, model: str = "briggs"):
    """Chance that a cytosine is a uracil, given its distances from its strand's 5' and 3' ends: it is
    single-stranded if its 5' end has a 5' overhang reaching it, or its 3' end a 3' overhang."""
    if model == "legacy":           # gargammel <= 1.1.4: a 5' overhang at every end, no 3' overhangs
        p = p_beyond(i_from5, lam)
    else:
        p5, p3 = 0.5 * p_beyond(i_from5, lam), 0.5 * p_beyond(i_from3, lam)
        p = p5 + p3 - p5 * p3
    return ds * p + dd * (1.0 - p)


def _overhang_dist(lam: float, L: int, model: str):
    """Probabilities of overhang lengths 0..L at one end."""
    o = np.arange(L + 1)
    g = lam * (1 - lam) ** o
    if model == "legacy":
        return g
    g = 0.5 * g
    g[0] += 0.5
    return g


def read_profile(L: int, nu: float, lam: float, dd: float, ds: float, model: str = "briggs"):
    """Expected C->T at each distance from a read's 5' end and G->A from its 3' end, for reads of
    length L from a double-stranded library, as deamSim produces them (``-damage`` for "briggs",
    gargammel 1.1.4's ``-damage`` = 1.1.5's ``-damagelegacy`` for "legacy"). Exact enumeration over the
    overhang lengths at both ends and the position of the first nick."""
    P = _overhang_dist(lam, L, model)
    ct, ga = np.zeros(L), np.zeros(L)
    i = np.arange(L)
    tot = 0.0
    for o5 in range(L - 1):
        for o3 in range(L - 1 - o5):              # the overhangs leave at least two double-stranded bases
            w = P[o5] * P[o3]
            if w < 1e-14:
                continue
            tot += w
            last = L - o3 - 1                     # last double-stranded base
            in5, in3 = i < o5, (L - i) <= o3
            ds_ = ~in5 & ~in3
            if model == "legacy":
                # a nick at each position with probability nu, looked for from the 5' end
                before = (1 - nu) ** (i + 1)      # no nick at or before i
            else:
                n = last - o5
                p_nick = nu / (n * nu + (1 - nu)) if n >= 1 else 0.0
                # the first nick lies after one of the bases o5..last-1, uniformly; "no nick" holds the rest
                before = np.clip(1 - (i - o5) * p_nick, 0, 1)
            c = np.where(in5, ds, np.where(ds_, dd * before, 0.0))
            g = np.where(in3, ds, np.where(ds_, dd * (1 - before), 0.0))
            ct += w * c
            ga += w * g
    return ct / tot, (ga / tot)[::-1]


def _residues(atoms):
    out = {}
    for a in atoms:
        out.setdefault((a[0], a[1]), []).append(a)
    return out


def to_uracil(res_atoms):
    """A cytosine nucleotide's atoms with the base deaminated to uracil (N4 -> O4 at 1.23 A)."""
    by = {a[3]: a for a in res_atoms}
    c4, n4 = np.asarray(by["C4"][4]), np.asarray(by["N4"][4])
    v = n4 - c4
    o4 = c4 + 1.23 * v / np.linalg.norm(v)
    out = []
    for c, n, r, name, xyz in res_atoms:
        out.append((c, n, "DU", "O4", o4) if name == "N4" else (c, n, "DU", name, xyz))
    return out


def molecule(seq: str, left: tuple, right: tuple, uracils: set) -> dict:
    """A damaged molecule. Strand I is ``seq`` (5'->3', left to right). ``left`` and ``right`` are
    (kind, length) for each end's overhang: kind "5" leaves the 5' strand protruding (strand I at the
    left end, strand II at the right), "3" the 3' strand. ``uracils``: set of (strand, distance from that
    strand's 5' end in the full-length helix) for cytosines that are uracil. Single strands are drawn
    continuing the B-form helix of the duplex they belonged to; in reality they are flexible."""
    n = len(seq)
    if left[1] + right[1] > n - 2:
        raise ValueError("the overhangs must leave at least two base pairs")
    m = build(seq, "B")
    gone = set()
    for k in range(left[1]):
        gone.add(("B", 2 * n - k) if left[0] == "5" else ("A", k + 1))
    for k in range(n - right[1], n):
        gone.add(("A", k + 1) if right[0] == "5" else ("B", 2 * n - k))
    present = [key for key in dict.fromkeys((a[0], a[1]) for a in m["atoms"]) if key not in gone]
    first = {}
    for c, num in present:
        first.setdefault(c, num)               # each strand's 5'-most remaining residue
    out = []
    for (chain, num), atoms in _residues(m["atoms"]).items():
        if (chain, num) in gone:
            continue
        k = num - 1 if chain == "A" else 2 * n - num
        dist = k if chain == "A" else n - 1 - k
        if first[chain] == num:                # a new 5' end is a hydroxyl: no phosphate
            atoms = [a for a in atoms if a[3] not in ("P", "OP1", "OP2")]
        if atoms[0][2][-1] == "C" and (chain, dist) in uracils:
            atoms = to_uracil(atoms)
        out.extend(atoms)
    return {"atoms": out, "sequence": seq, "left": list(left), "right": list(right), "frames": m["frames"]}


def single_stranded(n: int, left: tuple, right: tuple):
    """{(strand, distance from its 5' end): True} for the bases left single-stranded by the overhangs."""
    ss = set()
    for k in range(left[1]):                    # left end: strand I's 5' end / strand II's 3' end
        ss.add(("A", k) if left[0] == "5" else ("B", n - 1 - k))
    for k in range(n - right[1], n):            # right end: strand I's 3' end / strand II's 5' end
        ss.add(("B", n - 1 - k) if right[0] == "5" else ("A", k))
    return ss


def sample(seq: str, lam: float, dd: float, ds: float, rng=None, profile=None, model: str = "briggs") -> dict:
    """Draw one molecule: each end's overhang (5' or 3' with equal probability, geometric length), then
    each cytosine is a uracil with probability delta_s if single-stranded and delta_d if not. With a
    measured ``profile`` (C->T by distance from the 5' end) the cytosines' chances come from it instead."""
    rng = rng or np.random.default_rng()
    lam = float(np.clip(lam, 0.01, 0.99))
    dd, ds = float(np.clip(dd, 0, 1)), float(np.clip(ds, 0, 1))
    n = len(seq)
    for _ in range(1000):
        ends = []
        for _ in range(2):
            kind = "5" if (model == "legacy" or rng.random() < 0.5) else "3"
            ends.append((kind, int(rng.geometric(lam)) - 1))       # numpy counts trials: shift to 0-based
        if ends[0][1] + ends[1][1] <= n - 2:
            break
    else:
        ends = [("5", 0), ("5", 0)]
    left, right = ends
    ss = single_stranded(n, left, right)
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    strands = {"A": seq, "B": "".join(comp[b] for b in reversed(seq))}
    gone = {("B", n - 1 - k) for k in range(left[1])} if left[0] == "5" else {("A", k) for k in range(left[1])}
    gone |= {("A", k) for k in range(n - right[1], n)} if right[0] == "5" else {("B", n - 1 - k) for k in range(n - right[1], n)}
    ura, probs = set(), {"A": [], "B": []}
    for chain, s in strands.items():
        for i, b in enumerate(s):
            if b != "C" or (chain, i) in gone:
                probs[chain].append(None)
                continue
            if profile is not None:
                p = float(np.clip(profile[min(i, len(profile) - 1)], 0, 1))
            else:
                p = ds if (chain, i) in ss else dd
            probs[chain].append(p)
            if rng.random() < p:
                ura.add((chain, i))
    m = molecule(seq, left, right, ura)
    m.update({"uracils": sorted(ura), "p_uracil": probs, "single_stranded": sorted(ss),
              "model": "profile" if profile is not None else model})
    return m
