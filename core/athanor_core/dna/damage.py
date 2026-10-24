"""Ancient DNA damage in place: single-stranded overhangs and deaminated cytosines on a molecule.

Post-mortem, DNA breaks into short fragments whose ends are often single-stranded, and cytosines lose
their amino group (hydrolytic deamination), becoming uracils that a polymerase reads as thymine. Both
happen far more in single-stranded stretches, so an ancient molecule carries most of its uracils near
its 5' ends (Briggs et al., PNAS 104, 14616, 2007). In sequencing reads this appears as C->T at the
5' end and, from the complementary strand, G->A at the 3' end.

Two ways to put damage on a helix:

* the Briggs model: 5' overhang lengths are geometric with parameter lambda (the chance that an
  overhang reaches past position i is (1-lambda)^(i+1)); a cytosine deaminates with probability
  delta_s in an overhang and delta_d in the double-stranded part. CI checks these formulas, and the
  read-level profile they imply (with the nick parameter nu), against gargammel's deamSim simulator.
* a measured profile: the C->T frequency at each distance from the 5' end that the Ancient DNA Lab
  measured on real reads, used directly as each cytosine's chance of being a uracil.

A deaminated cytosine is drawn as uracil in place: the exocyclic N4 amine becomes the O4 carbonyl
(C4=O 1.23 A) on the same ring.
"""
from __future__ import annotations

import numpy as np

from .builder import build


def p_single(i, lam: float):
    """Chance that the base i positions from a 5' end (0 = the terminal base) lies in that end's
    single-stranded overhang, for geometric overhang lengths with parameter lambda."""
    return (1.0 - lam) ** (np.asarray(i, float) + 1.0)


def p_uracil(i, lam: float, dd: float, ds: float):
    """Chance that a cytosine i positions from its strand's 5' end is a uracil (Briggs model)."""
    p = p_single(i, lam)
    return ds * p + dd * (1.0 - p)


def read_profile(n: int, nu: float, lam: float, dd: float, ds: float):
    """C->T at each distance from a read's 5' end, as gargammel 1.1.4's deamSim produces it for long
    fragments (its "legacy" Briggs model): a 5' overhang on every molecule, deamination delta_s in it;
    in the double-stranded part, delta_d until the first nick (probability nu per position), after
    which the read copies the other strand, where that damage shows as G->A instead."""
    i = np.arange(n, dtype=float)
    p = p_single(i, lam)
    return ds * p + dd * (1.0 - p) * (1.0 - nu) ** (i + 1.0)


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
        if name == "N4":
            out.append((c, n, "DU", "O4", o4))
        else:
            out.append((c, n, "DU", name, xyz))
    return out


def molecule(seq: str, overhang5: int, overhang3: int, uracils: set) -> dict:
    """A damaged molecule: strand I is ``seq`` (5'->3'); its 5' end overhangs by ``overhang5`` bases
    (strand II missing there) and strand II's 5' end overhangs by ``overhang3`` bases at the other end
    (strand I missing there). ``uracils``: set of (strand, position from that strand's 5' end) for
    cytosines that are uracil. Single strands are drawn continuing the B-form helix of the duplex they
    belonged to; in reality they are flexible."""
    n = len(seq)
    if overhang5 + overhang3 > n - 2:
        raise ValueError("the overhangs must leave at least two base pairs")
    m = build(seq, "B")
    res = _residues(m["atoms"])
    out = []
    for (chain, num), atoms in res.items():
        k = num - 1 if chain == "A" else 2 * n - num          # base-pair index along strand I
        if chain == "A" and k >= n - overhang3:
            continue                                          # opposite strand II's 5' overhang
        if chain == "B" and k < overhang5:
            continue                                          # opposite strand I's 5' overhang
        dist = k if chain == "A" else n - 1 - k               # distance from this strand's own 5' end
        if atoms[0][2][-1] == "C" and (chain, dist) in uracils:
            atoms = to_uracil(atoms)
        out.extend(atoms)
    return {"atoms": out, "sequence": seq, "overhang5": overhang5, "overhang3": overhang3,
            "frames": m["frames"]}


def sample(seq: str, lam: float, dd: float, ds: float, rng=None, profile=None) -> dict:
    """Draw overhang lengths and uracils for one molecule. With ``profile`` (measured C->T by distance
    from the 5' end), each cytosine's chance comes from it and the overhangs are drawn from lambda."""
    rng = rng or np.random.default_rng()
    n = len(seq)
    for _ in range(1000):
        o5 = int(rng.geometric(lam)) - 1          # numpy's geometric counts trials: shift to 0-based
        o3 = int(rng.geometric(lam)) - 1
        if o5 + o3 <= n - 2:
            break
    else:
        o5 = o3 = 0
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    s2 = "".join(comp[b] for b in reversed(seq))           # strand II, 5'->3'
    ura = set()
    probs = {"A": [], "B": []}
    for chain, s, own in (("A", seq, o5), ("B", s2, o3)):
        present = n - (o3 if chain == "A" else o5)            # strand I lacks its 3' tail opposite strand II's overhang
        for i, b in enumerate(s[:present]):
            if b != "C":
                probs[chain].append(None)
                continue
            if profile is not None and i < len(profile):
                p = float(profile[i])
            elif profile is not None:
                p = float(profile[-1])
            else:
                p = ds if i < own else dd
            probs[chain].append(p)
            if rng.random() < p:
                ura.add((chain, i))
    m = molecule(seq, o5, o3, ura)
    m.update({"uracils": sorted(ura), "p_uracil": probs, "model": "profile" if profile is not None else "briggs"})
    return m
