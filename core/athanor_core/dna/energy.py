"""Why the helix holds together: base-pairing and base-stacking interaction energies.

A relaxed helix (OpenMM, with hydrogens) is cut into its bases. Each base is capped where its
glycosidic bond was: the sugar's C1' becomes a hydrogen 1.01 A from N9 (purines) or N1 (pyrimidines)
along the old bond, as in the S22 benchmark's nucleobase complexes (Jurecka et al., Phys. Chem. Chem.
Phys. 8, 1985, 2006). Then

* pairing energy of base pair i: E(base_I + base_II) - E(base_I) - E(base_II)
* stacking energy of step i: E(pair i + pair i+1) - E(pair i) - E(pair i+1)

computed in the gas phase with GFN2-xTB (fast; every pair and step) or DFT (B3LYP-D3(BJ)/def2-TZVP with
counterpoise correction; slow, one pair or step at a time). Both methods are checked in CI on the S22
set of noncovalent complexes against coupled-cluster CCSD(T)/CBS reference energies, including its
adenine-thymine Watson-Crick and stacked pairs. These are interaction energies of the bases alone in
vacuum: water, ions, the backbone and entropy all change what holds a real duplex together, so they
explain the pattern (G-C pairs bind more strongly than A-T; stacking depends on the step), not the
melting free energy.
"""
from __future__ import annotations

import numpy as np

from ..molecule import Molecule

HARTREE_KCAL = 627.509474
PURINE = {"A", "G"}
BASE_HEAVY = {
    "A": ["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"],
    "G": ["N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"],
    "C": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
    "T": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6"],
}
BASE_H = {"A": ["H8", "H61", "H62", "H2"], "G": ["H8", "H1", "H21", "H22"], "C": ["H41", "H42", "H5", "H6"],
          "T": ["H3", "H71", "H72", "H73", "H6"]}


def base_fragment(heavy: dict, hyd: dict, base: str):
    """Symbols and coordinates (A) of one capped nucleobase from a nucleotide's heavy atoms and hydrogens."""
    sym, xyz = [], []
    for n in BASE_HEAVY[base]:
        sym.append(n.lstrip("0123456789")[0])
        xyz.append(heavy[n])
    for n in BASE_H[base]:
        if n in hyd:
            sym.append("H")
            xyz.append(hyd[n])
    n_gly = "N9" if base in PURINE else "N1"
    v = heavy["C1'"] - heavy[n_gly]
    sym.append("H")
    xyz.append(heavy[n_gly] + 1.01 * v / np.linalg.norm(v))
    return sym, np.array(xyz)


def fragments(atoms, hydrogens, n: int):
    """[(strand I base, strand II base)] fragments for each base pair of a duplex of n bp."""
    by, hy = {}, {}
    names = {}
    for c, num, res, name, x in atoms:
        by.setdefault((c, num), {})[name] = np.asarray(x)
        names[(c, num)] = res[-1]
    for c, num, res, name, x in hydrogens:
        hy.setdefault((c, num), {})[name] = np.asarray(x)
    out = []
    for k in range(n):
        a, b = ("A", k + 1), ("B", 2 * n - k)
        out.append((base_fragment(by[a], hy.get(a, {}), names[a]), base_fragment(by[b], hy.get(b, {}), names[b])))
    return out


def _mol(*frags) -> Molecule:
    sym = [s for f in frags for s in f[0]]
    xyz = np.vstack([f[1] for f in frags])
    return Molecule(sym, xyz, 0, 1, [], "", [])


def interaction_gfn2(frags_a, frags_b, engine=None) -> float:
    """E(A+B) - E(A) - E(B) in kcal/mol with GFN2-xTB; frags_* are lists of (symbols, xyz)."""
    from ..engines.xtb import GFN2xTB
    eng = engine or GFN2xTB()
    e_ab = eng.compute(_mol(*frags_a, *frags_b)).energy
    e_a = eng.compute(_mol(*frags_a)).energy
    e_b = eng.compute(_mol(*frags_b)).energy
    return (e_ab - e_a - e_b) * HARTREE_KCAL


def interaction_dft(frags_a, frags_b, functional: str = "b3lyp-d3bj", basis: str = "def2-TZVP",
                    threads: int | None = None, memory: str = "4 GB") -> float:
    """Counterpoise-corrected interaction energy (kcal/mol) with Psi4."""
    import os

    import psi4

    from ..engines.dft import quiet_psi4
    quiet_psi4(psi4)
    psi4.set_memory(memory)
    psi4.set_num_threads(threads or max(1, (os.cpu_count() or 2) - 1))
    psi4.core.clean_options()
    psi4.set_options({"basis": basis, "scf_type": "df", "e_convergence": 1e-8, "d_convergence": 1e-8})

    def block(frags):
        return "\n".join(f"{s} {x:.8f} {y:.8f} {z:.8f}" for f in frags for s, (x, y, z) in zip(f[0], f[1]))
    geom = psi4.geometry(f"0 1\n{block(frags_a)}\n--\n0 1\n{block(frags_b)}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com")
    e = psi4.energy(functional, bsse_type="cp", molecule=geom)
    psi4.core.clean()
    return float(e) * HARTREE_KCAL


def helix_energies(atoms, hydrogens, n: int, method: str = "gfn2-xtb", on_step=None) -> dict:
    """Pairing energy of every base pair and stacking energy of every step (kcal/mol)."""
    fr = fragments(atoms, hydrogens, n)
    if method not in ("gfn2-xtb", "dft"):
        raise ValueError(f"unknown method {method!r}")
    calc = interaction_gfn2 if method == "gfn2-xtb" else interaction_dft
    pairs, steps = [], []
    total = 2 * n - 1
    for k, (a, b) in enumerate(fr):
        pairs.append(calc([a], [b]))
        if on_step:
            on_step({"task": "dna_energy", "step": k + 1, "steps": total, "stage": f"pairing energy of base pair {k + 1}"})
    for k in range(n - 1):
        steps.append(calc(list(fr[k]), list(fr[k + 1])))
        if on_step:
            on_step({"task": "dna_energy", "step": n + k + 1, "steps": total, "stage": f"stacking energy of step {k + 1}"})
    return {"pairs": pairs, "steps": steps, "method": method}
