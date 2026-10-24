"""DNA Structure: helix geometry and melting."""
import numpy as np

from athanor_core.dna import melting as M
from athanor_core.dna.geometry import (kabsch, next_frame, pair_frames, pair_parameters, rot_axis,
                                       step_parameters, step_transform)


def rand_frame(rng):
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return Q * np.sign(np.linalg.det(Q)), rng.normal(size=3) * 10


def test_step_parameters_invert_exactly():
    rng = np.random.default_rng(0)
    for _ in range(500):
        R1, o1 = rand_frame(rng)
        p = np.array([rng.normal(), rng.normal(), 3 + rng.normal(), 8 * rng.normal(), 12 * rng.normal(),
                      rng.uniform(-60, 60)])
        R2, o2 = next_frame(R1, o1, p)
        q, Rm, om = step_parameters(R1, o1, R2, o2)
        assert np.allclose(q, p, atol=1e-9)
        assert np.allclose(Rm, R1 @ step_transform(p)[1], atol=1e-12) and np.allclose(om, (o1 + o2) / 2)
        (Ri, oi), (Rii, oii) = pair_frames(R1, o1, p)
        pp, Rb, ob = pair_parameters(Ri, oi, Rii, oii)
        assert np.allclose(pp, p, atol=1e-9) and np.allclose(Rb, R1) and np.allclose(ob, o1)


def test_pure_motions():
    I, o = np.eye(3), np.zeros(3)
    # a pure twist about z, a pure rise, a pure roll (rotation about y)
    q, _, _ = step_parameters(I, o, rot_axis([0, 0, 1], 36), np.array([0, 0, 3.4]))
    assert np.allclose(q, [0, 0, 3.4, 0, 0, 36])
    q, _, _ = step_parameters(I, o, rot_axis([0, 1, 0], 10), o)
    assert np.allclose(q, [0, 0, 0, 0, 10, 0])
    q, _, _ = step_parameters(I, o, rot_axis([1, 0, 0], 10), o)
    assert np.allclose(q, [0, 0, 0, 10, 0, 0])
    # ten B-DNA steps of 36 degrees make one turn
    R, oo = I, o
    for _ in range(10):
        R, oo = next_frame(R, oo, [0, 0, 3.38, 0, 0, 36])
    assert np.allclose(R, I, atol=1e-9) and np.allclose(oo, [0, 0, 33.8])


def test_kabsch():
    rng = np.random.default_rng(2)
    P = rng.normal(size=(9, 3))
    R, t = rand_frame(rng)
    R2, t2, rmsd = kabsch(P, P @ R.T + t)
    assert np.allclose(R2, R) and np.allclose(t2, t) and rmsd < 1e-10


def test_melting_matches_the_formula_and_its_parts():
    r = M.tm("TGACGGAGCGATTCAGC", ct=2e-6, na=1.0, salt="santalucia1998")
    # 16 stacks (TG GA AC CG GG GA AG GC CG GA AT TT TC CA AG GC: -137.7 kcal/mol) + initiation with a
    # terminal A.T (2.3) and a terminal G.C (0.1), SantaLucia 1998 Table 2
    assert abs(r["dH"] - (-137.7 + 2.3 + 0.1)) < 1e-9
    assert abs(r["tm"] - (r["dH"] * 1000 / (r["dS"] + M.R * np.log(2e-6 / 4)) - 273.15)) < 1e-9
    assert abs(M.curve("TGACGGAGCGATTCAGC", [r["tm"]], 2e-6)[0] - 0.5) < 1e-9
    # self-complementary: symmetry term and ct (not ct/4)
    s = M.tm("CGCGAATTCGCG", ct=1e-4)
    assert s["self_complementary"] and abs(M.curve("CGCGAATTCGCG", [s["tm"]], 1e-4)[0] - 0.5) < 1e-9
    # lower salt melts lower; higher concentration melts higher
    assert M.tm("TGACGGAGCGATTCAGC", 2e-6, 0.1)["tm"] < r["tm"] < M.tm("TGACGGAGCGATTCAGC", 2e-4, 1.0)["tm"]


def test_zipper_frays_ends_first():
    s = "GCGCAAATTTAAATTTGCGC"
    t = M.tm(s, 1e-6)["tm"]
    z = M.zipper(s, t + 5, 1e-6)["closed"]
    assert z[len(s) // 2] > z[0] and z[len(s) // 2] > z[-1]
    cold = M.zipper(s, 20, 1e-6)
    assert cold["fraction_duplex"] > 0.99 and min(cold["closed"]) > 0.9
    hot = M.zipper(s, 95, 1e-6)
    assert hot["fraction_duplex"] < 0.01
