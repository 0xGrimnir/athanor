import numpy as np

from athanor_core.adna.damage import damage_profiles
from athanor_core.adna.truth import truth_by_reference_position


def test_cs_walk_maps_reference_positions_to_truth():
    # reference: ACGTACGTAC (10), truth: ACGAACGTTTAC -> sub T>A at ref 3, insertion TT after ref 7, ...
    truth = "ACGAACGTTTAC"
    paf = "\t".join(["truth", "12", "0", "12", "+", "ref", "10", "0", "10", "9", "12", "60",
                     "tp:A:P", "cs:Z::3*ta:4+tt:2"])
    m = truth_by_reference_position(paf, 10, truth)
    assert "".join(m) == "ACGAACGTAC"
    paf_del = "\t".join(["truth", "8", "0", "8", "+", "ref", "10", "0", "10", "8", "10", "60",
                         "tp:A:P", "cs:Z::4-ac:4"])
    assert "".join(truth_by_reference_position(paf_del, 10, "ACGTGTAC")) == "ACGT--GTAC"


def test_damage_profile_subtracts_background(tmp_path):
    ct = [0.30, 0.18, 0.10, 0.06] + [0.02] * 21
    (tmp_path / "5pCtoT_freq.txt").write_text("pos\t5pC>T\n" + "".join(f"{i+1}\t{v}\n" for i, v in enumerate(ct)))
    (tmp_path / "3pGtoA_freq.txt").write_text("pos\t3pG>A\n" + "".join(f"{i+1}\t{v}\n" for i, v in enumerate(ct)))
    d = damage_profiles(tmp_path)
    assert d["background_CtoT"] == 0.02
    np.testing.assert_allclose(d["delta5"][:4], [0.28, 0.16, 0.08, 0.04])
    assert d["delta5"][-1] == 0
