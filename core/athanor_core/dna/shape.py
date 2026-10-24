"""DNA shape from sequence: a NumPy implementation of Deep DNAshape.

Deep DNAshape (Li, Chiu & Rohs, Nature Commun. 15, 1243, 2024; https://github.com/JinsenLi/deepDNAshape,
BSD-3-Clause) predicts the shape a free DNA duplex takes in solution - minor-groove width, propeller
twist, roll, helix twist and the other base-pair and step parameters, and their fluctuations - from its
sequence. It was trained on Monte Carlo simulations of thousands of sequences and refined on
experimental structures. Each base (or each step) is a node in a chain; seven message-passing layers
let every position see one more base pair of flanking sequence per layer.

Athanor runs the published trained weights (exported from the TensorFlow checkpoints in CI, see
tools/dnadata/dnashape_export.py) with this NumPy implementation of the same network, so no deep
learning framework ships with the app. CI checks that it gives the same numbers as the original.
"""
from __future__ import annotations

import itertools
import json
from functools import lru_cache

import numpy as np

from .structure import DATA

SHAPE = DATA / "dnashape"
BASES = ("T", "G", "C", "A")                       # one-hot order used by Deep DNAshape
INTRA = ("Shear", "Stretch", "Stagger", "Buckle", "ProT", "Opening", "MGW", "EP")
INTER = ("Shift", "Slide", "Rise", "Tilt", "Roll", "HelT")
XAXIS = {"Shift", "Tilt", "Shear", "Buckle"}       # change sign when read on the other strand
COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
UNITS = {"MGW": "A", "EP": "kT/e", "Shear": "A", "Stretch": "A", "Stagger": "A", "Shift": "A", "Slide": "A",
         "Rise": "A"}
NAMES = {"MGW": "minor-groove width", "ProT": "propeller twist", "Roll": "roll", "HelT": "helix twist",
         "EP": "electrostatic potential in the minor groove", "Shear": "shear", "Stretch": "stretch",
         "Stagger": "stagger", "Buckle": "buckle", "Opening": "opening", "Shift": "shift", "Slide": "slide",
         "Rise": "rise", "Tilt": "tilt"}


def _onehot():
    one = {b: np.eye(4)[i] for i, b in enumerate(BASES)}
    one["N"] = np.full(4, 0.25)
    di = {}
    for i, (a, b) in enumerate(itertools.product(BASES, repeat=2)):
        v = np.zeros(16)
        v[i] = 1.0
        di[(a, b)] = v
    di[("N", "N")] = np.full(16, 1 / 16)
    for a in BASES:
        di[(a, "N")] = sum(di[(a, b)] for b in BASES) / 4
        di[("N", a)] = sum(di[(b, a)] for b in BASES) / 4
    return one, di


ONE, DI = _onehot()


@lru_cache(maxsize=None)
def index() -> dict:
    return json.loads((SHAPE / "index.json").read_text())


@lru_cache(maxsize=None)
def weights(feature: str) -> dict:
    """{name: array} for one feature, from the raw float32 file and its index."""
    meta = index()["features"][feature]
    raw = np.fromfile(SHAPE / meta["file"], dtype="<f4")
    out = {}
    for name, (off, shape) in meta["arrays"].items():
        n = int(np.prod(shape))
        out[name] = raw[off:off + n].reshape(shape).astype(np.float64)
    return out


def available() -> list[str]:
    try:
        return list(index()["features"])
    except FileNotFoundError:
        return []


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def forward(W: dict, x: np.ndarray, layers: int = 7) -> np.ndarray:
    """The network on one chain of nodes x [L, F]; returns the node outputs after each layer [L, layers+1]."""
    h = x @ W["selfconv/kernel"][0] + W["selfconv/bias"]
    out = [h.mean(1)]
    L = len(h)
    prev = np.r_[0, np.arange(L - 1)]
    nxt = np.r_[np.arange(1, L), L - 1]
    for k in range(layers):
        p = f"mp_layers/{k}/"
        hp, hn = h[prev], h[nxt]
        a = (hn @ W[p + "Weight1"] + hp @ W[p + "Weight2"] + W[p + "Bias1"]
             + hn @ W[p + "Weight1_term2"] + hp @ W[p + "Weight2_term2"])
        a = a * h + W[p + "BiasAll"]
        a = np.maximum(a, 0.0)
        a = (a - W[p + "bn/moving_mean"]) / np.sqrt(W[p + "bn/moving_variance"] + 1e-3) * W[p + "bn/gamma"] + W[p + "bn/beta"]
        K, R, b = W[p + "gru_layer/kernel"], W[p + "gru_layer/recurrent_kernel"], W[p + "gru_layer/bias"]
        xi = a @ K + b[0]
        hr = h @ R + b[1]
        n = h.shape[1]
        z = _sigmoid(xi[:, :n] + hr[:, :n])
        r = _sigmoid(xi[:, n:2 * n] + hr[:, n:2 * n])
        hh = np.tanh(xi[:, 2 * n:] + r * hr[:, 2 * n:])
        h = z * h + (1.0 - z) * hh
        out.append(h.mean(1))
    return np.stack(out, axis=1)


def _rescale(feature, y):
    p = index()["params"][feature]
    if p["method"] == "standard":
        return y * p["std"] + p["mean"]
    return y * p["percentile_range"] + p["median"]


def _encode(seq, inter):
    if inter:
        return np.array([DI[(seq[i], seq[i + 1])] for i in range(len(seq) - 1)])
    return np.array([ONE[b] for b in seq])


def predict(feature: str, seq: str, layer: int = 4) -> np.ndarray:
    """Deep DNAshape's prediction for every base pair (intra-base-pair features, MGW, EP) or every step
    (inter-base-pair features), with ``layer`` base pairs of flanking sequence considered (0-7; the
    authors' default is 4). The sequence is read on both strands and the two readings averaged, as the
    original does."""
    seq = "".join(seq.split()).upper().replace("U", "T")
    if any(b not in "ACGTN" for b in seq):
        raise ValueError("sequence must be A, C, G, T or N")
    base = feature.replace("-FL", "")
    inter = base in INTER
    W = weights(feature)
    s = "NN" + seq + "NN"
    r = "".join(COMP[b] for b in reversed(s))
    fwd = _rescale(feature, forward(W, _encode(s, inter)))
    rev = _rescale(feature, forward(W, _encode(r, inter)))
    if feature in XAXIS:
        rev = -rev
    y = (fwd + rev[::-1]) / 2
    return y[2:-2, layer]


def profile(seq: str, features=("MGW", "ProT", "Roll", "HelT"), layer: int = 4) -> dict:
    """Several features at once, as plain lists, for the app."""
    return {f: predict(f, seq, layer).tolist() for f in features if f in available()}
