"""Principal components of reference genotypes, and projection of new samples.

As in EIGENSOFT (Patterson, Price & Reich, PLoS Genet. 2, e190, 2006), each
SNP is centred by twice its reference allele frequency p and scaled by
sqrt(p(1 - p)). Components come from the eigendecomposition of the
reference samples' genetic relationship matrix. Samples with missing data
(ancient genomes, genotyping chips) are placed by least squares on the SNPs
they have, as smartpca's lsqproject does (Patterson in Lazaridis et al.,
Nature 513, 409, 2014).

With Z = U S V' for the reference samples, a reference sample's coordinates
are its row of U S; a new sample with normalised genotypes z on SNPs O gets
x = argmin ||z_O - V_O x||, which equals U S for a reference sample with no
missing data.
"""
from __future__ import annotations

import numpy as np


def normalise(dosage: np.ndarray, p: np.ndarray, scale: np.ndarray) -> np.ndarray:
    """dosage: alt-allele counts (0, 1, 2; NaN or negative = missing)."""
    d = np.asarray(dosage, dtype=float)
    d = np.where(d < 0, np.nan, d)
    return (d - 2 * p) / scale


def project(z: np.ndarray, loadings: np.ndarray, min_snps: int = 1000) -> tuple[np.ndarray, int]:
    """Least-squares coordinates of one sample from normalised genotypes z
    (NaN = missing) and SNP loadings V (S x k)."""
    ok = ~np.isnan(z)
    n = int(ok.sum())
    if n < min_snps:
        raise ValueError(f"only {n} SNPs overlap the PCA reference; at least {min_snps} are needed")
    V = loadings[ok].astype(float)
    A = V.T @ V
    b = V.T @ z[ok]
    return np.linalg.solve(A, b), n


def nearest(coords: np.ndarray, points: np.ndarray, weights: np.ndarray | None = None, k: int = 10):
    """Indices and distances of the k points nearest to coords (optionally
    weighting each axis)."""
    w = np.ones(points.shape[1]) if weights is None else weights
    d = np.sqrt((((points - coords[None, :]) ** 2) * w[None, :]).sum(1))
    idx = np.argsort(d)[:k]
    return idx, d[idx]
