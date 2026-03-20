"""
scRNAseq Hodge Pipeline — SPD Covariance & Log-Euclidean Operations
====================================================================
Per-donor covariance estimation (Ledoit-Wolf), matrix logarithm/exponential
on the SPD manifold, Log-Euclidean distances and means, and the
d_cov / d_corr / d_var decomposition.
"""
import logging
from typing import Dict, List, Tuple

import numpy as np
import scipy.linalg
from sklearn.covariance import LedoitWolf

from . import config

logger = logging.getLogger(__name__)


def estimate_spd(scores: np.ndarray) -> np.ndarray:
    """Ledoit-Wolf shrinkage covariance estimation."""
    n, k = scores.shape
    if n < 2:
        raise ValueError(f"Need >= 2 observations, got {n}")
    lw = LedoitWolf()
    lw.fit(scores)
    Sigma = lw.covariance_
    Sigma = (Sigma + Sigma.T) / 2.0
    return Sigma


def cov_to_corr(cov: np.ndarray, eps: float = 1e-12) -> Tuple[np.ndarray, np.ndarray]:
    """Convert covariance matrix to correlation matrix."""
    std_devs = np.sqrt(np.maximum(np.diag(cov), eps))
    D_inv = np.diag(1.0 / std_devs)
    corr = D_inv @ cov @ D_inv
    corr = (corr + corr.T) / 2.0
    np.fill_diagonal(corr, 1.0)
    return corr, std_devs


def spd_log(A: np.ndarray) -> np.ndarray:
    """Matrix logarithm for a symmetric positive-definite matrix."""
    A_sym = (A + A.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(A_sym)
    eigvals = np.maximum(eigvals, 1e-15)
    log_A = eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return (log_A + log_A.T) / 2.0


def spd_exp(A: np.ndarray) -> np.ndarray:
    """Matrix exponential for a symmetric matrix."""
    A_sym = (A + A.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(A_sym)
    exp_A = eigvecs @ np.diag(np.exp(eigvals)) @ eigvecs.T
    return (exp_A + exp_A.T) / 2.0


def log_euclidean_mean(matrices: List[np.ndarray]) -> np.ndarray:
    """Log-Euclidean Frechet mean of SPD matrices."""
    if not matrices:
        raise ValueError("Cannot compute mean of empty list")
    log_sum = spd_log(matrices[0]).copy()
    for M in matrices[1:]:
        log_sum += spd_log(M)
    log_mean = log_sum / len(matrices)
    return spd_exp(log_mean)


def log_euclidean_distance(A: np.ndarray, B: np.ndarray) -> float:
    """Log-Euclidean distance between two SPD matrices."""
    log_A = spd_log(A)
    log_B = spd_log(B)
    return float(np.linalg.norm(log_A - log_B, "fro"))


def decompose_distance(
    cov_A: np.ndarray,
    cov_B: np.ndarray,
    eps: float = 1e-12,
) -> Dict[str, float]:
    """Decompose Log-Euclidean covariance distance into d_cov, d_corr, d_var."""
    R_A, std_A = cov_to_corr(cov_A, eps=eps)
    R_B, std_B = cov_to_corr(cov_B, eps=eps)

    log_R_A = spd_log(R_A)
    log_R_B = spd_log(R_B)
    d_corr = float(np.linalg.norm(log_R_A - log_R_B, "fro"))
    d_var = float(np.linalg.norm(np.log(std_A) - np.log(std_B), 2))
    d_cov = log_euclidean_distance(cov_A, cov_B)

    return {"d_cov": d_cov, "d_corr": d_corr, "d_var": d_var}


def vectorize_upper_tri(matrix: np.ndarray) -> np.ndarray:
    """Extract upper triangle (including diagonal) as a flat vector."""
    k = matrix.shape[0]
    idx = np.triu_indices(k, k=0)
    return matrix[idx]


def batch_spd_log_gpu(
    matrices: np.ndarray,
    device: str = "cuda:0",
) -> np.ndarray:
    """Batch matrix logarithm on GPU via eigendecomposition."""
    import torch

    use_gpu = (
        config.USE_GPU
        and device.startswith("cuda")
        and torch.cuda.is_available()
    )
    actual_device = device if use_gpu else "cpu"

    M = torch.from_numpy(matrices.astype(np.float64)).to(actual_device)
    M = (M + M.transpose(-1, -2)) / 2.0
    eigvals, eigvecs = torch.linalg.eigh(M)
    eigvals = torch.clamp(eigvals, min=1e-15)
    log_eigvals = torch.log(eigvals)
    log_M = eigvecs @ torch.diag_embed(log_eigvals) @ eigvecs.transpose(-1, -2)
    log_M = (log_M + log_M.transpose(-1, -2)) / 2.0

    result = log_M.cpu().numpy()

    if use_gpu:
        del M, eigvals, eigvecs, log_eigvals, log_M
        torch.cuda.empty_cache()

    return result
