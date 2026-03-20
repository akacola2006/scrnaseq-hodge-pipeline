"""
scRNAseq Hodge Pipeline — GPU-Accelerated PCA Engine
=====================================================
Randomized PCA via torch.pca_lowrank on GPU, with automatic
VRAM-aware fallback to CPU.
"""
import logging
from typing import Dict, Tuple

import numpy as np

try:
    import torch
except ImportError:
    torch = None

from . import config

logger = logging.getLogger(__name__)

_VRAM_SAFETY_BYTES = 2 * 1024**3


def _estimate_gpu_bytes(n_cells: int, n_genes: int) -> float:
    return n_cells * n_genes * 4 * 2.5


def pca_gpu(
    residuals: np.ndarray,
    k: int,
    device: str = "cuda:0",
    n_oversamples: int = None,
    n_iter: int = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """GPU-accelerated randomized PCA.

    Returns
    -------
    scores : np.ndarray, shape (n_cells, k)
    loadings : np.ndarray, shape (n_genes, k)
    singular_values : np.ndarray, shape (k,)
    """
    if n_oversamples is None:
        n_oversamples = config.PCA_N_OVERSAMPLES
    if n_iter is None:
        n_iter = config.PCA_N_ITER

    n_cells, n_genes = residuals.shape
    logger.info("PCA: %d cells x %d genes, k=%d, device=%s", n_cells, n_genes, k, device)

    max_k = min(n_cells, n_genes) - 1
    if k >= min(n_cells, n_genes):
        k = max(max_k, 1)
        logger.warning("  k clamped to %d", k)

    use_gpu = (
        config.USE_GPU
        and device.startswith("cuda")
        and torch is not None
        and torch.cuda.is_available()
    )

    if use_gpu:
        try:
            free_vram = torch.cuda.mem_get_info(0)[0]
        except Exception:
            free_vram = 16 * 1024**3

        needed = _estimate_gpu_bytes(n_cells, n_genes)
        if needed > free_vram - _VRAM_SAFETY_BYTES:
            logger.warning("  Matrix too large for GPU. Falling back to CPU.")
            use_gpu = False

    actual_device = device if use_gpu else "cpu"
    logger.info("  Using device: %s", actual_device)

    A = torch.from_numpy(residuals.astype(np.float32))
    col_means = A.mean(dim=0, keepdim=True)
    A = A - col_means

    if use_gpu:
        try:
            A = A.to(device)
        except RuntimeError:
            logger.warning("  GPU transfer failed, falling back to CPU")
            A = A.to("cpu")
            use_gpu = False

    q = k + n_oversamples
    try:
        U, S, V = torch.pca_lowrank(A, q=q, niter=n_iter)
    except RuntimeError as exc:
        if use_gpu and "out of memory" in str(exc).lower():
            logger.warning("  GPU OOM during PCA, retrying on CPU")
            torch.cuda.empty_cache()
            A = A.cpu()
            U, S, V = torch.pca_lowrank(A, q=q, niter=n_iter)
        else:
            raise

    U = U[:, :k]
    S = S[:k]
    V = V[:, :k]

    scores = (U * S.unsqueeze(0)).cpu().numpy().astype(np.float64)
    loadings = V.cpu().numpy().astype(np.float64)
    singular_values = S.cpu().numpy().astype(np.float64)

    del A, U, S, V
    if use_gpu:
        torch.cuda.empty_cache()

    logger.info("  PCA complete: scores %s, loadings %s", scores.shape, loadings.shape)
    return scores, loadings, singular_values


def extract_donor_scores(
    scores: np.ndarray,
    donor_ids: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Split PC scores by donor."""
    donor_ids = np.asarray(donor_ids, dtype=str)
    unique_donors = np.unique(donor_ids)
    result = {}
    for d in unique_donors:
        mask = donor_ids == d
        result[d] = scores[mask]
    logger.info(
        "  Extracted scores for %d donors (cells range: %d - %d)",
        len(result),
        min(v.shape[0] for v in result.values()),
        max(v.shape[0] for v in result.values()),
    )
    return result
