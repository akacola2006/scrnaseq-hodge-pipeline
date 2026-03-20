"""
scRNAseq Hodge Pipeline — Lane B: Upstream PC & Gene Extraction
================================================================
Identifies upstream PCs via d_corr contribution decomposition
and extracts gene candidates from PC loadings.
"""
import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from . import config
from .spd import spd_log

logger = logging.getLogger(__name__)


def _standardize(x: np.ndarray) -> np.ndarray:
    """Z-score standardize a 1-D array."""
    s = x.std()
    if s < 1e-15:
        return np.zeros_like(x)
    return (x - x.mean()) / s


def identify_upstream_pcs(
    corr_matrices_by_window: Dict[int, np.ndarray],
    upstream_window: int,
    n_top: int = None,
) -> Dict[str, Any]:
    """Identify upstream PCs via d_corr contribution decomposition.

    For the upstream transition (window w* -> w*+1), computes element-wise
    change in log-correlation matrix and scores each PC on diagonal and row-sum.
    """
    if n_top is None:
        n_top = config.UPSTREAM_PC_TOP_N

    w_star = upstream_window
    w_next = w_star + 1

    if w_star not in corr_matrices_by_window:
        raise KeyError(f"Window {w_star} not found")
    if w_next not in corr_matrices_by_window:
        raise KeyError(f"Window {w_next} not found")

    R_w = corr_matrices_by_window[w_star]
    R_next = corr_matrices_by_window[w_next]
    k = R_w.shape[0]

    logger.info("Lane B: d_corr decomposition, window %d -> %d, k=%d", w_star, w_next, k)

    log_R_w = spd_log(R_w)
    log_R_next = spd_log(R_next)
    delta = log_R_next - log_R_w

    # Diagonal score: |Delta_ii|
    diag_scores = np.abs(np.diag(delta))

    # Row-sum score: sum_j |Delta_ij|
    row_scores = np.abs(delta).sum(axis=1)

    # Dual score: z(diag) + z(row)
    dual_scores = _standardize(diag_scores) + _standardize(row_scores)

    # Top N PCs
    top_idx = np.argsort(-dual_scores)[:n_top]

    logger.info("Top %d upstream PCs: %s", n_top, top_idx.tolist())

    return {
        "pc_indices": top_idx,
        "dual_scores": dual_scores,
        "diag_scores": diag_scores,
        "row_scores": row_scores,
        "delta_matrix": delta,
    }


def extract_gene_set(
    loadings: np.ndarray,
    pc_indices: np.ndarray,
    gene_names: List[str],
    top_n_genes: int = None,
) -> Tuple[List[str], np.ndarray]:
    """Extract top genes from upstream PC loadings.

    Parameters
    ----------
    loadings : np.ndarray, shape (n_genes, k)
        PCA loadings matrix.
    pc_indices : np.ndarray
        Indices of upstream PCs.
    gene_names : list of str
        Gene identifiers.
    top_n_genes : int, optional
        Top genes per PC. Defaults to config.GENE_LOADING_TOP_N.

    Returns
    -------
    gene_set : list of str
        Unique gene names from all upstream PCs.
    gene_indices : np.ndarray
        Indices into the original gene list.
    """
    if top_n_genes is None:
        top_n_genes = config.GENE_LOADING_TOP_N

    all_gene_indices = set()

    for pc_idx in pc_indices:
        if pc_idx >= loadings.shape[1]:
            logger.warning("PC index %d exceeds loadings shape; skipping", pc_idx)
            continue

        pc_loadings = np.abs(loadings[:, pc_idx])
        top_genes = np.argsort(-pc_loadings)[:top_n_genes]
        all_gene_indices.update(top_genes.tolist())

    gene_indices = np.array(sorted(all_gene_indices), dtype=np.int64)
    gene_set = [gene_names[i] for i in gene_indices]

    logger.info(
        "Gene set: %d unique genes from %d PCs x %d genes/PC",
        len(gene_set), len(pc_indices), top_n_genes,
    )

    return gene_set, gene_indices
