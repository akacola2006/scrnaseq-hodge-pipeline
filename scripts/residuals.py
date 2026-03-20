"""
scRNAseq Hodge Pipeline — GPU-Accelerated Regression Residuals
===============================================================
Computes residuals from the model: y ~ donor + nUMI + pct_mito
on log1p(CPM)-normalized expression, with GPU-batched least squares.

GPU strategy:
  - Design matrix X is built on CPU (small: n_cells x p, p ~ 50).
  - Expression matrix Y is transferred to GPU in gene-batches.
  - (X^T X)^{-1} X^T is precomputed once on GPU.
  - Falls back to CPU on OOM or when USE_GPU=False.
"""
import gc
import logging
from typing import List, Optional, Tuple

import anndata as ad
import numpy as np
import torch

from . import config

logger = logging.getLogger(__name__)

_GPU_GENE_BATCH = 4000


# ---------------------------------------------------------------------------
# 1. Normalization
# ---------------------------------------------------------------------------

def normalize_log1p_cpm(adata: ad.AnnData) -> np.ndarray:
    """Convert raw counts to log1p(CPM)."""
    X = np.asarray(adata.X, dtype=np.float32)
    totals = X.sum(axis=1, keepdims=True)
    totals = np.where(totals > 0, totals, 1.0)
    cpm = X / totals * 1e6
    return np.log1p(cpm).astype(np.float32)


# ---------------------------------------------------------------------------
# 2. Gene QC -- detection-rate filter
# ---------------------------------------------------------------------------

def filter_genes_by_detection(
    expr: np.ndarray,
    min_frac: float = None,
) -> np.ndarray:
    """Return boolean mask of genes detected in >= min_frac of cells."""
    if min_frac is None:
        min_frac = config.GENE_DETECTION_MIN_FRAC
    n_cells = expr.shape[0]
    detection_rate = (expr > 0).sum(axis=0) / n_cells
    return detection_rate >= min_frac


# ---------------------------------------------------------------------------
# 3. HVG selection
# ---------------------------------------------------------------------------

def select_hvg(
    expr_log: np.ndarray,
    n_hvg: Optional[int],
) -> np.ndarray:
    """Select highly variable genes by variance on log1p-CPM expression."""
    n_genes = expr_log.shape[1]
    if n_hvg is None or n_hvg >= n_genes:
        return np.ones(n_genes, dtype=bool)

    gene_var = expr_log.var(axis=0)
    threshold = np.partition(gene_var, -n_hvg)[-n_hvg]
    mask = gene_var >= threshold
    if mask.sum() > n_hvg:
        indices = np.where(mask)[0]
        sorted_idx = indices[np.argsort(-gene_var[indices])]
        mask = np.zeros(n_genes, dtype=bool)
        mask[sorted_idx[:n_hvg]] = True
    return mask


# ---------------------------------------------------------------------------
# 4. Design-matrix construction
# ---------------------------------------------------------------------------

def _build_design_matrix(obs) -> np.ndarray:
    """Build design matrix: intercept + donor (one-hot, drop-first) + nUMI + pct_mito."""
    donors = obs["donor_id"].astype(str).values
    unique_donors = sorted(set(donors))

    n_cells = len(donors)
    n_donor_cols = max(len(unique_donors) - 1, 0)
    donor_dummies = np.zeros((n_cells, n_donor_cols), dtype=np.float32)
    donor_map = {d: i for i, d in enumerate(unique_donors[1:])}
    for cell_idx, d in enumerate(donors):
        if d in donor_map:
            donor_dummies[cell_idx, donor_map[d]] = 1.0

    nUMI = obs["nUMI"].values.astype(np.float32).reshape(-1, 1)
    pct_mito = obs["pct_mito"].values.astype(np.float32).reshape(-1, 1)
    intercept = np.ones((n_cells, 1), dtype=np.float32)

    return np.hstack([intercept, donor_dummies, nUMI, pct_mito])


# ---------------------------------------------------------------------------
# 5. GPU-accelerated batch residuals
# ---------------------------------------------------------------------------

def compute_residuals_gpu(
    expr_log: np.ndarray,
    donor_ids: np.ndarray,
    nUMI: np.ndarray,
    pct_mito: np.ndarray,
    device: str = "cuda:0",
) -> np.ndarray:
    """GPU-accelerated batch OLS regression and residual computation."""
    import pandas as pd

    n_cells, n_genes = expr_log.shape

    obs_df = pd.DataFrame({
        "donor_id": donor_ids.astype(str),
        "nUMI": nUMI.astype(np.float32),
        "pct_mito": pct_mito.astype(np.float32),
    })
    X_design = _build_design_matrix(obs_df)
    p = X_design.shape[1]
    logger.info("  Design matrix: %d cells x %d covariates", n_cells, p)

    residuals = np.empty((n_cells, n_genes), dtype=np.float32)

    X_t = torch.from_numpy(X_design).to(device=device, dtype=torch.float32)
    XtX = X_t.T @ X_t
    XtX += torch.eye(p, device=device, dtype=torch.float32) * 1e-6
    XtX_inv = torch.linalg.inv(XtX)
    H_coeff = XtX_inv @ X_t.T

    gene_batch = _GPU_GENE_BATCH
    for start in range(0, n_genes, gene_batch):
        end = min(start + gene_batch, n_genes)
        Y_batch = torch.from_numpy(
            expr_log[:, start:end]
        ).to(device=device, dtype=torch.float32)

        beta = H_coeff @ Y_batch
        Y_hat = X_t @ beta
        R_batch = Y_batch - Y_hat

        residuals[:, start:end] = R_batch.cpu().numpy()

        del Y_batch, beta, Y_hat, R_batch

    del X_t, XtX, XtX_inv, H_coeff
    torch.cuda.empty_cache()

    return residuals


# ---------------------------------------------------------------------------
# 6. CPU fallback
# ---------------------------------------------------------------------------

def _cpu_residuals(X_design: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """CPU fallback: compute residuals via numpy least squares."""
    beta, _, _, _ = np.linalg.lstsq(X_design, Y, rcond=None)
    Y_hat = X_design @ beta
    return (Y - Y_hat).astype(np.float32)


# ---------------------------------------------------------------------------
# 7. Full pipeline entry point
# ---------------------------------------------------------------------------

def compute_residuals(
    adata: ad.AnnData,
    hvg_n: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Full residual-computation pipeline for one cell type.

    Steps:
      1. log1p(CPM) normalization
      2. Gene detection-rate filter
      3. Optional HVG selection
      4. Build design matrix (donor + nUMI + pct_mito)
      5. GPU-batched OLS residuals (CPU fallback on OOM / no GPU)

    Returns
    -------
    residuals : np.ndarray, shape (n_cells, n_genes_kept)
    gene_mask : np.ndarray of bool, shape (n_genes_original,)
    gene_names : list of str
    """
    device = config.GPU_DEVICE
    logger.info(
        "Computing residuals: %d cells x %d genes, device=%s, hvg_n=%s",
        adata.n_obs, adata.n_vars, device, hvg_n,
    )

    X_raw = np.asarray(adata.X, dtype=np.float32)
    expr_log = normalize_log1p_cpm(adata)
    logger.info("  Normalization done: log1p(CPM)")

    det_mask = filter_genes_by_detection(X_raw, min_frac=config.GENE_DETECTION_MIN_FRAC)
    n_det = int(det_mask.sum())
    logger.info(
        "  Gene detection filter: %d / %d pass (>= %.1f%% cells)",
        n_det, X_raw.shape[1], config.GENE_DETECTION_MIN_FRAC * 100,
    )

    if hvg_n is not None:
        hvg_mask_sub = select_hvg(expr_log[:, det_mask], n_hvg=hvg_n)
        gene_mask = np.zeros(adata.n_vars, dtype=bool)
        det_indices = np.where(det_mask)[0]
        gene_mask[det_indices[hvg_mask_sub]] = True
        logger.info(
            "  HVG filter: requested %d, kept %d",
            hvg_n, int(gene_mask.sum()),
        )
    else:
        gene_mask = det_mask

    n_genes_kept = int(gene_mask.sum())
    if n_genes_kept == 0:
        raise ValueError("No genes passed filtering.")

    expr_filtered = expr_log[:, gene_mask]
    gene_names = list(adata.var_names[gene_mask])

    use_gpu = (
        config.USE_GPU
        and device.startswith("cuda")
        and torch.cuda.is_available()
    )

    if use_gpu:
        try:
            logger.info("  Running GPU-batched OLS on %s (%d genes)", device, n_genes_kept)
            residuals = compute_residuals_gpu(
                expr_log=expr_filtered,
                donor_ids=adata.obs["donor_id"].values,
                nUMI=adata.obs["nUMI"].values,
                pct_mito=adata.obs["pct_mito"].values,
                device=device,
            )
            logger.info("  GPU residuals complete")
        except (torch.cuda.OutOfMemoryError, RuntimeError) as exc:
            logger.warning("  GPU OOM (%s), falling back to CPU", exc)
            torch.cuda.empty_cache()
            X_design = _build_design_matrix(adata.obs)
            residuals = _cpu_residuals(X_design, expr_filtered)
            logger.info("  CPU residuals complete (fallback)")
    else:
        logger.info("  Running CPU OLS (%d genes)", n_genes_kept)
        X_design = _build_design_matrix(adata.obs)
        residuals = _cpu_residuals(X_design, expr_filtered)
        logger.info("  CPU residuals complete")

    return residuals, gene_mask, gene_names
