"""
ALS IDS Stage-1 Rebuild -- Background Covariance Whitening
==========================================================
Implements ZCA whitening (Σ_ctrl^{-1/2}) to remove background covariance
before gene-level Hodge decomposition, following the IDS framework
(CellNavi-IDS, Kaneko et al. 2026).

The whitening matrix W is estimated once from control (PN) donors and
applied to ALL cells before per-donor log-correlation estimation.
This isolates perturbation-specific (ALS-specific) covariance changes
from background gene-gene correlation structure.

Key design decisions:
  - W is estimated from donor-split pseudobulk means (not raw cells)
    to match the donor-level granularity of the Hodge pipeline.
  - Ledoit-Wolf shrinkage ensures W is well-conditioned even when p > n.
  - ZCA form (W = V D^{-1/2} V^T) preserves gene-space interpretability.
  - Whitening is applied once; bootstrap iterations do NOT re-estimate W.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# =====================================================================
# 1. Donor Split Utilities
# =====================================================================

def split_indices_per_donor(
    donor_ids: np.ndarray,
    target_donors: Sequence[str],
    n_splits: int = 5,
) -> Dict[str, List[np.ndarray]]:
    """Split cell indices for each control donor into n_splits groups.

    Parameters
    ----------
    donor_ids : (n_cells,) array of donor IDs
    target_donors : donor IDs to split
    n_splits : number of splits per donor

    Returns
    -------
    {donor_id: [idx_split_0, idx_split_1, ...]}
    """
    result = {}
    for did in target_donors:
        mask = np.where(donor_ids == did)[0]
        if len(mask) == 0:
            continue
        # Deterministic split (no shuffle for reproducibility)
        splits = np.array_split(mask, min(n_splits, len(mask)))
        result[did] = [s for s in splits if len(s) > 0]
    return result


# =====================================================================
# 2. Build Control Whitener
# =====================================================================

def build_control_whitener(
    residuals: np.ndarray,
    donor_ids: np.ndarray,
    control_donor_ids: Sequence[str],
    gene_indices: np.ndarray,
    n_splits: int = 5,
    eigen_min: float = 1e-10,
    eps: float = 1e-12,
    max_cond: float = 1e6,
    reduce_p: Optional[int] = None,
    seed: Optional[int] = None,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Build ZCA whitening matrix from control donors.

    Pipeline:
      1. Extract control donor cells → residuals[:, gene_indices]
      2. Split each donor into n_splits equal parts
      3. Compute mean residual per split → one sample vector
      4. Estimate covariance via Ledoit-Wolf on split-mean samples
      5. Eigendecompose → clip eigenvalues → ZCA whitening matrix

    Parameters
    ----------
    residuals : (n_cells, n_all_genes) full residual matrix
    donor_ids : (n_cells,) donor ID per cell
    control_donor_ids : donor IDs classified as control/PN
    gene_indices : column indices of candidate genes
    n_splits : splits per donor for pseudobulk samples
    eigen_min : minimum eigenvalue (clipped)
    eps : numerical safety margin
    max_cond : maximum condition number before fallback
    reduce_p : if set, reduce gene count to this value for stability
    seed : random seed for gene selection in reduce_p

    Returns
    -------
    W : (p, p) whitening matrix
    meta : metadata dict
    """
    from sklearn.covariance import LedoitWolf

    p_orig = len(gene_indices)
    gene_idx_used = gene_indices.copy()

    # Optional gene reduction for p >> n stability
    if reduce_p is not None and reduce_p < p_orig:
        gene_idx_used = _reduce_genes_by_variance(
            residuals, donor_ids, control_donor_ids,
            gene_indices, reduce_p, seed=seed,
        )
        logger.info("  reduce_p: %d -> %d genes", p_orig, len(gene_idx_used))

    p = len(gene_idx_used)

    # Step 1-3: Build split-mean samples
    splits = split_indices_per_donor(donor_ids, control_donor_ids, n_splits)

    if len(splits) == 0:
        raise RuntimeError("No control donors found for whitening")

    samples = []
    for did, split_indices in splits.items():
        for s_idx in split_indices:
            split_resid = residuals[s_idx][:, gene_idx_used]
            samples.append(split_resid.mean(axis=0))

    X_ctrl = np.array(samples)  # (n_samples, p)
    n_samples = X_ctrl.shape[0]
    logger.info(
        "  Whitening: %d control donors, %d splits, %d samples, %d genes",
        len(splits), n_splits, n_samples, p,
    )

    if n_samples < 3:
        raise RuntimeError(
            f"Too few whitening samples: {n_samples} (need >= 3)"
        )

    # Step 4: Ledoit-Wolf covariance
    lw = LedoitWolf()
    lw.fit(X_ctrl)
    sigma = lw.covariance_
    shrinkage_alpha = float(lw.shrinkage_)

    # Step 5: Eigendecompose
    eigvals, eigvecs = np.linalg.eigh(sigma)
    eigvals_raw = eigvals.copy()
    n_clipped = int((eigvals < eigen_min).sum())
    eigvals = np.clip(eigvals, eigen_min, None)

    cond = float(eigvals.max() / eigvals.min())
    logger.info(
        "  Sigma: cond=%.2e, eigmin=%.2e, eigmax=%.2e, clipped=%d/%d, alpha=%.4f",
        cond, eigvals.min(), eigvals.max(), n_clipped, p, shrinkage_alpha,
    )

    if cond > max_cond:
        raise RuntimeError(
            f"Whitening matrix ill-conditioned: cond={cond:.2e} > {max_cond:.2e}. "
            f"Consider --whiten_fallback reduce_p or reducing gene set."
        )

    # ZCA: W = V @ diag(1/sqrt(eig)) @ V^T
    W = eigvecs @ np.diag(1.0 / np.sqrt(eigvals + eps)) @ eigvecs.T

    # Sanity check: control split covariance after whitening should be ~I
    X_whitened = X_ctrl @ W
    cov_whitened = np.cov(X_whitened, rowvar=False)
    frob_to_I = float(np.linalg.norm(cov_whitened - np.eye(p), 'fro'))
    frob_I_norm = frob_to_I / np.sqrt(p)
    logger.info("  Sanity: ||Cov_whitened - I||_F / sqrt(p) = %.4f", frob_I_norm)

    meta = {
        "n_ctrl_donors": len(splits),
        "n_ctrl_splits": n_splits,
        "n_samples": n_samples,
        "p": p,
        "p_original": p_orig,
        "cond": cond,
        "eigmin": float(eigvals.min()),
        "eigmax": float(eigvals.max()),
        "clipped_ratio": n_clipped / p if p > 0 else 0,
        "shrinkage_alpha": shrinkage_alpha,
        "whiten_control_mode": "custom",
        "fallback_used": reduce_p is not None and reduce_p < p_orig,
        "reduce_p": reduce_p,
        "seed": seed,
        "frob_whitened_to_I": frob_to_I,
        "frob_whitened_to_I_norm": frob_I_norm,
        "gene_indices_used": gene_idx_used.tolist(),
        "control_donor_ids": list(control_donor_ids),
    }

    return W, meta


# =====================================================================
# 3. Apply Whitener
# =====================================================================

def apply_whitener(
    residuals: np.ndarray,
    W: np.ndarray,
    gene_indices: np.ndarray,
) -> np.ndarray:
    """Apply whitening to residuals in-place (gene subset).

    Parameters
    ----------
    residuals : (n_cells, n_all_genes) full residual matrix
    W : (p, p) whitening matrix
    gene_indices : column indices of the p genes

    Returns
    -------
    whitened_residuals : (n_cells, n_all_genes) with whitened gene columns
    """
    p = len(gene_indices)
    assert W.shape == (p, p), f"W shape {W.shape} != ({p}, {p})"
    assert residuals.ndim == 2, f"residuals must be 2D, got {residuals.ndim}D"

    # Copy to avoid modifying the original
    out = residuals.copy()
    sub = out[:, gene_indices]  # (n_cells, p)

    # Check for NaN/Inf before
    n_nan_before = int(np.isnan(sub).sum())
    n_inf_before = int(np.isinf(sub).sum())
    if n_nan_before > 0 or n_inf_before > 0:
        logger.warning(
            "  apply_whitener: %d NaN, %d Inf in input",
            n_nan_before, n_inf_before,
        )

    # Apply: residuals @ W (NOT W @ residuals)
    whitened = sub @ W  # (n_cells, p)

    # Check for NaN/Inf after
    n_nan_after = int(np.isnan(whitened).sum())
    n_inf_after = int(np.isinf(whitened).sum())
    if n_nan_after > 0 or n_inf_after > 0:
        logger.warning(
            "  apply_whitener: %d NaN, %d Inf in output (PROBLEM)",
            n_nan_after, n_inf_after,
        )

    out[:, gene_indices] = whitened

    logger.info(
        "  apply_whitener: (%d, %d) @ (%d, %d) -> applied. "
        "NaN before=%d after=%d, Inf before=%d after=%d",
        sub.shape[0], sub.shape[1], W.shape[0], W.shape[1],
        n_nan_before, n_nan_after, n_inf_before, n_inf_after,
    )

    return out


# =====================================================================
# 4. Control Donor Identification
# =====================================================================

def identify_control_donors(
    pt_df: "pd.DataFrame",
    donor_ids: np.ndarray,
    mode: str = "pn_only",
    condition_col: Optional[str] = None,
    sample_info: Optional["pd.DataFrame"] = None,
) -> List[str]:
    """Identify control donors for whitening.

    Parameters
    ----------
    pt_df : DataFrame with donor_id, window, (optionally pt_value)
    donor_ids : all donor IDs in the dataset
    mode : "pn_only" or "pn_plus_early"
    condition_col : column name in sample_info with PN/ALS labels
    sample_info : sample metadata DataFrame

    Returns
    -------
    List of control donor IDs
    """
    import pandas as pd

    unique_donors = sorted(set(donor_ids))

    # Try to find condition labels
    donor_to_condition = {}

    if sample_info is not None:
        # Auto-detect condition column
        if condition_col is None:
            for col in ["subtype", "diagnosis", "condition", "Condition", "Group"]:
                if col in sample_info.columns:
                    condition_col = col
                    break

        if condition_col and condition_col in sample_info.columns:
            donor_to_condition = dict(
                zip(
                    sample_info["donor_id"].astype(str),
                    sample_info[condition_col],
                )
            )

    # PN donors
    pn_donors = [
        d for d in unique_donors
        if donor_to_condition.get(d, "").upper() in ("PN", "CONTROL", "NORMAL")
    ]
    logger.info("  Control donors (PN): %d / %d total", len(pn_donors), len(unique_donors))

    if mode == "pn_only":
        if len(pn_donors) == 0:
            logger.warning("  No PN donors found! Whitening will fail.")
        return pn_donors

    elif mode == "pn_plus_early":
        # Add early-window donors if PN count is low
        if len(pn_donors) >= 5:
            logger.info("  Sufficient PN donors (%d), not adding early-window", len(pn_donors))
            return pn_donors

        # Add window 0 donors
        w0_donors = pt_df[pt_df["window"] == 0]["donor_id"].astype(str).unique().tolist()
        early_donors = [d for d in w0_donors if d not in pn_donors]
        combined = pn_donors + early_donors
        logger.info(
            "  pn_plus_early: %d PN + %d early = %d total",
            len(pn_donors), len(early_donors), len(combined),
        )
        return combined

    else:
        raise ValueError(f"Unknown control mode: {mode}")


# =====================================================================
# 5. Gene Reduction for p >> n
# =====================================================================

def _reduce_genes_by_variance(
    residuals: np.ndarray,
    donor_ids: np.ndarray,
    control_donor_ids: Sequence[str],
    gene_indices: np.ndarray,
    target_p: int,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Reduce gene_indices to target_p by variance in control donors.

    Strategy: Keep top-variance genes among control cells.
    If target_p > len(gene_indices), return all.
    """
    if target_p >= len(gene_indices):
        return gene_indices

    # Get control cell mask
    ctrl_mask = np.isin(donor_ids, list(control_donor_ids))
    ctrl_resid = residuals[ctrl_mask][:, gene_indices]

    # Variance per gene
    var = ctrl_resid.var(axis=0)
    # Top target_p by variance
    top_idx = np.argsort(var)[-target_p:]
    top_idx = np.sort(top_idx)  # preserve original order

    return gene_indices[top_idx]


# =====================================================================
# 6. Summary Utilities
# =====================================================================

def summarize_cond(sigma: np.ndarray) -> Dict[str, float]:
    """Compute condition number and eigenvalue summary."""
    eigvals = np.linalg.eigvalsh(sigma)
    return {
        "cond": float(eigvals.max() / max(eigvals.min(), 1e-15)),
        "eigmin": float(eigvals.min()),
        "eigmax": float(eigvals.max()),
        "n_negative": int((eigvals < 0).sum()),
    }


def ensure_finite(arr: np.ndarray, name: str = "array") -> None:
    """Raise if array contains NaN or Inf."""
    if not np.all(np.isfinite(arr)):
        n_nan = int(np.isnan(arr).sum())
        n_inf = int(np.isinf(arr).sum())
        raise RuntimeError(
            f"{name} contains {n_nan} NaN and {n_inf} Inf values"
        )


def save_whitening_meta(meta: Dict[str, Any], output_dir: Path) -> None:
    """Save whitening metadata to JSON."""
    path = output_dir / "whitening_meta.json"
    with open(path, "w") as f:
        json.dump(meta, f, indent=2, default=str)
    logger.info("  Saved whitening_meta.json")


def save_control_ids(control_ids: List[str], output_dir: Path) -> None:
    """Save control donor IDs to JSON."""
    path = output_dir / "whiten_control_ids.json"
    with open(path, "w") as f:
        json.dump({"control_donor_ids": control_ids}, f, indent=2)
    logger.info("  Saved whiten_control_ids.json")
