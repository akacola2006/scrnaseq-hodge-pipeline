"""
scRNAseq Hodge Pipeline — 2-Axis Model (TRS x MSS)
=====================================================
Computes two orthogonal scores per celltype x window:

  TRS (Translation-Resource Score):
    Measures coordination among upstream (High-phi) genes.
    TRS_spd = mean upper-triangle of donor-level log-correlation among High genes.

  MSS (Morphogenesis-Structure Score):
    Measures structural gene program activity.
    MSS_pc1 = PC1 of a user-defined structural gene set (default: all Medium-phi genes),
    computed from residual expression.

The 2-axis trajectory (TRS, MSS) across windows reveals whether
disease progression follows a supply-demand pattern:
  - TRS leads (upstream genes coordinate first)
  - MSS follows (structural changes appear later)
"""
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import spearmanr
from sklearn.covariance import LedoitWolf
from sklearn.decomposition import PCA

from . import config
from .spd import spd_log

logger = logging.getLogger(__name__)


def compute_trs(
    residuals: np.ndarray,
    donor_ids: np.ndarray,
    gene_indices: np.ndarray,
    donor_window: Dict[str, int],
    n_windows: int = None,
) -> Dict[int, Dict[str, float]]:
    """Compute TRS (Translation-Resource Score) per window.

    TRS = mean upper-triangle of donor log-correlation among target genes.

    Parameters
    ----------
    residuals : np.ndarray, shape (n_cells, n_genes)
    donor_ids : np.ndarray, shape (n_cells,)
    gene_indices : np.ndarray
        Indices of High-phi genes in residual matrix.
    donor_window : dict
        {donor_id: window_index}.
    n_windows : int

    Returns
    -------
    dict of {window: {"mean": float, "std": float, "n_donors": int}}
    """
    if n_windows is None:
        n_windows = config.N_WINDOWS

    n_genes_sel = len(gene_indices)
    unique_donors = sorted(set(donor_ids))

    trs_per_window = {w: [] for w in range(n_windows)}

    for d_id in unique_donors:
        d_str = str(d_id)
        if d_str not in donor_window:
            continue
        w = donor_window[d_str]

        mask = donor_ids == d_id
        if mask.sum() < max(30, n_genes_sel // 5):
            continue

        X_d = residuals[mask][:, gene_indices]

        try:
            lw = LedoitWolf()
            lw.fit(X_d)
            cov = lw.covariance_

            d_std = np.sqrt(np.diag(cov))
            d_std[d_std < 1e-15] = 1.0
            corr = cov / np.outer(d_std, d_std)
            np.fill_diagonal(corr, 1.0)
            corr = np.clip(corr, -1.0, 1.0)

            # Ensure SPD
            eigvals, eigvecs = np.linalg.eigh(corr)
            eigvals = np.maximum(eigvals, 1e-10)
            corr = (eigvecs * eigvals) @ eigvecs.T
            corr = (corr + corr.T) / 2.0

            log_corr = spd_log(corr)
            triu_idx = np.triu_indices(n_genes_sel, k=1)
            trs_val = float(np.mean(log_corr[triu_idx]))
            trs_per_window[w].append(trs_val)
        except Exception:
            continue

    result = {}
    for w in range(n_windows):
        vals = trs_per_window[w]
        if vals:
            result[w] = {"mean": float(np.mean(vals)), "std": float(np.std(vals)),
                         "n_donors": len(vals)}
        else:
            result[w] = {"mean": float("nan"), "std": 0.0, "n_donors": 0}

    return result


def compute_mss(
    residuals: np.ndarray,
    donor_ids: np.ndarray,
    gene_indices: np.ndarray,
    donor_window: Dict[str, int],
    n_windows: int = None,
) -> Dict[int, Dict[str, float]]:
    """Compute MSS (Morphogenesis-Structure Score) per window.

    MSS = PC1 projection of structural gene set expression, averaged per window.

    Parameters
    ----------
    residuals : np.ndarray
    donor_ids : np.ndarray
    gene_indices : np.ndarray
        Indices of structural/Medium-phi genes.
    donor_window : dict
    n_windows : int

    Returns
    -------
    dict of {window: {"mean": float, "std": float, "n_donors": int}}
    """
    if n_windows is None:
        n_windows = config.N_WINDOWS

    unique_donors = sorted(set(donor_ids))

    # Compute PC1 across all cells for structural genes
    X_struct = residuals[:, gene_indices]
    pca = PCA(n_components=1)
    pc1_all = pca.fit_transform(X_struct).ravel()
    explained_var = float(pca.explained_variance_ratio_[0])
    logger.info("  MSS PC1 explained variance: %.1f%%", explained_var * 100)

    # Average PC1 per donor, then per window
    mss_per_window = {w: [] for w in range(n_windows)}

    for d_id in unique_donors:
        d_str = str(d_id)
        if d_str not in donor_window:
            continue
        w = donor_window[d_str]

        mask = donor_ids == d_id
        if mask.sum() < 10:
            continue

        donor_mss = float(np.mean(pc1_all[mask]))
        mss_per_window[w].append(donor_mss)

    result = {}
    for w in range(n_windows):
        vals = mss_per_window[w]
        if vals:
            result[w] = {"mean": float(np.mean(vals)), "std": float(np.std(vals)),
                         "n_donors": len(vals)}
        else:
            result[w] = {"mean": float("nan"), "std": 0.0, "n_donors": 0}

    return result


def run_two_axis_model(
    residuals: np.ndarray,
    donor_ids: np.ndarray,
    all_gene_names: List[str],
    hodge_gene_names: List[str],
    phi_scores: np.ndarray,
    classification: List[str],
    donor_window: Dict[str, int],
    n_windows: int = None,
    output_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run 2-axis model: TRS (High genes) x MSS (Medium genes).

    Parameters
    ----------
    residuals : np.ndarray, shape (n_cells, n_genes_kept)
        Full residual matrix for this cell type.
    donor_ids : np.ndarray, shape (n_cells,)
    all_gene_names : list of str
        Gene names corresponding to residual columns (length = n_genes_kept).
    hodge_gene_names : list of str
        Gene names from Hodge analysis (subset of all_gene_names).
    phi_scores : np.ndarray, shape (n_hodge_genes,)
    classification : list of str ("High", "Medium", "Low")
    donor_window : dict {donor_id_str: window_int}
    n_windows : int
    output_dir : Path, optional

    Returns
    -------
    dict with TRS and MSS per window, lead/lag analysis.
    """
    if n_windows is None:
        n_windows = config.N_WINDOWS

    # Build residual column index lookup
    resid_name_to_idx = {g: i for i, g in enumerate(all_gene_names)}

    # Separate High and Medium genes, map to residual column indices
    high_genes = [g for g, c in zip(hodge_gene_names, classification) if c == "High"]
    medium_genes = [g for g, c in zip(hodge_gene_names, classification) if c == "Medium"]

    logger.info("2-Axis Model: %d High genes (TRS), %d Medium genes (MSS)",
                len(high_genes), len(medium_genes))

    if len(high_genes) < 5:
        logger.warning("Too few High genes for TRS")
        return {"status": "SKIPPED", "reason": "too few High genes"}
    if len(medium_genes) < 5:
        logger.warning("Too few Medium genes for MSS")
        return {"status": "SKIPPED", "reason": "too few Medium genes"}

    # Map to residual matrix column indices
    high_idx = np.array([resid_name_to_idx[g] for g in high_genes if g in resid_name_to_idx])
    med_idx = np.array([resid_name_to_idx[g] for g in medium_genes if g in resid_name_to_idx])

    if len(high_idx) < 5 or len(med_idx) < 5:
        return {"status": "SKIPPED", "reason": "too few genes resolved in residual matrix"}

    # Compute TRS and MSS
    trs = compute_trs(residuals, donor_ids, high_idx, donor_window, n_windows)
    mss = compute_mss(residuals, donor_ids, med_idx, donor_window, n_windows)

    # Lead/lag analysis: does TRS peak before MSS?
    trs_means = [trs[w]["mean"] for w in range(n_windows)]
    mss_means = [mss[w]["mean"] for w in range(n_windows)]

    valid_trs = [(w, v) for w, v in enumerate(trs_means) if not np.isnan(v)]
    valid_mss = [(w, v) for w, v in enumerate(mss_means) if not np.isnan(v)]

    trs_peak_window = max(valid_trs, key=lambda x: abs(x[1]))[0] if valid_trs else None
    mss_peak_window = max(valid_mss, key=lambda x: abs(x[1]))[0] if valid_mss else None

    lead_lag = None
    if trs_peak_window is not None and mss_peak_window is not None:
        if trs_peak_window < mss_peak_window:
            lead_lag = "TRS_leads"
        elif trs_peak_window > mss_peak_window:
            lead_lag = "MSS_leads"
        else:
            lead_lag = "simultaneous"

    result = {
        "status": "OK",
        "n_high_genes": len(high_genes),
        "n_medium_genes": len(medium_genes),
        "trs_per_window": trs,
        "mss_per_window": mss,
        "trs_peak_window": trs_peak_window,
        "mss_peak_window": mss_peak_window,
        "lead_lag": lead_lag,
    }

    logger.info("  TRS peak: window %s, MSS peak: window %s -> %s",
                trs_peak_window, mss_peak_window, lead_lag)

    if output_dir is not None:
        import pandas as pd
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Trajectory table
        rows = []
        for w in range(n_windows):
            rows.append({
                "window": w,
                "trs_mean": trs[w]["mean"],
                "trs_std": trs[w]["std"],
                "trs_n": trs[w]["n_donors"],
                "mss_mean": mss[w]["mean"],
                "mss_std": mss[w]["std"],
                "mss_n": mss[w]["n_donors"],
            })
        pd.DataFrame(rows).to_csv(output_dir / "two_axis_trajectory.csv", index=False)

        summary = {k: v for k, v in result.items()
                   if k not in ("trs_per_window", "mss_per_window")}
        with open(output_dir / "two_axis_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)

        logger.info("2-Axis results saved to %s", output_dir)

    return result
