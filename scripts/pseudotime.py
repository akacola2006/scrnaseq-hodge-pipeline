"""
scRNAseq Hodge Pipeline — Pseudotime Construction
===================================================
PT-B: structure-based pseudotime from SPD covariance (primary).
Window assignment: partition donors into equal-size bins along PT axis.
"""
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.linalg
from scipy.spatial.distance import pdist, squareform

from . import config
from .spd import spd_log, vectorize_upper_tri

logger = logging.getLogger(__name__)


def build_pt_b(
    spd_matrices: Dict[str, Dict[str, np.ndarray]],
    sample_info: pd.DataFrame,
    celltype_set: Optional[List[str]] = None,
    k: int = None,
) -> Tuple[pd.DataFrame, Dict]:
    """Build PT-B from SPD covariance structure via diffusion map.

    Algorithm:
    1. For each donor x celltype, compute log(Sigma_d) and vectorize upper triangle.
    2. Per celltype, z-score standardize across donors.
    3. Concatenate celltype vectors per donor -> donor feature vector.
    4. Only use "complete case" donors (present in all celltypes).
    5. Compute pairwise Euclidean distance matrix.
    6. Diffusion map -> 2nd eigenvector as pseudotime.
    7. Sign convention: control donors have lowest mean PT.
    """
    if celltype_set is None:
        celltype_set = list(config.PT_CELLTYPE_SET)
    if k is None:
        k = config.PCA_K

    vec_len = k * (k + 1) // 2
    logger.info("Building PT-B with celltypes=%s, k=%d", celltype_set, k)

    # Determine complete-case donors
    donor_sets = []
    for ct in celltype_set:
        if ct not in spd_matrices:
            raise ValueError(f"Celltype '{ct}' not found in spd_matrices")
        donor_sets.append(set(spd_matrices[ct].keys()))
    complete_donors = sorted(set.intersection(*donor_sets))
    logger.info("Complete-case donors: %d", len(complete_donors))

    if len(complete_donors) < 5:
        raise ValueError(
            f"Only {len(complete_donors)} complete-case donors; need >= 5"
        )

    # Per-celltype: vectorize log-covariance, then z-score
    ct_features = {}
    for ct in celltype_set:
        vecs = np.empty((len(complete_donors), vec_len), dtype=np.float64)
        for i, donor in enumerate(complete_donors):
            log_sigma = spd_log(spd_matrices[ct][donor])
            vecs[i] = vectorize_upper_tri(log_sigma)

        mu = vecs.mean(axis=0)
        sigma = vecs.std(axis=0, ddof=1)
        sigma[sigma < 1e-12] = 1.0
        vecs_z = (vecs - mu) / sigma
        ct_features[ct] = vecs_z

    # Concatenate into donor feature vectors
    donor_features = np.hstack([ct_features[ct] for ct in celltype_set])
    total_dim = donor_features.shape[1]
    logger.info("Feature matrix: %d donors x %d dims", len(complete_donors), total_dim)

    # Pairwise Euclidean distance
    dist_vec = pdist(donor_features, metric="euclidean")
    dist_matrix = squareform(dist_vec)

    # Diffusion map
    eigvals, eigvecs = _diffusion_map(dist_matrix, n_components=5)

    ev_idx = config.PT_EIGENVECTOR
    pt_raw = eigvecs[:, ev_idx]

    # Sign convention
    cond_col = config.CONDITION_COLUMN
    info_map = sample_info.set_index("donor_id")[cond_col].to_dict() if cond_col in sample_info.columns else {}
    conditions = [info_map.get(d, "UNKNOWN") for d in complete_donors]

    pt_raw = _apply_sign_convention(
        pt_raw, conditions, convention=config.PT_SIGN_CONVENTION
    )

    pt_df = pd.DataFrame({
        "donor_id": complete_donors,
        "condition": conditions,
        "pt_b": pt_raw,
    }).sort_values("pt_b").reset_index(drop=True)

    metadata = {
        "n_complete_donors": len(complete_donors),
        "celltypes_used": celltype_set,
        "feature_dim": total_dim,
        "distance_matrix": dist_matrix,
        "eigenvalues": eigvals,
        "eigenvectors": eigvecs,
        "donor_ids_ordered": complete_donors,
    }

    logger.info(
        "PT-B built: %d donors, range [%.4f, %.4f]",
        len(pt_df), pt_df["pt_b"].min(), pt_df["pt_b"].max(),
    )
    return pt_df, metadata


def _diffusion_map(
    dist_matrix: np.ndarray,
    n_components: int = 5,
    knn_frac: float = 0.3,
) -> Tuple[np.ndarray, np.ndarray]:
    """Diffusion map on a distance matrix."""
    N = dist_matrix.shape[0]
    knn_k = max(3, int(N * knn_frac))

    sigmas = np.empty(N, dtype=np.float64)
    for i in range(N):
        sorted_dists = np.sort(dist_matrix[i])
        knn_dists = sorted_dists[1:knn_k + 1]
        sigmas[i] = np.median(knn_dists)

    sigmas = np.maximum(sigmas, 1e-10)

    bandwidth_matrix = np.outer(sigmas, sigmas)
    affinity = np.exp(-(dist_matrix ** 2) / bandwidth_matrix)

    row_sums = affinity.sum(axis=1)
    row_sums = np.maximum(row_sums, 1e-15)

    D_inv_sqrt = np.diag(1.0 / np.sqrt(row_sums))
    P_sym = D_inv_sqrt @ affinity @ D_inv_sqrt
    P_sym = (P_sym + P_sym.T) / 2.0

    n_eig = min(n_components, N - 1)
    eigenvalues, eigenvectors = scipy.linalg.eigh(
        P_sym,
        subset_by_index=[N - n_eig, N - 1],
    )

    sort_idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[sort_idx]
    eigenvectors = eigenvectors[:, sort_idx]

    eigenvectors = D_inv_sqrt @ eigenvectors

    logger.info("Diffusion map: top eigenvalues = %s",
                np.array2string(eigenvalues[:n_components], precision=4))

    return eigenvalues, eigenvectors


def _apply_sign_convention(
    pt_values: np.ndarray,
    conditions: List[str],
    convention: str = "control_low",
) -> np.ndarray:
    """Apply sign convention to pseudotime values.

    Supports:
    - "control_low": Control/reference donors have lowest mean PT.
    - "PN_low": Legacy alias for control_low (uses "PN" as control label).
    """
    pt = pt_values.copy()

    if convention in ("control_low", "PN_low"):
        control_label = config.CONTROL_LABEL
        # For legacy "PN_low" convention
        if convention == "PN_low":
            control_label = "PN"

        control_mask = np.array([c == control_label for c in conditions])
        if control_mask.any():
            control_mean = pt[control_mask].mean()
            overall_median = np.median(pt)
            if control_mean > overall_median:
                pt = -pt
                logger.info("PT sign flipped: control mean (%.4f) > median (%.4f)",
                           control_mean, overall_median)
            else:
                logger.info("PT sign OK: control mean (%.4f) <= median (%.4f)",
                           control_mean, overall_median)
        else:
            logger.warning("No control donors found (label='%s'); skipping sign convention",
                          control_label)
    else:
        raise ValueError(f"Unknown sign convention: {convention}")

    return pt


def assign_windows(
    pt_df: pd.DataFrame,
    n_windows: int = None,
    donor_subset: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Assign donors to equal-size windows based on PT ordering."""
    if n_windows is None:
        n_windows = config.N_WINDOWS

    df = pt_df.copy()

    if "pt_b" in df.columns:
        pt_col = "pt_b"
    elif "pt_a" in df.columns:
        pt_col = "pt_a"
    else:
        raise ValueError("pt_df must contain 'pt_b' or 'pt_a' column")

    if donor_subset is not None:
        subset_set = set(str(d) for d in donor_subset)
        df = df[df["donor_id"].isin(subset_set)].copy()

    if len(df) == 0:
        raise ValueError("No donors remain after subsetting")

    df = df.sort_values(pt_col).reset_index(drop=True)

    n_donors = len(df)
    if n_donors < n_windows:
        logger.warning("Fewer donors (%d) than windows (%d); adjusting", n_donors, n_windows)
        n_windows = n_donors

    window_labels = np.empty(n_donors, dtype=int)
    splits = np.array_split(np.arange(n_donors), n_windows)
    for w_idx, indices in enumerate(splits):
        window_labels[indices] = w_idx

    df["window"] = window_labels

    window_sizes = df.groupby("window").size()
    logger.info("Assigned %d donors to %d windows (sizes: %s)",
                n_donors, n_windows, window_sizes.to_dict())

    return df
