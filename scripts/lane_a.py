"""
scRNAseq Hodge Pipeline — Lane A: Cell-type Upstream Analysis
==============================================================
Steps:
1. Group donors into K windows by pseudotime.
2. Per window, per celltype: compute Log-Euclidean mean of SPD matrices.
3. Decompose distances into d_cov, d_corr, d_var per celltype.
4. Build influence flow -> Hodge decomposition -> celltype phi (upstream score).
5. Permutation test for gradient dominance.
"""
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import config
from .spd import log_euclidean_mean, decompose_distance
from .hodge import (
    build_incidence_matrices,
    hodge_decomposition,
    compute_spd_velocity,
    build_flow_pairwise,
    permutation_test_gradient,
)

logger = logging.getLogger(__name__)


def _compute_window_representatives(
    spd_matrices: Dict[str, Dict[str, np.ndarray]],
    pt_df: pd.DataFrame,
    celltypes: List[str],
    n_windows: int,
    min_donors_per_window: int,
) -> Dict[str, List[Optional[np.ndarray]]]:
    """Compute per-window, per-celltype SPD representative matrices."""
    if "window" not in pt_df.columns:
        raise ValueError("pt_df must contain 'window' column")

    donor_to_window = dict(zip(pt_df["donor_id"].astype(str), pt_df["window"].astype(int)))

    representatives = {}
    for ct in celltypes:
        ct_reps = []
        if ct not in spd_matrices:
            logger.warning("Celltype '%s' not in spd_matrices; filling with None", ct)
            representatives[ct] = [None] * n_windows
            continue

        ct_donors = spd_matrices[ct]
        for w in range(n_windows):
            window_matrices = [
                sigma for d_id, sigma in ct_donors.items()
                if donor_to_window.get(str(d_id)) == w
            ]

            if len(window_matrices) < min_donors_per_window:
                ct_reps.append(None)
            else:
                ct_reps.append(log_euclidean_mean(window_matrices))

        representatives[ct] = ct_reps

    return representatives


def run_lane_a(
    spd_matrices: Dict[str, Dict[str, np.ndarray]],
    pt_df: pd.DataFrame,
    celltypes: List[str],
    n_windows: int = None,
    min_donors_per_window: int = None,
    output_dir: Optional[Path] = None,
    run_permutation: bool = True,
) -> Dict[str, Any]:
    """Run full Lane A analysis.

    Returns dict with: phi_scores, gradient_fraction, valid_celltypes,
    hodge_result, distance_decomp, permutation_result.
    """
    if n_windows is None:
        n_windows = config.N_WINDOWS
    if min_donors_per_window is None:
        min_donors_per_window = config.MIN_DONORS_PER_WINDOW

    logger.info("Lane A: %d celltypes, %d windows", len(celltypes), n_windows)

    # Step 1: Window representatives
    representatives = _compute_window_representatives(
        spd_matrices, pt_df, celltypes, n_windows, min_donors_per_window,
    )

    # Step 2: Distance decomposition
    distance_decomp = {}
    for ct in celltypes:
        reps = representatives[ct]
        d_corrs = []
        pairs = []
        for w in range(n_windows - 1):
            if reps[w] is not None and reps[w + 1] is not None:
                dd = decompose_distance(reps[w], reps[w + 1])
                d_corrs.append(dd["d_corr"])
                pairs.append((w, w + 1))
        distance_decomp[ct] = {"d_corr": d_corrs, "window_pairs": pairs}

    # Determine valid celltypes
    valid_celltypes = sorted(
        ct for ct in celltypes if distance_decomp[ct]["d_corr"]
    )
    if len(valid_celltypes) < 2:
        logger.error("Fewer than 2 valid celltypes; cannot run Hodge")
        return {"phi_scores": {}, "gradient_fraction": 0.0,
                "valid_celltypes": valid_celltypes, "status": "FAILED"}

    # Step 3: Pairwise flow + Hodge
    d_corr_for_flow = {ct: distance_decomp[ct] for ct in valid_celltypes}
    flow = build_flow_pairwise(d_corr_for_flow, valid_celltypes)
    n_nodes = len(valid_celltypes)
    B0, B1 = build_incidence_matrices(n_nodes)
    hodge_result = hodge_decomposition(flow, B0, B1)

    phi = hodge_result["phi"]
    phi_scores = {ct: float(phi[i]) for i, ct in enumerate(valid_celltypes)}

    # Rank by phi
    ranked = sorted(phi_scores.items(), key=lambda x: -x[1])
    logger.info("Cell-type ranking by phi: %s",
                [(ct, f"{v:.4f}") for ct, v in ranked])

    # Step 4: Permutation test
    permutation_result = None
    if run_permutation:
        permutation_result = permutation_test_gradient(flow, B0, B1)

    result = {
        "phi_scores": phi_scores,
        "gradient_fraction": hodge_result["gradient_fraction"],
        "curl_fraction": hodge_result["curl_fraction"],
        "harmonic_fraction": hodge_result["harmonic_fraction"],
        "valid_celltypes": valid_celltypes,
        "hodge_result": hodge_result,
        "distance_decomp": distance_decomp,
        "permutation_result": permutation_result,
        "ranked_celltypes": [ct for ct, _ in ranked],
        "status": "OK",
    }

    # Save results
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        summary = {
            "phi_scores": phi_scores,
            "gradient_fraction": hodge_result["gradient_fraction"],
            "curl_fraction": hodge_result["curl_fraction"],
            "harmonic_fraction": hodge_result["harmonic_fraction"],
            "ranked_celltypes": [ct for ct, _ in ranked],
            "valid_celltypes": valid_celltypes,
        }
        if permutation_result:
            summary["permutation_p_value"] = permutation_result["p_value"]

        with open(output_dir / "lane_a_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)

        logger.info("Lane A results saved to %s", output_dir)

    return result
