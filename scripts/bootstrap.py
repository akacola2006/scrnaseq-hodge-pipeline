"""
scRNAseq Hodge Pipeline — Bootstrap Confidence Estimation
==========================================================
Donor-level cluster bootstrap: resample donors with replacement,
re-run Lane A, and collect per-iteration statistics over B iterations.
"""
import json
import logging
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import config
from .seed_utils import get_bootstrap_seed
from .log_utils import FailureLogger
from .spd import log_euclidean_mean, decompose_distance
from .hodge import (
    build_incidence_matrices,
    hodge_decomposition,
    build_flow_pairwise,
)

logger = logging.getLogger(__name__)


def _resample_donors(spd_matrices, pt_df, rng):
    """Resample donors with replacement (cluster bootstrap)."""
    donors = pt_df["donor_id"].unique()
    n_donors = len(donors)
    boot_indices = rng.choice(n_donors, size=n_donors, replace=True)

    rows = []
    for boot_idx, orig_idx in enumerate(boot_indices):
        orig_donor = donors[orig_idx]
        orig_row = pt_df[pt_df["donor_id"] == orig_donor].iloc[0].copy()
        orig_row["donor_id"] = f"{orig_donor}__b{boot_idx}"
        rows.append(orig_row)
    resampled_pt = pd.DataFrame(rows).reset_index(drop=True)

    resampled_spd = {}
    for ct, donor_dict in spd_matrices.items():
        resampled_spd[ct] = {}
        for boot_idx, orig_idx in enumerate(boot_indices):
            orig_donor = str(donors[orig_idx])
            if orig_donor in donor_dict:
                new_id = f"{orig_donor}__b{boot_idx}"
                resampled_spd[ct][new_id] = donor_dict[orig_donor]

    return resampled_spd, resampled_pt


def _run_lane_a_lite(spd_matrices, pt_df, celltypes, n_windows, min_donors_per_window):
    """Lightweight Lane A for bootstrap: no permutation test, no file I/O."""
    donor_to_window = dict(
        zip(pt_df["donor_id"].astype(str), pt_df["window"].astype(int))
    )

    representatives = {}
    for ct in celltypes:
        if ct not in spd_matrices:
            representatives[ct] = [None] * n_windows
            continue
        ct_donors = spd_matrices[ct]
        ct_reps = []
        for w in range(n_windows):
            window_matrices = [
                sigma for d_id, sigma in ct_donors.items()
                if donor_to_window.get(str(d_id)) == w
            ]
            if len(window_matrices) >= min_donors_per_window:
                ct_reps.append(log_euclidean_mean(window_matrices))
            else:
                ct_reps.append(None)
        representatives[ct] = ct_reps

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

    valid_celltypes = sorted(ct for ct in celltypes if distance_decomp[ct]["d_corr"])
    if len(valid_celltypes) < 2:
        return {"phi_scores": {}, "gradient_fraction": 0.0, "top_k": [], "status": "FAILED"}

    d_corr_for_flow = {ct: distance_decomp[ct] for ct in valid_celltypes}
    flow = build_flow_pairwise(d_corr_for_flow, valid_celltypes)
    n_nodes = len(valid_celltypes)
    B0, B1 = build_incidence_matrices(n_nodes)
    hodge_result = hodge_decomposition(flow, B0, B1)

    phi = hodge_result["phi"]
    phi_scores = {ct: float(phi[i]) for i, ct in enumerate(valid_celltypes)}

    ranked = sorted(phi_scores.items(), key=lambda x: -x[1])
    top_k = [ct for ct, _ in ranked[:config.BOOTSTRAP_TOP_K]]

    return {
        "phi_scores": phi_scores,
        "gradient_fraction": hodge_result["gradient_fraction"],
        "top_k": top_k,
        "status": "OK",
    }


def _bootstrap_single_iteration(spd_matrices, pt_df, celltypes, seed, n_windows,
                                 min_donors_per_window=None):
    if min_donors_per_window is None:
        min_donors_per_window = config.MIN_DONORS_PER_WINDOW
    try:
        rng = np.random.RandomState(seed)
        resampled_spd, resampled_pt = _resample_donors(spd_matrices, pt_df, rng)
        result = _run_lane_a_lite(
            resampled_spd, resampled_pt, celltypes, n_windows, min_donors_per_window,
        )
        result["seed"] = seed
        result["error"] = None
        return result
    except Exception as exc:
        return {
            "seed": seed, "phi_scores": {}, "gradient_fraction": 0.0,
            "top_k": [], "status": "FAILED",
            "error": f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
        }


def bootstrap_lane_a(
    spd_matrices=None, pt_df=None, celltypes=None,
    n_bootstrap=None, base_seed=None, n_workers=None,
    output_dir=None, failure_logger=None,
):
    """Run B bootstrap iterations of Lane A."""
    if n_bootstrap is None:
        n_bootstrap = config.BOOTSTRAP_N
    if base_seed is None:
        base_seed = config.BOOTSTRAP_BASE_SEED
    if n_workers is None:
        n_workers = config.BOOTSTRAP_N_WORKERS

    n_windows = pt_df["window"].nunique() if "window" in pt_df.columns else config.N_WINDOWS

    logger.info("Bootstrap: B=%d, %d celltypes, %d workers", n_bootstrap, len(celltypes), n_workers)

    iteration_seeds = [get_bootstrap_seed(base_seed, i) for i in range(n_bootstrap)]

    results = []
    if n_workers <= 1:
        for i, seed in enumerate(iteration_seeds):
            if (i + 1) % 10 == 0 or i == 0:
                logger.info("Bootstrap %d/%d", i + 1, n_bootstrap)
            res = _bootstrap_single_iteration(spd_matrices, pt_df, celltypes, seed, n_windows)
            results.append(res)
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(
                    _bootstrap_single_iteration,
                    spd_matrices, pt_df, celltypes, seed, n_windows,
                ): (i, seed) for i, seed in enumerate(iteration_seeds)
            }
            for future in as_completed(futures):
                i, seed = futures[future]
                try:
                    res = future.result()
                except Exception as exc:
                    res = {"seed": seed, "phi_scores": {}, "gradient_fraction": 0.0,
                           "top_k": [], "status": "FAILED", "error": str(exc)}
                results.append(res)

    successful = [r for r in results if r["status"] == "OK"]
    failed = [r for r in results if r["status"] == "FAILED"]
    n_success = len(successful)
    n_failed = len(failed)
    logger.info("Bootstrap: %d/%d successful", n_success, n_bootstrap)

    if n_success == 0:
        return {
            "phi_matrix": np.array([]), "phi_celltypes": celltypes,
            "gradient_fractions": np.array([]),
            "top3_reproducibility": {ct: 0.0 for ct in celltypes},
            "phi_ci": {ct: (np.nan, np.nan) for ct in celltypes},
            "gradient_fraction_ci": (np.nan, np.nan),
            "n_success": 0, "n_failed": n_failed,
            "failed_seeds": [r["seed"] for r in failed],
        }

    all_cts = sorted(set().union(*(r["phi_scores"].keys() for r in successful)))
    ct_to_idx = {ct: i for i, ct in enumerate(all_cts)}

    phi_matrix = np.full((n_success, len(all_cts)), np.nan, dtype=np.float64)
    gradient_fractions = np.empty(n_success, dtype=np.float64)

    for b, r in enumerate(successful):
        for ct, val in r["phi_scores"].items():
            phi_matrix[b, ct_to_idx[ct]] = val
        gradient_fractions[b] = r["gradient_fraction"]

    top_k_counts = {ct: 0 for ct in all_cts}
    for r in successful:
        for ct in r["top_k"]:
            if ct in top_k_counts:
                top_k_counts[ct] += 1
    top3_reproducibility = {ct: count / n_success for ct, count in top_k_counts.items()}

    phi_ci = {}
    for ct in all_cts:
        col = phi_matrix[:, ct_to_idx[ct]]
        valid = col[~np.isnan(col)]
        if len(valid) >= 2:
            phi_ci[ct] = (float(np.percentile(valid, 2.5)), float(np.percentile(valid, 97.5)))
        else:
            phi_ci[ct] = (np.nan, np.nan)

    gradient_fraction_ci = (
        float(np.percentile(gradient_fractions, 2.5)),
        float(np.percentile(gradient_fractions, 97.5)),
    )

    stability_classification = {}
    for ct, repro in top3_reproducibility.items():
        if repro >= config.BOOTSTRAP_STABLE_THRESHOLD:
            stability_classification[ct] = "stable"
        elif repro >= config.BOOTSTRAP_SUGGESTIVE_THRESHOLD:
            stability_classification[ct] = "suggestive"
        else:
            stability_classification[ct] = "unstable"

    result = {
        "phi_matrix": phi_matrix,
        "phi_celltypes": all_cts,
        "gradient_fractions": gradient_fractions,
        "top3_reproducibility": top3_reproducibility,
        "phi_ci": phi_ci,
        "gradient_fraction_ci": gradient_fraction_ci,
        "stability_classification": stability_classification,
        "n_success": n_success, "n_failed": n_failed,
        "failed_seeds": [r["seed"] for r in failed],
    }

    if output_dir is not None:
        _save_bootstrap_results(result, output_dir)

    return result


def _save_bootstrap_results(result, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    np.save(output_dir / "bootstrap_phi_matrix.npy", result["phi_matrix"])
    np.save(output_dir / "bootstrap_gradient_fractions.npy", result["gradient_fractions"])

    summary = {
        "phi_celltypes": result["phi_celltypes"],
        "top3_reproducibility": result["top3_reproducibility"],
        "phi_ci": {ct: list(ci) for ct, ci in result["phi_ci"].items()},
        "gradient_fraction_ci": list(result["gradient_fraction_ci"]),
        "stability_classification": result["stability_classification"],
        "n_success": result["n_success"],
        "n_failed": result["n_failed"],
    }
    with open(output_dir / "bootstrap_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    logger.info("Bootstrap results saved to %s", output_dir)
