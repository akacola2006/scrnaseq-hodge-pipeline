"""
IDS evaluation on Norman et al. 2019 Perturb-seq (paper Section 2.4).

Reproduces the Top-1 source identification accuracy reported in
Section 2.4 ("16.4% without whitening, 66.4% with whitening" on Norman).
This is the SAME IDS R2b pipeline as run_ids_replogle.py, but pointed at
the Norman test split (paper test cells, K_test = 235 perturbation
categories) and using dynamically-constructed candidate genes from the
perturbation labels.

Pipeline (faithfully ported from rebuild_faithful.py
ids_hodge_decomposition_v2, which produced the paper Table 1 numbers):

  1. Load norman_test.h5ad (34,624 cells, 235 categories).
  2. Build candidate gene set DYNAMICALLY from perturbation labels:
       for each category != 'control', split on '_' and keep tokens
       that match adata.var_names.
     (This matches NormanCellNaviDataset; the static candidate_genes.txt
     used by older scripts gives a different, smaller candidate pool.)
  3. Normalise per-cell: log1p(TP10K).
  4. Global whitening on all control cells, candidate columns only:
       Σ_ctrl = LedoitWolf(ctrl_expr)
       W = V D^{-1/2} V^T  (symmetric ZCA)
       X_w = X @ W
  5. Per perturbation group with n_pert >= 3 and >= 1 candidate target:
       Bootstrap n=100 (seed=42), control subsample size
         min(max(n_pert, 100), n_ctrl).
       Per bootstrap:
         lm_p = matrix_log(LedoitWolf(pert_expr_w[idx_p]))
         lm_c = matrix_log(LedoitWolf(ctrl_expr_w[idx_c]))
         delta = lm_p - lm_c
         phi   = analytical K_N Hodge potential of edge-weight flow
                 f(i,j) = |delta_ij| * sign(d_i - d_j),  d_i = ||delta_{i,:}||_2
                 div_i = -sum_{j != i} |delta_ij| * sign(d_i - d_j)
                 phi   = (div - mean(div)) / N
       Mean phi over bootstraps -> rank ascending (lowest phi = upstream).
  6. Report Top-1 / Top-5 / Top-10 hit rates and median rank.

Reference numbers (paper Table 1, Norman, run with the same code path):
    IDS R2b (whitened + Hodge) : ~66% Top-1, median rank 1.0
    d_corr  (whitened, no Hodge): ~57% Top-1, median rank 1.0
    Random                     :  ~1% Top-1

Data preparation
----------------
The Norman test split (norman_test.h5ad) is the paper test cells
extracted from Norman et al. 2019 (GSE133344) via scPerturb (see
benchmark/norman_perturb_seq/README.md for download instructions).
Place it under data/external/norman/processed/ or override the location
with the environment variable ALS_NORMAN_DATA_DIR.

Usage
-----
    export ALS_NORMAN_DATA_DIR=/path/to/norman_data
    python benchmark/norman_perturb_seq/run_ids_norman.py
"""

import gc
import json
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.covariance import LedoitWolf

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger(__name__)

# norman_test.h5ad lives directly inside DATA_DIR by default.
DATA_DIR = Path(os.environ.get("ALS_NORMAN_DATA_DIR", "data/external/norman"))
DATA_PATH = Path(os.environ.get(
    "ALS_NORMAN_TEST_H5AD",
    str(DATA_DIR / "processed" / "norman_test.h5ad"),
))
OUTPUT_DIR = Path(os.environ.get("ALS_NORMAN_OUTPUT_DIR", "data/external/norman_results"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Core R2b components (identical to run_ids_replogle.py)
# ---------------------------------------------------------------------------

def compute_log_matrix(X_norm):
    """LedoitWolf covariance -> correlation -> SPD guarantee -> matrix log."""
    lw = LedoitWolf()
    cov = lw.fit(X_norm).covariance_
    d = np.sqrt(np.maximum(np.diag(cov), 1e-30))
    mat = cov / np.outer(d, d)
    np.fill_diagonal(mat, 1.0)
    mat = (mat + mat.T) / 2
    eigvals, eigvecs = np.linalg.eigh(mat)
    eigvals = np.maximum(eigvals, 1e-15)
    log_mat = (eigvecs * np.log(eigvals)) @ eigvecs.T
    return (log_mat + log_mat.T) / 2


def hodge_edge_weight_phi(delta):
    """Analytical K_N Hodge potential for edge-weight flow.

    f(i, j) = |delta_ij| * sign(d_i - d_j),  d_i = ||delta_{i,:}||_2
    div_i   = -sum_{j != i} |delta_ij| * sign(d_i - d_j)
    phi     = (div - mean(div)) / N

    Sign convention
    ---------------
    This function returns the *raw* Hodge potential ``phi_raw``. With the
    flow definition above, divergence is most negative at sources, so
    ``phi_raw`` is **ascending = upstream** (lowest phi_raw = predicted
    source) -- this is the convention used internally throughout the
    pipeline (BV.61, BV.65, etc.).

    The manuscript uses a *source-oriented* presentation where larger
    values denote more upstream / source-like behaviour. To get that
    convention, negate the returned value::

        phi_raw       = hodge_edge_weight_phi(delta)   # ascending = upstream
        source_score  = -phi_raw                       # descending = upstream

    The ranking in this script is by ascending ``phi_raw`` (equivalently,
    descending ``source_score``); they differ only by sign.
    """
    d = np.linalg.norm(delta, axis=1)
    abs_delta = np.abs(delta).copy()
    np.fill_diagonal(abs_delta, 0)
    sign_d = np.sign(d[:, None] - d[None, :])
    div = -(abs_delta * sign_d).sum(axis=1)
    N = delta.shape[0]
    return (div - div.mean()) / N


def hodge_source_score(delta):
    """Source-oriented potential (manuscript convention).

    Returns ``-hodge_edge_weight_phi(delta)`` so that larger values
    indicate more upstream / source-like genes, matching the sign
    convention used in the manuscript figures and tables
    (e.g. "Oligo phi = 0.900" with positive = upstream). Internally
    identical to ranking by ascending ``hodge_edge_weight_phi``.
    """
    return -hodge_edge_weight_phi(delta)


def global_whitening(ctrl_expr):
    """Symmetric ZCA whitening from LedoitWolf control covariance."""
    lw = LedoitWolf()
    cov = lw.fit(ctrl_expr).covariance_
    eigvals, eigvecs = np.linalg.eigh(cov)
    eigvals = np.maximum(eigvals, 1e-10)
    cond = eigvals[-1] / eigvals[0]
    W = (eigvecs / np.sqrt(eigvals)) @ eigvecs.T
    logger.info(f"  Whitening: cond(Sigma) = {cond:.0f}")
    return W


# ---------------------------------------------------------------------------
# Norman test split loader + dynamic candidate construction
# ---------------------------------------------------------------------------

def load_norman_test():
    if not DATA_PATH.exists():
        raise FileNotFoundError(
            f"{DATA_PATH} not found.\n"
            "Download Norman et al. 2019 via scPerturb and extract the paper "
            "test split. See benchmark/norman_perturb_seq/README.md."
        )
    logger.info(f"Loading {DATA_PATH} ...")
    adata = sc.read_h5ad(DATA_PATH)
    logger.info(f"Loaded: {adata.shape}")
    return adata


# ---------------------------------------------------------------------------
# Main evaluation
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("IDS R2b on Norman Perturb-seq (paper Section 2.4)")
    print("  Pipeline: log1p TP10K -> whitening -> SPD log -> K_N Hodge phi")
    print("  Candidates: dynamic from perturbation categories")
    print("=" * 70)

    # 1. Load Norman test split
    adata = load_norman_test()

    pert_col = "perturbation"
    if pert_col not in adata.obs.columns:
        # Fall back to other common names.
        for col in ("gene", "target_gene", "guide_identity"):
            if col in adata.obs.columns:
                pert_col = col
                break
        else:
            raise ValueError("No perturbation column found in adata.obs")
    logger.info(f"Using perturbation column: {pert_col}")

    # 2. Perturbation categories (matches NormanCellNaviDataset).
    if hasattr(adata.obs[pert_col], "cat"):
        categories = list(adata.obs[pert_col].cat.categories)
    else:
        categories = sorted(adata.obs[pert_col].unique().tolist())
    logger.info(f"Perturbation categories: {len(categories)}")

    # 3. Dynamic candidate gene set: tokens of each category that are in var_names.
    gene_names = adata.var_names.tolist()
    gene_name_set = set(gene_names)
    candidate_genes = sorted({
        g for cat in categories if cat != "control"
        for g in cat.split("_") if g in gene_name_set
    })
    candidate_indices = np.array([gene_names.index(g) for g in candidate_genes])
    n_cand = len(candidate_genes)
    logger.info(f"Dynamic candidate genes: {n_cand}")

    # 4. Control / perturbed cell masks.
    adata.obs["is_control"] = adata.obs[pert_col] == "control"
    n_control = int(adata.obs["is_control"].sum())
    n_perturbed = int((~adata.obs["is_control"]).sum())
    logger.info(f"Control cells: {n_control}, Perturbed cells: {n_perturbed}")

    # 5. Build group_indices (skip control + n_pert < 3 + no candidate targets).
    groups = {}
    pert_values = adata.obs[pert_col].values
    for cat in categories:
        if cat == "control":
            continue
        mask = pert_values == cat
        n_cells = int(mask.sum())
        if n_cells < 3:
            continue
        targets = [g for g in cat.split("_") if g in gene_name_set]
        if not targets:
            continue
        groups[cat] = {"targets": targets, "n_cells": n_cells, "mask": mask}
    logger.info(f"Evaluable groups (>=3 cells, candidate targets): {len(groups)}")

    # 6. Normalise: log1p(TP10K).
    logger.info("Normalising (log1p TP10K)...")
    X = adata.X
    X_dense = X.toarray().astype(np.float64) if hasattr(X, "toarray") else X.astype(np.float64)
    row_sums = np.maximum(X_dense.sum(axis=1), 1e-12)
    X_norm = np.log1p(X_dense / row_sums[:, None] * 10000)
    del X_dense
    gc.collect()

    X_cand = X_norm[:, candidate_indices].copy()
    logger.info(f"Candidate expression matrix: {X_cand.shape}")
    control_mask = adata.obs["is_control"].values
    group_indices = {cat: np.where(info["mask"])[0] for cat, info in groups.items()}
    del adata
    gc.collect()

    # 7. Global whitening on candidate columns.
    ctrl_expr_cand = X_cand[control_mask]
    logger.info(f"Control expression for whitening: {ctrl_expr_cand.shape}")
    W = global_whitening(ctrl_expr_cand)
    logger.info("Applying whitening to all cells...")
    X_cand_w = X_cand @ W
    ctrl_expr_w = X_cand_w[control_mask]
    n_ctrl = ctrl_expr_w.shape[0]
    del X_cand, X_norm, ctrl_expr_cand
    gc.collect()

    # 8. Bootstrap evaluation.
    n_bootstrap = 100
    seed = 42
    rng = np.random.RandomState(seed)
    gene_to_idx = {g: i for i, g in enumerate(candidate_genes)}

    per_group = {}
    ranks_phi = []

    t_start = time.time()
    for i, (cat, info) in enumerate(groups.items()):
        if (i + 1) % 20 == 0 or i < 3:
            elapsed = time.time() - t_start
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta_min = (len(groups) - i - 1) / rate / 60 if rate > 0 else 999
            logger.info(f"  [{i + 1}/{len(groups)}] {rate:.2f}/s, ETA {eta_min:.0f} min")
            sys.stdout.flush()

        cell_idx = group_indices[cat]
        n_pert = len(cell_idx)
        pert_expr_w = X_cand_w[cell_idx]
        n_ctrl_sub = min(max(n_pert, 100), n_ctrl)
        valid_targets = [t for t in info["targets"] if t in gene_to_idx]
        if not valid_targets:
            continue

        boot_phi = []
        for _ in range(n_bootstrap):
            idx_p = rng.choice(n_pert, n_pert, replace=True)
            idx_c = rng.choice(n_ctrl, n_ctrl_sub, replace=True)
            lm_p = compute_log_matrix(pert_expr_w[idx_p])
            lm_c = compute_log_matrix(ctrl_expr_w[idx_c])
            delta = lm_p - lm_c
            boot_phi.append(hodge_edge_weight_phi(delta))

        mean_phi = np.mean(boot_phi, axis=0)
        # Ascending phi = upstream; lowest phi = predicted source.
        ranking = np.argsort(mean_phi)
        best_rank = min(
            int(np.where(ranking == gene_to_idx[t])[0][0]) + 1
            for t in valid_targets
        )
        per_group[cat] = {
            "n_cells": n_pert,
            "valid_targets": valid_targets,
            "best_rank": best_rank,
            "top1": best_rank == 1,
            "top5": best_rank <= 5,
            "top10": best_rank <= 10,
        }
        ranks_phi.append(best_rank)

    # 9. Summary.
    ranks_phi = np.asarray(ranks_phi)
    total = len(ranks_phi)
    top1 = int((ranks_phi == 1).sum())
    top5 = int((ranks_phi <= 5).sum())
    top10 = int((ranks_phi <= 10).sum())
    summary = {
        "total_groups": total,
        "top1_hits": top1,
        "top1_rate": float(top1 / total) if total else 0.0,
        "top5_hits": top5,
        "top5_rate": float(top5 / total) if total else 0.0,
        "top10_hits": top10,
        "top10_rate": float(top10 / total) if total else 0.0,
        "median_rank": float(np.median(ranks_phi)) if total else None,
        "mean_rank": float(np.mean(ranks_phi)) if total else None,
        "n_candidates": int(n_cand),
        "n_bootstrap": n_bootstrap,
        "seed": seed,
        "data_path": str(DATA_PATH),
    }

    print("\n" + "=" * 70)
    print(f"IDS R2b on Norman -- Top-1/5/10 ({total} groups, K_test = {n_cand})")
    print("=" * 70)
    print(f"  Top-1     : {top1}/{total} = {summary['top1_rate']:.1%}")
    print(f"  Top-5     : {top5}/{total} = {summary['top5_rate']:.1%}")
    print(f"  Top-10    : {top10}/{total} = {summary['top10_rate']:.1%}")
    print(f"  Median rk : {summary['median_rank']:.1f}")
    print(f"  Mean   rk : {summary['mean_rank']:.1f}")
    print()
    print("Reference (paper Section 2.4, Norman, same code path):")
    print("  IDS R2b (whitened + Hodge) : ~66% Top-1, median rank 1.0")
    print("  d_corr  (whitened, no Hodge): ~57% Top-1, median rank 1.0")

    out_json = OUTPUT_DIR / "ids_norman.json"
    with open(out_json, "w") as f:
        json.dump({"summary": summary, "per_group": per_group}, f, indent=2)
    logger.info(f"Saved summary to {out_json}")

    rows = []
    for cat, res in per_group.items():
        rows.append({
            "perturbation": cat,
            "n_cells": res["n_cells"],
            "targets": ",".join(res["valid_targets"]),
            "best_rank": res["best_rank"],
            "top1": res["top1"],
            "top5": res["top5"],
            "top10": res["top10"],
        })
    df = pd.DataFrame(rows)
    csv_path = OUTPUT_DIR / "ids_norman_per_group.csv"
    df.to_csv(csv_path, index=False)
    logger.info(f"Saved per-group CSV to {csv_path}")


if __name__ == "__main__":
    main()
