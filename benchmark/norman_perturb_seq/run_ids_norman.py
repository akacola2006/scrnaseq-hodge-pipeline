"""
IDS R2b evaluation on Replogle K562 essential Perturb-seq.
Faithful replication of the R2b pipeline from rebuild_faithful.py.

Key components faithfully ported from R2b:
  1. Global whitening: LedoitWolf Σ_ctrl → W = V D^{-1/2} V^T → X @ W
  2. SPD matrix log: LedoitWolf → corr → eigenvalue clamping → matrix log
  3. Edge-weight Hodge flow: f(i,j) = |Δ_ij| * sign(d_i - d_j)
  4. Ranking by ascending phi (lowest phi = upstream driver)
  5. Control subsample: min(max(n_pert, 100), n_ctrl) per bootstrap
  6. Bootstrap with replacement (n_bootstrap=100)
"""

import json
import gc
import logging
import time
import sys
import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path
from sklearn.covariance import LedoitWolf
from sklearn.metrics import accuracy_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger(__name__)

DATA_DIR = Path("D:/Projects/cell navi/cellnavi_ids/data/replogle")
OUTPUT_DIR = Path("D:/Projects/cell navi/benchmark/results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Core R2b components (from rebuild_faithful.py)
# ---------------------------------------------------------------------------

def compute_log_matrix(X_norm):
    """
    LedoitWolf covariance → correlation → SPD guarantee → matrix log.
    Exact port of _compute_log_matrix_v2 (use_precision=False, extract_idx=None).
    """
    lw = LedoitWolf()
    cov = lw.fit(X_norm).covariance_

    # Cov → Corr
    d = np.sqrt(np.maximum(np.diag(cov), 1e-30))
    mat = cov / np.outer(d, d)
    np.fill_diagonal(mat, 1.0)
    mat = (mat + mat.T) / 2

    # SPD guarantee + matrix log via eigendecomposition
    eigvals, eigvecs = np.linalg.eigh(mat)
    eigvals = np.maximum(eigvals, 1e-15)
    log_mat = (eigvecs * np.log(eigvals)) @ eigvecs.T
    return (log_mat + log_mat.T) / 2


def hodge_edge_weight(delta):
    """
    Hodge decomposition with edge-weight flow on K_N.
    Exact port of _hodge_edge_weight from rebuild_faithful.py.

    Flow: f(i,j) = |Δ_ij| * sign(d_i - d_j)
    Divergence: div_i = -sum_{j≠i} |Δ_ij| * sign(d_i - d_j)
    Potential: φ = (div - mean(div)) / N

    Returns phi (ascending = upstream driver) and gradient fraction.
    """
    N = delta.shape[0]
    d = np.linalg.norm(delta, axis=1)  # per-gene correlation-change magnitude

    abs_delta = np.abs(delta).copy()
    np.fill_diagonal(abs_delta, 0)
    sign_d = np.sign(d[:, None] - d[None, :])

    # Analytical divergence for edge_weight flow on K_N
    div = -(abs_delta * sign_d).sum(axis=1)

    # Hodge potential
    phi = (div - div.mean()) / N

    # Gradient fraction
    flow_mat = abs_delta * np.abs(sign_d)
    flow_energy = np.sum(np.triu(flow_mat ** 2, k=1))
    phi_diff = phi[None, :] - phi[:, None]
    grad_energy = np.sum(np.triu(phi_diff ** 2, k=1))
    gf = grad_energy / flow_energy if flow_energy > 1e-30 else 0.0

    return phi, float(gf)


def global_whitening(ctrl_expr):
    """
    Compute global whitening matrix from all control cells.
    Exact port of R2b whitening (rebuild_faithful.py L999-1009).

    W = V D^{-1/2} V^T  (symmetric ZCA whitening)
    """
    logger.info(f"  Computing global whitening matrix (n={ctrl_expr.shape[0]}, "
                f"p={ctrl_expr.shape[1]})...")
    lw = LedoitWolf()
    cov = lw.fit(ctrl_expr).covariance_
    eigvals, eigvecs = np.linalg.eigh(cov)
    eigvals = np.maximum(eigvals, 1e-10)
    cond = eigvals[-1] / eigvals[0]
    W = (eigvecs / np.sqrt(eigvals)) @ eigvecs.T
    logger.info(f"  Whitening done: cond(Σ)={cond:.0f}")
    return W


# ---------------------------------------------------------------------------
# Main evaluation
# ---------------------------------------------------------------------------

def load_replogle_singlecell():
    """Load Replogle K562 essential single-cell data."""
    sc_path = DATA_DIR / "K562_essential_raw_singlecell_01.h5ad"
    if not sc_path.exists():
        raise FileNotFoundError(
            f"{sc_path} not found. Download from Figshare first:\n"
            "https://ndownloader.figshare.com/files/35773219"
        )
    logger.info(f"Loading {sc_path} ...")
    adata = ad.read_h5ad(sc_path)
    logger.info(f"Loaded: {adata.shape}")
    return adata


def main():
    print("=" * 60)
    print("IDS R2b Replogle K562 External Replication")
    print("  (faithful port of rebuild_faithful.py R2b pipeline)")
    print("=" * 60)

    # 1. Load data
    adata = load_replogle_singlecell()

    # Gene names for candidate matching
    if 'gene_name' in adata.var.columns:
        gene_names = adata.var['gene_name'].tolist()
    else:
        gene_names = adata.var_names.tolist()

    # 2. Parse perturbation info
    logger.info("Parsing perturbation info...")
    pert_col = None
    for col in ['gene', 'perturbation', 'guide_identity', 'target_gene']:
        if col in adata.obs.columns:
            pert_col = col
            break
    if pert_col is None:
        raise ValueError("No perturbation column found")

    logger.info(f"Using perturbation column: {pert_col}")

    # Identify control and perturbed cells
    control_values = ['control', 'non-targeting', 'NT', 'non_targeting']
    adata.obs['is_control'] = adata.obs[pert_col].isin(control_values)

    n_control = adata.obs['is_control'].sum()
    n_perturbed = (~adata.obs['is_control']).sum()
    logger.info(f"Control cells: {n_control}, Perturbed cells: {n_perturbed}")

    # 3. Build candidate gene list
    pert_genes = adata.obs.loc[~adata.obs['is_control'], pert_col].unique()
    gene_set = set(gene_names)
    candidate_genes = sorted([g for g in pert_genes if g in gene_set])
    candidate_indices = np.array([gene_names.index(g) for g in candidate_genes])
    n_cand = len(candidate_genes)

    logger.info(f"Perturbation targets: {len(pert_genes)}")
    logger.info(f"Candidate genes (in expression matrix): {n_cand}")

    # 4. Normalize (TP10K + log1p)
    logger.info("Normalizing (log1p TP10K from raw counts)...")
    X = adata.X
    if hasattr(X, 'toarray'):
        X_dense = X.toarray().astype(np.float64)
    else:
        X_dense = X.astype(np.float64)

    row_sums = X_dense.sum(axis=1)
    row_sums = np.maximum(row_sums, 1e-12)
    X_norm = np.log1p(X_dense / row_sums[:, None] * 10000)
    del X_dense
    gc.collect()

    # Extract candidate gene columns only
    X_cand = X_norm[:, candidate_indices].copy()
    logger.info(f"Candidate expression matrix: {X_cand.shape}")

    # Build masks and indices
    control_mask = adata.obs['is_control'].values
    pert_col_values = adata.obs[pert_col].values

    # Build per-group cell indices
    group_indices = {}
    for i, val in enumerate(pert_col_values):
        if not control_mask[i]:
            if val not in group_indices:
                group_indices[val] = []
            group_indices[val].append(i)

    del adata
    gc.collect()

    # 5. Global whitening on candidate genes (from all control cells)
    ctrl_expr_cand = X_cand[control_mask]
    logger.info(f"Control expression for whitening: {ctrl_expr_cand.shape}")

    W = global_whitening(ctrl_expr_cand)

    # Apply whitening to ALL expression data (candidate genes)
    logger.info("Applying whitening to all cells...")
    X_cand_w = X_cand @ W
    ctrl_expr_w = X_cand_w[control_mask]
    logger.info(f"Whitened candidate matrix: {X_cand_w.shape}")

    del X_cand, X_norm, ctrl_expr_cand
    gc.collect()

    # 6. Evaluate per perturbation group
    n_ctrl = ctrl_expr_w.shape[0]
    n_bootstrap = 100
    seed = 42
    rng = np.random.RandomState(seed)

    unique_perts = sorted([g for g in pert_genes if g in candidate_genes])
    logger.info(f"\nEvaluating {len(unique_perts)} perturbation groups "
                f"(n_bootstrap={n_bootstrap})...")

    all_ranks = []
    all_gf = []
    per_group_results = {}
    failed = []

    t_start = time.time()

    for i, pert_gene in enumerate(unique_perts):
        if (i + 1) % 10 == 0 or i < 5:
            elapsed = time.time() - t_start
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (len(unique_perts) - i - 1) / rate / 60 if rate > 0 else 999
            logger.info(f"  Progress: {i+1}/{len(unique_perts)} "
                        f"({rate:.2f}/s, ETA {eta:.0f} min)")
            sys.stdout.flush()
            sys.stderr.flush()

        if pert_gene not in group_indices:
            failed.append({'gene': pert_gene, 'reason': 'no cells found'})
            continue

        cell_idx = group_indices[pert_gene]
        n_pert = len(cell_idx)
        if n_pert < 3:
            failed.append({'gene': pert_gene, 'reason': f'too few cells ({n_pert})'})
            continue

        pert_expr_w = X_cand_w[cell_idx]

        # R2b bootstrap: control subsample = min(max(n_pert, 100), n_ctrl)
        n_ctrl_sub = min(max(n_pert, 100), n_ctrl)

        try:
            boot_phi = []
            boot_gf = []

            for _ in range(n_bootstrap):
                # Bootstrap with replacement (matching R2b)
                idx_p = rng.choice(n_pert, n_pert, replace=True)
                idx_c = rng.choice(n_ctrl, n_ctrl_sub, replace=True)

                # SPD log-correlation matrices
                lm_p = compute_log_matrix(pert_expr_w[idx_p])
                lm_c = compute_log_matrix(ctrl_expr_w[idx_c])

                # Delta + Hodge edge-weight decomposition
                delta = lm_p - lm_c
                phi_b, gf_b = hodge_edge_weight(delta)

                boot_phi.append(phi_b)
                boot_gf.append(gf_b)

            # Average phi across bootstrap
            mean_phi = np.mean(boot_phi, axis=0)
            mean_gf = float(np.mean(boot_gf))

            # Rank by ASCENDING phi (lowest = upstream driver) — R2b convention
            ranking = np.argsort(mean_phi)

            true_idx = candidate_genes.index(pert_gene)
            rank = int(np.where(ranking == true_idx)[0][0]) + 1  # 1-indexed

            all_ranks.append(rank)
            all_gf.append(mean_gf)

            per_group_results[pert_gene] = {
                'n_cells': n_pert,
                'rank': rank,
                'gradient_fraction': mean_gf,
                'rank_1_gene': candidate_genes[ranking[0]],
            }

        except Exception as e:
            failed.append({'gene': pert_gene, 'reason': str(e)})

    logger.info(f"\nCompleted: {len(all_ranks)}/{len(unique_perts)} groups")
    logger.info(f"Failed: {len(failed)} groups")

    # 7. Compute metrics
    all_ranks = np.array(all_ranks)

    top1_acc = float(np.mean(all_ranks == 1))
    top5_acc = float(np.mean(all_ranks <= 5))
    top10_acc = float(np.mean(all_ranks <= 10))
    top50_acc = float(np.mean(all_ranks <= 50))
    top100_acc = float(np.mean(all_ranks <= 100))
    median_rank = float(np.median(all_ranks))
    mean_rank = float(np.mean(all_ranks))
    mean_gf = float(np.mean(all_gf))

    random_top1 = 1.0 / n_cand

    metrics = {
        'dataset': 'Replogle K562 essential',
        'pipeline': 'IDS R2b (faithful port from rebuild_faithful.py)',
        'components': {
            'global_whitening': True,
            'spd_matrix_log': True,
            'edge_weight_flow': True,
            'ascending_phi_ranking': True,
            'control_subsample': 'min(max(n_pert, 100), n_ctrl)',
            'bootstrap_replace': True,
        },
        'n_perturbation_groups': len(unique_perts),
        'n_evaluated': len(all_ranks),
        'n_candidate_genes': n_cand,
        'n_control_cells': int(n_control),
        'n_bootstrap': n_bootstrap,
        'top1_accuracy': top1_acc,
        'top5_accuracy': top5_acc,
        'top10_accuracy': top10_acc,
        'top50_accuracy': top50_acc,
        'top100_accuracy': top100_acc,
        'median_rank': median_rank,
        'mean_rank': mean_rank,
        'mean_gf': mean_gf,
        'random_top1': float(random_top1),
        'fold_over_random_top1': float(top1_acc / random_top1) if random_top1 > 0 else 0,
        'n_failed': len(failed),
        'norman_reference': {
            'top1_accuracy': 0.664,
            'n_candidate_genes': 102,
            'random_top1': 1/102,
            'median_rank': 1.0,
            'mean_gf': 0.308,
        }
    }

    # 8. Save results
    results_file = OUTPUT_DIR / "ids_replogle_r2b_results.json"
    with open(results_file, 'w') as f:
        json.dump(metrics, f, indent=2)

    # Save per-group detail
    detail_file = OUTPUT_DIR / "ids_replogle_r2b_per_group.json"
    with open(detail_file, 'w') as f:
        json.dump(per_group_results, f, indent=2)

    # Print summary
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY (IDS R2b on Replogle)")
    print("=" * 60)
    print(f"\nPipeline: IDS R2b (faithful port)")
    print(f"  Global whitening:    YES (LedoitWolf ZCA)")
    print(f"  SPD matrix log:      YES")
    print(f"  Edge-weight flow:    YES")
    print(f"  Ranking:             ascending phi")
    print(f"  Control subsample:   min(max(n_pert,100), n_ctrl)")
    print(f"  Bootstrap:           {n_bootstrap} (replace=True)")
    print(f"\nDataset: Replogle K562 essential Perturb-seq")
    print(f"Perturbation groups evaluated: {len(all_ranks)}/{len(unique_perts)}")
    print(f"Candidate genes: {n_cand}")
    print(f"Control cells: {n_control}")
    print(f"\nIDS R2b Performance (Replogle):")
    print(f"  Top-1 accuracy:  {top1_acc:.3f} ({top1_acc*100:.1f}%)")
    print(f"  Top-5 accuracy:  {top5_acc:.3f} ({top5_acc*100:.1f}%)")
    print(f"  Top-10 accuracy: {top10_acc:.3f} ({top10_acc*100:.1f}%)")
    print(f"  Top-50 accuracy: {top50_acc:.3f} ({top50_acc*100:.1f}%)")
    print(f"  Top-100 accuracy:{top100_acc:.3f} ({top100_acc*100:.1f}%)")
    print(f"  Median rank:     {median_rank:.1f} / {n_cand}")
    print(f"  Mean rank:       {mean_rank:.1f} / {n_cand}")
    print(f"  Mean GF:         {mean_gf:.3f}")
    print(f"  Random Top-1:    {random_top1:.4f} ({random_top1*100:.2f}%)")
    print(f"  Fold over random: {top1_acc/random_top1:.1f}x" if random_top1 > 0 else "")
    print(f"\nNorman R2b Reference:")
    print(f"  Top-1 accuracy:  66.4%")
    print(f"  Median rank:     1.0")
    print(f"  Mean GF:         0.308")
    print(f"  Candidates:      102")
    print(f"\nResults saved to {results_file}")
    print(f"Per-group detail: {detail_file}")

    # 9. Failed groups summary
    if failed:
        print(f"\nFailed groups ({len(failed)}):")
        for f_info in failed[:10]:
            print(f"  {f_info['gene']}: {f_info['reason']}")
        if len(failed) > 10:
            print(f"  ... and {len(failed)-10} more")


if __name__ == "__main__":
    main()
