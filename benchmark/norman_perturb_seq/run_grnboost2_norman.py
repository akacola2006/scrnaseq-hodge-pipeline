"""
GRNBoost2-equivalent benchmark on Norman Perturb-seq data.
Uses sklearn GradientBoostingRegressor (same algorithm as GRNBoost2)
since arboreto has dask compatibility issues.

Output: gene importance -> 3 ranking strategies -> Top-1/5/10 hit rates.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from collections import defaultdict
from sklearn.ensemble import GradientBoostingRegressor


def load_norman_data(data_dir):
    """Load Norman Perturb-seq data."""
    import anndata as ad
    adata = sc.read_h5ad(os.path.join(data_dir, "norman_raw_counts.h5ad"))
    with open(os.path.join(data_dir, "candidate_genes.txt")) as f:
        candidate_genes = [line.strip() for line in f if line.strip()]
    return adata, candidate_genes


def get_perturbation_groups(adata):
    """Extract perturbation groups with their target genes."""
    groups = {}
    for pert in adata.obs["perturbation"].unique():
        if pert == "control":
            continue
        # Parse targets
        if "+" in pert:
            targets = pert.split("+")
        elif "_" in pert:
            targets = pert.split("_")
        else:
            targets = [pert]
        targets = [t for t in targets if t != "NegCtrl" and len(t) > 1]
        if len(targets) > 0:
            n_cells = (adata.obs["perturbation"] == pert).sum()
            if n_cells >= 3:
                groups[pert] = {"targets": targets, "n_cells": n_cells}
    return groups


def run_grnboost2_equivalent(expr_matrix, gene_names, seed=42):
    """
    GRNBoost2-equivalent: for each gene, fit GBM with all other genes as features.
    Returns importance matrix (n_genes x n_genes).

    This is the same algorithm as GRNBoost2/GENIE3 but without dask.
    GRNBoost2 parameters from arboreto defaults.
    """
    n_genes = len(gene_names)
    importance_matrix = np.zeros((n_genes, n_genes))

    for j in range(n_genes):
        # Target gene j, features = all other genes
        y = expr_matrix[:, j]
        feature_idx = [i for i in range(n_genes) if i != j]
        X = expr_matrix[:, feature_idx]

        if y.std() < 1e-10:
            continue

        # GRNBoost2 uses stochastic GBM with early stopping
        gbm = GradientBoostingRegressor(
            n_estimators=100,
            max_features="sqrt",
            learning_rate=0.05,
            subsample=0.9,
            max_depth=3,
            random_state=seed,
        )
        gbm.fit(X, y)

        # Feature importances
        for k, fi in enumerate(feature_idx):
            importance_matrix[fi, j] = gbm.feature_importances_[k]

    return importance_matrix


def importance_to_rankings(importance_matrix, gene_names):
    """Convert importance matrix to gene rankings using 3 strategies."""
    rankings = {}

    # (a) row-sum: total outgoing regulation
    row_sums = importance_matrix.sum(axis=1)
    rankings["row_sum"] = np.argsort(-row_sums)

    # (b) max: maximum outgoing importance
    row_maxs = importance_matrix.max(axis=1)
    rankings["max"] = np.argsort(-row_maxs)

    # (c) in-degree: total incoming regulation
    col_sums = importance_matrix.sum(axis=0)
    rankings["in_degree"] = np.argsort(-col_sums)

    return rankings


def evaluate_rankings(rankings, candidate_genes, target_genes):
    """Evaluate Top-1, Top-5, Top-10 for each ranking strategy."""
    results = {}
    gene_to_idx = {g: i for i, g in enumerate(candidate_genes)}

    for strategy, ranked_indices in rankings.items():
        best_rank = float("inf")
        for tg in target_genes:
            if tg in gene_to_idx:
                idx = gene_to_idx[tg]
                rank = np.where(ranked_indices == idx)[0]
                if len(rank) > 0:
                    best_rank = min(best_rank, rank[0] + 1)

        if best_rank == float("inf"):
            best_rank = len(candidate_genes)

        results[strategy] = {
            "best_rank": int(best_rank),
            "top1": best_rank == 1,
            "top5": best_rank <= 5,
            "top10": best_rank <= 10,
        }

    return results


def main():
    print("=" * 60)
    print("GRNBoost2-equivalent Benchmark on Norman Perturb-seq")
    print("(sklearn GBM implementation, same algorithm as GRNBoost2)")
    print("=" * 60)

    data_dir = str(
        Path(__file__).parent.parent / "cellnavi_ids" / "data" / "norman" / "processed"
    )
    output_dir = str(Path(__file__).parent / "results")
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    print("\nLoading data...")
    adata, candidate_genes = load_norman_data(data_dir)
    n_candidates = len(candidate_genes)
    print(f"Candidates: {n_candidates} genes")

    # Normalize: log1p(TP10K)
    print("Normalizing (log1p TP10K)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    gene_names = adata.var_names.tolist()
    candidate_idx = [gene_names.index(g) for g in candidate_genes if g in gene_names]
    available_candidates = [gene_names[i] for i in candidate_idx]
    print(f"Available candidates in data: {len(available_candidates)}")

    # Get control cells
    control_mask = adata.obs["perturbation"] == "control"
    control_expr_full = adata[control_mask].X
    if hasattr(control_expr_full, "toarray"):
        control_expr_full = control_expr_full.toarray()

    # Subset to candidate genes
    control_expr = control_expr_full[:, candidate_idx]
    print(f"Control cells: {control_expr.shape[0]}")

    # Get perturbation groups
    groups = get_perturbation_groups(adata)
    valid_groups = {}
    for pert, info in groups.items():
        valid_targets = [t for t in info["targets"] if t in available_candidates]
        if valid_targets:
            valid_groups[pert] = {**info, "valid_targets": valid_targets}
    print(f"Valid perturbation groups: {len(valid_groups)}")

    # Run per group
    all_results = {}
    strategy_hits = defaultdict(
        lambda: {"top1": 0, "top5": 0, "top10": 0, "total": 0, "ranks": []}
    )

    for i, (pert, info) in enumerate(valid_groups.items()):
        print(
            f"\r[{i+1}/{len(valid_groups)}] {pert} (n={info['n_cells']})...",
            end="",
            flush=True,
        )

        pert_mask = adata.obs["perturbation"] == pert
        pert_expr_full = adata[pert_mask].X
        if hasattr(pert_expr_full, "toarray"):
            pert_expr_full = pert_expr_full.toarray()
        pert_expr = pert_expr_full[:, candidate_idx]

        try:
            # Subsample control
            n_ctrl = min(max(info["n_cells"], 100), control_expr.shape[0])
            rng = np.random.RandomState(42)
            ctrl_idx = rng.choice(control_expr.shape[0], n_ctrl, replace=False)
            ctrl_sub = control_expr[ctrl_idx]

            # Combine and run GBM
            combined = np.vstack([pert_expr, ctrl_sub])
            importance_matrix = run_grnboost2_equivalent(
                combined, available_candidates, seed=42
            )
            rankings = importance_to_rankings(importance_matrix, available_candidates)
            eval_result = evaluate_rankings(
                rankings, available_candidates, info["valid_targets"]
            )

            all_results[pert] = eval_result

            for strategy in eval_result:
                strategy_hits[strategy]["total"] += 1
                if eval_result[strategy]["top1"]:
                    strategy_hits[strategy]["top1"] += 1
                if eval_result[strategy]["top5"]:
                    strategy_hits[strategy]["top5"] += 1
                if eval_result[strategy]["top10"]:
                    strategy_hits[strategy]["top10"] += 1
                strategy_hits[strategy]["ranks"].append(
                    eval_result[strategy]["best_rank"]
                )

        except Exception as e:
            print(f"\n  ERROR on {pert}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print("\n")

    # Summary
    print("=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)

    summary = {}
    for strategy in ["row_sum", "max", "in_degree"]:
        hits = strategy_hits[strategy]
        total = hits["total"]
        ranks = hits["ranks"]
        result = {
            "total_groups": total,
            "top1_hits": hits["top1"],
            "top1_rate": hits["top1"] / total if total > 0 else 0,
            "top5_hits": hits["top5"],
            "top5_rate": hits["top5"] / total if total > 0 else 0,
            "top10_hits": hits["top10"],
            "top10_rate": hits["top10"] / total if total > 0 else 0,
            "median_rank": float(np.median(ranks)) if ranks else None,
            "mean_rank": float(np.mean(ranks)) if ranks else None,
        }
        summary[strategy] = result

        print(f"\nStrategy: {strategy}")
        print(f"  Total groups: {total}")
        print(f"  Top-1: {hits['top1']}/{total} = {result['top1_rate']:.1%}")
        print(f"  Top-5: {hits['top5']}/{total} = {result['top5_rate']:.1%}")
        print(f"  Top-10: {hits['top10']}/{total} = {result['top10_rate']:.1%}")
        print(f"  Median rank: {result['median_rank']}")
        print(f"  Mean rank: {result['mean_rank']:.1f}")

    # Comparison
    print("\n" + "=" * 60)
    print("COMPARISON WITH IDS")
    print("=" * 60)
    print(f"IDS R2b Top-1 (232-group): 65.9%")
    print(f"DE ranking Top-1 (232-group): 81.0%")
    print(f"d_corr Top-1 (232-group): 57.3%")
    best_strategy = max(summary, key=lambda s: summary[s]["top1_rate"])
    print(
        f"GRNBoost2 best Top-1 ({best_strategy}): {summary[best_strategy]['top1_rate']:.1%}"
    )

    # Save
    output_file = os.path.join(output_dir, "grnboost2_norman_results.json")
    with open(output_file, "w") as f:
        json.dump(
            {"summary": summary, "per_group": all_results},
            f,
            indent=2,
        )
    print(f"\nResults saved to {output_file}")

    # CSV
    rows = []
    for pert, eval_result in all_results.items():
        row = {
            "perturbation": pert,
            "targets": ",".join(valid_groups[pert]["valid_targets"]),
        }
        for strategy in ["row_sum", "max", "in_degree"]:
            if strategy in eval_result:
                row[f"{strategy}_rank"] = eval_result[strategy]["best_rank"]
                row[f"{strategy}_top1"] = eval_result[strategy]["top1"]
        rows.append(row)
    df_out = pd.DataFrame(rows)
    csv_file = os.path.join(output_dir, "grnboost2_norman_per_group.csv")
    df_out.to_csv(csv_file, index=False)
    print(f"Per-group results saved to {csv_file}")


if __name__ == "__main__":
    main()
