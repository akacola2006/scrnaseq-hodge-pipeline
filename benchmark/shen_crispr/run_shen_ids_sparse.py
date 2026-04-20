"""
run_track2_sparse.py
====================
P4 Track 2 on sparse k-NN graph.
Re-runs the CRISPR KO target gene phi percentile test
using k-NN Hodge instead of K_N Hodge.

Tests k = 17, 30, 50, and K_N (complete graph) as reference.

Usage:
    python run_track2_sparse.py
"""

import gzip
import json
import logging
import time
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import scipy.io
import scipy.sparse as sp
from scipy import stats as sp_stats
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import lsqr
from scipy.spatial.distance import pdist, squareform
from sklearn.covariance import LedoitWolf

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

# ── Paths ──
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "gse274058_extract"
OUTPUT_DIR = BASE_DIR / "p4_results"
OUTPUT_DIR.mkdir(exist_ok=True)

SAMPLES = [
    "GSM8442771_B1_1", "GSM8442772_B1_2", "GSM8442773_B2_1",
    "GSM8442774_B2_2", "GSM8442775_B3_1",
]

PERTURB = {
    "mSafe_gRNA": ("mSafe", "Ctrl(neg)"),
    "Trem2_gRNA": ("Trem2", "AD"),
    "Rraga_gRNA": ("Rraga", "Control"),
    "Fasn_gRNA": ("Fasn", "Control"),
    "Clu_gRNA": ("Clu", "AD"),
    "Dpp6_gRNA": ("Dpp6", "ALS"),
    "Tbk1_gRNA": ("Tbk1", "ALS"),
    "Flcn_gRNA": ("Flcn", "Control"),
    "Gfap_gRNA": ("Gfap", "Control"),
    "C9orf72_gRNA": ("C9orf72", "ALS"),
    "Cfap410_gRNA": ("Cfap410", "ALS"),
    "Stk39_gRNA": ("Stk39", "PD"),
    "Lrrk2_gRNA": ("Lrrk2", "PD"),
    "Ndufaf_gRNA": ("Ndufaf2", "AD"),
    "Sh3gl2_gRNA": ("Sh3gl2", "PD"),
    "Srf_gRNA": ("Srf", "Control"),
    "Rbfox_gRNA": ("Rbfox1", "AD"),
    "Olig2_gRNA": ("Olig2", "Control"),
}

N_HVG = 300
N_HVG_POOL = 5000
SEED = 42
MAX_CTRL_PER_LIB = 2000
MIN_KO_CELLS = 15
K_VALUES = [17, 30, 50]  # sparse k values; K_N added as special case


# =====================================================================
# Core functions (from track2_gene_hodge.py + sparse Hodge)
# =====================================================================

def estimate_correlation(expr: np.ndarray) -> np.ndarray:
    lw = LedoitWolf().fit(expr)
    cov = lw.covariance_
    std = np.sqrt(np.maximum(np.diag(cov), 1e-12))
    corr = cov / np.outer(std, std)
    corr = (corr + corr.T) / 2.0
    np.fill_diagonal(corr, 1.0)
    eigvals, eigvecs = np.linalg.eigh(corr)
    eigvals = np.maximum(eigvals, 1e-10)
    corr = (eigvecs * eigvals) @ eigvecs.T
    corr = (corr + corr.T) / 2.0
    np.fill_diagonal(corr, 1.0)
    return corr


def spd_log(A: np.ndarray) -> np.ndarray:
    A_sym = (A + A.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(A_sym)
    eigvals = np.maximum(eigvals, 1e-15)
    log_A = eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return (log_A + log_A.T) / 2.0


def build_knn_graph(delta: np.ndarray, k: int) -> np.ndarray:
    """Build symmetric k-NN graph. Returns upper-triangle edge list."""
    n = delta.shape[0]
    row_dist = squareform(pdist(delta, metric="euclidean"))
    adjacency = lil_matrix((n, n), dtype=np.float64)
    for i in range(n):
        dists = row_dist[i].copy()
        dists[i] = np.inf
        neighbours = np.argsort(dists)[:k]
        for j in neighbours:
            adjacency[i, j] = 1
            adjacency[j, i] = 1
    adj_csr = adjacency.tocsr()
    rows_arr, cols_arr = adj_csr.nonzero()
    mask = rows_arr < cols_arr
    return np.column_stack([rows_arr[mask], cols_arr[mask]])


def build_incidence(edges: np.ndarray, n_nodes: int) -> csr_matrix:
    n_edges = len(edges)
    rows = np.concatenate([np.arange(n_edges), np.arange(n_edges)])
    cols = np.concatenate([edges[:, 0], edges[:, 1]])
    vals = np.concatenate([-np.ones(n_edges), np.ones(n_edges)])
    return csr_matrix((vals, (rows, cols)), shape=(n_edges, n_nodes))


def sparse_hodge(edges: np.ndarray, flow: np.ndarray, n_nodes: int) -> Dict[str, Any]:
    B0 = build_incidence(edges, n_nodes)
    div = B0.T @ flow
    L = B0.T @ B0
    phi_raw, *_ = lsqr(L, div, atol=1e-12, btol=1e-12)
    phi = phi_raw - phi_raw.mean()
    grad = B0 @ phi
    flow_sq = float(np.dot(flow, flow))
    grad_sq = float(np.dot(grad, grad))
    gf = grad_sq / flow_sq if flow_sq > 1e-30 else 0.0
    return {"phi": phi, "gf": gf}


def hodge_gradient_kn(flow: np.ndarray, n_nodes: int) -> Dict[str, Any]:
    """Analytical Hodge on K_N (reference)."""
    edges = np.array(list(combinations(range(n_nodes), 2)), dtype=np.int32)
    div = np.zeros(n_nodes, dtype=np.float64)
    np.add.at(div, edges[:, 0], -flow)
    np.add.at(div, edges[:, 1], +flow)
    phi = (div - div.mean()) / n_nodes
    gradient = phi[edges[:, 1]] - phi[edges[:, 0]]
    flow_sq = np.dot(flow, flow)
    grad_sq = np.dot(gradient, gradient)
    gf = float(grad_sq / flow_sq) if flow_sq > 1e-30 else 0.0
    return {"phi": phi, "gf": gf}


def compute_delta(expr_ko: np.ndarray, expr_ctrl: np.ndarray) -> Optional[np.ndarray]:
    """Compute delta matrix (KO - Control in log-correlation space)."""
    try:
        corr_ctrl = estimate_correlation(expr_ctrl)
        corr_ko = estimate_correlation(expr_ko)
    except Exception as exc:
        logger.warning("Correlation estimation failed: %s", exc)
        return None
    return spd_log(corr_ko) - spd_log(corr_ctrl)


def run_hodge_on_delta(delta: np.ndarray, k: Optional[int] = None) -> Dict[str, Any]:
    """Run Hodge decomposition on delta. k=None means K_N."""
    n = delta.shape[0]
    d_corr = np.linalg.norm(delta, axis=1)  # d_mode=full always

    if k is None or k >= n - 1:
        # K_N: all edges
        n_edges = n * (n - 1) // 2
        flow = np.empty(n_edges, dtype=np.float64)
        idx = 0
        for i in range(n):
            js = np.arange(i + 1, n)
            w_ij = np.abs(delta[i, js])
            diff = d_corr[i] - d_corr[js]
            flow[idx:idx + len(js)] = w_ij * np.sign(diff)
            idx += len(js)
        result = hodge_gradient_kn(flow, n)
        result["n_edges"] = n_edges
        result["k"] = "KN"
    else:
        # k-NN sparse
        edges = build_knn_graph(delta, k)
        w_ij = np.abs(delta[edges[:, 0], edges[:, 1]])
        diff = d_corr[edges[:, 0]] - d_corr[edges[:, 1]]
        flow = w_ij * np.sign(diff)
        result = sparse_hodge(edges, flow, n)
        result["n_edges"] = len(edges)
        result["k"] = k

    return result


def get_target_percentile(phi: np.ndarray, target_idx: Optional[int]) -> Optional[float]:
    """Get target gene's percentile in phi ranking (higher = more upstream)."""
    if target_idx is None:
        return None
    n = len(phi)
    rank = int(np.sum(phi < phi[target_idx])) + 1  # 1 = lowest
    return 100.0 * rank / n


# =====================================================================
# Data Loading
# =====================================================================

def load_all_data():
    """Load expression data from all 5 libraries, identify perturbations."""
    logger.info("Loading data from %d libraries...", len(SAMPLES))

    # Read gene names from first library
    mdir = DATA_DIR / SAMPLES[0]
    with gzip.open(mdir / "features.tsv.gz", "rt") as f:
        features = [line.strip().split("\t") for line in f]
    gene_names = np.array([ft[1] for ft in features])
    n_total = len(gene_names)

    # Identify gRNA features
    grna_idx_map = {}
    for gn in PERTURB:
        idx = np.where(gene_names == gn)[0]
        if len(idx) == 1:
            grna_idx_map[gn] = idx[0]

    grna_set = set(grna_idx_map.values())
    myrf = np.where(gene_names == "Myrf_gRNA")[0]
    if len(myrf) > 0:
        grna_set.add(myrf[0])
    endo_mask = np.array([i not in grna_set for i in range(n_total)])
    endo_idx = np.where(endo_mask)[0]

    rng = np.random.default_rng(SEED)

    # Compute HVGs from first library
    mat0 = scipy.io.mmread(mdir / "matrix.mtx.gz").tocsc()
    sub = rng.choice(mat0.shape[1], size=min(2000, mat0.shape[1]), replace=False)
    X_sub = mat0[:, sub].T.toarray().astype(np.float64)[:, endo_mask]
    totals = X_sub.sum(axis=1, keepdims=True)
    totals = np.where(totals > 0, totals, 1.0)
    X_sub = np.log1p(X_sub / totals * 1e6)
    gene_var = X_sub.var(axis=0)

    hvg_local = np.argsort(gene_var)[-N_HVG:]
    hvg_local = np.sort(hvg_local)
    hvg_global = endo_idx[hvg_local]
    hvg_names = gene_names[hvg_global]

    # Ensure perturbation targets are included
    target_genes_to_add = {}
    for grna, (target, cat) in PERTURB.items():
        target_lower = target.lower()
        matches = [i for i, g in enumerate(gene_names[endo_idx]) if g.lower() == target_lower]
        if matches:
            global_idx = endo_idx[matches[0]]
            if global_idx not in hvg_global:
                target_genes_to_add[target] = global_idx

    if target_genes_to_add:
        extra = np.array(list(target_genes_to_add.values()))
        hvg_global = np.sort(np.concatenate([hvg_global, extra]))
        hvg_names = gene_names[hvg_global]
        logger.info("  Added %d target genes to HVG set, total: %d",
                     len(extra), len(hvg_global))

    # Load all libraries: identify KO cells per perturbation
    all_expr = {grna: [] for grna in PERTURB}
    all_ctrl = []

    for sample in SAMPLES:
        sdir = DATA_DIR / sample
        mat = scipy.io.mmread(sdir / "matrix.mtx.gz").tocsc()
        nc = mat.shape[1]

        # Read barcodes (not needed for KO identification)
        # Identify KO cells by gRNA expression
        ko_assignments = {}
        for grna, grna_global_idx in grna_idx_map.items():
            grna_expr = mat[grna_global_idx, :].toarray().ravel()
            ko_cells = np.where(grna_expr > 0)[0]
            if len(ko_cells) > 0:
                ko_assignments[grna] = ko_cells

        # All barcode-positive cells
        all_ko_cells = set()
        for cells in ko_assignments.values():
            all_ko_cells.update(cells)

        # Barcode-negative = control
        ctrl_cells = np.array([c for c in range(nc) if c not in all_ko_cells])
        if len(ctrl_cells) > MAX_CTRL_PER_LIB:
            ctrl_cells = rng.choice(ctrl_cells, size=MAX_CTRL_PER_LIB, replace=False)

        # Extract expression for HVG genes only
        X_ctrl = mat[hvg_global, :][:, ctrl_cells].T.toarray().astype(np.float64)
        totals_c = X_ctrl.sum(axis=1, keepdims=True)
        totals_c = np.where(totals_c > 0, totals_c, 1.0)
        X_ctrl = np.log1p(X_ctrl / totals_c * 1e6)
        all_ctrl.append(X_ctrl)

        for grna, ko_cells in ko_assignments.items():
            X_ko = mat[hvg_global, :][:, ko_cells].T.toarray().astype(np.float64)
            totals_k = X_ko.sum(axis=1, keepdims=True)
            totals_k = np.where(totals_k > 0, totals_k, 1.0)
            X_ko = np.log1p(X_ko / totals_k * 1e6)
            all_expr[grna].append(X_ko)

        del mat
        logger.info("  %s: %d ctrl, %d KO groups", sample, len(ctrl_cells), len(ko_assignments))

    # Concatenate
    ctrl_expr = np.vstack(all_ctrl)
    ko_expr_map = {}
    for grna in PERTURB:
        parts = all_expr[grna]
        if parts:
            ko_expr_map[grna] = np.vstack(parts)

    logger.info("Total ctrl: %d cells, %d genes", ctrl_expr.shape[0], ctrl_expr.shape[1])

    # Map target gene names to indices in hvg_names
    target_idx_map = {}
    for grna, (target, cat) in PERTURB.items():
        matches = [i for i, g in enumerate(hvg_names) if g.lower() == target.lower()]
        if matches:
            target_idx_map[grna] = matches[0]
        else:
            target_idx_map[grna] = None

    return hvg_names, ctrl_expr, ko_expr_map, target_idx_map


# =====================================================================
# Main
# =====================================================================

def main():
    t_start = time.time()
    logger.info("=" * 60)
    logger.info("P4 Track 2 Sparse: Causal Validation of k-NN Hodge phi")
    logger.info("=" * 60)

    hvg_names, ctrl_expr, ko_expr_map, target_idx_map = load_all_data()
    n_genes = len(hvg_names)

    all_results = []

    for grna, (target, category) in PERTURB.items():
        if grna not in ko_expr_map or ko_expr_map[grna].shape[0] < MIN_KO_CELLS:
            logger.info("  SKIP %s (n_ko=%d < %d)",
                         target, ko_expr_map.get(grna, np.empty((0,))).shape[0], MIN_KO_CELLS)
            continue

        ko_expr = ko_expr_map[grna]
        n_ko = ko_expr.shape[0]

        # Subsample control to ctrl_ratio × n_ko
        ctrl_ratio = 8
        n_ctrl_use = min(ctrl_ratio * n_ko, ctrl_expr.shape[0])
        rng = np.random.default_rng(SEED)
        ctrl_idx = rng.choice(ctrl_expr.shape[0], size=n_ctrl_use, replace=False)
        ctrl_sub = ctrl_expr[ctrl_idx]

        logger.info("--- %s (%s, n_ko=%d, n_ctrl=%d) ---", target, category, n_ko, n_ctrl_use)

        # Compute delta
        delta = compute_delta(ko_expr, ctrl_sub)
        if delta is None:
            continue

        target_gene_idx = target_idx_map.get(grna)

        # Run for each k value + K_N
        for k in K_VALUES + [None]:
            k_label = f"k={k}" if k is not None else "KN"
            hodge = run_hodge_on_delta(delta, k=k)
            target_pct = get_target_percentile(hodge["phi"], target_gene_idx)

            rec = {
                "gene": target,
                "category": category,
                "n_ko": n_ko,
                "n_genes": n_genes,
                "k": k if k is not None else "KN",
                "n_edges": hodge["n_edges"],
                "gf": hodge["gf"],
                "target_idx": target_gene_idx,
                "target_percentile": target_pct,
            }
            all_results.append(rec)

            if target_pct is not None:
                logger.info("  %s: pct=%.1f%%, GF=%.4f, |E|=%d",
                             k_label, target_pct, hodge["gf"], hodge["n_edges"])
            else:
                logger.info("  %s: target not in gene set, GF=%.4f", k_label, hodge["gf"])

    # ── Aggregate ──
    import pandas as pd
    df = pd.DataFrame(all_results)
    df.to_csv(OUTPUT_DIR / "track2_sparse_results.csv", index=False)

    # Summary per k
    summary = {}
    for k_val in K_VALUES + ["KN"]:
        sub = df[df["k"] == k_val]
        pcts = sub["target_percentile"].dropna().values

        if len(pcts) >= 3:
            # Wilcoxon: test if percentiles > 50
            stat, wilcox_p = sp_stats.wilcoxon(pcts - 50.0, alternative="greater")
            n_above = int(np.sum(pcts > 50))
            binom_p = float(sp_stats.binomtest(n_above, len(pcts), 0.5, alternative="greater").pvalue)

            summary[str(k_val)] = {
                "k": k_val,
                "n_perturbations": len(pcts),
                "mean_percentile": float(np.mean(pcts)),
                "std_percentile": float(np.std(pcts)),
                "median_percentile": float(np.median(pcts)),
                "n_above_50": n_above,
                "wilcoxon_p": float(wilcox_p),
                "binomial_p": binom_p,
                "mean_gf": float(sub["gf"].mean()),
            }
            logger.info("k=%s: mean_pct=%.1f%%, wilcox_p=%.4f, binom_p=%.4f, GF=%.4f",
                         k_val, np.mean(pcts), wilcox_p, binom_p, sub["gf"].mean())

    with open(OUTPUT_DIR / "track2_sparse_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Comparison table ──
    comparison_rows = []
    for grna, (target, cat) in PERTURB.items():
        row = {"gene": target, "category": cat}
        for k_val in K_VALUES + ["KN"]:
            sub = df[(df["gene"] == target) & (df["k"] == k_val)]
            if len(sub) > 0:
                row[f"pct_k{k_val}"] = sub.iloc[0]["target_percentile"]
                row[f"gf_k{k_val}"] = sub.iloc[0]["gf"]
        comparison_rows.append(row)

    comp_df = pd.DataFrame(comparison_rows)
    comp_df.to_csv(OUTPUT_DIR / "track2_sparse_vs_kn.csv", index=False)

    # ── Figure ──
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: mean target percentile by k
    k_labels = [str(k) for k in K_VALUES] + ["KN"]
    means = [summary.get(str(k), {}).get("mean_percentile", 0) for k in K_VALUES + ["KN"]]
    stds = [summary.get(str(k), {}).get("std_percentile", 0) for k in K_VALUES + ["KN"]]
    pvals = [summary.get(str(k), {}).get("wilcoxon_p", 1) for k in K_VALUES + ["KN"]]

    x = np.arange(len(k_labels))
    bars = axes[0].bar(x, means, yerr=stds, color=["steelblue"] * len(K_VALUES) + ["coral"],
                       edgecolor="black", linewidth=0.5, capsize=5)
    axes[0].axhline(50, color="grey", linestyle="--", linewidth=1, label="Null (50%)")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([f"k={k}" for k in K_VALUES] + ["K_N"])
    axes[0].set_ylabel("Mean Target Gene Percentile")
    axes[0].set_title("P4 Track 2: Sparse vs K_N")
    axes[0].legend()

    for i, p in enumerate(pvals):
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        axes[0].text(i, means[i] + stds[i] + 2, sig, ha="center", fontsize=11)

    # Right: GF by k
    gfs = [summary.get(str(k), {}).get("mean_gf", 0) for k in K_VALUES + ["KN"]]
    axes[1].bar(x, gfs, color=["steelblue"] * len(K_VALUES) + ["coral"],
                edgecolor="black", linewidth=0.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f"k={k}" for k in K_VALUES] + ["K_N"])
    axes[1].set_ylabel("Gradient Fraction")
    axes[1].set_title("GF: Sparse vs K_N")

    plt.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(OUTPUT_DIR / f"Fig_track2_sparse_comparison.{ext}",
                    dpi=200, bbox_inches="tight")
    plt.close(fig)

    dt = time.time() - t_start
    logger.info("=" * 60)
    logger.info("COMPLETE in %.1fs", dt)
    logger.info("=" * 60)

    # Print summary table
    print("\n" + "=" * 60)
    print("SUMMARY: Target Gene Percentile by Graph Type")
    print("=" * 60)
    for k_val in K_VALUES + ["KN"]:
        s = summary.get(str(k_val), {})
        print(f"  k={str(k_val):>3s}: mean={s.get('mean_percentile', 0):.1f}%, "
              f"wilcox_p={s.get('wilcoxon_p', 1):.4f}, "
              f"GF={s.get('mean_gf', 0):.4f}")
    print()


if __name__ == "__main__":
    main()
