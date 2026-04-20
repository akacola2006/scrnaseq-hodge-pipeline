"""
Step 8: Glioma pseudotime construction and validation
=====================================================
Constructs data-driven pseudotime for glioma bulk RNA-seq (723 samples)
using PT-A approach: expression PCA → diffusion map → DPT.

Then evaluates:
  1. Does pseudotime recover clinical grade ordering? (Spearman rho)
  2. AUC: Normal vs GBM separation by pseudotime
  3. Construct pseudotime-based windows (5 equal bins)
  4. Run Hodge decomposition on PT-based windows
  5. Compare PT-based phi to clinical-grade-based phi (Spearman rho)
  6. GSEA on PT-based phi to check if translation is upstream

Usage:
    python step8_pseudotime.py
"""
import json
import logging
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.covariance import LedoitWolf
from sklearn.metrics import roc_auc_score

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

BASE = Path(r"D:\Projects\IDS paper writing\glioma_hodge")
DATA_DIR = BASE / "data"
OUT_DIR = BASE / "results" / "pseudotime"


# ── Diffusion pseudotime ──────────────────────────────────────

def diffusion_pseudotime(X: np.ndarray, n_components: int = 30,
                         alpha: float = 1.0, seed: int = 42) -> np.ndarray:
    """
    Compute diffusion pseudotime from expression matrix.

    Steps:
      1. PCA to n_components
      2. Pairwise Euclidean distance in PCA space
      3. Gaussian kernel with adaptive bandwidth (median distance)
      4. Row-normalise → transition matrix
      5. Eigendecomposition → diffusion coordinates
      6. DPT = Euclidean distance from root in diffusion space

    Args:
        X: (n_samples, n_genes) expression matrix
        n_components: PCA dimensions
        alpha: diffusion map alpha parameter
        seed: random seed for PCA

    Returns:
        dpt: (n_samples,) pseudotime values
    """
    n_samples = X.shape[0]
    logger.info("  PCA: %d samples × %d genes → %d components",
                X.shape[0], X.shape[1], n_components)

    # Step 1: PCA
    pca = PCA(n_components=n_components, random_state=seed)
    X_pca = pca.fit_transform(X)
    var_explained = pca.explained_variance_ratio_.sum()
    logger.info("  PCA variance explained: %.1f%%", var_explained * 100)

    # Step 2: Pairwise distances
    from scipy.spatial.distance import pdist, squareform
    D = squareform(pdist(X_pca, metric='euclidean'))

    # Step 3: Gaussian kernel with median bandwidth
    sigma = np.median(D[D > 0])
    logger.info("  Kernel bandwidth (median distance): %.4f", sigma)
    K = np.exp(-D**2 / (2 * sigma**2))
    np.fill_diagonal(K, 0)  # no self-loops

    # Step 4: Density normalisation (alpha)
    if alpha > 0:
        q = K.sum(axis=1)
        q_outer = np.outer(q, q) ** alpha
        q_outer = np.maximum(q_outer, 1e-15)
        K = K / q_outer

    # Row normalise → transition matrix
    row_sums = K.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1e-15)
    T = K / row_sums

    # Step 5: Eigendecomposition
    # Use symmetric normalisation for numerical stability
    d_half = np.sqrt(K.sum(axis=1))
    d_half_inv = 1.0 / np.maximum(d_half, 1e-15)
    T_sym = (T * d_half.reshape(-1, 1)) * d_half_inv.reshape(1, -1)
    T_sym = (T_sym + T_sym.T) / 2  # ensure symmetry

    n_eig = min(20, n_samples - 1)
    eigvals, eigvecs = np.linalg.eigh(T_sym)

    # Sort descending
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Take top n_eig (excluding trivial eigenvalue ≈ 1)
    # Diffusion coordinates: psi_i = eigvec_i / eigvec_0
    diff_coords = eigvecs[:, 1:n_eig+1]  # skip first (trivial)
    diff_eigenvalues = eigvals[1:n_eig+1]

    logger.info("  Top 5 diffusion eigenvalues: %s",
                [f"{v:.4f}" for v in diff_eigenvalues[:5]])

    # Step 6: DPT from root
    # Root = sample with highest density (most central)
    density = K.sum(axis=1)
    root_idx = np.argmax(density)
    logger.info("  Root sample (highest density): index %d", root_idx)

    # DPT = L2 distance from root in diffusion coordinate space
    # Weight by eigenvalues for multiscale
    weighted_coords = diff_coords * diff_eigenvalues[np.newaxis, :]
    dpt = np.sqrt(np.sum((weighted_coords - weighted_coords[root_idx])**2, axis=1))

    return dpt


# ── SPD / Hodge utilities (from step2) ────────────────────────

def spd_log(A: np.ndarray) -> np.ndarray:
    A_sym = (A + A.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(A_sym)
    eigvals = np.maximum(eigvals, 1e-15)
    log_A = eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return (log_A + log_A.T) / 2.0


def compute_correlation_lw(X: np.ndarray) -> np.ndarray:
    lw = LedoitWolf()
    lw.fit(X)
    cov = lw.covariance_
    std = np.sqrt(np.diag(cov))
    std = np.maximum(std, 1e-15)
    corr = cov / np.outer(std, std)
    np.fill_diagonal(corr, 1.0)
    eigvals, eigvecs = np.linalg.eigh(corr)
    eigvals = np.maximum(eigvals, 1e-10)
    corr = (eigvecs * eigvals) @ eigvecs.T
    return (corr + corr.T) / 2.0


def hodge_gradient_kn(flow, n_nodes, edges):
    div = np.zeros(n_nodes, dtype=np.float64)
    np.add.at(div, edges[:, 0], -flow)
    np.add.at(div, edges[:, 1], +flow)
    phi = (div - div.mean()) / n_nodes
    gradient = phi[edges[:, 1]] - phi[edges[:, 0]]
    flow_sq = np.dot(flow, flow)
    grad_sq = np.dot(gradient, gradient)
    gf = float(grad_sq / flow_sq) if flow_sq > 1e-30 else 0.0
    return {"phi": phi, "gradient_fraction": gf}


def build_flow_edge_weight(delta):
    n_genes = delta.shape[0]
    gene_d = np.linalg.norm(delta, axis=1)
    n_edges = n_genes * (n_genes - 1) // 2
    flow = np.empty(n_edges, dtype=np.float64)
    idx = 0
    for i in range(n_genes):
        js = np.arange(i + 1, n_genes)
        delta_ij = np.abs(delta[i, js])
        diff = gene_d[i] - gene_d[js]
        flow[idx:idx + len(js)] = delta_ij * np.sign(diff)
        idx += len(js)
    flow_std = float(np.std(flow))
    if flow_std > 1e-15:
        flow = flow / flow_std
    return flow


# ── Main ──────────────────────────────────────────────────────

def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 8: Glioma pseudotime construction & validation")
    logger.info("=" * 60)

    # ── Load data ─────────────────────────────────────────────
    logger.info("Loading data...")
    expr = pd.read_parquet(DATA_DIR / "expression_gbmlgg_gtex.parquet")
    meta = pd.read_csv(DATA_DIR / "metadata_gbmlgg_gtex.csv", index_col=0)
    gene_names = list(expr.index)
    expr = expr.T  # → (samples × genes)

    common = expr.index.intersection(meta.index)
    expr = expr.loc[common]
    meta = meta.loc[common]
    logger.info("Samples: %d, Genes: %d", len(meta), expr.shape[1])

    # Exclude grade-unknown samples
    mask_known = meta["window_label"] != "LGG IDH-mut (grade unknown)"
    expr = expr.loc[mask_known]
    meta = meta.loc[mask_known]
    logger.info("After excluding grade-unknown: %d samples", len(meta))

    X = expr.values.astype(np.float64)

    # ── 1. Construct pseudotime ───────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("1. Constructing diffusion pseudotime (PT-A)")
    logger.info("=" * 60)

    dpt = diffusion_pseudotime(X, n_components=30, alpha=1.0, seed=42)
    meta = meta.copy()
    meta["pseudotime"] = dpt

    # ── 2. Compare pseudotime to clinical grade ───────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("2. Pseudotime vs clinical grade")
    logger.info("=" * 60)

    # Numeric grade: Normal=0, G2=1, G3=2, IDH-wt=3, GBM=4
    grade_map = meta["window"].values  # already 0-4
    rho, p_rho = stats.spearmanr(dpt, grade_map)
    logger.info("  Spearman rho(pseudotime, clinical grade): %.4f (p=%.4e)", rho, p_rho)

    # AUC: Normal (0) vs GBM (4)
    mask_normal = meta["window"] == 0
    mask_gbm = meta["window"] == 4
    mask_ng = mask_normal | mask_gbm
    y_ng = meta.loc[mask_ng, "window"].values
    y_binary = (y_ng == 4).astype(int)
    dpt_ng = dpt[mask_ng.values]
    auc_ng = roc_auc_score(y_binary, dpt_ng)
    # If AUC < 0.5, pseudotime is reversed relative to clinical grade
    if auc_ng < 0.5:
        auc_ng = 1 - auc_ng
        logger.info("  (Pseudotime direction reversed relative to grade)")
    logger.info("  AUC (Normal vs GBM): %.4f", auc_ng)

    # Per-window pseudotime statistics
    logger.info("")
    logger.info("  Per-window pseudotime (mean ± std):")
    for w in sorted(meta["window"].unique()):
        mask_w = meta["window"] == w
        pt_w = dpt[mask_w.values]
        label = meta.loc[mask_w, "window_label"].iloc[0]
        logger.info("    W%d (%s): %.4f ± %.4f (n=%d)",
                     w, label, pt_w.mean(), pt_w.std(), mask_w.sum())

    # ── 3. PT-based windows ───────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("3. Constructing pseudotime-based windows")
    logger.info("=" * 60)

    # 5 equal-sized bins by pseudotime quantile
    n_windows = 5
    pt_sorted_idx = np.argsort(dpt)
    n_total = len(dpt)
    window_size = n_total // n_windows
    pt_window = np.zeros(n_total, dtype=int)
    for w in range(n_windows):
        start = w * window_size
        end = (w + 1) * window_size if w < n_windows - 1 else n_total
        pt_window[pt_sorted_idx[start:end]] = w

    meta["pt_window"] = pt_window

    # Show composition of each PT window
    logger.info("  PT window composition (clinical labels):")
    for w in range(n_windows):
        mask_w = meta["pt_window"] == w
        counts = meta.loc[mask_w, "window_label"].value_counts()
        logger.info("    PT-W%d (n=%d):", w, mask_w.sum())
        for label, count in counts.items():
            logger.info("      %s: %d", label, count)

    # Cross-tabulation
    ct = pd.crosstab(meta["pt_window"], meta["window_label"])
    logger.info("")
    logger.info("  Cross-tabulation (PT window × clinical label):")
    logger.info("\n%s", ct.to_string())

    # ── 4. Hodge on PT-based windows ──────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("4. Hodge decomposition on PT-based windows")
    logger.info("=" * 60)

    n_genes = X.shape[1]
    edges_i, edges_j = np.triu_indices(n_genes, k=1)
    edges = np.stack([edges_i, edges_j], axis=1).astype(np.int32)

    # Compute per-PT-window correlation matrices
    pt_log_corr = {}
    for w in range(n_windows):
        mask_w = (meta["pt_window"] == w).values
        X_w = X[mask_w]
        logger.info("  PT-W%d: %d samples", w, X_w.shape[0])
        corr = compute_correlation_lw(X_w)
        pt_log_corr[w] = spd_log(corr)

    # Hodge on 4 consecutive transitions
    pt_phi_all = {}
    pt_results = {}
    for t in range(n_windows - 1):
        delta = pt_log_corr[t + 1] - pt_log_corr[t]
        flow = build_flow_edge_weight(delta)
        hodge = hodge_gradient_kn(flow, n_genes, edges)
        pt_phi_all[t] = hodge["phi"]
        pt_results[f"PT_T{t}_{t+1}"] = {
            "gradient_fraction": hodge["gradient_fraction"],
        }
        logger.info("  PT transition %d→%d: GF=%.4f", t, t+1,
                     hodge["gradient_fraction"])

    # Mean phi across PT transitions
    pt_phi_mean = np.mean([pt_phi_all[t] for t in range(n_windows - 1)], axis=0)

    # ── 5. Load clinical-grade phi and compare ────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("5. Comparing PT-based phi vs clinical-grade phi")
    logger.info("=" * 60)

    # Load clinical phi from step2 results
    clinical_phi_list = []
    for t_key in ["T0_1", "T1_2", "T2_3", "T3_4"]:
        phi_path = BASE / "results" / f"phi_{t_key}.npy"
        if phi_path.exists():
            clinical_phi_list.append(np.load(phi_path))
        else:
            logger.warning("  Missing: %s", phi_path)

    if len(clinical_phi_list) == 4:
        clinical_phi_mean = np.mean(clinical_phi_list, axis=0)

        rho_phi, p_phi = stats.spearmanr(pt_phi_mean, clinical_phi_mean)
        logger.info("  Spearman rho(PT phi, clinical phi): %.4f (p=%.4e)",
                     rho_phi, p_phi)

        # Per-transition comparison
        for t in range(4):
            r, p = stats.spearmanr(pt_phi_all[t], clinical_phi_list[t])
            logger.info("  Transition %d: rho=%.4f (p=%.4e)", t, r, p)
    else:
        rho_phi = np.nan
        p_phi = np.nan
        logger.warning("  Could not load all clinical phi files")

    # ── 6. Translation/ribosomal gene check ───────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("6. Translation gene position in PT-based phi")
    logger.info("=" * 60)

    # Ribosomal genes (RPL/RPS)
    ribo_mask = np.array([g.startswith("RPL") or g.startswith("RPS")
                          for g in gene_names])
    n_ribo = ribo_mask.sum()
    logger.info("  Ribosomal genes (RPL/RPS) in gene set: %d", n_ribo)

    if n_ribo > 0:
        # Percentile of ribosomal genes in mean PT phi
        ranks = stats.rankdata(pt_phi_mean)
        percentiles = ranks / len(ranks) * 100
        ribo_pctile = percentiles[ribo_mask].mean()
        logger.info("  Mean percentile of ribosomal genes (PT phi): %.1f",
                     ribo_pctile)

        # For comparison, in clinical phi
        if len(clinical_phi_list) == 4:
            ranks_clin = stats.rankdata(clinical_phi_mean)
            pctiles_clin = ranks_clin / len(ranks_clin) * 100
            ribo_pctile_clin = pctiles_clin[ribo_mask].mean()
            logger.info("  Mean percentile of ribosomal genes (clinical phi): %.1f",
                         ribo_pctile_clin)

        # Wilcoxon rank-sum: ribosomal vs non-ribosomal
        stat_w, p_w = stats.ranksums(pt_phi_mean[ribo_mask],
                                      pt_phi_mean[~ribo_mask])
        logger.info("  Wilcoxon rank-sum (ribo vs rest): z=%.3f, p=%.4e",
                     stat_w, p_w)

    # ── 7. Save results ───────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("7. Saving results")
    logger.info("=" * 60)

    results = {
        "pseudotime_vs_clinical_grade": {
            "spearman_rho": float(rho),
            "spearman_p": float(p_rho),
            "auc_normal_vs_gbm": float(auc_ng),
        },
        "phi_comparison": {
            "mean_phi_spearman_rho": float(rho_phi) if not np.isnan(rho_phi) else None,
            "mean_phi_spearman_p": float(p_phi) if not np.isnan(p_phi) else None,
        },
        "pt_hodge": pt_results,
        "ribosomal_genes": {
            "n_ribo": int(n_ribo),
            "mean_percentile_pt_phi": float(ribo_pctile) if n_ribo > 0 else None,
            "mean_percentile_clinical_phi": float(ribo_pctile_clin) if n_ribo > 0 and len(clinical_phi_list) == 4 else None,
            "wilcoxon_z": float(stat_w) if n_ribo > 0 else None,
            "wilcoxon_p": float(p_w) if n_ribo > 0 else None,
        },
        "n_samples": int(len(meta)),
        "n_genes": int(n_genes),
        "n_pt_windows": n_windows,
    }

    # Save
    with open(OUT_DIR / "pseudotime_results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Save pseudotime values
    pt_df = meta[["window", "window_label", "pseudotime", "pt_window"]].copy()
    pt_df.to_csv(OUT_DIR / "pseudotime_values.csv")

    # Save PT-based phi
    for t in range(n_windows - 1):
        np.save(OUT_DIR / f"phi_pt_T{t}_{t+1}.npy", pt_phi_all[t])
    np.save(OUT_DIR / "phi_pt_mean.npy", pt_phi_mean)

    logger.info("  Results saved to %s", OUT_DIR)

    # ── Summary ───────────────────────────────────────────────
    elapsed = time.time() - t0
    logger.info("")
    logger.info("=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)
    logger.info("  Pseudotime vs clinical grade: rho=%.4f (p=%.4e)", rho, p_rho)
    logger.info("  Normal vs GBM AUC: %.4f", auc_ng)
    logger.info("  PT phi vs clinical phi (mean): rho=%.4f", rho_phi)
    if n_ribo > 0:
        logger.info("  Ribosomal gene percentile (PT phi): %.1f", ribo_pctile)
        logger.info("  Ribosomal gene percentile (clinical phi): %.1f", ribo_pctile_clin)
    logger.info("")
    logger.info("Step 8 COMPLETE (%.1f s)", elapsed)


if __name__ == "__main__":
    main()
