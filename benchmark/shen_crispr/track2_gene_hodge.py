"""
P4 Track 2: Gene-Level Hodge Decomposition on Perturbation Data
=================================================================

Applies the IDS-SALS Hodge decomposition framework to CRISPR perturbation
data (GSE274058) to test whether the *method* captures perturbation-induced
gene hierarchy.

Design:
  - For each CRISPR target, define two conditions:
      Control (barcode-negative, w=0)  →  KO (target knockdown, w=1)
  - Estimate gene-gene correlation matrices per condition (pooled Ledoit-Wolf)
  - Compute delta in log-correlation space → gene-level d_corr → flow → Hodge
  - Test whether the KO target gene has high φ (upstream position)

Key difference from IDS-SALS:
  - IDS-SALS: per-donor correlation → Log-Euclidean mean per pseudotime window
  - Track 2:  pooled correlation per condition (no natural donor structure)
  - Bootstrap: cell resampling (not donor resampling)

The Hodge decomposition core (hodge_gradient_kn, flow construction) is
replicated here from the frozen gene_hodge.py to keep Track 2 self-contained.

Reference: gene_hodge.py from sals_analysis_frozen_20260211/scripts/
"""
import logging
from itertools import combinations
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from sklearn.covariance import LedoitWolf

logger = logging.getLogger(__name__)


# =====================================================================
# 1. Correlation Estimation (Pooled Ledoit-Wolf)
# =====================================================================

def estimate_correlation(
    expr: np.ndarray,
    eps: float = 1e-12,
) -> np.ndarray:
    """Estimate gene-gene correlation matrix using Ledoit-Wolf shrinkage.

    Parameters
    ----------
    expr : np.ndarray, shape (n_cells, n_genes)
        Gene expression matrix (log-normalised).
    eps : float
        Floor for diagonal entries.

    Returns
    -------
    corr : np.ndarray, shape (n_genes, n_genes)
        SPD correlation matrix.
    """
    lw = LedoitWolf()
    lw.fit(expr)
    cov = lw.covariance_

    # Cov → Corr
    std = np.sqrt(np.maximum(np.diag(cov), eps))
    corr = cov / np.outer(std, std)
    corr = (corr + corr.T) / 2.0
    np.fill_diagonal(corr, 1.0)

    # Ensure SPD (eigenvalue floor)
    eigvals, eigvecs = np.linalg.eigh(corr)
    eigvals = np.maximum(eigvals, 1e-10)
    corr = (eigvecs * eigvals) @ eigvecs.T
    corr = (corr + corr.T) / 2.0
    np.fill_diagonal(corr, 1.0)

    return corr


def spd_log(A: np.ndarray) -> np.ndarray:
    """Matrix logarithm for SPD matrix via eigendecomposition.

    Replicates spd.spd_log from frozen analysis.
    """
    A_sym = (A + A.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(A_sym)
    eigvals = np.maximum(eigvals, 1e-15)
    log_A = eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return (log_A + log_A.T) / 2.0


# =====================================================================
# 2. Flow Construction (from gene_hodge.py)
# =====================================================================

def build_flow_sign(gene_d_corr: np.ndarray, n_genes: int) -> np.ndarray:
    """Build pairwise flow: flow(i,j) = sign(d_corr_i - d_corr_j)."""
    n_edges = n_genes * (n_genes - 1) // 2
    flow = np.empty(n_edges, dtype=np.float64)
    idx = 0
    for i in range(n_genes):
        n_others = n_genes - i - 1
        if n_others == 0:
            break
        diff = gene_d_corr[i] - gene_d_corr[i + 1: n_genes]
        flow[idx: idx + n_others] = np.sign(diff)
        idx += n_others
    return flow


def build_flow_edge_weight(
    delta: np.ndarray,
    gene_d_corr: np.ndarray,
    n_genes: int,
) -> np.ndarray:
    """Build flow: flow(i,j) = |delta[i,j]| * sign(d_corr_i - d_corr_j).

    This is the "edge_weight" mode used in the IDS-SALS Oligo analysis.
    """
    n_edges = n_genes * (n_genes - 1) // 2
    flow = np.empty(n_edges, dtype=np.float64)
    idx = 0
    for i in range(n_genes):
        n_others = n_genes - i - 1
        if n_others == 0:
            break
        js = np.arange(i + 1, n_genes)
        delta_ij = delta[i, js]
        diff = gene_d_corr[i] - gene_d_corr[js]
        flow[idx: idx + n_others] = np.abs(delta_ij) * np.sign(diff)
        idx += n_others
    return flow


def standardize_flow(flow: np.ndarray) -> Tuple[np.ndarray, float]:
    """Standardize flow by dividing by std(flow)."""
    flow_std = float(np.std(flow))
    if flow_std < 1e-15:
        logger.warning("Flow std near zero (%.2e); skipping standardization", flow_std)
        return flow, flow_std
    return flow / flow_std, flow_std


# =====================================================================
# 3. Hodge Decomposition (K_N analytical)
# =====================================================================

def hodge_gradient_kn(
    flow: np.ndarray,
    n_nodes: int,
) -> Dict[str, Any]:
    """Gradient-only Hodge decomposition for complete graph K_n.

    L = nI - J,  L+ = (1/n)(I - J/n),  phi = L+ @ div.

    Replicates gene_hodge.hodge_gradient_kn from frozen analysis.
    """
    edges = np.array(list(combinations(range(n_nodes), 2)), dtype=np.int32)
    div = np.zeros(n_nodes, dtype=np.float64)
    np.add.at(div, edges[:, 0], -flow)
    np.add.at(div, edges[:, 1], +flow)

    phi = (div - div.mean()) / n_nodes
    gradient = phi[edges[:, 1]] - phi[edges[:, 0]]

    flow_sq = np.dot(flow, flow)
    grad_sq = np.dot(gradient, gradient)

    gradient_fraction = float(grad_sq / flow_sq) if flow_sq > 1e-30 else 0.0
    non_gradient_fraction = max(0.0, 1.0 - gradient_fraction)

    if flow_sq < 1e-30 or grad_sq < 1e-30:
        r_phi = 0.0
    else:
        r_phi = float(np.corrcoef(flow, gradient)[0, 1])
        if np.isnan(r_phi):
            r_phi = 0.0

    return {
        "phi": phi,
        "gradient_fraction": gradient_fraction,
        "non_gradient_fraction": non_gradient_fraction,
        "r_phi": r_phi,
    }


# =====================================================================
# 4. Permutation Test
# =====================================================================

def permutation_test(
    flow: np.ndarray,
    n_nodes: int,
    n_perm: int = 1000,
    seed: int = 42,
) -> Dict[str, Any]:
    """Null-distribution for gradient_fraction via flow shuffle + sign-flip."""
    rng = np.random.RandomState(seed)
    obs_gf = hodge_gradient_kn(flow, n_nodes)["gradient_fraction"]

    null_gf = np.empty(n_perm, dtype=np.float64)
    for p in range(n_perm):
        perm_flow = flow.copy()
        rng.shuffle(perm_flow)
        perm_flow *= rng.choice([-1.0, 1.0], size=len(perm_flow))
        null_gf[p] = hodge_gradient_kn(perm_flow, n_nodes)["gradient_fraction"]

    p_value = float((np.sum(null_gf >= obs_gf) + 1) / (n_perm + 1))
    return {
        "p_value": p_value,
        "observed_gradient_fraction": obs_gf,
        "null_mean": float(null_gf.mean()),
        "null_std": float(null_gf.std()),
    }


# =====================================================================
# 5. Gene Classification
# =====================================================================

def classify_genes(
    phi: np.ndarray,
    gene_names: List[str],
    high_sigma: float = 1.5,
) -> Dict[str, Any]:
    """Classify genes into High / Medium / Low by phi.

    High: phi > mean + high_sigma*std
    Medium: phi > mean
    Low: rest
    """
    mean_phi = float(phi.mean())
    std_phi = float(phi.std())
    high_thresh = mean_phi + high_sigma * std_phi

    high, medium, low = [], [], []
    for i, name in enumerate(gene_names):
        entry = {"gene": name, "phi": float(phi[i]), "rank": 0}
        if phi[i] > high_thresh:
            high.append(entry)
        elif phi[i] > mean_phi:
            medium.append(entry)
        else:
            low.append(entry)

    for lst in (high, medium, low):
        lst.sort(key=lambda x: -x["phi"])

    # Assign ranks
    all_sorted = sorted(enumerate(gene_names), key=lambda x: -phi[x[0]])
    rank_map = {name: r + 1 for r, (_, name) in enumerate(all_sorted)}
    for lst in (high, medium, low):
        for entry in lst:
            entry["rank"] = rank_map[entry["gene"]]

    return {
        "high": high, "medium": medium, "low": low,
        "thresholds": {"high": high_thresh, "medium": mean_phi},
        "stats": {
            "mean_phi": mean_phi, "std_phi": std_phi,
            "n_high": len(high), "n_medium": len(medium), "n_low": len(low),
        },
    }


# =====================================================================
# 6. Core: Single-Perturbation Hodge Analysis
# =====================================================================

def compute_perturbation_hodge(
    expr_ko: np.ndarray,
    expr_ctrl: np.ndarray,
    gene_names: List[str],
    flow_mode: str = "edge_weight",
) -> Optional[Dict[str, Any]]:
    """Run Hodge decomposition for one perturbation: Control → KO.

    Parameters
    ----------
    expr_ko : np.ndarray, shape (n_ko, n_genes)
        Expression matrix for KO cells (log-normalised, gene-subset).
    expr_ctrl : np.ndarray, shape (n_ctrl, n_genes)
        Expression matrix for Control cells.
    gene_names : list of str
        Gene names (length n_genes).
    flow_mode : str
        "sign" or "edge_weight" (matches IDS-SALS Oligo analysis).

    Returns
    -------
    dict with phi, gradient_fraction, classification, etc.
    None if estimation fails.
    """
    n_genes = len(gene_names)
    if n_genes < 10:
        logger.warning("Too few genes (%d) for Hodge analysis", n_genes)
        return None

    # 1. Estimate correlation matrices
    try:
        corr_ctrl = estimate_correlation(expr_ctrl)
        corr_ko = estimate_correlation(expr_ko)
    except Exception as exc:
        logger.warning("Correlation estimation failed: %s", exc)
        return None

    # 2. Log-correlation space delta
    log_ctrl = spd_log(corr_ctrl)
    log_ko = spd_log(corr_ko)
    delta = log_ko - log_ctrl  # (n_genes, n_genes)

    # 3. Per-gene d_corr (row-wise L2 norm of delta)
    gene_d_corr = np.linalg.norm(delta, axis=1)  # (n_genes,)

    # 4. Build flow
    if flow_mode == "edge_weight":
        flow = build_flow_edge_weight(delta, gene_d_corr, n_genes)
        flow, flow_std_raw = standardize_flow(flow)
    elif flow_mode == "sign":
        flow = build_flow_sign(gene_d_corr, n_genes)
        flow_std_raw = None
    else:
        raise ValueError(f"Unknown flow_mode: {flow_mode}")

    # 5. Hodge decomposition
    hodge = hodge_gradient_kn(flow, n_genes)

    # 6. Classification
    classification = classify_genes(hodge["phi"], gene_names)

    return {
        "phi": hodge["phi"],
        "flow": flow,
        "gene_d_corr": gene_d_corr,
        "gradient_fraction": hodge["gradient_fraction"],
        "non_gradient_fraction": hodge["non_gradient_fraction"],
        "r_phi": hodge["r_phi"],
        "flow_mode": flow_mode,
        "flow_std_raw": flow_std_raw,
        "classification": classification,
        "n_ko": expr_ko.shape[0],
        "n_ctrl": expr_ctrl.shape[0],
        "n_genes": n_genes,
    }


# =====================================================================
# 7. Cell-Level Bootstrap
# =====================================================================

def bootstrap_perturbation_hodge(
    expr_ko: np.ndarray,
    expr_ctrl: np.ndarray,
    gene_names: List[str],
    n_bootstrap: int = 100,
    flow_mode: str = "edge_weight",
    seed: int = 42,
) -> Dict[str, Any]:
    """Bootstrap Hodge via cell resampling within each condition.

    For each iteration:
      1. Resample KO cells with replacement
      2. Resample Control cells with replacement
      3. Re-estimate correlations and run Hodge
      4. Collect per-gene phi

    Returns CI, reproducibility, gradient_fraction CI.
    """
    rng = np.random.RandomState(seed)
    n_ko = expr_ko.shape[0]
    n_ctrl = expr_ctrl.shape[0]
    n_genes = len(gene_names)

    phi_matrix = []
    gf_list = []
    high_counts = np.zeros(n_genes, dtype=np.int64)
    n_success = 0
    n_failed = 0

    HIGH_SIGMA = 1.5  # Match IDS-SALS config

    for b in range(n_bootstrap):
        # Resample cells
        ko_idx = rng.randint(0, n_ko, size=n_ko)
        ctrl_idx = rng.randint(0, n_ctrl, size=n_ctrl)

        result = compute_perturbation_hodge(
            expr_ko[ko_idx],
            expr_ctrl[ctrl_idx],
            gene_names,
            flow_mode=flow_mode,
        )

        if result is None:
            n_failed += 1
            continue

        phi_b = result["phi"]
        phi_matrix.append(phi_b)
        gf_list.append(result["gradient_fraction"])

        # Classification for this iteration
        mean_b = phi_b.mean()
        std_b = phi_b.std()
        thresh_b = mean_b + HIGH_SIGMA * std_b if std_b > 1e-15 else mean_b
        high_counts += (phi_b > thresh_b).astype(np.int64)
        n_success += 1

        if (b + 1) % 20 == 0:
            logger.info("  Bootstrap %d / %d (success=%d)", b + 1, n_bootstrap, n_success)

    logger.info(
        "Bootstrap complete: %d / %d success, %d failed",
        n_success, n_bootstrap, n_failed,
    )

    if n_success < 10:
        return {
            "n_bootstrap": n_bootstrap,
            "n_success": n_success,
            "n_failed": n_failed,
            "insufficient": True,
        }

    phi_matrix = np.array(phi_matrix)  # (n_success, n_genes)
    gf_arr = np.array(gf_list)

    # 95% CI per gene
    phi_ci = {}
    for i, g in enumerate(gene_names):
        lo = float(np.percentile(phi_matrix[:, i], 2.5))
        hi = float(np.percentile(phi_matrix[:, i], 97.5))
        phi_ci[g] = [lo, hi]

    # High-tier reproducibility per gene
    high_repro = {}
    for i, g in enumerate(gene_names):
        high_repro[g] = float(high_counts[i] / n_success)

    n_stable_high = sum(1 for v in high_repro.values() if v >= 0.80)

    # Gradient fraction CI
    gf_ci = [float(np.percentile(gf_arr, 2.5)), float(np.percentile(gf_arr, 97.5))]

    # Mean phi per gene (point estimate from bootstrap)
    phi_mean = {}
    for i, g in enumerate(gene_names):
        phi_mean[g] = float(phi_matrix[:, i].mean())

    return {
        "n_bootstrap": n_bootstrap,
        "n_success": n_success,
        "n_failed": n_failed,
        "insufficient": False,
        "phi_ci": phi_ci,
        "phi_mean": phi_mean,
        "high_reproducibility": high_repro,
        "n_stable_high": n_stable_high,
        "gradient_fraction_ci": gf_ci,
        "gradient_fraction_mean": float(gf_arr.mean()),
    }


# =====================================================================
# 8. Target Gene Evaluation
# =====================================================================

def evaluate_target_gene(
    phi: np.ndarray,
    gene_names: List[str],
    target_gene: str,
    classification: Dict[str, Any],
    bootstrap_result: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Evaluate the KO target gene's position in the Hodge hierarchy.

    The key Track 2 question: does the target gene rank as "upstream"
    (high φ) when Hodge decomposition is applied to KO-vs-Control
    transcriptomic variation?

    Parameters
    ----------
    phi : np.ndarray
        Per-gene phi from the point estimate.
    gene_names : list of str
        Gene names (length n_genes).
    target_gene : str
        The CRISPR target gene name (in the gene_names space).
    classification : dict
        Gene classification from classify_genes().
    bootstrap_result : dict, optional
        Bootstrap results for CI and reproducibility.

    Returns
    -------
    dict with target gene evaluation metrics.
    """
    n_genes = len(gene_names)
    name_to_idx = {g: i for i, g in enumerate(gene_names)}

    if target_gene not in name_to_idx:
        return {
            "target_gene": target_gene,
            "in_gene_set": False,
            "reason": "Target gene not in HVG set",
        }

    idx = name_to_idx[target_gene]
    target_phi = float(phi[idx])

    # Rank (1 = highest phi)
    sorted_indices = np.argsort(-phi)
    rank = int(np.where(sorted_indices == idx)[0][0]) + 1
    percentile = 100.0 * (1.0 - rank / n_genes)

    # Z-score
    mean_phi = phi.mean()
    std_phi = phi.std()
    z_score = float((target_phi - mean_phi) / std_phi) if std_phi > 1e-15 else 0.0

    # Classification tier
    tier = "Low"
    for entry in classification["high"]:
        if entry["gene"] == target_gene:
            tier = "High"
            break
    if tier == "Low":
        for entry in classification["medium"]:
            if entry["gene"] == target_gene:
                tier = "Medium"
                break

    result = {
        "target_gene": target_gene,
        "in_gene_set": True,
        "phi": target_phi,
        "rank": rank,
        "n_genes": n_genes,
        "percentile": percentile,
        "z_score": z_score,
        "tier": tier,
    }

    # Bootstrap CI and reproducibility
    if bootstrap_result and not bootstrap_result.get("insufficient", True):
        phi_ci = bootstrap_result.get("phi_ci", {}).get(target_gene)
        high_repro = bootstrap_result.get("high_reproducibility", {}).get(target_gene)
        phi_mean_boot = bootstrap_result.get("phi_mean", {}).get(target_gene)

        if phi_ci is not None:
            result["phi_ci"] = phi_ci
            result["phi_ci_above_zero"] = phi_ci[0] > 0
        if high_repro is not None:
            result["high_reproducibility"] = high_repro
            result["bootstrap_stable_high"] = high_repro >= 0.80
        if phi_mean_boot is not None:
            result["phi_mean_bootstrap"] = phi_mean_boot

    return result


# =====================================================================
# 9. Aggregate Cross-Perturbation Analysis
# =====================================================================

def aggregate_track2_results(
    per_perturbation: Dict[str, Dict[str, Any]],
    target_evaluations: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """Aggregate Track 2 results across all perturbations.

    Tests:
    1. Mean target gene percentile > 50 (one-sample Wilcoxon signed-rank)
    2. Target gene tier distribution (how many High/Medium/Low)
    3. ALS targets vs Control targets comparison
    4. Gradient fraction across perturbations
    """
    from scipy import stats

    # Collect target gene percentiles
    percentiles = {}
    categories = {}
    for gene, eval_res in target_evaluations.items():
        if eval_res.get("in_gene_set", False):
            percentiles[gene] = eval_res["percentile"]
            # Find category from per_perturbation
            for pname, pres in per_perturbation.items():
                if pres.get("target_gene_mouse") == gene or pres.get("target_gene") == gene:
                    categories[gene] = pres.get("category", "Unknown")
                    break

    if not percentiles:
        return {"error": "No target genes found in gene sets"}

    pct_values = np.array(list(percentiles.values()))
    n_targets = len(pct_values)

    # One-sample Wilcoxon signed-rank: percentile > 50?
    if n_targets >= 5:
        stat, p_wilcoxon = stats.wilcoxon(pct_values - 50.0, alternative="greater")
    else:
        stat, p_wilcoxon = np.nan, np.nan

    # Binomial test: fraction above median > 0.5?
    n_above_median = int(np.sum(pct_values > 50.0))
    if n_targets >= 1:
        p_binom = float(stats.binomtest(n_above_median, n_targets, 0.5, alternative="greater").pvalue)
    else:
        p_binom = np.nan

    # Tier distribution
    tiers = {"High": 0, "Medium": 0, "Low": 0, "Not_in_set": 0}
    for gene, eval_res in target_evaluations.items():
        if eval_res.get("in_gene_set", False):
            tiers[eval_res.get("tier", "Low")] += 1
        else:
            tiers["Not_in_set"] += 1

    # Category-level analysis
    category_stats = {}
    for cat in set(categories.values()):
        cat_pcts = [pct for g, pct in percentiles.items() if categories.get(g) == cat]
        if cat_pcts:
            category_stats[cat] = {
                "n": len(cat_pcts),
                "mean_percentile": float(np.mean(cat_pcts)),
                "std_percentile": float(np.std(cat_pcts)) if len(cat_pcts) > 1 else 0.0,
                "percentiles": {g: pct for g, pct in percentiles.items() if categories.get(g) == cat},
            }

    # ALS vs Control comparison
    als_pcts = [pct for g, pct in percentiles.items() if categories.get(g) == "ALS"]
    ctrl_pcts = [pct for g, pct in percentiles.items() if categories.get(g) in ("Control", "Ctrl(neg)")]
    if len(als_pcts) >= 2 and len(ctrl_pcts) >= 2:
        mw_stat, p_als_vs_ctrl = stats.mannwhitneyu(
            als_pcts, ctrl_pcts, alternative="greater"
        )
    else:
        mw_stat, p_als_vs_ctrl = np.nan, np.nan

    # Gradient fraction summary
    gf_values = []
    for pres in per_perturbation.values():
        gf = pres.get("gradient_fraction")
        if gf is not None:
            gf_values.append(gf)

    gf_summary = {}
    if gf_values:
        gf_arr = np.array(gf_values)
        gf_summary = {
            "mean": float(gf_arr.mean()),
            "std": float(gf_arr.std()),
            "min": float(gf_arr.min()),
            "max": float(gf_arr.max()),
        }

    return {
        "n_targets_evaluated": n_targets,
        "n_targets_in_gene_set": n_targets,
        "mean_percentile": float(pct_values.mean()),
        "std_percentile": float(pct_values.std()) if n_targets > 1 else 0.0,
        "median_percentile": float(np.median(pct_values)),
        "n_above_median": n_above_median,
        "wilcoxon_p": float(p_wilcoxon) if not np.isnan(p_wilcoxon) else None,
        "binom_p": float(p_binom) if not np.isnan(p_binom) else None,
        "tier_distribution": tiers,
        "category_stats": category_stats,
        "als_vs_control_p": float(p_als_vs_ctrl) if not np.isnan(p_als_vs_ctrl) else None,
        "gradient_fraction_summary": gf_summary,
        "per_target_summary": {
            gene: {
                "percentile": pct,
                "category": categories.get(gene, "Unknown"),
                "tier": target_evaluations[gene].get("tier", "N/A"),
                "z_score": target_evaluations[gene].get("z_score", np.nan),
            }
            for gene, pct in percentiles.items()
        },
    }
