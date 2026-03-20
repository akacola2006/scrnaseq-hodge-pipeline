"""
scRNAseq Hodge Pipeline — Gene-Level Hodge Decomposition
==========================================================
Applies discrete Hodge decomposition to gene-level correlation dynamics
within the upstream celltype. Uses the analytical pseudoinverse of the
K_N graph Laplacian for efficiency with large gene sets.
"""
import json
import logging
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from . import config
from .spd import spd_log

logger = logging.getLogger(__name__)


# =====================================================================
# 1. Execution Gate
# =====================================================================

def check_execution_gate(
    bootstrap_summary: Dict[str, Any],
    upstream_celltype: str,
) -> Tuple[bool, str]:
    """Check whether gene Hodge should run for upstream_celltype."""
    top3_repro = bootstrap_summary.get("top3_reproducibility", {})
    phi_ci = bootstrap_summary.get("phi_ci", {})

    repro = top3_repro.get(upstream_celltype, 0.0)
    if repro < config.BOOTSTRAP_STABLE_THRESHOLD:
        return False, f"Bootstrap repro {repro:.1%} < {config.BOOTSTRAP_STABLE_THRESHOLD:.0%}"

    ci = phi_ci.get(upstream_celltype)
    if ci is None:
        return False, f"No phi CI for {upstream_celltype}"

    ci_lower, ci_upper = float(ci[0]), float(ci[1])
    if ci_lower <= 0.0:
        return False, f"phi CI crosses zero: [{ci_lower:.4f}, {ci_upper:.4f}]"

    return True, f"Gate PASSED: boot={repro:.1%}, phi CI=[{ci_lower:.4f}, {ci_upper:.4f}]"


# =====================================================================
# 2. Gene Set Resolution
# =====================================================================

def resolve_gene_indices(
    gene_set: List[str],
    all_gene_names: List[str],
) -> Tuple[np.ndarray, List[str]]:
    """Map gene set to column indices in the residual matrix."""
    name_to_idx = {g: i for i, g in enumerate(all_gene_names)}
    indices, resolved, n_miss = [], [], 0
    for g in gene_set:
        if g in name_to_idx:
            indices.append(name_to_idx[g])
            resolved.append(g)
        else:
            n_miss += 1
    if n_miss:
        logger.warning("%d / %d genes not found (dropped)", n_miss, len(gene_set))
    return np.array(indices, dtype=np.int64), resolved


# =====================================================================
# 3. Pre-compute Per-Donor log(Correlation) Matrices
# =====================================================================

def precompute_donor_log_corr(
    residuals: np.ndarray,
    donor_ids: np.ndarray,
    gene_indices: np.ndarray,
    pt_df: "pd.DataFrame",
) -> Tuple[Dict[str, np.ndarray], Dict[str, int]]:
    """Compute Ledoit-Wolf gene correlation -> log(corr) per donor."""
    from sklearn.covariance import LedoitWolf

    donor_to_window = dict(zip(
        pt_df["donor_id"].astype(str), pt_df["window"].astype(int),
    ))

    donor_log_corr = {}
    donor_window = {}
    unique_donors = sorted(set(donor_ids))

    for did in unique_donors:
        w = donor_to_window.get(did)
        if w is None:
            continue

        mask = donor_ids == did
        donor_resid = residuals[mask][:, gene_indices]
        n_cells = donor_resid.shape[0]
        n_genes_sel = len(gene_indices)

        if n_cells < max(3, n_genes_sel // 10):
            continue

        try:
            lw = LedoitWolf()
            lw.fit(donor_resid)
            cov = lw.covariance_

            std = np.sqrt(np.diag(cov))
            std = np.maximum(std, 1e-15)
            corr = cov / np.outer(std, std)
            np.fill_diagonal(corr, 1.0)

            eigvals, eigvecs = np.linalg.eigh(corr)
            eigvals = np.maximum(eigvals, 1e-10)
            corr = (eigvecs * eigvals) @ eigvecs.T
            corr = (corr + corr.T) / 2.0

            donor_log_corr[did] = spd_log(corr)
            donor_window[did] = w
        except Exception as exc:
            logger.warning("LW failed donor %s: %s", did, exc)

    logger.info("Pre-computed log(corr) for %d / %d donors",
                len(donor_log_corr), len(unique_donors))
    return donor_log_corr, donor_window


# =====================================================================
# 4. Core: Gene Hodge phi via K_N analytical pseudoinverse
# =====================================================================

def _window_mean_log(donor_log_corr, donor_window, donor_ids, target_window):
    """Log-Euclidean mean for a specific window."""
    logs = [
        donor_log_corr[d]
        for d in donor_ids
        if d in donor_log_corr and donor_window.get(d) == target_window
    ]
    if len(logs) < config.MIN_DONORS_PER_WINDOW:
        return None
    return np.mean(logs, axis=0)


def hodge_gradient_kn(flow: np.ndarray, n_genes: int) -> np.ndarray:
    """Analytical Hodge gradient on K_N.

    For the complete graph K_N:
        phi_i = (div_i - mean(div)) / N
    where div_i = sum_j f(i,j).
    """
    N = n_genes
    div = np.zeros(N, dtype=np.float64)

    e_idx = 0
    for i in range(N):
        for j in range(i + 1, N):
            div[i] += flow[e_idx]
            div[j] -= flow[e_idx]
            e_idx += 1

    phi = (div - div.mean()) / N
    return phi


def build_gene_flow(
    delta: np.ndarray,
    mode: str = "sign",
) -> np.ndarray:
    """Build gene-level edge flow from correlation change matrix Delta.

    Parameters
    ----------
    delta : np.ndarray, shape (N, N)
        Change in log-correlation between windows.
    mode : str
        Flow construction mode: "sign", "weighted", "edge_product", "edge_weight"
    """
    N = delta.shape[0]
    n_edges = N * (N - 1) // 2

    # Per-gene change magnitude
    d = np.linalg.norm(delta, axis=1)

    flow = np.empty(n_edges, dtype=np.float64)
    e_idx = 0
    for i in range(N):
        for j in range(i + 1, N):
            if mode == "sign":
                flow[e_idx] = np.sign(d[i] - d[j])
            elif mode == "weighted":
                flow[e_idx] = d[i] - d[j]
            elif mode == "edge_product":
                flow[e_idx] = delta[i, j] * (d[i] - d[j])
            elif mode == "edge_weight":
                flow[e_idx] = abs(delta[i, j]) * np.sign(d[i] - d[j])
            else:
                raise ValueError(f"Unknown flow mode: {mode}")
            e_idx += 1

    return flow


def run_gene_hodge(
    donor_log_corr: Dict[str, np.ndarray],
    donor_window: Dict[str, int],
    donor_ids: List[str],
    w_star: int,
    gene_names: List[str],
    flow_mode: str = None,
    n_perm: int = None,
) -> Dict[str, Any]:
    """Run gene-level Hodge decomposition for transition w* -> w*+1.

    Returns dict with: phi, gene_names, gradient_fraction, p_value, classification.
    """
    if flow_mode is None:
        flow_mode = config.GENE_HODGE_FLOW_MODE
    if n_perm is None:
        n_perm = config.GENE_HODGE_PERMUTATION_N

    N = len(gene_names)
    if N < config.GENE_HODGE_MIN_GENES:
        return {"status": "SKIPPED", "reason": f"Too few genes: {N}"}

    # Window means
    mean_w = _window_mean_log(donor_log_corr, donor_window, donor_ids, w_star)
    mean_w1 = _window_mean_log(donor_log_corr, donor_window, donor_ids, w_star + 1)

    if mean_w is None or mean_w1 is None:
        return {"status": "SKIPPED", "reason": "Insufficient donors in windows"}

    delta = mean_w1 - mean_w

    # Build flow and compute Hodge
    flow = build_gene_flow(delta, mode=flow_mode)
    phi = hodge_gradient_kn(flow, N)

    # Gradient fraction
    n_edges = N * (N - 1) // 2
    grad_flow = np.empty(n_edges, dtype=np.float64)
    e_idx = 0
    for i in range(N):
        for j in range(i + 1, N):
            grad_flow[e_idx] = phi[j] - phi[i]
            e_idx += 1

    f_sq = np.dot(flow, flow)
    g_sq = np.dot(grad_flow, grad_flow)
    gf = g_sq / f_sq if f_sq > 1e-30 else 0.0

    # Permutation test for GF
    rng = np.random.RandomState(42)
    null_gf = np.empty(n_perm, dtype=np.float64)
    for p in range(n_perm):
        perm_flow = flow.copy()
        rng.shuffle(perm_flow)
        signs = rng.choice([-1.0, 1.0], size=len(perm_flow))
        perm_flow *= signs

        perm_phi = hodge_gradient_kn(perm_flow, N)
        perm_grad = np.empty(n_edges, dtype=np.float64)
        e_idx = 0
        for i in range(N):
            for j in range(i + 1, N):
                perm_grad[e_idx] = perm_phi[j] - perm_phi[i]
                e_idx += 1

        perm_g_sq = np.dot(perm_grad, perm_grad)
        perm_f_sq = np.dot(perm_flow, perm_flow)
        null_gf[p] = perm_g_sq / perm_f_sq if perm_f_sq > 1e-30 else 0.0

    # Phipson & Smyth (2010) pseudocount correction
    p_value = (np.sum(null_gf >= gf) + 1) / (n_perm + 1)

    # Classify genes
    phi_mean = phi.mean()
    phi_std = phi.std()
    classification = []
    for i in range(N):
        if phi[i] > phi_mean + config.GENE_HODGE_HIGH_SIGMA * phi_std:
            classification.append("High")
        elif phi[i] > phi_mean + config.GENE_HODGE_MEDIUM_SIGMA * phi_std:
            classification.append("Medium")
        else:
            classification.append("Low")

    logger.info("Gene Hodge: N=%d, GF=%.4f, p=%.4f, High=%d, Med=%d, Low=%d",
                N, gf, p_value,
                classification.count("High"),
                classification.count("Medium"),
                classification.count("Low"))

    return {
        "status": "OK",
        "phi": phi,
        "gene_names": gene_names,
        "gradient_fraction": gf,
        "p_value": p_value,
        "classification": classification,
        "n_genes": N,
        "flow_mode": flow_mode,
        "transition": (w_star, w_star + 1),
    }
