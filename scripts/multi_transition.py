"""
scRNAseq Hodge Pipeline — Multi-Transition phi Integration
============================================================
Integrates phi estimates across all available pseudotime transitions,
instead of using only a single primary transition.

Three analyses:
  A. Per-transition phi profiles (n_genes x n_transitions matrix)
  B. Weighted multi-transition integration (4 weighting schemes)
  C. Concordance assessment (Spearman between transitions)
"""
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import spearmanr

from . import config
from .gene_hodge import hodge_gradient_kn, _window_mean_log, build_gene_flow

logger = logging.getLogger(__name__)


def run_single_transition(
    donor_log_corr: Dict[str, np.ndarray],
    donor_window: Dict[str, int],
    donor_ids: List[str],
    w_from: int,
    n_genes: int,
    flow_mode: str = "edge_weight",
) -> Optional[Dict[str, Any]]:
    """Run Hodge on a single transition w_from -> w_from+1."""
    mean_w = _window_mean_log(donor_log_corr, donor_window, donor_ids, w_from)
    mean_w1 = _window_mean_log(donor_log_corr, donor_window, donor_ids, w_from + 1)

    if mean_w is None or mean_w1 is None:
        return None

    delta = mean_w1 - mean_w
    flow = build_gene_flow(delta, mode=flow_mode)
    phi = hodge_gradient_kn(flow, n_genes)

    # Gradient fraction
    n_edges = n_genes * (n_genes - 1) // 2
    grad = np.empty(n_edges, dtype=np.float64)
    e = 0
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            grad[e] = phi[j] - phi[i]
            e += 1
    f_sq = np.dot(flow, flow)
    g_sq = np.dot(grad, grad)
    gf = g_sq / f_sq if f_sq > 1e-30 else 0.0

    return {"phi": phi, "gf": gf, "delta": delta}


def run_multi_transition(
    donor_log_corr: Dict[str, np.ndarray],
    donor_window: Dict[str, int],
    gene_names: List[str],
    n_windows: int = None,
    flow_mode: str = None,
    output_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run multi-transition phi integration.

    Parameters
    ----------
    donor_log_corr : dict
        {donor_id: log(corr) matrix}.
    donor_window : dict
        {donor_id: window_index}.
    gene_names : list of str
    n_windows : int
    flow_mode : str
    output_dir : Path, optional

    Returns
    -------
    dict with per-transition phi, integrated phi (4 schemes), concordance.
    """
    if n_windows is None:
        n_windows = config.N_WINDOWS
    if flow_mode is None:
        flow_mode = config.GENE_HODGE_FLOW_MODE

    N = len(gene_names)
    donor_ids = list(donor_log_corr.keys())
    logger.info("Multi-transition: N=%d genes, %d windows, %d possible transitions",
                N, n_windows, n_windows - 1)

    # ── A: Per-transition phi ─────────────────────────────────
    transitions = []
    phi_list = []
    gf_list = []

    for w in range(n_windows - 1):
        result = run_single_transition(
            donor_log_corr, donor_window, donor_ids, w, N, flow_mode,
        )
        if result is None:
            logger.warning("  w%d->w%d: insufficient donors, skipping", w, w + 1)
            continue

        transitions.append((w, w + 1))
        phi_list.append(result["phi"])
        gf_list.append(result["gf"])
        logger.info("  w%d->w%d: GF=%.4f", w, w + 1, result["gf"])

    if len(phi_list) == 0:
        logger.error("No valid transitions")
        return {"status": "FAILED", "reason": "no valid transitions"}

    phi_matrix = np.column_stack(phi_list)  # (N, n_transitions)
    n_trans = len(transitions)

    # ── B: Weighted integration (4 schemes) ───────────────────
    gf_arr = np.array(gf_list)

    # 1. Uniform: simple mean
    phi_uniform = phi_matrix.mean(axis=1)

    # 2. GF-weighted: weight by gradient fraction
    gf_weights = gf_arr / gf_arr.sum() if gf_arr.sum() > 0 else np.ones(n_trans) / n_trans
    phi_gf_weighted = phi_matrix @ gf_weights

    # 3. Rank-based: median rank across transitions
    from scipy.stats import rankdata
    rank_matrix = np.apply_along_axis(lambda x: rankdata(-x), 0, phi_matrix)
    phi_rank_median = -np.median(rank_matrix, axis=1)  # negative so higher = more upstream

    # 4. Max-GF: use only the transition with highest GF
    best_idx = np.argmax(gf_arr)
    phi_best = phi_matrix[:, best_idx]

    integrated = {
        "uniform": phi_uniform,
        "gf_weighted": phi_gf_weighted,
        "rank_median": phi_rank_median,
        "best_transition": phi_best,
    }

    # ── C: Concordance ────────────────────────────────────────
    concordance = np.ones((n_trans, n_trans), dtype=np.float64)
    for i in range(n_trans):
        for j in range(i + 1, n_trans):
            rho, _ = spearmanr(phi_list[i], phi_list[j])
            concordance[i, j] = rho
            concordance[j, i] = rho

    mean_concordance = float(concordance[np.triu_indices(n_trans, k=1)].mean()) if n_trans > 1 else 1.0

    logger.info("  Concordance: mean Spearman=%.3f across %d transition pairs",
                mean_concordance, n_trans * (n_trans - 1) // 2)

    # ── Assemble result ───────────────────────────────────────
    result = {
        "status": "OK",
        "n_genes": N,
        "n_transitions": n_trans,
        "transitions": transitions,
        "gf_per_transition": gf_list,
        "best_transition": transitions[best_idx],
        "best_gf": float(gf_arr[best_idx]),
        "mean_concordance": mean_concordance,
        "phi_matrix": phi_matrix,
        "integrated": integrated,
        "concordance_matrix": concordance,
    }

    # ── Save ──────────────────────────────────────────────────
    if output_dir is not None:
        import pandas as pd
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Gene-level results
        df = pd.DataFrame({"gene": gene_names})
        for i, (w_from, w_to) in enumerate(transitions):
            df[f"phi_w{w_from}_{w_to}"] = phi_list[i]
        df["phi_uniform"] = phi_uniform
        df["phi_gf_weighted"] = phi_gf_weighted
        df["phi_rank_median"] = phi_rank_median
        df = df.sort_values("phi_gf_weighted", ascending=False)
        df.to_csv(output_dir / "multi_transition_phi.csv", index=False)

        # Summary
        summary = {
            "n_genes": N,
            "n_transitions": n_trans,
            "transitions": transitions,
            "gf_per_transition": gf_list,
            "best_transition": transitions[best_idx],
            "best_gf": float(gf_arr[best_idx]),
            "mean_concordance": mean_concordance,
        }
        with open(output_dir / "multi_transition_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)

        logger.info("Multi-transition results saved to %s", output_dir)

    return result
