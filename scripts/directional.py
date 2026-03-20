"""
scRNAseq Hodge Pipeline — Directional Decomposition (Δ⁺ / Δ⁻ Split)
=====================================================================
Splits the correlation-change matrix Δ into gain (Δ⁺, correlations that
increased) and loss (Δ⁻, correlations that decreased), then runs Hodge
decomposition on each component separately.

This reveals whether upstream genes are driven by correlation gains
(coordinated activation) or losses (decoupling).
"""
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import spearmanr, skew, kurtosis

from . import config
from .gene_hodge import hodge_gradient_kn

logger = logging.getLogger(__name__)


def decompose_delta(delta: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Split symmetric Δ into gain (Δ>0) and loss (|Δ<0|) components."""
    delta_plus = np.maximum(delta, 0.0)
    delta_minus = np.abs(np.minimum(delta, 0.0))
    return delta_plus, delta_minus


def _build_flow_edge_weight(delta_component: np.ndarray) -> np.ndarray:
    """Build edge-weight flow: f(i,j) = |Δ_ij| * sign(d_i - d_j)."""
    n = delta_component.shape[0]
    d_corr = np.linalg.norm(delta_component, axis=1)
    n_edges = n * (n - 1) // 2
    flow = np.empty(n_edges, dtype=np.float64)
    idx = 0
    for i in range(n):
        js = np.arange(i + 1, n)
        if len(js) == 0:
            break
        w_ij = np.abs(delta_component[i, js])
        diff = d_corr[i] - d_corr[js]
        flow[idx:idx + len(js)] = w_ij * np.sign(diff)
        idx += len(js)
    return flow


def _hodge_with_gf(flow: np.ndarray, n: int) -> Dict[str, Any]:
    """Hodge gradient on K_N with gradient fraction."""
    phi = hodge_gradient_kn(flow, n)
    n_edges = n * (n - 1) // 2
    grad = np.empty(n_edges, dtype=np.float64)
    e = 0
    for i in range(n):
        for j in range(i + 1, n):
            grad[e] = phi[j] - phi[i]
            e += 1
    f_sq = np.dot(flow, flow)
    g_sq = np.dot(grad, grad)
    gf = g_sq / f_sq if f_sq > 1e-30 else 0.0
    return {"phi": phi, "gf": gf}


def _edge_weight_stats(flow: np.ndarray) -> Dict[str, float]:
    """CV, skewness, kurtosis of edge weights."""
    w = np.abs(flow)
    w = w[w > 0]
    if len(w) < 10:
        return {"cv": float("nan"), "skew": float("nan"), "kurt": float("nan")}
    return {
        "cv": float(np.std(w) / np.mean(w)),
        "skew": float(skew(w)),
        "kurt": float(kurtosis(w)),
    }


def run_directional_decomposition(
    delta: np.ndarray,
    gene_names: List[str],
    output_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run directional decomposition: full Δ, Δ⁺ (gain), Δ⁻ (loss).

    Parameters
    ----------
    delta : np.ndarray, shape (N, N)
        Correlation change matrix.
    gene_names : list of str
    output_dir : Path, optional

    Returns
    -------
    dict with phi scores and GF for each component, plus cross-correlations.
    """
    N = len(gene_names)
    logger.info("Directional decomposition: N=%d genes", N)

    delta_plus, delta_minus = decompose_delta(delta)

    # Energy fractions
    total_energy = np.sum(delta ** 2)
    plus_energy = np.sum(delta_plus ** 2)
    minus_energy = np.sum(delta_minus ** 2)
    plus_frac = plus_energy / total_energy if total_energy > 0 else 0
    minus_frac = minus_energy / total_energy if total_energy > 0 else 0

    logger.info("  Energy: Δ⁺=%.1f%%, Δ⁻=%.1f%%", plus_frac * 100, minus_frac * 100)

    results = {}
    for label, component in [("full", delta), ("gain", delta_plus), ("loss", delta_minus)]:
        flow = _build_flow_edge_weight(component)
        hodge = _hodge_with_gf(flow, N)
        stats = _edge_weight_stats(flow)

        results[label] = {
            "phi": hodge["phi"],
            "gf": hodge["gf"],
            "edge_stats": stats,
        }
        logger.info("  %s: GF=%.4f, CV=%.3f", label, hodge["gf"], stats.get("cv", 0))

    # Cross-correlations between components
    rho_gain_loss, p_gl = spearmanr(results["gain"]["phi"], results["loss"]["phi"])
    rho_full_gain, p_fg = spearmanr(results["full"]["phi"], results["gain"]["phi"])
    rho_full_loss, p_fl = spearmanr(results["full"]["phi"], results["loss"]["phi"])

    cross = {
        "gain_vs_loss": {"rho": float(rho_gain_loss), "p": float(p_gl)},
        "full_vs_gain": {"rho": float(rho_full_gain), "p": float(p_fg)},
        "full_vs_loss": {"rho": float(rho_full_loss), "p": float(p_fl)},
    }
    logger.info("  Spearman: full↔gain=%.3f, full↔loss=%.3f, gain↔loss=%.3f",
                rho_full_gain, rho_full_loss, rho_gain_loss)

    result = {
        "n_genes": N,
        "energy_fraction": {"gain": float(plus_frac), "loss": float(minus_frac)},
        "components": {k: {"gf": v["gf"], "edge_stats": v["edge_stats"]} for k, v in results.items()},
        "phi_full": results["full"]["phi"],
        "phi_gain": results["gain"]["phi"],
        "phi_loss": results["loss"]["phi"],
        "cross_correlations": cross,
    }

    if output_dir is not None:
        import pandas as pd
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        pd.DataFrame({
            "gene": gene_names,
            "phi_full": results["full"]["phi"],
            "phi_gain": results["gain"]["phi"],
            "phi_loss": results["loss"]["phi"],
        }).sort_values("phi_full", ascending=False).to_csv(
            output_dir / "directional_phi.csv", index=False,
        )

        summary = {k: v for k, v in result.items() if k not in ("phi_full", "phi_gain", "phi_loss")}
        with open(output_dir / "directional_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)

        logger.info("Directional results saved to %s", output_dir)

    return result
