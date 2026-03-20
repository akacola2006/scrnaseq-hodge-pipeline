"""
scRNAseq Hodge Pipeline — Sparse k-NN Hodge Decomposition (Dual-Mode)
=======================================================================
Hodge decomposition on k-NN sparsified gene graph, providing a dual-mode
analysis alongside the default complete graph K_N.

On K_N, gradient fraction is bounded by a structural baseline (~1/phi on
large N). Sparse graphs break this degeneracy and reveal local curl
structure invisible to K_N.

Two d_i variants:
  (a) "full"  — d_i from all columns of delta (same input as K_N, only graph changes)
  (b) "local" — d_i from connected edges only (graph + flow both change)

Analyses:
  - k sweep: GF, phi preservation (Spearman vs K_N phi) across k values
  - Curl structure: triangle curls, curl hub scores
  - Dual-mode comparison: K_N vs optimal sparse
"""
import json
import logging
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import lsqr
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr

from . import config
from .gene_hodge import hodge_gradient_kn, build_gene_flow

logger = logging.getLogger(__name__)


# =====================================================================
# 1. k-NN Graph Construction
# =====================================================================

def build_knn_graph(
    delta: np.ndarray,
    k: int,
) -> Tuple[np.ndarray, int, np.ndarray]:
    """Build symmetric k-NN graph from delta row distances.

    Parameters
    ----------
    delta : np.ndarray, shape (N, N)
        Change in log-correlation between windows.
    k : int
        Number of nearest neighbours per node.

    Returns
    -------
    edges : np.ndarray, shape (n_edges, 2)
        Upper-triangle edge list (i < j).
    n_components : int
        Number of connected components.
    row_dist : np.ndarray, shape (N, N)
        Pairwise row-distance matrix (reusable).
    """
    n = delta.shape[0]

    # Row-distance matrix
    row_dist = squareform(pdist(delta, metric="euclidean"))

    # k-NN per node (exclude self), symmetrised
    adjacency = lil_matrix((n, n), dtype=np.float64)
    for i in range(n):
        dists = row_dist[i].copy()
        dists[i] = np.inf
        neighbours = np.argsort(dists)[:k]
        for j in neighbours:
            adjacency[i, j] = 1
            adjacency[j, i] = 1

    adj_csr = adjacency.tocsr()

    # Connected components
    n_components, _ = connected_components(adj_csr, directed=False)

    # Extract upper-triangle edges
    rows_arr, cols_arr = adj_csr.nonzero()
    mask = rows_arr < cols_arr
    edges = np.column_stack([rows_arr[mask], cols_arr[mask]])

    logger.info("  k-NN graph: k=%d, |E|=%d (of %d possible), components=%d, avg_degree=%.1f",
                k, len(edges), n * (n - 1) // 2, n_components, 2 * len(edges) / n)

    return edges, n_components, row_dist


# =====================================================================
# 2. Flow Construction on Sparse Graph
# =====================================================================

def build_flow_sparse(
    delta: np.ndarray,
    edges: np.ndarray,
    d_mode: str = "full",
) -> Tuple[np.ndarray, np.ndarray]:
    """Build edge-weight flow on sparse graph.

    Parameters
    ----------
    delta : np.ndarray, shape (N, N)
        Correlation change matrix.
    edges : np.ndarray, shape (n_edges, 2)
        Edge list.
    d_mode : str
        "full" = d_i from all columns, "local" = d_i from connected edges only.

    Returns
    -------
    flow : np.ndarray, shape (n_edges,)
    d_corr : np.ndarray, shape (N,)
        Per-node correlation-change magnitude.
    """
    n = delta.shape[0]

    if d_mode == "full":
        d_corr = np.linalg.norm(delta, axis=1)
    else:
        d_sq = np.zeros(n, dtype=np.float64)
        for i, j in edges:
            val_sq = delta[i, j] ** 2
            d_sq[i] += val_sq
            d_sq[j] += val_sq
        d_corr = np.sqrt(d_sq)

    # Flow: f(i,j) = |delta_ij| * sign(d_i - d_j)
    w_ij = np.abs(delta[edges[:, 0], edges[:, 1]])
    diff = d_corr[edges[:, 0]] - d_corr[edges[:, 1]]
    flow = w_ij * np.sign(diff)

    return flow, d_corr


# =====================================================================
# 3. Sparse Hodge Decomposition
# =====================================================================

def _build_incidence(edges: np.ndarray, n_nodes: int) -> csr_matrix:
    """Build oriented incidence matrix B0: (n_edges, n_nodes)."""
    n_edges = len(edges)
    rows = np.concatenate([np.arange(n_edges), np.arange(n_edges)])
    cols = np.concatenate([edges[:, 0], edges[:, 1]])
    vals = np.concatenate([-np.ones(n_edges), np.ones(n_edges)])
    return csr_matrix((vals, (rows, cols)), shape=(n_edges, n_nodes))


def sparse_hodge_decomposition(
    edges: np.ndarray,
    flow: np.ndarray,
    n_nodes: int,
) -> Dict[str, Any]:
    """Hodge decomposition on arbitrary graph via LSQR.

    Decomposes f = gradient + residual, where gradient = B0 @ phi.

    Returns
    -------
    dict with: phi, gradient, residual, gf (gradient fraction), rf (residual fraction).
    """
    B0 = _build_incidence(edges, n_nodes)

    # Divergence: div = B0^T @ f
    div = B0.T @ flow

    # Graph Laplacian: L = B0^T @ B0
    L = B0.T @ B0

    # Solve L phi = div via LSQR
    phi_raw, *_ = lsqr(L, div, atol=1e-12, btol=1e-12)
    phi = phi_raw - phi_raw.mean()

    gradient = B0 @ phi
    residual = flow - gradient

    flow_sq = float(np.dot(flow, flow))
    grad_sq = float(np.dot(gradient, gradient))
    resid_sq = float(np.dot(residual, residual))

    gf = grad_sq / flow_sq if flow_sq > 1e-30 else 0.0
    rf = resid_sq / flow_sq if flow_sq > 1e-30 else 0.0

    logger.info("  Sparse Hodge: GF=%.4f, RF=%.4f", gf, rf)

    return {
        "phi": phi,
        "gradient": gradient,
        "residual": residual,
        "gf": gf,
        "rf": rf,
        "flow_norm_sq": flow_sq,
    }


# =====================================================================
# 4. Triangle / Curl Analysis
# =====================================================================

def find_triangles(edges: np.ndarray, n_nodes: int) -> np.ndarray:
    """Find all triangles in the graph. Returns (n_triangles, 3) with i<j<k."""
    adj = [set() for _ in range(n_nodes)]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)

    triangles = []
    for i, j in edges:
        common = adj[i] & adj[j]
        for k_node in common:
            if k_node > j:
                triangles.append((i, j, k_node))

    if triangles:
        return np.array(triangles, dtype=np.int32)
    return np.empty((0, 3), dtype=np.int32)


def compute_triangle_curls(
    edges: np.ndarray,
    flow: np.ndarray,
    triangles: np.ndarray,
) -> np.ndarray:
    """Compute curl (circulation) for each triangle."""
    if len(triangles) == 0:
        return np.array([], dtype=np.float64)

    edge_map = {(i, j): idx for idx, (i, j) in enumerate(edges)}

    curls = np.zeros(len(triangles), dtype=np.float64)
    for t_idx, (i, j, k) in enumerate(triangles):
        f_ij = flow[edge_map[(i, j)]]
        f_ik = flow[edge_map[(i, k)]]
        jk_key = (min(j, k), max(j, k))
        f_jk = flow[edge_map[jk_key]]
        curls[t_idx] = f_ij + f_jk - f_ik

    return curls


def compute_curl_hub_scores(
    triangle_curls: np.ndarray,
    triangles: np.ndarray,
    n_genes: int,
    percentile_threshold: float = 99.0,
) -> Dict[str, Any]:
    """Curl hub scores: genes participating in high-curl triangles."""
    if len(triangle_curls) == 0:
        return {
            "curl_hub": np.zeros(n_genes),
            "curl_energy": np.zeros(n_genes),
            "n_triangles_total": 0,
        }

    abs_curls = np.abs(triangle_curls)
    threshold = np.percentile(abs_curls, percentile_threshold)
    high_mask = abs_curls >= threshold

    curl_energy = np.zeros(n_genes, dtype=np.float64)
    curl_sq = triangle_curls.astype(np.float64) ** 2
    np.add.at(curl_energy, triangles[:, 0], curl_sq)
    np.add.at(curl_energy, triangles[:, 1], curl_sq)
    np.add.at(curl_energy, triangles[:, 2], curl_sq)

    high_indices = triangles[high_mask]
    n_high_per_gene = np.zeros(n_genes, dtype=np.int64)
    if len(high_indices) > 0:
        np.add.at(n_high_per_gene, high_indices[:, 0], 1)
        np.add.at(n_high_per_gene, high_indices[:, 1], 1)
        np.add.at(n_high_per_gene, high_indices[:, 2], 1)

    max_per_gene = max(n_high_per_gene.max(), 1)
    curl_hub = n_high_per_gene.astype(np.float64) / max_per_gene

    return {
        "curl_hub": curl_hub,
        "curl_energy": curl_energy,
        "n_high_triangles": n_high_per_gene,
        "n_triangles_total": len(triangle_curls),
        "n_high_total": int(high_mask.sum()),
        "threshold": float(threshold),
    }


# =====================================================================
# 5. Dual-Mode Analysis: K_N vs Sparse
# =====================================================================

def run_dual_mode(
    delta: np.ndarray,
    gene_names: List[str],
    k_values: List[int] = None,
    d_mode: str = "full",
    output_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run dual-mode Hodge analysis: K_N (complete graph) vs sparse k-NN.

    Parameters
    ----------
    delta : np.ndarray, shape (N, N)
        Correlation change matrix (log(corr_{w+1}) - log(corr_w)).
    gene_names : list of str
        Gene identifiers.
    k_values : list of int, optional
        k values for k-NN sweep. Default: [10, 20, 30, 50, 100].
    d_mode : str
        "full" or "local" for node strength computation.
    output_dir : Path, optional
        Save results here.

    Returns
    -------
    dict with K_N results, sparse sweep results, and comparison metrics.
    """
    if k_values is None:
        k_values = [10, 20, 30, 50, 100]

    N = delta.shape[0]
    logger.info("Dual-mode Hodge: N=%d genes, k_values=%s, d_mode=%s", N, k_values, d_mode)

    # ── K_N baseline ──────────────────────────────────────────
    logger.info("[K_N] Complete graph analysis...")
    flow_kn = build_gene_flow(delta, mode="edge_weight")
    phi_kn = hodge_gradient_kn(flow_kn, N)

    # K_N gradient fraction
    n_edges_kn = N * (N - 1) // 2
    grad_kn = np.empty(n_edges_kn, dtype=np.float64)
    e_idx = 0
    for i in range(N):
        for j in range(i + 1, N):
            grad_kn[e_idx] = phi_kn[j] - phi_kn[i]
            e_idx += 1

    f_sq_kn = np.dot(flow_kn, flow_kn)
    g_sq_kn = np.dot(grad_kn, grad_kn)
    gf_kn = g_sq_kn / f_sq_kn if f_sq_kn > 1e-30 else 0.0

    logger.info("[K_N] GF=%.4f, phi range=[%.4f, %.4f]", gf_kn, phi_kn.min(), phi_kn.max())

    kn_result = {
        "phi": phi_kn,
        "gf": gf_kn,
        "n_edges": n_edges_kn,
    }

    # ── Sparse k-NN sweep ─────────────────────────────────────
    sweep_results = []

    for k in k_values:
        if k >= N:
            logger.warning("  k=%d >= N=%d, skipping (would be complete graph)", k, N)
            continue

        logger.info("[Sparse] k=%d...", k)
        edges, n_comp, _ = build_knn_graph(delta, k)

        if len(edges) == 0:
            logger.warning("  k=%d: no edges, skipping", k)
            continue

        flow_sp, d_corr = build_flow_sparse(delta, edges, d_mode=d_mode)
        hodge_sp = sparse_hodge_decomposition(edges, flow_sp, N)

        # Phi preservation: Spearman correlation with K_N phi
        rho, pval = spearmanr(phi_kn, hodge_sp["phi"])

        # Curl analysis
        triangles = find_triangles(edges, N)
        if len(triangles) > 0:
            t_curls = compute_triangle_curls(edges, flow_sp, triangles)
            curl_hub = compute_curl_hub_scores(t_curls, triangles, N)
        else:
            t_curls = np.array([])
            curl_hub = {"n_triangles_total": 0, "curl_hub": np.zeros(N)}

        entry = {
            "k": k,
            "n_edges": len(edges),
            "n_components": n_comp,
            "gf": hodge_sp["gf"],
            "rf": hodge_sp["rf"],
            "phi": hodge_sp["phi"],
            "phi_spearman_vs_kn": float(rho),
            "phi_spearman_pval": float(pval),
            "n_triangles": curl_hub["n_triangles_total"],
            "curl_hub_scores": curl_hub["curl_hub"],
        }
        sweep_results.append(entry)

        logger.info("  k=%d: GF=%.4f, phi_rho=%.4f, triangles=%d",
                    k, hodge_sp["gf"], rho, curl_hub["n_triangles_total"])

    # ── Select optimal k ──────────────────────────────────────
    # Criterion: highest GF with phi_spearman >= 0.8
    candidates = [r for r in sweep_results if r["phi_spearman_vs_kn"] >= 0.8]
    if candidates:
        optimal = max(candidates, key=lambda r: r["gf"])
    elif sweep_results:
        optimal = max(sweep_results, key=lambda r: r["phi_spearman_vs_kn"])
    else:
        optimal = None

    # ── Assemble result ───────────────────────────────────────
    result = {
        "kn": kn_result,
        "sweep": sweep_results,
        "optimal_k": optimal["k"] if optimal else None,
        "optimal_gf": optimal["gf"] if optimal else None,
        "optimal_phi_rho": optimal["phi_spearman_vs_kn"] if optimal else None,
        "n_genes": N,
        "d_mode": d_mode,
    }

    # ── Save ──────────────────────────────────────────────────
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Summary JSON
        summary = {
            "n_genes": N,
            "d_mode": d_mode,
            "kn_gf": gf_kn,
            "optimal_k": result["optimal_k"],
            "optimal_gf": result["optimal_gf"],
            "optimal_phi_rho": result["optimal_phi_rho"],
            "sweep": [
                {
                    "k": r["k"],
                    "n_edges": r["n_edges"],
                    "n_components": r["n_components"],
                    "gf": r["gf"],
                    "rf": r["rf"],
                    "phi_spearman_vs_kn": r["phi_spearman_vs_kn"],
                    "n_triangles": r["n_triangles"],
                }
                for r in sweep_results
            ],
        }
        with open(output_dir / "dual_mode_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        # Gene-level comparison (K_N phi vs optimal sparse phi)
        if optimal:
            import pandas as pd
            comparison = pd.DataFrame({
                "gene": gene_names,
                "phi_kn": phi_kn,
                "phi_sparse": optimal["phi"],
                "curl_hub_sparse": optimal["curl_hub_scores"],
            })
            comparison = comparison.sort_values("phi_kn", ascending=False)
            comparison.to_csv(output_dir / "dual_mode_gene_comparison.csv", index=False)

        logger.info("Dual-mode results saved to %s", output_dir)

    return result
