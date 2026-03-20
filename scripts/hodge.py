"""
scRNAseq Hodge Pipeline — Hodge Decomposition
===============================================
Discrete Hodge decomposition on the complete simplicial complex for
celltype-level flow analysis. Decomposes edge flow into gradient,
curl, and harmonic components.
"""
import logging
from itertools import combinations
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import nnls

from . import config
from .spd import spd_log

logger = logging.getLogger(__name__)


def _edge_index(i: int, j: int, n: int) -> int:
    """Map an ordered edge (i, j) with i < j to its linear index."""
    return i * n - i * (i + 1) // 2 + (j - i - 1)


def build_incidence_matrices(n_nodes: int) -> Tuple[np.ndarray, np.ndarray]:
    """Build B0 (node-edge) and B1 (edge-triangle) incidence matrices
    for the complete graph on n_nodes vertices."""
    n_edges = n_nodes * (n_nodes - 1) // 2
    n_triangles = n_nodes * (n_nodes - 1) * (n_nodes - 2) // 6

    B0 = np.zeros((n_edges, n_nodes), dtype=np.float64)
    for i, j in combinations(range(n_nodes), 2):
        e_idx = _edge_index(i, j, n_nodes)
        B0[e_idx, i] = -1.0
        B0[e_idx, j] = +1.0

    B1 = np.zeros((n_triangles, n_edges), dtype=np.float64)
    t_idx = 0
    for i, j, k in combinations(range(n_nodes), 3):
        e_ij = _edge_index(i, j, n_nodes)
        e_ik = _edge_index(i, k, n_nodes)
        e_jk = _edge_index(j, k, n_nodes)
        B1[t_idx, e_ij] = +1.0
        B1[t_idx, e_ik] = -1.0
        B1[t_idx, e_jk] = +1.0
        t_idx += 1

    logger.info("Built incidence matrices: %d nodes, %d edges, %d triangles",
                n_nodes, n_edges, n_triangles)
    return B0, B1


def hodge_decomposition(
    flow: np.ndarray,
    B0: np.ndarray,
    B1: np.ndarray,
) -> Dict[str, Any]:
    """Full Hodge decomposition of edge flow on a simplicial complex.

    Returns dict with: phi, gradient, curl_coeffs, curl, harmonic,
    gradient_fraction, curl_fraction, harmonic_fraction, r_phi.
    """
    phi, _, _, _ = np.linalg.lstsq(B0, flow, rcond=None)
    gradient = B0 @ phi

    if B1.shape[0] > 0:
        curl_coeffs, _, _, _ = np.linalg.lstsq(B1.T, flow, rcond=None)
        curl = B1.T @ curl_coeffs
    else:
        curl_coeffs = np.array([], dtype=np.float64)
        curl = np.zeros_like(flow)

    harmonic = flow - gradient - curl

    f_sq = np.dot(flow, flow)
    grad_sq = np.dot(gradient, gradient)
    curl_sq = np.dot(curl, curl)
    harm_sq = np.dot(harmonic, harmonic)

    total_sq = grad_sq + curl_sq + harm_sq
    if total_sq < 1e-30:
        gradient_fraction = curl_fraction = harmonic_fraction = 0.0
    else:
        gradient_fraction = float(grad_sq / total_sq)
        curl_fraction = float(curl_sq / total_sq)
        harmonic_fraction = float(harm_sq / total_sq)

    if f_sq < 1e-30 or grad_sq < 1e-30:
        r_phi = 0.0
    else:
        r_phi = float(np.corrcoef(flow, gradient)[0, 1])
        if np.isnan(r_phi):
            r_phi = 0.0

    logger.info("Hodge: gradient=%.3f, curl=%.3f, harmonic=%.3f, R(phi)=%.3f",
                gradient_fraction, curl_fraction, harmonic_fraction, r_phi)

    return {
        "phi": phi,
        "gradient": gradient,
        "curl_coeffs": curl_coeffs,
        "curl": curl,
        "harmonic": harmonic,
        "gradient_fraction": gradient_fraction,
        "curl_fraction": curl_fraction,
        "harmonic_fraction": harmonic_fraction,
        "r_phi": r_phi,
    }


def compute_spd_velocity(spd_windows: List[np.ndarray]) -> List[np.ndarray]:
    """Compute SPD velocity as consecutive differences in log space."""
    if len(spd_windows) < 2:
        raise ValueError(f"Need >= 2 windows, got {len(spd_windows)}")
    log_spds = [spd_log(S) for S in spd_windows]
    velocities = []
    for w in range(len(log_spds) - 1):
        v = log_spds[w + 1] - log_spds[w]
        v = (v + v.T) / 2.0
        velocities.append(v)
    return velocities


def build_flow_pairwise(
    d_corr_data: Dict[str, Dict],
    celltypes: List[str],
) -> np.ndarray:
    """Build edge flow from pairwise window-level d_corr comparisons."""
    n = len(celltypes)
    n_edges = n * (n - 1) // 2
    flow = np.zeros(n_edges, dtype=np.float64)

    ct_lookup = {}
    for ct in celltypes:
        data = d_corr_data.get(ct, {"d_corr": [], "window_pairs": []})
        pairs = data.get("window_pairs", [])
        d_corrs = data.get("d_corr", [])
        ct_lookup[ct] = {tuple(p): v for p, v in zip(pairs, d_corrs)}

    for i, j in combinations(range(n), 2):
        ct_i, ct_j = celltypes[i], celltypes[j]
        lookup_i = ct_lookup[ct_i]
        lookup_j = ct_lookup[ct_j]

        common_pairs = set(lookup_i.keys()) & set(lookup_j.keys())
        if not common_pairs:
            continue

        wins_i = wins_j = 0
        for p in common_pairs:
            if lookup_i[p] > lookup_j[p]:
                wins_i += 1
            elif lookup_j[p] > lookup_i[p]:
                wins_j += 1

        n_common = len(common_pairs)
        e_idx = _edge_index(i, j, n)
        flow[e_idx] = (wins_i - wins_j) / n_common

    logger.info("Pairwise flow: %d edges, range [%.4f, %.4f]",
                n_edges, flow.min(), flow.max())
    return flow


def permutation_test_gradient(
    flow: np.ndarray,
    B0: np.ndarray,
    B1: np.ndarray,
    n_perm: int = None,
    seed: int = 42,
) -> Dict[str, Any]:
    """Permutation test for gradient dominance."""
    if n_perm is None:
        n_perm = config.PERMUTATION_N

    rng = np.random.RandomState(seed)

    obs_result = hodge_decomposition(flow, B0, B1)
    obs_gf = obs_result["gradient_fraction"]
    obs_r = obs_result["r_phi"]

    null_gf = np.empty(n_perm, dtype=np.float64)
    null_r = np.empty(n_perm, dtype=np.float64)

    for p in range(n_perm):
        perm_flow = flow.copy()
        rng.shuffle(perm_flow)
        signs = rng.choice([-1.0, 1.0], size=len(perm_flow))
        perm_flow = perm_flow * signs

        perm_result = hodge_decomposition(perm_flow, B0, B1)
        null_gf[p] = perm_result["gradient_fraction"]
        null_r[p] = perm_result["r_phi"]

    p_value = float(np.mean(null_gf >= obs_gf))

    logger.info("Permutation test (%d): obs GF=%.4f, p=%.4f", n_perm, obs_gf, p_value)

    return {
        "p_value": p_value,
        "observed_gradient_fraction": obs_gf,
        "null_distribution": null_gf,
        "observed_r_phi": obs_r,
        "null_r_phi": null_r,
    }
