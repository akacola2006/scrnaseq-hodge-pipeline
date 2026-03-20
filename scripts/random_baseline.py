"""
scRNAseq Hodge Pipeline — Random Matrix Baseline
==================================================
Computes the null-expectation gradient fraction (GF) from random matrices
to establish whether observed GF represents genuine signal or structural
artifact of the complete graph K_N.

On K_N with N nodes, random edge-weight flows produce GF ≈ 1/φ (golden
ratio inverse ≈ 0.618) for sign-based flow, and GF ≈ 0.425 for
edge-weight flow. This is a geometric property of the graph, not signal.

The baseline is computed by:
1. Generating random symmetric matrices (same size as delta)
2. Running the same flow construction + Hodge decomposition
3. Building a null distribution of GF values
4. Comparing observed GF against this null
"""
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from scipy.stats import percentileofscore

from . import config
from .gene_hodge import hodge_gradient_kn, build_gene_flow

logger = logging.getLogger(__name__)


def _random_symmetric_matrix(n: int, rng: np.random.RandomState) -> np.ndarray:
    """Generate a random symmetric matrix with same statistical properties."""
    A = rng.randn(n, n)
    return (A + A.T) / 2.0


def compute_null_gf(
    n_genes: int,
    flow_mode: str = None,
    n_replicates: int = 200,
    seed: int = 42,
) -> Dict[str, Any]:
    """Compute null distribution of gradient fraction for random matrices.

    Parameters
    ----------
    n_genes : int
        Number of genes (determines graph size K_N).
    flow_mode : str
        Flow construction mode. Default: config.GENE_HODGE_FLOW_MODE.
    n_replicates : int
        Number of random matrices to generate.
    seed : int

    Returns
    -------
    dict with null_gf distribution, mean, std, percentiles.
    """
    if flow_mode is None:
        flow_mode = config.GENE_HODGE_FLOW_MODE

    rng = np.random.RandomState(seed)
    null_gfs = np.empty(n_replicates, dtype=np.float64)

    logger.info("Computing null GF baseline: N=%d, mode=%s, %d replicates",
                n_genes, flow_mode, n_replicates)

    n_edges = n_genes * (n_genes - 1) // 2

    for rep in range(n_replicates):
        delta_rand = _random_symmetric_matrix(n_genes, rng)
        flow = build_gene_flow(delta_rand, mode=flow_mode)
        phi = hodge_gradient_kn(flow, n_genes)

        # Compute GF
        grad = np.empty(n_edges, dtype=np.float64)
        e = 0
        for i in range(n_genes):
            for j in range(i + 1, n_genes):
                grad[e] = phi[j] - phi[i]
                e += 1
        f_sq = np.dot(flow, flow)
        g_sq = np.dot(grad, grad)
        null_gfs[rep] = g_sq / f_sq if f_sq > 1e-30 else 0.0

    result = {
        "null_gf": null_gfs,
        "mean": float(null_gfs.mean()),
        "std": float(null_gfs.std()),
        "median": float(np.median(null_gfs)),
        "p5": float(np.percentile(null_gfs, 5)),
        "p95": float(np.percentile(null_gfs, 95)),
        "n_genes": n_genes,
        "flow_mode": flow_mode,
        "n_replicates": n_replicates,
    }

    logger.info("  Null GF: mean=%.4f, std=%.4f, [p5=%.4f, p95=%.4f]",
                result["mean"], result["std"], result["p5"], result["p95"])

    return result


def compare_to_baseline(
    observed_gf: float,
    null_result: Dict[str, Any],
) -> Dict[str, Any]:
    """Compare observed GF to null baseline.

    Returns
    -------
    dict with:
        signal_above_null: bool
        observed_gf, null_mean, excess_gf, z_score, percentile_rank
    """
    null_gfs = null_result["null_gf"]
    null_mean = null_result["mean"]
    null_std = null_result["std"]

    excess = observed_gf - null_mean
    z_score = excess / null_std if null_std > 1e-10 else 0.0
    pct_rank = percentileofscore(null_gfs, observed_gf)

    signal = observed_gf > null_result["p95"]

    result = {
        "observed_gf": observed_gf,
        "null_mean": null_mean,
        "null_std": null_std,
        "excess_gf": excess,
        "z_score": z_score,
        "percentile_rank": pct_rank,
        "signal_above_null": signal,
        "interpretation": (
            f"GF={observed_gf:.4f} is {'ABOVE' if signal else 'WITHIN'} "
            f"null baseline (mean={null_mean:.4f}, p95={null_result['p95']:.4f}). "
            f"{'Genuine gradient signal detected.' if signal else 'No significant gradient signal beyond structural baseline.'}"
        ),
    }

    logger.info("  Baseline comparison: %s", result["interpretation"])

    return result


def run_random_baseline(
    observed_gf: float,
    n_genes: int,
    flow_mode: str = None,
    n_replicates: int = 200,
    output_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Full random matrix baseline analysis.

    Parameters
    ----------
    observed_gf : float
        The GF from the actual gene Hodge analysis.
    n_genes : int
    flow_mode : str
    n_replicates : int
    output_dir : Path, optional

    Returns
    -------
    dict with null distribution and comparison results.
    """
    null_result = compute_null_gf(n_genes, flow_mode, n_replicates)
    comparison = compare_to_baseline(observed_gf, null_result)

    result = {
        "null": {k: v for k, v in null_result.items() if k != "null_gf"},
        "comparison": comparison,
    }

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        np.save(output_dir / "null_gf_distribution.npy", null_result["null_gf"])

        with open(output_dir / "random_baseline_summary.json", "w") as f:
            json.dump(result, f, indent=2, default=str)

        logger.info("Random baseline results saved to %s", output_dir)

    return result
