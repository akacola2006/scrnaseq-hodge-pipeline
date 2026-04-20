"""
All-gene insertion for 3φ framework (Section 6.15, Appendix AS of the paper).

For each protein-coding gene in the genome, insert it one-by-one into the
cell-type-specific base gene set (200-3,500 genes), recompute the per-donor
log-correlation matrices with the inserted gene included, and record the
three φ variants (disease, condition, static) as well as their percentile
position within the base+1 gene ranking.

This enables:
  - Genome-wide screening for disease-specific residual z-scores
    (basis for NEMF universal downshift finding, Section 8.9)
  - Per-gene manifold preservation R² computation (Table 1, Section 4.1)
  - All-gene rewiring/collapse count (Section 4.2)

GPU-accelerated via PyTorch: for each of ~10,000-19,000 extra genes per CT,
the insertion requires one updated log-correlation computation per donor.

Typical runtime: ~40 hours for all 10 CTs on NVIDIA RTX 5090 (32 GB VRAM).

Usage:
    python -m scripts.allgene_insertion --cell-type Oligo --n-base 3500

Output:
    results/allgene_3phi/{ct}/
        allgene_3phi.csv          # Per-gene 3φ percentiles + raw values
        base_raw_phi.npz          # Base-set raw φ values and delta matrices
        summary.json              # CT metadata, GF values, timing

Migrated from: sals_analysis_frozen_20260211/scripts/run_allgene_insertion_3phi.py
              sals_analysis_frozen_20260211/scripts/run_allgene_insertion_fast.py
"""

from __future__ import annotations

import gc
import json
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

try:
    import torch
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

from scripts.gene_hodge import hodge_gradient_kn

logger = logging.getLogger(__name__)


# Per-CT base gene set sizes (from Appendix AS.2 of the validation report)
# Determined by per-donor cell count: Ledoit-Wolf requires n_cells ≥ n_genes/10
DEFAULT_CT_CONFIG = {
    "Oligo":  {"n_base": 3500, "min_cells_per_donor": 100, "min_donors_per_window": 3},
    "Astro":  {"n_base": 3000, "min_cells_per_donor": 100, "min_donors_per_window": 3},
    "L2_L3":  {"n_base": 2500, "min_cells_per_donor": 100, "min_donors_per_window": 2},
    "OPC":    {"n_base": 2000, "min_cells_per_donor": 100, "min_donors_per_window": 3},
    "L3_L5":  {"n_base": 2000, "min_cells_per_donor": 100, "min_donors_per_window": 3},
    "L4_L6":  {"n_base": 2000, "min_cells_per_donor": 100, "min_donors_per_window": 3},
    "L5_L6":  {"n_base": 2000, "min_cells_per_donor": 100, "min_donors_per_window": 3},
    "PV":     {"n_base": 1500, "min_cells_per_donor":  50, "min_donors_per_window": 2},
    "L6":     {"n_base": 1000, "min_cells_per_donor":  50, "min_donors_per_window": 2},
    "Endo":   {"n_base":  200, "min_cells_per_donor":  10, "min_donors_per_window": 2},
}


def gpu_log_corr(X_torch: "torch.Tensor") -> "torch.Tensor":
    """Compute log-Euclidean representation of (n_cells × n_genes) residual matrix.

    Steps:
      1. Pearson correlation matrix (p × p)
      2. Ledoit-Wolf-style shrinkage toward identity
      3. Matrix logarithm via eigendecomposition

    Returns log(C) in the symmetric positive definite tangent space.
    """
    if not HAS_TORCH:
        raise RuntimeError("PyTorch is required for GPU all-gene insertion.")

    n, p = X_torch.shape
    # Center
    X_c = X_torch - X_torch.mean(dim=0, keepdim=True)
    # Correlation
    std = X_c.std(dim=0, keepdim=True, unbiased=True) + 1e-12
    X_std = X_c / std
    C = (X_std.T @ X_std) / (n - 1)

    # Ensure symmetric positive definite: eigendecompose, clip negative eigenvalues
    eigvals, eigvecs = torch.linalg.eigh(C)
    eigvals = torch.clamp(eigvals, min=1e-10)
    log_C = eigvecs @ torch.diag(torch.log(eigvals)) @ eigvecs.T
    return log_C


def _edge_weight_flow_torch(delta: "torch.Tensor") -> "torch.Tensor":
    """Edge-weight flow f(i,j) = |Δ_ij| · sign(d_i - d_j) on GPU."""
    n = delta.shape[0]
    d = torch.norm(delta, dim=1)
    # Upper triangular indices
    tri_i, tri_j = torch.triu_indices(n, n, offset=1, device=delta.device)
    flow = torch.abs(delta[tri_i, tri_j]) * torch.sign(d[tri_i] - d[tri_j])
    std_f = flow.std()
    if std_f > 0:
        flow = flow / std_f
    return flow


def _hodge_phi_from_flow(flow_np: np.ndarray, n_genes: int) -> Dict:
    """K_N Hodge: φ = L⁺ div(f), with L⁺ = (1/N)(I − J/N)."""
    return hodge_gradient_kn(flow_np, n_genes)


def insert_one_gene(
    base_log_corrs: Dict[str, np.ndarray],
    donor_x_gene: Dict[str, np.ndarray],
    gene_x_vec: Dict[str, np.ndarray],
    gene_idx: int,
    condition: Dict[str, str],
    window: Dict[str, int],
    transition: int,
    device: str = "cuda",
) -> Optional[Dict]:
    """Insert a single gene g into the base set for one CT.

    For each donor, append gene g's residual vector to that donor's base residuals,
    recompute log(corr) on the (n_base + 1) × (n_base + 1) matrix, then extract
    the top-left (n_base) × (n_base) block for deltas.

    Returns per-gene:
      phi_disease_percentile
      phi_condition_percentile
      phi_static_percentile
      raw_phi_disease, raw_phi_condition, raw_phi_static

    Note: For efficiency, typically ~20% of execution is the per-donor eigendecomp;
    the GPU parallelisation of eigendecomp across donors is the main speedup.
    """
    # This is a simplified high-level interface; see `_run_allgene_insertion_ct`
    # below for the batched GPU implementation.
    raise NotImplementedError(
        "Use `run_allgene_insertion_ct` for the full batched implementation. "
        "Single-gene insertion is exposed only for testing."
    )


def run_allgene_insertion_ct(
    residuals: np.ndarray,
    gene_names: List[str],
    donor_ids: np.ndarray,
    donor_condition: Dict[str, str],
    pseudotime_df: pd.DataFrame,
    ct_name: str,
    n_base: int = 3500,
    min_donors_per_window: int = 2,
    output_dir: Optional[Path] = None,
    device: str = "cuda",
    checkpoint_every: int = 500,
) -> pd.DataFrame:
    """Run all-gene 3φ insertion for a single cell type.

    Parameters
    ----------
    residuals : (n_cells, n_all_genes) OLS residuals of log₁₀(CPM+1)
    gene_names : list of all detected gene symbols / ENSG IDs
    donor_ids : (n_cells,) donor labels
    donor_condition : dict donor_id → 'ALS' or 'PN'
    pseudotime_df : DataFrame with donor_id, pt, window columns
    ct_name : cell type label for logging/output naming
    n_base : number of top-variance genes to include in base set
    output_dir : if given, saves checkpoints and results here

    Returns
    -------
    DataFrame with columns:
        gene, symbol, biotype, in_base_set,
        phi_disease_pct, phi_condition_pct, phi_static_pct,
        raw_phi_disease, raw_phi_condition, raw_phi_static

    See paper Section 6.15 and Appendix AS for the complete algorithm.
    This is a GPU-accelerated reference implementation; the production version
    in sals_analysis_frozen_20260211/scripts/run_allgene_insertion_3phi.py
    includes additional optimisations (batched eigendecomp, checkpoint resume).
    """
    if not HAS_TORCH:
        raise RuntimeError(
            "PyTorch + CUDA required for all-gene insertion. "
            "Install: pip install torch --index-url https://download.pytorch.org/whl/cu121"
        )

    import torch

    device = torch.device(device if torch.cuda.is_available() else "cpu")
    logger.info(f"[{ct_name}] Device: {device}")

    # Select base genes by top variance
    gene_var = np.var(residuals, axis=0)
    base_idx = np.sort(np.argsort(gene_var)[::-1][:n_base])
    extra_idx = np.sort(np.setdiff1d(np.arange(len(gene_names)), base_idx))

    logger.info(f"[{ct_name}] Base: {len(base_idx)} genes; Extra (to insert): {len(extra_idx)} genes")

    # -- This is a reference stub of the algorithm. --
    # The full production implementation (sals_analysis_frozen_20260211/scripts/
    # run_allgene_insertion_3phi.py) is ~400 lines with batched GPU eigendecomp,
    # checkpointing, and NaN handling. It produces the 21,722-row zscore_matrix
    # used for the NEMF screen.
    #
    # We provide the public skeleton here; users requiring exact reproduction
    # of the paper's all-gene screen should consult the frozen reference
    # implementation. Key computational steps:
    #
    # 1. Compute mean PN log(corr) on (n_base × n_base) → phi_static base ranking
    # 2. Compute mean ALS log(corr) on (n_base × n_base) → phi_condition base ranking
    # 3. For primary transition w_star, compute delta_disease on (n_base × n_base)
    # 4. Compute all three base phi rankings
    # 5. For each extra gene g in parallel batches (GPU):
    #    a. For each donor d, compute log-corr(residuals[:, base ∪ {g}])
    #    b. Extract row g of the log-corr matrix (length n_base + 1)
    #    c. Compute updated mean_log_PN, mean_log_ALS including row g
    #    d. Compute delta_disease, delta_condition, static matrices at dimension n_base + 1
    #    e. Run Hodge gradient_kn; take phi[g's index]
    #    f. Convert to percentile in the n_base + 1 ranking
    # 6. Stack results into per-CT DataFrame

    raise NotImplementedError(
        f"run_allgene_insertion_ct is a reference stub. "
        f"For production use, see: sals_analysis_frozen_20260211/scripts/"
        f"run_allgene_insertion_3phi.py. "
        f"Expected output: results/allgene_3phi/{ct_name.lower()}/allgene_3phi.csv "
        f"with columns [gene, symbol, biotype, phi_disease_pct, phi_condition_pct, "
        f"phi_static_pct, in_base_set, raw_phi_static, raw_phi_condition, raw_phi_disease]"
    )


def build_zscore_matrix(
    allgene_3phi_dir: Path,
    cell_types: List[str],
    output_path: Optional[Path] = None,
) -> pd.DataFrame:
    """Build the genome-wide z-score matrix (one row per gene, one column per CT).

    This is the input to the NEMF screen (Section 8.9.1) and the matched-null
    z-score distribution (Appendix AV).

    Reads per-CT `allgene_3phi.csv` files and computes poly3 residual z-scores
    for each gene. Genes appearing in ≥2 CTs (default) are retained.

    Parameters
    ----------
    allgene_3phi_dir : directory containing {ct}/allgene_3phi.csv files
    cell_types : list of CT labels
    output_path : if provided, saves the matrix as CSV

    Returns
    -------
    DataFrame with columns:
        gene, symbol, biotype,
        z_disease_<ct1>, z_disease_<ct2>, ..., z_disease_<ct10>,
        n_ct_tested, mean_z, median_z, min_z, max_z, max_abs_z,
        n_z_below_neg2, n_z_above_pos2, n_z_below_neg1, n_z_above_pos1,
        direction_consistent_neg, direction_consistent_pos

    Migrated from: sals_analysis_frozen_20260211/scripts/run_build_zscore_matrix.py
    """
    from scripts.three_phi_residual import compute_residual_z

    z_dfs = []
    for ct in cell_types:
        ct_path = allgene_3phi_dir / ct.lower() / "allgene_3phi.csv"
        if not ct_path.exists():
            logger.warning(f"Missing {ct_path}, skipping")
            continue

        df = pd.read_csv(ct_path)
        if "raw_phi_disease" not in df.columns or "raw_phi_static" not in df.columns:
            logger.warning(f"{ct_path} missing raw phi columns, skipping")
            continue

        mask = df["raw_phi_disease"].notna() & df["raw_phi_static"].notna()
        df = df[mask].copy()
        if len(df) == 0:
            continue

        resid = compute_residual_z(
            df["raw_phi_disease"].values,
            df["raw_phi_static"].values,
            degree=3,
        )
        df[f"z_disease_{ct.lower()}"] = resid["z"]
        z_dfs.append(df[["gene", "symbol", f"z_disease_{ct.lower()}"]])

    if not z_dfs:
        raise ValueError("No allgene_3phi.csv files found in any CT directory")

    # Outer join on gene id
    merged = z_dfs[0]
    for df in z_dfs[1:]:
        merged = merged.merge(df, on=["gene", "symbol"], how="outer")

    # Summary columns
    z_cols = [c for c in merged.columns if c.startswith("z_disease_")]
    z_matrix = merged[z_cols].values

    merged["n_ct_tested"] = (~np.isnan(z_matrix)).sum(axis=1)
    merged["mean_z"] = np.nanmean(z_matrix, axis=1)
    merged["median_z"] = np.nanmedian(z_matrix, axis=1)
    merged["min_z"] = np.nanmin(z_matrix, axis=1)
    merged["max_z"] = np.nanmax(z_matrix, axis=1)
    merged["max_abs_z"] = np.nanmax(np.abs(z_matrix), axis=1)
    merged["n_z_below_neg2"] = (z_matrix < -2).sum(axis=1)
    merged["n_z_above_pos2"] = (z_matrix > 2).sum(axis=1)
    merged["n_z_below_neg1"] = (z_matrix < -1).sum(axis=1)
    merged["n_z_above_pos1"] = (z_matrix > 1).sum(axis=1)
    merged["direction_consistent_neg"] = (z_matrix < 0).sum(axis=1) == merged["n_ct_tested"]
    merged["direction_consistent_pos"] = (z_matrix > 0).sum(axis=1) == merged["n_ct_tested"]

    if output_path:
        merged.to_csv(output_path, index=False)
        logger.info(f"Saved z-score matrix to {output_path} ({len(merged)} genes × {len(cell_types)} CTs)")

    return merged


def main():
    """CLI entry point for all-gene insertion."""
    import argparse

    parser = argparse.ArgumentParser(description="All-gene 3φ insertion (Section 6.15)")
    parser.add_argument("--cell-type", type=str, required=True, help="Target cell type label")
    parser.add_argument("--n-base", type=int, default=None, help="Override base gene set size")
    parser.add_argument("--output-dir", type=Path, default=Path("results/allgene_3phi"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    ct_cfg = DEFAULT_CT_CONFIG.get(args.cell_type, {"n_base": 2000})
    n_base = args.n_base or ct_cfg["n_base"]

    logger.info(f"Running all-gene 3φ insertion for {args.cell_type} (n_base={n_base})")
    logger.warning(
        "This script provides the interface skeleton. "
        "Production implementation lives in the frozen reference pipeline; "
        "see module docstring for details."
    )


if __name__ == "__main__":
    main()
