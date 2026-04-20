"""
NEMF genome-wide screen and matched-null z-score analysis (Section 8.9, Appendix AV).

Inputs the per-CT `allgene_3phi.csv` files produced by `allgene_insertion.py`
(or the precomputed `zscore_matrix_wide.csv` built by
`allgene_insertion.build_zscore_matrix`), and performs:

  1. Per-gene matched-null z-score using a static-phi-matched window
     of neighbouring genes as the null reference (Section 6.15).
  2. Genome-wide multiple-testing correction (Bonferroni + binomial null)
     for the "z < -2 in ≥K cell types" criterion (Appendix AV).
  3. Cross-CT ρ correlation analysis: for each gene, compute Spearman ρ
     between its 10-dimensional z-vector and NEMF's z-vector across CTs,
     identifying the co-collapse axis (positive ρ) and compensatory axis
     (negative ρ; Section 8.9.5, Figure 6B).
  4. Pathway-level matched-null aggregation via one-sample t-test
     (Translation/Ribosome 9-10/10 CTs, Chaperone 4/10, ATP synthase 3/10;
     Section 4.2, Figure 5).

Output:
    results/nemf_screen/
        matched_null_z.csv              # Per-gene matched-null z per CT
        nemf_bonferroni.csv             # Individual-CT Bonferroni correction
        genome_wide_summary.json        # Binomial null FDR computation
        cross_ct_correlation.csv        # Cross-CT ρ(gene_z, NEMF_z)
        pathway_test_results.csv        # Pathway-level one-sample t-test

Migrated from: sals_analysis_frozen_20260211/scripts/run_3phi_verification.py
              sals_analysis_frozen_20260211/SALS_VALIDATION_REPORT.md (Appendix AV)
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_1samp, norm, binom

logger = logging.getLogger(__name__)


def matched_null_z_static_window(
    df: pd.DataFrame,
    gene_symbol: str,
    target_col: str = "residual_z",
    static_col: str = "raw_phi_static",
    window: float = 0.002,
    expand_window: float = 0.003,
    min_null_genes: int = 20,
    absolute_min: int = 5,
) -> tuple[float, int]:
    """Compute matched-null z-score using genes at similar static φ position.

    For a target gene, find all other genes whose raw_phi_static is within ±`window`
    (expand to `expand_window` if fewer than `min_null_genes` genes match). Treat
    these as the null reference. Compute:

        matched_z = (target_z - null_mean) / null_std

    This corrects for position-dependent variance structure (high-φ genes tend
    to have narrower residual distributions than low-φ genes).

    Parameters
    ----------
    df : DataFrame with symbol, target_col, static_col columns
    gene_symbol : target gene symbol
    target_col : column containing the residual or raw z-score to standardise
    static_col : column containing the static φ coordinate for matching
    window : initial window width (absolute difference in static φ)
    expand_window : expanded window if initial yields too few null genes
    min_null_genes : preferred minimum null gene count
    absolute_min : absolute minimum null gene count (below which return NaN)

    Returns
    -------
    matched_z : float (NaN if insufficient null)
    n_null    : number of genes used in null
    """
    row = df[df["symbol"] == gene_symbol]
    if len(row) == 0:
        return np.nan, 0

    r = row.iloc[0]
    if pd.isna(r.get(target_col, np.nan)) or pd.isna(r.get(static_col, np.nan)):
        return np.nan, 0

    sta = r[static_col]
    mask = (
        (np.abs(df[static_col] - sta) < window)
        & (df["symbol"] != gene_symbol)
        & df[target_col].notna()
    )
    null = df[mask]

    if len(null) < min_null_genes:
        mask = (
            (np.abs(df[static_col] - sta) < expand_window)
            & (df["symbol"] != gene_symbol)
            & df[target_col].notna()
        )
        null = df[mask]

    if len(null) < absolute_min:
        return np.nan, len(null)

    null_mean = null[target_col].mean()
    null_std = null[target_col].std()
    if null_std == 0:
        return np.nan, len(null)

    z = (r[target_col] - null_mean) / null_std
    return z, len(null)


def pathway_matched_null_test(
    df: pd.DataFrame,
    gene_list: Sequence[str],
    target_col: str = "residual_z",
    static_col: str = "raw_phi_static",
    window: float = 0.002,
    min_genes: int = 3,
) -> Dict:
    """Test whether a gene set is upstream-shifted via one-sample t-test on matched-null z.

    For each gene in the set, compute its matched-null z. Then perform a
    one-sample t-test against H₀: mean z = 0.

    A positive mean z with small p indicates upstream-specific enrichment
    beyond what static topology would predict.

    Used for Section 4.2 pathway hierarchy:
      - Translation/ribosome: 10/10 CTs p < 0.05
      - Chaperone: 4/10 CTs
      - ATP synthase: 3/10 CTs
      - etc.
    """
    z_values = []
    for g in gene_list:
        z, n = matched_null_z_static_window(df, g, target_col, static_col, window)
        if not np.isnan(z):
            z_values.append(z)

    if len(z_values) < min_genes:
        return {
            "n_genes_found": len(z_values),
            "mean_z": np.nan,
            "median_z": np.nan,
            "p_one_sample_t": np.nan,
            "fraction_above_pos2": np.nan,
            "fraction_below_neg2": np.nan,
        }

    t_stat, p_val = ttest_1samp(z_values, 0)
    return {
        "n_genes_found": len(z_values),
        "mean_z": float(np.mean(z_values)),
        "median_z": float(np.median(z_values)),
        "p_one_sample_t": float(p_val),
        "fraction_above_pos2": float(np.mean(np.array(z_values) > 2)),
        "fraction_below_neg2": float(np.mean(np.array(z_values) < -2)),
    }


def bonferroni_correction(
    z_values: Dict[str, float],
    n_tests: int,
) -> Dict[str, Dict[str, float]]:
    """Per-CT Bonferroni correction across the ~8,874 protein-coding genes.

    Under N(0,1) assumption, two-sided p-value is converted to Bonferroni-corrected p.

    For NEMF L3_L5: z = -6.95 → raw p = 3.7 × 10⁻¹² → Bonferroni p_adj = 3.3 × 10⁻⁸
    (Section 8.9.1, Appendix AV.3)
    """
    result = {}
    for ct, z in z_values.items():
        if np.isnan(z):
            continue
        raw_p = 2 * norm.sf(np.abs(z))  # two-sided
        bonf_p = min(1.0, raw_p * n_tests)
        result[ct] = {
            "z": float(z),
            "raw_p": float(raw_p),
            "bonferroni_p": float(bonf_p),
            "significant_genome_wide": bool(bonf_p < 0.05),
        }
    return result


def binomial_cross_ct_test(
    n_cts_below_threshold: int,
    n_cts_total: int,
    per_ct_p: float = 0.0228,
    n_genome_tests: int = 8874,
) -> Dict:
    """Binomial test for "z < -2 in ≥K cell types" criterion (Appendix AV.4).

    Under the null hypothesis that z-scores are independent N(0,1) across CTs:
      - P(z < -2 in a single CT) ≈ 0.0228
      - P(z < -2 in ≥K / N CTs) = Σ_{k=K}^{N} binom(N, k) · p^k · (1-p)^(N-k)
      - Expected false positives genome-wide = that probability × n_genome_tests

    For NEMF (7/10 CTs): P = 3.62×10⁻¹⁰ → expected FP = 3.2×10⁻⁶ → FDR < 10⁻⁵

    Returns
    -------
    dict:
        p_binomial : probability of hitting threshold by chance
        expected_fp : expected false positives across genome
        fdr : effective FDR (expected FP / 1 discovery)
    """
    p_by_chance = 1 - binom.cdf(n_cts_below_threshold - 1, n_cts_total, per_ct_p)
    expected_fp = p_by_chance * n_genome_tests
    return {
        "p_binomial": float(p_by_chance),
        "expected_false_positives": float(expected_fp),
        "fdr": float(expected_fp),  # FP / 1 discovery (NEMF is the sole hit)
    }


def cross_ct_correlation_with_target(
    zscore_matrix: pd.DataFrame,
    target_gene: str = "NEMF",
    z_cols: Optional[List[str]] = None,
    min_n_cts: int = 5,
) -> pd.DataFrame:
    """Compute cross-CT Spearman ρ between each gene's z-vector and the target's.

    For Figure 6B (Section 8.9.5):
      - Positive ρ = co-collapse axis (LSM14A ρ=+0.94, DROSHA ρ=+0.89, FXR1 ρ=+0.89)
      - Negative ρ = compensatory axis (EIF4EBP2 ρ=-0.99, MRPL39 ρ=-0.89, COX6B1 ρ=-0.86)

    Parameters
    ----------
    zscore_matrix : output of `allgene_insertion.build_zscore_matrix`
    target_gene : reference gene (default NEMF)
    z_cols : list of z_disease_<ct> columns to use (default: all)
    min_n_cts : minimum number of non-NaN CTs required per gene

    Returns
    -------
    DataFrame:
        gene, symbol, cross_ct_rho_<target>, cross_ct_p_<target>,
        n_cts_common (with target), axis_class ('co-collapse', 'compensatory', 'neutral')
    """
    if z_cols is None:
        z_cols = [c for c in zscore_matrix.columns if c.startswith("z_disease_")]

    target_row = zscore_matrix[zscore_matrix["symbol"] == target_gene]
    if len(target_row) == 0:
        raise ValueError(f"Target gene {target_gene} not found in z-score matrix")

    target_z_vec = target_row.iloc[0][z_cols].values.astype(float)

    results = []
    for _, row in zscore_matrix.iterrows():
        gene_z_vec = row[z_cols].values.astype(float)
        valid = ~np.isnan(gene_z_vec) & ~np.isnan(target_z_vec)
        if valid.sum() < min_n_cts:
            continue
        rho, p_val = spearmanr(gene_z_vec[valid], target_z_vec[valid])

        if rho > 0.3 and p_val < 0.1:
            axis = "co-collapse"
        elif rho < -0.3 and p_val < 0.1:
            axis = "compensatory"
        else:
            axis = "neutral"

        results.append({
            "gene": row["gene"],
            "symbol": row["symbol"],
            f"cross_ct_rho_{target_gene}": float(rho),
            f"cross_ct_p_{target_gene}": float(p_val),
            "n_cts_common": int(valid.sum()),
            "axis_class": axis,
        })

    return pd.DataFrame(results)


def identify_universal_downshift_genes(
    zscore_matrix: pd.DataFrame,
    z_cols: Optional[List[str]] = None,
    z_threshold: float = -2.0,
    min_ct_threshold: int = 7,
    direction_consistent: bool = True,
) -> pd.DataFrame:
    """Identify genes satisfying 'z < threshold in ≥K of 10 CTs'.

    NEMF screen criterion (Section 8.9.1):
      - z < -2 in ≥7/10 CTs
      - (Optional) direction-consistent (all 10 CTs have z < 0)

    Under the paper's configuration (threshold=-2, min_ct=7), NEMF is the
    sole gene in the genome satisfying this criterion.
    """
    if z_cols is None:
        z_cols = [c for c in zscore_matrix.columns if c.startswith("z_disease_")]

    z_matrix = zscore_matrix[z_cols].values
    n_below = (z_matrix < z_threshold).sum(axis=1)

    hit_mask = n_below >= min_ct_threshold
    if direction_consistent:
        n_ct_tested = (~np.isnan(z_matrix)).sum(axis=1)
        n_negative = (z_matrix < 0).sum(axis=1)
        hit_mask &= (n_negative == n_ct_tested)

    hits = zscore_matrix[hit_mask].copy()
    hits[f"n_z_below_{z_threshold}"] = n_below[hit_mask]
    hits = hits.sort_values("mean_z" if "mean_z" in hits.columns else f"n_z_below_{z_threshold}")
    return hits


# -------------------------------------------------------------------------
# CLI entry point
# -------------------------------------------------------------------------
def main():
    """Execute the NEMF genome-wide screen on precomputed z-score matrix."""
    import argparse

    parser = argparse.ArgumentParser(description="NEMF screen (Section 8.9)")
    parser.add_argument(
        "--zscore-matrix", type=Path,
        default=Path("results/allgene_3phi/verification/zscore_matrix_wide.csv"),
        help="Path to z-score matrix CSV (wide format: one row per gene, one col per CT)"
    )
    parser.add_argument("--target-gene", type=str, default="NEMF")
    parser.add_argument("--z-threshold", type=float, default=-2.0)
    parser.add_argument("--min-ct", type=int, default=7)
    parser.add_argument("--output-dir", type=Path, default=Path("results/nemf_screen"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading z-score matrix: {args.zscore_matrix}")
    zmat = pd.read_csv(args.zscore_matrix)

    # 1. Universal downshift genes
    logger.info(f"Identifying genes with z < {args.z_threshold} in ≥{args.min_ct}/N CTs...")
    hits = identify_universal_downshift_genes(
        zmat, z_threshold=args.z_threshold, min_ct_threshold=args.min_ct,
    )
    hits.to_csv(args.output_dir / "universal_downshift_hits.csv", index=False)
    logger.info(f"Found {len(hits)} genes meeting criterion; top hit: {hits['symbol'].iloc[0] if len(hits) else 'none'}")

    # 2. Per-CT Bonferroni for target gene
    z_cols = [c for c in zmat.columns if c.startswith("z_disease_")]
    target_row = zmat[zmat["symbol"] == args.target_gene]
    if len(target_row) > 0:
        target_z_dict = {
            c.replace("z_disease_", ""): target_row.iloc[0][c]
            for c in z_cols
        }
        n_tests = len(zmat)
        bonf = bonferroni_correction(target_z_dict, n_tests)
        pd.DataFrame(bonf).T.to_csv(args.output_dir / f"{args.target_gene}_bonferroni.csv")
        logger.info(f"{args.target_gene} Bonferroni: see {args.output_dir}/{args.target_gene}_bonferroni.csv")

    # 3. Binomial null
    n_cts = len(z_cols)
    binom_result = binomial_cross_ct_test(args.min_ct, n_cts, n_genome_tests=len(zmat))
    with open(args.output_dir / "binomial_null.json", "w") as f:
        json.dump(binom_result, f, indent=2)
    logger.info(
        f"Binomial null: P = {binom_result['p_binomial']:.2e}, "
        f"expected FP = {binom_result['expected_false_positives']:.2e}"
    )

    # 4. Cross-CT correlation with target
    logger.info(f"Computing cross-CT ρ vs {args.target_gene}...")
    cross_ct = cross_ct_correlation_with_target(zmat, target_gene=args.target_gene)
    cross_ct.to_csv(args.output_dir / f"cross_ct_rho_vs_{args.target_gene}.csv", index=False)

    logger.info(f"\nSaved NEMF screen results to {args.output_dir}")


if __name__ == "__main__":
    main()
