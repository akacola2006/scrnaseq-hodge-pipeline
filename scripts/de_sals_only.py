"""
Differential expression (SALS-only DESeq2) and GSEA-style rank test (Section 6.16).

Implements the DE pipeline used throughout Section 4.2 and Section 8.9 for
per-gene expression comparisons between disease (SALS) and pathological normal
(PN) donors, with sex as a covariate.

This is the "SALS-only" version that excludes C9ALS donors so that the DE
results correspond to the same cohort used in the main Hodge pipeline.
A parallel ALS-combined analysis (SALS + C9ALS versus PN) can be run as a
sensitivity check by passing disease_labels=["SALS", "C9ALS"].

Pipeline:
  1. Per-CT pseudobulk: sum raw UMI counts across cells within each (donor, CT)
  2. Filter: ≥ 20 cells/donor/CT, ≥ 1 count in ≥ 3 donors
  3. PyDESeq2 Wald test on design `~ condition + sex`
  4. Apply apeglm shrinkage for log₂FC
  5. BH-FDR for multiple testing

Gene class analyses (Section 4.2):
  - RPL+RPS (cytoplasmic ribosomal proteins)
  - Translation_full (RPL+RPS+EIF+EEF)
  - Class mean log₂FC per CT
  - One-sided Mann-Whitney U rank test vs all other genes (GSEA-style)
  - Top-decile enrichment (fraction of class in top 10 %)

Expected outputs match paper Section 4.2 Table:
    | CT     | RPL+RPS mean log₂FC | Translation log₂FC | Rank-test p   |
    | Oligo  | -0.06               | 0.00               | 0.75          |
    | L4/6   | +0.24               | +0.28              | 6×10⁻¹⁴       |
    | PV     | +0.27               | +0.27              | 2×10⁻¹⁵       |
    | Pooled | +0.14               | +0.18              | 5.6×10⁻²¹     |

Migrated from:
    sals_analysis_frozen_20260211/scripts/run_de_sals_only.py
    sals_analysis_frozen_20260211/scripts/run_de_gsea_rank_test.py
"""

from __future__ import annotations

import gc
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import mannwhitneyu

logger = logging.getLogger(__name__)

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    HAS_PYDESEQ2 = True
except ImportError:
    HAS_PYDESEQ2 = False


# Gene class definitions (Section 6.16)
# Uses HGNC gene-symbol prefix matching after DESeq2 low-expression filtering.
GENE_CLASS_PREFIXES = {
    "RPL_RPS": ["RPL", "RPS"],
    "Translation_full": ["RPL", "RPS", "EIF", "EEF"],
    "EIF": ["EIF"],
    "EEF": ["EEF"],
}


def build_pseudobulk(
    adata,
    donor_ids: np.ndarray,
    min_cells_per_donor: int = 20,
    min_counts_per_gene: int = 1,
    min_donors_per_gene: int = 3,
) -> tuple[pd.DataFrame, pd.Series]:
    """Build pseudobulk count matrix (donors × genes) from single-cell AnnData.

    Returns
    -------
    counts_df : DataFrame of summed raw UMI counts (donors × genes)
    donors_kept : Series of donor IDs after filtering
    """
    unique_donors, counts_per_donor = np.unique(donor_ids, return_counts=True)
    donors_kept = unique_donors[counts_per_donor >= min_cells_per_donor]
    logger.info(f"  Donors with ≥ {min_cells_per_donor} cells: {len(donors_kept)}")

    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X
    # AnnData X should be raw counts; verify
    if X.max() < 50 and X.dtype.kind == "f":
        logger.warning("  adata.X appears to be log-normalised; pseudobulk needs raw counts!")

    pseudobulk_rows = []
    for d in donors_kept:
        mask = donor_ids == d
        summed = X[mask].sum(axis=0)
        pseudobulk_rows.append(summed)

    counts = pd.DataFrame(
        np.array(pseudobulk_rows),
        index=donors_kept,
        columns=adata.var_names,
    )

    # Filter genes: ≥ min_counts in ≥ min_donors
    expressed = (counts >= min_counts_per_gene).sum(axis=0) >= min_donors_per_gene
    counts = counts.loc[:, expressed]
    logger.info(f"  Genes after low-expression filter: {counts.shape[1]}")

    return counts.astype(int), pd.Series(donors_kept)


def run_deseq2_sals_vs_pn(
    counts_df: pd.DataFrame,
    donor_condition: Dict[str, str],
    donor_sex: Dict[str, str],
    disease_labels: List[str] = None,
    control_label: str = "PN",
    min_donors_per_group: int = 4,
) -> Optional[pd.DataFrame]:
    """Run PyDESeq2 Wald test with ~condition + sex design.

    Parameters
    ----------
    counts_df : donors × genes raw count matrix
    donor_condition : dict donor_id → condition label
    donor_sex : dict donor_id → 'M' or 'F'
    disease_labels : list of disease group labels (default ["SALS"])
    control_label : reference group label

    Returns
    -------
    DataFrame with columns:
        gene, log2FoldChange, pvalue, padj, baseMean
        (apeglm-shrunken LFC applied)
    """
    if not HAS_PYDESEQ2:
        raise RuntimeError(
            "PyDESeq2 not installed. Run: pip install pydeseq2==0.5.0"
        )

    if disease_labels is None:
        disease_labels = ["SALS"]

    # Filter donors to disease_labels ∪ {control_label}
    kept_donors = [
        d for d in counts_df.index
        if donor_condition.get(str(d)) in disease_labels + [control_label]
    ]
    counts_df = counts_df.loc[kept_donors]

    # Collapse multi-disease labels to "ALS" (SALS+C9ALS sensitivity) or "SALS"
    def _condition_label(d):
        cond = donor_condition.get(str(d))
        if cond in disease_labels:
            return "DISEASE"
        if cond == control_label:
            return "CONTROL"
        return None

    metadata = pd.DataFrame({
        "donor_id": counts_df.index,
        "condition": [_condition_label(d) for d in counts_df.index],
        "sex": [donor_sex.get(str(d), "M") for d in counts_df.index],
    }).set_index("donor_id")

    # Keep only donors with valid condition
    metadata = metadata[metadata["condition"].notna()]
    counts_df = counts_df.loc[metadata.index]

    n_disease = (metadata["condition"] == "DISEASE").sum()
    n_control = (metadata["condition"] == "CONTROL").sum()
    if n_disease < min_donors_per_group or n_control < min_donors_per_group:
        logger.warning(
            f"  Insufficient donors: {n_disease} DISEASE, {n_control} CONTROL "
            f"(min: {min_donors_per_group})"
        )
        return None

    # DESeq2
    try:
        inf = DefaultInference(n_cpus=8)
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design_factors=["condition", "sex"],
            refit_cooks=True,
            inference=inf,
        )
        dds.deseq2()
        stats = DeseqStats(
            dds,
            contrast=["condition", "DISEASE", "CONTROL"],
            inference=inf,
        )
        stats.summary()
        stats.lfc_shrink(coeff="condition_DISEASE_vs_CONTROL")

        result = stats.results_df.reset_index().rename(columns={"index": "gene"})
        # BH-FDR already applied by DESeq2; apeglm shrinkage in log2FoldChange
        return result
    except Exception as e:
        logger.error(f"  DESeq2 failed: {e}")
        return None


def classify_genes_by_prefix(
    gene_symbols: pd.Series,
    class_prefixes: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, pd.Index]:
    """Return dict of gene-class name → Index of rows matching class prefix."""
    class_prefixes = class_prefixes or GENE_CLASS_PREFIXES
    class_indices = {}
    for class_name, prefixes in class_prefixes.items():
        pattern = "|".join([f"^{p}" for p in prefixes])
        mask = gene_symbols.str.match(pattern, na=False)
        class_indices[class_name] = gene_symbols.index[mask]
    return class_indices


def gsea_rank_test(
    de_df: pd.DataFrame,
    gene_class_indices: Dict[str, pd.Index],
    lfc_col: str = "log2FoldChange",
    alternative: str = "greater",
) -> Dict[str, Dict]:
    """Pineda 2024-style GSEA rank test via one-sided Mann-Whitney U.

    For each gene class, test whether its log₂FC rank distribution is
    systematically shifted upward (alternative="greater") vs all other genes.

    Matches paper Section 4.2:
      - L4_L6 Translation rank-test p = 6 × 10⁻¹⁴
      - PV Translation rank-test p = 2 × 10⁻¹⁵
      - Pooled neurons p = 5.6 × 10⁻²¹

    Returns
    -------
    dict: {class_name: {n_class, n_other, mean_lfc, median_lfc, mwu_u, mwu_p,
                        top_decile_fraction, top_decile_enrichment_fold}}
    """
    result = {}
    for class_name, idx in gene_class_indices.items():
        class_lfc = de_df.loc[idx, lfc_col].dropna().values
        other_lfc = de_df.loc[~de_df.index.isin(idx), lfc_col].dropna().values

        if len(class_lfc) < 3 or len(other_lfc) < 3:
            continue

        u_stat, p_val = mannwhitneyu(class_lfc, other_lfc, alternative=alternative)

        # Top-decile enrichment
        all_lfc = np.concatenate([class_lfc, other_lfc])
        top10_threshold = np.percentile(all_lfc, 90)
        class_in_top10 = (class_lfc > top10_threshold).mean()

        result[class_name] = {
            "n_class": int(len(class_lfc)),
            "n_other": int(len(other_lfc)),
            "mean_lfc": float(class_lfc.mean()),
            "median_lfc": float(np.median(class_lfc)),
            "mwu_u": float(u_stat),
            "mwu_p": float(p_val),
            "top_decile_fraction": float(class_in_top10),
            "top_decile_enrichment_fold": float(class_in_top10 / 0.10),
        }
    return result


def run_pooled_neuron_analysis(
    per_ct_de: Dict[str, pd.DataFrame],
    neuronal_cts: List[str] = None,
) -> Dict:
    """Pool neuronal cell types (L2/3, L4/6, L5/6, PV) for pan-neuronal test.

    Concatenates per-CT DE tables, then runs the GSEA rank test on the pooled
    set. Matches Pineda 2024 Fig. 2C methodology.

    Section 4.2 expects:
      - n = 543 Translation genes
      - Pooled rank-test p = 5.6 × 10⁻²¹
      - Top-decile = 21.4 % (vs expected 10 %)
    """
    if neuronal_cts is None:
        neuronal_cts = ["L2_L3", "L4_L6", "L5_L6", "PV"]

    pooled_rows = []
    for ct in neuronal_cts:
        if ct not in per_ct_de:
            continue
        df = per_ct_de[ct].copy()
        df["ct"] = ct
        pooled_rows.append(df)

    if not pooled_rows:
        return {}

    pooled = pd.concat(pooled_rows, ignore_index=True)
    pooled_gene_class = classify_genes_by_prefix(pooled["gene"])
    rank_test = gsea_rank_test(pooled, pooled_gene_class)

    return {
        "pooled_cts": neuronal_cts,
        "n_total_rows": len(pooled),
        "rank_test": rank_test,
    }


# -------------------------------------------------------------------------
# CLI entry point
# -------------------------------------------------------------------------
def main():
    """CLI for SALS-only DE + GSEA rank test across all cell types."""
    import argparse

    parser = argparse.ArgumentParser(description="SALS-only DESeq2 + GSEA rank test (Section 6.16)")
    parser.add_argument(
        "--cell-types", nargs="+",
        default=["Oligo", "OPC", "Astro", "L2_L3", "L4_L6", "L5_L6", "PV"],
    )
    parser.add_argument("--include-c9als", action="store_true",
                        help="Sensitivity analysis: SALS + C9ALS vs PN")
    parser.add_argument("--min-cells-per-donor", type=int, default=20)
    parser.add_argument("--output-dir", type=Path, default=Path("results/de_sals_only"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    from scripts.data_loader import load_all_samples, get_cells_for_celltype, get_sample_info

    sample_info = get_sample_info()
    donor_condition = dict(zip(sample_info["donor_id"].astype(str), sample_info["condition"]))
    donor_sex = (
        dict(zip(sample_info["donor_id"].astype(str), sample_info["sex"]))
        if "sex" in sample_info.columns
        else {}
    )

    disease_labels = ["SALS", "C9ALS"] if args.include_c9als else ["SALS"]
    logger.info(f"Disease labels: {disease_labels}")

    logger.info("Loading all samples...")
    all_adata = load_all_samples()

    per_ct_de = {}
    class_summary_rows = []
    for ct in args.cell_types:
        logger.info(f"\n=== {ct} ===")
        try:
            adata = get_cells_for_celltype(all_adata, ct, min_cells_per_donor=args.min_cells_per_donor)
            donor_ids = adata.obs["donor_id"].values.astype(str)

            pseudobulk, donors_kept = build_pseudobulk(
                adata, donor_ids, min_cells_per_donor=args.min_cells_per_donor,
            )
            del adata
            gc.collect()

            de_df = run_deseq2_sals_vs_pn(
                pseudobulk, donor_condition, donor_sex,
                disease_labels=disease_labels,
            )
            if de_df is None:
                continue

            de_df.to_csv(args.output_dir / f"{ct.lower()}_de_sals.csv", index=False)
            per_ct_de[ct] = de_df

            # Gene class analysis
            gene_classes = classify_genes_by_prefix(de_df["gene"])
            rank_test = gsea_rank_test(de_df, gene_classes)

            for class_name, stats in rank_test.items():
                class_summary_rows.append({
                    "ct": ct,
                    "class": class_name,
                    **stats,
                })
            logger.info(f"  Translation_full n={rank_test.get('Translation_full', {}).get('n_class', 0)}, "
                        f"mean LFC={rank_test.get('Translation_full', {}).get('mean_lfc', 0):.3f}, "
                        f"rank p={rank_test.get('Translation_full', {}).get('mwu_p', np.nan):.2e}")

        except Exception as e:
            logger.error(f"{ct} FAILED: {e}")
            import traceback
            traceback.print_exc()

    # Per-class summary
    if class_summary_rows:
        class_df = pd.DataFrame(class_summary_rows)
        class_df.to_csv(args.output_dir / "class_summary.csv", index=False)

    # Pooled neuron analysis
    pooled = run_pooled_neuron_analysis(per_ct_de)
    with open(args.output_dir / "pooled_neuron_analysis.json", "w") as f:
        json.dump(pooled, f, indent=2, default=str)

    logger.info(f"\nSaved DE + GSEA results to {args.output_dir}")
    if pooled and "rank_test" in pooled:
        trans = pooled["rank_test"].get("Translation_full", {})
        if trans:
            logger.info(
                f"Pooled neurons Translation: n={trans['n_class']}, "
                f"mean LFC={trans['mean_lfc']:.3f}, "
                f"rank p={trans['mwu_p']:.2e}, "
                f"top-10% enrichment={trans['top_decile_enrichment_fold']:.2f}x"
            )


if __name__ == "__main__":
    main()
