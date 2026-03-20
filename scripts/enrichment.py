"""
scRNAseq Hodge Pipeline — Functional Module Enrichment
========================================================
Fisher's exact test enrichment of gene Hodge results against
pre-defined functional gene modules (e.g., Synaptic, Oxidative_Stress,
Myelination, etc.).

Also provides gene annotation lookups (ENSEMBL ID <-> gene symbol,
chromosome, neuron/glia classification, TDP-43 probability, etc.).
"""
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from . import config

logger = logging.getLogger(__name__)


# =====================================================================
# 1. Functional Modules
# =====================================================================

_MODULES_CACHE: Optional[Dict[str, List[str]]] = None


def load_functional_modules() -> Dict[str, List[str]]:
    """Load functional gene modules from JSON.

    Returns dict of {module_name: [gene_symbol, ...]}.
    """
    global _MODULES_CACHE
    if _MODULES_CACHE is not None:
        return _MODULES_CACHE

    path = config.FUNCTIONAL_MODULES_PATH
    if path is None or not path.exists():
        logger.warning("Functional modules file not found: %s", path)
        return {}

    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    # Support both flat format and nested {"modules": {...}} format
    if "modules" in data:
        modules = data["modules"]
    else:
        modules = {k: v for k, v in data.items() if isinstance(v, list)}

    _MODULES_CACHE = modules
    logger.info("Loaded %d functional modules (%s)",
                len(modules),
                {k: len(v) for k, v in list(modules.items())[:5]})
    return modules


def get_module_summary() -> Dict[str, int]:
    """Return module name -> gene count summary."""
    modules = load_functional_modules()
    return {name: len(genes) for name, genes in modules.items()}


# =====================================================================
# 2. Gene Annotation
# =====================================================================

_ANNOTATION_CACHE: Optional[pd.DataFrame] = None


def load_gene_annotation() -> Optional[pd.DataFrame]:
    """Load gene annotation CSV.

    Expected columns: gene_id, gene_symbol, chromosome, gene_type,
    neuron_or_glia, neuron, glia, tdp43_prob, mirna_score, etc.
    """
    global _ANNOTATION_CACHE
    if _ANNOTATION_CACHE is not None:
        return _ANNOTATION_CACHE

    path = config.GENE_ANNOTATION_PATH
    if path is None or not path.exists():
        logger.warning("Gene annotation file not found: %s", path)
        return None

    annot = pd.read_csv(path)
    if "gene_id" in annot.columns:
        annot = annot.set_index("gene_id")

    _ANNOTATION_CACHE = annot
    logger.info("Loaded gene annotation: %d genes, columns: %s",
                len(annot), list(annot.columns))
    return annot


def ensembl_to_symbol(ensembl_ids: List[str]) -> Dict[str, str]:
    """Map ENSEMBL IDs to gene symbols."""
    annot = load_gene_annotation()
    if annot is None or "gene_symbol" not in annot.columns:
        return {}

    result = {}
    for eid in ensembl_ids:
        bare = eid.split(".")[0]
        if bare in annot.index:
            result[eid] = annot.loc[bare, "gene_symbol"]
    return result


def symbol_to_ensembl(symbols: List[str]) -> Dict[str, str]:
    """Map gene symbols to ENSEMBL IDs."""
    annot = load_gene_annotation()
    if annot is None or "gene_symbol" not in annot.columns:
        return {}

    sym_to_ens = {}
    for eid, row in annot.iterrows():
        sym = row.get("gene_symbol", "")
        if sym:
            sym_to_ens[sym] = eid
    return {s: sym_to_ens[s] for s in symbols if s in sym_to_ens}


def get_gene_info(gene_ids: List[str]) -> pd.DataFrame:
    """Get annotation info for a list of genes (ENSEMBL or symbol).

    Returns DataFrame with available annotation columns.
    """
    annot = load_gene_annotation()
    if annot is None:
        return pd.DataFrame({"gene": gene_ids})

    rows = []
    for gid in gene_ids:
        bare = gid.split(".")[0]
        if bare in annot.index:
            row = annot.loc[bare].to_dict()
            row["query"] = gid
            rows.append(row)
        elif "gene_symbol" in annot.columns:
            match = annot[annot["gene_symbol"] == gid]
            if len(match) > 0:
                row = match.iloc[0].to_dict()
                row["query"] = gid
                rows.append(row)
            else:
                rows.append({"query": gid})
        else:
            rows.append({"query": gid})

    return pd.DataFrame(rows)


# =====================================================================
# 3. Fisher Enrichment Test
# =====================================================================

def fisher_enrichment(
    target_genes: List[str],
    background_genes: List[str],
    modules: Optional[Dict[str, List[str]]] = None,
    alpha: float = 0.05,
    use_symbols: bool = True,
) -> pd.DataFrame:
    """Fisher's exact test enrichment of target genes against functional modules.

    Parameters
    ----------
    target_genes : list of str
        Gene set to test (e.g., High-phi genes from gene Hodge).
    background_genes : list of str
        Background gene set (e.g., all genes in the analysis).
    modules : dict, optional
        {module_name: [gene, ...]}. Defaults to loaded functional_modules.json.
    alpha : float
        FDR significance threshold.
    use_symbols : bool
        If True and genes are ENSEMBL IDs, convert to symbols for matching.

    Returns
    -------
    pd.DataFrame
        Columns: module, n_target, n_background, n_module, n_overlap,
        odds_ratio, p_value, fdr, significant, overlap_genes.
        Sorted by p_value.
    """
    if modules is None:
        modules = load_functional_modules()

    if not modules:
        logger.warning("No functional modules available for enrichment")
        return pd.DataFrame()

    # Convert ENSEMBL to symbols if needed
    if use_symbols:
        annot = load_gene_annotation()
        if annot is not None and "gene_symbol" in annot.columns:
            # Check if target genes look like ENSEMBL IDs
            if target_genes and target_genes[0].startswith("ENSG"):
                ens_map = ensembl_to_symbol(target_genes)
                target_genes = [ens_map.get(g, g) for g in target_genes]
                ens_map_bg = ensembl_to_symbol(background_genes)
                background_genes = [ens_map_bg.get(g, g) for g in background_genes]

    target_set = set(target_genes)
    bg_set = set(background_genes)
    N = len(bg_set)

    results = []
    for module_name, module_genes in modules.items():
        module_set = set(module_genes) & bg_set  # Restrict to background
        n_module = len(module_set)

        if n_module == 0:
            continue

        # 2x2 contingency table
        overlap = target_set & module_set
        n_overlap = len(overlap)
        n_target = len(target_set & bg_set)

        a = n_overlap                          # target AND module
        b = n_target - n_overlap               # target AND NOT module
        c = n_module - n_overlap               # NOT target AND module
        d = N - n_target - n_module + n_overlap  # NOT target AND NOT module

        # Fisher's exact test (one-sided: over-representation)
        if a + b == 0 or a + c == 0:
            continue

        table = np.array([[a, b], [c, d]])
        odds_ratio, p_value = fisher_exact(table, alternative="greater")

        results.append({
            "module": module_name,
            "n_target": n_target,
            "n_background": N,
            "n_module": n_module,
            "n_overlap": n_overlap,
            "odds_ratio": odds_ratio,
            "p_value": p_value,
            "overlap_genes": ", ".join(sorted(overlap)) if n_overlap <= 20
                             else f"{', '.join(sorted(overlap)[:20])}... (+{n_overlap - 20})",
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)

    # FDR correction
    reject, fdr_pvals, _, _ = multipletests(df["p_value"].values, alpha=alpha, method="fdr_bh")
    df["fdr"] = fdr_pvals
    df["significant"] = reject

    df = df.sort_values("p_value").reset_index(drop=True)

    n_sig = df["significant"].sum()
    logger.info("Enrichment: %d modules tested, %d significant (FDR < %.2f)",
                len(df), n_sig, alpha)

    return df


# =====================================================================
# 4. Annotate Gene Hodge Results
# =====================================================================

def annotate_gene_hodge_results(
    gene_results: pd.DataFrame,
    add_modules: bool = True,
) -> pd.DataFrame:
    """Annotate gene Hodge results with gene info and module membership.

    Parameters
    ----------
    gene_results : pd.DataFrame
        Must have 'gene' column. Typically also has 'phi', 'classification'.
    add_modules : bool
        If True, add boolean columns for each functional module.

    Returns
    -------
    pd.DataFrame
        Annotated DataFrame with additional columns.
    """
    df = gene_results.copy()

    # Add gene annotation
    annot = load_gene_annotation()
    if annot is not None and "gene_symbol" in annot.columns:
        gene_info = get_gene_info(df["gene"].tolist())

        # Map relevant columns
        for col in ["gene_symbol", "chromosome", "gene_type",
                     "neuron_or_glia", "tdp43_prob", "mirna_score"]:
            if col in gene_info.columns:
                col_map = dict(zip(gene_info["query"], gene_info[col]))
                df[col] = df["gene"].map(col_map)

    # Add module membership
    if add_modules:
        modules = load_functional_modules()
        # Convert gene names to symbols if needed
        genes = df["gene"].tolist()
        if genes and genes[0].startswith("ENSG"):
            ens_map = ensembl_to_symbol(genes)
            gene_symbols = [ens_map.get(g, g) for g in genes]
        else:
            gene_symbols = genes

        for module_name, module_genes in modules.items():
            module_set = set(module_genes)
            df[f"module_{module_name}"] = [
                sym in module_set for sym in gene_symbols
            ]

    return df


def run_enrichment_analysis(
    gene_hodge_result: Dict[str, Any],
    all_gene_names: List[str],
    output_dir: Optional[Path] = None,
) -> Dict[str, pd.DataFrame]:
    """Run full enrichment analysis on gene Hodge results.

    Tests enrichment for:
    - High-tier genes
    - High + Medium-tier genes
    - All genes (as sanity check)

    Returns dict of {test_name: enrichment_DataFrame}.
    """
    if gene_hodge_result.get("status") != "OK":
        logger.warning("Gene Hodge not OK; skipping enrichment")
        return {}

    gene_names = gene_hodge_result["gene_names"]
    classification = gene_hodge_result["classification"]

    results = {}

    # High tier
    high_genes = [g for g, c in zip(gene_names, classification) if c == "High"]
    if high_genes:
        results["high"] = fisher_enrichment(high_genes, all_gene_names)
        logger.info("High-tier enrichment: %d genes", len(high_genes))

    # High + Medium tier
    high_med_genes = [g for g, c in zip(gene_names, classification) if c in ("High", "Medium")]
    if high_med_genes:
        results["high_medium"] = fisher_enrichment(high_med_genes, all_gene_names)

    # Save
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for name, df in results.items():
            df.to_csv(output_dir / f"enrichment_{name}.csv", index=False)

        # Summary
        summary = {}
        for name, df in results.items():
            sig = df[df["significant"]]
            summary[name] = {
                "n_tested": len(df),
                "n_significant": len(sig),
                "top_modules": sig[["module", "odds_ratio", "fdr"]].head(10).to_dict("records")
                    if len(sig) > 0 else [],
            }
        with open(output_dir / "enrichment_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)

        logger.info("Enrichment results saved to %s", output_dir)

    return results
