"""
P4 Validation — Run Track 2: Gene-Level Hodge on Perturbation Data
===================================================================
Applies IDS-SALS Hodge decomposition to CRISPR perturbation data.
For each target gene's KO, computes gene-gene correlation dynamics
(Control → KO) and tests whether the target gene is ranked as
"upstream" (high φ) by the Hodge method.

Ultra-low-memory implementation: processes one library at a time.

Usage (on the analysis machine):
  python p4_scripts/run_track2.py

Required files:
  - GSE274058 scRNA-seq data (gse274058_extract/ with 5 libraries)

Outputs:
  - p4_results/track2_results.json         (full results)
  - p4_results/track2_per_perturbation/    (per-perturbation Hodge results)
  - P4_TRACK2_RESULTS_REPORT.md            (human-readable report)
"""
import logging
import json
import time
import gc
import gzip
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger(__name__)

# ── Adjust these paths for your machine ──
# On Cowork VM:
DATA_DIR = Path("/sessions/great-affectionate-newton/gse274058_extract")
OUTPUT_DIR = Path("/sessions/great-affectionate-newton/p4_validation/results")

# On analysis machine (uncomment and adjust):
# DATA_DIR = Path(r"D:\Projects\scRNAseq\p4_validation\gse274058_extract")
# OUTPUT_DIR = Path(r"D:\Projects\scRNAseq\p4_validation\results")

SAMPLES = [
    "GSM8442771_B1_1", "GSM8442772_B1_2", "GSM8442773_B2_1",
    "GSM8442774_B2_2", "GSM8442775_B3_1",
]

PERTURB = {
    "mSafe_gRNA": ("mSafe", "Ctrl(neg)"),
    "Trem2_gRNA": ("Trem2", "AD"),
    "Rraga_gRNA": ("Rraga", "Control"),
    "Fasn_gRNA": ("Fasn", "Control"),
    "Clu_gRNA": ("Clu", "AD"),
    "Dpp6_gRNA": ("Dpp6", "ALS"),
    "Tbk1_gRNA": ("Tbk1", "ALS"),
    "Flcn_gRNA": ("Flcn", "Control"),
    "Gfap_gRNA": ("Gfap", "Control"),
    "C9orf72_gRNA": ("C9orf72", "ALS"),
    "Cfap410_gRNA": ("Cfap410", "ALS"),
    "Stk39_gRNA": ("Stk39", "PD"),
    "Lrrk2_gRNA": ("Lrrk2", "PD"),
    "Ndufaf_gRNA": ("Ndufaf2", "AD"),
    "Sh3gl2_gRNA": ("Sh3gl2", "PD"),
    "Srf_gRNA": ("Srf", "Control"),
    "Rbfox_gRNA": ("Rbfox1", "AD"),
    "Olig2_gRNA": ("Olig2", "Control"),
}

CATEGORY_MAP = {v[0]: v[1] for v in PERTURB.values()}

# ── Track 2 Parameters ──
N_HVG = 300           # Gene set size for Hodge (smaller than Track 3's 5000)
                       # Ledoit-Wolf needs n_cells > n_genes / 3 for good estimates
                       # Smallest KO group: 12 cells → max ~36 genes for ideal LW
                       # With LW shrinkage, 300 works even for small n (n >= 15)
N_HVG_POOL = 5000     # Initial HVG pool (same as Track 3)
SEED = 42
MAX_CTRL_PER_LIB = 2000
FLOW_MODE = "edge_weight"  # Matching IDS-SALS Oligo analysis
N_BOOTSTRAP = 50      # Cell-level bootstrap iterations (reduced for speed)
N_PERM = 500           # Permutation test iterations
MIN_KO_CELLS = 15     # Minimum KO cells for Hodge analysis
HIGH_SIGMA = 1.5       # φ classification threshold (matching IDS-SALS)


def main():
    import sys
    scripts_dir = str(Path(__file__).parent)
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)

    from track2_gene_hodge import (
        compute_perturbation_hodge,
        bootstrap_perturbation_hodge,
        permutation_test,
        evaluate_target_gene,
        aggregate_track2_results,
        hodge_gradient_kn,
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    per_pert_dir = OUTPUT_DIR / "track2_per_perturbation"
    per_pert_dir.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)

    t0 = time.time()

    # ════════════════════════════════════════════════════════
    # Phase 1: Load gene names & identify HVGs
    # ════════════════════════════════════════════════════════
    logger.info("Phase 1: Loading gene names and computing HVGs...")

    mdir = DATA_DIR / SAMPLES[0]
    with gzip.open(mdir / "features.tsv.gz", "rt") as f:
        features = [l.strip().split("\t") for l in f]
    gene_names = np.array([ft[1] for ft in features])
    n_genes_total = len(gene_names)
    logger.info(f"  Total features: {n_genes_total}")

    # Identify gRNA indices
    grna_idx_map = {}
    for gn in PERTURB:
        idx = np.where(gene_names == gn)[0]
        if len(idx) == 1:
            grna_idx_map[gn] = idx[0]

    # Endogenous gene mask (exclude gRNA features)
    grna_set = set(grna_idx_map.values())
    myrf = np.where(gene_names == "Myrf_gRNA")[0]
    if len(myrf) > 0:
        grna_set.add(myrf[0])
    endo_mask = np.array([i not in grna_set for i in range(n_genes_total)])
    endo_idx = np.where(endo_mask)[0]

    # HVG from first library subsample
    mat = scipy.io.mmread(mdir / "matrix.mtx.gz").tocsc()
    nc = mat.shape[1]
    sub_cells = rng.choice(nc, size=min(2000, nc), replace=False)
    X_sub = mat[:, sub_cells].T.toarray().astype(np.float64)
    X_sub = X_sub[:, endo_mask]
    totals = X_sub.sum(axis=1, keepdims=True)
    totals = np.where(totals > 0, totals, 1.0)
    X_sub = np.log1p(X_sub / totals * 1e6)
    gene_var = X_sub.var(axis=0)

    # Get top N_HVG_POOL variance genes, then select top N_HVG for Hodge
    hvg_pool_local = np.argsort(gene_var)[-N_HVG_POOL:]
    hvg_pool_local = np.sort(hvg_pool_local)

    # For Track 2, use only top N_HVG genes
    # But first check which perturbation targets are in the HVG set
    hvg_local = np.argsort(gene_var)[-N_HVG:]
    hvg_local = np.sort(hvg_local)
    hvg_global = endo_idx[hvg_local]
    hvg_gene_names = gene_names[hvg_global]

    # Ensure perturbation target genes are included in the gene set
    target_genes_mouse = [v[0] for v in PERTURB.values()]
    targets_in_hvg = [g for g in target_genes_mouse if g in hvg_gene_names]
    targets_not_in_hvg = [g for g in target_genes_mouse if g not in hvg_gene_names]

    if targets_not_in_hvg:
        logger.info(f"  Adding {len(targets_not_in_hvg)} target genes not in top-{N_HVG} HVGs")
        for tg in targets_not_in_hvg:
            tg_global_idx = np.where(gene_names == tg)[0]
            if len(tg_global_idx) == 1 and endo_mask[tg_global_idx[0]]:
                hvg_global = np.append(hvg_global, tg_global_idx[0])
        hvg_global = np.unique(np.sort(hvg_global))
        hvg_gene_names = gene_names[hvg_global]

    logger.info(f"  HVG gene set: {len(hvg_global)} genes (top-{N_HVG} + targets)")
    logger.info(f"  Targets in set: {len([g for g in target_genes_mouse if g in hvg_gene_names])}/{len(set(target_genes_mouse))}")

    del mat, X_sub, gene_var
    gc.collect()

    # ════════════════════════════════════════════════════════
    # Phase 2: Load scRNA-seq (library-by-library)
    # ════════════════════════════════════════════════════════
    logger.info("Phase 2: Loading scRNA-seq data...")

    perturb_cells = {v[0]: [] for v in PERTURB.values()}
    ctrl_cells = []

    for lib in SAMPLES:
        logger.info(f"  Loading {lib}...")
        mdir = DATA_DIR / lib
        mat = scipy.io.mmread(mdir / "matrix.mtx.gz").tocsc()
        nc = mat.shape[1]
        mat_csr = mat.T.tocsr()

        # QC
        genes_pc = np.diff(mat_csr.indptr)
        mito_m = np.array([g.startswith("mt-") for g in gene_names])
        if mito_m.any():
            mito_c = np.asarray(mat_csr[:, mito_m].sum(axis=1)).ravel()
            tot_c = np.asarray(mat_csr.sum(axis=1)).ravel()
            pctm = np.where(tot_c > 0, mito_c / tot_c, 0)
        else:
            pctm = np.zeros(nc)
        qc = (genes_pc >= 200) & (genes_pc <= 7000) & (pctm <= 0.15)

        # Classify by barcode
        n_bc = np.zeros(nc, dtype=int)
        grna_positive = {}
        for gn, gidx in grna_idx_map.items():
            col = np.asarray(mat[gidx, :].todense()).ravel()
            pos = col > 0
            grna_positive[gn] = pos
            n_bc[pos] += 1

        # Collect perturbation cells (HVG expression, log-normalised)
        for gn, gidx in grna_idx_map.items():
            short = PERTURB[gn][0]
            single = grna_positive[gn] & (n_bc == 1) & qc
            if single.sum() > 0:
                idx = np.where(single)[0]
                X_small = mat_csr[idx][:, hvg_global].toarray().astype(np.float64)
                t = X_small.sum(axis=1, keepdims=True)
                t = np.where(t > 0, t, 1.0)
                X_small = np.log1p(X_small / t * 1e6).astype(np.float32)
                perturb_cells[short].append(X_small)

        # Collect control cells (barcode-negative)
        neg = (n_bc == 0) & qc
        neg_idx = np.where(neg)[0]
        if len(neg_idx) > MAX_CTRL_PER_LIB:
            neg_idx = rng.choice(neg_idx, MAX_CTRL_PER_LIB, replace=False)
        X_ctrl = mat_csr[neg_idx][:, hvg_global].toarray().astype(np.float64)
        t = X_ctrl.sum(axis=1, keepdims=True)
        t = np.where(t > 0, t, 1.0)
        X_ctrl = np.log1p(X_ctrl / t * 1e6).astype(np.float32)
        ctrl_cells.append(X_ctrl)

        del mat, mat_csr, grna_positive
        gc.collect()
        logger.info(f"    {lib}: {qc.sum()} QC pass, {neg.sum()} neg")

    X_ctrl_all = np.vstack(ctrl_cells).astype(np.float64)
    del ctrl_cells
    gc.collect()
    logger.info(f"  Control pool: {X_ctrl_all.shape[0]} cells × {X_ctrl_all.shape[1]} genes")

    # ════════════════════════════════════════════════════════
    # Phase 3: Hodge Decomposition per Perturbation
    # ════════════════════════════════════════════════════════
    logger.info("Phase 3: Running Hodge decomposition per perturbation...")

    gene_name_list = list(hvg_gene_names)
    per_perturbation = {}
    target_evaluations = {}

    CTRL_RATIO = 8  # Use up to 8× KO count as controls

    for gene in sorted(perturb_cells.keys()):
        pieces = perturb_cells[gene]
        if not pieces:
            logger.info(f"  {gene}: no cells, skipping")
            per_perturbation[gene] = {
                "error": "no cells", "n_ko": 0,
                "category": CATEGORY_MAP.get(gene, "Unknown"),
            }
            continue

        X_ko = np.vstack(pieces).astype(np.float64)
        n_ko = X_ko.shape[0]
        category = CATEGORY_MAP.get(gene, "Unknown")

        if n_ko < MIN_KO_CELLS:
            logger.info(f"  {gene}: only {n_ko} KO cells < {MIN_KO_CELLS}, skipping")
            per_perturbation[gene] = {
                "error": f"too few cells ({n_ko})",
                "n_ko": n_ko,
                "category": category,
            }
            continue

        # Matched control subset
        n_ctrl_use = min(n_ko * CTRL_RATIO, X_ctrl_all.shape[0])
        ctrl_idx = rng.choice(X_ctrl_all.shape[0], size=n_ctrl_use, replace=False)
        X_ctrl_matched = X_ctrl_all[ctrl_idx]

        logger.info(f"  {gene} ({category}): {n_ko} KO, {n_ctrl_use} Ctrl, {len(gene_name_list)} genes")

        # ── Point estimate ──
        result = compute_perturbation_hodge(
            X_ko, X_ctrl_matched, gene_name_list,
            flow_mode=FLOW_MODE,
        )

        if result is None:
            logger.warning(f"  {gene}: Hodge computation failed")
            per_perturbation[gene] = {
                "error": "Hodge failed",
                "n_ko": n_ko,
                "category": category,
            }
            continue

        # ── Permutation test ──
        logger.info(f"    Permutation test ({N_PERM})...")
        try:
            perm_result = permutation_test(
                result["flow"], len(gene_name_list),
                n_perm=N_PERM, seed=SEED,
            )
        except Exception as exc:
            logger.warning(f"    Permutation test failed: {exc}")
            perm_result = {"p_value": None, "observed_gradient_fraction": None}

        # ── Bootstrap ──
        logger.info(f"    Bootstrap ({N_BOOTSTRAP})...")
        boot_result = bootstrap_perturbation_hodge(
            X_ko, X_ctrl_matched, gene_name_list,
            n_bootstrap=N_BOOTSTRAP,
            flow_mode=FLOW_MODE,
            seed=SEED + hash(gene) % 1000,
        )

        # ── Target gene evaluation ──
        target_eval = evaluate_target_gene(
            result["phi"], gene_name_list, gene,
            result["classification"],
            bootstrap_result=boot_result,
        )

        # Store results
        per_perturbation[gene] = {
            "target_gene": gene,
            "category": category,
            "n_ko": n_ko,
            "n_ctrl": n_ctrl_use,
            "n_genes": result["n_genes"],
            "gradient_fraction": result["gradient_fraction"],
            "non_gradient_fraction": result["non_gradient_fraction"],
            "r_phi": result["r_phi"],
            "flow_mode": result["flow_mode"],
            "permutation_p": perm_result.get("p_value"),
            "permutation_null_mean": perm_result.get("null_mean"),
            "target_phi": target_eval.get("phi"),
            "target_rank": target_eval.get("rank"),
            "target_percentile": target_eval.get("percentile"),
            "target_z_score": target_eval.get("z_score"),
            "target_tier": target_eval.get("tier"),
            "target_in_gene_set": target_eval.get("in_gene_set"),
            "bootstrap_n_success": boot_result.get("n_success", 0),
            "bootstrap_n_stable_high": boot_result.get("n_stable_high", 0),
            "bootstrap_gf_ci": boot_result.get("gradient_fraction_ci"),
            "classification_stats": result["classification"]["stats"],
        }

        target_evaluations[gene] = target_eval

        # Save per-perturbation detail
        detail = {
            **per_perturbation[gene],
            "gene_phi_top20": [
                {"gene": e["gene"], "phi": e["phi"], "rank": e["rank"]}
                for e in result["classification"]["high"][:20]
            ],
            "target_bootstrap_ci": boot_result.get("phi_ci", {}).get(gene),
            "target_high_repro": boot_result.get("high_reproducibility", {}).get(gene),
        }
        with open(per_pert_dir / f"{gene}_hodge.json", "w") as f:
            json.dump(detail, f, indent=2, default=str)

        pct_str = f"{target_eval['percentile']:.1f}" if 'percentile' in target_eval else "N/A"
        logger.info(
            f"    → gf={result['gradient_fraction']:.4f} "
            f"r_phi={result['r_phi']:.4f} "
            f"perm_p={perm_result.get('p_value', 'N/A')} | "
            f"target: rank={target_eval.get('rank', 'N/A')}/{result['n_genes']} "
            f"pct={pct_str} "
            f"tier={target_eval.get('tier', 'N/A')}"
        )

        del X_ko, X_ctrl_matched
        gc.collect()

    # ════════════════════════════════════════════════════════
    # Phase 4: Aggregate Analysis
    # ════════════════════════════════════════════════════════
    logger.info("Phase 4: Aggregating results...")

    aggregate = aggregate_track2_results(per_perturbation, target_evaluations)

    # ════════════════════════════════════════════════════════
    # Phase 5: Save Results
    # ════════════════════════════════════════════════════════
    logger.info("Phase 5: Saving results...")

    full_results = {
        "track": "Track 2",
        "description": "Gene-level Hodge decomposition on perturbation data",
        "parameters": {
            "n_hvg": N_HVG,
            "n_hvg_pool": N_HVG_POOL,
            "flow_mode": FLOW_MODE,
            "n_bootstrap": N_BOOTSTRAP,
            "n_perm": N_PERM,
            "min_ko_cells": MIN_KO_CELLS,
            "high_sigma": HIGH_SIGMA,
            "ctrl_ratio": CTRL_RATIO,
        },
        "gene_set_size": len(gene_name_list),
        "gene_names": gene_name_list,
        "per_perturbation": per_perturbation,
        "target_evaluations": {
            g: {k: v for k, v in e.items()}
            for g, e in target_evaluations.items()
        },
        "aggregate": aggregate,
    }

    with open(OUTPUT_DIR / "track2_results.json", "w") as f:
        json.dump(full_results, f, indent=2, default=str)
    logger.info(f"  Saved: {OUTPUT_DIR / 'track2_results.json'}")

    # ════════════════════════════════════════════════════════
    # Phase 6: Generate Report
    # ════════════════════════════════════════════════════════
    logger.info("Phase 6: Generating report...")
    generate_report(full_results, per_perturbation, target_evaluations, aggregate)

    elapsed = time.time() - t0
    logger.info(f"Track 2 complete in {elapsed:.0f}s")


def generate_report(full_results, per_perturbation, target_evaluations, aggregate):
    """Generate human-readable markdown report."""
    report_path = OUTPUT_DIR.parent / "P4_TRACK2_RESULTS_REPORT.md"

    lines = []
    lines.append("# P4 Track 2: Gene-Level Hodge on Perturbation Data — Results Report")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append("**Question**: When the IDS-SALS Hodge decomposition method is applied to")
    lines.append("CRISPR perturbation data (Control → KO transcriptomic variation), does the")
    lines.append("KO target gene rank as \"upstream\" (high φ)?")
    lines.append("")

    agg = aggregate
    mean_pct = agg.get("mean_percentile", 0)
    n_eval = agg.get("n_targets_evaluated", 0)
    wp = agg.get("wilcoxon_p")
    bp = agg.get("binom_p")

    if wp is not None and wp < 0.05:
        answer = "**Yes** — target genes are significantly enriched in the upstream (high-φ) positions."
    elif mean_pct > 60:
        answer = "**Trend** — target genes show a tendency toward upstream positions, but not statistically significant."
    else:
        answer = "**No** — target genes do not systematically rank as upstream in the Hodge hierarchy."

    lines.append(f"**Answer**: {answer}")
    lines.append("")
    lines.append(f"- Mean target percentile: {mean_pct:.1f}% (50% = random expectation)")
    lines.append(f"- Wilcoxon signed-rank p: {wp}")
    lines.append(f"- Binomial test p: {bp}")
    tiers = agg.get("tier_distribution", {})
    lines.append(f"- Target tier distribution: {tiers.get('High',0)} High / {tiers.get('Medium',0)} Medium / {tiers.get('Low',0)} Low")
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## 1. Per-Perturbation Results")
    lines.append("")
    lines.append("| Gene | Category | n_KO | n_genes | grad_frac | perm_p | target_rank | target_pct | tier |")
    lines.append("|---|---|---|---|---|---|---|---|---|")

    for gene in sorted(per_perturbation.keys()):
        p = per_perturbation[gene]
        if "error" in p:
            lines.append(f"| {gene} | {p.get('category','')} | {p.get('n_ko',0)} | — | — | — | — | — | SKIP |")
            continue
        tpct = p.get('target_percentile')
        tpct_s = f"{tpct:.1f}%" if isinstance(tpct, (int, float)) else "N/A"
        tgf = p.get('gradient_fraction')
        tgf_s = f"{tgf:.4f}" if isinstance(tgf, (int, float)) else "N/A"
        lines.append(
            f"| {gene} | {p.get('category','')} "
            f"| {p.get('n_ko',0)} "
            f"| {p.get('n_genes',0)} "
            f"| {tgf_s} "
            f"| {p.get('permutation_p','N/A')} "
            f"| {p.get('target_rank','N/A')}/{p.get('n_genes',0)} "
            f"| {tpct_s} "
            f"| {p.get('target_tier','N/A')} |"
        )

    lines.append("")
    lines.append("## 2. Category Summary")
    lines.append("")
    cat_stats = agg.get("category_stats", {})
    if cat_stats:
        lines.append("| Category | n | mean pct | std pct |")
        lines.append("|---|---|---|---|")
        for cat in ["ALS", "PD", "AD", "Control", "Ctrl(neg)"]:
            if cat in cat_stats:
                cs = cat_stats[cat]
                lines.append(f"| {cat} | {cs['n']} | {cs['mean_percentile']:.1f}% | {cs['std_percentile']:.1f} |")

    lines.append("")
    lines.append("## 3. Gradient Fraction Summary")
    lines.append("")
    gf = agg.get("gradient_fraction_summary", {})
    if gf:
        lines.append(f"- Mean gradient fraction: {gf.get('mean', 0):.4f} ± {gf.get('std', 0):.4f}")
        lines.append(f"- Range: [{gf.get('min', 0):.4f}, {gf.get('max', 0):.4f}]")
        lines.append("")
        lines.append("For reference, the IDS-SALS Oligo edge_weight analysis had")
        lines.append("gradient_fraction = 0.398, r_phi = 0.629.")

    lines.append("")
    lines.append("## 4. Interpretation")
    lines.append("")
    lines.append("### 4.1 What Track 2 Tests")
    lines.append("")
    lines.append("Unlike Track 3 (which compared *human φ values* with *mouse perturbation")
    lines.append("effects*), Track 2 applies the Hodge method *de novo* to perturbation data.")
    lines.append("This tests whether the *method itself* — not the specific IDS-SALS φ values")
    lines.append("— captures meaningful perturbation-induced gene hierarchy.")
    lines.append("")
    lines.append("### 4.2 Design Choices")
    lines.append("")
    lines.append("- **Pooled Ledoit-Wolf** (not pseudo-donors): Perturbation data has no")
    lines.append("  natural donor structure. Pooled estimation with cell-bootstrap provides")
    lines.append("  uncertainty quantification.")
    lines.append("- **edge_weight flow mode**: Matches the IDS-SALS Oligo analysis.")
    lines.append(f"- **{full_results.get('gene_set_size', N_HVG)} genes**: Top HVGs + perturbation targets.")
    lines.append("  Smaller than IDS-SALS (922) due to small KO group sizes.")
    lines.append("")
    lines.append("### 4.3 Key Caveats")
    lines.append("")
    lines.append("- Small KO groups (12-75 cells) limit correlation estimation precision")
    lines.append("- Mouse hippocampus mixed cell types ≠ human motor cortex Oligodendrocytes")
    lines.append("- Acute CRISPR knockdown ≠ chronic disease progression")
    lines.append("- The \"upstream\" interpretation differs: in IDS-SALS, φ reflects position")
    lines.append("  along disease pseudotime; here, φ reflects position along KO→Control axis")
    lines.append("")

    params = full_results.get("parameters", {})
    lines.append("## 5. Technical Parameters")
    lines.append("")
    for k, v in params.items():
        lines.append(f"- {k}: {v}")

    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append(f"*Generated: 2026-03-17 | Pipeline: P4 Track 2 | Data: GSE274058*")
    lines.append("")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    logger.info(f"  Report: {report_path}")


if __name__ == "__main__":
    main()
