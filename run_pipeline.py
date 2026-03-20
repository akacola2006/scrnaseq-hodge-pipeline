"""
scRNAseq Hodge Decomposition Pipeline — Main Runner
=====================================================
Orchestrates the full analysis pipeline from h5ad files to gene-level
upstream identification via discrete Hodge decomposition.

Usage:
    python run_pipeline.py                    # Run full pipeline
    python run_pipeline.py --step residuals   # Run specific step
    python run_pipeline.py --step lane_a      # Run Lane A only
    python run_pipeline.py --list-steps       # List available steps
"""
import argparse
import json
import logging
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from scripts import config
from scripts.seed_utils import set_global_seed
from scripts.log_utils import FailureLogger

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("pipeline")

STEPS = [
    "validate",
    "load_data",
    "residuals",
    "pca",
    "spd",
    "pseudotime",
    "lane_a",
    "bootstrap",
    "lane_b",
    "gene_hodge",
    "enrichment",
]


def validate_setup():
    """Validate configuration and data availability."""
    logger.info("=" * 60)
    logger.info("Validating pipeline setup...")
    logger.info("=" * 60)

    issues = config.validate_config()
    if issues:
        for issue in issues:
            logger.error("CONFIG ERROR: %s", issue)
        raise RuntimeError(
            "Configuration validation failed. Please fix the issues above "
            "and edit project_config.yaml accordingly."
        )

    # Check h5ad files exist
    h5ad_files = list(config.H5AD_DIR.glob("*.h5ad"))
    logger.info("Found %d h5ad files in %s", len(h5ad_files), config.H5AD_DIR)
    if not h5ad_files:
        raise RuntimeError(
            f"No h5ad files found in {config.H5AD_DIR}. "
            "Please place your scRNAseq h5ad files there."
        )

    # Check sample info
    sample_info = pd.read_csv(config.SAMPLE_INFO_PATH)
    logger.info("Sample info: %d samples, columns: %s",
                len(sample_info), list(sample_info.columns))

    # Check conditions
    if config.CONDITION_COLUMN in sample_info.columns:
        conditions = sample_info[config.CONDITION_COLUMN].value_counts()
        logger.info("Conditions: %s", conditions.to_dict())
    else:
        logger.warning("Condition column '%s' not found in sample_info.csv",
                       config.CONDITION_COLUMN)

    logger.info("Cell types to analyze: %s", config.CELLTYPE_LIST)
    logger.info("Pseudotime cell types: %s", config.PT_CELLTYPE_SET)
    logger.info("Control label: %s", config.CONTROL_LABEL)
    logger.info("GPU enabled: %s", config.USE_GPU)

    # Check GPU
    try:
        import torch
        if config.USE_GPU and torch.cuda.is_available():
            logger.info("GPU: %s (%.1f GB)",
                       torch.cuda.get_device_name(0),
                       torch.cuda.get_device_properties(0).total_mem / 1e9)
        elif config.USE_GPU:
            logger.warning("GPU requested but CUDA not available. Will use CPU.")
    except ImportError:
        if config.USE_GPU:
            logger.warning("PyTorch not installed. Will use CPU.")

    logger.info("Validation PASSED")
    return sample_info


def run_full_pipeline(start_step: str = None):
    """Run the full analysis pipeline."""
    t0 = time.time()

    # Setup
    set_global_seed(config.GLOBAL_SEED)
    run_id = config.make_run_id()
    dirs = config.setup_run_dirs(run_id)
    failure_logger = FailureLogger(dirs["logs"])

    logger.info("=" * 60)
    logger.info("scRNAseq Hodge Decomposition Pipeline")
    logger.info("Run ID: %s", run_id)
    logger.info("Output: %s", dirs["root"])
    logger.info("=" * 60)

    # Determine which steps to run
    if start_step:
        if start_step not in STEPS:
            logger.error("Unknown step: %s. Available: %s", start_step, STEPS)
            sys.exit(1)
        step_idx = STEPS.index(start_step)
        steps_to_run = STEPS[step_idx:]
        logger.info("Starting from step: %s", start_step)
    else:
        steps_to_run = STEPS

    # ── Step 1: Validate ──────────────────────────────────────
    if "validate" in steps_to_run:
        sample_info = validate_setup()
    else:
        sample_info = pd.read_csv(config.SAMPLE_INFO_PATH)
        sample_info["donor_id"] = sample_info["donor_id"].astype(str)

    # ── Step 2: Load Data ─────────────────────────────────────
    if "load_data" in steps_to_run:
        logger.info("\n[STEP 2/10] Loading h5ad data...")
        from scripts.data_loader import load_all_samples
        adata_dict = load_all_samples()
        logger.info("Loaded %d samples", len(adata_dict))
    else:
        adata_dict = None  # Will need to reload if needed

    # ── Steps 3-5: Residuals, PCA, SPD per celltype ──────────
    celltypes = config.CELLTYPE_LIST
    all_residuals = {}
    all_gene_masks = {}
    all_gene_names = {}
    all_loadings = {}
    spd_matrices = {}

    need_processing = any(s in steps_to_run for s in ["residuals", "pca", "spd"])

    if need_processing:
        if adata_dict is None:
            from scripts.data_loader import load_all_samples
            adata_dict = load_all_samples()

        from scripts.data_loader import get_cells_for_celltype
        from scripts.residuals import compute_residuals
        from scripts.pca_engine import pca_gpu, extract_donor_scores
        from scripts.spd import estimate_spd

        for ct in celltypes:
            logger.info("\n[PROCESSING] Cell type: %s", ct)

            try:
                # Get cells
                adata_ct = get_cells_for_celltype(adata_dict, ct)
                logger.info("  %s: %d cells, %d donors",
                           ct, adata_ct.n_obs,
                           adata_ct.obs["donor_id"].nunique())

                # Residuals
                if "residuals" in steps_to_run:
                    logger.info("  [Step 3] Computing residuals...")
                    resid, gene_mask, gene_names = compute_residuals(
                        adata_ct, hvg_n=config.HVG_N
                    )
                    all_residuals[ct] = resid
                    all_gene_masks[ct] = gene_mask
                    all_gene_names[ct] = gene_names
                    logger.info("  Residuals: %s", resid.shape)

                # PCA
                if "pca" in steps_to_run and ct in all_residuals:
                    logger.info("  [Step 4] Running PCA (k=%d)...", config.PCA_K)
                    scores, loadings, svals = pca_gpu(
                        all_residuals[ct], k=config.PCA_K,
                        device=config.GPU_DEVICE,
                    )
                    all_loadings[ct] = loadings
                    donor_scores = extract_donor_scores(
                        scores, adata_ct.obs["donor_id"].values
                    )
                    logger.info("  PCA: %d components, %d donors",
                               scores.shape[1], len(donor_scores))

                    # SPD estimation
                    if "spd" in steps_to_run:
                        logger.info("  [Step 5] Estimating SPD covariance...")
                        spd_matrices[ct] = {}
                        for donor_id, ds in donor_scores.items():
                            if ds.shape[0] >= 2:
                                spd_matrices[ct][donor_id] = estimate_spd(ds)
                        logger.info("  SPD: %d donors", len(spd_matrices[ct]))

            except Exception as e:
                logger.error("  FAILED for %s: %s", ct, e)
                failure_logger.log("processing", e, celltype=ct)
                continue

    # ── Step 6: Pseudotime ────────────────────────────────────
    pt_df = None
    if "pseudotime" in steps_to_run and spd_matrices:
        logger.info("\n[STEP 6/10] Building pseudotime (PT-B)...")
        from scripts.pseudotime import build_pt_b, assign_windows

        try:
            pt_df, pt_metadata = build_pt_b(
                spd_matrices, sample_info,
                celltype_set=config.PT_CELLTYPE_SET,
                k=config.PCA_K,
            )
            pt_df = assign_windows(pt_df, n_windows=config.N_WINDOWS)

            # Save PT
            pt_df.to_csv(dirs["pt"] / "pseudotime_windows.csv", index=False)
            logger.info("Pseudotime saved: %d donors, %d windows",
                       len(pt_df), pt_df["window"].nunique())
        except Exception as e:
            logger.error("Pseudotime FAILED: %s", e)
            failure_logger.log("pseudotime", e)

    # ── Step 7: Lane A ────────────────────────────────────────
    lane_a_result = None
    if "lane_a" in steps_to_run and pt_df is not None and spd_matrices:
        logger.info("\n[STEP 7/10] Running Lane A (cell-type upstream analysis)...")
        from scripts.lane_a import run_lane_a

        try:
            lane_a_result = run_lane_a(
                spd_matrices, pt_df, celltypes,
                output_dir=dirs["laneA"],
            )
            logger.info("Lane A: upstream ranking = %s",
                       lane_a_result.get("ranked_celltypes", []))
        except Exception as e:
            logger.error("Lane A FAILED: %s", e)
            failure_logger.log("lane_a", e)

    # ── Step 8: Bootstrap ─────────────────────────────────────
    bootstrap_result = None
    if "bootstrap" in steps_to_run and pt_df is not None and spd_matrices:
        logger.info("\n[STEP 8/10] Running bootstrap (B=%d)...", config.BOOTSTRAP_N)
        from scripts.bootstrap import bootstrap_lane_a

        try:
            bootstrap_result = bootstrap_lane_a(
                spd_matrices=spd_matrices,
                pt_df=pt_df,
                celltypes=celltypes,
                output_dir=dirs["bootstrap"],
                failure_logger=failure_logger,
            )
            logger.info("Bootstrap: top-3 reproducibility = %s",
                       {ct: f"{v:.2f}" for ct, v in
                        sorted(bootstrap_result["top3_reproducibility"].items(),
                               key=lambda x: -x[1])[:5]})
        except Exception as e:
            logger.error("Bootstrap FAILED: %s", e)
            failure_logger.log("bootstrap", e)

    # ── Step 9: Lane B ────────────────────────────────────────
    if "lane_b" in steps_to_run and lane_a_result and bootstrap_result:
        logger.info("\n[STEP 9/10] Running Lane B (gene-level analysis)...")
        from scripts.lane_b import identify_upstream_pcs, extract_gene_set
        from scripts.gene_hodge import check_execution_gate

        # Identify upstream celltype
        ranked = lane_a_result.get("ranked_celltypes", [])
        if not ranked:
            logger.warning("No ranked celltypes from Lane A")
        else:
            upstream_ct = ranked[0]
            logger.info("Upstream celltype: %s", upstream_ct)

            # Check execution gate
            gate_pass, gate_msg = check_execution_gate(bootstrap_result, upstream_ct)
            logger.info("Gate check: %s", gate_msg)

            if gate_pass and upstream_ct in all_loadings:
                # Build correlation matrices by window for Lane B
                from scripts.spd import cov_to_corr
                if upstream_ct in spd_matrices:
                    donor_to_window = dict(zip(
                        pt_df["donor_id"].astype(str),
                        pt_df["window"].astype(int),
                    ))

                    # Compute window-level correlation matrices from SPD
                    from scripts.spd import log_euclidean_mean
                    k = config.PCA_K
                    corr_by_window = {}
                    for w in range(config.N_WINDOWS):
                        window_mats = [
                            sigma for did, sigma in spd_matrices[upstream_ct].items()
                            if donor_to_window.get(str(did)) == w
                        ]
                        if len(window_mats) >= config.MIN_DONORS_PER_WINDOW:
                            mean_cov = log_euclidean_mean(window_mats)
                            corr, _ = cov_to_corr(mean_cov)
                            corr_by_window[w] = corr

                    # Find upstream window (max d_corr transition)
                    d_corr_data = lane_a_result.get("distance_decomp", {}).get(upstream_ct, {})
                    d_corrs = d_corr_data.get("d_corr", [])
                    w_pairs = d_corr_data.get("window_pairs", [])

                    if d_corrs:
                        max_idx = np.argmax(d_corrs)
                        upstream_window = w_pairs[max_idx][0]
                    else:
                        upstream_window = 0

                    try:
                        lb_result = identify_upstream_pcs(
                            corr_by_window, upstream_window,
                        )
                        gene_set, gene_indices = extract_gene_set(
                            all_loadings[upstream_ct],
                            lb_result["pc_indices"],
                            all_gene_names.get(upstream_ct, []),
                        )

                        # Save Lane B results
                        lb_summary = {
                            "upstream_celltype": upstream_ct,
                            "upstream_window": int(upstream_window),
                            "n_upstream_pcs": len(lb_result["pc_indices"]),
                            "pc_indices": lb_result["pc_indices"].tolist(),
                            "n_genes": len(gene_set),
                        }
                        with open(dirs["laneB"] / "lane_b_summary.json", "w") as f:
                            json.dump(lb_summary, f, indent=2)

                        # Save gene set
                        pd.DataFrame({
                            "gene": gene_set,
                            "index": gene_indices.tolist(),
                        }).to_csv(dirs["laneB"] / "upstream_gene_set.csv", index=False)

                        logger.info("Lane B: %d upstream PCs, %d genes",
                                   len(lb_result["pc_indices"]), len(gene_set))

                    except Exception as e:
                        logger.error("Lane B FAILED: %s", e)
                        failure_logger.log("lane_b", e, celltype=upstream_ct)

    # ── Step 10: Gene Hodge ───────────────────────────────────
    if "gene_hodge" in steps_to_run:
        gene_set_path = dirs["laneB"] / "upstream_gene_set.csv"
        if gene_set_path.exists() and lane_a_result:
            logger.info("\n[STEP 10/10] Running gene-level Hodge decomposition...")
            from scripts.gene_hodge import (
                resolve_gene_indices, precompute_donor_log_corr, run_gene_hodge,
            )

            ranked = lane_a_result.get("ranked_celltypes", [])
            if ranked:
                upstream_ct = ranked[0]
                gene_set_df = pd.read_csv(gene_set_path)
                gene_set = gene_set_df["gene"].tolist()

                if upstream_ct in all_residuals and upstream_ct in all_gene_names:
                    gene_indices, resolved_genes = resolve_gene_indices(
                        gene_set, all_gene_names[upstream_ct],
                    )

                    if len(resolved_genes) >= config.GENE_HODGE_MIN_GENES:
                        # Need adata for donor_ids
                        from scripts.data_loader import get_cells_for_celltype
                        if adata_dict is None:
                            from scripts.data_loader import load_all_samples
                            adata_dict = load_all_samples()
                        adata_ct = get_cells_for_celltype(adata_dict, upstream_ct)

                        donor_log_corr, donor_window = precompute_donor_log_corr(
                            all_residuals[upstream_ct],
                            adata_ct.obs["donor_id"].values.astype(str),
                            gene_indices,
                            pt_df,
                        )

                        # Find upstream window
                        d_corr_data = lane_a_result.get("distance_decomp", {}).get(upstream_ct, {})
                        d_corrs = d_corr_data.get("d_corr", [])
                        w_pairs = d_corr_data.get("window_pairs", [])
                        if d_corrs:
                            max_idx = np.argmax(d_corrs)
                            w_star = w_pairs[max_idx][0]
                        else:
                            w_star = 0

                        try:
                            gh_result = run_gene_hodge(
                                donor_log_corr, donor_window,
                                list(donor_log_corr.keys()),
                                w_star, resolved_genes,
                            )

                            if gh_result["status"] == "OK":
                                # Save results
                                gene_results = pd.DataFrame({
                                    "gene": gh_result["gene_names"],
                                    "phi": gh_result["phi"],
                                    "classification": gh_result["classification"],
                                })
                                gene_results = gene_results.sort_values("phi", ascending=False)
                                gene_results.to_csv(
                                    dirs["gene_hodge"] / "gene_phi_scores.csv",
                                    index=False,
                                )

                                gh_summary = {
                                    "upstream_celltype": upstream_ct,
                                    "n_genes": gh_result["n_genes"],
                                    "gradient_fraction": gh_result["gradient_fraction"],
                                    "p_value": gh_result["p_value"],
                                    "n_high": gh_result["classification"].count("High"),
                                    "n_medium": gh_result["classification"].count("Medium"),
                                    "n_low": gh_result["classification"].count("Low"),
                                    "flow_mode": gh_result["flow_mode"],
                                    "transition": gh_result["transition"],
                                }
                                with open(dirs["gene_hodge"] / "gene_hodge_summary.json", "w") as f:
                                    json.dump(gh_summary, f, indent=2, default=str)

                                logger.info(
                                    "Gene Hodge: %d genes, GF=%.4f, p=%.4f, "
                                    "%d High / %d Medium / %d Low",
                                    gh_result["n_genes"],
                                    gh_result["gradient_fraction"],
                                    gh_result["p_value"],
                                    gh_summary["n_high"],
                                    gh_summary["n_medium"],
                                    gh_summary["n_low"],
                                )
                            else:
                                logger.warning("Gene Hodge skipped: %s",
                                             gh_result.get("reason", "unknown"))

                        except Exception as e:
                            logger.error("Gene Hodge FAILED: %s", e)
                            failure_logger.log("gene_hodge", e, celltype=upstream_ct)

    # ── Step 11: Enrichment Analysis ─────────────────────────
    if "enrichment" in steps_to_run:
        gene_scores_path = dirs["gene_hodge"] / "gene_phi_scores.csv"
        if gene_scores_path.exists():
            logger.info("\n[STEP 11] Running functional enrichment analysis...")
            from scripts.enrichment import (
                run_enrichment_analysis, annotate_gene_hodge_results,
                load_functional_modules,
            )

            try:
                modules = load_functional_modules()
                if modules:
                    gene_scores_df = pd.read_csv(gene_scores_path)

                    # Build gh_result-like dict for enrichment
                    gh_result_for_enrich = {
                        "status": "OK",
                        "gene_names": gene_scores_df["gene"].tolist(),
                        "classification": gene_scores_df["classification"].tolist(),
                    }

                    # Background = all genes in analysis
                    bg_genes = []
                    for ct_genes in all_gene_names.values():
                        bg_genes.extend(ct_genes)
                    bg_genes = list(set(bg_genes)) if bg_genes else gh_result_for_enrich["gene_names"]

                    enrichment_results = run_enrichment_analysis(
                        gh_result_for_enrich, bg_genes,
                        output_dir=dirs["gene_hodge"] / "enrichment",
                    )

                    # Annotate gene results
                    annotated = annotate_gene_hodge_results(gene_scores_df)
                    annotated.to_csv(
                        dirs["gene_hodge"] / "gene_phi_scores_annotated.csv",
                        index=False,
                    )
                    logger.info("Enrichment and annotation complete")

                    for test_name, enrich_df in enrichment_results.items():
                        sig = enrich_df[enrich_df["significant"]]
                        if len(sig) > 0:
                            logger.info("  %s: %d significant modules — %s",
                                       test_name, len(sig),
                                       sig[["module", "fdr"]].head(5).to_dict("records"))

            except Exception as e:
                logger.error("Enrichment FAILED: %s", e)
                failure_logger.log("enrichment", e)

    # ── Summary ───────────────────────────────────────────────
    elapsed = time.time() - t0
    logger.info("\n" + "=" * 60)
    logger.info("Pipeline complete in %.1f seconds", elapsed)
    logger.info("Results: %s", dirs["root"])
    logger.info("=" * 60)

    # Save run summary
    summary = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "elapsed_seconds": elapsed,
        "project_name": config.PROJECT_NAME,
        "n_celltypes": len(celltypes),
        "celltypes": celltypes,
        "n_samples": len(pd.read_csv(config.SAMPLE_INFO_PATH)),
        "pca_k": config.PCA_K,
        "n_windows": config.N_WINDOWS,
        "bootstrap_n": config.BOOTSTRAP_N,
    }
    if lane_a_result and lane_a_result.get("status") == "OK":
        summary["upstream_ranking"] = lane_a_result.get("ranked_celltypes", [])
        summary["gradient_fraction"] = lane_a_result.get("gradient_fraction", 0)
    if bootstrap_result:
        summary["top3_reproducibility"] = bootstrap_result.get("top3_reproducibility", {})

    with open(dirs["root"] / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)


def main():
    parser = argparse.ArgumentParser(
        description="scRNAseq Hodge Decomposition Pipeline"
    )
    parser.add_argument(
        "--step", type=str, default=None,
        help="Start from a specific step (validate, load_data, residuals, "
             "pca, spd, pseudotime, lane_a, bootstrap, lane_b, gene_hodge)"
    )
    parser.add_argument(
        "--list-steps", action="store_true",
        help="List available pipeline steps"
    )
    args = parser.parse_args()

    if args.list_steps:
        print("Available pipeline steps:")
        for i, step in enumerate(STEPS, 1):
            print(f"  {i:2d}. {step}")
        return

    run_full_pipeline(start_step=args.step)


if __name__ == "__main__":
    main()
