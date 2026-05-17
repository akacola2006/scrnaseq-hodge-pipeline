# Centralised path resolution for Mendelian-randomization benchmark scripts.
#
# Each path is read from an environment variable (highest priority) or falls
# back to the repository-relative default under data/external/<subdir>/. None
# of the underlying datasets are redistributed by this repository; see
# benchmark/mendelian_randomization/README.md for the original sources.
#
# Environment variables (override any default):
#   ALS_BRYOIS_EQTL_DIR    -> Bryois et al. brain-region cis-eQTL directory
#   ALS_GWAS_FILE          -> van Rheenen 2021 ALS GWAS summary stats file
#   ALS_FROZEN_DIR         -> Frozen SALS analysis directory (stable-High, dark-grey lists)
#   ALS_MR_RESULTS_DIR     -> Base directory for MR result output
#   ALS_GLIOMA_DIR         -> Glioma working directory (for cross-disease scripts)

get_env_path <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}

DATA_DIR        <- get_env_path("ALS_BRYOIS_EQTL_DIR", "data/external/bryois_eqtl")
ALS_GWAS_FILE   <- get_env_path("ALS_GWAS_FILE",
                                "data/external/als_gwas/GCST90027164_buildGRCh37.tsv.gz")
FROZEN_DIR      <- get_env_path("ALS_FROZEN_DIR",
                                "data/external/sals_analysis_frozen_20260211")
STABLE_HIGH_FILE <- file.path(FROZEN_DIR,
                              "results/track_b/sals_upstream_gene_list/stable_high_genes.csv")
DARKGREY_FILE    <- file.path(FROZEN_DIR,
                              "results/track_b/laneB/darkgrey_identity/darkgrey_genes.tsv")
MR_RESULTS_BASE  <- get_env_path("ALS_MR_RESULTS_DIR", "data/external/mr_results")
