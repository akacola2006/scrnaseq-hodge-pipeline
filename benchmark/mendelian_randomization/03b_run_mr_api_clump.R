#!/usr/bin/env Rscript
# ============================================================
# 03b_run_mr_api_clump.R
# MR with IEU OpenGWAS API LD clumping (JWT authenticated)
# ============================================================

user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))

# Load JWT from .Renviron
readRenviron(file.path(Sys.getenv("USERPROFILE"), ".Renviron"))

library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(ggplot2)

cat(sprintf("JWT set: %s...\n", substr(Sys.getenv("OPENGWAS_JWT"), 1, 20)))

# --- Paths ---
DATA_DIR      <- "D:/Projects/MR/data/bryois_eqtl"
ALS_GWAS_FILE <- "D:/Projects/MR/data/als_gwas/GCST90027164_buildGRCh37.tsv.gz"
OUTPUT_DIR    <- "D:/Projects/MR/results_clumped"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

STABLE_HIGH_FILE <- "D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/sals_upstream_gene_list/stable_high_genes.csv"
DARKGREY_FILE    <- "D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/laneB/darkgrey_identity/darkgrey_genes.tsv"

# ============================================================
# Reusable functions (same as 02_run_mr.R)
# ============================================================
load_snp_pos <- function() {
  fread(file.path(DATA_DIR, "snp_pos.txt.gz"))
}

load_celltype_eqtl <- function(cell_type, snp_pos, gene_symbols = NULL) {
  cat(sprintf("\n--- Loading %s eQTL ---\n", cell_type))
  all_chr <- lapply(1:22, function(chr) {
    f <- file.path(DATA_DIR, cell_type, sprintf("%s.%d.gz", cell_type, chr))
    if (!file.exists(f)) return(NULL)
    fread(f, header = FALSE)
  })
  eqtl <- rbindlist(all_chr[!sapply(all_chr, is.null)])
  setnames(eqtl, c("gene_id", "snp_id", "dist_tss", "pvalue", "beta"))
  eqtl[, gene_symbol := sub("_ENSG.*", "", gene_id)]

  if (!is.null(gene_symbols))
    eqtl <- eqtl[gene_symbol %in% gene_symbols]

  eqtl_m <- merge(eqtl, snp_pos, by.x = "snp_id", by.y = "SNP")
  eqtl_m[, se := abs(beta) / qnorm(pvalue / 2, lower.tail = FALSE)]
  eqtl_m[!is.finite(se) | se <= 0, se := NA]

  eqtl_sig <- eqtl_m[pvalue < 5e-8 & !is.na(se)]
  if (nrow(eqtl_sig) < 3)
    eqtl_sig <- eqtl_m[pvalue < 5e-6 & !is.na(se)]

  eqtl_sig[, F_stat := (beta / se)^2]
  eqtl_sig <- eqtl_sig[F_stat > 10]
  cat(sprintf("  %d instruments (F>10), %d genes\n",
              nrow(eqtl_sig), uniqueN(eqtl_sig$gene_symbol)))
  return(eqtl_sig)
}

to_exposure_format <- function(eqtl_sig, label) {
  eqtl_sig[, chr_num := as.integer(sub("chr", "", sub(":.*", "", SNP_id_hg19)))]
  eqtl_sig[, pos_hg19 := as.integer(sub(".*:", "", SNP_id_hg19))]

  exposure <- data.frame(
    SNP = eqtl_sig$snp_id, beta.exposure = eqtl_sig$beta,
    se.exposure = eqtl_sig$se, pval.exposure = eqtl_sig$pvalue,
    effect_allele.exposure = eqtl_sig$effect_allele,
    other_allele.exposure = eqtl_sig$other_allele,
    eaf.exposure = eqtl_sig$MAF,
    exposure = label, id.exposure = label,
    gene.exposure = eqtl_sig$gene_symbol,
    chr.exposure = eqtl_sig$chr_num,
    pos.exposure = eqtl_sig$pos_hg19,
    mr_keep.exposure = TRUE, stringsAsFactors = FALSE
  ) %>%
    group_by(SNP) %>%
    slice_min(pval.exposure, n = 1, with_ties = FALSE) %>%
    ungroup() %>% as.data.frame()

  return(exposure)
}

format_als_outcome <- function(als_gwas, snps) {
  als_filt <- als_gwas[rsid %in% snps]
  if (nrow(als_filt) == 0) return(NULL)
  data.frame(
    SNP = als_filt$rsid, beta.outcome = als_filt$beta,
    se.outcome = als_filt$standard_error,
    pval.outcome = als_filt$p_value,
    effect_allele.outcome = toupper(als_filt$effect_allele),
    other_allele.outcome = toupper(als_filt$other_allele),
    eaf.outcome = als_filt$effect_allele_frequency,
    outcome = "ALS", id.outcome = "ALS",
    mr_keep.outcome = TRUE,
    samplesize.outcome = als_filt$N_effective,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Core MR with API clumping
# ============================================================
run_mr <- function(exposure_dat, label, als_gwas) {
  cat(sprintf("\n========== %s ==========\n", label))
  rdir <- file.path(OUTPUT_DIR, label)
  dir.create(rdir, showWarnings = FALSE, recursive = TRUE)

  n0 <- nrow(exposure_dat)
  cat(sprintf("  Pre-clump: %d SNPs\n", n0))

  # API clumping (now with JWT)
  exposure_c <- tryCatch({
    clump_data(exposure_dat, clump_r2 = 0.001, clump_kb = 10000)
  }, error = function(e) {
    cat(sprintf("  Clump error: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(exposure_c) || nrow(exposure_c) == 0) {
    cat("  FAIL: Clumping returned nothing.\n")
    return(NULL)
  }
  cat(sprintf("  Post-clump: %d independent instruments\n", nrow(exposure_c)))

  if (nrow(exposure_c) < 3) {
    cat("  SKIP: <3 instruments.\n")
    return(NULL)
  }

  # Outcome
  outcome <- format_als_outcome(als_gwas, exposure_c$SNP)
  if (is.null(outcome)) { cat("  SKIP: no outcome.\n"); return(NULL) }
  cat(sprintf("  Outcome matched: %d\n", nrow(outcome)))

  # Harmonize
  dat <- harmonise_data(exposure_c, outcome)
  dat <- dat[dat$mr_keep, ]
  cat(sprintf("  Harmonized: %d\n", nrow(dat)))
  if (nrow(dat) < 3) { cat("  SKIP: <3 harmonized.\n"); return(NULL) }

  # MR
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                  "mr_weighted_median", "mr_weighted_mode"))
  pleio <- mr_pleiotropy_test(dat)
  hetero <- mr_heterogeneity(dat)
  single <- mr_singlesnp(dat)
  loo <- mr_leaveoneout(dat)

  cat(sprintf("\n  --- %s (CLUMPED) ---\n", label))
  print(res[, c("method", "nsnp", "b", "se", "pval")])
  if (nrow(pleio) > 0) cat(sprintf("  Egger intercept p = %.4f\n", pleio$pval[1]))
  ivw_h <- hetero[hetero$method == "Inverse variance weighted", ]
  if (nrow(ivw_h) > 0) cat(sprintf("  Q_pval = %.4f\n", ivw_h$Q_pval[1]))

  # Plots
  tryCatch({
    ggsave(file.path(rdir, "scatter.png"), mr_scatter_plot(res, dat)[[1]], width=8, height=6)
    ggsave(file.path(rdir, "forest.png"), mr_forest_plot(single)[[1]], width=8, height=8)
    ggsave(file.path(rdir, "funnel.png"), mr_funnel_plot(single)[[1]], width=8, height=6)
    ggsave(file.path(rdir, "loo.png"), mr_leaveoneout_plot(loo)[[1]], width=8, height=8)
  }, error = function(e) cat(sprintf("  Plot error: %s\n", conditionMessage(e))))

  obj <- list(mr_results = res, dat = dat, pleiotropy = pleio,
              heterogeneity = hetero, single = single, loo = loo,
              n = nrow(dat), label = label)
  saveRDS(obj, file.path(rdir, "result.rds"))
  write.csv(res, file.path(rdir, "mr_results.csv"), row.names = FALSE)
  return(obj)
}

# ============================================================
# Main
# ============================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("MR with API LD Clumping (JWT Authenticated)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

snp_pos <- load_snp_pos()
cat(sprintf("SNP positions: %d\n", nrow(snp_pos)))

als <- fread(ALS_GWAS_FILE)
cat(sprintf("ALS GWAS: %d SNPs\n", nrow(als)))

stable_high <- fread(STABLE_HIGH_FILE)
darkgrey <- fread(DARKGREY_FILE)
dg_col <- intersect(names(darkgrey), c("symbol", "gene", "gene_symbol"))
dg_symbols <- if (length(dg_col) > 0) darkgrey[[dg_col[1]]] else darkgrey[[1]]

rtk_genes <- c("NTRK2", "IGF1R", "INSR", "FYN", "DOCK3", "ROR1", "LRP4",
                "ADAM10", "ADAM12", "MEGF9", "NEO1", "EPS8", "EPS15",
                "RASA1", "RASGRF1", "RASGRF2", "RASGRP3", "MAP2K5",
                "MAPK10", "PRKCA")

R <- list()

# Track A
oligo <- load_celltype_eqtl("Oligodendrocytes", snp_pos)
R[["A_Oligo"]] <- run_mr(to_exposure_format(oligo, "Oligo_eQTL"), "A_Oligo", als)

# Track B
oligo_sh <- load_celltype_eqtl("Oligodendrocytes", snp_pos, stable_high$symbol)
if (nrow(oligo_sh) >= 3)
  R[["B_StableHigh"]] <- run_mr(to_exposure_format(oligo_sh, "StableHigh_eQTL"), "B_StableHigh", als)

# Track D
oligo_dg <- load_celltype_eqtl("Oligodendrocytes", snp_pos, dg_symbols)
if (nrow(oligo_dg) >= 3)
  R[["D_Darkgrey"]] <- run_mr(to_exposure_format(oligo_dg, "Darkgrey_eQTL"), "D_Darkgrey", als)

# Track E
oligo_rtk <- load_celltype_eqtl("Oligodendrocytes", snp_pos, rtk_genes)
if (nrow(oligo_rtk) >= 3)
  R[["E_RTK"]] <- run_mr(to_exposure_format(oligo_rtk, "RTK_eQTL"), "E_RTK", als)

# Summary
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("FINAL SUMMARY (LD CLUMPED)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

sdf <- do.call(rbind, lapply(names(R), function(nm) {
  r <- R[[nm]]
  if (is.null(r)) return(NULL)
  ivw <- r$mr_results[r$mr_results$method == "Inverse variance weighted", ]
  if (nrow(ivw) == 0) return(NULL)
  data.frame(track = nm, n = ivw$nsnp, beta = round(ivw$b, 4),
             se = round(ivw$se, 4), pval = signif(ivw$pval, 3),
             OR = round(exp(ivw$b), 3),
             CI_lo = round(exp(ivw$b - 1.96*ivw$se), 3),
             CI_hi = round(exp(ivw$b + 1.96*ivw$se), 3),
             egger_p = signif(r$pleiotropy$pval[1], 3),
             stringsAsFactors = FALSE)
}))

if (!is.null(sdf)) {
  print(sdf, row.names = FALSE)
  write.csv(sdf, file.path(OUTPUT_DIR, "summary_clumped.csv"), row.names = FALSE)

  p <- ggplot(sdf, aes(x = reorder(track, OR), y = OR)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = "MR: Oligo eQTL -> ALS (LD Clumped, r2<0.001)",
         x = "", y = "Odds Ratio (95% CI)") +
    theme_minimal(base_size = 14)
  ggsave(file.path(OUTPUT_DIR, "forest_clumped.png"), p, width = 10, height = 5)
}

cat("\nDone.\n")
