#!/usr/bin/env Rscript
# ============================================================
# 04_relaxed_threshold_mr.R
# MR with relaxed p-value threshold (p < 5e-6)
# + Power analysis
# ============================================================

user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))
readRenviron(file.path(Sys.getenv("USERPROFILE"), ".Renviron"))

library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(ggplot2)

DATA_DIR      <- "D:/Projects/MR/data/bryois_eqtl"
ALS_GWAS_FILE <- "D:/Projects/MR/data/als_gwas/GCST90027164_buildGRCh37.tsv.gz"
OUTPUT_DIR    <- "D:/Projects/MR/results_relaxed"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

STABLE_HIGH_FILE <- "D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/sals_upstream_gene_list/stable_high_genes.csv"
DARKGREY_FILE    <- "D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/laneB/darkgrey_identity/darkgrey_genes.tsv"

# ============================================================
# Load functions
# ============================================================
load_snp_pos <- function() fread(file.path(DATA_DIR, "snp_pos.txt.gz"))

load_celltype_eqtl <- function(cell_type, snp_pos, gene_symbols = NULL,
                                pval_thresh = 5e-6) {
  cat(sprintf("\n--- %s (p < %.0e) ---\n", cell_type, pval_thresh))
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

  eqtl_sig <- eqtl_m[pvalue < pval_thresh & !is.na(se)]
  eqtl_sig[, F_stat := (beta / se)^2]
  eqtl_sig <- eqtl_sig[F_stat > 10]

  cat(sprintf("  %d instruments, %d genes, median F=%.1f\n",
              nrow(eqtl_sig), uniqueN(eqtl_sig$gene_symbol),
              median(eqtl_sig$F_stat, na.rm = TRUE)))
  return(eqtl_sig)
}

to_exposure <- function(eqtl_sig, label) {
  eqtl_sig[, chr_num := as.integer(sub("chr", "", sub(":.*", "", SNP_id_hg19)))]
  eqtl_sig[, pos_hg19 := as.integer(sub(".*:", "", SNP_id_hg19))]
  data.frame(
    SNP = eqtl_sig$snp_id, beta.exposure = eqtl_sig$beta,
    se.exposure = eqtl_sig$se, pval.exposure = eqtl_sig$pvalue,
    effect_allele.exposure = eqtl_sig$effect_allele,
    other_allele.exposure = eqtl_sig$other_allele,
    eaf.exposure = eqtl_sig$MAF,
    exposure = label, id.exposure = label,
    gene.exposure = eqtl_sig$gene_symbol,
    chr.exposure = eqtl_sig$chr_num, pos.exposure = eqtl_sig$pos_hg19,
    mr_keep.exposure = TRUE, stringsAsFactors = FALSE
  ) %>% group_by(SNP) %>%
    slice_min(pval.exposure, n = 1, with_ties = FALSE) %>%
    ungroup() %>% as.data.frame()
}

format_outcome <- function(als, snps) {
  af <- als[rsid %in% snps]
  if (nrow(af) == 0) return(NULL)
  data.frame(
    SNP = af$rsid, beta.outcome = af$beta,
    se.outcome = af$standard_error, pval.outcome = af$p_value,
    effect_allele.outcome = toupper(af$effect_allele),
    other_allele.outcome = toupper(af$other_allele),
    eaf.outcome = af$effect_allele_frequency,
    outcome = "ALS", id.outcome = "ALS",
    mr_keep.outcome = TRUE, samplesize.outcome = af$N_effective,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# MR with API clumping
# ============================================================
run_mr <- function(exp_dat, label, als) {
  cat(sprintf("\n===== %s =====\n", label))
  rdir <- file.path(OUTPUT_DIR, label)
  dir.create(rdir, showWarnings = FALSE, recursive = TRUE)

  n0 <- nrow(exp_dat)
  exp_c <- tryCatch(
    clump_data(exp_dat, clump_r2 = 0.001, clump_kb = 10000),
    error = function(e) { cat(sprintf("  Clump err: %s\n", conditionMessage(e))); NULL }
  )
  if (is.null(exp_c) || nrow(exp_c) < 3) {
    cat(sprintf("  SKIP: %d -> %d IVs (need >=3)\n", n0,
                ifelse(is.null(exp_c), 0, nrow(exp_c))))
    return(NULL)
  }
  cat(sprintf("  Clumped: %d -> %d\n", n0, nrow(exp_c)))

  # Mean F-statistic for reporting
  mean_F <- mean((exp_c$beta.exposure / exp_c$se.exposure)^2)
  cat(sprintf("  Mean F-statistic: %.1f\n", mean_F))

  out <- format_outcome(als, exp_c$SNP)
  if (is.null(out)) { cat("  No outcome.\n"); return(NULL) }

  dat <- harmonise_data(exp_c, out)
  dat <- dat[dat$mr_keep, ]
  cat(sprintf("  Harmonized: %d\n", nrow(dat)))
  if (nrow(dat) < 3) { cat("  <3 harmonized.\n"); return(NULL) }

  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                  "mr_weighted_median", "mr_weighted_mode"))
  pleio <- mr_pleiotropy_test(dat)
  hetero <- mr_heterogeneity(dat)
  single <- mr_singlesnp(dat)
  loo <- mr_leaveoneout(dat)

  cat(sprintf("\n  --- %s ---\n", label))
  print(res[, c("method", "nsnp", "b", "se", "pval")])
  if (nrow(pleio) > 0) cat(sprintf("  Egger intercept p = %.4f\n", pleio$pval[1]))
  ivw_h <- hetero[hetero$method == "Inverse variance weighted", ]
  if (nrow(ivw_h) > 0) cat(sprintf("  Q_pval = %.4f\n", ivw_h$Q_pval[1]))

  tryCatch({
    ggsave(file.path(rdir, "scatter.png"), mr_scatter_plot(res, dat)[[1]], width=8, height=6)
    ggsave(file.path(rdir, "forest.png"), mr_forest_plot(single)[[1]], width=8, height=8)
    ggsave(file.path(rdir, "funnel.png"), mr_funnel_plot(single)[[1]], width=8, height=6)
    ggsave(file.path(rdir, "loo.png"), mr_leaveoneout_plot(loo)[[1]], width=8, height=8)
  }, error = function(e) cat(sprintf("  Plot err: %s\n", conditionMessage(e))))

  obj <- list(mr_results = res, dat = dat, pleiotropy = pleio,
              heterogeneity = hetero, single = single, loo = loo,
              n = nrow(dat), label = label, mean_F = mean_F)
  saveRDS(obj, file.path(rdir, "result.rds"))
  write.csv(res, file.path(rdir, "mr_results.csv"), row.names = FALSE)
  return(obj)
}

# ============================================================
# Power analysis
# ============================================================
mr_power <- function(n_iv, var_explained, n_outcome, alpha = 0.05) {
  # Approximate power for IVW MR
  # Based on Brion et al. 2013 / Burgess 2014
  # n_outcome = effective sample size of outcome GWAS
  # var_explained = variance in exposure explained by IVs (R2)
  # Returns: minimum detectable OR at 80% power
  ncp_per_unit <- n_outcome * var_explained
  se_approx <- 1 / sqrt(ncp_per_unit)

  # OR detectable at 80% power
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta  <- qnorm(0.80)
  min_beta <- (z_alpha + z_beta) * se_approx
  min_OR <- exp(min_beta)

  list(
    n_iv = n_iv,
    R2 = var_explained,
    n_outcome = n_outcome,
    se_approx = se_approx,
    min_beta = min_beta,
    min_OR = min_OR,
    min_OR_protective = exp(-min_beta)
  )
}

# ============================================================
# Main
# ============================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("MR with Relaxed p-threshold (p < 5e-6) + Power Analysis\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

snp_pos <- load_snp_pos()
als <- fread(ALS_GWAS_FILE)
cat(sprintf("ALS GWAS: %d SNPs, median N_eff = %d\n",
            nrow(als), median(als$N_effective, na.rm = TRUE)))

stable_high <- fread(STABLE_HIGH_FILE)
darkgrey <- fread(DARKGREY_FILE)
dg_col <- intersect(names(darkgrey), c("symbol", "gene", "gene_symbol"))
dg_symbols <- if (length(dg_col) > 0) darkgrey[[dg_col[1]]] else darkgrey[[1]]

rtk_genes <- c("NTRK2", "IGF1R", "INSR", "FYN", "DOCK3", "ROR1", "LRP4",
                "ADAM10", "ADAM12", "MEGF9", "NEO1", "EPS8", "EPS15",
                "RASA1", "RASGRF1", "RASGRF2", "RASGRP3", "MAP2K5",
                "MAPK10", "PRKCA")

R <- list()

# --- p < 5e-6 threshold ---
PTHRESH <- 5e-6

# Track A: Oligo all
oligo <- load_celltype_eqtl("Oligodendrocytes", snp_pos, pval_thresh = PTHRESH)
R[["A_Oligo"]] <- run_mr(to_exposure(oligo, "Oligo_eQTL"), "A_Oligo_5e6", als)

# Track B: Stable-High
oligo_sh <- load_celltype_eqtl("Oligodendrocytes", snp_pos,
                                gene_symbols = stable_high$symbol,
                                pval_thresh = PTHRESH)
if (nrow(oligo_sh) >= 3)
  R[["B_StableHigh"]] <- run_mr(to_exposure(oligo_sh, "StableHigh_eQTL"),
                                 "B_StableHigh_5e6", als)

# Track D: Darkgrey
oligo_dg <- load_celltype_eqtl("Oligodendrocytes", snp_pos,
                                gene_symbols = dg_symbols,
                                pval_thresh = PTHRESH)
if (nrow(oligo_dg) >= 3)
  R[["D_Darkgrey"]] <- run_mr(to_exposure(oligo_dg, "Darkgrey_eQTL"),
                                "D_Darkgrey_5e6", als)

# Track E: RTK
oligo_rtk <- load_celltype_eqtl("Oligodendrocytes", snp_pos,
                                  gene_symbols = rtk_genes,
                                  pval_thresh = PTHRESH)
if (nrow(oligo_rtk) >= 3)
  R[["E_RTK"]] <- run_mr(to_exposure(oligo_rtk, "RTK_eQTL"),
                           "E_RTK_5e6", als)

# ============================================================
# Power Analysis
# ============================================================
cat("\n\n", paste(rep("=", 70), collapse = ""), "\n")
cat("POWER ANALYSIS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

n_outcome_median <- median(als$N_effective, na.rm = TRUE)
cat(sprintf("Outcome N_effective (median): %d\n\n", n_outcome_median))

power_results <- list()
for (nm in names(R)) {
  r <- R[[nm]]
  if (is.null(r)) next

  # Estimate R2 from F-statistics
  n_iv <- r$n
  mean_F <- r$mean_F
  # R2 approx = n_iv * F / (n_iv * F + n_exposure)
  # For eQTL: n_exposure ~ 196 (Bryois, ~196 individuals)
  n_exposure <- 196
  R2_approx <- (n_iv * mean_F) / (n_iv * mean_F + n_exposure)

  pw <- mr_power(n_iv, R2_approx, n_outcome_median)
  pw$track <- nm
  pw$mean_F <- mean_F
  power_results[[nm]] <- pw

  cat(sprintf("[%s] n_IV=%d, mean_F=%.1f, R2=%.4f, min_detectable_OR=%.3f (protective: %.3f)\n",
              nm, n_iv, mean_F, R2_approx, pw$min_OR, pw$min_OR_protective))
}

# ============================================================
# Summary
# ============================================================
cat("\n\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SUMMARY (p < 5e-6, LD CLUMPED)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

sdf <- do.call(rbind, lapply(names(R), function(nm) {
  r <- R[[nm]]
  if (is.null(r)) return(NULL)
  ivw <- r$mr_results[r$mr_results$method == "Inverse variance weighted", ]
  if (nrow(ivw) == 0) return(NULL)
  pw <- power_results[[nm]]
  data.frame(
    track = nm, n_iv = ivw$nsnp, mean_F = round(r$mean_F, 1),
    beta = round(ivw$b, 4), se = round(ivw$se, 4),
    pval = signif(ivw$pval, 3),
    OR = round(exp(ivw$b), 3),
    CI_lo = round(exp(ivw$b - 1.96 * ivw$se), 3),
    CI_hi = round(exp(ivw$b + 1.96 * ivw$se), 3),
    egger_p = signif(r$pleiotropy$pval[1], 3),
    min_OR_80pwr = round(pw$min_OR, 3),
    stringsAsFactors = FALSE
  )
}))

if (!is.null(sdf) && nrow(sdf) > 0) {
  print(sdf, row.names = FALSE)
  write.csv(sdf, file.path(OUTPUT_DIR, "summary_relaxed.csv"), row.names = FALSE)

  p <- ggplot(sdf, aes(x = reorder(track, OR), y = OR)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = "MR: Oligo eQTL -> ALS (p<5e-6, LD clumped)",
         subtitle = sprintf("Labels show min detectable OR at 80%% power"),
         x = "", y = "Odds Ratio (95% CI)") +
    geom_text(aes(label = sprintf("min OR=%.2f", min_OR_80pwr)),
              hjust = -0.1, size = 3) +
    theme_minimal(base_size = 14)
  ggsave(file.path(OUTPUT_DIR, "forest_relaxed.png"), p, width = 11, height = 5)
}

# Also run at p < 1e-5 for comparison
cat("\n\n--- Additional: p < 1e-5 for Track A ---\n")
oligo_1e5 <- load_celltype_eqtl("Oligodendrocytes", snp_pos, pval_thresh = 1e-5)
if (nrow(oligo_1e5) >= 3) {
  exp_1e5 <- to_exposure(oligo_1e5, "Oligo_eQTL_1e5")
  r_1e5 <- run_mr(exp_1e5, "A_Oligo_1e5", als)
  if (!is.null(r_1e5)) {
    ivw_1e5 <- r_1e5$mr_results[r_1e5$mr_results$method == "Inverse variance weighted", ]
    cat(sprintf("\n  Track A (p<1e-5): n=%d, beta=%.4f, p=%.3e, OR=%.3f\n",
                ivw_1e5$nsnp, ivw_1e5$b, ivw_1e5$pval, exp(ivw_1e5$b)))
  }
}

cat("\nDone.\n")
