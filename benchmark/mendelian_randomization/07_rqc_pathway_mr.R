#!/usr/bin/env Rscript
# ============================================================
# 07_rqc_pathway_mr.R
# Track F: RQC (Ribosome-associated Quality Control) pathway
# MR + Coloc + Power Analysis + Per-gene breakdown
# ============================================================

user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))
readRenviron(file.path(Sys.getenv("USERPROFILE"), ".Renviron"))

library(TwoSampleMR)
library(ieugwasr)
library(coloc)
library(data.table)
library(dplyr)
library(ggplot2)

DATA_DIR      <- "D:/Projects/MR/data/bryois_eqtl"
ALS_GWAS_FILE <- "D:/Projects/MR/data/als_gwas/GCST90027164_buildGRCh37.tsv.gz"
OUTPUT_DIR    <- "D:/Projects/MR/results_rqc"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PTHRESH <- 5e-6

# ============================================================
# RQC gene sets
# ============================================================
rqc_core <- c("LTN1", "NEMF", "TCF25", "ZNF598")
rqc_related <- c("GCN1", "EIF2AK4", "VCP", "UFD1", "NPL4",
                  "ABCE1", "RNF25", "CNOT4", "PELO", "HBS1L")
rqc_all <- c(rqc_core, rqc_related)

cat("RQC gene set:\n")
cat(sprintf("  Core (%d): %s\n", length(rqc_core), paste(rqc_core, collapse = ", ")))
cat(sprintf("  Related (%d): %s\n", length(rqc_related), paste(rqc_related, collapse = ", ")))
cat(sprintf("  Total: %d\n\n", length(rqc_all)))

# ============================================================
# Load data
# ============================================================
cat("Loading SNP positions...\n")
snp_pos <- fread(file.path(DATA_DIR, "snp_pos.txt.gz"))
snp_pos[, chr := as.integer(sub("chr", "", sub(":.*", "", SNP_id_hg19)))]
snp_pos[, pos := as.integer(sub(".*:", "", SNP_id_hg19))]

cat("Loading ALS GWAS...\n")
als <- fread(ALS_GWAS_FILE)

cat("Loading Oligo eQTL (all chr)...\n")
all_chr <- lapply(1:22, function(chr) {
  f <- file.path(DATA_DIR, "Oligodendrocytes",
                 sprintf("Oligodendrocytes.%d.gz", chr))
  if (!file.exists(f)) return(NULL)
  fread(f, header = FALSE)
})
eqtl_raw <- rbindlist(all_chr[!sapply(all_chr, is.null)])
setnames(eqtl_raw, c("gene_id", "snp_id", "dist_tss", "pvalue", "beta"))
eqtl_raw[, gene_symbol := sub("_ENSG.*", "", gene_id)]

# Merge positions
eqtl_full <- merge(eqtl_raw, snp_pos[, .(SNP, chr, pos, MAF, effect_allele, other_allele)],
                   by.x = "snp_id", by.y = "SNP")
eqtl_full[, se := abs(beta) / qnorm(pvalue / 2, lower.tail = FALSE)]
eqtl_full[!is.finite(se) | se <= 0, se := NA]

cat(sprintf("Total eQTL: %d pairs\n\n", nrow(eqtl_full)))

# ============================================================
# Check which RQC genes have eQTL data
# ============================================================
cat("=== RQC Gene Coverage in Bryois Oligo eQTL ===\n")
rqc_coverage <- data.frame(
  gene = rqc_all,
  category = ifelse(rqc_all %in% rqc_core, "Core RQC", "Related"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(rqc_coverage))) {
  g <- rqc_coverage$gene[i]
  g_data <- eqtl_full[gene_symbol == g]
  rqc_coverage$n_snps[i] <- nrow(g_data)
  rqc_coverage$n_sig_5e6[i] <- sum(g_data$pvalue < 5e-6, na.rm = TRUE)
  rqc_coverage$n_sig_5e8[i] <- sum(g_data$pvalue < 5e-8, na.rm = TRUE)
  rqc_coverage$min_p[i] <- ifelse(nrow(g_data) > 0, min(g_data$pvalue, na.rm = TRUE), NA)
}

cat("\n")
print(rqc_coverage, row.names = FALSE)
write.csv(rqc_coverage, file.path(OUTPUT_DIR, "rqc_gene_coverage.csv"), row.names = FALSE)

# ============================================================
# Prepare exposure (RQC pathway combined)
# ============================================================
cat("\n=== Preparing RQC exposure (p < 5e-6) ===\n")

rqc_eqtl <- eqtl_full[gene_symbol %in% rqc_all & pvalue < PTHRESH & !is.na(se)]
rqc_eqtl[, F_stat := (beta / se)^2]
rqc_eqtl <- rqc_eqtl[F_stat > 10]

cat(sprintf("RQC instruments (p<5e-6, F>10): %d from %d genes\n",
            nrow(rqc_eqtl), uniqueN(rqc_eqtl$gene_symbol)))
cat(sprintf("Genes with instruments: %s\n",
            paste(unique(rqc_eqtl$gene_symbol), collapse = ", ")))

if (nrow(rqc_eqtl) == 0) {
  cat("\nWARNING: No RQC instruments at p<5e-6. Trying p<1e-4...\n")
  rqc_eqtl <- eqtl_full[gene_symbol %in% rqc_all & pvalue < 1e-4 & !is.na(se)]
  rqc_eqtl[, F_stat := (beta / se)^2]
  rqc_eqtl <- rqc_eqtl[F_stat > 10]
  cat(sprintf("RQC instruments (p<1e-4, F>10): %d from %d genes\n",
              nrow(rqc_eqtl), uniqueN(rqc_eqtl$gene_symbol)))
  PTHRESH_USED <- 1e-4
} else {
  PTHRESH_USED <- PTHRESH
}

# Format exposure
if (nrow(rqc_eqtl) > 0) {
  rqc_eqtl[, chr_num := chr]
  rqc_eqtl[, pos_hg19 := pos]

  rqc_exp <- data.frame(
    SNP = rqc_eqtl$snp_id, beta.exposure = rqc_eqtl$beta,
    se.exposure = rqc_eqtl$se, pval.exposure = rqc_eqtl$pvalue,
    effect_allele.exposure = rqc_eqtl$effect_allele,
    other_allele.exposure = rqc_eqtl$other_allele,
    eaf.exposure = rqc_eqtl$MAF,
    exposure = "RQC_eQTL", id.exposure = "RQC_eQTL",
    gene.exposure = rqc_eqtl$gene_symbol,
    chr.exposure = rqc_eqtl$chr_num,
    pos.exposure = rqc_eqtl$pos_hg19,
    mr_keep.exposure = TRUE, stringsAsFactors = FALSE
  ) %>%
    group_by(SNP) %>%
    slice_min(pval.exposure, n = 1, with_ties = FALSE) %>%
    ungroup() %>% as.data.frame()

  cat(sprintf("Unique instruments: %d\n", nrow(rqc_exp)))
}

# ============================================================
# MR: RQC pathway combined (Track F)
# ============================================================
run_track_f <- function(exp_dat, label, als_data) {
  cat(sprintf("\n========== %s ==========\n", label))

  rdir <- file.path(OUTPUT_DIR, label)
  dir.create(rdir, showWarnings = FALSE, recursive = TRUE)

  n0 <- nrow(exp_dat)
  cat(sprintf("  Pre-clump: %d\n", n0))

  # Clump
  exp_c <- tryCatch(
    clump_data(exp_dat, clump_r2 = 0.001, clump_kb = 10000),
    error = function(e) {
      cat(sprintf("  Clump err: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(exp_c) || nrow(exp_c) < 3) {
    n_c <- ifelse(is.null(exp_c), 0, nrow(exp_c))
    cat(sprintf("  %d -> %d IVs. ", n0, n_c))
    if (n_c > 0) {
      cat("Running with available IVs (Wald ratio if n=1, or limited MR).\n")
      exp_c_final <- exp_c
    } else {
      cat("SKIP.\n")
      return(NULL)
    }
  } else {
    exp_c_final <- exp_c
  }

  cat(sprintf("  Clumped: %d -> %d\n", n0, nrow(exp_c_final)))
  mean_F <- mean((exp_c_final$beta.exposure / exp_c_final$se.exposure)^2)
  cat(sprintf("  Mean F: %.1f\n", mean_F))

  # Outcome
  als_filt <- als_data[rsid %in% exp_c_final$SNP]
  if (nrow(als_filt) == 0) { cat("  No outcome.\n"); return(NULL) }

  out <- data.frame(
    SNP = als_filt$rsid, beta.outcome = als_filt$beta,
    se.outcome = als_filt$standard_error, pval.outcome = als_filt$p_value,
    effect_allele.outcome = toupper(als_filt$effect_allele),
    other_allele.outcome = toupper(als_filt$other_allele),
    eaf.outcome = als_filt$effect_allele_frequency,
    outcome = "ALS", id.outcome = "ALS", mr_keep.outcome = TRUE,
    samplesize.outcome = als_filt$N_effective, stringsAsFactors = FALSE
  )

  dat <- harmonise_data(exp_c_final, out)
  dat <- dat[dat$mr_keep, ]
  cat(sprintf("  Harmonized: %d\n", nrow(dat)))

  if (nrow(dat) == 0) return(NULL)

  # MR (use Wald ratio if only 1 SNP)
  if (nrow(dat) == 1) {
    res <- mr(dat, method_list = "mr_wald_ratio")
  } else if (nrow(dat) == 2) {
    res <- mr(dat, method_list = c("mr_ivw"))
  } else {
    res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                    "mr_weighted_median", "mr_weighted_mode"))
  }

  # Sensitivity
  pleio <- if (nrow(dat) >= 3) mr_pleiotropy_test(dat) else data.frame(pval = NA)
  hetero <- if (nrow(dat) >= 3) mr_heterogeneity(dat) else data.frame(Q_pval = NA)
  single <- if (nrow(dat) >= 2) mr_singlesnp(dat) else NULL
  loo <- if (nrow(dat) >= 3) mr_leaveoneout(dat) else NULL

  cat(sprintf("\n  --- %s ---\n", label))
  print(res[, c("method", "nsnp", "b", "se", "pval")])

  # Plots
  if (!is.null(single) && nrow(dat) >= 3) {
    tryCatch({
      ggsave(file.path(rdir, "scatter.png"), mr_scatter_plot(res, dat)[[1]], width=8, height=6)
      ggsave(file.path(rdir, "forest.png"), mr_forest_plot(single)[[1]], width=8, height=8)
      ggsave(file.path(rdir, "funnel.png"), mr_funnel_plot(single)[[1]], width=8, height=6)
      if (!is.null(loo))
        ggsave(file.path(rdir, "loo.png"), mr_leaveoneout_plot(loo)[[1]], width=8, height=8)
    }, error = function(e) cat(sprintf("  Plot err: %s\n", conditionMessage(e))))
  }

  obj <- list(mr_results = res, dat = dat, pleiotropy = pleio,
              heterogeneity = hetero, single = single, loo = loo,
              n = nrow(dat), label = label, mean_F = mean_F)
  saveRDS(obj, file.path(rdir, "result.rds"))
  write.csv(res, file.path(rdir, "mr_results.csv"), row.names = FALSE)
  return(obj)
}

# ============================================================
# Per-gene MR (individual gene analysis)
# ============================================================
run_per_gene_mr <- function(gene_sym, eqtl_data, als_data, pthresh) {
  g_eqtl <- eqtl_data[gene_symbol == gene_sym & pvalue < pthresh & !is.na(se)]
  g_eqtl[, F_stat := (beta / se)^2]
  g_eqtl <- g_eqtl[F_stat > 10]

  if (nrow(g_eqtl) == 0)
    return(data.frame(gene = gene_sym, n_iv_pre = 0, n_iv_post = 0,
                      beta = NA, se = NA, pval = NA, method = NA,
                      stringsAsFactors = FALSE))

  exp <- data.frame(
    SNP = g_eqtl$snp_id, beta.exposure = g_eqtl$beta,
    se.exposure = g_eqtl$se, pval.exposure = g_eqtl$pvalue,
    effect_allele.exposure = g_eqtl$effect_allele,
    other_allele.exposure = g_eqtl$other_allele,
    eaf.exposure = g_eqtl$MAF,
    exposure = gene_sym, id.exposure = gene_sym,
    gene.exposure = gene_sym,
    mr_keep.exposure = TRUE, stringsAsFactors = FALSE
  ) %>% group_by(SNP) %>%
    slice_min(pval.exposure, n = 1, with_ties = FALSE) %>%
    ungroup() %>% as.data.frame()

  n_pre <- nrow(exp)

  # Clump
  exp_c <- tryCatch(
    clump_data(exp, clump_r2 = 0.001, clump_kb = 10000),
    error = function(e) NULL
  )
  if (is.null(exp_c) || nrow(exp_c) == 0)
    return(data.frame(gene = gene_sym, n_iv_pre = n_pre, n_iv_post = 0,
                      beta = NA, se = NA, pval = NA, method = NA,
                      stringsAsFactors = FALSE))

  # Outcome
  als_f <- als_data[rsid %in% exp_c$SNP]
  if (nrow(als_f) == 0)
    return(data.frame(gene = gene_sym, n_iv_pre = n_pre, n_iv_post = nrow(exp_c),
                      beta = NA, se = NA, pval = NA, method = "no_outcome",
                      stringsAsFactors = FALSE))

  out <- data.frame(
    SNP = als_f$rsid, beta.outcome = als_f$beta,
    se.outcome = als_f$standard_error, pval.outcome = als_f$p_value,
    effect_allele.outcome = toupper(als_f$effect_allele),
    other_allele.outcome = toupper(als_f$other_allele),
    eaf.outcome = als_f$effect_allele_frequency,
    outcome = "ALS", id.outcome = "ALS", mr_keep.outcome = TRUE,
    stringsAsFactors = FALSE
  )

  dat <- tryCatch(harmonise_data(exp_c, out), error = function(e) NULL)
  if (is.null(dat) || nrow(dat[dat$mr_keep, ]) == 0)
    return(data.frame(gene = gene_sym, n_iv_pre = n_pre, n_iv_post = nrow(exp_c),
                      beta = NA, se = NA, pval = NA, method = "harmonize_fail",
                      stringsAsFactors = FALSE))

  dat <- dat[dat$mr_keep, ]

  if (nrow(dat) == 1) {
    res <- mr(dat, method_list = "mr_wald_ratio")
  } else if (nrow(dat) >= 3) {
    res <- mr(dat, method_list = "mr_ivw")
  } else {
    res <- mr(dat, method_list = "mr_ivw")
  }

  ivw <- res[1, ]
  data.frame(gene = gene_sym, n_iv_pre = n_pre, n_iv_post = nrow(dat),
             beta = ivw$b, se = ivw$se, pval = ivw$pval,
             method = ivw$method, stringsAsFactors = FALSE)
}

# ============================================================
# Colocalization for RQC genes
# ============================================================
run_coloc_gene <- function(gene_sym, eqtl_data, als_data,
                            window = 500000, n_eqtl = 196) {
  g_eqtl <- eqtl_data[gene_symbol == gene_sym & !is.na(se)]
  if (nrow(g_eqtl) < 10) return(NULL)

  lead <- g_eqtl[which.min(pvalue)]
  g_win <- g_eqtl[chr == lead$chr & pos >= (lead$pos - window) & pos <= (lead$pos + window)]
  als_win <- als_data[chromosome == lead$chr &
                        base_pair_location >= (lead$pos - window) &
                        base_pair_location <= (lead$pos + window)]

  common <- intersect(g_win$snp_id, als_win$rsid)
  if (length(common) < 10) return(NULL)

  e_sub <- g_win[snp_id %in% common][!duplicated(snp_id)][order(snp_id)]
  a_sub <- als_win[rsid %in% common][!duplicated(rsid)]
  a_sub <- a_sub[match(e_sub$snp_id, rsid)]

  keep <- !is.na(e_sub$se) & !is.na(a_sub$standard_error) &
          e_sub$se > 0 & a_sub$standard_error > 0
  e_sub <- e_sub[keep]; a_sub <- a_sub[keep]
  if (nrow(e_sub) < 10) return(NULL)

  d1 <- list(beta = e_sub$beta, varbeta = e_sub$se^2, snp = e_sub$snp_id,
             position = e_sub$pos, type = "quant", N = n_eqtl,
             MAF = e_sub$MAF, sdY = 1)
  d2 <- list(beta = a_sub$beta, varbeta = a_sub$standard_error^2,
             snp = a_sub$rsid, position = a_sub$base_pair_location,
             type = "cc", N = median(a_sub$N_effective, na.rm = TRUE), s = 0.5)

  res <- tryCatch(coloc.abf(d1, d2), error = function(e) NULL)
  if (is.null(res)) return(NULL)

  pp <- res$summary
  list(gene = gene_sym, n_snps = nrow(e_sub), lead_snp = lead$snp_id,
       lead_pval = lead$pvalue, chr = lead$chr, pos = lead$pos,
       PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"],
       PP.H2 = pp["PP.H2.abf"], PP.H3 = pp["PP.H3.abf"],
       PP.H4 = pp["PP.H4.abf"])
}

# ============================================================
# MAIN EXECUTION
# ============================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Track F: RQC Pathway -> ALS Risk\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# --- Combined pathway MR ---
if (nrow(rqc_eqtl) > 0) {
  res_combined <- run_track_f(rqc_exp, "F_RQC_combined", als)
} else {
  cat("No instruments for combined analysis.\n")
  res_combined <- NULL
}

# --- Per-gene MR ---
cat("\n=== Per-Gene MR ===\n")
per_gene_results <- lapply(rqc_all, function(g) {
  cat(sprintf("  %s...\n", g))
  run_per_gene_mr(g, eqtl_full, als, PTHRESH_USED)
})
per_gene_df <- do.call(rbind, per_gene_results)
per_gene_df$category <- ifelse(per_gene_df$gene %in% rqc_core, "Core RQC", "Related")
per_gene_df$OR <- exp(per_gene_df$beta)

cat("\n=== Per-Gene Results ===\n")
print(per_gene_df[, c("gene", "category", "n_iv_pre", "n_iv_post",
                       "beta", "pval", "OR", "method")], row.names = FALSE)
write.csv(per_gene_df, file.path(OUTPUT_DIR, "per_gene_mr.csv"), row.names = FALSE)

# --- Colocalization ---
cat("\n=== Colocalization (RQC genes) ===\n")
coloc_list <- list()
for (g in rqc_all) {
  cat(sprintf("  coloc %s...\n", g))
  r <- run_coloc_gene(g, eqtl_full, als)
  if (!is.null(r)) coloc_list[[g]] <- r
}

if (length(coloc_list) > 0) {
  coloc_df <- do.call(rbind, lapply(coloc_list, function(x) {
    data.frame(gene = x$gene, n_snps = x$n_snps, lead_snp = x$lead_snp,
               lead_pval = x$lead_pval, chr = x$chr,
               PP.H3 = round(x$PP.H3, 4), PP.H4 = round(x$PP.H4, 4),
               stringsAsFactors = FALSE)
  }))
  coloc_df <- coloc_df[order(-coloc_df$PP.H4), ]
  cat("\n")
  print(coloc_df, row.names = FALSE)
  write.csv(coloc_df, file.path(OUTPUT_DIR, "coloc_rqc.csv"), row.names = FALSE)
} else {
  cat("  No genes with sufficient data for coloc.\n")
  coloc_df <- NULL
}

# --- Power analysis ---
cat("\n=== Power Analysis ===\n")
if (!is.null(res_combined) && res_combined$n >= 1) {
  n_out <- median(als$N_effective, na.rm = TRUE)
  n_iv <- res_combined$n
  mean_F <- res_combined$mean_F
  n_exp <- 196
  R2 <- (n_iv * mean_F) / (n_iv * mean_F + n_exp)
  se_approx <- 1 / sqrt(n_out * R2)
  min_beta <- (qnorm(0.975) + qnorm(0.80)) * se_approx
  cat(sprintf("  N_IV: %d, R2: %.4f, N_outcome: %d\n", n_iv, R2, n_out))
  cat(sprintf("  Min detectable beta (80%% power): %.4f\n", min_beta))
  cat(sprintf("  Min detectable OR: %.3f (risk) / %.3f (protective)\n",
              exp(min_beta), exp(-min_beta)))
}

# --- Comparison with Track B ---
cat("\n\n", paste(rep("=", 70), collapse = ""), "\n")
cat("COMPARISON: Track B (Stable-High) vs Track F (RQC)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Track B (Stable-High 135 ribosomal genes):\n")
cat("  N IVs: 15, IVW OR: 0.988 (0.972-1.005), p=0.159\n")
cat("  Egger intercept p=0.791, Q p=0.493\n")
cat("  Direction: All 4 methods protective\n\n")

if (!is.null(res_combined)) {
  ivw_f <- res_combined$mr_results[res_combined$mr_results$method == "Inverse variance weighted" |
                                    res_combined$mr_results$method == "Wald ratio", ]
  if (nrow(ivw_f) > 0) {
    cat(sprintf("Track F (RQC %d genes):\n", length(rqc_all)))
    cat(sprintf("  N IVs: %d, IVW OR: %.3f (%.3f-%.3f), p=%.3f\n",
                ivw_f$nsnp, exp(ivw_f$b),
                exp(ivw_f$b - 1.96 * ivw_f$se),
                exp(ivw_f$b + 1.96 * ivw_f$se),
                ivw_f$pval))
  }
}

# Gene overlap
cat("\n=== Gene Overlap ===\n")
stable_high <- fread("D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/sals_upstream_gene_list/stable_high_genes.csv")
overlap <- intersect(rqc_all, stable_high$symbol)
cat(sprintf("RQC genes in Stable-High: %d / %d\n", length(overlap), length(rqc_all)))
if (length(overlap) > 0) cat(sprintf("  Overlap: %s\n", paste(overlap, collapse = ", ")))

cat("\nTrack F analysis complete.\n")
