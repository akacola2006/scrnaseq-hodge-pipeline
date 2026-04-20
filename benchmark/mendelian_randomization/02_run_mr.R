#!/usr/bin/env Rscript
# ============================================================
# 02_run_mr.R
# Mendelian Randomization: Brain cell-type eQTL -> ALS risk
#
# Data sources:
#   Exposure: Bryois et al. 2022 Nat Neurosci (Zenodo 7276971)
#   Outcome:  van Rheenen et al. 2021 (GCST90027164)
#
# Tracks:
#   A: Oligo all eQTL -> ALS
#   B: Stable-High 135 genes only
#   C: All 8 cell types comparison (when data available)
#   D: Darkgrey 327 genes (structural source)
#   E: RTK/Trophic pathway genes (input pathway)
# ============================================================

user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)

# --- Configuration ---
DATA_DIR      <- "D:/Projects/MR/data/bryois_eqtl"
ALS_GWAS_FILE <- "D:/Projects/MR/data/als_gwas/GCST90027164_buildGRCh37.tsv.gz"
OUTPUT_DIR    <- "D:/Projects/MR/results"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

STABLE_HIGH_FILE <- "D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/sals_upstream_gene_list/stable_high_genes.csv"
DARKGREY_FILE    <- "D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/laneB/darkgrey_identity/darkgrey_genes.tsv"

# ============================================================
# Load SNP position file (with alleles and MAF)
# Columns: SNP, SNP_id_hg38, SNP_id_hg19, effect_allele, other_allele, MAF
# ============================================================
load_snp_pos <- function() {
  cat("Loading SNP positions...\n")
  snp_pos <- fread(file.path(DATA_DIR, "snp_pos.txt.gz"))
  cat(sprintf("  %d SNPs loaded\n", nrow(snp_pos)))
  return(snp_pos)
}

# ============================================================
# Load Bryois eQTL for one cell type
# File format: gene_symbol_ENSGID, rsID, dist_tss, pvalue, beta
# ============================================================
load_celltype_eqtl <- function(cell_type, snp_pos, gene_symbols = NULL) {
  cat(sprintf("\n--- Loading %s eQTL ---\n", cell_type))

  all_chr <- lapply(1:22, function(chr) {
    f <- file.path(DATA_DIR, cell_type, sprintf("%s.%d.gz", cell_type, chr))
    if (!file.exists(f)) return(NULL)
    fread(f, header = FALSE)
  })
  eqtl <- rbindlist(all_chr[!sapply(all_chr, is.null)])
  setnames(eqtl, c("gene_id", "snp_id", "dist_tss", "pvalue", "beta"))

  # Parse gene_id: "SYMBOL_ENSGID" format
  eqtl[, gene_symbol := sub("_ENSG.*", "", gene_id)]
  eqtl[, ensg_id := sub(".*_(ENSG[0-9]+)$", "\\1", gene_id)]

  cat(sprintf("  Raw: %d pairs, %d genes, %d SNPs\n",
              nrow(eqtl), uniqueN(eqtl$gene_id), uniqueN(eqtl$snp_id)))

  # Filter by gene symbols if provided
  if (!is.null(gene_symbols)) {
    eqtl <- eqtl[gene_symbol %in% gene_symbols]
    cat(sprintf("  After gene filter (%d targets): %d pairs, %d genes matched\n",
                length(gene_symbols), nrow(eqtl), uniqueN(eqtl$gene_symbol)))
  }

  # Merge with SNP positions (alleles + MAF)
  eqtl_m <- merge(eqtl, snp_pos, by.x = "snp_id", by.y = "SNP")
  cat(sprintf("  After SNP merge: %d pairs\n", nrow(eqtl_m)))

  # Compute SE from beta and p-value
  eqtl_m[, se := abs(beta) / qnorm(pvalue / 2, lower.tail = FALSE)]
  eqtl_m[!is.finite(se) | se <= 0, se := NA]

  # Filter genome-wide significant
  eqtl_sig <- eqtl_m[pvalue < 5e-8 & !is.na(se)]
  cat(sprintf("  Genome-wide significant (p<5e-8): %d\n", nrow(eqtl_sig)))

  if (nrow(eqtl_sig) < 3) {
    cat("  Trying suggestive threshold p<5e-6...\n")
    eqtl_sig <- eqtl_m[pvalue < 5e-6 & !is.na(se)]
    cat(sprintf("  Suggestive (p<5e-6): %d\n", nrow(eqtl_sig)))
  }

  # F-statistic filter
  eqtl_sig[, F_stat := (beta / se)^2]
  n_weak <- sum(eqtl_sig$F_stat <= 10, na.rm = TRUE)
  eqtl_sig <- eqtl_sig[F_stat > 10]
  cat(sprintf("  After F>10: %d instruments (removed %d weak)\n",
              nrow(eqtl_sig), n_weak))

  if (nrow(eqtl_sig) > 0) {
    cat(sprintf("  Median F: %.1f | Unique genes: %d\n",
                median(eqtl_sig$F_stat), uniqueN(eqtl_sig$gene_symbol)))
  }

  return(eqtl_sig)
}

# ============================================================
# Convert to TwoSampleMR exposure format
# ============================================================
to_exposure_format <- function(eqtl_sig, label) {
  # Extract chr and pos from SNP_id_hg19 (format: "chr1:234313")
  eqtl_sig[, chr_num := as.integer(sub("chr", "", sub(":.*", "", SNP_id_hg19)))]
  eqtl_sig[, pos_hg19 := as.integer(sub(".*:", "", SNP_id_hg19))]

  exposure <- data.frame(
    SNP = eqtl_sig$snp_id,
    beta.exposure = eqtl_sig$beta,
    se.exposure = eqtl_sig$se,
    pval.exposure = eqtl_sig$pvalue,
    effect_allele.exposure = eqtl_sig$effect_allele,
    other_allele.exposure = eqtl_sig$other_allele,
    eaf.exposure = eqtl_sig$MAF,
    exposure = label,
    id.exposure = label,
    gene.exposure = eqtl_sig$gene_symbol,
    chr.exposure = eqtl_sig$chr_num,
    pos.exposure = eqtl_sig$pos_hg19,
    mr_keep.exposure = TRUE,
    stringsAsFactors = FALSE
  )

  # For multiple genes per SNP, keep the strongest association
  exposure <- exposure %>%
    group_by(SNP) %>%
    slice_min(pval.exposure, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    as.data.frame()

  cat(sprintf("  Exposure formatted: %d unique instruments\n", nrow(exposure)))
  return(exposure)
}

# ============================================================
# Load ALS GWAS as local outcome
# ============================================================
load_als_gwas_local <- function() {
  cat("\nLoading local ALS GWAS...\n")
  if (!file.exists(ALS_GWAS_FILE) || file.size(ALS_GWAS_FILE) < 1000) {
    cat("  File not available. Will use IEU OpenGWAS API.\n")
    return(NULL)
  }
  als <- fread(ALS_GWAS_FILE, nrows = 5)
  cat(sprintf("  Columns: %s\n", paste(names(als), collapse = ", ")))
  print(head(als, 2))
  als_full <- fread(ALS_GWAS_FILE)
  cat(sprintf("  %d total SNPs\n", nrow(als_full)))
  return(als_full)
}

# ============================================================
# Format local ALS GWAS as outcome
# Known columns: rsid, chromosome, base_pair_location,
#   effect_allele, other_allele, effect_allele_frequency,
#   beta, standard_error, p_value, N_effective
# ============================================================
format_als_outcome <- function(als_gwas, snps) {
  als_filt <- als_gwas[rsid %in% snps]
  cat(sprintf("  ALS outcome: matched %d of %d SNPs\n", nrow(als_filt), length(snps)))

  if (nrow(als_filt) == 0) return(NULL)

  # Alleles need to be uppercase to match Bryois format
  outcome <- data.frame(
    SNP = als_filt$rsid,
    beta.outcome = als_filt$beta,
    se.outcome = als_filt$standard_error,
    pval.outcome = als_filt$p_value,
    effect_allele.outcome = toupper(als_filt$effect_allele),
    other_allele.outcome = toupper(als_filt$other_allele),
    eaf.outcome = als_filt$effect_allele_frequency,
    outcome = "ALS_vanRheenen2021",
    id.outcome = "ALS",
    mr_keep.outcome = TRUE,
    samplesize.outcome = als_filt$N_effective,
    stringsAsFactors = FALSE
  )

  return(outcome)
}

# ============================================================
# Core MR pipeline
# ============================================================
run_mr_analysis <- function(exposure_dat, label, als_gwas = NULL,
                            output_dir = OUTPUT_DIR) {
  cat(sprintf("\n========== MR: %s ==========\n", label))

  result_dir <- file.path(output_dir, label)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  # LD clumping via IEU API
  cat("  LD clumping...\n")
  n_pre <- nrow(exposure_dat)
  exposure_clumped <- tryCatch({
    clump_data(exposure_dat, clump_r2 = 0.001, clump_kb = 10000)
  }, error = function(e) {
    cat(sprintf("  Clumping API error: %s\n", conditionMessage(e)))
    cat("  Falling back to p-value sorting (no LD pruning)...\n")
    exposure_dat[order(exposure_dat$pval.exposure), ]
  })
  cat(sprintf("  Clumped: %d -> %d\n", n_pre, nrow(exposure_clumped)))

  if (nrow(exposure_clumped) < 3) {
    cat(sprintf("  SKIP: <3 instruments for %s\n", label))
    saveRDS(list(status = "insufficient_instruments", n = nrow(exposure_clumped)),
            file.path(result_dir, "result.rds"))
    return(NULL)
  }

  # Outcome data
  cat("  Getting outcome data...\n")
  if (!is.null(als_gwas)) {
    outcome_dat <- format_als_outcome(als_gwas, exposure_clumped$SNP)
  } else {
    outcome_dat <- tryCatch({
      extract_outcome_data(snps = exposure_clumped$SNP,
                            outcomes = "ebi-a-GCST90027164")
    }, error = function(e) {
      cat(sprintf("  Outcome error: %s\n", conditionMessage(e)))
      return(NULL)
    })
  }

  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    cat("  SKIP: No outcome data.\n")
    return(NULL)
  }
  cat(sprintf("  Outcome SNPs: %d\n", nrow(outcome_dat)))

  # Harmonize
  dat <- harmonise_data(exposure_clumped, outcome_dat)
  dat <- dat[dat$mr_keep == TRUE, ]
  cat(sprintf("  Harmonized: %d SNPs\n", nrow(dat)))

  if (nrow(dat) < 3) {
    cat(sprintf("  SKIP: <3 harmonized for %s\n", label))
    return(NULL)
  }

  # MR (4 methods)
  cat("  Running MR...\n")
  mr_res <- mr(dat, method_list = c(
    "mr_ivw", "mr_egger_regression",
    "mr_weighted_median", "mr_weighted_mode"))

  # Sensitivity
  pleio <- mr_pleiotropy_test(dat)
  hetero <- mr_heterogeneity(dat)
  single <- mr_singlesnp(dat)
  loo <- mr_leaveoneout(dat)

  # Print results
  cat(sprintf("\n  --- %s Results ---\n", label))
  print(mr_res[, c("method", "nsnp", "b", "se", "pval")])
  if (nrow(pleio) > 0) cat(sprintf("  Egger intercept p = %.4f\n", pleio$pval[1]))
  ivw_h <- hetero[hetero$method == "Inverse variance weighted", ]
  if (nrow(ivw_h) > 0) cat(sprintf("  IVW Q_pval = %.4f\n", ivw_h$Q_pval[1]))

  # Plots
  tryCatch({
    p1 <- mr_scatter_plot(mr_res, dat)
    ggsave(file.path(result_dir, "scatter.png"), p1[[1]], width = 8, height = 6)
    p2 <- mr_forest_plot(single)
    ggsave(file.path(result_dir, "forest.png"), p2[[1]], width = 8, height = 8)
    p3 <- mr_funnel_plot(single)
    ggsave(file.path(result_dir, "funnel.png"), p3[[1]], width = 8, height = 6)
    p4 <- mr_leaveoneout_plot(loo)
    ggsave(file.path(result_dir, "leaveoneout.png"), p4[[1]], width = 8, height = 8)
    cat("  Plots saved.\n")
  }, error = function(e) cat(sprintf("  Plot error: %s\n", conditionMessage(e))))

  # Save
  result <- list(mr_results = mr_res, dat = dat, pleiotropy = pleio,
                 heterogeneity = hetero, single_snp = single, loo = loo,
                 n_instruments = nrow(dat), label = label)
  saveRDS(result, file.path(result_dir, "result.rds"))
  write.csv(mr_res, file.path(result_dir, "mr_results.csv"), row.names = FALSE)
  return(result)
}

# ============================================================
# Main
# ============================================================
main <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("Oligodendrocyte-ALS Mendelian Randomization\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  snp_pos <- load_snp_pos()

  # Try loading local ALS GWAS
  als_gwas <- load_als_gwas_local()

  # Load gene lists
  cat("\nLoading gene lists...\n")
  stable_high <- fread(STABLE_HIGH_FILE)
  darkgrey    <- fread(DARKGREY_FILE)
  cat(sprintf("  Stable-High: %d genes (symbols: %s...)\n",
              nrow(stable_high), paste(head(stable_high$symbol, 5), collapse = ", ")))
  cat(sprintf("  Darkgrey: %d genes\n", nrow(darkgrey)))

  # RTK/Trophic pathway genes (from darkgrey input analysis)
  rtk_genes <- c("NTRK2", "IGF1R", "INSR", "FYN", "DOCK3", "ROR1", "LRP4",
                  "ADAM10", "ADAM12", "MEGF9", "NEO1", "EPS8", "EPS15",
                  "RASA1", "RASGRF1", "RASGRF2", "RASGRP3", "MAP2K5",
                  "MAPK10", "PRKCA")

  all_results <- list()

  # === Track A: Oligo all eQTL ===
  cat("\n### TRACK A: Oligodendrocyte All eQTL -> ALS ###\n")
  oligo_eqtl <- load_celltype_eqtl("Oligodendrocytes", snp_pos)
  if (nrow(oligo_eqtl) >= 3) {
    oligo_exp <- to_exposure_format(oligo_eqtl, "Oligo_eQTL")
    all_results[["TrackA_Oligo"]] <- run_mr_analysis(oligo_exp, "TrackA_Oligo", als_gwas)
  }

  # === Track B: Stable-High genes ===
  cat("\n### TRACK B: Stable-High 135 Genes -> ALS ###\n")
  sh_symbols <- stable_high$symbol
  oligo_sh <- load_celltype_eqtl("Oligodendrocytes", snp_pos, gene_symbols = sh_symbols)
  if (nrow(oligo_sh) >= 3) {
    sh_exp <- to_exposure_format(oligo_sh, "Oligo_StableHigh_eQTL")
    all_results[["TrackB_StableHigh"]] <- run_mr_analysis(sh_exp, "TrackB_StableHigh", als_gwas)
  }

  # === Track D: Darkgrey genes ===
  cat("\n### TRACK D: Darkgrey 327 Genes -> ALS ###\n")
  dg_col <- intersect(names(darkgrey), c("symbol", "gene", "gene_symbol"))
  if (length(dg_col) > 0) {
    dg_symbols <- darkgrey[[dg_col[1]]]
  } else {
    dg_symbols <- darkgrey[[1]]
  }
  oligo_dg <- load_celltype_eqtl("Oligodendrocytes", snp_pos, gene_symbols = dg_symbols)
  if (nrow(oligo_dg) >= 3) {
    dg_exp <- to_exposure_format(oligo_dg, "Oligo_Darkgrey_eQTL")
    all_results[["TrackD_Darkgrey"]] <- run_mr_analysis(dg_exp, "TrackD_Darkgrey", als_gwas)
  }

  # === Track E: RTK/Trophic pathway ===
  cat("\n### TRACK E: RTK/Trophic Pathway -> ALS ###\n")
  oligo_rtk <- load_celltype_eqtl("Oligodendrocytes", snp_pos, gene_symbols = rtk_genes)
  if (nrow(oligo_rtk) >= 3) {
    rtk_exp <- to_exposure_format(oligo_rtk, "Oligo_RTK_eQTL")
    all_results[["TrackE_RTK"]] <- run_mr_analysis(rtk_exp, "TrackE_RTK", als_gwas)
  }

  # === Track C: All cell types (if data available) ===
  cat("\n### TRACK C: All Cell Types ###\n")
  cell_types <- c("Oligodendrocytes", "OPCs_COPs", "Astrocytes",
                  "Excitatory_neurons", "Inhibitory_neurons",
                  "Microglia", "Endothelial", "Pericytes")

  for (ct in cell_types) {
    ct_dir <- file.path(DATA_DIR, ct)
    if (!dir.exists(ct_dir)) {
      cat(sprintf("  %s: data not downloaded, skipping\n", ct))
      next
    }
    ct_eqtl <- load_celltype_eqtl(ct, snp_pos)
    if (nrow(ct_eqtl) >= 3) {
      ct_exp <- to_exposure_format(ct_eqtl, paste0(ct, "_eQTL"))
      all_results[[paste0("TrackC_", ct)]] <- run_mr_analysis(
        ct_exp, paste0("TrackC_", ct), als_gwas)
    }
  }

  # === Summary ===
  cat("\n\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SUMMARY\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  summary_df <- do.call(rbind, lapply(names(all_results), function(name) {
    r <- all_results[[name]]
    if (is.null(r)) return(NULL)
    ivw <- r$mr_results[r$mr_results$method == "Inverse variance weighted", ]
    if (nrow(ivw) == 0) return(NULL)
    data.frame(
      track = name, method = "IVW", n_snp = ivw$nsnp,
      beta = round(ivw$b, 4), se = round(ivw$se, 4),
      pval = signif(ivw$pval, 3),
      OR = round(exp(ivw$b), 3),
      OR_CI_low = round(exp(ivw$b - 1.96 * ivw$se), 3),
      OR_CI_high = round(exp(ivw$b + 1.96 * ivw$se), 3),
      egger_p = signif(r$pleiotropy$pval[1], 3),
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(summary_df) && nrow(summary_df) > 0) {
    cat("\n=== IVW Results ===\n")
    print(summary_df, row.names = FALSE)
    write.csv(summary_df, file.path(OUTPUT_DIR, "summary_all_tracks.csv"), row.names = FALSE)

    # Comparison forest plot
    p <- ggplot(summary_df, aes(x = reorder(track, beta), y = OR)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = OR_CI_low, ymax = OR_CI_high), width = 0.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = "Cell-type-specific MR: eQTL -> ALS risk",
           x = "", y = "Odds Ratio (95% CI)") +
      theme_minimal(base_size = 14)
    ggsave(file.path(OUTPUT_DIR, "comparison_forest.png"), p, width = 10, height = 6)
    cat("  Comparison forest plot saved.\n")
  }

  cat("\nPipeline complete.\n")
  return(all_results)
}

# Run
results <- main()
