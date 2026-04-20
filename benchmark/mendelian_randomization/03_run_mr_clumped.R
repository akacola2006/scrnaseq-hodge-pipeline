#!/usr/bin/env Rscript
# ============================================================
# 03_run_mr_clumped.R
# MR with LOCAL LD clumping (PLINK + 1000G EUR)
# Re-runs Track A-E with proper clumping
# ============================================================

user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))

library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(ggplot2)

# --- Paths ---
DATA_DIR      <- "D:/Projects/MR/data/bryois_eqtl"
ALS_GWAS_FILE <- "D:/Projects/MR/data/als_gwas/GCST90027164_buildGRCh37.tsv.gz"
OUTPUT_DIR    <- "D:/Projects/MR/results_clumped"
PLINK_PATH    <- "D:/Projects/MR/tools/plink.exe"
REF_DIR       <- "D:/Projects/MR/data/ld_reference"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

STABLE_HIGH_FILE <- "D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/sals_upstream_gene_list/stable_high_genes.csv"
DARKGREY_FILE    <- "D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/laneB/darkgrey_identity/darkgrey_genes.tsv"

# --- Detect EUR reference bfile ---
detect_bfile <- function() {
  candidates <- c(
    file.path(REF_DIR, "EUR"),
    file.path(REF_DIR, "1kg.v3", "EUR"),
    file.path(REF_DIR, "1kg_v3", "EUR")
  )
  for (c in candidates) {
    if (file.exists(paste0(c, ".bed"))) {
      cat(sprintf("  EUR reference found: %s\n", c))
      return(c)
    }
  }
  # Search recursively
  bed_files <- list.files(REF_DIR, pattern = "EUR\\.bed$", recursive = TRUE, full.names = TRUE)
  if (length(bed_files) > 0) {
    bfile <- sub("\\.bed$", "", bed_files[1])
    cat(sprintf("  EUR reference found: %s\n", bfile))
    return(bfile)
  }
  stop("EUR reference panel not found in ", REF_DIR)
}

# ============================================================
# Local LD clumping
# ============================================================
local_clump <- function(exposure_dat, bfile, plink_path = PLINK_PATH) {
  cat("  Local LD clumping (r2<0.001, kb=10000)...\n")

  # Prepare data for ieugwasr::ld_clump_local
  clump_input <- data.frame(
    rsid = exposure_dat$SNP,
    pval = exposure_dat$pval.exposure,
    id = exposure_dat$id.exposure,
    stringsAsFactors = FALSE
  )

  clumped_rsids <- tryCatch({
    result <- ld_clump(
      dat = clump_input,
      clump_r2 = 0.001,
      clump_kb = 10000,
      plink_bin = plink_path,
      bfile = bfile
    )
    result$rsid
  }, error = function(e) {
    cat(sprintf("  Local clump error: %s\n", conditionMessage(e)))
    cat("  Trying ld_clump_local...\n")
    tryCatch({
      result <- ld_clump_local(
        dat = clump_input,
        clump_r2 = 0.001,
        clump_kb = 10000,
        plink_bin = plink_path,
        bfile = bfile
      )
      result$rsid
    }, error = function(e2) {
      cat(sprintf("  ld_clump_local also failed: %s\n", conditionMessage(e2)))
      NULL
    })
  })

  if (is.null(clumped_rsids) || length(clumped_rsids) == 0) {
    cat("  WARNING: Clumping failed. Returning original data.\n")
    return(exposure_dat)
  }

  exposure_clumped <- exposure_dat[exposure_dat$SNP %in% clumped_rsids, ]
  cat(sprintf("  Clumped: %d -> %d independent instruments\n",
              nrow(exposure_dat), nrow(exposure_clumped)))
  return(exposure_clumped)
}

# ============================================================
# Load SNP positions
# ============================================================
load_snp_pos <- function() {
  cat("Loading SNP positions...\n")
  snp_pos <- fread(file.path(DATA_DIR, "snp_pos.txt.gz"))
  cat(sprintf("  %d SNPs\n", nrow(snp_pos)))
  return(snp_pos)
}

# ============================================================
# Load eQTL for one cell type
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

  eqtl[, gene_symbol := sub("_ENSG.*", "", gene_id)]

  cat(sprintf("  Raw: %d pairs, %d genes\n", nrow(eqtl), uniqueN(eqtl$gene_id)))

  if (!is.null(gene_symbols)) {
    eqtl <- eqtl[gene_symbol %in% gene_symbols]
    cat(sprintf("  After gene filter: %d pairs, %d genes\n",
                nrow(eqtl), uniqueN(eqtl$gene_symbol)))
  }

  eqtl_m <- merge(eqtl, snp_pos, by.x = "snp_id", by.y = "SNP")
  eqtl_m[, se := abs(beta) / qnorm(pvalue / 2, lower.tail = FALSE)]
  eqtl_m[!is.finite(se) | se <= 0, se := NA]

  eqtl_sig <- eqtl_m[pvalue < 5e-8 & !is.na(se)]
  cat(sprintf("  GW-significant: %d\n", nrow(eqtl_sig)))

  if (nrow(eqtl_sig) < 3) {
    eqtl_sig <- eqtl_m[pvalue < 5e-6 & !is.na(se)]
    cat(sprintf("  Suggestive (p<5e-6): %d\n", nrow(eqtl_sig)))
  }

  eqtl_sig[, F_stat := (beta / se)^2]
  eqtl_sig <- eqtl_sig[F_stat > 10]
  cat(sprintf("  After F>10: %d (median F=%.1f)\n",
              nrow(eqtl_sig), median(eqtl_sig$F_stat)))

  return(eqtl_sig)
}

# ============================================================
# To exposure format
# ============================================================
to_exposure_format <- function(eqtl_sig, label) {
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

  exposure <- exposure %>%
    group_by(SNP) %>%
    slice_min(pval.exposure, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    as.data.frame()

  return(exposure)
}

# ============================================================
# Format ALS outcome
# ============================================================
format_als_outcome <- function(als_gwas, snps) {
  als_filt <- als_gwas[rsid %in% snps]
  cat(sprintf("  ALS outcome: %d / %d SNPs matched\n", nrow(als_filt), length(snps)))
  if (nrow(als_filt) == 0) return(NULL)

  data.frame(
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
}

# ============================================================
# MR analysis
# ============================================================
run_mr_analysis <- function(exposure_dat, label, als_gwas, bfile,
                            output_dir = OUTPUT_DIR) {
  cat(sprintf("\n========== MR: %s ==========\n", label))
  result_dir <- file.path(output_dir, label)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  # LOCAL LD clumping
  exposure_clumped <- local_clump(exposure_dat, bfile)

  if (nrow(exposure_clumped) < 3) {
    cat(sprintf("  SKIP: <3 instruments for %s\n", label))
    saveRDS(list(status = "insufficient", n = nrow(exposure_clumped)),
            file.path(result_dir, "result.rds"))
    return(NULL)
  }

  # Outcome
  outcome_dat <- format_als_outcome(als_gwas, exposure_clumped$SNP)
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    cat("  SKIP: No outcome data.\n")
    return(NULL)
  }

  # Harmonize
  dat <- harmonise_data(exposure_clumped, outcome_dat)
  dat <- dat[dat$mr_keep == TRUE, ]
  cat(sprintf("  Harmonized: %d SNPs\n", nrow(dat)))

  if (nrow(dat) < 3) {
    cat("  SKIP: <3 harmonized.\n")
    return(NULL)
  }

  # MR
  mr_res <- mr(dat, method_list = c(
    "mr_ivw", "mr_egger_regression",
    "mr_weighted_median", "mr_weighted_mode"))

  # Sensitivity
  pleio <- mr_pleiotropy_test(dat)
  hetero <- mr_heterogeneity(dat)
  single <- mr_singlesnp(dat)
  loo <- mr_leaveoneout(dat)

  # Print
  cat(sprintf("\n  --- %s Results (CLUMPED) ---\n", label))
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
  }, error = function(e) cat(sprintf("  Plot error: %s\n", conditionMessage(e))))

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
  cat("MR with LOCAL LD Clumping (PLINK + 1000G EUR)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # Verify PLINK
  if (!file.exists(PLINK_PATH)) stop("PLINK not found: ", PLINK_PATH)
  cat(sprintf("PLINK: %s\n", PLINK_PATH))

  # Detect reference
  bfile <- detect_bfile()

  # Load data
  snp_pos <- load_snp_pos()
  als_gwas <- fread(ALS_GWAS_FILE)
  cat(sprintf("ALS GWAS: %d SNPs\n", nrow(als_gwas)))

  stable_high <- fread(STABLE_HIGH_FILE)
  darkgrey    <- fread(DARKGREY_FILE)

  rtk_genes <- c("NTRK2", "IGF1R", "INSR", "FYN", "DOCK3", "ROR1", "LRP4",
                  "ADAM10", "ADAM12", "MEGF9", "NEO1", "EPS8", "EPS15",
                  "RASA1", "RASGRF1", "RASGRF2", "RASGRP3", "MAP2K5",
                  "MAPK10", "PRKCA")

  all_results <- list()

  # Track A
  cat("\n### TRACK A: Oligo All eQTL ###\n")
  oligo_eqtl <- load_celltype_eqtl("Oligodendrocytes", snp_pos)
  oligo_exp <- to_exposure_format(oligo_eqtl, "Oligo_eQTL")
  all_results[["TrackA_Oligo"]] <- run_mr_analysis(oligo_exp, "TrackA_Oligo", als_gwas, bfile)

  # Track B
  cat("\n### TRACK B: Stable-High 135 ###\n")
  oligo_sh <- load_celltype_eqtl("Oligodendrocytes", snp_pos, gene_symbols = stable_high$symbol)
  if (nrow(oligo_sh) >= 3) {
    sh_exp <- to_exposure_format(oligo_sh, "Oligo_StableHigh_eQTL")
    all_results[["TrackB_StableHigh"]] <- run_mr_analysis(sh_exp, "TrackB_StableHigh", als_gwas, bfile)
  }

  # Track D
  cat("\n### TRACK D: Darkgrey 327 ###\n")
  dg_col <- intersect(names(darkgrey), c("symbol", "gene", "gene_symbol"))
  dg_symbols <- if (length(dg_col) > 0) darkgrey[[dg_col[1]]] else darkgrey[[1]]
  oligo_dg <- load_celltype_eqtl("Oligodendrocytes", snp_pos, gene_symbols = dg_symbols)
  if (nrow(oligo_dg) >= 3) {
    dg_exp <- to_exposure_format(oligo_dg, "Oligo_Darkgrey_eQTL")
    all_results[["TrackD_Darkgrey"]] <- run_mr_analysis(dg_exp, "TrackD_Darkgrey", als_gwas, bfile)
  }

  # Track E
  cat("\n### TRACK E: RTK/Trophic ###\n")
  oligo_rtk <- load_celltype_eqtl("Oligodendrocytes", snp_pos, gene_symbols = rtk_genes)
  if (nrow(oligo_rtk) >= 3) {
    rtk_exp <- to_exposure_format(oligo_rtk, "Oligo_RTK_eQTL")
    all_results[["TrackE_RTK"]] <- run_mr_analysis(rtk_exp, "TrackE_RTK", als_gwas, bfile)
  }

  # Summary
  cat("\n\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SUMMARY (WITH LD CLUMPING)\n")
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
    cat("\n=== IVW Results (Clumped) ===\n")
    print(summary_df, row.names = FALSE)
    write.csv(summary_df, file.path(OUTPUT_DIR, "summary_clumped.csv"), row.names = FALSE)

    p <- ggplot(summary_df, aes(x = reorder(track, OR), y = OR)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = OR_CI_low, ymax = OR_CI_high), width = 0.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = "MR: Oligo eQTL -> ALS (LD clumped)",
           x = "", y = "Odds Ratio (95% CI)") +
      theme_minimal(base_size = 14)
    ggsave(file.path(OUTPUT_DIR, "comparison_forest_clumped.png"), p, width = 10, height = 6)
  }

  cat("\nPipeline complete (with LD clumping).\n")
}

main()
