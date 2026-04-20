#!/usr/bin/env Rscript
# ============================================================
# 05_smr_heidi.R
# SMR/HEIDI-like analysis using coloc
# (True SMR requires the SMR binary; we use coloc as the
#  Bayesian equivalent for eQTL-GWAS colocalization)
#
# For each gene with a significant eQTL in oligodendrocytes,
# test whether the eQTL and ALS GWAS signal colocalize
# (share a causal variant) at that locus.
# ============================================================

user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))

library(coloc)
library(data.table)
library(dplyr)

DATA_DIR      <- "D:/Projects/MR/data/bryois_eqtl"
ALS_GWAS_FILE <- "D:/Projects/MR/data/als_gwas/GCST90027164_buildGRCh37.tsv.gz"
OUTPUT_DIR    <- "D:/Projects/MR/results_coloc"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

STABLE_HIGH_FILE <- "D:/Projects/MR/sals_analysis_frozen_20260211/results/track_b/sals_upstream_gene_list/stable_high_genes.csv"

# ============================================================
# Load data
# ============================================================
cat("Loading SNP positions...\n")
snp_pos <- fread(file.path(DATA_DIR, "snp_pos.txt.gz"))
# Extract chr and pos from SNP_id_hg19
snp_pos[, chr := as.integer(sub("chr", "", sub(":.*", "", SNP_id_hg19)))]
snp_pos[, pos := as.integer(sub(".*:", "", SNP_id_hg19))]

cat("Loading ALS GWAS...\n")
als <- fread(ALS_GWAS_FILE)
setkey(als, rsid)

cat("Loading Oligo eQTL (all chr)...\n")
all_chr <- lapply(1:22, function(chr) {
  f <- file.path(DATA_DIR, "Oligodendrocytes",
                 sprintf("Oligodendrocytes.%d.gz", chr))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, header = FALSE)
  setnames(dt, c("gene_id", "snp_id", "dist_tss", "pvalue", "beta"))
  dt
})
eqtl <- rbindlist(all_chr[!sapply(all_chr, is.null)])
eqtl[, gene_symbol := sub("_ENSG.*", "", gene_id)]
cat(sprintf("  %d SNP-gene pairs, %d genes\n", nrow(eqtl), uniqueN(eqtl$gene_id)))

# Merge with positions
eqtl_pos <- merge(eqtl, snp_pos[, .(SNP, chr, pos, MAF)],
                  by.x = "snp_id", by.y = "SNP")

# Compute SE
eqtl_pos[, se := abs(beta) / qnorm(pvalue / 2, lower.tail = FALSE)]
eqtl_pos[!is.finite(se) | se <= 0, se := NA]

# Load gene lists
stable_high <- fread(STABLE_HIGH_FILE)
cat(sprintf("Stable-High genes: %d\n", nrow(stable_high)))

# ============================================================
# Colocalization per gene
# ============================================================
run_coloc_for_gene <- function(gene_sym, eqtl_data, als_data,
                                window = 500000, n_eqtl = 196) {
  # Get eQTL data for this gene
  gene_eqtl <- eqtl_data[gene_symbol == gene_sym & !is.na(se)]

  if (nrow(gene_eqtl) < 10) return(NULL)

  # Find the lead SNP
  lead <- gene_eqtl[which.min(pvalue)]
  lead_chr <- lead$chr
  lead_pos <- lead$pos

  # Define window around lead SNP
  gene_eqtl_window <- gene_eqtl[chr == lead_chr &
                                   pos >= (lead_pos - window) &
                                   pos <= (lead_pos + window)]

  if (nrow(gene_eqtl_window) < 10) return(NULL)

  # Get ALS data in same window
  als_window <- als_data[chromosome == lead_chr &
                           base_pair_location >= (lead_pos - window) &
                           base_pair_location <= (lead_pos + window)]

  # Find overlapping SNPs
  common_snps <- intersect(gene_eqtl_window$snp_id, als_window$rsid)
  if (length(common_snps) < 10) return(NULL)

  eqtl_sub <- gene_eqtl_window[snp_id %in% common_snps]
  als_sub  <- als_window[rsid %in% common_snps]

  # Deduplicate
  eqtl_sub <- eqtl_sub[!duplicated(snp_id)]
  als_sub  <- als_sub[!duplicated(rsid)]

  # Align
  eqtl_sub <- eqtl_sub[order(snp_id)]
  als_sub  <- als_sub[match(eqtl_sub$snp_id, rsid)]

  # Remove NAs
  keep <- !is.na(eqtl_sub$se) & !is.na(als_sub$standard_error) &
          eqtl_sub$se > 0 & als_sub$standard_error > 0
  eqtl_sub <- eqtl_sub[keep]
  als_sub  <- als_sub[keep]

  if (nrow(eqtl_sub) < 10) return(NULL)

  # coloc input
  dataset1 <- list(
    beta = eqtl_sub$beta,
    varbeta = eqtl_sub$se^2,
    snp = eqtl_sub$snp_id,
    position = eqtl_sub$pos,
    type = "quant",
    N = n_eqtl,
    MAF = eqtl_sub$MAF,
    sdY = 1
  )

  # ALS: case-control
  n_als <- median(als_sub$N_effective, na.rm = TRUE)
  dataset2 <- list(
    beta = als_sub$beta,
    varbeta = als_sub$standard_error^2,
    snp = als_sub$rsid,
    position = als_sub$base_pair_location,
    type = "cc",
    N = n_als,
    s = 0.5  # approximate case fraction
  )

  # Run coloc
  result <- tryCatch(
    coloc.abf(dataset1, dataset2),
    error = function(e) {
      cat(sprintf("    coloc error for %s: %s\n", gene_sym, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(result)) return(NULL)

  # Extract posterior probabilities
  pp <- result$summary
  list(
    gene = gene_sym,
    n_snps = nrow(eqtl_sub),
    lead_snp = lead$snp_id,
    lead_pval = lead$pvalue,
    chr = lead_chr,
    pos = lead_pos,
    PP.H0 = pp["PP.H0.abf"],
    PP.H1 = pp["PP.H1.abf"],
    PP.H2 = pp["PP.H2.abf"],
    PP.H3 = pp["PP.H3.abf"],
    PP.H4 = pp["PP.H4.abf"]
  )
}

# ============================================================
# Run coloc for top genes
# ============================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Colocalization Analysis (coloc.abf)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Get genes with at least one significant eQTL
gene_min_p <- eqtl_pos[, .(min_p = min(pvalue, na.rm = TRUE),
                             n_sig = sum(pvalue < 5e-6, na.rm = TRUE)),
                         by = gene_symbol]
gene_min_p <- gene_min_p[n_sig >= 1][order(min_p)]

cat(sprintf("Genes with >= 1 eQTL at p<5e-6: %d\n", nrow(gene_min_p)))
cat(sprintf("Of which in Stable-High: %d\n",
            sum(gene_min_p$gene_symbol %in% stable_high$symbol)))

# Run for top 100 genes (or all with sig eQTLs)
top_genes <- head(gene_min_p$gene_symbol, 100)
cat(sprintf("\nRunning coloc for top %d genes...\n\n", length(top_genes)))

coloc_results <- list()
for (i in seq_along(top_genes)) {
  g <- top_genes[i]
  if (i %% 10 == 0) cat(sprintf("  [%d/%d] %s\n", i, length(top_genes), g))
  res <- run_coloc_for_gene(g, eqtl_pos, als)
  if (!is.null(res)) coloc_results[[g]] <- res
}

cat(sprintf("\nCompleted: %d / %d genes with valid coloc results\n",
            length(coloc_results), length(top_genes)))

# ============================================================
# Results
# ============================================================
if (length(coloc_results) > 0) {
  coloc_df <- do.call(rbind, lapply(coloc_results, function(x) {
    data.frame(
      gene = x$gene, n_snps = x$n_snps,
      lead_snp = x$lead_snp, lead_pval = x$lead_pval,
      chr = x$chr, pos = x$pos,
      PP.H0 = round(x$PP.H0, 4),
      PP.H1 = round(x$PP.H1, 4),
      PP.H2 = round(x$PP.H2, 4),
      PP.H3 = round(x$PP.H3, 4),
      PP.H4 = round(x$PP.H4, 4),
      in_stable_high = x$gene %in% stable_high$symbol,
      stringsAsFactors = FALSE
    )
  }))
  coloc_df <- coloc_df[order(-coloc_df$PP.H4), ]

  cat("\n=== Top Colocalization Results (by PP.H4) ===\n")
  cat("PP.H4 > 0.5 suggests shared causal variant\n")
  cat("PP.H4 > 0.8 is strong evidence\n\n")

  # Show top 20
  print(head(coloc_df[, c("gene", "n_snps", "lead_snp", "PP.H3", "PP.H4",
                            "in_stable_high")], 20), row.names = FALSE)

  # Highlight genes with PP.H4 > 0.5
  strong <- coloc_df[coloc_df$PP.H4 > 0.5, ]
  if (nrow(strong) > 0) {
    cat(sprintf("\n=== Genes with PP.H4 > 0.5 (%d) ===\n", nrow(strong)))
    print(strong[, c("gene", "chr", "lead_snp", "PP.H4", "in_stable_high")],
          row.names = FALSE)
  } else {
    cat("\nNo genes with PP.H4 > 0.5\n")
  }

  suggestive <- coloc_df[coloc_df$PP.H4 > 0.3, ]
  if (nrow(suggestive) > 0) {
    cat(sprintf("\n=== Genes with PP.H4 > 0.3 (suggestive, %d) ===\n", nrow(suggestive)))
    print(suggestive[, c("gene", "chr", "lead_snp", "PP.H4", "in_stable_high")],
          row.names = FALSE)
  }

  # Save
  write.csv(coloc_df, file.path(OUTPUT_DIR, "coloc_results.csv"), row.names = FALSE)
  saveRDS(coloc_results, file.path(OUTPUT_DIR, "coloc_full.rds"))

  # Summary stats
  cat(sprintf("\n=== Summary ===\n"))
  cat(sprintf("Total genes tested: %d\n", nrow(coloc_df)))
  cat(sprintf("PP.H4 > 0.8 (strong): %d\n", sum(coloc_df$PP.H4 > 0.8)))
  cat(sprintf("PP.H4 > 0.5 (moderate): %d\n", sum(coloc_df$PP.H4 > 0.5)))
  cat(sprintf("PP.H4 > 0.3 (suggestive): %d\n", sum(coloc_df$PP.H4 > 0.3)))
  cat(sprintf("PP.H3 > 0.5 (distinct signals): %d\n", sum(coloc_df$PP.H3 > 0.5)))
  cat(sprintf("Median PP.H4: %.4f\n", median(coloc_df$PP.H4)))
}

cat("\nColoc analysis complete.\n")
