#!/usr/bin/env Rscript
# ============================================================
# 08_mam_pathway_mr.R
# Track G: MAM (Mitochondria-Associated Membrane) pathway
# Oligodendrocyte eQTL (Bryois 2022) -> ALS GWAS (van Rheenen 2021)
#
# Gene set rationale:
#   Motivated by all-gene phi insertion analysis (Appendix AL, AN):
#   - STX17/SPTLC1/ATG14/EIF2AK3/SNW1 are the ONLY 5 genes
#     in the top 25% across ALL 5 cell types simultaneously
#   - All 5 are MAM core functions or MAM-associated stress responses
#   - STX17 is the sole gene in top 15% across ALL 5 CTs
#   - GDI2 is the convergent gene: exome P=9.4e-4 AND phi top 5%
#
# Sub-tracks:
#   G1_PhiTop5   - 5 CT-universal phi-upstream genes
#   G2_Tethering - MAM structural tethering
#   G3_Autophagy - MAM-localised autophagy
#   G4_Lipid     - MAM lipid synthesis / ceramide
#   G5_All       - All MAM gene sets combined
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
OUTPUT_DIR    <- "D:/Projects/MR/results_mam"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PTHRESH <- 5e-8
PTHRESH_RELAX <- 5e-6

# ============================================================
# MAM gene sets
# ============================================================

# G1: CT-universal phi-upstream genes (sole 5 genes top 25% all 5 CTs)
#     From all-gene insertion analysis (Appendix AL / AN)
mam_phi_top5 <- c(
  "STX17",    # Autophagosome-lysosome SNARE (Oligo 1.4%, L4/6 2.3%, Astro 4.5%, L2/3 1.3%, OPC 6.6%)
  "SPTLC1",   # Ceramide synthesis control / ORMDL complex (GoF -> juvenile ALS)
  "ATG14",    # PI3K complex, autophagy initiation at MAM
  "EIF2AK3",  # PERK: MAM-localised ER stress sensor -> eIF2alpha-P -> global translation repression
  "SNW1"      # Splicing factor / ALS genetic risk
)

# G2: MAM structural tethering
mam_tethering <- c(
  "VAPB",     # MAM tethering (VAPB-PTPIP51 contact), ALS8 locus
  "RMDN3",    # PTPIP51: VAPB interaction partner
  "MFN2",     # Mitofusin-2: ER-mito tethering
  "PACS2",    # MAM organiser
  "SIGMAR1",  # Sigma-1 receptor: MAM chaperone, binds TDP-43 / FUS
  "HSPA9"     # GRP75/mortalin: IP3R-VDAC bridge at MAM
)

# G3: MAM-localised autophagy
mam_autophagy <- c(
  "STX17",    # (also in G1) Autophagosome fusion SNARE
  "ATG14",    # (also in G1) PI3K/Beclin complex
  "SQSTM1",   # p62 autophagy receptor (Oligo 2.8%, Astro 5.1%)
  "TBK1",     # Autophagy kinase (OPTN/SQSTM1 phosphorylation)
  "OPTN",     # Optineurin: autophagy receptor
  "BECN1",    # Beclin-1: autophagy initiation
  "PIK3C3",   # VPS34: PI3P production for autophagosome
  "RUBCN",    # Rubicon: autophagosome-lysosome fusion regulator
  "RAB7A"     # Late endosome/autophagosome maturation
)

# G4: MAM lipid synthesis and ceramide metabolism
mam_lipid <- c(
  "SPTLC1",   # (also in G1) Serine palmitoyltransferase subunit 1
  "SPTLC2",   # Serine palmitoyltransferase catalytic subunit
  "CERT1",    # COL4A3BP: ceramide transfer protein at MAM
  "ACSL4",    # Acyl-CoA synthetase: arachidonic acid -> ER-mito lipid transfer
  "CERS2",    # Ceramide synthase 2: myelin-enriched ceramide
  "DEGS1",    # Sphingolipid delta desaturase
  "GBA"       # Glucocerebrosidase: sphingolipid recycling
)

# G5: Ca2+ transfer at MAM
mam_ca2 <- c(
  "ITPR2",    # IP3R2: ER Ca2+ release channel at MAM (Oligo 3.1%)
  "VDAC1",    # Voltage-dependent anion channel: mitochondrial Ca2+ entry
  "MCU",      # Mitochondrial Ca2+ uniporter
  "MICU1",    # MCU gatekeeper (neuron upstream)
  "MICU3"     # MCU enhancer (neuron upstream, L4/6 0.4%)
)

# Exome convergence (GDI2 + USP10)
mam_exome_convergent <- c(
  "GDI2",     # GDP dissociation inhibitor 2: Rho GTPase GDI (exome P=9.4e-4, phi Astro 0.7%)
  "USP10"     # Deubiquitinase (exome FDR<0.05, Oligo phi 4.8%)
)

# Combined unique gene list
mam_all <- unique(c(mam_phi_top5, mam_tethering, mam_autophagy,
                    mam_lipid, mam_ca2, mam_exome_convergent))

cat("MAM gene sets:\n")
cat(sprintf("  G1 PhiTop5    (%2d): %s\n", length(mam_phi_top5), paste(mam_phi_top5, collapse=", ")))
cat(sprintf("  G2 Tethering  (%2d): %s\n", length(mam_tethering), paste(mam_tethering, collapse=", ")))
cat(sprintf("  G3 Autophagy  (%2d): %s\n", length(mam_autophagy), paste(mam_autophagy, collapse=", ")))
cat(sprintf("  G4 Lipid      (%2d): %s\n", length(mam_lipid), paste(mam_lipid, collapse=", ")))
cat(sprintf("  G5 Ca2+       (%2d): %s\n", length(mam_ca2), paste(mam_ca2, collapse=", ")))
cat(sprintf("  Exome conv.   (%2d): %s\n", length(mam_exome_convergent), paste(mam_exome_convergent, collapse=", ")))
cat(sprintf("  G_All (union) (%2d): %s\n\n", length(mam_all), paste(mam_all, collapse=", ")))

# ============================================================
# Load data
# ============================================================
cat("Loading SNP positions...\n")
snp_pos <- fread(file.path(DATA_DIR, "snp_pos.txt.gz"))
snp_pos[, chr := as.integer(sub("chr", "", sub(":.*", "", SNP_id_hg19)))]
snp_pos[, pos := as.integer(sub(".*:", "", SNP_id_hg19))]

cat("Loading ALS GWAS...\n")
als <- fread(ALS_GWAS_FILE)
cat(sprintf("  ALS GWAS: %d variants\n", nrow(als)))

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

eqtl_full <- merge(eqtl_raw,
                   snp_pos[, .(SNP, chr, pos, MAF, effect_allele, other_allele)],
                   by.x = "snp_id", by.y = "SNP")
eqtl_full[, se := abs(beta) / qnorm(pvalue / 2, lower.tail = FALSE)]
eqtl_full[!is.finite(se) | se <= 0, se := NA]
cat(sprintf("  Total eQTL: %d pairs\n\n", nrow(eqtl_full)))

# ============================================================
# Gene coverage check
# ============================================================
cat("=== MAM Gene Coverage in Bryois Oligo eQTL ===\n")

coverage <- data.frame(
  gene      = mam_all,
  in_g1     = mam_all %in% mam_phi_top5,
  in_g2     = mam_all %in% mam_tethering,
  in_g3     = mam_all %in% mam_autophagy,
  in_g4     = mam_all %in% mam_lipid,
  in_g5     = mam_all %in% mam_ca2,
  in_exome  = mam_all %in% mam_exome_convergent,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(coverage))) {
  g     <- coverage$gene[i]
  g_dat <- eqtl_full[gene_symbol == g]
  coverage$n_snps[i]    <- nrow(g_dat)
  coverage$n_sig_5e8[i] <- sum(g_dat$pvalue < 5e-8, na.rm = TRUE)
  coverage$n_sig_5e6[i] <- sum(g_dat$pvalue < 5e-6, na.rm = TRUE)
  coverage$min_p[i]     <- if (nrow(g_dat) > 0) min(g_dat$pvalue, na.rm = TRUE) else NA
}

cat("\n")
print(coverage[, c("gene","in_g1","n_snps","n_sig_5e8","n_sig_5e6","min_p")], row.names = FALSE)
write.csv(coverage, file.path(OUTPUT_DIR, "mam_gene_coverage.csv"), row.names = FALSE)

# ============================================================
# Shared MR function (identical to Track F pattern)
# ============================================================
format_exposure <- function(eqtl_sub, label) {
  eqtl_sub[, F_stat := (beta / se)^2]
  eqtl_sub <- eqtl_sub[!is.na(se) & F_stat > 10]
  if (nrow(eqtl_sub) == 0) return(NULL)
  data.frame(
    SNP                    = eqtl_sub$snp_id,
    beta.exposure          = eqtl_sub$beta,
    se.exposure            = eqtl_sub$se,
    pval.exposure          = eqtl_sub$pvalue,
    effect_allele.exposure = eqtl_sub$effect_allele,
    other_allele.exposure  = eqtl_sub$other_allele,
    eaf.exposure           = eqtl_sub$MAF,
    exposure               = label, id.exposure = label,
    gene.exposure          = eqtl_sub$gene_symbol,
    chr.exposure           = eqtl_sub$chr,
    pos.exposure           = eqtl_sub$pos,
    mr_keep.exposure       = TRUE,
    stringsAsFactors = FALSE
  ) %>%
    group_by(SNP) %>%
    slice_min(pval.exposure, n = 1, with_ties = FALSE) %>%
    ungroup() %>% as.data.frame()
}

format_outcome <- function(als_filt) {
  data.frame(
    SNP                   = als_filt$rsid,
    beta.outcome          = als_filt$beta,
    se.outcome            = als_filt$standard_error,
    pval.outcome          = als_filt$p_value,
    effect_allele.outcome = toupper(als_filt$effect_allele),
    other_allele.outcome  = toupper(als_filt$other_allele),
    eaf.outcome           = als_filt$effect_allele_frequency,
    outcome               = "ALS", id.outcome = "ALS",
    mr_keep.outcome       = TRUE,
    samplesize.outcome    = als_filt$N_effective,
    stringsAsFactors = FALSE
  )
}

run_mam_track <- function(gene_set, label, eqtl_data, als_data,
                           pthresh = PTHRESH, pthresh_relax = PTHRESH_RELAX) {
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse="")))
  cat(sprintf("Track %s  |  %d genes\n", label, length(gene_set)))
  cat(sprintf("%s\n", paste(rep("=", 60), collapse="")))

  rdir <- file.path(OUTPUT_DIR, label)
  dir.create(rdir, showWarnings = FALSE, recursive = TRUE)

  # Filter eQTLs
  sub <- eqtl_data[gene_symbol %in% gene_set & pvalue < pthresh & !is.na(se)]
  if (nrow(sub) == 0) {
    cat(sprintf("  No instruments at p<%s. Relaxing to p<%s...\n",
                format(pthresh, scientific=TRUE), format(pthresh_relax, scientific=TRUE)))
    sub <- eqtl_data[gene_symbol %in% gene_set & pvalue < pthresh_relax & !is.na(se)]
    pthresh_used <- pthresh_relax
  } else {
    pthresh_used <- pthresh
  }

  cat(sprintf("  Detected genes: %s\n",
              paste(unique(sub$gene_symbol), collapse=", ")))
  cat(sprintf("  Pre-clump instruments (p<%s): %d from %d genes\n",
              format(pthresh_used, scientific=TRUE), nrow(sub),
              uniqueN(sub$gene_symbol)))

  exp_dat <- format_exposure(sub, label)
  if (is.null(exp_dat) || nrow(exp_dat) == 0) {
    cat("  No valid instruments (F>10). SKIP.\n")
    return(NULL)
  }

  # Clump
  exp_c <- tryCatch(
    clump_data(exp_dat, clump_r2 = 0.001, clump_kb = 10000),
    error = function(e) {
      cat(sprintf("  Clump error: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(exp_c) || nrow(exp_c) == 0) {
    cat("  Clumping returned 0 instruments. SKIP.\n")
    return(NULL)
  }
  cat(sprintf("  Post-clump: %d independent instruments\n", nrow(exp_c)))
  cat(sprintf("  Genes with IVs: %s\n",
              paste(unique(exp_c$gene.exposure), collapse=", ")))
  mean_F <- mean((exp_c$beta.exposure / exp_c$se.exposure)^2)
  cat(sprintf("  Mean F-statistic: %.1f\n", mean_F))

  # Outcome
  als_filt <- als_data[rsid %in% exp_c$SNP]
  if (nrow(als_filt) == 0) { cat("  No ALS GWAS variants matched.\n"); return(NULL) }
  out <- format_outcome(als_filt)

  dat <- harmonise_data(exp_c, out)
  dat <- dat[dat$mr_keep, ]
  cat(sprintf("  Harmonised: %d\n", nrow(dat)))
  if (nrow(dat) == 0) return(NULL)

  # MR
  if (nrow(dat) == 1) {
    res <- mr(dat, method_list = "mr_wald_ratio")
  } else if (nrow(dat) == 2) {
    res <- mr(dat, method_list = "mr_ivw")
  } else {
    res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                    "mr_weighted_median", "mr_weighted_mode"))
  }

  # Sensitivity
  pleio  <- if (nrow(dat) >= 3) mr_pleiotropy_test(dat) else data.frame(pval = NA)
  hetero <- if (nrow(dat) >= 3) mr_heterogeneity(dat) else data.frame(Q = NA, Q_pval = NA)
  single <- if (nrow(dat) >= 2) mr_singlesnp(dat) else NULL
  loo    <- if (nrow(dat) >= 3) mr_leaveoneout(dat) else NULL

  cat(sprintf("\n--- %s Results ---\n", label))
  print(res[, c("method", "nsnp", "b", "se", "pval")])

  if (nrow(dat) >= 3) {
    ivw_row <- res[grepl("Inverse variance", res$method), ]
    if (nrow(ivw_row) > 0) {
      cat(sprintf("\n  IVW OR: %.3f (%.3f-%.3f), p=%.4f\n",
                  exp(ivw_row$b[1]),
                  exp(ivw_row$b[1] - 1.96 * ivw_row$se[1]),
                  exp(ivw_row$b[1] + 1.96 * ivw_row$se[1]),
                  ivw_row$pval[1]))
      cat(sprintf("  Egger intercept p: %.4f\n",
                  ifelse(is.na(pleio$pval[1]), NA, round(pleio$pval[1], 4))))
      ivw_het <- hetero[grepl("Inverse variance", hetero$method), ]
      if (nrow(ivw_het) > 0)
        cat(sprintf("  Cochran Q p: %.4f\n", round(ivw_het$Q_pval[1], 4)))
    }
  } else if (nrow(dat) == 1) {
    wr <- res[res$method == "Wald ratio", ]
    cat(sprintf("\n  Wald ratio OR: %.3f (%.3f-%.3f), p=%.4f\n",
                exp(wr$b), exp(wr$b - 1.96 * wr$se),
                exp(wr$b + 1.96 * wr$se), wr$pval))
  }

  # Power
  n_out  <- median(als_data$N_effective, na.rm = TRUE)
  n_exp  <- 196
  R2     <- (nrow(dat) * mean_F) / (nrow(dat) * mean_F + n_exp)
  se_pwr <- 1 / sqrt(n_out * R2)
  min_b  <- (qnorm(0.975) + qnorm(0.80)) * se_pwr
  cat(sprintf("  Min detectable OR (80%% power): %.3f\n", exp(min_b)))

  # Plots
  if (!is.null(single) && nrow(dat) >= 3) {
    tryCatch({
      ggsave(file.path(rdir, "scatter.png"), mr_scatter_plot(res, dat)[[1]], width=8, height=6)
      ggsave(file.path(rdir, "forest.png"),  mr_forest_plot(single)[[1]],   width=8, height=max(6, nrow(dat)*0.4))
      ggsave(file.path(rdir, "funnel.png"),  mr_funnel_plot(single)[[1]],   width=8, height=6)
      if (!is.null(loo))
        ggsave(file.path(rdir, "loo.png"), mr_leaveoneout_plot(loo)[[1]], width=8, height=max(6, nrow(dat)*0.4))
    }, error = function(e) cat(sprintf("  Plot err: %s\n", conditionMessage(e))))
  }

  obj <- list(mr_results=res, dat=dat, pleiotropy=pleio, heterogeneity=hetero,
              single=single, loo=loo, n=nrow(dat), label=label, mean_F=mean_F,
              genes_with_iv=unique(exp_c$gene.exposure), pthresh=pthresh_used)
  saveRDS(obj, file.path(rdir, "result.rds"))
  write.csv(res, file.path(rdir, "mr_results.csv"), row.names = FALSE)
  return(obj)
}

# ============================================================
# Per-gene MR (individual gene Wald ratio or IVW)
# ============================================================
run_per_gene <- function(gene_sym, eqtl_data, als_data,
                          pthresh = PTHRESH, relax = PTHRESH_RELAX) {
  g <- eqtl_data[gene_symbol == gene_sym]
  if (nrow(g) == 0)
    return(data.frame(gene=gene_sym, detected=FALSE, n_iv=0, OR=NA, OR_lo=NA, OR_hi=NA,
                      pval=NA, method=NA, stringsAsFactors=FALSE))

  g_sig <- g[pvalue < pthresh & !is.na(se)]
  pt_used <- pthresh
  if (nrow(g_sig) == 0) {
    g_sig <- g[pvalue < relax & !is.na(se)]
    pt_used <- relax
  }
  g_sig[, F_stat := (beta / se)^2]
  g_sig <- g_sig[F_stat > 10]

  if (nrow(g_sig) == 0)
    return(data.frame(gene=gene_sym, detected=TRUE, n_iv=0, OR=NA, OR_lo=NA, OR_hi=NA,
                      pval=NA, method=paste0("no_iv_at_", format(relax, scientific=TRUE)),
                      stringsAsFactors=FALSE))

  exp <- data.frame(
    SNP=g_sig$snp_id, beta.exposure=g_sig$beta, se.exposure=g_sig$se,
    pval.exposure=g_sig$pvalue, effect_allele.exposure=g_sig$effect_allele,
    other_allele.exposure=g_sig$other_allele, eaf.exposure=g_sig$MAF,
    exposure=gene_sym, id.exposure=gene_sym, gene.exposure=gene_sym,
    mr_keep.exposure=TRUE, stringsAsFactors=FALSE
  ) %>% group_by(SNP) %>% slice_min(pval.exposure, n=1, with_ties=FALSE) %>%
    ungroup() %>% as.data.frame()

  exp_c <- tryCatch(clump_data(exp, clump_r2=0.001, clump_kb=10000), error=function(e) NULL)
  if (is.null(exp_c) || nrow(exp_c) == 0)
    return(data.frame(gene=gene_sym, detected=TRUE, n_iv=0, OR=NA, OR_lo=NA, OR_hi=NA,
                      pval=NA, method="clump_fail", stringsAsFactors=FALSE))

  als_f <- als_data[rsid %in% exp_c$SNP]
  if (nrow(als_f) == 0)
    return(data.frame(gene=gene_sym, detected=TRUE, n_iv=nrow(exp_c), OR=NA, OR_lo=NA, OR_hi=NA,
                      pval=NA, method="no_outcome", stringsAsFactors=FALSE))

  out <- format_outcome(als_f)
  dat <- tryCatch(harmonise_data(exp_c, out), error=function(e) NULL)
  if (is.null(dat)) return(data.frame(gene=gene_sym, detected=TRUE, n_iv=nrow(exp_c),
                                       OR=NA, OR_lo=NA, OR_hi=NA, pval=NA, method="harmonize_fail",
                                       stringsAsFactors=FALSE))
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) == 0) return(data.frame(gene=gene_sym, detected=TRUE, n_iv=nrow(exp_c),
                                         OR=NA, OR_lo=NA, OR_hi=NA, pval=NA, method="no_harmonized",
                                         stringsAsFactors=FALSE))

  method_list <- if (nrow(dat) == 1) "mr_wald_ratio" else "mr_ivw"
  res <- tryCatch(mr(dat, method_list=method_list), error=function(e) NULL)
  if (is.null(res)) return(data.frame(gene=gene_sym, detected=TRUE, n_iv=nrow(dat),
                                       OR=NA, OR_lo=NA, OR_hi=NA, pval=NA, method="mr_fail",
                                       stringsAsFactors=FALSE))

  data.frame(gene=gene_sym, detected=TRUE, n_iv=nrow(dat),
             OR     = exp(res$b[1]),
             OR_lo  = exp(res$b[1] - 1.96 * res$se[1]),
             OR_hi  = exp(res$b[1] + 1.96 * res$se[1]),
             pval   = res$pval[1],
             method = res$method[1],
             pthresh_used = format(pt_used, scientific=TRUE),
             stringsAsFactors=FALSE)
}

# ============================================================
# Colocalization
# ============================================================
run_coloc <- function(gene_sym, eqtl_data, als_data,
                       window=500000, n_eqtl=196) {
  g <- eqtl_data[gene_symbol == gene_sym & !is.na(se)]
  if (nrow(g) < 10) return(NULL)
  lead <- g[which.min(pvalue)]
  g_w  <- g[chr == lead$chr & pos >= (lead$pos - window) & pos <= (lead$pos + window)]
  a_w  <- als_data[chromosome == lead$chr &
                     base_pair_location >= (lead$pos - window) &
                     base_pair_location <= (lead$pos + window)]
  common <- intersect(g_w$snp_id, a_w$rsid)
  if (length(common) < 10) return(NULL)
  e <- g_w[snp_id %in% common][!duplicated(snp_id)][order(snp_id)]
  a <- a_w[rsid %in% common][!duplicated(rsid)]
  a <- a[match(e$snp_id, rsid)]
  keep <- !is.na(e$se) & !is.na(a$standard_error) & e$se > 0 & a$standard_error > 0
  e <- e[keep]; a <- a[keep]
  if (nrow(e) < 10) return(NULL)
  d1 <- list(beta=e$beta, varbeta=e$se^2, snp=e$snp_id, position=e$pos,
             type="quant", N=n_eqtl, MAF=e$MAF, sdY=1)
  d2 <- list(beta=a$beta, varbeta=a$standard_error^2, snp=a$rsid,
             position=a$base_pair_location, type="cc",
             N=median(a$N_effective, na.rm=TRUE), s=0.5)
  res <- tryCatch(coloc.abf(d1, d2), error=function(e) NULL)
  if (is.null(res)) return(NULL)
  pp <- res$summary
  list(gene=gene_sym, n_snps=nrow(e), lead_snp=lead$snp_id, lead_pval=lead$pvalue,
       chr=lead$chr, pos=lead$pos,
       PP.H3=pp["PP.H3.abf"], PP.H4=pp["PP.H4.abf"])
}

# ============================================================
# MAIN EXECUTION
# ============================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("Track G: MAM Pathway -> ALS Risk\n")
cat("Oligodendrocyte eQTL (Bryois 2022) x van Rheenen 2021 ALS GWAS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

results <- list()

results[["G1_PhiTop5"]] <- run_mam_track(
  mam_phi_top5, "G1_PhiTop5", eqtl_full, als)

results[["G2_Tethering"]] <- run_mam_track(
  mam_tethering, "G2_Tethering", eqtl_full, als)

results[["G3_Autophagy"]] <- run_mam_track(
  mam_autophagy, "G3_Autophagy", eqtl_full, als)

results[["G4_Lipid"]] <- run_mam_track(
  mam_lipid, "G4_Lipid", eqtl_full, als)

results[["G5_Ca2"]] <- run_mam_track(
  mam_ca2, "G5_Ca2", eqtl_full, als)

results[["G_All"]] <- run_mam_track(
  mam_all, "G_All", eqtl_full, als)

# ============================================================
# Per-gene MR
# ============================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("Per-Gene MR (all MAM genes)\n")
cat(paste(rep("=", 70), collapse=""), "\n")

per_gene <- lapply(mam_all, function(g) {
  cat(sprintf("  %-12s ... ", g))
  r <- run_per_gene(g, eqtl_full, als)
  cat(sprintf("n_iv=%d  OR=%.3f  p=%.4f\n",
              r$n_iv, ifelse(is.na(r$OR), NA, r$OR),
              ifelse(is.na(r$pval), NA, r$pval)))
  r
})
per_gene_df <- bind_rows(per_gene)
per_gene_df$set <- with(per_gene_df, {
  s <- ifelse(gene %in% mam_phi_top5, "G1_PhiTop5", "")
  s <- ifelse(gene %in% mam_tethering, paste(s, "G2"), s)
  s <- ifelse(gene %in% mam_autophagy, paste(s, "G3"), s)
  s <- ifelse(gene %in% mam_lipid,     paste(s, "G4"), s)
  s <- ifelse(gene %in% mam_ca2,       paste(s, "G5"), s)
  s <- ifelse(gene %in% mam_exome_convergent, paste(s, "Exome"), s)
  trimws(s)
})

cat("\n=== Per-Gene Results ===\n")
print(per_gene_df[, c("gene","set","n_iv","OR","OR_lo","OR_hi","pval","method")],
      row.names=FALSE)
write.csv(per_gene_df, file.path(OUTPUT_DIR, "mam_per_gene_mr.csv"), row.names=FALSE)

# ============================================================
# Colocalization
# ============================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("Colocalization (MAM genes)\n")
cat(paste(rep("=", 70), collapse=""), "\n")

coloc_list <- list()
for (g in mam_all) {
  cat(sprintf("  coloc %-12s ... ", g))
  r <- run_coloc(g, eqtl_full, als)
  if (!is.null(r)) {
    coloc_list[[g]] <- r
    cat(sprintf("PP.H4=%.3f\n", r$PP.H4))
  } else {
    cat("insufficient data\n")
  }
}

if (length(coloc_list) > 0) {
  coloc_df <- do.call(rbind, lapply(coloc_list, function(x)
    data.frame(gene=x$gene, n_snps=x$n_snps, lead_snp=x$lead_snp,
               lead_pval=x$lead_pval, chr=x$chr,
               PP.H3=round(x$PP.H3, 4), PP.H4=round(x$PP.H4, 4),
               stringsAsFactors=FALSE)))
  coloc_df <- coloc_df[order(-coloc_df$PP.H4), ]
  cat("\n=== Top Colocalization Results ===\n")
  print(coloc_df, row.names=FALSE)
  write.csv(coloc_df, file.path(OUTPUT_DIR, "mam_coloc.csv"), row.names=FALSE)
} else {
  cat("No genes with sufficient data for colocalization.\n")
}

# ============================================================
# Summary comparison across all tracks
# ============================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("SUMMARY: Track G sub-tracks\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

sum_rows <- list()
for (nm in names(results)) {
  r <- results[[nm]]
  if (is.null(r)) {
    sum_rows[[nm]] <- data.frame(Track=nm, N_IV=0, OR=NA, OR_lo=NA, OR_hi=NA,
                                  p=NA, Egger_p=NA, Q_p=NA, stringsAsFactors=FALSE)
    next
  }
  ivw <- r$mr_results[grepl("Inverse variance|Wald ratio", r$mr_results$method), ]
  if (nrow(ivw) == 0) next
  ivw <- ivw[1, ]
  pleio_p <- tryCatch(r$pleiotropy$pval[1], error=function(e) NA)
  het <- r$heterogeneity[grepl("Inverse variance", r$heterogeneity$method), ]
  Q_p <- if (nrow(het) > 0) het$Q_pval[1] else NA
  sum_rows[[nm]] <- data.frame(
    Track = nm, N_IV = ivw$nsnp,
    OR    = round(exp(ivw$b), 3),
    OR_lo = round(exp(ivw$b - 1.96 * ivw$se), 3),
    OR_hi = round(exp(ivw$b + 1.96 * ivw$se), 3),
    p     = round(ivw$pval, 4),
    Egger_p = round(pleio_p, 4),
    Q_p   = round(Q_p, 4),
    stringsAsFactors=FALSE)
}

if (length(sum_rows) > 0) {
  sum_df <- do.call(rbind, sum_rows)
  cat("\n")
  print(sum_df, row.names=FALSE)
  write.csv(sum_df, file.path(OUTPUT_DIR, "mam_summary.csv"), row.names=FALSE)
}

# Comparison with existing tracks
cat("\n--- Comparison with existing tracks (from MR_VALIDATION_REPORT.md) ---\n")
cat("Track A | All Oligo (10,462g)   | N=576 | IVW OR=1.000 (0.998-1.002) | p=0.766\n")
cat("Track B | Stable-High (135 ribo)| N=15  | IVW OR=0.988 (0.972-1.005) | p=0.159\n")
cat("Track D | Darkgrey (327 Rho/Ras)| N=74  | IVW OR=1.001 (0.995-1.008) | p=0.660\n")
cat("Track F | RQC (HBS1L)           | N=1   | Wald OR=0.957 (0.913-1.004)| p=0.072\n")

cat(sprintf("\nTrack G results saved to: %s\n", OUTPUT_DIR))
cat("Track G analysis complete.\n")
