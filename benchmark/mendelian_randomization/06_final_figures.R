#!/usr/bin/env Rscript
# ============================================================
# 06_final_figures.R
# Publication-quality figures for MR + Coloc results
# ============================================================

user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))

library(data.table)
library(dplyr)
library(ggplot2)

OUTPUT_DIR <- "D:/Projects/MR/results_final"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Figure 1: Multi-track MR Forest Plot (publication quality)
# ============================================================
cat("=== Figure 1: MR Forest Plot ===\n")

# Combine p<5e-8 and p<5e-6 results
mr_data <- data.frame(
  track = c(
    "A: Oligo All (p<5e-8)", "A: Oligo All (p<5e-6)",
    "B: Stable-High 135 (p<5e-8)", "B: Stable-High 135 (p<5e-6)",
    "D: Darkgrey 327 (p<5e-8)", "D: Darkgrey 327 (p<5e-6)",
    "E: RTK/Trophic (p<5e-8)", "E: RTK/Trophic (p<5e-6)"
  ),
  n_iv = c(431, 576, 9, 15, 53, 74, 6, 8),
  beta = c(-0.0004, -0.0003, -0.0096, -0.0118, 0.0012, 0.0015, -0.0052, -0.0031),
  se =   c(0.0011, 0.0011, 0.0093, 0.0084, 0.0038, 0.0034, 0.0134, 0.0118),
  pval = c(0.741, 0.766, 0.304, 0.159, 0.758, 0.660, 0.699, 0.795),
  threshold = c("5e-8", "5e-6", "5e-8", "5e-6", "5e-8", "5e-6", "5e-8", "5e-6"),
  gene_set = c("All", "All", "Stable-High", "Stable-High",
               "Darkgrey", "Darkgrey", "RTK/Trophic", "RTK/Trophic"),
  stringsAsFactors = FALSE
)

mr_data$OR <- exp(mr_data$beta)
mr_data$CI_lo <- exp(mr_data$beta - 1.96 * mr_data$se)
mr_data$CI_hi <- exp(mr_data$beta + 1.96 * mr_data$se)
mr_data$label <- sprintf("%d IVs", mr_data$n_iv)
mr_data$pval_label <- ifelse(mr_data$pval < 0.001,
                              sprintf("p=%.0e", mr_data$pval),
                              sprintf("p=%.3f", mr_data$pval))

# Order for display
mr_data$track <- factor(mr_data$track,
  levels = rev(c(
    "A: Oligo All (p<5e-8)", "A: Oligo All (p<5e-6)",
    "B: Stable-High 135 (p<5e-8)", "B: Stable-High 135 (p<5e-6)",
    "D: Darkgrey 327 (p<5e-8)", "D: Darkgrey 327 (p<5e-6)",
    "E: RTK/Trophic (p<5e-8)", "E: RTK/Trophic (p<5e-6)"
  ))
)

fig1 <- ggplot(mr_data, aes(x = track, y = OR, color = threshold, shape = gene_set)) +
  geom_point(size = 3, position = position_dodge(width = 0)) +
  geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = label), hjust = -0.3, size = 2.8, color = "grey30",
            position = position_dodge(width = 0)) +
  coord_flip(ylim = c(0.94, 1.06)) +
  scale_color_manual(values = c("5e-8" = "#2166AC", "5e-6" = "#B2182B"),
                     name = "eQTL threshold") +
  scale_shape_manual(values = c("All" = 16, "Stable-High" = 17,
                                 "Darkgrey" = 15, "RTK/Trophic" = 18),
                     name = "Gene set") +
  labs(
    title = "Mendelian Randomization: Oligodendrocyte eQTL \u2192 ALS Risk",
    subtitle = "IVW estimates with 95% CI | Bryois et al. 2022 \u00d7 van Rheenen et al. 2021",
    x = "", y = "Odds Ratio (95% CI)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

ggsave(file.path(OUTPUT_DIR, "Fig_MR_forest.png"), fig1, width = 11, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Fig_MR_forest.pdf"), fig1, width = 11, height = 6)
cat("  Saved Fig_MR_forest\n")

# ============================================================
# Figure 2: Track B Detail (4-method comparison)
# ============================================================
cat("=== Figure 2: Track B Detail ===\n")

trackB_methods <- data.frame(
  method = c("IVW", "MR-Egger", "Weighted Median", "Weighted Mode"),
  threshold = rep(c("p<5e-8", "p<5e-6"), each = 4),
  beta = c(-0.0096, -0.0064, -0.0159, -0.0260,   # p<5e-8
           -0.0118, -0.0066, -0.0157, -0.0166),   # p<5e-6
  se = c(0.0093, 0.0271, 0.0125, 0.0239,
         0.0084, 0.0211, 0.0120, 0.0337),
  stringsAsFactors = FALSE
)
trackB_methods$OR <- exp(trackB_methods$beta)
trackB_methods$CI_lo <- exp(trackB_methods$beta - 1.96 * trackB_methods$se)
trackB_methods$CI_hi <- exp(trackB_methods$beta + 1.96 * trackB_methods$se)
trackB_methods$method <- factor(trackB_methods$method,
  levels = c("IVW", "MR-Egger", "Weighted Median", "Weighted Mode"))

fig2 <- ggplot(trackB_methods, aes(x = method, y = OR, color = threshold)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi), width = 0.3,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("p<5e-8" = "#2166AC", "p<5e-6" = "#B2182B")) +
  coord_flip(ylim = c(0.92, 1.08)) +
  labs(
    title = "Track B: Stable-High Genes (Ribosomal) \u2192 ALS Risk",
    subtitle = "All methods show consistent protective direction (OR < 1)\nEgger intercept p = 0.79 (no pleiotropy) | Q p = 0.49 (no heterogeneity)",
    x = "", y = "Odds Ratio (95% CI)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 13, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "Fig_TrackB_detail.png"), fig2, width = 9, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Fig_TrackB_detail.pdf"), fig2, width = 9, height = 5)
cat("  Saved Fig_TrackB_detail\n")

# ============================================================
# Figure 3: Colocalization Heatmap
# ============================================================
cat("=== Figure 3: Coloc Results ===\n")

coloc_df <- fread("D:/Projects/MR/results_coloc/coloc_results.csv")

# Top 15 genes by PP.H4
top_coloc <- head(coloc_df[order(-PP.H4)], 15)

# Reshape for heatmap
coloc_long <- melt(top_coloc,
  id.vars = c("gene", "n_snps", "lead_snp", "lead_pval", "chr", "pos", "in_stable_high"),
  measure.vars = c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"),
  variable.name = "hypothesis", value.name = "PP"
)

coloc_long$hypothesis <- factor(coloc_long$hypothesis,
  levels = c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"),
  labels = c("H0: Neither", "H1: eQTL only", "H2: GWAS only",
             "H3: Both, distinct", "H4: Shared causal")
)

coloc_long$gene_label <- sprintf("%s (chr%d)", coloc_long$gene, coloc_long$chr)
coloc_long$gene_label <- factor(coloc_long$gene_label,
  levels = rev(unique(coloc_long$gene_label[order(coloc_long$PP[coloc_long$hypothesis == "H4: Shared causal"])])))

fig3 <- ggplot(coloc_long, aes(x = hypothesis, y = gene_label, fill = PP)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", PP)), size = 2.8) +
  scale_fill_gradient2(low = "white", mid = "#FEE08B", high = "#D73027",
                       midpoint = 0.4, limits = c(0, 1),
                       name = "Posterior\nProbability") +
  labs(
    title = "Colocalization: Oligodendrocyte eQTL \u00d7 ALS GWAS",
    subtitle = "Top 15 genes by PP.H4 | KANSL1 (PP.H4=0.79) is the strongest candidate",
    x = "", y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(size = 13, face = "bold"),
    panel.grid = element_blank()
  )

ggsave(file.path(OUTPUT_DIR, "Fig_coloc_heatmap.png"), fig3, width = 10, height = 7, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Fig_coloc_heatmap.pdf"), fig3, width = 10, height = 7)
cat("  Saved Fig_coloc_heatmap\n")

# ============================================================
# Figure 4: Power Analysis
# ============================================================
cat("=== Figure 4: Power Curve ===\n")

# Power as function of outcome sample size
or_range <- seq(1.001, 1.05, by = 0.001)
n_range <- c(50000, 80713, 150000, 300000, 500000)

power_grid <- expand.grid(OR = or_range, N = n_range)
power_grid$beta <- log(power_grid$OR)

# For Track B: 15 IVs, R2 ~ 0.72
R2_trackB <- 0.7248
power_grid$ncp <- power_grid$N * R2_trackB * power_grid$beta^2
power_grid$power <- pnorm(sqrt(power_grid$ncp) - qnorm(0.975))
power_grid$N_label <- sprintf("N=%.0fk", power_grid$N / 1000)

fig4 <- ggplot(power_grid, aes(x = OR, y = power, color = N_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 1.012, linetype = "dotted", color = "red", alpha = 0.5) +
  annotate("text", x = 1.014, y = 0.05, label = "OR=1.012\n(current limit)",
           size = 3, color = "red") +
  annotate("text", x = 1.04, y = 0.83, label = "80% power", size = 3, color = "grey40") +
  scale_color_manual(
    values = c("N=50k" = "#4DAF4A", "N=81k" = "#377EB8",
               "N=150k" = "#FF7F00", "N=300k" = "#E41A1C", "N=500k" = "#984EA3"),
    name = "ALS GWAS\nSample Size"
  ) +
  labs(
    title = "Statistical Power: Stable-High eQTL \u2192 ALS (Track B)",
    subtitle = "15 IVs, R\u00b2=0.72 | Current N_eff=80,713 can detect OR\u22651.012",
    x = "Odds Ratio (effect size)", y = "Power"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    legend.position = c(0.85, 0.3)
  )

ggsave(file.path(OUTPUT_DIR, "Fig_power.png"), fig4, width = 9, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Fig_power.pdf"), fig4, width = 9, height = 6)
cat("  Saved Fig_power\n")

# ============================================================
# Summary table for paper
# ============================================================
cat("\n=== Summary Table ===\n")

summary_table <- data.frame(
  Track = c("A: All Oligo eQTL", "B: Stable-High (ribosomal)",
            "D: Darkgrey (Rho/Ras GTPase)", "E: RTK/Trophic pathway"),
  Gene_set_size = c("10,462 genes", "135 genes", "327 genes", "20 genes"),
  N_IV_5e8 = c(431, 9, 53, 6),
  N_IV_5e6 = c(576, 15, 74, 8),
  IVW_OR_5e6 = c("1.000", "0.988", "1.001", "0.997"),
  IVW_CI_5e6 = c("0.998-1.002", "0.972-1.005", "0.995-1.008", "0.974-1.020"),
  IVW_p_5e6 = c(0.766, 0.159, 0.660, 0.795),
  Egger_p = c(0.122, 0.791, 0.239, 0.356),
  Q_p = c("<0.001", "0.493", "0.002", "0.040"),
  Direction_consistent = c("Yes", "Yes (all 4 methods)", "No", "Yes"),
  stringsAsFactors = FALSE
)

cat("\n")
print(summary_table, row.names = FALSE)

write.csv(summary_table, file.path(OUTPUT_DIR, "Table_MR_summary.csv"), row.names = FALSE)

# KANSL1 coloc summary
cat("\n=== KANSL1 Colocalization ===\n")
kansl1 <- coloc_df[coloc_df$gene == "KANSL1", ]
cat(sprintf("Gene: KANSL1 (chr%d:%d)\n", kansl1$chr, kansl1$pos))
cat(sprintf("Lead SNP: %s (eQTL p = %.2e)\n", kansl1$lead_snp, kansl1$lead_pval))
cat(sprintf("N SNPs in window: %d\n", kansl1$n_snps))
cat(sprintf("PP.H4 (shared causal): %.3f\n", kansl1$PP.H4))
cat(sprintf("PP.H3 (distinct causal): %.3f\n", kansl1$PP.H3))
cat(sprintf("In stable-High: %s\n", kansl1$in_stable_high))

cat("\nAll figures saved to:", OUTPUT_DIR, "\n")
