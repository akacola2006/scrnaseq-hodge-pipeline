user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))

install.packages("R.utils", lib = user_lib, repos = "https://cran.r-project.org", quiet = TRUE)
library(data.table)

DATA_DIR <- "D:/Projects/MR/data/bryois_eqtl"

# 1. SNP position file
cat("=== SNP Position File ===\n")
snp_pos <- fread(file.path(DATA_DIR, "snp_pos.txt.gz"), nrows = 10)
cat("Columns:", paste(names(snp_pos), collapse = ", "), "\n")
cat("Ncol:", ncol(snp_pos), "\n")
print(head(snp_pos, 5))

# 2. eQTL file (chr 22 - smallest)
cat("\n=== Oligodendrocytes eQTL (chr 22, first 10 rows) ===\n")
eqtl <- fread(file.path(DATA_DIR, "Oligodendrocytes/Oligodendrocytes.22.gz"), nrows = 10)
cat("Columns:", paste(names(eqtl), collapse = ", "), "\n")
cat("Ncol:", ncol(eqtl), "\n")
print(head(eqtl, 5))

# 3. Count total rows in chr 22
cat("\n=== Row count chr 22 ===\n")
eqtl_full <- fread(file.path(DATA_DIR, "Oligodendrocytes/Oligodendrocytes.22.gz"))
cat("Total rows chr22:", nrow(eqtl_full), "\n")
cat("Unique genes:", length(unique(eqtl_full$V1)), "\n")
cat("Unique SNPs:", length(unique(eqtl_full$V2)), "\n")
cat("Min p-value:", min(eqtl_full$V4), "\n")
cat("N significant (p<5e-8):", sum(eqtl_full$V4 < 5e-8), "\n")

# 4. Check gene ID format
cat("\nSample gene IDs:", paste(head(unique(eqtl_full$V1), 5), collapse = ", "), "\n")

# 5. Check ALS GWAS
als_file <- "D:/Projects/MR/data/als_gwas/GCST90027164_buildGRCh37.tsv.gz"
if (file.exists(als_file) && file.size(als_file) > 1000) {
  cat("\n=== ALS GWAS ===\n")
  als <- fread(als_file, nrows = 5)
  cat("Columns:", paste(names(als), collapse = ", "), "\n")
  print(head(als, 3))
} else {
  cat("\n=== ALS GWAS file not yet downloaded or incomplete ===\n")
}
