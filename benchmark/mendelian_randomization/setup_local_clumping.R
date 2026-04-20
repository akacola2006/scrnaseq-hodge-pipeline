#!/usr/bin/env Rscript
# ============================================================
# setup_local_clumping.R
# Download PLINK and 1000G EUR reference for local LD clumping
# ============================================================

user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))

library(ieugwasr)

PLINK_DIR <- "D:/Projects/MR/tools"
REF_DIR   <- "D:/Projects/MR/data/ld_reference"
dir.create(PLINK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Check if genetics.binaRies package has PLINK
# ============================================================
cat("Checking for PLINK...\n")
plink_path <- tryCatch({
  genetics.binaRies::get_plink_binary()
}, error = function(e) {
  cat("  genetics.binaRies not available\n")
  NULL
})

if (!is.null(plink_path) && file.exists(plink_path)) {
  cat(sprintf("  PLINK found: %s\n", plink_path))
} else {
  cat("  PLINK not found via R package. Will download manually.\n")
  # Download PLINK 1.9 for Windows
  plink_url <- "https://s3.amazonaws.com/plink1-assets/plink_win64_20231018.zip"
  plink_zip <- file.path(PLINK_DIR, "plink.zip")

  cat(sprintf("  Downloading PLINK from %s...\n", plink_url))
  download.file(plink_url, plink_zip, mode = "wb", quiet = TRUE)
  unzip(plink_zip, exdir = PLINK_DIR)
  plink_path <- file.path(PLINK_DIR, "plink.exe")

  if (file.exists(plink_path)) {
    cat(sprintf("  PLINK installed: %s\n", plink_path))
  } else {
    # Try alternate name
    plink_files <- list.files(PLINK_DIR, pattern = "plink", full.names = TRUE)
    cat(sprintf("  Files in tools dir: %s\n", paste(plink_files, collapse = ", ")))
  }
}

# ============================================================
# 2. Download 1000 Genomes EUR reference panel
# ============================================================
cat("\nChecking for 1000G EUR reference panel...\n")

# MRCIEU hosts pre-formatted files for ieugwasr
ref_url <- "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
ref_file <- file.path(REF_DIR, "1kg.v3.tgz")

ref_bim <- file.path(REF_DIR, "EUR.bim")
if (file.exists(ref_bim)) {
  cat("  Reference panel already exists.\n")
} else {
  cat(sprintf("  Downloading from %s (~1.5GB, this takes a while)...\n", ref_url))
  download.file(ref_url, ref_file, mode = "wb", quiet = FALSE)

  cat("  Extracting...\n")
  untar(ref_file, exdir = REF_DIR)

  # Check what was extracted
  ref_files <- list.files(REF_DIR, pattern = "EUR|1kg", full.names = TRUE)
  cat(sprintf("  Reference files: %s\n", paste(basename(ref_files), collapse = ", ")))
}

# ============================================================
# 3. Test local clumping
# ============================================================
cat("\n=== Configuration Summary ===\n")
cat(sprintf("PLINK path: %s\n", plink_path))
cat(sprintf("Reference dir: %s\n", REF_DIR))
cat(sprintf("Reference files:\n"))
ref_all <- list.files(REF_DIR, recursive = FALSE)
cat(paste(" ", ref_all, collapse = "\n"), "\n")

cat("\nSetup complete. Use these paths in 02_run_mr.R for local clumping.\n")
