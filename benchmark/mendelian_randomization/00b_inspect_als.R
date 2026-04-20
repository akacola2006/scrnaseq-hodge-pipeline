user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))
library(data.table)

als <- fread("D:/Projects/MR/data/als_gwas/GCST90027164_buildGRCh37.tsv.gz", nrows = 5)
cat("Columns:", paste(names(als), collapse = ", "), "\n")
cat("Ncol:", ncol(als), "\n\n")
print(als)
