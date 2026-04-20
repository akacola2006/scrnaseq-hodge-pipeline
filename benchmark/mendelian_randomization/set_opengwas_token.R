user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
.libPaths(c(user_lib, .libPaths()))

# Set OpenGWAS JWT token
token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJha2Fjb2xhMjAwNkB5YWhvby5jby5qcCIsImlhdCI6MTc3NDMxNTE0MywiZXhwIjoxNzc1NTI0NzQzfQ.IrgBuXh8Els1P6zl-Zc-csBG7jNvYMW6kZCuaOLKWcwmXKfbOR3aoYeJq5D74iRPLJREPXyLOrOYxEVhQAZjwusit-2SFkCStP4f_9o2u3sarBdnWNC-fAWPiAx5gV1pje3g35XZeVEY5opB1detO-vWGT_Tzlsw4OyI1vJjwD1oU8pNNSuUfz_ZfjUdY7O_rGu1LoLRTrqVWLPXFNwA7wPppX7T0bLr76NOIil23wXmjYkoXrH4vUkckj9N96fZqhaLSz87ftoyOgdfAr0DarzWht6iLDEk4pCbz-Y7vg3BULOk5n9rFQlwd-d4QHHisRO3Qdr4Vixg_uWzowK2bA"

Sys.setenv(OPENGWAS_JWT = token)

# Test the token
library(ieugwasr)
cat("Testing OpenGWAS API with JWT...\n")
result <- tryCatch({
  gwasinfo("ebi-a-GCST90027164")
}, error = function(e) {
  cat(sprintf("Error: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(result)) {
  cat("SUCCESS! Token works.\n")
  cat(sprintf("  Study: %s\n", result$trait))
  cat(sprintf("  Sample size: %d\n", result$sample_size))

  # Also write to .Renviron for persistence
  renviron <- file.path(Sys.getenv("USERPROFILE"), ".Renviron")
  lines <- if (file.exists(renviron)) readLines(renviron) else character(0)
  # Remove old token if exists
  lines <- lines[!grepl("^OPENGWAS_JWT=", lines)]
  lines <- c(lines, paste0("OPENGWAS_JWT=", token))
  writeLines(lines, renviron)
  cat(sprintf("  Token saved to %s\n", renviron))
} else {
  cat("FAILED. Token may be invalid.\n")
}
