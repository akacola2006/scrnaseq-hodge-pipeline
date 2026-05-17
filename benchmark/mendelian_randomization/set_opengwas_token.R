# Configure OpenGWAS JWT for ieugwasr.
#
# The token is read from the OPENGWAS_JWT environment variable (or .Renviron).
# Obtain a token from https://api.opengwas.io/profile/ and either:
#   (1) export OPENGWAS_JWT="<your token>" before launching R, or
#   (2) add OPENGWAS_JWT=<your token> to your ~/.Renviron.
#
# Never commit a JWT to a public repository: any committed token must be
# revoked at https://api.opengwas.io/profile/ before re-issuing a new one.

user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", "4.4")
if (nzchar(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

token <- Sys.getenv("OPENGWAS_JWT")
if (!nzchar(token)) {
  stop(
    "OPENGWAS_JWT is not set. Obtain a JWT from https://api.opengwas.io/profile/ ",
    "and export it as the OPENGWAS_JWT environment variable (or add it to .Renviron)."
  )
}

# ieugwasr reads OPENGWAS_JWT from the environment; the explicit Sys.setenv
# call below is kept for clarity and to support shells that did not export it.
Sys.setenv(OPENGWAS_JWT = token)

library(ieugwasr)
cat("Testing OpenGWAS API with JWT from OPENGWAS_JWT...\n")
result <- tryCatch(
  gwasinfo("ebi-a-GCST90027164"),
  error = function(e) {
    cat(sprintf("Error: %s\n", conditionMessage(e)))
    NULL
  }
)

if (!is.null(result)) {
  cat("SUCCESS! Token works.\n")
  cat(sprintf("  Study: %s\n", result$trait))
  cat(sprintf("  Sample size: %d\n", result$sample_size))
} else {
  cat("FAILED. Token may be invalid or expired -- check https://api.opengwas.io/profile/.\n")
}
