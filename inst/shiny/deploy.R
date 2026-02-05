# deploy.R — Create a self-contained deployment bundle and deploy to shinyapps.io
#
# Usage:
#   source(system.file("shiny", "deploy.R", package = "evalinger"))
#
# Or from the package source directory:
#   source("inst/shiny/deploy.R")
#
# This script:
#   1. Creates a temporary deployment directory
#   2. Copies app.R and all evalinger R source files into it
#   3. Deploys to shinyapps.io via rsconnect
#
# Prerequisites:
#   install.packages(c("rsconnect", "shiny", "bslib"))
#   rsconnect::setAccountInfo(name = "...", token = "...", secret = "...")

library(rsconnect)

# --- Locate source files ---
# Try installed package first, fall back to source tree
pkg_installed <- requireNamespace("evalinger", quietly = TRUE)

if (pkg_installed) {
  app_dir <- system.file("shiny", package = "evalinger")
  # R source files are in the installed package's R/ directory
  # We need the raw .R files, not the installed lazy-load database
  # So we look for the source tree
  pkg_path <- find.package("evalinger")
  r_source <- file.path(pkg_path, "R")
  if (!dir.exists(r_source) || length(list.files(r_source, pattern = "\\.R$")) == 0) {
    # Installed via install_github — R/ contains lazy-load, not source .R files
    # Fall back to finding source tree
    r_source <- NULL
  }
} else {
  app_dir <- NULL
  r_source <- NULL
}

# Try source tree (for running from cloned repo)
if (is.null(r_source)) {
  # Look relative to this script's location
  candidates <- c(
    file.path(getwd(), "R"),                           # if run from package root
    file.path(getwd(), "..", "..", "R"),                # if run from inst/shiny/
    normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..", "R"),
                  mustWork = FALSE)
  )
  r_source <- Find(function(d) dir.exists(d) && length(list.files(d, pattern = "\\.R$")) > 0,
                    candidates)
}

if (is.null(app_dir)) {
  candidates <- c(
    file.path(getwd(), "inst", "shiny"),
    normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "."),
                  mustWork = FALSE)
  )
  app_dir <- Find(function(d) file.exists(file.path(d, "app.R")), candidates)
}

if (is.null(r_source)) stop("Cannot find evalinger R/ source files. Run from the package root or install the package from source.")
if (is.null(app_dir))  stop("Cannot find app.R. Run from the package root or install the package.")

cat("Using R source files from:", r_source, "\n")
cat("Using app.R from:", app_dir, "\n")

# --- Create deployment bundle ---
deploy_dir <- file.path(tempdir(), "evalinger-deploy")
if (dir.exists(deploy_dir)) unlink(deploy_dir, recursive = TRUE)
dir.create(deploy_dir)
dir.create(file.path(deploy_dir, "R"))

# Copy app.R
file.copy(file.path(app_dir, "app.R"), file.path(deploy_dir, "app.R"))

# Copy all R source files
r_files <- list.files(r_source, pattern = "\\.R$", full.names = TRUE)
file.copy(r_files, file.path(deploy_dir, "R"))

cat("Deployment bundle created at:", deploy_dir, "\n")
cat("Contents:\n")
cat("  app.R\n")
cat(paste0("  R/", basename(r_files), "\n"))

# --- Deploy ---
cat("\nDeploying to shinyapps.io...\n")
rsconnect::deployApp(
  appDir   = deploy_dir,
  appName  = "evalinger",
  appTitle = "evalinger: E-Values for Clinical Trials",
  forceUpdate = TRUE
)
