# R dependency installer for python_scmega
# Target runtime: R 4.5.x (e.g. 4.5.3 on Ubuntu 24.04 + CRAN noble-cran40)

userlib <- path.expand("~/R/library")
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))

cat("R version:", R.version.string, "\n")
cat("Using .libPaths():\n")
print(.libPaths())

# Preflight checks for system-level dependencies.
has_openblas <- nzchar(system("ldconfig -p 2>/dev/null | grep -F 'libopenblas.so.0' | head -n 1", intern = TRUE))
has_gsl <- nzchar(Sys.which("gsl-config"))
has_cmake <- nzchar(Sys.which("cmake"))
if (!has_openblas || !has_gsl || !has_cmake) {
  cat("\nMissing system dependencies detected.\n")
  if (!has_openblas) cat("- Missing: libopenblas.so.0 (apt: libopenblas0 libopenblas-dev)\n")
  if (!has_gsl) cat("- Missing: gsl-config (apt: libgsl-dev)\n")
  if (!has_cmake) cat("- Missing: cmake (apt: cmake)\n")
  cat("\nPlease run in WSL:\n")
  cat("  sudo apt update && sudo apt install -y libopenblas0 libopenblas-dev libgsl-dev cmake libuv1-dev\n\n")
  stop("System dependency check failed. Install OS packages first, then rerun this script.")
}

# Bioconductor version mapping (keep aligned with R major.minor).
r_ver <- getRversion()
if (r_ver >= "4.5.0" && r_ver < "4.6.0") {
  target_bioc_version <- "3.22"
} else {
  stop(
    sprintf(
      "Unsupported R version for this script: %s. Please update the Bioconductor mapping first.",
      as.character(r_ver)
    )
  )
}
cat("Target Bioconductor version:", target_bioc_version, "\n")

# Cleanup stale lock dirs from interrupted installs.
locks <- Sys.glob(file.path(userlib, "00LOCK*"))
if (length(locks) > 0) {
  unlink(locks, recursive = TRUE, force = TRUE)
}

# Ensure BiocManager is available from CRAN.
install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = userlib)

# Pin Bioconductor repository set for this R runtime.
BiocManager::install(version = target_bioc_version, ask = FALSE, update = FALSE, lib = userlib)
options(repos = BiocManager::repositories(version = target_bioc_version))

# Stage 1: root dependencies that commonly trigger cascades.
install.packages(c("fs", "igraph"), repos = "https://cloud.r-project.org", lib = userlib)
BiocManager::install("TFMPvalue", ask = FALSE, update = FALSE, lib = userlib)

# Stage 2: CRAN visualization chain.
install.packages(
  c("sass", "bslib", "rmarkdown", "htmlwidgets", "shiny", "miniUI", "plotly", "DT"),
  repos = "https://cloud.r-project.org",
  lib = userlib
)

# Stage 3: Bioconductor chain + targets.
bioc_pkgs <- c(
  "JASPAR2020",
  "bluster",
  "DirichletMultinomial",
  "scran",
  "TFBSTools",
  "motifmatchr",
  "chromVAR",
  "scDblFinder"
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE, lib = userlib)

# Final verification.
check_pkgs <- c("TFBSTools", "JASPAR2020", "motifmatchr", "chromVAR", "scDblFinder")
status <- sapply(check_pkgs, function(p) requireNamespace(p, quietly = TRUE))

cat("\nInstall check:\n")
for (p in check_pkgs) {
  cat(p, "=>", status[[p]], "\n")
}

if (!all(status)) {
  failed <- names(status)[!status]
  stop(sprintf("Installation incomplete. Missing packages: %s", paste(failed, collapse = ", ")))
}

cat("\nAll required R packages are installed.\n")
