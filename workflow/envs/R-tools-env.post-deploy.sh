#!/usr/bin/env bash
set -euo pipefail

CODA4MICROBIOME_VERSION="${CODA4MICROBIOME_VERSION:-0.2.4}"
CODA4MICROBIOME_REPOS="${CODA4MICROBIOME_REPOS:-https://cloud.r-project.org}"

echo "Installing coda4microbiome ${CODA4MICROBIOME_VERSION} from CRAN..."
Rscript --no-environ - "${CODA4MICROBIOME_VERSION}" "${CODA4MICROBIOME_REPOS}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
coda_version <- args[[1]]
repos <- args[[2]]

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("The R package 'remotes' is required to install coda4microbiome.", call. = FALSE)
}

installed <- requireNamespace("coda4microbiome", quietly = TRUE)
if (!installed || as.character(utils::packageVersion("coda4microbiome")) != coda_version) {
  remotes::install_version(
    "coda4microbiome",
    version = coda_version,
    repos = repos,
    dependencies = FALSE,
    upgrade = "never"
  )
}

stopifnot(
  requireNamespace("coda4microbiome", quietly = TRUE),
  as.character(utils::packageVersion("coda4microbiome")) == coda_version,
  is.function(getExportedValue("coda4microbiome", "coda_coxnet"))
)
RSCRIPT
