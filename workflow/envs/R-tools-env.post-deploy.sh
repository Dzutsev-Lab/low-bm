#!/usr/bin/env bash
set -euo pipefail

ANCOMBC_GIT_URL="${ANCOMBC_GIT_URL:-https://github.com/FrederickHuangLin/ANCOMBC.git}"
ANCOMBC_GIT_REF="${ANCOMBC_GIT_REF:-4595750750e354dfa61645f4a3f1f6c53645f683}"
CODA4MICROBIOME_VERSION="${CODA4MICROBIOME_VERSION:-0.2.4}"
CODA4MICROBIOME_REPOS="${CODA4MICROBIOME_REPOS:-https://cloud.r-project.org}"

if [[ ! "${ANCOMBC_GIT_REF}" =~ ^[0-9a-fA-F]{7,40}$ ]]; then
    echo "ANCOMBC_GIT_REF must be a fixed Git commit SHA, not a branch or moving tag: ${ANCOMBC_GIT_REF}" >&2
    exit 2
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

ancombc_dir="${tmpdir}/ANCOMBC"
echo "Cloning ANCOMBC..."
git clone "${ANCOMBC_GIT_URL}" "${ancombc_dir}"
git -C "${ancombc_dir}" checkout --detach "${ANCOMBC_GIT_REF}"

echo "Installing upstream ANCOMBC with quadprog trend optimization..."
R CMD INSTALL --no-test-load "${ancombc_dir}"

Rscript --no-environ -e 'stopifnot(length(find.package("ANCOMBC", quiet = TRUE)) == 1); stopifnot(requireNamespace("microbiome", quietly = TRUE)); stopifnot(requireNamespace("quadprog", quietly = TRUE)); stopifnot(utils::packageVersion("ANCOMBC") >= "2.13.2"); stopifnot(!"CVXR" %in% names(getNamespaceImports("ANCOMBC")))'

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
