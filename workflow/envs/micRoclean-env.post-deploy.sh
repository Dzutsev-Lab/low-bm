#!/usr/bin/env bash
set -euo pipefail

find_source_file() {
    local candidate

    # 1. Explicit override, useful for HPC/containers/custom launchers.
    if [[ -n "${MICROCLEAN_SOURCE_ENV:-}" ]]; then
        if [[ -f "${MICROCLEAN_SOURCE_ENV}" ]]; then
            printf '%s\n' "${MICROCLEAN_SOURCE_ENV}"
            return 0
        fi
        echo "MICROCLEAN_SOURCE_ENV is set but does not exist: ${MICROCLEAN_SOURCE_ENV}" >&2
        return 1
    fi

    # 2. Normal Snakemake execution: run from the workflow checkout.
    candidate="${PWD}/workflow/envs/micRoclean-source.env"
    if [[ -f "${candidate}" ]]; then
        printf '%s\n' "${candidate}"
        return 0
    fi

    # 3. Direct/manual execution from the original workflow/envs directory.
    candidate="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/micRoclean-source.env"
    if [[ -f "${candidate}" ]]; then
        printf '%s\n' "${candidate}"
        return 0
    fi

    return 1
}

if source_file="$(find_source_file)"; then
    set -a
    # shellcheck disable=SC1090
    source "${source_file}"
    set +a
fi

export R_REMOTES_NO_ERRORS_FROM_WARNINGS=false
export TORCH_VERIFY_LOAD=FALSE
export TORCH_LOAD=0

#Installing SCRuB (unused dependency for micRoclean)
SCRUB_GIT_URL="${SCRUB_GIT_URL:-https://github.com/Shenhav-and-Korem-labs/SCRuB.git}"
: "${SCRUB_GIT_REF:?Set SCRUB_GIT_REF to the fixed SCRuB commit SHA.}"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

scrub_dir="${tmpdir}/SCRuB"
echo "Cloning SCRuB..."
git clone "${SCRUB_GIT_URL}" "${scrub_dir}"
git -C "${scrub_dir}" checkout --detach "${SCRUB_GIT_REF}"

echo "Installing SCRuB without test-loading torch..."
R CMD INSTALL --no-test-load "${scrub_dir}"

Rscript --no-environ -e 'stopifnot(length(find.package("SCRuB", quiet = TRUE)) == 1)'

ANCOMBC_GIT_URL="${ANCOMBC_GIT_URL:-https://github.com/FrederickHuangLin/ANCOMBC.git}"
: "${ANCOMBC_GIT_REF:?Set ANCOMBC_GIT_REF to the fixed upstream ANCOMBC commit SHA.}"

if [[ ! "${ANCOMBC_GIT_REF}" =~ ^[0-9a-fA-F]{7,40}$ ]]; then
    echo "ANCOMBC_GIT_REF must be a fixed Git commit SHA, not a branch or moving tag: ${ANCOMBC_GIT_REF}" >&2
    exit 2
fi

ancombc_dir="${tmpdir}/ANCOMBC"
echo "Cloning ANCOMBC..."
git clone "${ANCOMBC_GIT_URL}" "${ancombc_dir}"
git -C "${ancombc_dir}" checkout --detach "${ANCOMBC_GIT_REF}"

echo "Installing upstream ANCOMBC with quadprog trend optimization..."
R CMD INSTALL --no-test-load "${ancombc_dir}"

Rscript --no-environ -e 'stopifnot(length(find.package("ANCOMBC", quiet = TRUE)) == 1); stopifnot(requireNamespace("microbiome", quietly = TRUE)); stopifnot(requireNamespace("quadprog", quietly = TRUE)); stopifnot(utils::packageVersion("ANCOMBC") >= "2.13.2"); stopifnot(!"CVXR" %in% names(getNamespaceImports("ANCOMBC")))'

: "${MICROCLEAN_GIT_URL:?Set MICROCLEAN_GIT_URL, MICROCLEAN_SOURCE_ENV, or create workflow/envs/micRoclean-source.env.}"
: "${MICROCLEAN_GIT_REF:?Set MICROCLEAN_GIT_REF to the fixed patched micRoclean commit SHA.}"

if [[ ! "${MICROCLEAN_GIT_REF}" =~ ^[0-9a-fA-F]{7,40}$ ]]; then
    echo "MICROCLEAN_GIT_REF must be a fixed Git commit SHA, not a branch or moving tag: ${MICROCLEAN_GIT_REF}" >&2
    exit 2
fi

microclean_dir="${tmpdir}/micRoclean"
echo "Cloning micRoclean..."
git clone "${MICROCLEAN_GIT_URL}" "${microclean_dir}"
git -C "${microclean_dir}" checkout --detach "${MICROCLEAN_GIT_REF}"

package_dir="${microclean_dir}"
if [[ -n "${MICROCLEAN_GIT_SUBDIR:-}" ]]; then
    package_dir="${microclean_dir}/${MICROCLEAN_GIT_SUBDIR}"
fi
if [[ ! -f "${package_dir}/DESCRIPTION" || ! -f "${package_dir}/NAMESPACE" ]]; then
    echo "micRoclean package source is missing DESCRIPTION or NAMESPACE: ${package_dir}" >&2
    exit 2
fi

echo "Patching micRoclean namespace imports..."
Rscript --no-environ - "${package_dir}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
package_dir <- args[[1]]
required_imports <- c("stringr", "tibble")

namespace_path <- file.path(package_dir, "NAMESPACE")
namespace_lines <- readLines(namespace_path, warn = FALSE)
for (pkg in required_imports) {
  import_line <- paste0("import(", pkg, ")")
  if (!any(namespace_lines == import_line)) {
    import_idx <- grep("^import\\(", namespace_lines)
    insert_after <- if (length(import_idx) > 0) max(import_idx) else length(namespace_lines)
    namespace_lines <- append(namespace_lines, import_line, after = insert_after)
  }
}
writeLines(namespace_lines, namespace_path)

description_path <- file.path(package_dir, "DESCRIPTION")
description <- read.dcf(description_path)
imports <- if ("Imports" %in% colnames(description)) description[1, "Imports"] else ""
has_import <- function(pkg) {
  grepl(paste0("(^|[,\\n[:space:]])", pkg, "([,[:space:]\\n]|$)"), imports)
}
missing_imports <- required_imports[!vapply(required_imports, has_import, logical(1))]
if (length(missing_imports) > 0) {
  import_suffix <- paste(missing_imports, collapse = ",\n    ")
  new_imports <- if (nzchar(trimws(imports))) paste0(imports, ",\n    ", import_suffix) else import_suffix
  if ("Imports" %in% colnames(description)) {
    description[1, "Imports"] <- new_imports
  } else {
    description <- cbind(description, Imports = new_imports)
  }
  write.dcf(description, description_path)
}
RSCRIPT

echo "Installing patched micRoclean without test-loading torch..."
R CMD INSTALL --no-test-load "${package_dir}"

Rscript --no-environ -e 'ns <- asNamespace("micRoclean"); stopifnot(length(find.package("micRoclean", quiet = TRUE)) == 1); stopifnot(is.function(get("str_extract", envir = ns, inherits = TRUE))); stopifnot(is.function(get("add_column", envir = ns, inherits = TRUE)))'
