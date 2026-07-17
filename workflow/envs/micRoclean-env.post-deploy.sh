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

: "${MICROCLEAN_GIT_URL:?Set MICROCLEAN_GIT_URL, MICROCLEAN_SOURCE_ENV, or create workflow/envs/micRoclean-source.env.}"
: "${MICROCLEAN_GIT_REF:?Set MICROCLEAN_GIT_REF to the fixed patched micRoclean commit SHA.}"

if [[ ! "${MICROCLEAN_GIT_REF}" =~ ^[0-9a-fA-F]{7,40}$ ]]; then
    echo "MICROCLEAN_GIT_REF must be a fixed Git commit SHA, not a branch or moving tag: ${MICROCLEAN_GIT_REF}" >&2
    exit 2
fi

export R_REMOTES_NO_ERRORS_FROM_WARNINGS=false

Rscript --no-environ -e '
repo <- Sys.getenv("MICROCLEAN_GIT_URL")
ref <- Sys.getenv("MICROCLEAN_GIT_REF")
subdir <- Sys.getenv("MICROCLEAN_GIT_SUBDIR", unset = NA)

Sys.setenv(
  TORCH_LOAD = "0",
  TORCH_VERIFY_LOAD = "FALSE",
  R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true"
)

args <- list(
  url = repo,
  ref = ref,
  upgrade = "never",
  dependencies = FALSE,
  INSTALL_opts = "--no-test-load"
)
if (!is.na(subdir) && nzchar(subdir)) {
  args$subdir <- subdir
}

do.call(remotes::install_git, args)
'

Rscript --no-environ -e 'stopifnot(length(find.package("micRoclean", quiet = TRUE)) == 1)'

