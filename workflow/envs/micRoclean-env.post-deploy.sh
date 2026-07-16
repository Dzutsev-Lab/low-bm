#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_file="${script_dir}/micRoclean-source.env"

if [[ -f "${source_file}" ]]; then
    set -a
    # shellcheck disable=SC1090
    source "${source_file}"
    set +a
fi

: "${MICROCLEAN_GIT_URL:?Set MICROCLEAN_GIT_URL or create workflow/envs/micRoclean-source.env from micRoclean-source.env.example.}"
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

args <- list(
  url = repo,
  ref = ref,
  upgrade = "never",
  dependencies = FALSE
)
if (!is.na(subdir) && nzchar(subdir)) {
  args$subdir <- subdir
}

do.call(remotes::install_git, args)
'

Rscript --no-environ -e 'stopifnot(requireNamespace("micRoclean", quietly = TRUE))'

