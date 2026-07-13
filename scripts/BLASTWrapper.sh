#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/BLASTWrapper.sh --analysis-config analysis_config.yaml
  bash scripts/BLASTWrapper.sh <blast_analysis_dir> <reference_db> <comparison> [comparison ...]

Config mode reads blast_confirmation from analysis_config.yaml and writes BLAST hits
inside <io_dir>/BlastAnalysis/<comparison>/<taxon>/. In config mode,
blast_confirmation.candidate_comparisons is preferred over the legacy
blast_confirmation.DA_comparisons.
USAGE
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required command not found: $cmd" >&2
    exit 1
  fi
}

load_config() {
  local analysis_config="$1"
  local config_text

  config_text="$(
    Rscript -e '
args <- commandArgs(trailingOnly = TRUE)
source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
cfg <- load_yaml_config(args[1])
project <- cfg$project %||% list()
blast <- cfg$blast_confirmation
if (is.null(blast)) {
  stop("analysis_config.yaml is missing blast_confirmation.", call. = FALSE)
}
is_missing <- function(x) {
  is.null(x) || length(x) == 0 || all(is.na(x)) || all(!nzchar(trimws(as.character(x))))
}
io_dir <- analysis_output_dir(project, blast, section_keys = c("io_dir", "output_dir", "out_dir"))
comparisons <- config_value(blast, "candidate_comparisons") %||% config_value(blast, "DA_comparisons")
required <- c("reference_db")
missing <- required[vapply(required, function(k) is_missing(config_value(blast, k)), logical(1))]
if (is_missing(io_dir)) {
  missing <- c("project.output_dir or blast_confirmation.io_dir", missing)
}
if (is_missing(comparisons)) {
  missing <- c("blast_confirmation.candidate_comparisons or blast_confirmation.DA_comparisons", missing)
}
if (length(missing) > 0) {
  stop("blast_confirmation is missing required field(s): ", paste(missing, collapse = ", "), call. = FALSE)
}
emit <- function(key, value) {
  if (is_missing(value)) value <- ""
  cat(key, "\t", paste(as.character(value), collapse = ","), "\n", sep = "")
}
emit("TRIAL_DIR", file.path(io_dir, "BlastAnalysis"))
emit("REFERENCE_DB", config_value(blast, "reference_db"))
emit("REF_BASE_DIR", config_value(blast, "ref_base_dir") %||% "Ref_Data")
emit("NUM_THREADS", config_value(blast, "num_threads"))
emit("MAX_TARGET_SEQS", config_value(blast, "max_target_seqs") %||% 5)
for (comparison in comparisons) {
  emit("COMPARISON", comparison)
}
' "$analysis_config"
  )"

  COMPARISONS=()
  while IFS=$'\t' read -r key value; do
    case "$key" in
      TRIAL_DIR) TRIAL_DIR="$value" ;;
      REFERENCE_DB) REFERENCE_DB="$value" ;;
      REF_BASE_DIR) REF_BASE_DIR="$value" ;;
      NUM_THREADS) NUM_THREADS="$value" ;;
      MAX_TARGET_SEQS) MAX_TARGET_SEQS="$value" ;;
      COMPARISON) COMPARISONS+=("$value") ;;
    esac
  done <<< "$config_text"
}

blast_db_exists() {
  [ -f "${BLAST_DB_BASE}.nsq" ] || [ -f "${BLAST_DB_BASE}.ndb" ]
}

if [ "$#" -eq 0 ]; then
  usage
  exit 1
fi

TRIAL_DIR=""
REFERENCE_DB=""
REF_BASE_DIR="Ref_Data"
NUM_THREADS=""
MAX_TARGET_SEQS=5
COMPARISONS=()
OUTFMT='6 qseqid sseqid pident length qcovs evalue bitscore staxids'

if [ "${1:-}" = "--analysis-config" ]; then
  if [ "$#" -ne 2 ]; then
    usage
    exit 1
  fi
  load_config "$2"
else
  if [ "$#" -lt 3 ]; then
    usage
    exit 1
  fi
  TRIAL_DIR="$1"
  REFERENCE_DB="$2"
  shift 2
  COMPARISONS=("$@")
fi

if [ -z "$TRIAL_DIR" ] || [ -z "$REFERENCE_DB" ] || [ "${#COMPARISONS[@]}" -eq 0 ]; then
  echo "Missing required BLAST input. Check analysis_config.yaml or legacy arguments." >&2
  exit 1
fi

if [ -z "${NUM_THREADS:-}" ]; then
  NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
fi

REFERENCE_FASTA="${REF_BASE_DIR}/${REFERENCE_DB}.fasta"
TAXID_MAP="${REF_BASE_DIR}/${REFERENCE_DB}/seqid2taxid.map"
BLAST_DB_BASE="${REF_BASE_DIR}/${REFERENCE_DB}/blast_format_db/${REFERENCE_DB}_blast"

if [ ! -f "$REFERENCE_FASTA" ]; then
  echo "Missing reference FASTA: $REFERENCE_FASTA" >&2
  exit 1
fi

if [ ! -f "$TAXID_MAP" ]; then
  echo "Missing BLAST taxid map: $TAXID_MAP" >&2
  exit 1
fi

require_command blastn

if ! blast_db_exists; then
  require_command makeblastdb
  echo "BLAST database not found. Building database at: ${BLAST_DB_BASE}"
  mkdir -p "$(dirname "$BLAST_DB_BASE")"
  makeblastdb \
    -in "$REFERENCE_FASTA" \
    -dbtype nucl \
    -taxid_map "$TAXID_MAP" \
    -out "$BLAST_DB_BASE"
else
  echo "BLAST database already exists at: ${BLAST_DB_BASE}"
fi

for comparison in "${COMPARISONS[@]}"; do
  comp_dir="${TRIAL_DIR}/${comparison}"

  if [ ! -d "$comp_dir" ]; then
    echo "Skipping missing comparison directory: $comp_dir"
    continue
  fi

  echo "Processing comparison: $comparison"

  shopt -s nullglob
  taxon_dirs=("$comp_dir"/*/)
  shopt -u nullglob

  if [ "${#taxon_dirs[@]}" -eq 0 ]; then
    echo "No taxon directories found in: $comp_dir"
    continue
  fi

  for taxon_path in "${taxon_dirs[@]}"; do
    taxon_name="$(basename "$taxon_path")"
    fasta_file="${taxon_path}/${taxon_name}_ASV.fasta"

    if [ ! -f "$fasta_file" ]; then
      echo "Skipping $taxon_name: missing FASTA file $fasta_file"
      continue
    fi

    echo "  BLASTing taxon: $taxon_name"

    blastn \
      -query "$fasta_file" \
      -db "$BLAST_DB_BASE" \
      -task megablast \
      -num_threads "$NUM_THREADS" \
      -max_target_seqs "$MAX_TARGET_SEQS" \
      -outfmt "$OUTFMT" \
      > "${taxon_path}/${taxon_name}_blast_hits.tsv"
  done
done

echo "BLAST complete"
