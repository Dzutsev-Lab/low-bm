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

  require_command python3
  config_text="$(
    python3 - "$analysis_config" <<'PY'
import sys
from pathlib import Path

try:
    import yaml
except ModuleNotFoundError:
    raise SystemExit(
        "Required Python package not found: pyyaml. "
        "Recreate the low-bm-bio-tools analysis environment from workflow/envs/bio-tools-env.yaml."
    )


def is_missing(value):
    if value is None:
        return True
    if isinstance(value, list):
        return len(value) == 0 or all(is_missing(item) for item in value)
    return str(value).strip().lower() in {"", "null", "none", "~"}


def as_list(value):
    if isinstance(value, list):
        return [str(item) for item in value if not is_missing(item)]
    if is_missing(value):
        return []
    return [str(value)]


def first_present(values, keys):
    for key in keys:
        value = values.get(key)
        if not is_missing(value):
            return value
    return None


def emit(key, value):
    if is_missing(value):
        value = ""
    print(f"{key}\t{value}")


config_path = Path(sys.argv[1])
with config_path.open() as handle:
    cfg = yaml.safe_load(handle) or {}

project = cfg.get("project") or {}
blast = cfg.get("blast_confirmation")
if blast is None:
    raise SystemExit("analysis_config.yaml is missing blast_confirmation.")
if not isinstance(project, dict) or not isinstance(blast, dict):
    raise SystemExit("analysis_config.yaml project and blast_confirmation must be mappings.")

io_dir = first_present(blast, ("io_dir", "output_dir", "out_dir")) or project.get("output_dir")
comparisons = blast.get("candidate_comparisons")
if is_missing(comparisons):
    comparisons = blast.get("DA_comparisons")

missing = []
if is_missing(io_dir):
    missing.append("project.output_dir or blast_confirmation.io_dir")
if is_missing(blast.get("reference_db")):
    missing.append("reference_db")
if is_missing(comparisons):
    missing.append("blast_confirmation.candidate_comparisons or blast_confirmation.DA_comparisons")
if missing:
    raise SystemExit("blast_confirmation is missing required field(s): " + ", ".join(missing))

emit("TRIAL_DIR", str(Path(str(io_dir)) / "BlastAnalysis"))
emit("REFERENCE_DB", blast.get("reference_db"))
emit("REF_BASE_DIR", blast.get("ref_base_dir") or "Ref_Data")
emit("NUM_THREADS", blast.get("num_threads"))
emit("MAX_TARGET_SEQS", blast.get("max_target_seqs") or 5)
for comparison in as_list(comparisons):
    emit("COMPARISON", comparison)
PY
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
