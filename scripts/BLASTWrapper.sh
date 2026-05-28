#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: run_genus_blasts.sh <genus_fasta_dir> <blast_db_base>"
  echo "BLAST hits output to <genus_fasta_dir>"
  exit 1
fi

TRIAL_DIR="$1"
REFERENCE_DB="$2"
#Now that other input variables are set, can discard from input to collect comparisons
shift 2

COMPARISONS=("$@")

REF_BASE_DIR="Ref_Data"
REFERENCE_FASTA="${REF_BASE_DIR}/${REFERENCE_DB}.fasta"
TAXID_MAP="${REF_BASE_DIR}/${REFERENCE_DB}/seqid2taxid.map"
BLAST_DB_BASE="${REF_BASE_DIR}/${REFERENCE_DB}/blast_format_db/${REFERENCE_DB}_blast"
NUM_THREADS="${SLURM_CPUS_PER_TASK}"
MAX_TARGET_SEQS=5
OUTFMT='6 qseqid sseqid pident length qcovs evalue bitscore staxids'

#Build BLAST database only if missing
if [ ! -f "${REF_BASE_DIR}/${REFERENCE_DB}/blast_format_db/${REFERENCE_DB}_blast.nsq" ]; then
    echo "BLAST database not found. Building database at: ${BLAST_DB_BASE}"
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
    genus_dirs=("$comp_dir"/*/)
    shopt -u nullglob


    if [ "${#genus_dirs[@]}" -eq 0 ]; then
        echo "No genus directories found in: $comp_dir"
        continue
    fi

    for genus_path in "${genus_dirs[@]}"; do
        genus_name="$(basename "$genus_path")"
        fasta_file="${genus_path}/${genus_name}_ASV.fasta"

        if [ ! -f "$fasta_file" ]; then
            echo "Skipping $genus_name: missing FASTA file $fasta_file"
            continue
        fi

        echo "  BLASTing genus: $genus_name"

        blastn \
            -query "$fasta_file" \
            -db "$BLAST_DB_BASE" \
            -task megablast \
            -num_threads "$NUM_THREADS" \
            -max_target_seqs "$MAX_TARGET_SEQS" \
            -outfmt "$OUTFMT" \
            > "${genus_path}/${genus_name}_blast_hits.tsv"
    done
done

echo "BLAST Complete"