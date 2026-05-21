#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=25GB
#SBATCH --output=SLURM_stdout/slurm-%A_%a.out
#SBATCH --error=SLURM_stderr/slurm-%A_%a.err
#SBATCH --time 5:00:00

source myconda
mamba activate low-bm-base
cd /data/$USER/low-bm

line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" experiment_batch_configs.tsv)
IFS=$'\t' read -r trialID trial_descript exp_dir metadata <<< "$line"
RUN_CONFIG_DIR="experiment_batch_configs/"
mkdir -p "$RUN_CONFIG_DIR"
RUN_CONFIG_FILE="${RUN_CONFIG_DIR}/${trialID}_runconfig.yaml"

cat > "$RUN_CONFIG_FILE" <<EOF
trialID: "$trialID"
trial_descript: "$trial_descript"
exp_dir: "$exp_dir"
metadata: "$metadata"
EOF

LOG_DIR="snakemake_logs/${trialID}"
mkdir -p "$LOG_DIR"

snakemake --use-conda --cores "${SLURM_CPUS_PER_TASK}" all\
    --configfile "$RUN_CONFIG_FILE" \
    --rerun-incomplete \
    2> "$LOG_DIR/snakemake.out"