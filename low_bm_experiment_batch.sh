#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=25GB
#SBATCH --output=SLURM_stdout/slurm-%A_%a.out
#SBATCH --error=SLURM_stderr/slurm-%A_%a.err
#SBATCH --time 6:00:00

source myconda
mamba activate low-bm-base
cd /data/$USER/low-bm

BASE_CONFIG="${BASE_CONFIG:-config.yaml}"
BATCH_TABLE="${BATCH_TABLE:-experiment_batch_configs.tsv}"
RUN_CONFIG_DIR="${RUN_CONFIG_DIR:-experiment_batch_configs}"
mkdir -p "$RUN_CONFIG_DIR"

mapfile -t batch_info < <(python3 - "$BATCH_TABLE" "$SLURM_ARRAY_TASK_ID" "$RUN_CONFIG_DIR" <<'PY'
import csv
import os
import sys

batch_table, slurm_task_id, run_config_dir = sys.argv[1], int(sys.argv[2]), sys.argv[3]
row_index = slurm_task_id - 1
canonical_cols = [
    "trialID",
    "trial_descript",
    "exp_dir",
    "metadata",
    "batch_label",
    "include_processing",
    "include_analysis",
]
optional_config_cols = [
    "process_umis",
]

def truthy(value):
    return str(value).strip().lower() not in {"", "0", "false", "f", "no", "n"}

with open(batch_table, newline="") as handle:
    rows = list(csv.reader(handle, delimiter="\t"))

if not rows:
    raise SystemExit(f"Batch table is empty: {batch_table}")

has_header = rows[0][:4] == canonical_cols[:4] or rows[0][0] == "trialID"

if has_header:
    with open(batch_table, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        records = list(reader)
    if row_index < 0 or row_index >= len(records):
        raise SystemExit(f"SLURM_ARRAY_TASK_ID {slurm_task_id} is outside 1-{len(records)} batch rows.")
    fieldnames = reader.fieldnames or []
    selected_raw_row = rows[row_index + 1]
    if len(selected_raw_row) != len(fieldnames):
        raise SystemExit(
            f"Batch row {slurm_task_id} has {len(selected_raw_row)} tab-separated field(s), "
            f"but the header has {len(fieldnames)} column(s). "
            "Check for a missing value or an extra/missing tab."
        )
    missing_cols = [key for key in canonical_cols if key not in fieldnames]
    if missing_cols:
        suspicious = [
            field for field in fieldnames
            if any(key in field for key in canonical_cols + optional_config_cols)
            and field not in canonical_cols
            and field not in optional_config_cols
        ]
        hint = ""
        if suspicious:
            hint = (
                " Suspicious header value(s): "
                + ", ".join(repr(field) for field in suspicious)
                + ". Check that columns are separated by tabs, not spaces."
            )
        raise SystemExit(
            f"Missing required batch table column(s): {', '.join(missing_cols)}."
            + hint
        )
    row = {key: (records[row_index].get(key, "") or "").strip() for key in canonical_cols}
    for key in optional_config_cols:
        if key in fieldnames:
            row[key] = (records[row_index].get(key, "") or "").strip()
else:
    if row_index < 0 or row_index >= len(rows):
        raise SystemExit(f"SLURM_ARRAY_TASK_ID {slurm_task_id} is outside 1-{len(rows)} batch rows.")
    values = [field.strip() for field in rows[row_index]]
    if len(values) < 4:
        raise SystemExit("Legacy batch rows must contain trialID, trial_descript, exp_dir, metadata.")
    row = dict(zip(canonical_cols[:4], values[:4]))
    row["batch_label"] = row["trial_descript"]
    row["include_processing"] = "true"
    row["include_analysis"] = "true"

for required in canonical_cols[:4]:
    if not row.get(required):
        raise SystemExit(f"Missing required batch table value: {required}")

if not truthy(row.get("include_processing", "true")):
    print("__SKIP__")
    print(row["trialID"])
    raise SystemExit(0)

os.makedirs(run_config_dir, exist_ok=True)
run_config_file = os.path.join(run_config_dir, f"{row['trialID']}_runconfig.yaml")

def yaml_quote(value):
    value = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{value}"'

with open(run_config_file, "w", newline="\n") as out:
    config_cols = canonical_cols + [
        key for key in optional_config_cols
        if row.get(key, "") != ""
    ]
    for key in config_cols:
        default = "true" if key in {"include_processing", "include_analysis"} else ""
        out.write(f"{key}: {yaml_quote(row.get(key, default))}\n")

print(run_config_file)
print(row["trialID"])
PY
)

RUN_CONFIG_FILE="${batch_info[0]}"
trialID="${batch_info[1]}"

if [[ "$RUN_CONFIG_FILE" == "__SKIP__" ]]; then
    echo "Skipping ${trialID}: include_processing is false in ${BATCH_TABLE}."
    exit 0
fi

LOG_DIR="snakemake_logs/${trialID}"
mkdir -p "$LOG_DIR"

snakemake --use-conda --cores "${SLURM_CPUS_PER_TASK}" all\
    --configfile "$BASE_CONFIG" \
    --configfile "$RUN_CONFIG_FILE" \
    --rerun-incomplete \
    2> "$LOG_DIR/snakemake.out"
