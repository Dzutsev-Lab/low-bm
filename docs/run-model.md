# Run Model

This repository is moving from a SLURM-array batch launcher to a
baseline-style master-job launcher.

## Old Model

`low_bm_experiment_batch.sh` is submitted as a SLURM array. Each array task:

1. Selects one row from `experiment_batch_configs.tsv`.
2. Writes a per-row config file.
3. Runs that batch's whole Snakemake DAG inside the array allocation.

This isolates batches, but all rules in a batch compete inside one fixed
allocation.

## New Model

`low-bm run` is one batch. In SLURM mode it writes provenance files, submits one
master job, and that master job runs Snakemake with `profiles/slurm`. Snakemake
then submits rule-specific jobs through the Snakemake 9 SLURM executor plugin.

`low-bm batch submit` keeps the batch-table convenience by submitting one
independent master job per row in the table.

The SLURM profile still delegates rule execution from the master job to SLURM.
To reduce scheduler overhead for repeated sample-level prep work, it submits
`norm_fastq`, `umi_selection`, `umi_dedup`, and `no_umi_count_summary` through
Snakemake's SLURM array-job support rather than many independent submissions.
The profile caps each array submission at 100 tasks; if a rule has more ready
sample jobs than that, Snakemake should split them across multiple arrays.

## Baseline Files Worth Comparing

Review these files in `OpenOmics/baseline` for the closest structural analogs:

- `baseline`: CLI command dispatch and user-facing arguments.
- `src/run.py`: output initialization, config construction, dry-run, and launch flow.
- `src/run.sh`: master-job script that lets Snakemake submit rule-level jobs.
- `workflow/Snakefile`: generated-config assumptions and included rule files.
- `workflow/rules/hooks.smk`: `RUNNING`, `COMPLETED`, `FAILED`, and job-summary hooks.
- `config/cluster.json`: baseline's older per-rule resource map.
- `docs/usage/run.md`: user-facing explanation of the run command.

The concept to compare is array task runs whole DAG versus master job runs
Snakemake and Snakemake launches rule jobs.
