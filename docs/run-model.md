# Run Model

This repository is moving from a SLURM-array batch launcher to a
baseline-style master-job launcher.

## Old Model

`low_bm_experiment_batch.sh` is submitted as a SLURM array. Each array task:

1. Selects one row from the legacy root `experiment_batch_configs.tsv`.
2. Writes a per-row config file.
3. Runs that batch's whole Snakemake DAG inside the array allocation.

This isolates batches, but all rules in a batch compete inside one fixed
allocation.

## New Model

`low-bm run` is one batch. In SLURM mode it writes provenance files, submits one
master job, and that master job runs Snakemake with `profiles/slurm`. Snakemake
then submits rule-specific jobs through the Snakemake 9 SLURM executor plugin.

`low-bm batch submit` keeps the batch-table convenience by submitting one
independent master job per row in the table. By default it reads the ignored
local batch table at `config/local/batch.tsv`, which should be initialized from
the tracked template at `config/templates/batch.tsv`.

For `low-bm batch submit --mode slurm`, each row now uses an isolated
Snakemake working directory under `.low-bm/snakemake-workdirs/`. The directory
name is derived from the row's `<trialID>_<trial_descript>` output name, with
unsafe path characters replaced. The submitted Snakemake command uses absolute
paths for `--snakefile`, `--directory`, `--profile`, and all config files, so
repo-relative config values still point at this checkout while the Snakemake
lock lives in the row-specific workdir. Rule conda environments are shared
through `.low-bm/snakemake-conda` by default, preventing each isolated row from
rebuilding the same environments.

Use `--shared-workdir` with `batch submit` only for debugging or legacy
compatibility. `low-bm run` keeps the shared checkout-level Snakemake workdir
by default, but `--isolated-workdir` enables the same per-trial pathing for a
single run.

Config loading is controlled by the launcher rather than by a `configfile:`
directive in the Snakefile. The processing config stack is:

1. `config/local/processing.yaml`, initialized from
   `config/templates/processing.yaml`.
2. Any `--extra-configfile` overrides, such as an ignored
   `config/local/processing-overrides.yaml`.
3. The generated per-row run config under `experiment_batch_configs/`.

Snakemake itself is launched through a project-local runner environment created
by `low-bm setup runner`. This runner layer contains Snakemake and the SLURM
executor plugin, while rule-level bioinformatics tools remain in
Snakemake-managed conda environments. See `docs/portability.md` for the
portability rationale.

The default processing runner calls `.low-bm/runner/env/bin/snakemake`
directly and prepends `.low-bm/runner/env/bin` to `PATH` inside local runs and
master jobs. This avoids concurrent `mamba run` lock contention while keeping
Snakemake and its executor plugin pinned to the project-local runner prefix.

When rule conda environments may need to be created, prepare them serially
before launching many SLURM rows:

```bash
./low-bm batch prepare-envs
```

This uses the same batch table, config stack, isolated workdir resolver, and
shared `.low-bm/snakemake-conda` prefix as `batch submit`, but it runs locally
with `--conda-create-envs-only`.

## Unlocking Stale Locks

Do not use `--nolock` as the normal fix for concurrent batch jobs. If a
Snakemake master job dies and leaves a stale lock behind:

1. Confirm no relevant Snakemake or master SLURM jobs are still running.
2. Unlock only the affected isolated workdir, for example:

   ```bash
   ./low-bm batch unlock --trial-id 072826.6
   ```

3. Resubmit the failed row or batch.

Use `./low-bm batch unlock --all` only when every row in the batch table is
known to be stopped. `low-bm run --unlock` remains available for the single-run
path and unlocks whichever workdir that run resolves to, shared by default or
isolated when `--isolated-workdir` is supplied.

Before launching many isolated SLURM rows concurrently, make sure shared BWA
reference indexes such as `<reference>.bwt` already exist. Isolated Snakemake
workdirs prevent lock contention between batches, which also means they no
longer coordinate first-time creation of shared reference index files.

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

Post-processing analyses are handled outside the processing DAG through
`low-bm analysis`; see `docs/analysis-run-model.md`. Multi-batch compilation
and result-level meta-analyses live under `low-bm meta`; see
`docs/meta-run-model.md`.

The processing DAG's canonical `<trialID>_physeq.RData` endpoint is
pre-micRoclean and Kraken-annotated. Optional micRoclean decontamination now
runs through `low-bm meta decontaminate-phyloseq` to support deliberate pre/post
comparison endpoints.
