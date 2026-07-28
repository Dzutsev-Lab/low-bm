# HPC YAML Environment Validation

This workflow now uses repo-owned Snakemake conda environment files in
`workflow/envs/` instead of absolute paths to existing HPC environments.

The YAML files are recipes for recreating environments. They do not copy locally
installed packages from the current HPC envs. The patched `micRoclean` install is
handled by `workflow/envs/micRoclean-env.post-deploy.sh`.

## 1. Export Current HPC Envs For Audit

Run this on the HPC before building fresh Snakemake-managed envs:

```bash
mkdir -p env_exports
for env in low-bm-base bio-tools-env AmpUMI-env R-tools-env micRoclean-env kraken-env; do
  conda env export -n "$env" --no-builds > "env_exports/${env}.full.yaml"
  conda list -n "$env" --explicit > "env_exports/${env}.linux-64.pin.txt"
done
```

Use these exports to compare against `workflow/envs/*.yaml` if environment
creation fails or a package is missing.

## 2. Confirm Patched micRoclean Source Pins

This is required only if you plan to run optional
`low-bm meta decontaminate-phyloseq`. This repository tracks
`workflow/envs/micRoclean-source.env` as part of the reproducibility record. It
should contain the patched package Git URL plus fixed commit SHA values:

```bash
MICROCLEAN_GIT_URL=git@github.com:your-org/micRoclean.git
MICROCLEAN_GIT_REF=0123456789abcdef0123456789abcdef01234567
SCRUB_GIT_URL=https://github.com/Shenhav-and-Korem-labs/SCRuB.git
SCRUB_GIT_REF=0123456789abcdef0123456789abcdef01234567
```

The post-deploy script intentionally installs with `dependencies = FALSE`.
Missing `micRoclean` dependencies should be added to `micRoclean-env.yaml`
rather than installed opportunistically.

## 3. Create A Small-Batch Test Config

Initialize ignored local configs from the tracked templates:

```bash
mkdir -p config/local
cp config/templates/processing.yaml config/local/processing.yaml
cp config/templates/processing-overrides.yaml config/local/processing-overrides.yaml
cp config/templates/batch.tsv config/local/batch.tsv
```

Edit `config/local/processing.yaml` for the HPC paths and reference data. Edit
`config/local/processing-overrides.yaml` when you want isolated validation
outputs, and edit `config/local/batch.tsv` down to a one-row small batch. Use a
distinct `trialID` so the outputs do not collide with trusted runs.

The launcher will layer configs in this order:

```bash
config/local/processing.yaml
config/local/processing-overrides.yaml
experiment_batch_configs/<trialID>_runconfig.yaml
```

## 4. Dry Run And Build Envs

Create and validate the project-local runner first:

```bash
./low-bm setup runner
./low-bm doctor runner --mode slurm
```

Use `low-bm` for the normal validation path:

```bash
./low-bm batch submit \
  --batch-table config/local/batch.tsv \
  --configfile config/local/processing.yaml \
  --extra-configfile config/local/processing-overrides.yaml \
  --dry-run \
  --mode slurm
```

SLURM batch submission uses isolated Snakemake workdirs by default, rooted at
`.low-bm/snakemake-workdirs`, and shares rule conda environments through
`.low-bm/snakemake-conda`. The dry-run output should show a different
`--directory` value for each batch-table row.

If the shared Snakemake rule environments are not already built, create them
serially before submitting many rows:

```bash
./low-bm batch prepare-envs \
  --batch-table config/local/batch.tsv \
  --configfile config/local/processing.yaml \
  --extra-configfile config/local/processing-overrides.yaml
```

This runs locally with `--conda-create-envs-only` and avoids concurrent
environment creation races. It also avoids the outer `mamba run` wrapper used
by older launcher versions.

If this step reports `bad interpreter: .low-bm/runner/env/bin/python` or
`Error running conda info`, the runner env's conda entrypoint is using a
relative interpreter path that breaks once Snakemake runs from an isolated
workdir. Current `low-bm` launchers export `.low-bm/runner/conda-shell.sh`
through `BASH_ENV`, so Snakemake's bash probes call conda through absolute
runner-prefix paths. After updating, validate with `./low-bm doctor runner` and
rerun `batch prepare-envs`; do not delete mamba or conda lock files unless you
have confirmed no related processes are still running.

If you need to call Snakemake directly after the row config has been generated,
keep `--configfile` as one flag followed by all config layers, then end the
variable-length config list with `--` before the target:

```bash
snakemake --profile profiles/slurm --dry-run \
  --configfile config/local/processing.yaml \
    config/local/processing-overrides.yaml \
    experiment_batch_configs/<trialID>_runconfig.yaml \
  -- all
```

## 5. Run The Small Batch

```bash
./low-bm batch submit \
  --batch-table config/local/batch.tsv \
  --configfile config/local/processing.yaml \
  --extra-configfile config/local/processing-overrides.yaml \
  --mode slurm
```

For a single direct run, provide the batch row fields and the same config stack:

```bash
./low-bm run \
  --trial-id <trialID> \
  --trial-descript <trial_descript> \
  --exp-dir <exp_dir> \
  --metadata <metadata> \
  --configfile config/local/processing.yaml \
  --extra-configfile config/local/processing-overrides.yaml \
  --mode slurm
```

`low-bm run` keeps the shared checkout-level Snakemake workdir by default for
compatibility. Add `--isolated-workdir` when a single direct run should use the
same isolated workdir and shared conda-prefix layout as batch SLURM rows.

If a master job fails and leaves a stale isolated lock, first confirm that no
matching Snakemake or master SLURM jobs are still running, then unlock only the
affected row:

```bash
./low-bm batch unlock --trial-id <trialID>
```

Use `./low-bm batch unlock --all` only after confirming every row in the batch
table is stopped. Do not use `--nolock` as the routine fix for lock errors.
Do not delete mamba or conda lock files until you have confirmed no mamba,
conda, Snakemake, or master jobs are still running.

Expected terminal files:

- `<out_root>/<trialID>_<trial_descript>/<trialID>_physeq.RData`
- `<out_root>/<trialID>_<trial_descript>/<trialID>_ASV.fasta`
- `<out_root>/<trialID>_<trial_descript>/effective_config.yaml`
- clean rule logs under `<out_root>/<trialID>_<trial_descript>/Logs`

The phyloseq and ASV FASTA endpoints are pre-micRoclean, Kraken-annotated
processing outputs. Optional micRoclean decontamination is run later with
`low-bm meta decontaminate-phyloseq`.

## 6. Containerization Follow-Up

After the YAML-managed run succeeds, generate an Apptainer definition on the HPC:

```bash
snakemake --containerize apptainer > low-bm.def
```

Build and test the image using the same small-batch config, bind-mounting the
project directory plus external data/reference locations. Keep `Exp_Data`,
`IP_Data`, `Exp_Output`, and `Ref_Data` outside the image.
