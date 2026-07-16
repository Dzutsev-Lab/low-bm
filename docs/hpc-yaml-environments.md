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

## 2. Configure Patched micRoclean

Create an untracked source config from the example:

```bash
cp workflow/envs/micRoclean-source.env.example workflow/envs/micRoclean-source.env
```

Edit `workflow/envs/micRoclean-source.env` with the patched package Git URL and
fixed commit SHA:

```bash
MICROCLEAN_GIT_URL=git@github.com:your-org/micRoclean.git
MICROCLEAN_GIT_REF=0123456789abcdef0123456789abcdef01234567
```

The post-deploy script intentionally installs with `dependencies = FALSE`.
Missing `micRoclean` dependencies should be added to `micRoclean-env.yaml`
rather than installed opportunistically.

## 3. Create A Small-Batch Test Config

Copy the isolated root overrides:

```bash
cp config.hpc-yamltest.example.yaml config.hpc-yamltest.yaml
```

Create a one-row run config from an existing small batch, or use the batch
wrapper to generate one from a one-row `experiment_batch_configs.tsv`. Use a
distinct `trialID` so the outputs do not collide with trusted runs.

The SLURM wrapper accepts two optional variables for this validation path:

```bash
EXTRA_CONFIGFILES=config.hpc-yamltest.yaml
SNAKEMAKE_DEPLOY_ARGS="--sdm conda"
```

If the HPC Snakemake is older, leave `SNAKEMAKE_DEPLOY_ARGS` unset; the wrapper
defaults to `--use-conda`.

## 4. Dry Run And Build Envs

Use Snakemake's storage-deployment syntax when available:

```bash
snakemake -n --sdm conda --cores 8 all \
  --configfile config.yaml \
  --configfile config.hpc-yamltest.yaml \
  --configfile experiment_batch_configs/yamltest_runconfig.yaml

snakemake --sdm conda --conda-create-envs-only --cores 1 all \
  --configfile config.yaml \
  --configfile config.hpc-yamltest.yaml \
  --configfile experiment_batch_configs/yamltest_runconfig.yaml
```

If the HPC Snakemake version does not recognize `--sdm conda`, use:

```bash
snakemake -n --use-conda --cores 8 all \
  --configfile config.yaml \
  --configfile config.hpc-yamltest.yaml \
  --configfile experiment_batch_configs/yamltest_runconfig.yaml

snakemake --use-conda --conda-create-envs-only --cores 1 all \
  --configfile config.yaml \
  --configfile config.hpc-yamltest.yaml \
  --configfile experiment_batch_configs/yamltest_runconfig.yaml
```

## 5. Run The Small Batch

```bash
snakemake --sdm conda --cores 8 all \
  --configfile config.yaml \
  --configfile config.hpc-yamltest.yaml \
  --configfile experiment_batch_configs/yamltest_runconfig.yaml \
  --rerun-incomplete
```

Older Snakemake:

```bash
snakemake --use-conda --cores 8 all \
  --configfile config.yaml \
  --configfile config.hpc-yamltest.yaml \
  --configfile experiment_batch_configs/yamltest_runconfig.yaml \
  --rerun-incomplete
```

Expected terminal files:

- `<out_root>/<trialID>_<trial_descript>/<trialID>_physeq.RData`
- `<out_root>/<trialID>_<trial_descript>/<trialID>_ASV.fasta`
- `<out_root>/<trialID>_<trial_descript>/effective_config.yaml`
- clean rule logs under `<out_root>/<trialID>_<trial_descript>/Logs`

## 6. Containerization Follow-Up

After the YAML-managed run succeeds, generate an Apptainer definition on the HPC:

```bash
snakemake --containerize apptainer > low-bm.def
```

Build and test the image using the same small-batch config, bind-mounting the
project directory plus external data/reference locations. Keep `Exp_Data`,
`IP_Data`, `Exp_Output`, and `Ref_Data` outside the image.
