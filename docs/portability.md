# Portability Notes

`low-bm` now treats runtime portability as layered rather than monolithic. That
keeps the implementation easier to debug today and leaves a clean path toward
containerization later.

## Three Runtime Layers

1. **Launcher layer**: the `low-bm` Python CLI. It reads batch tables, prepares
   config stacks, writes provenance, and submits or runs Snakemake.
2. **Runner layer**: the environment that runs Snakemake itself, including the
   SLURM executor plugin. V1 creates this at `.low-bm/runner/env`.
3. **Rule environment layer**: the Snakemake-managed conda environments in
   `workflow/envs/`. These contain bioinformatics and R tooling for individual
   rules.

Keeping these layers separate is the portability trick. The runner can be made
reproducible without collapsing every rule dependency into one giant env, and a
future container runner can replace only the runner layer first.

## Batch Workdirs And Locks

Concurrent `low-bm batch submit --mode slurm` runs use isolated Snakemake
working directories by default:

```bash
./low-bm batch submit \
  --workdir-root .low-bm/snakemake-workdirs \
  --snakemake-conda-prefix .low-bm/snakemake-conda
```

Each batch-table row gets its own lock scope under `--workdir-root`, while
`--snakemake-conda-prefix` keeps rule environments shared across rows. The
workflow still writes biological outputs to the configured `IP_Data` and
`Exp_Output` roots; the isolated workdir is for Snakemake metadata, locks, and
runtime bookkeeping.

Before launching many rows whose rule environments have not been built yet, run
the environment preparation step serially:

```bash
./low-bm batch prepare-envs
```

This creates the Snakemake-managed conda environments in the shared
`.low-bm/snakemake-conda` prefix without submitting rule jobs to SLURM.

Shared BWA reference indexes are prepared outside the processing DAG:

```bash
./low-bm references prepare-bwa-indexes
./low-bm references check
```

The processing DAG requires the existing `.amb`, `.ann`, `.bwt`, `.pac`, and
`.sa` sidecars for each configured host, viral, and bacterial reference. This
keeps isolated batch masters from racing over shared reference files.

Use `--shared-workdir` on `batch submit` to recover the legacy checkout-level
lock behavior for debugging. Use `--isolated-workdir` on `low-bm run` when a
single run should use the same isolated layout.

Unlock stale isolated locks locally with:

```bash
./low-bm batch unlock --trial-id <trialID>
```

Only unlock after confirming no matching Snakemake or master jobs are still
running.

## Why A Project-Local Runner Prefix?

Global conda env names such as `low-bm-runner` are convenient, but they are also
easy to mutate accidentally or reuse across unrelated checkouts. A prefix such
as `.low-bm/runner/env` belongs to this repository copy, so two pipeline
versions can carry different runner environments side by side.

Create the runner with:

```bash
./low-bm setup runner
```

Validate it with:

```bash
./low-bm doctor runner --mode local
./low-bm doctor runner --mode slurm
./low-bm doctor runner --rule-env-smoke-test
```

## Why Direct Runner Entrypoints?

Shell activation depends on startup files, shell type, and site-specific module
state. `mamba run --prefix ...` also creates a process lock under the user's
global mamba cache, which can fail when many master jobs start at once on a
shared filesystem. The processing launcher therefore runs the pinned entrypoint
directly:

```bash
.low-bm/runner/env/bin/snakemake ...
```

The launcher calls that Snakemake executable by absolute path, then gives
Snakemake a conda-base-first environment and passes
`--conda-base-path <real-conda-base>`. The conda base is detected during
`low-bm setup runner` with `<manager> info --base` and recorded beside the
runner metadata. If the HPC needs a specific base path, provide it explicitly:

```bash
./low-bm setup runner --conda-base-prefix /path/to/conda/base
```

This keeps the boundary simple: the runner env starts Snakemake, and the real
conda base activates Snakemake-managed rule environments under
`.low-bm/snakemake-conda`. The launcher does not rewrite
`.low-bm/runner/env/bin/activate`, and normal runs do not depend on generated
activation shims. If `doctor runner` reports missing conda base activation
support, rerun `setup runner` with `--conda-base-prefix`; rebuilding the runner
env is only needed when `doctor runner` reports missing runner executables. Use
`doctor runner --rule-env-smoke-test` when you need to confirm that Snakemake
rule jobs use `.low-bm/snakemake-conda` instead of the runner Python.

`--activate-command` remains available as an advanced submitted-SLURM fallback,
but the direct runner path is the default lock-safe processing path.

## How This Leads To Containers

The inner Snakemake command is still built independently from the runner. Today
the host runner wraps it like this:

```bash
.low-bm/runner/env/bin/snakemake --conda-base-path /path/to/conda/base ...
```

A later container runner can wrap the same inner command like this:

```bash
apptainer exec low-bm-runner.sif snakemake ...
```

That keeps containerization as a new runner implementation rather than a rewrite
of batch parsing, config stacking, provenance, or SLURM submission logic.
