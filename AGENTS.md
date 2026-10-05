# AGENTS.md

## Repository Map

- `Snakefile` is the processing DAG. It should stay focused on processing 
  workflow graph, paths, rule inputs/outputs, and Snakemake integration.
- `./low-bm` is the stable user-facing launcher for setup, processing runs,
  reference checks, analysis, meta-analysis, and runner diagnostics.
- `src/low_bm/` contains the Python CLI/orchestration layer. It prepares config
  stacks, provenance, dry-runs, local execution, and SLURM master-job commands.
- `scripts/` contains R, Python, and shell analysis/workflow scripts. Shared R
  helpers live under `scripts/Rhelpers/`.
- `workflow/envs/` contains repo-owned conda environment recipes and pinned
  environment records.
- `config/templates/` contains tracked config templates. `config/local/` is
  ignored local configuration.

Start with these docs before changing related behavior:

- Processing, batch submission, config layering, and unlock behavior:
  `docs/run-model.md`
- Post-processing analysis commands and ANCOM-BC2 conventions:
  `docs/analysis-run-model.md`
- Multi-batch meta-analysis and optional micRoclean flow:
  `docs/meta-run-model.md`
- Runner, rule environment, lock, and portability model:
  `docs/portability.md`
- HPC environment validation and small-batch dry-run workflow:
  `docs/hpc-yaml-environments.md`

## Operating Rules

- Prefer `./low-bm` over direct script or Snakemake invocation unless debugging
  the underlying layer is the task.
- Preserve config layering order: base processing config, extra override files,
  then generated per-row run config.
- Keep the runtime layers separate: the launcher is Python, the runner env runs
  Snakemake, and rule envs are Snakemake-managed conda environments.
- Treat scientific defaults and comparison ordering as intentional. Do not
  change analysis formulas, `factor_levels`, `ordered_levels`, decontamination
  defaults, reference selection, or taxonomy behavior unless the task asks for
  that scientific change.
- Treat this checkout as a tool-development environment only. Do not run real
  biological analyses, production processing, meta-analysis, reference builds,
  or other data-bearing workflows locally, even when explicitly requested.
  Local validation must use unit tests, dry-runs, command/configuration checks,
  mocks, fixtures, or synthetic non-sensitive data only. Do not access, copy,
  persist, or export real sample-level data, patient data, credentials, or
  analysis outputs here. Full analyses and validation involving real data must
  be performed in a separate execution environment after changes are pushed to
  the shared remote repository. If requested validation requires real data or
  a substantive analysis, stop and report that boundary instead of running it.
- Do not launch real SLURM jobs, build large reference indexes, or run full
  biological workflows locally. Do not remove Snakemake/conda locks unless
  explicitly asked and the documented unlock conditions are satisfied. Prefer
  dry-runs and diagnostics while developing.
- Do not use `--nolock` as a routine fix. Follow the documented unlock flow and
  confirm related jobs are stopped before unlocking.
- Never commit secrets, personal absolute paths, HPC credentials, patient data,
  or sample-level data exports.
- Leave ignored local/run-state paths uncommitted unless explicitly requested:
  `.low-bm/`, `config/local/`, `experiment_batch_configs/`, `snakemake_logs/`,
  `analysis_logs/`, `meta_logs/`, `reference_logs/`, `SLURM_stdout/`,
  `SLURM_stderr/`, `Exp_Data/`, `IP_Data/`, `Exp_Output/`, and `Ref_Data/`.

## Validation

- Python CLI or helper changes: run
  `python -m unittest discover -s tests -p 'test_*.py'`.
- R helper or analysis-script changes: run `Rscript tests/test_helpers.R` when
  the required R packages are available.
- Shell wrapper changes: run the targeted Python tests when present, such as
  `python -m unittest tests/test_blast_wrapper.py`.
- Workflow/config changes: use a `./low-bm ... --dry-run` path with local
  configuration and synthetic or fixture inputs only. Do not execute real
  biological workflows or analyses in this development checkout.
- Runner/environment changes: use `./low-bm doctor runner --mode local`; add
  `--mode slurm` or `--rule-env-smoke-test` only when relevant and available.
- Documentation-only changes normally need review of the diff and linked docs;
  do not run large workflow checks just to validate prose.

## Code Review Rules

- Flag accidental writes to ignored run-state, reference, log, or output
  directories.
- Flag changes that bypass `./low-bm` command construction, provenance writing,
  config stack ordering, isolated workdir behavior, or runner/rule env
  separation.
- Flag scientific changes that lack an explicit rationale, especially comparison
  ordering, ANCOM-BC2 settings, decontamination controls, host/reference
  selection, ASV filtering, taxonomy reconciliation, or metadata joins.
- Flag real-data execution, substantive analysis commands, or analysis outputs
  in this checkout; prefer dry-runs, mocks, fixtures, synthetic data, and
  isolated tests.
