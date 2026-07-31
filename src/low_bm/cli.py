"""Command-line entry point for low-bm workflow execution.

This module is the small orchestration layer between a human-friendly run
command and Snakemake itself. It does not know the workflow DAG; it prepares the
right config files, selects a Snakemake profile, records provenance, and either
runs Snakemake directly or submits a lightweight SLURM master job.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .batch import (
    BatchRow,
    VALID_HOSTS,
    load_simple_run_config,
    normalize_host,
    read_batch_table,
    write_run_config,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BASE_CONFIG = "config/local/processing.yaml"
DEFAULT_BATCH_TABLE = "config/local/batch.tsv"
DEFAULT_ANALYSIS_CONFIG = "config/local/analysis.yaml"
DEFAULT_META_CONFIG = "config/local/meta.yaml"
TEMPLATE_BASE_CONFIG = "config/templates/processing.yaml"
TEMPLATE_BATCH_TABLE = "config/templates/batch.tsv"
TEMPLATE_ANALYSIS_CONFIG = "config/templates/analysis.yaml"
TEMPLATE_META_CONFIG = "config/templates/meta.yaml"
DEFAULT_RUN_CONFIG_DIR = "experiment_batch_configs"
DEFAULT_LOG_ROOT = "snakemake_logs"
DEFAULT_ANALYSIS_LOG_ROOT = "analysis_logs"
DEFAULT_META_LOG_ROOT = "meta_logs"
DEFAULT_REFERENCE_LOG_ROOT = "reference_logs"
DEFAULT_ANALYSIS_ENV_ROOT = ".low-bm/analysis/envs"
DEFAULT_RUNNER_PREFIX = ".low-bm/runner/env"
DEFAULT_RUNNER_ENV_FILE = "workflow/envs/runner-env.yaml"
DEFAULT_SNAKEMAKE_WORKDIR_ROOT = ".low-bm/snakemake-workdirs"
DEFAULT_SNAKEMAKE_CONDA_PREFIX = ".low-bm/snakemake-conda"
ENV_MANAGERS = ("mamba", "conda", "micromamba")
GIT_SHA_RE = re.compile(r"^[0-9a-fA-F]{7,40}$")
REQUIRED_PROCESSING_CONFIG_KEYS = ("in_root", "ip_root", "out_root", "script_dir")
ANALYSIS_ENV_MODES = ("managed", "named", "prefix", "direct")
BWA_INDEX_EXTENSIONS = (".amb", ".ann", ".bwt", ".pac", ".sa")
BWA_INDEX_MANIFEST_SCHEMA_VERSION = 1


@dataclass
class RunSpec:
    """Resolved inputs for one Snakemake invocation.

    Argument parsing leaves values in a CLI-shaped Namespace. RunSpec is the
    normalized, execution-shaped object used after row config generation and
    default resolution are complete.
    """

    configfiles: list[Path]
    target: str
    profile: Path
    mode: str
    log_dir: Path
    trial_id: str
    job_name: str
    dry_run: bool
    unlock: bool
    extra_snakemake_args: list[str]
    trial_name: str | None = None
    workdir: Path | None = None
    snakefile: Path | None = None
    snakemake_conda_prefix: Path | None = None


@dataclass(frozen=True)
class AnalysisStepSpec:
    """One analysis operation exposed through the low-bm CLI."""

    name: str
    argv_template: tuple[str, ...]
    required_section: str | None
    env_file: Path | None


@dataclass(frozen=True)
class AnalysisCommandSpec:
    """Resolved command for one analysis step before and after env wrapping."""

    step: str
    command: list[str]
    execution_command: list[str]
    env_name: str | None
    env_prefix: Path | None
    env_file: Path | None


@dataclass(frozen=True)
class AnalysisRunSpec:
    """Resolved inputs for one explicit analysis invocation."""

    analysis_config: Path
    requested_steps: list[str]
    expanded_steps: list[str]
    log_dir: Path
    dry_run: bool
    env_mode: str
    manager: EnvManagerSpec | None
    commands: list[AnalysisCommandSpec]


@dataclass(frozen=True)
class EnvManagerSpec:
    """The conda-compatible executable used to create and run the host runner."""

    name: str
    executable: str


@dataclass(frozen=True)
class HostRunnerSpec:
    """Resolved host runner state.

    The runner env is the controller environment: it runs Snakemake and the
    SLURM executor plugin. Rule environments remain separate Snakemake-managed
    conda envs under workflow/envs/. Keeping those layers separate is the first
    portability lesson here: the scheduler/controller can be made reproducible
    without merging every bioinformatics tool into one environment.
    """

    prefix: Path
    manager: EnvManagerSpec
    metadata_path: Path
    bin_dir: Path
    conda_base_prefix: Path


@dataclass(frozen=True)
class BwaReferenceSpec:
    """One configured reference FASTA and its processing config key."""

    key: str
    path: Path


R_TOOLS_ENV = Path("workflow/envs/R-tools-env.yaml")
BIO_TOOLS_ENV = Path("workflow/envs/bio-tools-env.yaml")
MICROCLEAN_ENV = Path("workflow/envs/micRoclean-env.yaml")
R_TOOLS_ENV_PREFIX_ENVVAR = "LOW_BM_R_TOOLS_PREFIX"
BIO_TOOLS_ENV_PREFIX_ENVVAR = "LOW_BM_BIO_TOOLS_PREFIX"
MICROCLEAN_ENV_PREFIX_ENVVAR = "LOW_BM_MICROCLEAN_PREFIX"

ANALYSIS_STEP_REGISTRY: dict[str, AnalysisStepSpec] = {
    "ordination": AnalysisStepSpec(
        name="ordination",
        argv_template=("Rscript", "scripts/OrdinationAnalysis.R", "--analysis-config", "{analysis_config}"),
        required_section="ordination",
        env_file=R_TOOLS_ENV,
    ),
    "abundance-barplots": AnalysisStepSpec(
        name="abundance-barplots",
        argv_template=("Rscript", "scripts/AbundanceBarPlots.R", "--analysis-config", "{analysis_config}"),
        required_section="abundance_barplots",
        env_file=R_TOOLS_ENV,
    ),
    "differential-abundance": AnalysisStepSpec(
        name="differential-abundance",
        argv_template=("Rscript", "scripts/DiffAbundAnalysis.R", "--analysis-config", "{analysis_config}"),
        required_section="differential_abundance",
        env_file=R_TOOLS_ENV,
    ),
    "lefse": AnalysisStepSpec(
        name="lefse",
        argv_template=("Rscript", "scripts/LEfSeAnalysis.R", "--analysis-config", "{analysis_config}"),
        required_section="lefse_analysis",
        env_file=R_TOOLS_ENV,
    ),
    "heatmap-violin": AnalysisStepSpec(
        name="heatmap-violin",
        argv_template=("Rscript", "scripts/HeatMapsViolinPlots.R", "--analysis-config", "{analysis_config}"),
        required_section="heatmap_violin",
        env_file=R_TOOLS_ENV,
    ),
    "xgboost": AnalysisStepSpec(
        name="xgboost",
        argv_template=("Rscript", "scripts/XGBoostClassification.R", "--analysis-config", "{analysis_config}"),
        required_section="xgboost_classification",
        env_file=R_TOOLS_ENV,
    ),
    "survival": AnalysisStepSpec(
        name="survival",
        argv_template=("Rscript", "scripts/SurvivalAnalysis.R", "--analysis-config", "{analysis_config}"),
        required_section="survival_analysis",
        env_file=R_TOOLS_ENV,
    ),
    "blast-candidates": AnalysisStepSpec(
        name="blast-candidates",
        argv_template=("Rscript", "scripts/BLASTCandidatePreparation.R", "--analysis-config", "{analysis_config}"),
        required_section="blast_confirmation",
        env_file=R_TOOLS_ENV,
    ),
    "blast-search": AnalysisStepSpec(
        name="blast-search",
        argv_template=("bash", "scripts/BLASTWrapper.sh", "--analysis-config", "{analysis_config}"),
        required_section="blast_confirmation",
        env_file=BIO_TOOLS_ENV,
    ),
    "blast-plots": AnalysisStepSpec(
        name="blast-plots",
        argv_template=("Rscript", "scripts/BLASTPlotting.R", "--analysis-config", "{analysis_config}"),
        required_section="blast_confirmation",
        env_file=R_TOOLS_ENV,
    ),
}

ANALYSIS_COMPOSITES: dict[str, tuple[str, ...]] = {
    "blast-confirmation": ("blast-candidates", "blast-search", "blast-plots"),
}

ANALYSIS_STEP_CHOICES = tuple(sorted(set(ANALYSIS_STEP_REGISTRY) | set(ANALYSIS_COMPOSITES)))
ANALYSIS_ENDPOINT_STEPS = frozenset(
    set(ANALYSIS_STEP_REGISTRY) - {"blast-search", "blast-plots"}
)

META_STEP_REGISTRY: dict[str, AnalysisStepSpec] = {
    "compile-phyloseq": AnalysisStepSpec(
        name="compile-phyloseq",
        argv_template=("Rscript", "scripts/PhyloseqCompiler.R", "--analysis-config", "{analysis_config}"),
        required_section="meta_compile",
        env_file=R_TOOLS_ENV,
    ),
    "decontaminate-phyloseq": AnalysisStepSpec(
        name="decontaminate-phyloseq",
        argv_template=("Rscript", "scripts/PhyloseqDecontamination.R", "--analysis-config", "{analysis_config}"),
        required_section="meta_decontamination",
        env_file=MICROCLEAN_ENV,
    ),
    "differential-abundance": AnalysisStepSpec(
        name="differential-abundance",
        argv_template=("Rscript", "scripts/DiffAbundMetaAnalysis.R", "--analysis-config", "{analysis_config}"),
        required_section="meta_differential_abundance",
        env_file=R_TOOLS_ENV,
    ),
}

META_STEP_CHOICES = tuple(sorted(META_STEP_REGISTRY))


def main(argv: list[str] | None = None) -> int:
    """Parse command-line arguments and dispatch to the selected subcommand."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level CLI with one-batch and batch-table entry points."""
    parser = argparse.ArgumentParser(
        prog="low-bm",
        description="Run and submit the low-biomass Snakemake workflow.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # `low-bm run` is the primitive operation: one batch, one Snakemake call.
    run_parser = subparsers.add_parser("run", help="Run or submit one batch.")
    add_run_arguments(run_parser)
    run_parser.set_defaults(func=run_command)

    # `low-bm batch submit` is a convenience wrapper around repeated `run`
    # calls. Each table row still becomes an independent master job.
    batch_parser = subparsers.add_parser("batch", help="Batch-table utilities.")
    batch_subparsers = batch_parser.add_subparsers(dest="batch_command", required=True)
    submit_parser = batch_subparsers.add_parser(
        "submit",
        help="Submit one independent master job per batch-table row.",
    )
    add_batch_submit_arguments(submit_parser)
    submit_parser.set_defaults(func=batch_submit_command)

    prepare_envs_parser = batch_subparsers.add_parser(
        "prepare-envs",
        help="Create Snakemake rule conda environments before concurrent batch runs.",
    )
    add_batch_prepare_envs_arguments(prepare_envs_parser)
    prepare_envs_parser.set_defaults(func=batch_prepare_envs_command)

    unlock_parser = batch_subparsers.add_parser(
        "unlock",
        help="Unlock stale Snakemake locks for batch-table workdir(s).",
    )
    add_batch_unlock_arguments(unlock_parser)
    unlock_parser.set_defaults(func=batch_unlock_command)

    references_parser = subparsers.add_parser("references", help="Check and prepare shared reference assets.")
    references_subparsers = references_parser.add_subparsers(dest="references_command", required=True)
    references_check_parser = references_subparsers.add_parser(
        "check",
        help="Validate configured BWA reference indexes.",
    )
    add_references_check_arguments(references_check_parser)
    references_check_parser.set_defaults(func=references_check_command)

    references_prepare_parser = references_subparsers.add_parser(
        "prepare-bwa-indexes",
        help="Create or refresh configured BWA reference indexes.",
    )
    add_references_prepare_arguments(references_prepare_parser)
    references_prepare_parser.set_defaults(func=references_prepare_bwa_indexes_command)

    analysis_parser = subparsers.add_parser("analysis", help="Run explicit post-processing analyses.")
    analysis_subparsers = analysis_parser.add_subparsers(dest="analysis_command", required=True)
    analysis_init_parser = analysis_subparsers.add_parser(
        "init",
        help="Create an editable local analysis config from the template.",
    )
    add_analysis_init_arguments(analysis_init_parser)
    analysis_init_parser.set_defaults(func=analysis_init_command)

    analysis_validate_parser = analysis_subparsers.add_parser(
        "validate",
        help="Validate an analysis config and optional selected step sections.",
    )
    add_analysis_validate_arguments(analysis_validate_parser)
    analysis_validate_parser.set_defaults(func=analysis_validate_command)

    analysis_run_parser = analysis_subparsers.add_parser(
        "run",
        help="Run one or more explicit analysis steps.",
    )
    add_analysis_run_arguments(analysis_run_parser)
    analysis_run_parser.set_defaults(func=analysis_run_command)

    meta_parser = subparsers.add_parser("meta", help="Prepare and run multi-batch meta-analysis steps.")
    meta_subparsers = meta_parser.add_subparsers(dest="meta_command", required=True)
    meta_init_parser = meta_subparsers.add_parser(
        "init",
        help="Create an editable local meta-analysis config from the template.",
    )
    add_meta_init_arguments(meta_init_parser)
    meta_init_parser.set_defaults(func=meta_init_command)

    meta_validate_parser = meta_subparsers.add_parser(
        "validate",
        help="Validate a meta-analysis config.",
    )
    add_meta_validate_arguments(meta_validate_parser)
    meta_validate_parser.set_defaults(func=meta_validate_command)

    meta_compile_parser = meta_subparsers.add_parser(
        "compile-phyloseq",
        help="Compile multiple processed batches into one phyloseq/ASV endpoint.",
    )
    add_meta_step_arguments(meta_compile_parser)
    meta_compile_parser.set_defaults(func=meta_compile_phyloseq_command)

    meta_decontam_parser = meta_subparsers.add_parser(
        "decontaminate-phyloseq",
        help="Run optional micRoclean decontamination on one phyloseq endpoint.",
    )
    add_meta_step_arguments(meta_decontam_parser)
    meta_decontam_parser.set_defaults(func=meta_decontaminate_phyloseq_command)

    meta_da_parser = meta_subparsers.add_parser(
        "differential-abundance",
        help="Run differential-abundance meta-analysis across batch DA outputs.",
    )
    add_meta_step_arguments(meta_da_parser)
    meta_da_parser.set_defaults(func=meta_differential_abundance_command)

    setup_parser = subparsers.add_parser("setup", help="Create project-local runtime assets.")
    setup_subparsers = setup_parser.add_subparsers(dest="setup_command", required=True)
    setup_runner_parser = setup_subparsers.add_parser(
        "runner",
        help="Create or validate the project-local Snakemake runner environment.",
    )
    add_setup_runner_arguments(setup_runner_parser)
    setup_runner_parser.set_defaults(func=setup_runner_command)

    doctor_parser = subparsers.add_parser("doctor", help="Check runtime portability setup.")
    doctor_subparsers = doctor_parser.add_subparsers(dest="doctor_command", required=True)
    doctor_runner_parser = doctor_subparsers.add_parser(
        "runner",
        help="Validate the project-local Snakemake runner environment.",
    )
    add_doctor_runner_arguments(doctor_runner_parser)
    doctor_runner_parser.set_defaults(func=doctor_runner_command)

    return parser


def add_processing_config_stack_arguments(parser: argparse.ArgumentParser) -> None:
    """Add processing config-layering options."""
    # Snakemake applies later --configfile entries as overrides, so the row
    # config is appended after the base config in build_run_spec().
    parser.add_argument(
        "--configfile",
        action="append",
        default=[],
        help=(
            "Base Snakemake config file. May be repeated. "
            f"Defaults to {DEFAULT_BASE_CONFIG}."
        ),
    )
    parser.add_argument(
        "--extra-configfile",
        action="append",
        default=[],
        help="Additional config override file, layered after --configfile.",
    )


def add_common_config_arguments(parser: argparse.ArgumentParser) -> None:
    """Add config-layering options shared by single-row and table runs."""
    add_processing_config_stack_arguments(parser)
    parser.add_argument(
        "--run-config-dir",
        default=DEFAULT_RUN_CONFIG_DIR,
        help="Directory for generated per-batch run config files.",
    )


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the single-batch run path."""
    add_common_config_arguments(parser)
    add_runner_arguments(parser)
    # A run can either consume a row config already written from a batch table,
    # or generate that same row config from direct command-line fields.
    parser.add_argument(
        "--row-config",
        help="Per-row config file generated from the batch table.",
    )
    parser.add_argument("--trial-id", help="Trial ID for direct one-batch runs.")
    parser.add_argument("--trial-descript", help="Trial description for direct one-batch runs.")
    parser.add_argument("--exp-dir", help="Experiment directory for direct one-batch runs.")
    parser.add_argument("--metadata", help="Metadata path for direct one-batch runs.")
    parser.add_argument(
        "--host",
        choices=VALID_HOSTS,
        help="Host reference selection for direct one-batch runs.",
    )
    parser.add_argument(
        "--process-umis",
        choices=["true", "false"],
        help="Optional process_umis override for direct one-batch runs.",
    )
    parser.add_argument(
        "--mode",
        choices=["local", "slurm"],
        default="local",
        help="Run locally now or submit a SLURM master job.",
    )
    parser.add_argument("--profile", help="Snakemake profile path. Defaults to profiles/<mode>.")
    parser.add_argument("--target", default="all", help="Snakemake target rule or file.")
    parser.add_argument("--log-root", default=DEFAULT_LOG_ROOT, help="Root for master logs.")
    parser.add_argument("--job-name", help="Name for the SLURM master job.")
    parser.add_argument("--dry-run", action="store_true", help="Run a Snakemake dry-run.")
    parser.add_argument(
        "--unlock",
        action="store_true",
        help=(
            "Run snakemake --unlock for the resolved Snakemake workdir. "
            "Use only after confirming no matching Snakemake/master jobs are running."
        ),
    )
    parser.add_argument(
        "--isolated-workdir",
        action="store_true",
        help=(
            "Run Snakemake in an isolated per-trial workdir under --workdir-root. "
            "The default for low-bm run remains the shared checkout workdir."
        ),
    )
    add_snakemake_workdir_arguments(parser)
    parser.add_argument(
        "--print-command",
        action="store_true",
        help="Print the resolved Snakemake or sbatch command before running.",
    )
    parser.add_argument(
        "--snakemake-arg",
        action="append",
        default=[],
        help="Extra raw argument passed to Snakemake. May be repeated.",
    )
    add_master_job_arguments(parser)


def add_batch_submit_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for submitting all rows in a canonical batch table."""
    add_common_config_arguments(parser)
    add_runner_arguments(parser)
    parser.add_argument(
        "--batch-table",
        default=DEFAULT_BATCH_TABLE,
        help=f"Canonical batch TSV. Defaults to {DEFAULT_BATCH_TABLE}.",
    )
    parser.add_argument(
        "--mode",
        choices=["slurm", "local"],
        default="slurm",
        help="Submit SLURM master jobs or run each batch locally in sequence.",
    )
    parser.add_argument("--profile", help="Snakemake profile path. Defaults to profiles/<mode>.")
    parser.add_argument("--target", default="all", help="Snakemake target rule or file.")
    parser.add_argument("--log-root", default=DEFAULT_LOG_ROOT, help="Root for master logs.")
    parser.add_argument(
        "--shared-workdir",
        action="store_true",
        help=(
            "Use the checkout-level Snakemake workdir instead of the default "
            "isolated per-row workdirs for SLURM batch submissions."
        ),
    )
    add_snakemake_workdir_arguments(parser)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write row configs and print what would be submitted or run.",
    )
    parser.add_argument(
        "--print-command",
        action="store_true",
        help="Print each resolved Snakemake or sbatch command.",
    )
    parser.add_argument(
        "--snakemake-arg",
        action="append",
        default=[],
        help="Extra raw argument passed to Snakemake. May be repeated.",
    )
    add_master_job_arguments(parser)


def add_batch_unlock_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for unlocking batch-row Snakemake workdirs."""
    add_common_config_arguments(parser)
    add_runner_arguments(parser)
    parser.add_argument(
        "--batch-table",
        default=DEFAULT_BATCH_TABLE,
        help=f"Canonical batch TSV. Defaults to {DEFAULT_BATCH_TABLE}.",
    )
    parser.add_argument(
        "--trial-id",
        action="append",
        default=[],
        help="Trial ID to unlock. May be repeated. Mutually exclusive with --all.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Unlock every row in the batch table.",
    )
    parser.add_argument(
        "--mode",
        choices=["slurm", "local"],
        default="slurm",
        help="Profile mode to resolve for the unlock command. Defaults to slurm.",
    )
    parser.add_argument("--profile", help="Snakemake profile path. Defaults to profiles/<mode>.")
    parser.add_argument("--target", default="all", help=argparse.SUPPRESS)
    parser.add_argument("--log-root", default=DEFAULT_LOG_ROOT, help="Root for unlock logs.")
    parser.add_argument(
        "--shared-workdir",
        action="store_true",
        help="Unlock the shared checkout-level Snakemake workdir instead of isolated row workdirs.",
    )
    add_snakemake_workdir_arguments(parser)
    parser.add_argument("--dry-run", action="store_true", help="Print unlock command(s) without running them.")
    parser.add_argument("--print-command", action="store_true", help="Print each resolved unlock command.")
    parser.add_argument(
        "--snakemake-arg",
        action="append",
        default=[],
        help="Extra raw argument passed to Snakemake. May be repeated.",
    )


def add_batch_prepare_envs_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for serialized Snakemake rule-environment creation."""
    add_common_config_arguments(parser)
    add_runner_arguments(parser)
    parser.add_argument(
        "--batch-table",
        default=DEFAULT_BATCH_TABLE,
        help=f"Canonical batch TSV. Defaults to {DEFAULT_BATCH_TABLE}.",
    )
    parser.add_argument(
        "--trial-id",
        action="append",
        default=[],
        help="Trial ID whose rule environments should be prepared. May be repeated. Defaults to all rows.",
    )
    parser.add_argument("--profile", help="Snakemake profile path. Defaults to profiles/local.")
    parser.add_argument("--target", default="all", help="Snakemake target rule or file.")
    parser.add_argument("--log-root", default=DEFAULT_LOG_ROOT, help="Root for environment-prep logs.")
    parser.add_argument(
        "--shared-workdir",
        action="store_true",
        help="Prepare environments from the shared checkout-level Snakemake workdir.",
    )
    add_snakemake_workdir_arguments(parser)
    parser.add_argument("--dry-run", action="store_true", help="Print preparation command(s) without running them.")
    parser.add_argument("--print-command", action="store_true", help="Print each resolved preparation command.")
    parser.add_argument(
        "--snakemake-arg",
        action="append",
        default=[],
        help="Extra raw argument passed to Snakemake. May be repeated.",
    )


def add_snakemake_workdir_arguments(parser: argparse.ArgumentParser) -> None:
    """Add Snakemake workdir and rule-env cache options."""
    parser.add_argument(
        "--workdir-root",
        default=DEFAULT_SNAKEMAKE_WORKDIR_ROOT,
        help=(
            "Root for isolated per-trial Snakemake workdirs. "
            f"Defaults to {DEFAULT_SNAKEMAKE_WORKDIR_ROOT}."
        ),
    )
    parser.add_argument(
        "--snakemake-conda-prefix",
        default=DEFAULT_SNAKEMAKE_CONDA_PREFIX,
        help=(
            "Shared Snakemake rule-environment prefix used with isolated workdirs. "
            f"Defaults to {DEFAULT_SNAKEMAKE_CONDA_PREFIX}."
        ),
    )


def add_reference_scope_arguments(parser: argparse.ArgumentParser) -> None:
    """Add options selecting which configured BWA references are considered."""
    parser.add_argument(
        "--all-configured-hosts",
        action="store_true",
        help=(
            "Include both human_ref and mouse_ref instead of only a config-selected "
            "host reference. Use this when preparing indexes for mixed-host batch tables."
        ),
    )


def add_references_check_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for validating configured reference indexes."""
    add_processing_config_stack_arguments(parser)
    add_reference_scope_arguments(parser)


def add_references_prepare_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for preparing configured BWA reference indexes."""
    add_processing_config_stack_arguments(parser)
    add_reference_scope_arguments(parser)
    parser.add_argument("--log-root", default=DEFAULT_REFERENCE_LOG_ROOT, help="Root for reference preparation logs.")
    parser.add_argument("--log-dir", help="Exact reference preparation log directory. Defaults under --log-root.")
    parser.add_argument("--dry-run", action="store_true", help="Print work to do without running BWA or writing manifests.")
    parser.add_argument("--print-command", action="store_true", help="Print each resolved BWA/index command.")
    parser.add_argument("--force", action="store_true", help="Rebuild indexes even when all sidecars already exist.")
    parser.add_argument(
        "--env-mode",
        choices=ANALYSIS_ENV_MODES,
        default="managed",
        help="Create/use a project-local bio-tools env, use a named env, use an explicit prefix, or run directly.",
    )
    parser.add_argument(
        "--manager",
        default="auto",
        help="Conda-compatible manager for --env-mode managed, named, or prefix. Defaults to auto.",
    )
    parser.add_argument(
        "--analysis-env-root",
        default=DEFAULT_ANALYSIS_ENV_ROOT,
        help=(
            "Root for project-managed conda prefixes when --env-mode managed is used. "
            f"Defaults to {DEFAULT_ANALYSIS_ENV_ROOT}."
        ),
    )
    parser.add_argument(
        "--bio-env-prefix",
        default=os.environ.get(BIO_TOOLS_ENV_PREFIX_ENVVAR),
        help=(
            "Conda environment prefix for BWA/bioinformatics tools when --env-mode "
            f"prefix is used. Defaults to ${BIO_TOOLS_ENV_PREFIX_ENVVAR}."
        ),
    )


def add_analysis_config_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--analysis-config",
        default=DEFAULT_ANALYSIS_CONFIG,
        help=(
            "Analysis YAML config. "
            f"Defaults to {DEFAULT_ANALYSIS_CONFIG}."
        ),
    )


def add_analysis_init_arguments(parser: argparse.ArgumentParser) -> None:
    add_analysis_config_argument(parser)
    parser.add_argument(
        "--template",
        default=TEMPLATE_ANALYSIS_CONFIG,
        help=f"Template to copy. Defaults to {TEMPLATE_ANALYSIS_CONFIG}.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing analysis config.",
    )


def add_analysis_validate_arguments(parser: argparse.ArgumentParser) -> None:
    add_analysis_config_argument(parser)
    parser.add_argument(
        "steps",
        nargs="*",
        choices=ANALYSIS_STEP_CHOICES,
        help="Optional analysis steps whose config sections should be checked.",
    )


def add_analysis_run_arguments(parser: argparse.ArgumentParser) -> None:
    add_analysis_config_argument(parser)
    parser.add_argument(
        "steps",
        nargs="+",
        choices=ANALYSIS_STEP_CHOICES,
        help="Analysis step(s) to run in the order provided.",
    )
    parser.add_argument("--log-root", default=DEFAULT_ANALYSIS_LOG_ROOT, help="Root for analysis logs.")
    parser.add_argument("--log-dir", help="Exact analysis log directory. Defaults under --log-root.")
    parser.add_argument("--dry-run", action="store_true", help="Validate and print commands without running scripts.")
    parser.add_argument("--print-command", action="store_true", help="Print each resolved analysis command.")
    parser.add_argument(
        "--env-mode",
        choices=ANALYSIS_ENV_MODES,
        default="managed",
        help=(
            "Create/use project-local envs from workflow/envs/*.yaml, use named "
            "envs, use explicit prefixes, or run directly in the current environment."
        ),
    )
    parser.add_argument(
        "--manager",
        default="auto",
        help="Conda-compatible manager for --env-mode managed, named, or prefix. Defaults to auto.",
    )
    add_analysis_env_prefix_arguments(parser)


def add_meta_config_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--analysis-config",
        default=DEFAULT_META_CONFIG,
        help=(
            "Meta-analysis YAML config. "
            f"Defaults to {DEFAULT_META_CONFIG}."
        ),
    )


def add_meta_init_arguments(parser: argparse.ArgumentParser) -> None:
    add_meta_config_argument(parser)
    parser.add_argument(
        "--template",
        default=TEMPLATE_META_CONFIG,
        help=f"Template to copy. Defaults to {TEMPLATE_META_CONFIG}.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing meta-analysis config.",
    )


def add_meta_validate_arguments(parser: argparse.ArgumentParser) -> None:
    add_meta_config_argument(parser)
    parser.add_argument(
        "steps",
        nargs="*",
        choices=META_STEP_CHOICES,
        help="Optional meta-analysis steps whose config sections should be checked.",
    )


def add_meta_step_arguments(parser: argparse.ArgumentParser) -> None:
    add_meta_config_argument(parser)
    parser.add_argument("--log-root", default=DEFAULT_META_LOG_ROOT, help="Root for meta-analysis logs.")
    parser.add_argument("--log-dir", help="Exact meta-analysis log directory. Defaults under --log-root.")
    parser.add_argument("--dry-run", action="store_true", help="Validate and print the command without running it.")
    parser.add_argument("--print-command", action="store_true", help="Print the resolved meta-analysis command.")
    parser.add_argument(
        "--env-mode",
        choices=ANALYSIS_ENV_MODES,
        default="managed",
        help=(
            "Create/use project-local envs from workflow/envs/*.yaml, use named "
            "envs, use explicit prefixes, or run directly in the current environment."
        ),
    )
    parser.add_argument(
        "--manager",
        default="auto",
        help="Conda-compatible manager for --env-mode managed, named, or prefix. Defaults to auto.",
    )
    add_analysis_env_prefix_arguments(parser)


def add_analysis_env_prefix_arguments(parser: argparse.ArgumentParser) -> None:
    """Add explicit analysis environment prefix options."""
    parser.add_argument(
        "--analysis-env-root",
        default=os.environ.get("LOW_BM_ANALYSIS_ENV_ROOT", DEFAULT_ANALYSIS_ENV_ROOT),
        help=(
            "Root for project-managed analysis conda prefixes when --env-mode "
            f"managed is used. Defaults to {DEFAULT_ANALYSIS_ENV_ROOT}."
        ),
    )
    parser.add_argument(
        "--r-env-prefix",
        default=os.environ.get(R_TOOLS_ENV_PREFIX_ENVVAR),
        help=(
            "Conda environment prefix for R-heavy analysis steps when "
            f"--env-mode prefix is used. Defaults to ${R_TOOLS_ENV_PREFIX_ENVVAR}."
        ),
    )
    parser.add_argument(
        "--bio-env-prefix",
        default=os.environ.get(BIO_TOOLS_ENV_PREFIX_ENVVAR),
        help=(
            "Conda environment prefix for BLAST/bioinformatics analysis steps "
            f"when --env-mode prefix is used. Defaults to ${BIO_TOOLS_ENV_PREFIX_ENVVAR}."
        ),
    )
    parser.add_argument(
        "--microclean-env-prefix",
        default=os.environ.get(MICROCLEAN_ENV_PREFIX_ENVVAR),
        help=(
            "Conda environment prefix for micRoclean meta steps when "
            f"--env-mode prefix is used. Defaults to ${MICROCLEAN_ENV_PREFIX_ENVVAR}."
        ),
    )


def add_master_job_arguments(parser: argparse.ArgumentParser) -> None:
    """Add options that affect only the outer SLURM master job."""
    parser.add_argument(
        "--activate-command",
        default=os.environ.get("LOW_BM_ACTIVATE_COMMAND", ""),
        help=(
            "Advanced fallback shell command evaluated in the master job before "
            "Snakemake. Normal runs should use `low-bm setup runner` instead."
        ),
    )
    parser.add_argument("--master-cpus", default="4", help="CPUs for the SLURM master job.")
    parser.add_argument("--master-mem", default="8G", help="Memory for the SLURM master job.")
    parser.add_argument(
        "--master-time",
        default="1-00:00:00",
        help="Walltime for the SLURM master job.",
    )
    parser.add_argument("--master-partition", help="Optional partition for the master job.")
    parser.add_argument(
        "--master-extra-sbatch",
        action="append",
        default=[],
        help="Extra raw #SBATCH directive body, such as '--mail-type=END,FAIL'.",
    )


def add_runner_arguments(parser: argparse.ArgumentParser) -> None:
    """Add V1 host-runner options shared by run entry points."""
    parser.add_argument(
        "--runner",
        choices=["host"],
        default="host",
        help="Runner implementation. V1 supports host; container support will plug in here later.",
    )
    parser.add_argument(
        "--runner-prefix",
        default=DEFAULT_RUNNER_PREFIX,
        help=f"Project-local runner environment prefix. Defaults to {DEFAULT_RUNNER_PREFIX}.",
    )
    parser.add_argument(
        "--manager",
        choices=("auto",) + ENV_MANAGERS,
        default="auto",
        help="Conda-compatible manager used to run the runner env.",
    )
    parser.add_argument(
        "--conda-base-prefix",
        default=os.environ.get("LOW_BM_CONDA_BASE_PREFIX"),
        help=(
            "Real conda base prefix used by Snakemake to activate rule environments. "
            "Defaults to the value detected during `low-bm setup runner`."
        ),
    )


def add_setup_runner_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for creating the project-local runner env."""
    parser.add_argument(
        "--runner-prefix",
        default=DEFAULT_RUNNER_PREFIX,
        help=f"Runner environment prefix to create. Defaults to {DEFAULT_RUNNER_PREFIX}.",
    )
    parser.add_argument(
        "--manager",
        choices=("auto",) + ENV_MANAGERS,
        default="auto",
        help="Conda-compatible manager used to create the runner env.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove and recreate an existing runner environment prefix.",
    )
    parser.add_argument(
        "--conda-base-prefix",
        default=os.environ.get("LOW_BM_CONDA_BASE_PREFIX"),
        help=(
            "Real conda base prefix to store for Snakemake rule-env activation. "
            "Defaults to `<manager> info --base`."
        ),
    )


def add_doctor_runner_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for checking the project-local runner env."""
    parser.add_argument(
        "--mode",
        choices=["local", "slurm"],
        default="local",
        help="Runtime mode to validate.",
    )
    parser.add_argument(
        "--runner",
        choices=["host"],
        default="host",
        help="Runner implementation to validate. V1 supports host.",
    )
    parser.add_argument(
        "--runner-prefix",
        default=DEFAULT_RUNNER_PREFIX,
        help=f"Runner environment prefix to validate. Defaults to {DEFAULT_RUNNER_PREFIX}.",
    )
    parser.add_argument(
        "--manager",
        choices=("auto",) + ENV_MANAGERS,
        default="auto",
        help="Conda-compatible manager used to run the runner env.",
    )
    parser.add_argument(
        "--conda-base-prefix",
        default=os.environ.get("LOW_BM_CONDA_BASE_PREFIX"),
        help=(
            "Real conda base prefix used by Snakemake to activate rule environments. "
            "Defaults to the value recorded by `low-bm setup runner`."
        ),
    )
    parser.add_argument(
        "--rule-env-smoke-test",
        action="store_true",
        help=(
            "Run a tiny Snakemake conda-rule smoke test and confirm the rule Python "
            "comes from the Snakemake rule-env cache, not the runner env."
        ),
    )


def run_command(args: argparse.Namespace) -> int:
    """Run or submit one batch after resolving its per-row config file."""
    row_config = resolve_row_config(args)
    spec = build_run_spec(args, row_config)
    return execute_run_spec(spec, args)


def batch_submit_command(args: argparse.Namespace) -> int:
    """Submit one independent run for every row in the batch table."""
    batch_table = Path(args.batch_table)
    ensure_batch_table_exists(batch_table)
    rows = read_batch_table(batch_table)
    if not rows:
        print(f"No rows found in {args.batch_table}.")
        return 0
    validate_unique_batch_rows(rows)

    exit_code = 0
    for row in rows:
        row_config = write_run_config(row, args.run_config_dir)
        if should_use_batch_isolated_workdirs(args):
            require_shared_reference_indexes(
                resolve_configfiles(args.configfile, args.extra_configfile),
                row_config,
            )
        # Reuse the same run-building path as `low-bm run` so batch-table and
        # direct one-batch launches cannot drift into different behaviors.
        row_args = argparse.Namespace(**vars(args))
        row_args.row_config = str(row_config)
        row_args.trial_id = None
        row_args.trial_descript = None
        row_args.exp_dir = None
        row_args.metadata = None
        row_args.host = None
        row_args.process_umis = None
        row_args.unlock = False
        row_args.isolated_workdir = should_use_batch_isolated_workdirs(args)
        row_args.job_name = None
        spec = build_run_spec(row_args, row_config)
        print(f"[{row.trialID}] prepared {row_config}")
        if args.dry_run:
            # In batch dry-run mode we intentionally avoid writing the master
            # script or invoking Snakemake; this reports the command shape only.
            command = build_execution_command(spec, row_args)
            if args.mode == "slurm":
                script = master_script_path(spec)
                print(f"[{row.trialID}] would submit master job script: {script}")
                print(f"[{row.trialID}] {shlex.join(command)}")
            else:
                print(f"[{row.trialID}] would run: {shlex.join(command)}")
            continue
        result = execute_run_spec(spec, row_args)
        if result != 0:
            exit_code = result
            if args.mode == "local":
                break
    return exit_code


def batch_unlock_command(args: argparse.Namespace) -> int:
    """Unlock one or more batch-row Snakemake workdirs locally."""
    if bool(args.all) == bool(args.trial_id):
        raise SystemExit("Choose exactly one of --trial-id or --all for batch unlock.")

    batch_table = Path(args.batch_table)
    ensure_batch_table_exists(batch_table)
    rows = read_batch_table(batch_table)
    validate_unique_batch_rows(rows)
    selected_rows = select_batch_unlock_rows(rows, args.trial_id, args.all)

    exit_code = 0
    for row in selected_rows:
        row_config = write_run_config(row, args.run_config_dir)
        row_args = argparse.Namespace(**vars(args))
        row_args.row_config = str(row_config)
        row_args.trial_id = None
        row_args.trial_descript = None
        row_args.exp_dir = None
        row_args.metadata = None
        row_args.host = None
        row_args.process_umis = None
        row_args.job_name = None
        row_args.unlock = True
        row_args.dry_run = False
        row_args.isolated_workdir = not args.shared_workdir
        row_args.print_command = args.print_command or args.dry_run
        spec = build_run_spec(row_args, row_config)
        if args.dry_run:
            command = build_execution_command(spec, row_args)
            print(f"[{row.trialID}] would unlock: {shlex.join(command)}")
            continue
        result = execute_run_spec(spec, row_args)
        if result != 0:
            exit_code = result
        else:
            target = spec.workdir if spec.workdir else REPO_ROOT
            print(f"[{row.trialID}] unlocked Snakemake workdir: {target}")
    return exit_code


def batch_prepare_envs_command(args: argparse.Namespace) -> int:
    """Prepare Snakemake rule conda environments sequentially for batch rows."""
    batch_table = Path(args.batch_table)
    ensure_batch_table_exists(batch_table)
    rows = read_batch_table(batch_table)
    if not rows:
        print(f"No rows found in {args.batch_table}.")
        return 0
    validate_unique_batch_rows(rows)
    selected_rows = select_batch_rows_by_trial_ids(rows, args.trial_id)

    exit_code = 0
    for row in selected_rows:
        row_config = write_run_config(row, args.run_config_dir)
        row_args = argparse.Namespace(**vars(args))
        row_args.row_config = str(row_config)
        row_args.trial_id = None
        row_args.trial_descript = None
        row_args.exp_dir = None
        row_args.metadata = None
        row_args.host = None
        row_args.process_umis = None
        row_args.job_name = f"low-bm-prepare-envs-{row.trialID}"
        row_args.mode = "local"
        row_args.profile = args.profile or str(REPO_ROOT / "profiles" / "local")
        row_args.unlock = False
        row_args.dry_run = False
        row_args.isolated_workdir = not args.shared_workdir
        row_args.print_command = args.print_command or args.dry_run
        row_args.snakemake_arg = list(args.snakemake_arg or []) + ["--conda-create-envs-only"]
        spec = build_run_spec(row_args, row_config)
        command = build_execution_command(spec, row_args)
        if args.dry_run:
            print(f"[{row.trialID}] would prepare envs: {shlex.join(command)}")
            continue
        result = execute_run_spec(spec, row_args)
        if result != 0:
            exit_code = result
            break
        print(f"[{row.trialID}] prepared Snakemake rule environments.")
    return exit_code


def references_check_command(args: argparse.Namespace) -> int:
    """Validate configured shared BWA indexes without modifying them."""
    refs = resolve_bwa_references_for_args(args)
    warnings = bwa_reference_manifest_warnings(refs)
    for warning in warnings:
        print(f"Warning: {warning}", file=sys.stderr)

    failures = missing_bwa_reference_inputs(refs)
    if failures:
        raise SystemExit(render_bwa_reference_failure(failures))

    for ref in refs:
        print(f"[ok] {ref.key}: {ref.path}")
    print(f"BWA reference check passed for {len(refs)} reference(s).")
    return 0


def references_prepare_bwa_indexes_command(args: argparse.Namespace) -> int:
    """Prepare configured shared BWA indexes outside the processing DAG."""
    refs = resolve_bwa_references_for_args(args)
    missing_fastas = [f"{ref.key}: missing FASTA {ref.path}" for ref in refs if not ref.path.exists()]
    if missing_fastas:
        raise SystemExit("Cannot prepare BWA indexes for missing reference FASTA(s):\n" + "\n".join(f"- {item}" for item in missing_fastas))

    created_utc = datetime.now(timezone.utc)
    log_dir = resolve_reference_log_dir(args, created_utc)
    if not args.dry_run:
        log_dir.mkdir(parents=True, exist_ok=True)

    env_mode, manager, env_prefix = resolve_reference_bio_env(args, validate_exists=not args.dry_run)
    if not args.dry_run:
        env_result = ensure_reference_bio_env(env_mode, manager, env_prefix, log_dir)
        if env_result != 0:
            return env_result

    exit_code = 0
    bwa_version: str | None = None
    for ref in refs:
        missing_sidecars = missing_bwa_index_sidecars(ref.path)
        manifest_status, manifest_detail = bwa_index_manifest_status(ref)
        needs_rebuild = bool(args.force or missing_sidecars or manifest_status == "stale_checksum")
        command = build_reference_bwa_command(["bwa", "index", str(ref.path)], env_mode, manager, env_prefix)

        if needs_rebuild:
            if args.print_command or args.dry_run:
                print(f"[{ref.key}] {shlex.join(command)}")
            if args.dry_run:
                continue

            log_file = reference_step_log_file(log_dir, ref)
            print(f"[run] {ref.key}: indexing {ref.path}")
            with log_file.open("a") as log:
                log.write(f"\n# {datetime.now(timezone.utc).isoformat()}\n")
                log.write(f"$ {shlex.join(command)}\n")
                log.flush()
                completed = subprocess.run(
                    command,
                    cwd=REPO_ROOT,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            if completed.returncode != 0:
                print(f"BWA indexing failed for {ref.key} with exit code {completed.returncode}. See {log_file}.", file=sys.stderr)
                exit_code = completed.returncode
                continue

            missing_after = missing_bwa_index_sidecars(ref.path)
            if missing_after:
                print(
                    f"BWA indexing completed for {ref.key}, but required sidecar(s) are still missing: "
                    + ", ".join(str(path) for path in missing_after),
                    file=sys.stderr,
                )
                exit_code = 1
                continue

            if bwa_version is None:
                bwa_version = detect_bwa_version(env_mode, manager, env_prefix)
            write_bwa_index_manifest(ref, bwa_version)
            print(f"[ok] {ref.key}: indexed {ref.path}")
            continue

        if manifest_status == "fresh":
            print(f"[skip] {ref.key}: complete BWA index already exists")
            continue

        if args.dry_run:
            print(f"[{ref.key}] would write manifest: {bwa_index_manifest_path(ref.path)} ({manifest_detail})")
            continue

        if bwa_version is None:
            bwa_version = detect_bwa_version(env_mode, manager, env_prefix)
        write_bwa_index_manifest(ref, bwa_version)
        print(f"[manifest] {ref.key}: wrote {bwa_index_manifest_path(ref.path)}")

    return exit_code


def resolve_bwa_references_for_args(args: argparse.Namespace) -> list[BwaReferenceSpec]:
    configfiles = resolve_configfiles(args.configfile, args.extra_configfile)
    validate_processing_config_stack(configfiles)
    values = merge_top_level_yaml_scalar_values(configfiles)
    return resolve_bwa_references(values, all_configured_hosts=bool(getattr(args, "all_configured_hosts", False)))


def resolve_bwa_references(values: dict[str, str], all_configured_hosts: bool = False) -> list[BwaReferenceSpec]:
    refs: list[BwaReferenceSpec] = []
    missing_keys: list[str] = []
    seen_paths: set[Path] = set()
    for key in selected_bwa_reference_keys(values, all_configured_hosts):
        value = clean_yaml_scalar(values.get(key))
        if not value:
            missing_keys.append(key)
            continue
        path = repo_path(Path(value)).resolve()
        if path in seen_paths:
            continue
        refs.append(BwaReferenceSpec(key=key, path=path))
        seen_paths.add(path)
    if missing_keys:
        raise SystemExit("Processing config is missing BWA reference value(s): " + ", ".join(missing_keys))
    return refs


def selected_bwa_reference_keys(values: dict[str, str], all_configured_hosts: bool = False) -> list[str]:
    if all_configured_hosts:
        host_keys = ["human_ref", "mouse_ref"]
    else:
        host = clean_yaml_scalar(values.get("host", "")).lower()
        if host == "human":
            host_keys = ["human_ref"]
        elif host == "mouse":
            host_keys = ["mouse_ref"]
        else:
            raise SystemExit(
                "Config key 'host' must be either 'human' or 'mouse' for BWA reference resolution. "
                "Use --all-configured-hosts when validating or preparing both host references."
            )
    return [*host_keys, "viral_ref", "bact16s_ref"]


def bwa_index_sidecars(reference: Path) -> list[Path]:
    return [Path(f"{reference}{extension}") for extension in BWA_INDEX_EXTENSIONS]


def missing_bwa_index_sidecars(reference: Path) -> list[Path]:
    return [path for path in bwa_index_sidecars(reference) if not path.exists()]


def missing_bwa_reference_inputs(refs: list[BwaReferenceSpec]) -> list[str]:
    failures: list[str] = []
    for ref in refs:
        if not ref.path.exists():
            failures.append(f"{ref.key}: missing FASTA {ref.path}")
            continue
        for sidecar in missing_bwa_index_sidecars(ref.path):
            failures.append(f"{ref.key}: missing BWA sidecar {sidecar}")
    return failures


def render_bwa_reference_failure(failures: list[str]) -> str:
    rendered = "\n".join(f"- {failure}" for failure in failures)
    return (
        "Missing BWA reference index files:\n"
        f"{rendered}\n"
        "Run `low-bm references prepare-bwa-indexes` before processing batches that use these references. "
        "Add `--all-configured-hosts` when preparing indexes for mixed-host batch tables."
    )


def bwa_index_manifest_path(reference: Path) -> Path:
    return Path(f"{reference}.low-bm-bwa-index.json")


def bwa_index_manifest_status(ref: BwaReferenceSpec) -> tuple[str, str]:
    manifest = bwa_index_manifest_path(ref.path)
    if not manifest.exists():
        return "missing", f"missing manifest {manifest}"
    try:
        values = json.loads(manifest.read_text())
    except json.JSONDecodeError:
        return "invalid", f"invalid manifest JSON {manifest}"
    if values.get("schema_version") != BWA_INDEX_MANIFEST_SCHEMA_VERSION:
        return "invalid", f"unsupported manifest schema in {manifest}"
    if values.get("reference_sha256") != hash_file(ref.path):
        return "stale_checksum", f"manifest checksum does not match {ref.path}"
    expected_sidecars = [str(path) for path in bwa_index_sidecars(ref.path)]
    if values.get("sidecars") != expected_sidecars:
        return "invalid", f"manifest sidecar list is stale for {ref.path}"
    if not values.get("bwa_version"):
        return "invalid", f"manifest is missing bwa_version for {ref.path}"
    return "fresh", str(manifest)


def bwa_reference_manifest_warnings(refs: list[BwaReferenceSpec]) -> list[str]:
    warnings: list[str] = []
    for ref in refs:
        if not ref.path.exists() or missing_bwa_index_sidecars(ref.path):
            continue
        status, detail = bwa_index_manifest_status(ref)
        if status != "fresh":
            warnings.append(f"{ref.key}: {detail}")
    return warnings


def write_bwa_index_manifest(ref: BwaReferenceSpec, bwa_version: str) -> None:
    manifest = {
        "schema_version": BWA_INDEX_MANIFEST_SCHEMA_VERSION,
        "reference_path": str(ref.path),
        "reference_sha256": hash_file(ref.path),
        "bwa_version": bwa_version,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sidecars": [str(path) for path in bwa_index_sidecars(ref.path)],
    }
    bwa_index_manifest_path(ref.path).write_text(json.dumps(manifest, indent=2) + "\n")


def resolve_reference_log_dir(args: argparse.Namespace, created_utc: datetime) -> Path:
    if getattr(args, "log_dir", None):
        return repo_path(Path(args.log_dir))
    timestamp = created_utc.strftime("%Y%m%dT%H%M%SZ")
    return repo_path(Path(getattr(args, "log_root", DEFAULT_REFERENCE_LOG_ROOT))) / f"{timestamp}_bwa-indexes"


def reference_step_log_file(log_dir: Path, ref: BwaReferenceSpec) -> Path:
    path_hash = hashlib.sha256(str(ref.path).encode()).hexdigest()[:8]
    return log_dir / f"bwa-index.{sanitize_path_slug(ref.key)}.{path_hash}.log"


def resolve_reference_bio_env(
    args: argparse.Namespace,
    validate_exists: bool,
) -> tuple[str, EnvManagerSpec | None, Path | None]:
    env_mode = getattr(args, "env_mode", "managed")
    if env_mode not in ANALYSIS_ENV_MODES:
        raise SystemExit(f"Unsupported reference env mode: {env_mode}")
    manager = resolve_env_manager(getattr(args, "manager", "auto")) if env_mode in {"managed", "named", "prefix"} else None
    if env_mode == "managed":
        env_root = repo_path(Path(getattr(args, "analysis_env_root", DEFAULT_ANALYSIS_ENV_ROOT))).resolve()
        return env_mode, manager, managed_analysis_env_prefix(BIO_TOOLS_ENV, env_root)
    if env_mode == "prefix":
        value = getattr(args, "bio_env_prefix", None)
        if not value:
            raise SystemExit(f"--env-mode prefix requires --bio-env-prefix or ${BIO_TOOLS_ENV_PREFIX_ENVVAR}.")
        prefix = repo_path(Path(value)).resolve()
        if validate_exists and not prefix.is_dir():
            raise SystemExit(f"Missing bio-tools environment prefix: {prefix}")
        return env_mode, manager, prefix
    return env_mode, manager, None


def ensure_reference_bio_env(
    env_mode: str,
    manager: EnvManagerSpec | None,
    env_prefix: Path | None,
    log_dir: Path,
) -> int:
    if env_mode != "managed" or env_prefix is None:
        return 0
    if env_prefix.is_dir():
        return 0
    if manager is None:
        raise SystemExit("--env-mode managed requires a conda-compatible manager.")
    env_prefix.parent.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "reference_env_setup.log"
    command = build_conda_env_create_command(env_prefix, repo_path(BIO_TOOLS_ENV), manager)
    print(f"Creating reference bio-tools environment at {env_prefix}.")
    with log_file.open("a") as log:
        log.write(f"\n# {datetime.now(timezone.utc).isoformat()}\n")
        log.write(f"$ {shlex.join(command)}\n")
        log.flush()
        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        print(f"Reference bio-tools environment setup failed with exit code {completed.returncode}. See {log_file}.", file=sys.stderr)
        return completed.returncode
    write_analysis_env_metadata(env_prefix, BIO_TOOLS_ENV, manager)
    return 0


def build_reference_bwa_command(
    command: list[str],
    env_mode: str,
    manager: EnvManagerSpec | None,
    env_prefix: Path | None,
) -> list[str]:
    return wrap_analysis_command(
        command,
        analysis_env_name(BIO_TOOLS_ENV),
        env_prefix,
        env_mode,
        manager,
    )


def detect_bwa_version(env_mode: str, manager: EnvManagerSpec | None, env_prefix: Path | None) -> str:
    command = build_reference_bwa_command(["bwa"], env_mode, manager, env_prefix)
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        if stripped.lower().startswith("version:"):
            return stripped.split(":", 1)[1].strip()
    return "unknown"


def analysis_init_command(args: argparse.Namespace) -> int:
    """Create an editable local analysis config from the tracked template."""
    destination_arg = Path(args.analysis_config)
    destination = repo_path(destination_arg)
    template_arg = Path(args.template)
    template = repo_path(template_arg)

    if not template.exists():
        raise SystemExit(f"Missing analysis template: {template_arg}")
    if destination.exists() and not args.force:
        raise SystemExit(
            f"Analysis config already exists: {destination_arg}. "
            "Use --force to overwrite it."
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(template, destination)
    print(f"Wrote analysis config: {destination_arg}")
    return 0


def analysis_validate_command(args: argparse.Namespace) -> int:
    """Validate analysis config shape without running analysis scripts."""
    analysis_config = resolve_analysis_config(args.analysis_config)
    expanded_steps = expand_analysis_steps(list(args.steps or []))
    validate_analysis_config(analysis_config, expanded_steps)
    if expanded_steps:
        print("Analysis config is valid for step(s): " + ", ".join(expanded_steps))
    else:
        print(f"Analysis config is valid: {analysis_config}")
    return 0


def analysis_run_command(args: argparse.Namespace) -> int:
    """Run selected post-processing analyses in explicit order."""
    spec = build_analysis_run_spec(args)
    spec.log_dir.mkdir(parents=True, exist_ok=True)
    write_analysis_provenance(spec)

    if args.print_command or spec.dry_run:
        for command_spec in spec.commands:
            print(shlex.join(command_spec.execution_command))

    if spec.dry_run:
        print(f"Analysis dry-run complete. See {spec.log_dir}.")
        return 0

    env_returncode = ensure_managed_analysis_envs(spec, "Analysis")
    if env_returncode != 0:
        return env_returncode

    for index, command_spec in enumerate(spec.commands, start=1):
        returncode = run_analysis_step(command_spec, spec.log_dir, index)
        if returncode != 0:
            print(
                f"Analysis step '{command_spec.step}' exited with {returncode}. "
                f"See {analysis_step_log_file(spec.log_dir, index, command_spec.step)}.",
                file=sys.stderr,
            )
            return returncode
    print(f"Analysis complete. See {spec.log_dir}.")
    return 0


def meta_init_command(args: argparse.Namespace) -> int:
    """Create an editable local meta-analysis config from the tracked template."""
    destination_arg = Path(args.analysis_config)
    destination = repo_path(destination_arg)
    template_arg = Path(args.template)
    template = repo_path(template_arg)

    if not template.exists():
        raise SystemExit(f"Missing meta-analysis template: {template_arg}")
    if destination.exists() and not args.force:
        raise SystemExit(
            f"Meta-analysis config already exists: {destination_arg}. "
            "Use --force to overwrite it."
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(template, destination)
    print(f"Wrote meta-analysis config: {destination_arg}")
    return 0


def meta_validate_command(args: argparse.Namespace) -> int:
    """Validate meta-analysis config shape without running scripts."""
    meta_config = resolve_meta_config(args.analysis_config)
    steps = list(args.steps or [])
    validate_meta_config(meta_config, steps)
    if steps:
        print("Meta-analysis config is valid for step(s): " + ", ".join(steps))
    else:
        print(f"Meta-analysis config is valid: {meta_config}")
    return 0


def meta_compile_phyloseq_command(args: argparse.Namespace) -> int:
    """Run multi-batch phyloseq and ASV FASTA compilation."""
    return execute_meta_step(args, "compile-phyloseq")


def meta_decontaminate_phyloseq_command(args: argparse.Namespace) -> int:
    """Run optional micRoclean decontamination on one phyloseq endpoint."""
    return execute_meta_step(args, "decontaminate-phyloseq")


def meta_differential_abundance_command(args: argparse.Namespace) -> int:
    """Run differential-abundance meta-analysis."""
    return execute_meta_step(args, "differential-abundance")


def execute_meta_step(args: argparse.Namespace, step: str) -> int:
    """Run one meta-analysis step with provenance and logging."""
    spec = build_meta_run_spec(args, step)
    spec.log_dir.mkdir(parents=True, exist_ok=True)
    write_meta_provenance(spec)

    command_spec = spec.commands[0]
    if args.print_command or spec.dry_run:
        print(shlex.join(command_spec.execution_command))

    if spec.dry_run:
        print(f"Meta-analysis dry-run complete. See {spec.log_dir}.")
        return 0

    env_returncode = ensure_managed_analysis_envs(spec, "Meta-analysis")
    if env_returncode != 0:
        return env_returncode

    returncode = run_analysis_step(command_spec, spec.log_dir, 1)
    if returncode != 0:
        print(
            f"Meta-analysis step '{command_spec.step}' exited with {returncode}. "
            f"See {analysis_step_log_file(spec.log_dir, 1, command_spec.step)}.",
            file=sys.stderr,
        )
        return returncode
    print(f"Meta-analysis complete. See {spec.log_dir}.")
    return 0


def setup_runner_command(args: argparse.Namespace) -> int:
    """Create or validate the project-local Snakemake runner environment."""
    prefix_arg = Path(args.runner_prefix)
    prefix = repo_path(prefix_arg)
    env_file = REPO_ROOT / DEFAULT_RUNNER_ENV_FILE
    if not env_file.exists():
        raise SystemExit(f"Missing runner environment file: {DEFAULT_RUNNER_ENV_FILE}")

    metadata_path = runner_metadata_path(prefix_arg)
    manager = resolve_env_manager(args.manager, metadata_path=metadata_path)
    conda_base_prefix = resolve_conda_base_prefix(
        getattr(args, "conda_base_prefix", None),
        metadata_path,
        manager,
    )
    if path_exists_for_repo_run(prefix):
        write_runner_metadata(prefix_arg, env_file, manager, conda_base_prefix, preserve_created=True)
        if not args.force:
            resolved_prefix = prefix.resolve()
            runner = HostRunnerSpec(
                prefix=resolved_prefix,
                manager=manager,
                metadata_path=metadata_path,
                bin_dir=runner_bin_dir(resolved_prefix),
                conda_base_prefix=conda_base_prefix,
            )
            print(f"Runner env already exists: {prefix}")
            return check_runner_command(runner, ["snakemake", "--version"], "snakemake")
        remove_runner_prefix(prefix)

    command = build_runner_create_command(prefix_arg, env_file, manager)
    print(f"Creating runner env at {prefix_arg} with {manager.name}.")
    completed = subprocess.run(command, cwd=REPO_ROOT, check=False)
    if completed.returncode != 0:
        return completed.returncode

    write_runner_metadata(prefix_arg, env_file, manager, conda_base_prefix, preserve_created=False)
    print(f"Wrote runner metadata: {metadata_path}")
    resolved_prefix = prefix.resolve()
    return check_runner_command(
        HostRunnerSpec(
            prefix=resolved_prefix,
            manager=manager,
            metadata_path=metadata_path,
            bin_dir=runner_bin_dir(resolved_prefix),
            conda_base_prefix=conda_base_prefix,
        ),
        ["snakemake", "--version"],
        "snakemake",
    )


def doctor_runner_command(args: argparse.Namespace) -> int:
    """Check whether the host runner is ready for the selected execution mode."""
    runner = resolve_host_runner(args)
    checks: list[tuple[str, bool, str]] = []

    snakemake = run_runner_command(runner, ["snakemake", "--version"])
    checks.append(("snakemake", snakemake.returncode == 0, runner_check_detail(snakemake)))

    conda_base_errors = conda_base_activation_errors(runner.conda_base_prefix)
    checks.append(
        (
            "conda base activation",
            not conda_base_errors,
            str(runner.conda_base_prefix) if not conda_base_errors else "; ".join(conda_base_errors),
        )
    )
    if getattr(args, "rule_env_smoke_test", False):
        smoke_ok, smoke_detail = run_rule_env_smoke_test(runner)
        checks.append(("rule env python smoke test", smoke_ok, smoke_detail))

    if args.mode == "slurm":
        plugin = run_runner_command(
            runner,
            ["python", "-c", "import snakemake_executor_plugin_slurm"],
        )
        checks.append(("snakemake-executor-plugin-slurm", plugin.returncode == 0, runner_check_detail(plugin)))
        sbatch = run_runner_command(
            runner,
            [
                "python",
                "-c",
                "import shutil, sys; sys.exit(0 if shutil.which('sbatch') else 1)",
            ],
        )
        checks.append(("sbatch visible from runner", sbatch.returncode == 0, runner_check_detail(sbatch)))

    source_errors = validate_microclean_source(REPO_ROOT / "workflow/envs/micRoclean-source.env")
    checks.append(("micRoclean source pins", not source_errors, "; ".join(source_errors)))

    ok = True
    for label, passed, detail in checks:
        marker = "ok" if passed else "fail"
        print(f"[{marker}] {label}" + (f": {detail}" if detail else ""))
        ok = ok and passed
    return 0 if ok else 1


def resolve_row_config(args: argparse.Namespace) -> Path:
    """Return an existing row config or create one from direct run fields."""
    if args.row_config:
        return Path(args.row_config)
    required = {
        "trial_id": args.trial_id,
        "trial_descript": args.trial_descript,
        "exp_dir": args.exp_dir,
        "metadata": args.metadata,
        "host": getattr(args, "host", None),
    }
    missing = [name.replace("_", "-") for name, value in required.items() if not value]
    if missing:
        raise SystemExit(
            "Provide --row-config or direct batch fields: "
            + ", ".join(f"--{name}" for name in missing)
        )
    try:
        host = normalize_host(args.host, context="direct run --host")
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    row = BatchRow(
        row_number=1,
        trialID=args.trial_id,
        trial_descript=args.trial_descript,
        exp_dir=args.exp_dir,
        metadata=args.metadata,
        host=host,
        process_umis=args.process_umis,
    )
    return write_run_config(row, args.run_config_dir)


def build_run_spec(args: argparse.Namespace, row_config: Path) -> RunSpec:
    """Resolve defaults and config layering for one Snakemake invocation."""
    if not row_config.exists():
        raise SystemExit(f"Missing row config: {row_config}")
    row_values = load_simple_run_config(row_config)
    require_row_config_host(row_values, row_config)
    trial_id = row_values.get("trialID") or args.trial_id or "run"
    trial_description = row_values.get("trial_descript") or getattr(args, "trial_descript", None) or ""
    trial_name = build_trial_name(trial_id, trial_description)
    isolated_workdir = should_use_isolated_workdir(args)
    profile = resolve_profile_path(args, isolated_workdir)
    # Config order is deliberate: base project defaults, optional user
    # overrides, then the generated per-row values for this batch.
    configfiles = resolve_configfiles(args.configfile, args.extra_configfile)
    validate_processing_config_stack(configfiles)
    configfiles.append(row_config)
    log_dir = Path(args.log_root) / trial_id
    job_name = args.job_name or f"low-bm-{trial_id}"
    workdir = resolve_isolated_snakemake_workdir(args, trial_name) if isolated_workdir else None
    snakefile = (REPO_ROOT / "Snakefile").resolve() if isolated_workdir else None
    conda_prefix = (
        repo_path(Path(getattr(args, "snakemake_conda_prefix", DEFAULT_SNAKEMAKE_CONDA_PREFIX))).resolve()
        if isolated_workdir
        else None
    )
    return RunSpec(
        configfiles=configfiles,
        target=args.target,
        profile=profile,
        mode=args.mode,
        log_dir=log_dir,
        trial_id=trial_id,
        job_name=job_name,
        dry_run=args.dry_run,
        unlock=args.unlock,
        extra_snakemake_args=list(args.snakemake_arg or []),
        trial_name=trial_name,
        workdir=workdir,
        snakefile=snakefile,
        snakemake_conda_prefix=conda_prefix,
    )


def require_row_config_host(row_values: dict[str, str], row_config: Path) -> None:
    """Fail fast when a per-batch row config does not choose a valid host."""
    host = row_values.get("host", "")
    if not host:
        raise SystemExit(
            f"Row config is missing required key 'host': {row_config}. "
            "Add a host column to the batch table or pass --host for direct runs."
        )
    try:
        normalize_host(host, context=f"row config {row_config} key 'host'")
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc


def build_trial_name(trial_id: str, trial_description: str) -> str:
    """Return the output trial name used by the Snakefile."""
    if trial_description:
        return f"{trial_id}_{trial_description}"
    return trial_id


def should_use_batch_isolated_workdirs(args: argparse.Namespace) -> bool:
    """Return the batch-submit default: isolated only for SLURM unless opted out."""
    return getattr(args, "mode", "slurm") == "slurm" and not bool(getattr(args, "shared_workdir", False))


def should_use_isolated_workdir(args: argparse.Namespace) -> bool:
    """Return whether this resolved Snakemake invocation should use --directory."""
    if bool(getattr(args, "shared_workdir", False)):
        return False
    return bool(getattr(args, "isolated_workdir", False))


def resolve_profile_path(args: argparse.Namespace, isolated_workdir: bool) -> Path:
    """Resolve the profile path, using absolute paths for isolated commands."""
    profile = Path(args.profile) if getattr(args, "profile", None) else REPO_ROOT / "profiles" / args.mode
    if isolated_workdir:
        return repo_path(profile).resolve()
    return profile


def resolve_isolated_snakemake_workdir(args: argparse.Namespace, trial_name: str) -> Path:
    """Resolve the per-trial Snakemake workdir under --workdir-root."""
    root = repo_path(Path(getattr(args, "workdir_root", DEFAULT_SNAKEMAKE_WORKDIR_ROOT))).resolve()
    return root / sanitize_path_slug(trial_name)


def sanitize_path_slug(value: str) -> str:
    """Keep trial-derived workdir names stable and shell/path friendly."""
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._-")
    slug = re.sub(r"_+", "_", slug)
    return slug or "run"


def validate_unique_batch_rows(rows: list[BatchRow]) -> None:
    """Reject batch rows that would collide in output or generated config paths."""
    duplicate_names = duplicate_values(row.trial_name for row in rows)
    if duplicate_names:
        raise SystemExit(
            "Duplicate batch output trial name(s): "
            + ", ".join(duplicate_names)
            + ". Each trialID + trial_descript pair must resolve to a unique output directory."
        )
    duplicate_ids = duplicate_values(row.trialID for row in rows)
    if duplicate_ids:
        raise SystemExit(
            "Duplicate trialID value(s): "
            + ", ".join(duplicate_ids)
            + ". Generated row configs are keyed by trialID, so trialID values must be unique."
        )
    duplicate_workdirs = duplicate_values(sanitize_path_slug(row.trial_name) for row in rows)
    if duplicate_workdirs:
        raise SystemExit(
            "Duplicate isolated Snakemake workdir slug(s): "
            + ", ".join(duplicate_workdirs)
            + ". Adjust trialID/trial_descript values so sanitized workdir names are unique."
        )


def duplicate_values(values) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def select_batch_unlock_rows(
    rows: list[BatchRow],
    trial_ids: list[str],
    unlock_all: bool,
) -> list[BatchRow]:
    """Return rows selected by `batch unlock` or raise for missing IDs."""
    if unlock_all:
        return rows
    return select_batch_rows_by_trial_ids(rows, trial_ids)


def select_batch_rows_by_trial_ids(rows: list[BatchRow], trial_ids: list[str]) -> list[BatchRow]:
    """Return all rows or the rows matching requested trial IDs."""
    if not trial_ids:
        return rows
    requested = set(trial_ids)
    selected = [row for row in rows if row.trialID in requested]
    found = {row.trialID for row in selected}
    missing = sorted(requested - found)
    if missing:
        raise SystemExit(
            "Requested trial ID(s) not found in batch table: " + ", ".join(missing)
        )
    return selected


def require_shared_reference_indexes(configfiles: list[Path], row_config: Path | None = None) -> None:
    """Fail when isolated batch runs are missing shared BWA indexes."""
    effective_configfiles = [*configfiles]
    if row_config is not None:
        effective_configfiles.append(row_config)
    refs = resolve_bwa_references(
        merge_top_level_yaml_scalar_values(effective_configfiles),
        all_configured_hosts=False,
    )
    failures = missing_bwa_reference_inputs(refs)
    if not failures:
        return
    raise SystemExit(
        render_bwa_reference_failure(failures)
        + "\nIsolated batch workdirs do not create or coordinate shared reference indexes."
    )


def merge_top_level_yaml_scalar_values(configfiles: list[Path]) -> dict[str, str]:
    values: dict[str, str] = {}
    for configfile in configfiles:
        values.update(top_level_yaml_scalar_values(repo_path(configfile)))
    return values


def resolve_configfiles(configfiles: list[str], extra_configfiles: list[str]) -> list[Path]:
    """Resolve and validate the processing config stack before row overrides."""
    resolved = [Path(path) for path in configfiles]
    if not resolved:
        default_config = Path(DEFAULT_BASE_CONFIG)
        if not path_exists_for_repo_run(default_config):
            raise SystemExit(
                f"Missing default processing config: {DEFAULT_BASE_CONFIG}. "
                f"Initialize it from {TEMPLATE_BASE_CONFIG}, edit it for this environment, "
                "or pass --configfile."
            )
        resolved = [default_config]

    for configfile in resolved:
        if not path_exists_for_repo_run(configfile):
            raise SystemExit(f"Missing config file: {configfile}")

    for extra_configfile in extra_configfiles:
        path = Path(extra_configfile)
        if not path_exists_for_repo_run(path):
            raise SystemExit(f"Missing extra config file: {path}")
        resolved.append(path)
    return resolved


def validate_processing_config_stack(configfiles: list[Path]) -> None:
    """Catch missing base processing configs before submitting a master job.

    Snakemake does not merge config files until it starts parsing the Snakefile.
    If the base processing layer is accidentally omitted, users otherwise see a
    remote master-job failure like `KeyError: in_root`. Checking the top-level
    keys here keeps the portability layer honest: the submitted job should fail
    only after the local launcher has proven its config stack is structurally
    complete.
    """
    seen_keys: set[str] = set()
    for configfile in configfiles:
        seen_keys.update(top_level_yaml_keys(repo_path(configfile)))

    missing = [key for key in REQUIRED_PROCESSING_CONFIG_KEYS if key not in seen_keys]
    if missing:
        raise SystemExit(
            "Processing config stack is missing required key(s): "
            + ", ".join(missing)
            + ". Make sure config/local/processing.yaml is included with --configfile "
            "and put run-specific path overrides under --extra-configfile."
        )


def top_level_yaml_keys(path: Path) -> set[str]:
    """Return simple top-level YAML keys without adding a PyYAML dependency."""
    keys: set[str] = set()
    for line in path.read_text().splitlines():
        if not line or line.startswith((" ", "\t", "#")):
            continue
        key, sep, _value = line.partition(":")
        if sep:
            keys.add(key.strip())
    return keys


def top_level_yaml_scalar_values(path: Path) -> dict[str, str]:
    """Return simple top-level scalar key/value pairs without PyYAML."""
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith((" ", "\t", "#")):
            continue
        key, sep, value = line.partition(":")
        if not sep:
            continue
        cleaned = clean_yaml_scalar(value)
        if cleaned:
            values[key.strip()] = cleaned
    return values


def top_level_yaml_section_values(path: Path, section: str) -> dict[str, str]:
    """Return simple scalar key/value pairs from one top-level YAML section."""
    values: dict[str, str] = {}
    in_section = False
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = len(line) - len(line.lstrip(" "))
        key, sep, value = stripped.partition(":")
        if indent == 0:
            in_section = sep == ":" and key.strip() == section
            continue
        if not in_section or sep != ":":
            continue
        if indent <= 2:
            values[key.strip()] = value.strip().split("#", 1)[0].strip()
    return values


def is_missing_yaml_scalar(value: str | None) -> bool:
    if value is None:
        return True
    value = clean_yaml_scalar(value)
    return not value or value.lower() in {"null", "none", "~"}


def clean_yaml_scalar(value: str | None) -> str:
    """Strip comments and simple YAML quotes from a scalar value."""
    if value is None:
        return ""
    value = value.split("#", 1)[0].strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        value = value[1:-1]
    return value


def ensure_batch_table_exists(batch_table: Path) -> None:
    """Raise a setup-oriented error when the default local batch table is absent."""
    if path_exists_for_repo_run(batch_table):
        return
    if str(batch_table) == DEFAULT_BATCH_TABLE:
        raise SystemExit(
            f"Missing default batch table: {DEFAULT_BATCH_TABLE}. "
            f"Initialize it from {TEMPLATE_BATCH_TABLE}, edit it for this experiment, "
            "or pass --batch-table."
        )
    raise SystemExit(f"Missing batch table: {batch_table}")


def resolve_analysis_config(path: str | Path) -> Path:
    """Return an existing analysis config or raise with setup guidance."""
    analysis_config = Path(path)
    if path_exists_for_repo_run(analysis_config):
        return analysis_config
    if str(analysis_config) == DEFAULT_ANALYSIS_CONFIG:
        raise SystemExit(
            f"Missing default analysis config: {DEFAULT_ANALYSIS_CONFIG}. "
            f"Initialize it with `low-bm analysis init`, copy {TEMPLATE_ANALYSIS_CONFIG}, "
            "or pass --analysis-config."
        )
    raise SystemExit(f"Missing analysis config: {analysis_config}")


def resolve_meta_config(path: str | Path) -> Path:
    """Return an existing meta-analysis config or raise with setup guidance."""
    meta_config = Path(path)
    if path_exists_for_repo_run(meta_config):
        return meta_config
    if str(meta_config) == DEFAULT_META_CONFIG:
        raise SystemExit(
            f"Missing default meta-analysis config: {DEFAULT_META_CONFIG}. "
            f"Initialize it with `low-bm meta init`, copy {TEMPLATE_META_CONFIG}, "
            "or pass --analysis-config."
        )
    raise SystemExit(f"Missing meta-analysis config: {meta_config}")


def expand_analysis_steps(steps: list[str]) -> list[str]:
    """Expand composite analysis aliases while preserving user order."""
    expanded: list[str] = []
    for step in steps:
        if step in ANALYSIS_COMPOSITES:
            expanded.extend(ANALYSIS_COMPOSITES[step])
        elif step in ANALYSIS_STEP_REGISTRY:
            expanded.append(step)
        else:
            raise SystemExit(
                f"Unsupported analysis step: {step}. "
                f"Supported steps: {', '.join(ANALYSIS_STEP_CHOICES)}"
            )
    return expanded


def validate_analysis_config(analysis_config: Path, expanded_steps: list[str]) -> None:
    """Check the config has the selected top-level analysis sections."""
    config_path = repo_path(analysis_config)
    keys = top_level_yaml_keys(config_path)
    if "project" not in keys:
        raise SystemExit(f"Analysis config is missing required top-level section: project")
    project_values = top_level_yaml_section_values(config_path, "project")

    missing_sections: list[str] = []
    for step in expanded_steps:
        step_spec = ANALYSIS_STEP_REGISTRY[step]
        if step_spec.required_section and step_spec.required_section not in keys:
            missing_sections.append(f"{step_spec.required_section} (for {step})")
        validate_analysis_step_assets(step_spec)

    if missing_sections:
        raise SystemExit(
            "Analysis config is missing section(s) required by selected step(s): "
            + ", ".join(missing_sections)
        )

    endpoint_steps = [step for step in expanded_steps if step in ANALYSIS_ENDPOINT_STEPS]
    if endpoint_steps and is_missing_yaml_scalar(project_values.get("compiled_physeq")):
        raise SystemExit(
            "Standard analysis steps require project.compiled_physeq in the analysis config. "
            "For multiple processing batches, run `low-bm meta compile-phyloseq` first, "
            "then point project.compiled_physeq at the compiled endpoint."
        )


def validate_meta_config(meta_config: Path, steps: list[str]) -> None:
    """Check the meta-analysis config has selected top-level sections."""
    config_path = repo_path(meta_config)
    keys = top_level_yaml_keys(config_path)
    if "project" not in keys:
        raise SystemExit("Meta-analysis config is missing required top-level section: project")

    missing_sections: list[str] = []
    for step in steps:
        step_spec = META_STEP_REGISTRY[step]
        if step_spec.required_section and step_spec.required_section not in keys:
            missing_sections.append(f"{step_spec.required_section} (for {step})")
        validate_analysis_step_assets(step_spec)

    if missing_sections:
        raise SystemExit(
            "Meta-analysis config is missing section(s) required by selected step(s): "
            + ", ".join(missing_sections)
        )


def validate_analysis_step_assets(step_spec: AnalysisStepSpec) -> None:
    """Fail early when the wrapper points at a missing local script or env file."""
    argv = list(step_spec.argv_template)
    script_paths = [part for part in argv if part.startswith("scripts/")]
    for script_path in script_paths:
        if not repo_path(Path(script_path)).exists():
            raise SystemExit(f"Missing analysis script for {step_spec.name}: {script_path}")
    if step_spec.env_file is not None and not repo_path(step_spec.env_file).exists():
        raise SystemExit(f"Missing analysis environment file for {step_spec.name}: {step_spec.env_file}")


def build_analysis_run_spec(
    args: argparse.Namespace,
    created_utc: datetime | None = None,
) -> AnalysisRunSpec:
    """Resolve config, steps, logs, env wrapping, and commands for analysis."""
    analysis_config = resolve_analysis_config(args.analysis_config)
    requested_steps = list(args.steps)
    expanded_steps = expand_analysis_steps(requested_steps)
    validate_analysis_config(analysis_config, expanded_steps)

    env_mode = getattr(args, "env_mode", "named")
    if env_mode not in ANALYSIS_ENV_MODES:
        raise SystemExit(f"Unsupported analysis env mode: {env_mode}")

    manager = (
        resolve_env_manager(getattr(args, "manager", "auto"))
        if env_mode in {"managed", "named", "prefix"}
        else None
    )
    if env_mode == "managed":
        env_prefixes = resolve_managed_analysis_env_prefixes(
            args,
            expanded_steps,
            ANALYSIS_STEP_REGISTRY,
        )
    elif env_mode == "prefix":
        env_prefixes = resolve_analysis_env_prefixes(
            args,
            expanded_steps,
            ANALYSIS_STEP_REGISTRY,
            validate_exists=not bool(getattr(args, "dry_run", False)),
        )
    else:
        env_prefixes = {}
    timestamp = created_utc or datetime.now(timezone.utc)
    log_dir = resolve_analysis_log_dir(args, expanded_steps, timestamp)

    commands = [
        build_analysis_command_spec(
            step=step,
            analysis_config=analysis_config,
            env_mode=env_mode,
            manager=manager,
            env_prefixes=env_prefixes,
            registry=ANALYSIS_STEP_REGISTRY,
        )
        for step in expanded_steps
    ]
    return AnalysisRunSpec(
        analysis_config=analysis_config,
        requested_steps=requested_steps,
        expanded_steps=expanded_steps,
        log_dir=log_dir,
        dry_run=bool(getattr(args, "dry_run", False)),
        env_mode=env_mode,
        manager=manager,
        commands=commands,
    )


def build_meta_run_spec(
    args: argparse.Namespace,
    step: str,
    created_utc: datetime | None = None,
) -> AnalysisRunSpec:
    """Resolve config, logs, env wrapping, and command for one meta step."""
    meta_config = resolve_meta_config(args.analysis_config)
    validate_meta_config(meta_config, [step])

    env_mode = getattr(args, "env_mode", "named")
    if env_mode not in ANALYSIS_ENV_MODES:
        raise SystemExit(f"Unsupported meta-analysis env mode: {env_mode}")

    manager = (
        resolve_env_manager(getattr(args, "manager", "auto"))
        if env_mode in {"managed", "named", "prefix"}
        else None
    )
    if env_mode == "managed":
        env_prefixes = resolve_managed_analysis_env_prefixes(
            args,
            [step],
            META_STEP_REGISTRY,
        )
    elif env_mode == "prefix":
        env_prefixes = resolve_analysis_env_prefixes(
            args,
            [step],
            META_STEP_REGISTRY,
            validate_exists=not bool(getattr(args, "dry_run", False)),
        )
    else:
        env_prefixes = {}
    timestamp = created_utc or datetime.now(timezone.utc)
    log_dir = resolve_analysis_log_dir(args, [step], timestamp)
    command = build_analysis_command_spec(
        step=step,
        analysis_config=meta_config,
        env_mode=env_mode,
        manager=manager,
        env_prefixes=env_prefixes,
        registry=META_STEP_REGISTRY,
    )
    return AnalysisRunSpec(
        analysis_config=meta_config,
        requested_steps=[step],
        expanded_steps=[step],
        log_dir=log_dir,
        dry_run=bool(getattr(args, "dry_run", False)),
        env_mode=env_mode,
        manager=manager,
        commands=[command],
    )


def resolve_analysis_log_dir(
    args: argparse.Namespace,
    expanded_steps: list[str],
    created_utc: datetime,
) -> Path:
    if getattr(args, "log_dir", None):
        return Path(args.log_dir)
    step_slug = "-".join(expanded_steps)
    step_slug = re.sub(r"[^A-Za-z0-9._-]+", "_", step_slug) or "analysis"
    timestamp = created_utc.strftime("%Y%m%dT%H%M%SZ")
    return Path(getattr(args, "log_root", DEFAULT_ANALYSIS_LOG_ROOT)) / f"{timestamp}_{step_slug}"


def build_analysis_command_spec(
    step: str,
    analysis_config: Path,
    env_mode: str,
    manager: EnvManagerSpec | None,
    env_prefixes: dict[Path, Path] | None = None,
    registry: dict[str, AnalysisStepSpec] = ANALYSIS_STEP_REGISTRY,
) -> AnalysisCommandSpec:
    step_spec = registry[step]
    command = [
        token.format(analysis_config=str(analysis_config))
        for token in step_spec.argv_template
    ]
    env_name = analysis_env_name(step_spec.env_file) if step_spec.env_file is not None else None
    env_prefix = None
    if step_spec.env_file is not None and env_prefixes is not None:
        env_prefix = env_prefixes.get(step_spec.env_file)
    execution_command = wrap_analysis_command(command, env_name, env_prefix, env_mode, manager)
    return AnalysisCommandSpec(
        step=step,
        command=command,
        execution_command=execution_command,
        env_name=env_name,
        env_prefix=env_prefix,
        env_file=step_spec.env_file,
    )


def wrap_analysis_command(
    command: list[str],
    env_name: str | None,
    env_prefix: Path | None,
    env_mode: str,
    manager: EnvManagerSpec | None,
) -> list[str]:
    if env_mode == "direct" or (env_name is None and env_prefix is None):
        return command
    if manager is None:
        raise SystemExit(f"--env-mode {env_mode} requires a conda-compatible manager.")
    if env_mode in {"managed", "prefix"}:
        if env_prefix is None:
            raise SystemExit(f"--env-mode {env_mode} requires an environment prefix for each selected step.")
        return [
            manager.executable,
            "run",
            "--prefix",
            str(env_prefix),
            *command,
        ]
    return [
        manager.executable,
        "run",
        "--name",
        env_name,
        *command,
    ]


def resolve_analysis_env_prefixes(
    args: argparse.Namespace,
    steps: list[str],
    registry: dict[str, AnalysisStepSpec],
    validate_exists: bool,
) -> dict[Path, Path]:
    """Resolve explicit conda prefixes needed by the selected analysis steps."""
    prefixes: dict[Path, Path] = {}
    for env_file in selected_analysis_env_files(steps, registry):
        option_name, envvar_name, label = analysis_env_prefix_source(env_file)
        value = getattr(args, option_name, None)
        if not value:
            raise SystemExit(
                f"--env-mode prefix requires {label} for steps that use {analysis_env_name(env_file)} "
                f"(or set ${envvar_name})."
            )
        prefix = repo_path(Path(value)).resolve()
        if validate_exists and not prefix.is_dir():
            raise SystemExit(f"Missing analysis environment prefix for {analysis_env_name(env_file)}: {prefix}")
        prefixes[env_file] = prefix
    return prefixes


def resolve_managed_analysis_env_prefixes(
    args: argparse.Namespace,
    steps: list[str],
    registry: dict[str, AnalysisStepSpec],
) -> dict[Path, Path]:
    """Resolve project-local prefixes derived from selected env YAMLs."""
    env_root = repo_path(Path(getattr(args, "analysis_env_root", DEFAULT_ANALYSIS_ENV_ROOT))).resolve()
    return {
        env_file: managed_analysis_env_prefix(env_file, env_root)
        for env_file in selected_analysis_env_files(steps, registry)
    }


def managed_analysis_env_prefix(env_file: Path, env_root: Path) -> Path:
    """Return a deterministic prefix for one analysis environment YAML."""
    env_name = re.sub(r"[^A-Za-z0-9._-]+", "-", analysis_env_name(env_file)).strip("-")
    env_hash = hash_file(repo_path(env_file))[:12]
    return env_root / f"{env_name}-{env_hash}"


def selected_analysis_env_files(
    steps: list[str],
    registry: dict[str, AnalysisStepSpec],
) -> list[Path]:
    """Return unique environment YAML paths used by selected steps."""
    env_files: list[Path] = []
    seen: set[Path] = set()
    for step in steps:
        env_file = registry[step].env_file
        if env_file is not None and env_file not in seen:
            env_files.append(env_file)
            seen.add(env_file)
    return env_files


def analysis_env_prefix_source(env_file: Path) -> tuple[str, str, str]:
    """Map an analysis env file to its CLI/env-var prefix source."""
    if env_file == R_TOOLS_ENV:
        return "r_env_prefix", R_TOOLS_ENV_PREFIX_ENVVAR, "--r-env-prefix"
    if env_file == BIO_TOOLS_ENV:
        return "bio_env_prefix", BIO_TOOLS_ENV_PREFIX_ENVVAR, "--bio-env-prefix"
    if env_file == MICROCLEAN_ENV:
        return "microclean_env_prefix", MICROCLEAN_ENV_PREFIX_ENVVAR, "--microclean-env-prefix"
    raise SystemExit(f"No prefix option is registered for analysis environment file: {env_file}")


def ensure_managed_analysis_envs(spec: AnalysisRunSpec, label: str) -> int:
    """Create any missing project-managed analysis environments."""
    if spec.env_mode != "managed":
        return 0
    if spec.manager is None:
        raise SystemExit("--env-mode managed requires a conda-compatible manager.")

    required_envs = required_managed_analysis_envs(spec)
    if not required_envs:
        return 0

    log_file = spec.log_dir / "analysis_env_setup.log"
    with log_file.open("a") as log:
        for env_file, prefix in required_envs.items():
            if prefix.is_dir():
                log.write(f"# Reusing {analysis_env_name(env_file)} at {prefix}\n")
                continue

            prefix.parent.mkdir(parents=True, exist_ok=True)
            env_path = repo_path(env_file)
            command = build_conda_env_create_command(prefix, env_path, spec.manager)
            print(f"Creating {label.lower()} environment {analysis_env_name(env_file)} at {prefix}.")
            log.write(f"\n# {datetime.now(timezone.utc).isoformat()}\n")
            log.write(f"$ {shlex.join(command)}\n")
            log.flush()
            completed = subprocess.run(
                command,
                cwd=REPO_ROOT,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
            if completed.returncode != 0:
                print(
                    f"{label} environment setup failed for {analysis_env_name(env_file)} "
                    f"with exit code {completed.returncode}. See {log_file}.",
                    file=sys.stderr,
                )
                return completed.returncode
            post_deploy_returncode = run_analysis_env_post_deploy(prefix, env_file, spec.manager, log)
            if post_deploy_returncode != 0:
                print(
                    f"{label} environment post-deploy failed for {analysis_env_name(env_file)} "
                    f"with exit code {post_deploy_returncode}. See {log_file}.",
                    file=sys.stderr,
                )
                return post_deploy_returncode
            write_analysis_env_metadata(prefix, env_file, spec.manager)
    return 0


def analysis_env_post_deploy_script(env_file: Path) -> Path | None:
    """Return a Snakemake-style env post-deploy script when the env has one."""
    script = repo_path(env_file).with_suffix(".post-deploy.sh")
    return script if script.exists() else None


def run_analysis_env_post_deploy(
    prefix: Path,
    env_file: Path,
    manager: EnvManagerSpec,
    log,
) -> int:
    """Run a managed analysis env's post-deploy script inside the new prefix."""
    script = analysis_env_post_deploy_script(env_file)
    if script is None:
        return 0

    command = [
        manager.executable,
        "run",
        "--prefix",
        str(prefix),
        "bash",
        str(script),
    ]
    log.write(f"$ {shlex.join(command)}\n")
    log.flush()
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return completed.returncode


def required_managed_analysis_envs(spec: AnalysisRunSpec) -> dict[Path, Path]:
    """Return unique env YAML to prefix mappings from an analysis run spec."""
    required: dict[Path, Path] = {}
    for command_spec in spec.commands:
        if command_spec.env_file is not None and command_spec.env_prefix is not None:
            required[command_spec.env_file] = command_spec.env_prefix
    return required


def write_analysis_env_metadata(prefix: Path, env_file: Path, manager: EnvManagerSpec) -> None:
    """Write a small sidecar showing which YAML produced a managed env."""
    metadata_path = prefix.parent / f"{prefix.name}.json"
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "env_prefix": str(prefix),
        "env_file": str(env_file),
        "env_file_sha256": hash_file(repo_path(env_file)),
        "manager": manager.name,
        "manager_executable": manager.executable,
        "low_bm_version": __version__,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")


def analysis_env_name(env_file: Path) -> str:
    """Read the top-level conda environment name from an env YAML file."""
    path = repo_path(env_file)
    for line in path.read_text().splitlines():
        if line.startswith("name:"):
            value = line.split(":", 1)[1].strip().strip('"').strip("'")
            if value:
                return value
    raise SystemExit(f"Missing top-level name in analysis environment file: {env_file}")


def write_analysis_provenance(spec: AnalysisRunSpec) -> None:
    """Record enough metadata to reconstruct an analysis invocation."""
    command_entries = [
        {
            "step": command_spec.step,
            "command": command_spec.command,
            "execution_command": command_spec.execution_command,
            "env_name": command_spec.env_name,
            "env_prefix": str(command_spec.env_prefix) if command_spec.env_prefix else None,
            "env_file": str(command_spec.env_file) if command_spec.env_file else None,
            "env_file_sha256": hash_file(repo_path(command_spec.env_file)) if command_spec.env_file else None,
        }
        for command_spec in spec.commands
    ]
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "low_bm_version": __version__,
        "analysis_config": str(spec.analysis_config),
        "requested_steps": spec.requested_steps,
        "expanded_steps": spec.expanded_steps,
        "dry_run": spec.dry_run,
        "env_mode": spec.env_mode,
        "manager": spec.manager.name if spec.manager else None,
        "manager_executable": spec.manager.executable if spec.manager else None,
        "git_commit": git_commit(),
        "commands": command_entries,
    }
    (spec.log_dir / "analysis_metadata.json").write_text(json.dumps(provenance, indent=2) + "\n")
    command_lines = []
    for command_spec in spec.commands:
        command_lines.append(f"# {command_spec.step}")
        command_lines.append(shlex.join(command_spec.execution_command))
    (spec.log_dir / "analysis_commands.txt").write_text("\n".join(command_lines) + "\n")


def write_meta_provenance(spec: AnalysisRunSpec) -> None:
    """Record enough metadata to reconstruct a meta-analysis invocation."""
    command_entries = [
        {
            "step": command_spec.step,
            "command": command_spec.command,
            "execution_command": command_spec.execution_command,
            "env_name": command_spec.env_name,
            "env_prefix": str(command_spec.env_prefix) if command_spec.env_prefix else None,
            "env_file": str(command_spec.env_file) if command_spec.env_file else None,
            "env_file_sha256": hash_file(repo_path(command_spec.env_file)) if command_spec.env_file else None,
        }
        for command_spec in spec.commands
    ]
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "low_bm_version": __version__,
        "meta_config": str(spec.analysis_config),
        "requested_steps": spec.requested_steps,
        "expanded_steps": spec.expanded_steps,
        "dry_run": spec.dry_run,
        "env_mode": spec.env_mode,
        "manager": spec.manager.name if spec.manager else None,
        "manager_executable": spec.manager.executable if spec.manager else None,
        "git_commit": git_commit(),
        "commands": command_entries,
    }
    (spec.log_dir / "meta_metadata.json").write_text(json.dumps(provenance, indent=2) + "\n")
    command_lines = []
    for command_spec in spec.commands:
        command_lines.append(f"# {command_spec.step}")
        command_lines.append(shlex.join(command_spec.execution_command))
    (spec.log_dir / "meta_commands.txt").write_text("\n".join(command_lines) + "\n")


def analysis_step_log_file(log_dir: Path, index: int, step: str) -> Path:
    safe_step = re.sub(r"[^A-Za-z0-9._-]+", "_", step)
    return log_dir / f"{index:02d}_{safe_step}.log"


def run_analysis_step(command_spec: AnalysisCommandSpec, log_dir: Path, index: int) -> int:
    """Run one analysis command and capture stdout/stderr in its step log."""
    log_file = analysis_step_log_file(log_dir, index, command_spec.step)
    with log_file.open("a") as log:
        log.write(f"\n# {datetime.now(timezone.utc).isoformat()}\n")
        log.write(f"$ {shlex.join(command_spec.execution_command)}\n")
        log.flush()
        completed = subprocess.run(
            command_spec.execution_command,
            cwd=REPO_ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return completed.returncode


def path_exists_for_repo_run(path: Path) -> bool:
    """Return true when a path exists as given or relative to the repo root."""
    return path.exists() or (not path.is_absolute() and (REPO_ROOT / path).exists())


def repo_path(path: Path) -> Path:
    """Resolve a path as written by the user or relative to the repository."""
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def runner_metadata_path(prefix: Path) -> Path:
    """Return the metadata file beside the project-local runner prefix."""
    return prefix.parent / "runner.json"


def resolve_env_manager(manager_name: str = "auto", metadata_path: Path | None = None) -> EnvManagerSpec:
    """Find a conda-compatible manager.

    We prefer existing tools over bootstrapping micromamba in V1. That keeps the
    first portability step easy to understand: the repo owns the runner prefix,
    while the user or HPC module system provides the package manager.
    """
    if manager_name != "auto":
        executable = shutil.which(manager_name)
        if not executable:
            raise SystemExit(
                f"Could not find requested environment manager: {manager_name}. "
                "Load it on PATH or rerun with --manager auto."
            )
        return EnvManagerSpec(name=manager_name, executable=executable)

    metadata = load_runner_metadata(metadata_path) if metadata_path else {}
    recorded_executable = metadata.get("manager_executable")
    recorded_name = metadata.get("manager")
    if recorded_executable and Path(recorded_executable).exists():
        return EnvManagerSpec(name=str(recorded_name or Path(recorded_executable).name), executable=str(recorded_executable))

    for candidate in ENV_MANAGERS:
        executable = shutil.which(candidate)
        if executable:
            return EnvManagerSpec(name=candidate, executable=executable)
    raise SystemExit(
        "Could not find mamba, conda, or micromamba on PATH. "
        "Load one of those tools, then run `low-bm setup runner`."
    )


def load_runner_metadata(metadata_path: Path | None) -> dict[str, object]:
    if not metadata_path:
        return {}
    path = repo_path(metadata_path)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}


def write_runner_metadata(
    prefix_arg: Path,
    env_file: Path,
    manager: EnvManagerSpec,
    conda_base_prefix: Path,
    *,
    preserve_created: bool,
) -> None:
    """Write runner provenance, preserving creation time for existing envs."""
    metadata_path = runner_metadata_path(prefix_arg)
    existing = load_runner_metadata(metadata_path) if preserve_created else {}
    metadata = {
        "created_utc": existing.get("created_utc") or datetime.now(timezone.utc).isoformat(),
        "runner_prefix": str(prefix_arg),
        "env_file": DEFAULT_RUNNER_ENV_FILE,
        "env_file_sha256": hash_file(env_file),
        "manager": manager.name,
        "manager_executable": manager.executable,
        "conda_base_prefix": str(conda_base_prefix),
        "low_bm_version": __version__,
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")


def resolve_path_prefix(value: str | Path) -> Path:
    """Resolve user-supplied prefix paths without treating absolute paths as repo-relative."""
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = repo_path(path)
    return path.resolve()


def resolve_conda_base_prefix(
    override: str | None,
    metadata_path: Path | None,
    manager: EnvManagerSpec,
) -> Path:
    """Resolve the real conda base Snakemake should use for rule env activation."""
    if override:
        return resolve_path_prefix(override)

    metadata = load_runner_metadata(metadata_path)
    recorded = metadata.get("conda_base_prefix")
    if recorded:
        return resolve_path_prefix(str(recorded))

    return detect_conda_base_prefix(manager)


def detect_conda_base_prefix(manager: EnvManagerSpec) -> Path:
    """Ask the configured conda-compatible manager for its base prefix."""
    completed = subprocess.run(
        [manager.executable, "info", "--base"],
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode == 0 and completed.stdout.strip():
        return resolve_path_prefix(completed.stdout.strip().splitlines()[-1])

    json_completed = subprocess.run(
        [manager.executable, "info", "--json"],
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if json_completed.returncode == 0 and json_completed.stdout.strip():
        try:
            payload = json.loads(json_completed.stdout)
        except json.JSONDecodeError:
            payload = {}
        for key in ("root_prefix", "conda_prefix", "base_prefix"):
            value = payload.get(key)
            if value:
                return resolve_path_prefix(str(value))

    detail = completed.stderr.strip() or json_completed.stderr.strip() or "manager did not report a base prefix"
    raise SystemExit(
        f"Could not detect conda base prefix from {manager.executable}: {detail}. "
        "Rerun with --conda-base-prefix /path/to/conda/base."
    )


def conda_base_activation_errors(prefix: Path) -> list[str]:
    """Return missing pieces needed by Snakemake's conda activation path."""
    errors = []
    if not prefix.is_dir():
        errors.append(f"missing conda base directory {prefix}")
        return errors
    activate = prefix / "bin" / "activate"
    conda = prefix / "bin" / "conda"
    if not activate.is_file():
        errors.append(f"missing conda activation script {activate}")
    if not conda.is_file():
        errors.append(f"missing conda executable {conda}")
    elif not os.access(conda, os.X_OK):
        errors.append(f"conda executable is not executable {conda}")
    return errors


def build_runner_create_command(prefix: Path, env_file: Path, manager: EnvManagerSpec) -> list[str]:
    """Build the command that creates the project-local runner environment."""
    return build_conda_env_create_command(prefix, env_file, manager)


def build_conda_env_create_command(prefix: Path, env_file: Path, manager: EnvManagerSpec) -> list[str]:
    """Build a conda-compatible environment create command."""
    if is_micromamba_manager(manager):
        return [
            manager.executable,
            "create",
            "--yes",
            "--prefix",
            str(prefix),
            "--file",
            str(env_file),
        ]
    return [
        manager.executable,
        "env",
        "create",
        "--yes",
        "--prefix",
        str(prefix),
        "--file",
        str(env_file),
    ]


def is_micromamba_manager(manager: EnvManagerSpec) -> bool:
    """Return true for micromamba even when the executable was passed by path."""
    return manager.name == "micromamba" or Path(manager.executable).name == "micromamba"


def remove_runner_prefix(prefix: Path) -> None:
    """Remove an existing runner prefix after a few guardrails."""
    target = repo_path(prefix).resolve()
    blocked = {REPO_ROOT.resolve(), Path.home().resolve(), Path("/").resolve()}
    if target in blocked or ".low-bm" not in target.parts:
        raise SystemExit(f"Refusing to remove unsafe runner prefix: {target}")
    shutil.rmtree(target)


def runner_bin_dir(prefix: Path) -> Path:
    """Return the executable directory for the project-local runner prefix."""
    return prefix / ("Scripts" if os.name == "nt" else "bin")


def required_runner_executables(runner: HostRunnerSpec) -> list[Path]:
    """Return runner executables needed for direct processing launches."""
    return [runner.bin_dir / name for name in ("snakemake", "python")]


def validate_host_runner_executables(runner: HostRunnerSpec) -> None:
    """Fail early when a runner prefix is missing direct entrypoints."""
    missing = [path for path in required_runner_executables(runner) if not path.exists()]
    if missing:
        rendered = ", ".join(str(path) for path in missing)
        raise SystemExit(
            f"Missing runner executable(s): {rendered}. "
            "Run `low-bm setup runner --force` to recreate the runner prefix."
        )


def runner_process_env(runner: HostRunnerSpec) -> dict[str, str]:
    """Return an environment where runner-prefix tools win PATH lookup."""
    env = os.environ.copy()
    env["PATH"] = str(runner.bin_dir) + os.pathsep + env.get("PATH", "")
    env.pop("PYTHONTZPATH", None)
    return env


def snakemake_process_env(runner: HostRunnerSpec) -> dict[str, str]:
    """Return an environment for Snakemake without leaking runner conda."""
    env = os.environ.copy()
    for key in list(env):
        if key.startswith("CONDA_"):
            env.pop(key, None)
    env.pop("PYTHONTZPATH", None)
    path_parts = [part for part in env.get("PATH", "").split(os.pathsep) if part]
    filtered_parts = [part for part in path_parts if part != str(runner.bin_dir)]
    conda_paths = [runner.conda_base_prefix / "bin", runner.conda_base_prefix / "condabin"]
    env["PATH"] = os.pathsep.join([str(path) for path in conda_paths] + filtered_parts)
    return env


def execution_environment(args: argparse.Namespace) -> dict[str, str] | None:
    """Return subprocess env overrides for default direct-runner execution."""
    if getattr(args, "activate_command", ""):
        return None
    runner = resolve_host_runner(args)
    return snakemake_process_env(runner)


def resolve_host_runner(args: argparse.Namespace) -> HostRunnerSpec:
    """Resolve the V1 host runner.

    A prefix path is more reproducible than a global env name because it is tied
    to this checkout. Two projects can carry different runner prefixes without
    fighting over a shared `low-bm-runner` name in the user's conda install.
    """
    if getattr(args, "runner", "host") != "host":
        raise SystemExit("V1 only supports --runner host. Container runner support will be added later.")
    prefix_arg = Path(getattr(args, "runner_prefix", DEFAULT_RUNNER_PREFIX))
    prefix = repo_path(prefix_arg)
    metadata_path = runner_metadata_path(prefix_arg)
    if not prefix.is_dir():
        raise SystemExit(
            f"Missing runner env: {prefix_arg}. "
            "Run `low-bm setup runner` before running the workflow, "
            "or use --activate-command for an advanced submitted-SLURM fallback."
        )
    manager = resolve_env_manager(getattr(args, "manager", "auto"), metadata_path=metadata_path)
    resolved_prefix = prefix.resolve()
    runner = HostRunnerSpec(
        prefix=resolved_prefix,
        manager=manager,
        metadata_path=metadata_path,
        bin_dir=runner_bin_dir(resolved_prefix),
        conda_base_prefix=resolve_conda_base_prefix(
            getattr(args, "conda_base_prefix", None),
            metadata_path,
            manager,
        ),
    )
    validate_host_runner_executables(runner)
    return runner


def wrap_runner_command(runner: HostRunnerSpec, command: list[str]) -> list[str]:
    """Run known runner commands directly from the project-local prefix.

    Direct entrypoints avoid concurrent `mamba run` processes contending for
    the user's global mamba lock directory while preserving the pinned runner
    environment.
    """
    if not command:
        return command
    executable = command[0]
    if executable in {"snakemake", "python"}:
        return [str(runner.bin_dir / executable), *command[1:]]
    return command


def wrap_snakemake_command(runner: HostRunnerSpec, command: list[str]) -> list[str]:
    """Wrap a Snakemake argv list in the active runner implementation."""
    wrapped = wrap_runner_command(runner, command)
    if "--conda-base-path" not in wrapped:
        wrapped = [wrapped[0], "--conda-base-path", str(runner.conda_base_prefix), *wrapped[1:]]
    return wrapped


def run_runner_command(runner: HostRunnerSpec, command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        wrap_runner_command(runner, command),
        cwd=REPO_ROOT,
        env=runner_process_env(runner),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def run_rule_env_smoke_test(runner: HostRunnerSpec) -> tuple[bool, str]:
    """Run a tiny conda-backed rule and verify it does not use runner Python."""
    with tempfile.TemporaryDirectory(prefix="low-bm-rule-env-smoke-") as tmp:
        root = Path(tmp)
        workdir = root / "workdir"
        profile = root / "profile"
        conda_prefix = root / "snakemake-conda"
        profile.mkdir()
        workdir.mkdir()
        conda_prefix.mkdir()
        (root / "env.yaml").write_text(
            "channels:\n"
            "  - conda-forge\n"
            "dependencies:\n"
            "  - python=3.11\n"
        )
        (profile / "config.yaml").write_text(
            "executor: local\n"
            "software-deployment-method:\n"
            "  - conda\n"
            "cores: 1\n"
        )
        snakefile = root / "Snakefile"
        snakefile.write_text(
            "rule all:\n"
            "    input:\n"
            "        'python_path.txt'\n\n"
            "rule write_python:\n"
            "    output:\n"
            "        'python_path.txt'\n"
            "    conda:\n"
            "        'env.yaml'\n"
            "    shell:\n"
            "        \"python -c 'import sys; open(\\\"{output}\\\", \\\"w\\\").write(sys.executable)'\"\n"
        )
        command = wrap_snakemake_command(
            runner,
            [
                "snakemake",
                "--snakefile",
                str(snakefile),
                "--directory",
                str(workdir),
                "--profile",
                str(profile),
                "--conda-prefix",
                str(conda_prefix),
            ],
        )
        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=snakemake_process_env(runner),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if completed.returncode != 0:
            return False, completed.stdout.strip()
        python_path = (workdir / "python_path.txt").read_text().strip()
        if str(runner.prefix) in python_path:
            return False, f"rule used runner Python: {python_path}"
        if str(conda_prefix) not in python_path:
            return False, f"rule Python was outside smoke conda prefix: {python_path}"
        return True, python_path


def check_runner_command(runner: HostRunnerSpec, command: list[str], label: str) -> int:
    completed = run_runner_command(runner, command)
    if completed.returncode == 0:
        detail = completed.stdout.strip()
        print(f"[ok] {label}" + (f": {detail}" if detail else ""))
    else:
        print(f"[fail] {label}", file=sys.stderr)
        if completed.stderr:
            print(completed.stderr, file=sys.stderr, end="")
    return completed.returncode


def runner_check_detail(completed: subprocess.CompletedProcess[str]) -> str:
    """Return useful doctor output for both passing and failing checks."""
    if completed.returncode == 0:
        return completed.stdout.strip()
    return completed.stderr.strip() or completed.stdout.strip()


def validate_microclean_source(path: Path) -> list[str]:
    """Validate fixed source pins used by micRoclean's post-deploy step."""
    if not path.exists():
        return [f"missing {path}"]
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or "=" not in stripped:
            continue
        key, value = stripped.split("=", 1)
        values[key.strip()] = value.strip().strip('"').strip("'")

    errors = []
    for key in ("MICROCLEAN_GIT_REF", "SCRUB_GIT_REF"):
        value = values.get(key, "")
        if not GIT_SHA_RE.match(value):
            errors.append(f"{key} must be a fixed Git SHA")
    for key in ("MICROCLEAN_GIT_URL", "SCRUB_GIT_URL"):
        if not values.get(key):
            errors.append(f"{key} is required")
    return errors


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_snakemake_command(spec: RunSpec) -> list[str]:
    """Build the Snakemake argv list without executing it."""
    command = ["snakemake"]
    if spec.snakefile:
        command.extend(["--snakefile", str(spec.snakefile)])
    if spec.workdir:
        command.extend(["--directory", str(spec.workdir)])
    command.extend(["--profile", str(spec.profile)])
    if spec.snakemake_conda_prefix:
        command.extend(["--conda-prefix", str(spec.snakemake_conda_prefix)])
    command.extend(["--config", f"low_bm_repo_root={REPO_ROOT.resolve()}"])
    if spec.dry_run:
        command.append("--dry-run")
    if spec.unlock:
        command.append("--unlock")
    command.extend(spec.extra_snakemake_args)
    if spec.configfiles:
        # Snakemake 9 parses --configfile as one option followed by one or
        # more files. Keep the whole config stack under a single flag so later
        # files override earlier files in the documented order.
        command.append("--configfile")
        command.extend(format_snakemake_path(configfile, absolute=bool(spec.workdir)) for configfile in spec.configfiles)
    if not spec.unlock and spec.target:
        # The explicit "--" ends --configfile's variable-length file list.
        # Without it, a target like "all" can be swallowed as another config
        # path; placing the target before --configfile prevented Snakemake from
        # applying the config stack in Snakemake 9.23.
        command.extend(["--", spec.target])
    return command


def format_snakemake_path(path: Path, absolute: bool) -> str:
    """Render a command path, anchoring repo-relative paths when isolated."""
    if absolute:
        return str(repo_path(path).resolve())
    return str(path)


def build_execution_command(spec: RunSpec, args: argparse.Namespace) -> list[str]:
    """Build the command that should actually run on this system.

    Snakemake remains the inner command because it expresses workflow intent.
    The runner wrapper expresses portability: V1 uses a host conda prefix, while
    a later container runner can wrap this same inner command with apptainer exec.
    """
    command = build_snakemake_command(spec)
    activate_command = getattr(args, "activate_command", "")
    if activate_command:
        if spec.mode != "slurm" or spec.dry_run or spec.unlock:
            raise SystemExit(
                "--activate-command is only supported for submitted SLURM master jobs. "
                "Run `low-bm setup runner` for dry-runs, unlocks, and local execution."
            )
        return command
    runner = resolve_host_runner(args)
    return wrap_snakemake_command(runner, command)


def execute_run_spec(spec: RunSpec, args: argparse.Namespace) -> int:
    """Record provenance, then run locally or submit a SLURM master job."""
    spec.log_dir.mkdir(parents=True, exist_ok=True)
    prepare_snakemake_runtime_dirs(spec)
    command = build_execution_command(spec, args)
    write_provenance(spec, args, command)

    if args.print_command:
        print(shlex.join(command))

    # Dry-runs and unlocks should happen immediately in the current shell; they
    # are Snakemake control operations, not cluster work to submit.
    if spec.mode == "local" or spec.dry_run or spec.unlock:
        return run_local_command(command, spec, args)
    return submit_master_job(command, spec, args)


def prepare_snakemake_runtime_dirs(spec: RunSpec) -> None:
    """Create isolated workdir and shared conda cache roots before execution."""
    if spec.workdir:
        spec.workdir.mkdir(parents=True, exist_ok=True)
    if spec.snakemake_conda_prefix and not spec.unlock:
        spec.snakemake_conda_prefix.mkdir(parents=True, exist_ok=True)


def run_local_command(command: list[str], spec: RunSpec, args: argparse.Namespace) -> int:
    """Execute Snakemake directly and tee all output into the run log."""
    log_file = spec.log_dir / ("snakemake.dryrun.log" if spec.dry_run else "snakemake.log")
    with log_file.open("a") as log:
        log.write(f"\n# {datetime.now(timezone.utc).isoformat()}\n")
        log.write(f"$ {shlex.join(command)}\n")
        log.flush()
        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=execution_environment(args),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        print(f"Snakemake exited with {completed.returncode}. See {log_file}.", file=sys.stderr)
    return completed.returncode


def submit_master_job(command: list[str], spec: RunSpec, args: argparse.Namespace) -> int:
    """Write and submit the SLURM master job that will launch Snakemake."""
    script = write_master_script(command, spec, args)
    sbatch_command = ["sbatch", "--parsable", str(script)]
    if args.print_command:
        print(shlex.join(sbatch_command))
    completed = subprocess.run(
        sbatch_command,
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        if completed.stderr:
            print(completed.stderr, file=sys.stderr, end="")
        return completed.returncode
    job_id = completed.stdout.strip()
    (spec.log_dir / "master_job_id.txt").write_text(job_id + "\n")
    print(f"[{spec.trial_id}] submitted master job {job_id}")
    return 0


def master_script_path(spec: RunSpec) -> Path:
    """Return the stable path for the generated master-job script."""
    return spec.log_dir / "master_job.sh"


def write_master_script(command: list[str], spec: RunSpec, args: argparse.Namespace) -> Path:
    """Create the lightweight SLURM script that runs Snakemake once.

    The resources here are for the master process only. Rule-level resources
    are managed by the Snakemake profile, especially profiles/slurm/config.yaml.
    """
    script = master_script_path(spec)
    script.parent.mkdir(parents=True, exist_ok=True)
    output_log = spec.log_dir / "master-%j.log"
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={spec.job_name}",
        f"#SBATCH --cpus-per-task={args.master_cpus}",
        f"#SBATCH --mem={args.master_mem}",
        f"#SBATCH --time={args.master_time}",
        f"#SBATCH --output={output_log}",
        f"#SBATCH --error={output_log}",
    ]
    if args.master_partition:
        lines.append(f"#SBATCH --partition={args.master_partition}")
    for directive in args.master_extra_sbatch or []:
        lines.append(f"#SBATCH {directive}")
    lines.extend(
        [
            "set -euo pipefail",
            f"cd {shlex.quote(str(REPO_ROOT))}",
        ]
    )
    activate_command = getattr(args, "activate_command", "")
    if activate_command:
        # Keep environment activation configurable because HPC shell setup tends
        # to be site- and account-specific.
        lines.append(activate_command)
    else:
        prefix = repo_path(Path(args.runner_prefix)).resolve()
        bin_dir = runner_bin_dir(prefix)
        metadata_path = runner_metadata_path(Path(args.runner_prefix))
        manager = resolve_env_manager(getattr(args, "manager", "auto"), metadata_path=metadata_path)
        conda_base = resolve_conda_base_prefix(getattr(args, "conda_base_prefix", None), metadata_path, manager)
        # A submitted master job starts later, possibly on another node. This
        # guard catches deleted or unmounted project-local runner envs before a
        # confusing Snakemake failure scrolls by.
        missing_message = f"Missing runner env: {prefix}. Run low-bm setup runner first."
        conda_message = f"Missing conda base activation support: {conda_base}. Rerun low-bm setup runner --conda-base-prefix."
        lines.extend(
            [
                f"if [[ ! -d {shlex.quote(str(prefix))} ]]; then",
                f"  echo {shlex.quote(missing_message)} >&2",
                "  exit 2",
                "fi",
                f"if [[ ! -x {shlex.quote(str(bin_dir / 'snakemake'))} ]]; then",
                (
                    "  echo "
                    + shlex.quote(
                        "Runner env is missing snakemake. Run low-bm setup runner --force."
                    )
                    + " >&2"
                ),
                "  exit 2",
                "fi",
                (
                    f"if [[ ! -f {shlex.quote(str(conda_base / 'bin' / 'activate'))} "
                    f"|| ! -x {shlex.quote(str(conda_base / 'bin' / 'conda'))} ]]; then"
                ),
                f"  echo {shlex.quote(conda_message)} >&2",
                "  exit 2",
                "fi",
                "for var_name in ${!CONDA_@}; do unset \"$var_name\"; done",
                f"export PATH={shlex.quote(str(conda_base / 'bin'))}:{shlex.quote(str(conda_base / 'condabin'))}:$PATH",
                "unset PYTHONTZPATH",
            ]
        )
    lines.extend(
        [
            f"echo '[low-bm] starting {spec.trial_id} at '$(date -Is)",
            "echo '[low-bm] command follows:'",
            f"printf '%s\\n' {shlex.quote(shlex.join(command))}",
            shlex.join(command),
            f"echo '[low-bm] completed {spec.trial_id} at '$(date -Is)",
            "",
        ]
    )
    script.write_text("\n".join(lines))
    script.chmod(0o755)
    return script


def provenance_runner_bin_dir(args: argparse.Namespace) -> str | None:
    """Return the direct runner bin directory recorded for processing runs."""
    if getattr(args, "activate_command", ""):
        return None
    prefix = repo_path(Path(getattr(args, "runner_prefix", DEFAULT_RUNNER_PREFIX))).resolve()
    return str(runner_bin_dir(prefix))


def write_provenance(spec: RunSpec, args: argparse.Namespace, execution_command: list[str]) -> None:
    """Write enough launch metadata to reconstruct what was submitted."""
    snakemake_command = build_snakemake_command(spec)
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "low_bm_version": __version__,
        "trial_id": spec.trial_id,
        "trial_name": spec.trial_name,
        "mode": spec.mode,
        "target": spec.target,
        "profile": str(spec.profile),
        "configfiles": [str(path) for path in spec.configfiles],
        "isolated_workdir": bool(spec.workdir),
        "snakemake_workdir": str(spec.workdir) if spec.workdir else None,
        "snakefile": str(spec.snakefile) if spec.snakefile else None,
        "snakemake_conda_prefix": (
            str(spec.snakemake_conda_prefix) if spec.snakemake_conda_prefix else None
        ),
        "snakemake_command": snakemake_command,
        "execution_command": execution_command,
        "runner": getattr(args, "runner", "host"),
        "runner_prefix": getattr(args, "runner_prefix", None),
        "runner_bin_dir": provenance_runner_bin_dir(args),
        "manager": getattr(args, "manager", None),
        "git_commit": git_commit(),
        "activate_command_set": bool(getattr(args, "activate_command", "")),
    }
    (spec.log_dir / "run_metadata.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (spec.log_dir / "snakemake_command.txt").write_text(shlex.join(snakemake_command) + "\n")
    (spec.log_dir / "execution_command.txt").write_text(shlex.join(execution_command) + "\n")


def git_commit() -> str | None:
    """Return the current git commit when the checkout is available."""
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()


if __name__ == "__main__":
    raise SystemExit(main())
