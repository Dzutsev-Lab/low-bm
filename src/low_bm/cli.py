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
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .batch import (
    BatchRow,
    load_simple_run_config,
    read_batch_table,
    write_run_config,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BASE_CONFIG = "config/local/processing.yaml"
DEFAULT_BATCH_TABLE = "config/local/batch.tsv"
TEMPLATE_BASE_CONFIG = "config/templates/processing.yaml"
TEMPLATE_BATCH_TABLE = "config/templates/batch.tsv"
DEFAULT_RUN_CONFIG_DIR = "experiment_batch_configs"
DEFAULT_LOG_ROOT = "snakemake_logs"
DEFAULT_RUNNER_PREFIX = ".low-bm/runner/env"
DEFAULT_RUNNER_ENV_FILE = "workflow/envs/runner-env.yaml"
ENV_MANAGERS = ("mamba", "conda", "micromamba")
GIT_SHA_RE = re.compile(r"^[0-9a-fA-F]{7,40}$")


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


def add_common_config_arguments(parser: argparse.ArgumentParser) -> None:
    """Add config-layering options shared by single-row and table runs."""
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
    parser.add_argument("--unlock", action="store_true", help="Unlock the Snakemake working directory.")
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

    exit_code = 0
    for row in rows:
        row_config = write_run_config(row, args.run_config_dir)
        # Reuse the same run-building path as `low-bm run` so batch-table and
        # direct one-batch launches cannot drift into different behaviors.
        row_args = argparse.Namespace(**vars(args))
        row_args.row_config = str(row_config)
        row_args.trial_id = None
        row_args.trial_descript = None
        row_args.exp_dir = None
        row_args.metadata = None
        row_args.process_umis = None
        row_args.unlock = False
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


def setup_runner_command(args: argparse.Namespace) -> int:
    """Create or validate the project-local Snakemake runner environment."""
    prefix = Path(args.runner_prefix)
    env_file = REPO_ROOT / DEFAULT_RUNNER_ENV_FILE
    if not env_file.exists():
        raise SystemExit(f"Missing runner environment file: {DEFAULT_RUNNER_ENV_FILE}")

    manager = resolve_env_manager(args.manager, metadata_path=runner_metadata_path(prefix))
    if path_exists_for_repo_run(prefix):
        if not args.force:
            runner = HostRunnerSpec(
                prefix=prefix,
                manager=manager,
                metadata_path=runner_metadata_path(prefix),
            )
            print(f"Runner env already exists: {prefix}")
            return check_runner_command(runner, ["snakemake", "--version"], "snakemake")
        remove_runner_prefix(prefix)

    command = build_runner_create_command(prefix, env_file, manager)
    print(f"Creating runner env at {prefix} with {manager.name}.")
    completed = subprocess.run(command, cwd=REPO_ROOT, check=False)
    if completed.returncode != 0:
        return completed.returncode

    metadata_path = runner_metadata_path(prefix)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "runner_prefix": str(prefix),
        "env_file": DEFAULT_RUNNER_ENV_FILE,
        "env_file_sha256": hash_file(env_file),
        "manager": manager.name,
        "manager_executable": manager.executable,
        "low_bm_version": __version__,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote runner metadata: {metadata_path}")
    return check_runner_command(
        HostRunnerSpec(prefix=prefix, manager=manager, metadata_path=metadata_path),
        ["snakemake", "--version"],
        "snakemake",
    )


def doctor_runner_command(args: argparse.Namespace) -> int:
    """Check whether the host runner is ready for the selected execution mode."""
    runner = resolve_host_runner(args)
    checks: list[tuple[str, bool, str]] = []

    snakemake = run_runner_command(runner, ["snakemake", "--version"])
    checks.append(("snakemake", snakemake.returncode == 0, runner_check_detail(snakemake)))

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
    }
    missing = [name.replace("_", "-") for name, value in required.items() if not value]
    if missing:
        raise SystemExit(
            "Provide --row-config or direct batch fields: "
            + ", ".join(f"--{name}" for name in missing)
        )
    row = BatchRow(
        row_number=1,
        trialID=args.trial_id,
        trial_descript=args.trial_descript,
        exp_dir=args.exp_dir,
        metadata=args.metadata,
        process_umis=args.process_umis,
    )
    return write_run_config(row, args.run_config_dir)


def build_run_spec(args: argparse.Namespace, row_config: Path) -> RunSpec:
    """Resolve defaults and config layering for one Snakemake invocation."""
    if not row_config.exists():
        raise SystemExit(f"Missing row config: {row_config}")
    row_values = load_simple_run_config(row_config)
    trial_id = row_values.get("trialID") or args.trial_id or "run"
    profile = Path(args.profile) if args.profile else REPO_ROOT / "profiles" / args.mode
    # Config order is deliberate: base project defaults, optional user
    # overrides, then the generated per-row values for this batch.
    configfiles = resolve_configfiles(args.configfile, args.extra_configfile)
    configfiles.append(row_config)
    log_dir = Path(args.log_root) / trial_id
    job_name = args.job_name or f"low-bm-{trial_id}"
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
    )


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


def build_runner_create_command(prefix: Path, env_file: Path, manager: EnvManagerSpec) -> list[str]:
    """Build the command that creates the project-local runner environment."""
    if manager.name == "micromamba":
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


def remove_runner_prefix(prefix: Path) -> None:
    """Remove an existing runner prefix after a few guardrails."""
    target = repo_path(prefix).resolve()
    blocked = {REPO_ROOT.resolve(), Path.home().resolve(), Path("/").resolve()}
    if target in blocked or ".low-bm" not in target.parts:
        raise SystemExit(f"Refusing to remove unsafe runner prefix: {target}")
    shutil.rmtree(target)


def resolve_host_runner(args: argparse.Namespace) -> HostRunnerSpec:
    """Resolve the V1 host runner.

    A prefix path is more reproducible than a global env name because it is tied
    to this checkout. Two projects can carry different runner prefixes without
    fighting over a shared `low-bm-runner` name in the user's conda install.
    """
    if getattr(args, "runner", "host") != "host":
        raise SystemExit("V1 only supports --runner host. Container runner support will be added later.")
    prefix = Path(getattr(args, "runner_prefix", DEFAULT_RUNNER_PREFIX))
    metadata_path = runner_metadata_path(prefix)
    if not repo_path(prefix).is_dir():
        raise SystemExit(
            f"Missing runner env: {prefix}. "
            "Run `low-bm setup runner` before running the workflow, "
            "or use --activate-command for an advanced submitted-SLURM fallback."
        )
    manager = resolve_env_manager(getattr(args, "manager", "auto"), metadata_path=metadata_path)
    return HostRunnerSpec(prefix=prefix, manager=manager, metadata_path=metadata_path)


def wrap_runner_command(runner: HostRunnerSpec, command: list[str]) -> list[str]:
    """Run a command inside the project-local runner prefix without shell activation.

    `conda run -p` style invocation is a portability bridge: it avoids relying on
    interactive shell startup files, `source activate`, or a site-specific module
    being reloaded inside every generated master job.
    """
    return [
        runner.manager.executable,
        "run",
        "--prefix",
        str(runner.prefix),
        *command,
    ]


def wrap_snakemake_command(runner: HostRunnerSpec, command: list[str]) -> list[str]:
    """Wrap a Snakemake argv list in the active runner implementation."""
    return wrap_runner_command(runner, command)


def run_runner_command(runner: HostRunnerSpec, command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        wrap_runner_command(runner, command),
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


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
    command = ["snakemake", "--profile", str(spec.profile)]
    if spec.dry_run:
        command.append("--dry-run")
    if spec.unlock:
        command.append("--unlock")
    if not spec.unlock and spec.target:
        # Snakemake 9 parses --configfile as FILE [FILE ...]. Keep positional
        # targets before configfile entries so targets like "all" are not
        # interpreted as additional config files.
        command.append(spec.target)
    for configfile in spec.configfiles:
        command.extend(["--configfile", str(configfile)])
    command.extend(spec.extra_snakemake_args)
    return command


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
    command = build_execution_command(spec, args)
    write_provenance(spec, args, command)

    if args.print_command:
        print(shlex.join(command))

    # Dry-runs and unlocks should happen immediately in the current shell; they
    # are Snakemake control operations, not cluster work to submit.
    if spec.mode == "local" or spec.dry_run or spec.unlock:
        return run_local_command(command, spec)
    return submit_master_job(command, spec, args)


def run_local_command(command: list[str], spec: RunSpec) -> int:
    """Execute Snakemake directly and tee all output into the run log."""
    log_file = spec.log_dir / ("snakemake.dryrun.log" if spec.dry_run else "snakemake.log")
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
        prefix = Path(args.runner_prefix)
        # A submitted master job starts later, possibly on another node. This
        # guard catches deleted or unmounted project-local runner envs before a
        # confusing Snakemake failure scrolls by.
        missing_message = f"Missing runner env: {prefix}. Run low-bm setup runner first."
        lines.extend(
            [
                f"if [[ ! -d {shlex.quote(str(prefix))} ]]; then",
                f"  echo {shlex.quote(missing_message)} >&2",
                "  exit 2",
                "fi",
            ]
        )
    lines.extend(
        [
            f"echo '[low-bm] starting {spec.trial_id} at '$(date -Is)",
            shlex.join(command),
            f"echo '[low-bm] completed {spec.trial_id} at '$(date -Is)",
            "",
        ]
    )
    script.write_text("\n".join(lines))
    script.chmod(0o755)
    return script


def write_provenance(spec: RunSpec, args: argparse.Namespace, execution_command: list[str]) -> None:
    """Write enough launch metadata to reconstruct what was submitted."""
    snakemake_command = build_snakemake_command(spec)
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "low_bm_version": __version__,
        "trial_id": spec.trial_id,
        "mode": spec.mode,
        "target": spec.target,
        "profile": str(spec.profile),
        "configfiles": [str(path) for path in spec.configfiles],
        "snakemake_command": snakemake_command,
        "execution_command": execution_command,
        "runner": getattr(args, "runner", "host"),
        "runner_prefix": getattr(args, "runner_prefix", None),
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
