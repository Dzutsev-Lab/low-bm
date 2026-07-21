"""Command-line entry point for low-bm workflow execution.

This module is the small orchestration layer between a human-friendly run
command and Snakemake itself. It does not know the workflow DAG; it prepares the
right config files, selects a Snakemake profile, records provenance, and either
runs Snakemake directly or submits a lightweight SLURM master job.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
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
DEFAULT_BASE_CONFIG = "config.yaml"
DEFAULT_RUN_CONFIG_DIR = "experiment_batch_configs"
DEFAULT_LOG_ROOT = "snakemake_logs"


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

    return parser


def add_common_config_arguments(parser: argparse.ArgumentParser) -> None:
    """Add config-layering options shared by single-row and table runs."""
    # Snakemake applies later --configfile entries as overrides, so the row
    # config is appended after the base config in build_run_spec().
    parser.add_argument(
        "--configfile",
        action="append",
        default=[],
        help="Base Snakemake config file. May be repeated. Defaults to config.yaml.",
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
    parser.add_argument(
        "--batch-table",
        default="experiment_batch_configs.tsv",
        help="Canonical batch TSV.",
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
            "Shell command evaluated in the master job before Snakemake, "
            "for example: 'source myconda && mamba activate low-bm-runner'."
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


def run_command(args: argparse.Namespace) -> int:
    """Run or submit one batch after resolving its per-row config file."""
    row_config = resolve_row_config(args)
    spec = build_run_spec(args, row_config)
    return execute_run_spec(spec, args)


def batch_submit_command(args: argparse.Namespace) -> int:
    """Submit one independent run for every row in the batch table."""
    rows = read_batch_table(args.batch_table)
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
            command = build_snakemake_command(spec)
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
    row_values = load_simple_run_config(row_config)
    trial_id = row_values.get("trialID") or args.trial_id or "run"
    profile = Path(args.profile) if args.profile else REPO_ROOT / "profiles" / args.mode
    # Config order is deliberate: base project defaults, optional user
    # overrides, then the generated per-row values for this batch.
    configfiles = [Path(p) for p in args.configfile] or [Path(DEFAULT_BASE_CONFIG)]
    configfiles.extend(Path(p) for p in args.extra_configfile)
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


def execute_run_spec(spec: RunSpec, args: argparse.Namespace) -> int:
    """Record provenance, then run locally or submit a SLURM master job."""
    spec.log_dir.mkdir(parents=True, exist_ok=True)
    write_provenance(spec, args)
    command = build_snakemake_command(spec)

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
    if args.activate_command:
        # Keep environment activation configurable because HPC shell setup tends
        # to be site- and account-specific.
        lines.append(args.activate_command)
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


def write_provenance(spec: RunSpec, args: argparse.Namespace) -> None:
    """Write enough launch metadata to reconstruct what was submitted."""
    command = build_snakemake_command(spec)
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "low_bm_version": __version__,
        "trial_id": spec.trial_id,
        "mode": spec.mode,
        "target": spec.target,
        "profile": str(spec.profile),
        "configfiles": [str(path) for path in spec.configfiles],
        "snakemake_command": command,
        "git_commit": git_commit(),
        "activate_command_set": bool(args.activate_command),
    }
    (spec.log_dir / "run_metadata.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (spec.log_dir / "snakemake_command.txt").write_text(shlex.join(command) + "\n")


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
