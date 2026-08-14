#!/usr/bin/env python3

from __future__ import annotations

import io
import os
import shlex
import tempfile
import unittest
import json
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from low_bm.batch import read_batch_table, write_run_config
import low_bm.cli as cli
from low_bm.cli import (
    EnvManagerSpec,
    batch_submit_command,
    build_execution_command,
    build_run_spec,
    build_runner_create_command,
    build_snakemake_command,
    references_check_command,
    references_prepare_bwa_indexes_command,
    resolve_env_manager,
    validate_microclean_source,
    write_master_script,
)


class BatchParsingTests(unittest.TestCase):
    def test_headered_table_parses_all_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "\t".join(
                    [
                        "trialID",
                        "trial_descript",
                        "exp_dir",
                        "metadata",
                        "host",
                        "process_umis",
                    ]
                )
                + "\n"
                + "\t".join(
                    [
                        "010126.1",
                        "batch1",
                        "Batch1",
                        "metadata/batch1.xlsx",
                        "Human",
                        "false",
                    ]
                )
                + "\n"
                + "\t".join(
                    [
                        "010126.2",
                        "batch2",
                        "Batch2",
                        "metadata/batch2.xlsx",
                        "mouse",
                        "",
                    ]
                )
                + "\n"
            )

            rows = read_batch_table(table)
            self.assertEqual([row.trialID for row in rows], ["010126.1", "010126.2"])
            self.assertEqual(rows[0].trial_name, "010126.1_batch1")
            self.assertEqual([row.host for row in rows], ["human", "mouse"])
            self.assertEqual(rows[0].process_umis, "false")

    def test_unheadered_table_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text("010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\n")
            with self.assertRaisesRegex(ValueError, "Missing required batch table column"):
                read_batch_table(table)

    def test_run_config_quotes_scalars(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\tprocess_umis\n"
                "010126.1\tbatch with spaces\tBatch1\tmetadata/batch 1.xlsx\thuman\tfalse\n"
            )
            row = read_batch_table(table)[0]
            run_config = write_run_config(row, Path(tmp) / "configs")
            text = run_config.read_text()
            self.assertIn('trialID: "010126.1"', text)
            self.assertIn('trial_descript: "batch with spaces"', text)
            self.assertIn('host: "human"', text)
            self.assertIn('process_umis: "false"', text)
            self.assertNotIn("batch_label", text)
            self.assertNotIn("include_processing", text)
            self.assertNotIn("include_analysis", text)

    def test_batch_table_requires_host_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\tprocess_umis\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\tfalse\n"
            )
            with self.assertRaisesRegex(ValueError, "Missing required batch table column.*host"):
                read_batch_table(table)

    def test_batch_table_rejects_invalid_host(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\trat\n"
            )
            with self.assertRaisesRegex(ValueError, "Invalid host value"):
                read_batch_table(table)


class CommandConstructionTests(unittest.TestCase):
    def test_snakemake_command_preserves_config_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\nhost: "human"\n')
            base_config = Path(tmp) / "processing.yaml"
            write_minimal_processing_config(base_config)
            override_config = Path(tmp) / "override.yaml"
            override_config.write_text("out_root: Exp_Output_validation/\n")
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[str(override_config)],
                row_config=str(run_config),
                trial_id=None,
                mode="local",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                job_name=None,
                dry_run=True,
                unlock=False,
                snakemake_arg=["--quiet"],
            )
            spec = build_run_spec(args, run_config)
            command = build_snakemake_command(spec)
            self.assertEqual(command[0:2], ["snakemake", "--profile"])
            self.assertIn("--dry-run", command)
            self.assertEqual(command.count("--configfile"), 1)
            self.assertEqual(
                command[command.index("--configfile") + 1 : command.index("--")],
                [str(base_config), str(override_config), str(run_config)],
            )
            self.assertLess(command.index("--quiet"), command.index("--configfile"))
            self.assertEqual(command[-2:], ["--", "all"])

    def test_processing_config_stack_requires_base_layer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\nhost: "human"\n')
            override_only = Path(tmp) / "processing-overrides.yaml"
            override_only.write_text("out_root: Exp_Output_validation/\n")
            args = Args(
                configfile=[str(override_only)],
                extra_configfile=[],
                row_config=str(run_config),
                trial_id=None,
                mode="local",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                job_name=None,
                dry_run=True,
                unlock=False,
                snakemake_arg=[],
            )
            with self.assertRaisesRegex(SystemExit, "config/local/processing.yaml"):
                build_run_spec(args, run_config)

    def test_row_config_requires_host_even_when_processing_config_has_host(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\n')
            base_config = Path(tmp) / "processing.yaml"
            write_processing_config_with_refs(base_config, Path(tmp) / "refs")
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[],
                row_config=str(run_config),
                trial_id=None,
                mode="local",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                job_name=None,
                dry_run=True,
                unlock=False,
                snakemake_arg=[],
            )
            with self.assertRaisesRegex(SystemExit, "missing required key 'host'"):
                build_run_spec(args, run_config)

    def test_direct_run_writes_required_host(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = Args(
                row_config=None,
                trial_id="010126.1",
                trial_descript="batch1",
                exp_dir="Batch1",
                metadata="metadata/batch1.xlsx",
                host="mouse",
                process_umis=None,
                run_config_dir=str(Path(tmp) / "configs"),
            )
            run_config = cli.resolve_row_config(args)
            self.assertIn('host: "mouse"', run_config.read_text())

    def test_direct_run_requires_host(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = Args(
                row_config=None,
                trial_id="010126.1",
                trial_descript="batch1",
                exp_dir="Batch1",
                metadata="metadata/batch1.xlsx",
                host=None,
                process_umis=None,
                run_config_dir=str(Path(tmp) / "configs"),
            )
            with self.assertRaisesRegex(SystemExit, "--host"):
                cli.resolve_row_config(args)

    def test_direct_run_host_choices_are_limited(self) -> None:
        with self.assertRaises(SystemExit):
            cli.build_parser().parse_args(
                [
                    "run",
                    "--trial-id",
                    "010126.1",
                    "--trial-descript",
                    "batch1",
                    "--exp-dir",
                    "Batch1",
                    "--metadata",
                    "metadata/batch1.xlsx",
                    "--host",
                    "rat",
                ]
            )

    def test_batch_submit_dry_run_reuses_run_builder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\thuman\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            for ref_name in ("human.fa", "viral.fa", "bacteria.fa"):
                write_bwa_sidecars(refs / ref_name)
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[str(base_config)],
                extra_configfile=[],
                mode="slurm",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                dry_run=True,
                print_command=False,
                snakemake_arg=[],
                activate_command="",
                master_cpus="4",
                master_mem="8G",
                master_time="1-00:00:00",
                master_partition=None,
                master_extra_sbatch=[],
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
            )
            self.assertEqual(batch_submit_command(args), 0)

    def test_batch_submit_dry_run_emits_unique_isolated_directories(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tbatch one\tBatch1\tmetadata/batch1.xlsx\thuman\n"
                "010126.2\tbatch two\tBatch2\tmetadata/batch2.xlsx\thuman\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            for ref_name in ("human.fa", "viral.fa", "bacteria.fa"):
                write_bwa_sidecars(refs / ref_name)
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            workdir_root = Path(tmp) / "workdirs"
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[str(base_config)],
                extra_configfile=[],
                mode="slurm",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                shared_workdir=False,
                workdir_root=str(workdir_root),
                snakemake_conda_prefix=str(Path(tmp) / "snakemake-conda"),
                dry_run=True,
                print_command=False,
                snakemake_arg=[],
                activate_command="",
                master_cpus="4",
                master_mem="8G",
                master_time="1-00:00:00",
                master_partition=None,
                master_extra_sbatch=[],
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
            )

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(batch_submit_command(args), 0)
            text = stdout.getvalue()
            self.assertIn("--directory", text)
            self.assertIn(str((workdir_root / "010126.1_batch_one").resolve()), text)
            self.assertIn(str((workdir_root / "010126.2_batch_two").resolve()), text)
            self.assertIn("--conda-prefix", text)

    def test_isolated_snakemake_command_uses_absolute_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\ntrial_descript: "batch with spaces"\nhost: "human"\n')
            base_config = Path(tmp) / "processing.yaml"
            write_minimal_processing_config(base_config)
            override_config = Path(tmp) / "override.yaml"
            override_config.write_text("out_root: Exp_Output_validation/\n")
            workdir_root = Path(tmp) / "workdirs"
            conda_prefix = Path(tmp) / "snakemake-conda"
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[str(override_config)],
                row_config=str(run_config),
                trial_id=None,
                trial_descript=None,
                mode="slurm",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                job_name=None,
                dry_run=True,
                unlock=False,
                snakemake_arg=[],
                isolated_workdir=True,
                workdir_root=str(workdir_root),
                snakemake_conda_prefix=str(conda_prefix),
            )
            spec = build_run_spec(args, run_config)
            command = build_snakemake_command(spec)

            self.assertEqual(command[0:2], ["snakemake", "--snakefile"])
            self.assertEqual(command[command.index("--snakefile") + 1], str((REPO_ROOT / "Snakefile").resolve()))
            self.assertEqual(command[command.index("--directory") + 1], str((workdir_root / "010126.1_batch_with_spaces").resolve()))
            self.assertEqual(command[command.index("--profile") + 1], str((REPO_ROOT / "profiles" / "slurm").resolve()))
            self.assertEqual(command[command.index("--conda-prefix") + 1], str(conda_prefix.resolve()))
            config_paths = command[command.index("--configfile") + 1 : command.index("--")]
            self.assertEqual(
                config_paths,
                [str(base_config.resolve()), str(override_config.resolve()), str(run_config.resolve())],
            )

    def test_shared_workdir_preserves_legacy_command_shape(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\ntrial_descript: "batch1"\nhost: "human"\n')
            base_config = Path(tmp) / "processing.yaml"
            write_minimal_processing_config(base_config)
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[],
                row_config=str(run_config),
                trial_id=None,
                trial_descript=None,
                mode="slurm",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                job_name=None,
                dry_run=True,
                unlock=False,
                snakemake_arg=[],
                shared_workdir=True,
                isolated_workdir=True,
                workdir_root=str(Path(tmp) / "workdirs"),
                snakemake_conda_prefix=str(Path(tmp) / "snakemake-conda"),
            )
            command = build_snakemake_command(build_run_spec(args, run_config))
            self.assertEqual(command[0:2], ["snakemake", "--profile"])
            self.assertNotIn("--snakefile", command)
            self.assertNotIn("--directory", command)
            self.assertNotIn("--conda-prefix", command)

    def test_batch_unlock_dry_run_builds_local_unlock_without_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tbatch one\tBatch1\tmetadata/batch1.xlsx\thuman\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            write_minimal_processing_config(base_config)
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            workdir_root = Path(tmp) / "workdirs"
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[str(base_config)],
                extra_configfile=[],
                trial_id=["010126.1"],
                all=False,
                mode="slurm",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                shared_workdir=False,
                workdir_root=str(workdir_root),
                snakemake_conda_prefix=str(Path(tmp) / "snakemake-conda"),
                dry_run=True,
                print_command=False,
                snakemake_arg=[],
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
            )

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(cli.batch_unlock_command(args), 0)
            text = stdout.getvalue()
            self.assertIn("--unlock", text)
            self.assertIn("--directory", text)
            self.assertIn(str((workdir_root / "010126.1_batch_one").resolve()), text)
            self.assertNotIn("-- all", text)

    def test_batch_submit_rejects_duplicate_output_trial_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\thuman\n"
                "010126.1\tbatch1\tBatch2\tmetadata/batch2.xlsx\thuman\n"
            )
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[],
                extra_configfile=[],
                mode="slurm",
                shared_workdir=False,
            )
            with self.assertRaisesRegex(SystemExit, "Duplicate batch output trial name"):
                batch_submit_command(args)

    def test_isolated_batch_submit_rejects_missing_shared_bwa_indexes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\thuman\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[str(base_config)],
                extra_configfile=[],
                mode="slurm",
                shared_workdir=False,
            )

            with self.assertRaisesRegex(SystemExit, "Missing BWA reference index files"):
                batch_submit_command(args)

    def test_isolated_batch_submit_allows_existing_shared_bwa_indexes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\thuman\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            for ref_name in ("human.fa", "viral.fa", "bacteria.fa"):
                write_bwa_sidecars(refs / ref_name)
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[str(base_config)],
                extra_configfile=[],
                mode="slurm",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                shared_workdir=False,
                workdir_root=str(Path(tmp) / "workdirs"),
                snakemake_conda_prefix=str(Path(tmp) / "snakemake-conda"),
                dry_run=True,
                print_command=False,
                snakemake_arg=[],
                activate_command="",
                master_cpus="4",
                master_mem="8G",
                master_time="1-00:00:00",
                master_partition=None,
                master_extra_sbatch=[],
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
            )

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(batch_submit_command(args), 0)
            self.assertIn("would submit master job script", stdout.getvalue())

    def test_isolated_batch_submit_uses_row_host_for_bwa_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tmouse_batch\tBatch1\tmetadata/batch1.xlsx\tmouse\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            for ref_name in ("human.fa", "viral.fa", "bacteria.fa"):
                write_bwa_sidecars(refs / ref_name)
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[str(base_config)],
                extra_configfile=[],
                mode="slurm",
                shared_workdir=False,
            )

            with self.assertRaisesRegex(SystemExit, "mouse_ref: missing BWA sidecar"):
                batch_submit_command(args)

    def test_isolated_batch_submit_mixed_hosts_requires_both_host_indexes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\thuman_batch\tBatch1\tmetadata/batch1.xlsx\thuman\n"
                "010126.2\tmouse_batch\tBatch2\tmetadata/batch2.xlsx\tmouse\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            for ref_name in ("human.fa", "mouse.fa", "viral.fa", "bacteria.fa"):
                write_bwa_sidecars(refs / ref_name)
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[str(base_config)],
                extra_configfile=[],
                mode="slurm",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                shared_workdir=False,
                workdir_root=str(Path(tmp) / "workdirs"),
                snakemake_conda_prefix=str(Path(tmp) / "snakemake-conda"),
                dry_run=True,
                print_command=False,
                snakemake_arg=[],
                activate_command="",
                master_cpus="4",
                master_mem="8G",
                master_time="1-00:00:00",
                master_partition=None,
                master_extra_sbatch=[],
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
            )

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(batch_submit_command(args), 0)
            self.assertIn("[010126.1] would submit master job script", stdout.getvalue())
            self.assertIn("[010126.2] would submit master job script", stdout.getvalue())

    def test_shared_workdir_batch_submit_skips_shared_bwa_index_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\thuman\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            write_processing_config_with_refs(base_config, Path(tmp) / "refs")
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[str(base_config)],
                extra_configfile=[],
                mode="slurm",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                shared_workdir=True,
                workdir_root=str(Path(tmp) / "workdirs"),
                snakemake_conda_prefix=str(Path(tmp) / "snakemake-conda"),
                dry_run=True,
                print_command=False,
                snakemake_arg=[],
                activate_command="",
                master_cpus="4",
                master_mem="8G",
                master_time="1-00:00:00",
                master_partition=None,
                master_extra_sbatch=[],
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
            )

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(batch_submit_command(args), 0)
            text = stdout.getvalue()
            self.assertIn("would submit master job script", text)
            self.assertNotIn("--directory", text)

    def test_references_check_succeeds_with_sidecars_and_warns_without_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            for ref_name in ("human.fa", "viral.fa", "bacteria.fa"):
                write_bwa_sidecars(refs / ref_name)
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[],
                all_configured_hosts=False,
            )

            stdout = io.StringIO()
            stderr = io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(stderr):
                self.assertEqual(references_check_command(args), 0)

            self.assertIn("BWA reference check passed for 3 reference(s)", stdout.getvalue())
            self.assertIn("missing manifest", stderr.getvalue())

    def test_references_check_fails_for_missing_sidecars(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[],
                all_configured_hosts=False,
            )

            with self.assertRaisesRegex(SystemExit, "missing BWA sidecar"):
                references_check_command(args)

    def test_references_check_warns_on_stale_manifest_without_failing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            for ref_name in ("human.fa", "viral.fa", "bacteria.fa"):
                reference = refs / ref_name
                write_bwa_sidecars(reference)
                write_bwa_manifest(reference, reference_sha256="stale")
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[],
                all_configured_hosts=False,
            )

            stdout = io.StringIO()
            stderr = io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(stderr):
                self.assertEqual(references_check_command(args), 0)

            self.assertIn("manifest checksum does not match", stderr.getvalue())

    def test_prepare_bwa_indexes_dry_run_prints_missing_and_stale_refs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            for ref_name in ("human.fa", "viral.fa"):
                reference = refs / ref_name
                write_bwa_sidecars(reference)
                write_bwa_manifest(reference, reference_sha256="stale")
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[],
                all_configured_hosts=False,
                log_root=str(Path(tmp) / "reference_logs"),
                log_dir=None,
                dry_run=True,
                print_command=False,
                force=False,
                env_mode="direct",
                manager="auto",
                analysis_env_root=str(Path(tmp) / "analysis_envs"),
                bio_env_prefix=None,
            )

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(references_prepare_bwa_indexes_command(args), 0)
            text = stdout.getvalue()
            self.assertEqual(text.count("bwa index"), 3)
            self.assertIn(str(refs / "human.fa"), text)
            self.assertIn(str(refs / "viral.fa"), text)
            self.assertIn(str(refs / "bacteria.fa"), text)

    def test_prepare_bwa_indexes_skips_complete_refs_unless_forced(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base_config = Path(tmp) / "processing.yaml"
            refs = Path(tmp) / "refs"
            write_processing_config_with_refs(base_config, refs)
            for ref_name in ("human.fa", "viral.fa", "bacteria.fa"):
                reference = refs / ref_name
                write_bwa_sidecars(reference)
                write_bwa_manifest(reference)
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[],
                all_configured_hosts=False,
                log_root=str(Path(tmp) / "reference_logs"),
                log_dir=None,
                dry_run=True,
                print_command=False,
                force=False,
                env_mode="direct",
                manager="auto",
                analysis_env_root=str(Path(tmp) / "analysis_envs"),
                bio_env_prefix=None,
            )

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(references_prepare_bwa_indexes_command(args), 0)
            self.assertNotIn("bwa index", stdout.getvalue())

            args.force = True
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(references_prepare_bwa_indexes_command(args), 0)
            self.assertEqual(stdout.getvalue().count("bwa index"), 3)

    def test_execution_command_uses_direct_runner_snakemake(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\nhost: "human"\n')
            base_config = Path(tmp) / "processing.yaml"
            write_minimal_processing_config(base_config)
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[],
                row_config=str(run_config),
                trial_id=None,
                mode="local",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                job_name=None,
                dry_run=True,
                unlock=False,
                snakemake_arg=["--quiet"],
                activate_command="",
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
            )
            spec = build_run_spec(args, run_config)
            command = build_execution_command(spec, args)
            conda_base = Path(tmp) / "conda-base"
            self.assertEqual(command[0], str(runner_prefix.resolve() / "bin" / "snakemake"))
            self.assertNotIn("run", command[0:4])
            self.assertEqual(
                command[1:3],
                ["--conda-base-path", str(conda_base.resolve())],
            )
            self.assertEqual(command.count("--configfile"), 1)
            self.assertIn("--config", command)
            self.assertIn(f"low_bm_repo_root={REPO_ROOT.resolve()}", command)
            self.assertEqual(command[-2:], ["--", "all"])

    def test_snakemake_process_env_uses_conda_base_not_runner_bin(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            runner = cli.resolve_host_runner(
                Args(runner="host", runner_prefix=str(runner_prefix), manager="auto")
            )

            with patch.dict(
                os.environ,
                {
                    "PATH": "/usr/bin",
                    "PYTHONTZPATH": ".low-bm/runner/env/share/zoneinfo",
                },
                clear=True,
            ):
                env = cli.snakemake_process_env(runner)

            path_parts = env["PATH"].split(os.pathsep)
            self.assertEqual(path_parts[0], str(runner.conda_base_prefix / "bin"))
            self.assertEqual(path_parts[1], str(runner.conda_base_prefix / "condabin"))
            self.assertNotIn(str(runner.bin_dir), path_parts)
            self.assertNotIn("PYTHONTZPATH", env)
            self.assertNotIn("BASH_ENV", env)
            self.assertFalse(any(key.startswith("CONDA_") for key in env))

    def test_conda_base_prefix_override_wins(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            override = Path(tmp) / "override-conda"
            make_fake_conda_base(override)
            runner = cli.resolve_host_runner(
                Args(
                    runner="host",
                    runner_prefix=str(runner_prefix),
                    manager="auto",
                    conda_base_prefix=str(override),
                )
            )

            self.assertEqual(runner.conda_base_prefix, override.resolve())

    def test_missing_runner_prefix_has_setup_hint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\nhost: "human"\n')
            base_config = Path(tmp) / "processing.yaml"
            write_minimal_processing_config(base_config)
            args = Args(
                configfile=[str(base_config)],
                extra_configfile=[],
                row_config=str(run_config),
                trial_id=None,
                mode="local",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                job_name=None,
                dry_run=True,
                unlock=False,
                snakemake_arg=[],
                activate_command="",
                runner="host",
                runner_prefix=str(Path(tmp) / "missing" / "env"),
                manager="auto",
            )
            spec = build_run_spec(args, run_config)
            with self.assertRaisesRegex(SystemExit, "low-bm setup runner"):
                build_execution_command(spec, args)

    def test_master_script_uses_runner_without_activation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            log_dir = Path(tmp) / "logs"
            spec = cli.RunSpec(
                configfiles=[Path("processing.yaml"), Path("row.yaml")],
                target="all",
                profile=Path("profiles/slurm"),
                mode="slurm",
                log_dir=log_dir,
                trial_id="010126.1",
                job_name="low-bm-010126.1",
                dry_run=False,
                unlock=False,
                extra_snakemake_args=[],
            )
            runner_prefix, manager = write_fake_runner(Path(tmp))
            conda_base = Path(tmp) / "conda-base"
            command = [
                str(runner_prefix.resolve() / "bin" / "snakemake"),
                "--conda-base-path",
                str(conda_base.resolve()),
                "all",
            ]
            args = Args(
                activate_command="",
                master_cpus="4",
                master_mem="8G",
                master_time="1-00:00:00",
                master_partition=None,
                master_extra_sbatch=[],
                runner_prefix=str(runner_prefix),
            )
            script = write_master_script(command, spec, args)
            text = script.read_text()
            self.assertNotIn(str(manager), text)
            self.assertIn(str(runner_prefix.resolve()), text)
            conda_base = (Path(tmp) / "conda-base").resolve()
            self.assertIn(f"export PATH={conda_base / 'bin'}:{conda_base / 'condabin'}:$PATH", text)
            self.assertIn("for var_name in ${!CONDA_@}; do unset \"$var_name\"; done", text)
            self.assertNotIn(f"export PATH={runner_prefix.resolve() / 'bin'}:$PATH", text)
            self.assertNotIn("BASH_ENV", text)
            self.assertNotIn("conda-base-shim", text)
            self.assertIn("unset PYTHONTZPATH", text)
            self.assertIn("Missing runner env", text)
            self.assertNotIn("activate low-bm-runner", text)
            self.assertFalse((runner_prefix / "bin" / "activate").exists())

    def test_batch_prepare_envs_dry_run_builds_serial_local_commands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\thost\n"
                "010126.1\tbatch one\tBatch1\tmetadata/batch1.xlsx\thuman\n"
                "010126.2\tbatch two\tBatch2\tmetadata/batch2.xlsx\thuman\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            write_minimal_processing_config(base_config)
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[str(base_config)],
                extra_configfile=[],
                trial_id=[],
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                shared_workdir=False,
                workdir_root=str(Path(tmp) / "workdirs"),
                snakemake_conda_prefix=str(Path(tmp) / "snakemake-conda"),
                dry_run=True,
                print_command=False,
                snakemake_arg=[],
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
            )

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(cli.batch_prepare_envs_command(args), 0)
            text = stdout.getvalue()
            self.assertIn("[010126.1] would prepare envs:", text)
            self.assertIn("[010126.2] would prepare envs:", text)
            self.assertEqual(text.count("--conda-create-envs-only"), 2)
            self.assertEqual(text.count(str(runner_prefix.resolve() / "bin" / "snakemake")), 2)
            self.assertNotIn("sbatch", text)

    def test_snakefile_uses_launcher_repo_root_config(self) -> None:
        text = (REPO_ROOT / "Snakefile").read_text()
        self.assertIn('config.get("low_bm_repo_root"', text)
        self.assertNotIn("workflow.source_path", text)

    def test_missing_default_processing_config_has_setup_hint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\nhost: "human"\n')
            missing_default = Path(tmp) / "config" / "local" / "processing.yaml"
            args = Args(
                configfile=[],
                extra_configfile=[],
                row_config=str(run_config),
                trial_id=None,
                mode="local",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                job_name=None,
                dry_run=True,
                unlock=False,
                snakemake_arg=[],
            )
            old_default = cli.DEFAULT_BASE_CONFIG
            try:
                cli.DEFAULT_BASE_CONFIG = str(missing_default)
                with self.assertRaisesRegex(
                    SystemExit,
                    r"Missing default processing config: .*config.*/local.*/processing\.yaml.*config/templates/processing\.yaml",
                ):
                    build_run_spec(args, run_config)
            finally:
                cli.DEFAULT_BASE_CONFIG = old_default

    def test_missing_default_batch_table_has_setup_hint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            missing_default = Path(tmp) / "config" / "local" / "batch.tsv"
            args = Args(
                batch_table=str(missing_default),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=[],
                extra_configfile=[],
                mode="slurm",
                profile=None,
                target="all",
                log_root=str(Path(tmp) / "logs"),
                dry_run=True,
                print_command=False,
                snakemake_arg=[],
                activate_command="",
                master_cpus="4",
                master_mem="8G",
                master_time="1-00:00:00",
                master_partition=None,
                master_extra_sbatch=[],
            )
            old_default = cli.DEFAULT_BATCH_TABLE
            try:
                cli.DEFAULT_BATCH_TABLE = str(missing_default)
                with self.assertRaisesRegex(
                    SystemExit,
                    r"Missing default batch table: .*config.*/local.*/batch\.tsv.*config/templates/batch\.tsv",
                ):
                    batch_submit_command(args)
            finally:
                cli.DEFAULT_BATCH_TABLE = old_default


class AnalysisCommandTests(unittest.TestCase):
    def test_compile_phyloseq_is_meta_not_standard_analysis(self) -> None:
        self.assertNotIn("compile-phyloseq", cli.ANALYSIS_STEP_REGISTRY)
        self.assertNotIn("compile-phyloseq", cli.ANALYSIS_STEP_CHOICES)
        self.assertIn("compile-phyloseq", cli.META_STEP_REGISTRY)
        self.assertIn("decontaminate-phyloseq", cli.META_STEP_REGISTRY)

    def test_analysis_run_defaults_to_managed_envs(self) -> None:
        args = cli.build_parser().parse_args(
            [
                "analysis",
                "run",
                "abundance-barplots",
                "--analysis-config",
                "config/local/analysis.yaml",
            ]
        )
        self.assertEqual(args.env_mode, "managed")
        self.assertEqual(args.analysis_env_root, cli.DEFAULT_ANALYSIS_ENV_ROOT)

    def test_analysis_command_construction_covers_registered_steps(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            analysis_config = Path(tmp) / "analysis.yaml"
            write_minimal_analysis_config(analysis_config)
            args = Args(
                analysis_config=str(analysis_config),
                steps=list(cli.ANALYSIS_STEP_REGISTRY),
                log_root=str(Path(tmp) / "analysis_logs"),
                log_dir=None,
                dry_run=True,
                env_mode="direct",
                manager="auto",
            )
            spec = cli.build_analysis_run_spec(
                args,
                created_utc=datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc),
            )
            self.assertEqual(spec.expanded_steps, list(cli.ANALYSIS_STEP_REGISTRY))
            self.assertEqual(spec.log_dir, Path(tmp) / "analysis_logs" / ("20260102T030405Z_" + "-".join(spec.expanded_steps)))
            for command_spec in spec.commands:
                self.assertEqual(command_spec.command, command_spec.execution_command)
                self.assertIn("--analysis-config", command_spec.command)
                self.assertEqual(command_spec.command[-1], str(analysis_config))

    def test_analysis_managed_env_wraps_r_steps_with_project_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            analysis_config = Path(tmp) / "analysis.yaml"
            write_minimal_analysis_config(analysis_config)
            env_root = Path(tmp) / "analysis_envs"
            args = Args(
                analysis_config=str(analysis_config),
                steps=["abundance-barplots"],
                log_root=str(Path(tmp) / "analysis_logs"),
                log_dir=None,
                dry_run=True,
                env_mode="managed",
                manager="/bin/echo",
                analysis_env_root=str(env_root),
                r_env_prefix=None,
                bio_env_prefix=None,
            )
            spec = cli.build_analysis_run_spec(
                args,
                created_utc=datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc),
            )
            expected_prefix = cli.managed_analysis_env_prefix(cli.R_TOOLS_ENV, env_root.resolve())
            self.assertEqual(spec.commands[0].env_prefix, expected_prefix)
            self.assertEqual(
                spec.commands[0].execution_command,
                [
                    "/bin/echo",
                    "run",
                    "--prefix",
                    str(expected_prefix),
                    "Rscript",
                    "scripts/AbundanceBarPlots.R",
                    "--analysis-config",
                    str(analysis_config),
                ],
            )

    def test_analysis_prefix_env_wraps_r_steps_with_prefix(self) -> None:
        analysis_config = Path("config/local/analysis.yaml")
        r_prefix = Path("/data/taylorng/conda/envs/low-bm-r-tools")
        command_spec = cli.build_analysis_command_spec(
            step="abundance-barplots",
            analysis_config=analysis_config,
            env_mode="prefix",
            manager=EnvManagerSpec(name="mamba", executable="/data/taylorng/conda/bin/mamba"),
            env_prefixes={cli.R_TOOLS_ENV: r_prefix},
        )
        self.assertEqual(command_spec.env_prefix, r_prefix)
        self.assertEqual(
            command_spec.execution_command,
            [
                "/data/taylorng/conda/bin/mamba",
                "run",
                "--prefix",
                str(r_prefix),
                "Rscript",
                "scripts/AbundanceBarPlots.R",
                "--analysis-config",
                str(analysis_config),
            ],
        )

    def test_analysis_prefix_env_requires_prefix_for_selected_env(self) -> None:
        args = Args(r_env_prefix=None, bio_env_prefix=None)
        with self.assertRaisesRegex(SystemExit, "--r-env-prefix"):
            cli.resolve_analysis_env_prefixes(
                args,
                ["abundance-barplots"],
                cli.ANALYSIS_STEP_REGISTRY,
                validate_exists=False,
            )

    def test_analysis_prefix_env_resolves_blast_confirmation_envs(self) -> None:
        args = Args(
            r_env_prefix="/data/taylorng/conda/envs/low-bm-r-tools",
            bio_env_prefix="/data/taylorng/conda/envs/low-bm-bio-tools",
            microclean_env_prefix=None,
        )
        prefixes = cli.resolve_analysis_env_prefixes(
            args,
            ["blast-candidates", "blast-search", "blast-plots"],
            cli.ANALYSIS_STEP_REGISTRY,
            validate_exists=False,
        )
        self.assertEqual(prefixes[cli.R_TOOLS_ENV], Path(args.r_env_prefix))
        self.assertEqual(prefixes[cli.BIO_TOOLS_ENV], Path(args.bio_env_prefix))

    def test_analysis_composite_expands_in_order(self) -> None:
        self.assertEqual(
            cli.expand_analysis_steps(["differential-abundance", "blast-confirmation", "survival"]),
            [
                "differential-abundance",
                "blast-candidates",
                "blast-search",
                "blast-plots",
                "survival",
            ],
        )

    def test_analysis_run_dry_run_writes_provenance_without_execution(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            analysis_config = Path(tmp) / "analysis.yaml"
            write_minimal_analysis_config(analysis_config)
            log_dir = Path(tmp) / "analysis_run"
            args = Args(
                analysis_config=str(analysis_config),
                steps=["abundance-barplots", "ordination"],
                log_root=str(Path(tmp) / "analysis_logs"),
                log_dir=str(log_dir),
                dry_run=True,
                print_command=False,
                env_mode="direct",
                manager="auto",
            )
            self.assertEqual(cli.analysis_run_command(args), 0)
            metadata = json.loads((log_dir / "analysis_metadata.json").read_text())
            self.assertEqual(metadata["expanded_steps"], ["abundance-barplots", "ordination"])
            self.assertTrue(metadata["dry_run"])
            commands_text = (log_dir / "analysis_commands.txt").read_text()
            self.assertLess(commands_text.index("# abundance-barplots"), commands_text.index("# ordination"))
            self.assertFalse((log_dir / "01_abundance-barplots.log").exists())

    def test_analysis_missing_config_has_setup_hint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            missing_config = Path(tmp) / "config" / "local" / "analysis.yaml"
            old_default = cli.DEFAULT_ANALYSIS_CONFIG
            try:
                cli.DEFAULT_ANALYSIS_CONFIG = str(missing_config)
                with self.assertRaisesRegex(SystemExit, "low-bm analysis init"):
                    cli.resolve_analysis_config(missing_config)
            finally:
                cli.DEFAULT_ANALYSIS_CONFIG = old_default

    def test_analysis_missing_selected_section_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            analysis_config = Path(tmp) / "analysis.yaml"
            analysis_config.write_text("project:\n  trialID: test\n  compiled_physeq: physeq.RData\n")
            with self.assertRaisesRegex(SystemExit, "ordination"):
                cli.validate_analysis_config(analysis_config, ["ordination"])

    def test_analysis_requires_compiled_physeq_endpoint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            analysis_config = Path(tmp) / "analysis.yaml"
            analysis_config.write_text(
                "project:\n"
                "  trialID: test\n"
                "  compiled_physeq: null\n"
                "ordination: {}\n"
            )
            with self.assertRaisesRegex(SystemExit, "project.compiled_physeq"):
                cli.validate_analysis_config(analysis_config, ["ordination"])


class MetaCommandTests(unittest.TestCase):
    def test_meta_compile_command_construction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            meta_config = Path(tmp) / "meta.yaml"
            write_minimal_meta_config(meta_config)
            args = Args(
                analysis_config=str(meta_config),
                log_root=str(Path(tmp) / "meta_logs"),
                log_dir=None,
                dry_run=True,
                env_mode="direct",
                manager="auto",
            )
            spec = cli.build_meta_run_spec(
                args,
                "compile-phyloseq",
                created_utc=datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc),
            )
            self.assertEqual(spec.expanded_steps, ["compile-phyloseq"])
            self.assertEqual(
                spec.commands[0].command,
                ["Rscript", "scripts/PhyloseqCompiler.R", "--analysis-config", str(meta_config)],
            )
            self.assertEqual(spec.commands[0].command, spec.commands[0].execution_command)
            self.assertEqual(spec.log_dir, Path(tmp) / "meta_logs" / "20260102T030405Z_compile-phyloseq")

    def test_meta_decontaminate_phyloseq_command_construction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            meta_config = Path(tmp) / "meta.yaml"
            write_minimal_meta_config(meta_config)
            args = Args(
                analysis_config=str(meta_config),
                log_root=str(Path(tmp) / "meta_logs"),
                log_dir=None,
                dry_run=True,
                env_mode="direct",
                manager="auto",
            )
            spec = cli.build_meta_run_spec(
                args,
                "decontaminate-phyloseq",
                created_utc=datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc),
            )
            self.assertEqual(spec.expanded_steps, ["decontaminate-phyloseq"])
            self.assertEqual(
                spec.commands[0].command,
                ["Rscript", "scripts/PhyloseqDecontamination.R", "--analysis-config", str(meta_config)],
            )
            self.assertEqual(spec.commands[0].command, spec.commands[0].execution_command)
            self.assertEqual(spec.log_dir, Path(tmp) / "meta_logs" / "20260102T030405Z_decontaminate-phyloseq")

    def test_meta_decontaminate_prefix_env_uses_microclean_prefix(self) -> None:
        args = Args(
            r_env_prefix=None,
            bio_env_prefix=None,
            microclean_env_prefix="/data/taylorng/conda/envs/low-bm-microclean",
        )
        prefixes = cli.resolve_analysis_env_prefixes(
            args,
            ["decontaminate-phyloseq"],
            cli.META_STEP_REGISTRY,
            validate_exists=False,
        )
        self.assertEqual(prefixes[cli.MICROCLEAN_ENV], Path(args.microclean_env_prefix))

    def test_meta_dry_run_writes_provenance_without_execution(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            meta_config = Path(tmp) / "meta.yaml"
            write_minimal_meta_config(meta_config)
            log_dir = Path(tmp) / "meta_run"
            args = Args(
                analysis_config=str(meta_config),
                log_root=str(Path(tmp) / "meta_logs"),
                log_dir=str(log_dir),
                dry_run=True,
                print_command=False,
                env_mode="direct",
                manager="auto",
            )
            self.assertEqual(cli.meta_compile_phyloseq_command(args), 0)
            metadata = json.loads((log_dir / "meta_metadata.json").read_text())
            self.assertEqual(metadata["expanded_steps"], ["compile-phyloseq"])
            self.assertTrue(metadata["dry_run"])
            self.assertIn("# compile-phyloseq", (log_dir / "meta_commands.txt").read_text())
            self.assertFalse((log_dir / "01_compile-phyloseq.log").exists())

    def test_meta_missing_config_has_setup_hint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            missing_config = Path(tmp) / "config" / "local" / "meta.yaml"
            old_default = cli.DEFAULT_META_CONFIG
            try:
                cli.DEFAULT_META_CONFIG = str(missing_config)
                with self.assertRaisesRegex(SystemExit, "low-bm meta init"):
                    cli.resolve_meta_config(missing_config)
            finally:
                cli.DEFAULT_META_CONFIG = old_default


class ProfileConfigurationTests(unittest.TestCase):
    def test_slurm_profile_does_not_array_batch_level_sample_prep(self) -> None:
        profile_text = (REPO_ROOT / "profiles" / "slurm" / "config.yaml").read_text()
        self.assertNotIn("slurm-array-jobs", profile_text)
        self.assertNotIn("slurm-array-limit", profile_text)

    def test_profiles_do_not_configure_bwa_index_rule(self) -> None:
        for profile in ("slurm", "local"):
            profile_text = (REPO_ROOT / "profiles" / profile / "config.yaml").read_text()
            self.assertNotIn("index_ref_bwa", profile_text)

    def test_snakefile_consumes_bwa_indexes_without_building_them(self) -> None:
        text = (REPO_ROOT / "Snakefile").read_text()
        self.assertNotIn("rule index_ref_bwa:", text)
        self.assertIn("def bwa_index_files(reference):", text)
        self.assertIn("reference_indexes = lambda wc: bwa_index_files", text)
        self.assertIn("reference_indexes = bwa_index_files(BACT16S_REF)", text)

    def test_snakefile_uses_manifest_and_batch_level_sample_prep(self) -> None:
        text = (REPO_ROOT / "Snakefile").read_text()
        self.assertIn("rule validate_fastqs:", text)
        self.assertIn("FASTQ_MANIFEST = ", text)
        self.assertIn("NORM_FASTQ_DONE", text)
        self.assertIn("UMI_SELECTION_DONE", text)
        self.assertIn("UMI_DEDUP_DONE", text)
        self.assertIn("count_summary_done = count_summary_done()", text)
        self.assertIn("sample_prep_done = sample_prep_done()", text)


class RunnerSetupTests(unittest.TestCase):
    def test_build_runner_create_command_uses_manager_specific_syntax(self) -> None:
        prefix = Path(".low-bm/runner/env")
        env_file = Path("workflow/envs/runner-env.yaml")
        conda = EnvManagerSpec(name="conda", executable="/opt/conda/bin/conda")
        micro = EnvManagerSpec(name="micromamba", executable="/usr/bin/micromamba")
        self.assertEqual(
            build_runner_create_command(prefix, env_file, conda),
            [
                "/opt/conda/bin/conda",
                "env",
                "create",
                "--yes",
                "--prefix",
                str(prefix),
                "--file",
                str(env_file),
            ],
        )
        self.assertEqual(
            build_runner_create_command(prefix, env_file, micro)[0:2],
            ["/usr/bin/micromamba", "create"],
        )

    def test_doctor_runner_checks_real_conda_base(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            args = Args(
                mode="local",
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
            )

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                self.assertEqual(cli.doctor_runner_command(args), 0)

            self.assertIn("[ok] conda base activation", stdout.getvalue())

    def test_doctor_runner_optional_rule_env_smoke_test(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            runner_prefix, _manager = write_fake_runner(Path(tmp))
            args = Args(
                mode="local",
                runner="host",
                runner_prefix=str(runner_prefix),
                manager="auto",
                rule_env_smoke_test=True,
            )

            stdout = io.StringIO()
            with patch.object(cli, "run_rule_env_smoke_test", return_value=(True, "rule-env/bin/python")):
                with redirect_stdout(stdout):
                    self.assertEqual(cli.doctor_runner_command(args), 0)

            self.assertIn("[ok] rule env python smoke test: rule-env/bin/python", stdout.getvalue())

    def test_setup_runner_records_detected_conda_base(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env_file = root / "runner-env.yaml"
            env_file.write_text("name: runner\n")
            conda_base = root / "conda-base"
            make_fake_conda_base(conda_base)
            manager = root / "bin" / "mamba"
            manager.parent.mkdir(parents=True)
            manager.write_text(
                "#!/usr/bin/env bash\n"
                "if [[ \"$1\" == \"info\" && \"$2\" == \"--base\" ]]; then\n"
                f"  printf '%s\\n' {shlex.quote(str(conda_base))}\n"
                "  exit 0\n"
                "fi\n"
                "if [[ \"$1\" == \"env\" && \"$2\" == \"create\" ]]; then\n"
                "  mkdir -p \"$5/bin\"\n"
                "  printf '#!/usr/bin/env bash\\nexit 0\\n' > \"$5/bin/snakemake\"\n"
                "  printf '#!/usr/bin/env bash\\nexit 0\\n' > \"$5/bin/python\"\n"
                "  printf '#!/usr/bin/env bash\\nexit 0\\n' > \"$5/bin/conda\"\n"
                "  chmod +x \"$5/bin/snakemake\" \"$5/bin/python\" \"$5/bin/conda\"\n"
                "  exit 0\n"
                "fi\n"
                "exit 0\n"
            )
            manager.chmod(0o755)

            old_env_file = cli.DEFAULT_RUNNER_ENV_FILE
            try:
                cli.DEFAULT_RUNNER_ENV_FILE = str(env_file.relative_to(REPO_ROOT)) if env_file.is_relative_to(REPO_ROOT) else str(env_file)
                args = Args(
                    runner_prefix=str(root / "runner" / "env"),
                    manager=str(manager),
                    force=False,
                    conda_base_prefix=None,
                )
                with patch.object(cli, "REPO_ROOT", root):
                    self.assertEqual(cli.setup_runner_command(args), 0)
            finally:
                cli.DEFAULT_RUNNER_ENV_FILE = old_env_file

            metadata = json.loads(cli.runner_metadata_path(root / "runner" / "env").read_text())
            self.assertEqual(metadata["conda_base_prefix"], str(conda_base.resolve()))

    def test_resolve_manager_uses_runner_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            runner_prefix, manager = write_fake_runner(Path(tmp))
            resolved = resolve_env_manager("auto", cli.runner_metadata_path(runner_prefix))
            self.assertEqual(resolved.name, "mamba")
            self.assertEqual(resolved.executable, str(manager))

    def test_managed_microclean_prefix_includes_post_deploy_and_source_recipe(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env_dir = root / "workflow" / "envs"
            env_dir.mkdir(parents=True)
            (env_dir / "micRoclean-env.yaml").write_text("name: low-bm-microclean\n")
            post_deploy = env_dir / "micRoclean-env.post-deploy.sh"
            post_deploy.write_text("#!/usr/bin/env bash\necho v1\n")
            source = env_dir / "micRoclean-source.env"
            source.write_text("MICROCLEAN_GIT_REF=aaaaaaaa\nSCRUB_GIT_REF=bbbbbbbb\n")
            env_root = root / ".low-bm" / "analysis" / "envs"

            with patch.object(cli, "REPO_ROOT", root):
                first_prefix = cli.managed_analysis_env_prefix(cli.MICROCLEAN_ENV, env_root)
                post_deploy.write_text("#!/usr/bin/env bash\necho v2\n")
                post_deploy_prefix = cli.managed_analysis_env_prefix(cli.MICROCLEAN_ENV, env_root)
                source.write_text("MICROCLEAN_GIT_REF=cccccccc\nSCRUB_GIT_REF=bbbbbbbb\n")
                source_prefix = cli.managed_analysis_env_prefix(cli.MICROCLEAN_ENV, env_root)

            self.assertNotEqual(first_prefix, post_deploy_prefix)
            self.assertNotEqual(post_deploy_prefix, source_prefix)

    def test_managed_r_tools_prefix_includes_post_deploy_recipe(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env_dir = root / "workflow" / "envs"
            env_dir.mkdir(parents=True)
            (env_dir / "R-tools-env.yaml").write_text("name: low-bm-r-tools\n")
            post_deploy = env_dir / "R-tools-env.post-deploy.sh"
            post_deploy.write_text("#!/usr/bin/env bash\necho v1\n")
            env_root = root / ".low-bm" / "analysis" / "envs"

            with patch.object(cli, "REPO_ROOT", root):
                first_prefix = cli.managed_analysis_env_prefix(cli.R_TOOLS_ENV, env_root)
                post_deploy.write_text("#!/usr/bin/env bash\necho v2\n")
                post_deploy_prefix = cli.managed_analysis_env_prefix(cli.R_TOOLS_ENV, env_root)

            self.assertNotEqual(first_prefix, post_deploy_prefix)

    def test_r_tools_post_deploy_installs_pinned_coda4microbiome(self) -> None:
        script = (REPO_ROOT / "workflow/envs/R-tools-env.post-deploy.sh").read_text()
        self.assertIn('CODA4MICROBIOME_VERSION="${CODA4MICROBIOME_VERSION:-0.2.4}"', script)
        self.assertIn('remotes::install_version(', script)
        self.assertIn('"coda4microbiome"', script)
        self.assertIn('ANCOMBC_GIT_REF="${ANCOMBC_GIT_REF:-4595750750e354dfa61645f4a3f1f6c53645f683}"', script)
        self.assertIn("Installing upstream ANCOMBC with quadprog trend optimization", script)
        self.assertIn('requireNamespace("microbiome", quietly = TRUE)', script)
        self.assertIn('requireNamespace("quadprog", quietly = TRUE)', script)
        self.assertIn('!"CVXR" %in% names(getNamespaceImports("ANCOMBC"))', script)

    def test_r_tools_env_includes_upstream_ancombc_runtime_dependencies(self) -> None:
        env = (REPO_ROOT / "workflow/envs/R-tools-env.yaml").read_text()
        self.assertIn("  - bioconductor-microbiome\n", env)
        self.assertIn("  - r-quadprog\n", env)

    def test_microclean_post_deploy_patches_tibble_namespace_import(self) -> None:
        script = (REPO_ROOT / "workflow/envs/micRoclean-env.post-deploy.sh").read_text()
        self.assertIn('required_imports <- c("stringr", "tibble")', script)
        self.assertIn('get("add_column", envir = ns, inherits = TRUE)', script)
        self.assertIn("ANCOMBC_GIT_REF", script)
        self.assertIn("Installing upstream ANCOMBC with quadprog trend optimization", script)
        self.assertIn('requireNamespace("microbiome", quietly = TRUE)', script)
        self.assertIn('requireNamespace("quadprog", quietly = TRUE)', script)
        self.assertIn('!"CVXR" %in% names(getNamespaceImports("ANCOMBC"))', script)

    def test_microclean_env_includes_upstream_ancombc_runtime_dependencies(self) -> None:
        env = (REPO_ROOT / "workflow/envs/micRoclean-env.yaml").read_text()
        self.assertIn("  - bioconductor-microbiome\n", env)
        self.assertIn("  - r-quadprog\n", env)

    def test_microclean_source_requires_fixed_shas(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "micRoclean-source.env"
            source.write_text(
                "MICROCLEAN_GIT_URL=https://example.org/micro.git\n"
                "MICROCLEAN_GIT_REF=8d823ed8ad11f2d4da91ebbf6b75710f9af082d7\n"
                "SCRUB_GIT_URL=https://example.org/scrub.git\n"
                "SCRUB_GIT_REF=fcbb8524190f0b27b7ad52cde232c8c4f59810e0\n"
                "ANCOMBC_GIT_URL=https://example.org/ancombc.git\n"
                "ANCOMBC_GIT_REF=4595750750e354dfa61645f4a3f1f6c53645f683\n"
            )
            self.assertEqual(validate_microclean_source(source), [])
            source.write_text(
                "MICROCLEAN_GIT_URL=https://example.org/micro.git\n"
                "MICROCLEAN_GIT_REF=main\n"
                "SCRUB_GIT_URL=https://example.org/scrub.git\n"
                "SCRUB_GIT_REF=fcbb8524190f0b27b7ad52cde232c8c4f59810e0\n"
                "ANCOMBC_GIT_URL=https://example.org/ancombc.git\n"
                "ANCOMBC_GIT_REF=4595750750e354dfa61645f4a3f1f6c53645f683\n"
            )
            self.assertIn("MICROCLEAN_GIT_REF must be a fixed Git SHA", validate_microclean_source(source))


def write_fake_runner(root: Path) -> tuple[Path, Path]:
    runner_prefix = root / "runner" / "env"
    runner_prefix.mkdir(parents=True)
    runner_bin = runner_prefix / "bin"
    runner_bin.mkdir(parents=True)
    for executable in ("snakemake", "python", "conda"):
        path = runner_bin / executable
        path.write_text("#!/usr/bin/env bash\nexit 0\n")
        path.chmod(0o755)
    conda_base = root / "conda-base"
    make_fake_conda_base(conda_base)
    manager = root / "bin" / "mamba"
    manager.parent.mkdir(parents=True)
    manager.write_text(
        "#!/usr/bin/env bash\n"
        "if [[ \"$1\" == \"info\" && \"$2\" == \"--base\" ]]; then\n"
        f"  printf '%s\\n' {shlex.quote(str(conda_base))}\n"
        "  exit 0\n"
        "fi\n"
        "exit 0\n"
    )
    manager.chmod(0o755)
    metadata = {
        "manager": "mamba",
        "manager_executable": str(manager),
        "conda_base_prefix": str(conda_base.resolve()),
    }
    cli.runner_metadata_path(runner_prefix).write_text(json.dumps(metadata) + "\n")
    return runner_prefix, manager


def make_fake_conda_base(prefix: Path) -> None:
    bin_dir = prefix / "bin"
    bin_dir.mkdir(parents=True)
    activate = bin_dir / "activate"
    activate.write_text("#!/usr/bin/env bash\nreturn 0 2>/dev/null || exit 0\n")
    activate.chmod(0o755)
    conda = bin_dir / "conda"
    conda.write_text("#!/usr/bin/env bash\nexit 0\n")
    conda.chmod(0o755)


def write_minimal_processing_config(path: Path) -> None:
    path.write_text(
        "in_root: Exp_Data/\n"
        "ip_root: IP_Data/\n"
        "out_root: Exp_Output/\n"
        "script_dir: scripts\n"
    )


def write_processing_config_with_refs(path: Path, refs: Path) -> None:
    refs.mkdir(parents=True)
    human = refs / "human.fa"
    mouse = refs / "mouse.fa"
    viral = refs / "viral.fa"
    bacteria = refs / "bacteria.fa"
    for fasta in (human, mouse, viral, bacteria):
        fasta.write_text(">ref\nACGT\n")
    path.write_text(
        "in_root: Exp_Data/\n"
        "ip_root: IP_Data/\n"
        "out_root: Exp_Output/\n"
        "script_dir: scripts\n"
        "host: human\n"
        f"human_ref: {human}\n"
        f"mouse_ref: {mouse}\n"
        f"viral_ref: {viral}\n"
        f"bact16s_ref: {bacteria}\n"
    )


def write_bwa_sidecars(reference: Path) -> None:
    for extension in cli.BWA_INDEX_EXTENSIONS:
        Path(f"{reference}{extension}").write_text("index\n")


def write_bwa_manifest(reference: Path, reference_sha256: str | None = None) -> None:
    manifest = {
        "schema_version": cli.BWA_INDEX_MANIFEST_SCHEMA_VERSION,
        "reference_path": str(reference),
        "reference_sha256": reference_sha256 or cli.hash_file(reference),
        "bwa_version": "0.7.17-test",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sidecars": [str(Path(f"{reference}{extension}")) for extension in cli.BWA_INDEX_EXTENSIONS],
    }
    cli.bwa_index_manifest_path(reference).write_text(json.dumps(manifest, indent=2) + "\n")


def write_minimal_analysis_config(path: Path) -> None:
    path.write_text(
        "project:\n"
        "  trialID: test\n"
        "  output_dir: Exp_Output/test\n"
        "  compiled_physeq: Exp_Output/test/test_physeq.RData\n"
        "ordination: {}\n"
        "abundance_barplots:\n"
        "  plots:\n"
        "    - name: SampleType\n"
        "      x: SampleID\n"
        "differential_abundance:\n"
        "  comparisons:\n"
        "    - name: TestComparison\n"
        "lefse_analysis: {}\n"
        "heatmap_violin: {}\n"
        "xgboost_classification:\n"
        "  models:\n"
        "    - name: TestClassifier\n"
        "      target:\n"
        "        column: SampleType\n"
        "        negative: [Control]\n"
        "        positive: [Treatment]\n"
        "survival_analysis:\n"
        "  analyses:\n"
        "    - name: TestSurvival\n"
        "blast_confirmation:\n"
        "  candidate_comparisons: [Test]\n"
    )


def write_minimal_meta_config(path: Path) -> None:
    path.write_text(
        "project:\n"
        "  base_dir: Exp_Output\n"
        "  output_dir: Exp_Output/meta\n"
        "meta_compile:\n"
        "  batch_table: config/local/batch.tsv\n"
        "  output_physeq: CompPhyseq.RData\n"
        "  output_asv_fasta: MergedASV.fasta\n"
        "meta_decontamination:\n"
        "  input_physeq: Exp_Output/test/test_physeq.RData\n"
        "  output_dir: Exp_Output/test_microclean\n"
        "  output_physeq: DecontamPhyseq.RData\n"
        "meta_differential_abundance:\n"
        "  batch1_name: 010126.1_batch1\n"
        "  batch2_name: 010126.2_batch2\n"
        "  DA_method: ANCOMBC\n"
        "  comparison: PatientSample\n"
    )


class Args:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


if __name__ == "__main__":
    unittest.main()
