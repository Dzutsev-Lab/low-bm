#!/usr/bin/env python3

from __future__ import annotations

import io
import tempfile
import unittest
import json
from contextlib import redirect_stdout
from datetime import datetime, timezone
from pathlib import Path

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
                        "",
                    ]
                )
                + "\n"
            )

            rows = read_batch_table(table)
            self.assertEqual([row.trialID for row in rows], ["010126.1", "010126.2"])
            self.assertEqual(rows[0].trial_name, "010126.1_batch1")
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
                "trialID\ttrial_descript\texp_dir\tmetadata\tprocess_umis\n"
                "010126.1\tbatch with spaces\tBatch1\tmetadata/batch 1.xlsx\tfalse\n"
            )
            row = read_batch_table(table)[0]
            run_config = write_run_config(row, Path(tmp) / "configs")
            text = run_config.read_text()
            self.assertIn('trialID: "010126.1"', text)
            self.assertIn('trial_descript: "batch with spaces"', text)
            self.assertIn('process_umis: "false"', text)
            self.assertNotIn("batch_label", text)
            self.assertNotIn("include_processing", text)
            self.assertNotIn("include_analysis", text)


class CommandConstructionTests(unittest.TestCase):
    def test_snakemake_command_preserves_config_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\n')
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
            run_config.write_text('trialID: "010126.1"\n')
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

    def test_batch_submit_dry_run_reuses_run_builder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\n"
            )
            base_config = Path(tmp) / "processing.yaml"
            write_minimal_processing_config(base_config)
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
                "trialID\ttrial_descript\texp_dir\tmetadata\n"
                "010126.1\tbatch one\tBatch1\tmetadata/batch1.xlsx\n"
                "010126.2\tbatch two\tBatch2\tmetadata/batch2.xlsx\n"
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
            run_config.write_text('trialID: "010126.1"\ntrial_descript: "batch with spaces"\n')
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
            run_config.write_text('trialID: "010126.1"\ntrial_descript: "batch1"\n')
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
                "trialID\ttrial_descript\texp_dir\tmetadata\n"
                "010126.1\tbatch one\tBatch1\tmetadata/batch1.xlsx\n"
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
                "trialID\ttrial_descript\texp_dir\tmetadata\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\n"
                "010126.1\tbatch1\tBatch2\tmetadata/batch2.xlsx\n"
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

    def test_execution_command_wraps_snakemake_in_runner_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\n')
            base_config = Path(tmp) / "processing.yaml"
            write_minimal_processing_config(base_config)
            runner_prefix, manager = write_fake_runner(Path(tmp))
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
            self.assertEqual(command[0:4], [str(manager), "run", "--prefix", str(runner_prefix.resolve())])
            self.assertIn("snakemake", command)
            self.assertEqual(command.count("--configfile"), 1)
            self.assertEqual(command[-2:], ["--", "all"])

    def test_missing_runner_prefix_has_setup_hint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\n')
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
            command = [str(manager), "run", "--prefix", str(runner_prefix), "snakemake", "all"]
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
            self.assertIn(str(manager), text)
            self.assertIn(str(runner_prefix), text)
            self.assertIn("Missing runner env", text)
            self.assertNotIn("activate low-bm-runner", text)

    def test_missing_default_processing_config_has_setup_hint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_config = Path(tmp) / "010126.1_runconfig.yaml"
            run_config.write_text('trialID: "010126.1"\n')
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
    def test_slurm_profile_caps_sample_prep_arrays(self) -> None:
        profile_text = (REPO_ROOT / "profiles" / "slurm" / "config.yaml").read_text()
        self.assertIn(
            "slurm-array-jobs: norm_fastq,umi_selection,umi_dedup,no_umi_count_summary",
            profile_text,
        )
        self.assertRegex(profile_text, r"(?m)^slurm-array-limit:\s*100$")


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

    def test_resolve_manager_uses_runner_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            runner_prefix, manager = write_fake_runner(Path(tmp))
            resolved = resolve_env_manager("auto", cli.runner_metadata_path(runner_prefix))
            self.assertEqual(resolved.name, "mamba")
            self.assertEqual(resolved.executable, str(manager))

    def test_microclean_source_requires_fixed_shas(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "micRoclean-source.env"
            source.write_text(
                "MICROCLEAN_GIT_URL=https://example.org/micro.git\n"
                "MICROCLEAN_GIT_REF=8d823ed8ad11f2d4da91ebbf6b75710f9af082d7\n"
                "SCRUB_GIT_URL=https://example.org/scrub.git\n"
                "SCRUB_GIT_REF=fcbb8524190f0b27b7ad52cde232c8c4f59810e0\n"
            )
            self.assertEqual(validate_microclean_source(source), [])
            source.write_text(
                "MICROCLEAN_GIT_URL=https://example.org/micro.git\n"
                "MICROCLEAN_GIT_REF=main\n"
                "SCRUB_GIT_URL=https://example.org/scrub.git\n"
                "SCRUB_GIT_REF=fcbb8524190f0b27b7ad52cde232c8c4f59810e0\n"
            )
            self.assertIn("MICROCLEAN_GIT_REF must be a fixed Git SHA", validate_microclean_source(source))


def write_fake_runner(root: Path) -> tuple[Path, Path]:
    runner_prefix = root / "runner" / "env"
    runner_prefix.mkdir(parents=True)
    manager = root / "bin" / "mamba"
    manager.parent.mkdir(parents=True)
    manager.write_text("#!/usr/bin/env bash\nexit 0\n")
    manager.chmod(0o755)
    metadata = {
        "manager": "mamba",
        "manager_executable": str(manager),
    }
    cli.runner_metadata_path(runner_prefix).write_text(json.dumps(metadata) + "\n")
    return runner_prefix, manager


def write_minimal_processing_config(path: Path) -> None:
    path.write_text(
        "in_root: Exp_Data/\n"
        "ip_root: IP_Data/\n"
        "out_root: Exp_Output/\n"
        "script_dir: scripts\n"
    )


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
        "xgboost_classification: {}\n"
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
