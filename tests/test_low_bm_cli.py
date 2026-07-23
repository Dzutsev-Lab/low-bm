#!/usr/bin/env python3

from __future__ import annotations

import tempfile
import unittest
import json
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
            self.assertLess(command.index("all"), command.index("--configfile"))
            self.assertEqual(
                [
                    command[i + 1]
                    for i, token in enumerate(command)
                    if token == "--configfile"
                ],
                [str(base_config), str(override_config), str(run_config)],
            )
            self.assertEqual(command[-1], "--quiet")

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
            self.assertEqual(command[0:4], [str(manager), "run", "--prefix", str(runner_prefix)])
            self.assertIn("snakemake", command)
            self.assertLess(command.index("all"), command.index("--configfile"))

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


class Args:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


if __name__ == "__main__":
    unittest.main()
