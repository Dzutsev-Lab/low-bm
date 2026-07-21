#!/usr/bin/env python3

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from low_bm.batch import read_batch_table, write_run_config
from low_bm.cli import batch_submit_command, build_run_spec, build_snakemake_command


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
            args = Args(
                configfile=["config.yaml"],
                extra_configfile=["override.yaml"],
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
            self.assertEqual(
                [
                    command[i + 1]
                    for i, token in enumerate(command)
                    if token == "--configfile"
                ],
                ["config.yaml", "override.yaml", str(run_config)],
            )
            self.assertEqual(command[-2:], ["--quiet", "all"])

    def test_batch_submit_dry_run_reuses_run_builder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "batch.tsv"
            table.write_text(
                "trialID\ttrial_descript\texp_dir\tmetadata\n"
                "010126.1\tbatch1\tBatch1\tmetadata/batch1.xlsx\n"
            )
            args = Args(
                batch_table=str(table),
                run_config_dir=str(Path(tmp) / "configs"),
                configfile=["config.yaml"],
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
            self.assertEqual(batch_submit_command(args), 0)


class ProfileConfigurationTests(unittest.TestCase):
    def test_slurm_profile_caps_sample_prep_arrays(self) -> None:
        profile_text = (REPO_ROOT / "profiles" / "slurm" / "config.yaml").read_text()
        self.assertIn(
            "slurm-array-jobs: norm_fastq,umi_selection,umi_dedup,no_umi_count_summary",
            profile_text,
        )
        self.assertRegex(profile_text, r"(?m)^slurm-array-limit:\s*100$")


class Args:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


if __name__ == "__main__":
    unittest.main()
