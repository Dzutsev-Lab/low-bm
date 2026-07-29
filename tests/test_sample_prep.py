from __future__ import annotations

import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import SamplePrep  # noqa: E402


class SamplePrepManifestTests(unittest.TestCase):
    def test_umi_selection_parser_normalizes_poly_g_threshold(self) -> None:
        parser = SamplePrep.build_parser()
        args = parser.parse_args(
            [
                "umi-selection",
                "--manifest",
                "manifest.tsv",
                "--norm-dir",
                "norm",
                "--out-dir",
                "umi",
                "--sample-log-dir",
                "logs",
                "--umi-selection-script",
                "UMISelection.py",
                "--r2-primer-motif",
                "ACGT",
                "--poly-G-threshold",
                "0.75",
                "--umi-len",
                "12",
                "--max-offset",
                "2",
                "--done",
                "umi.done",
            ]
        )

        self.assertEqual(args.poly_g_threshold, 0.75)
        self.assertNotIn("poly_G_threshold", vars(args))

    def test_validate_writes_manifest_and_sample_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw"
            raw.mkdir()
            write_fastq(raw / "SampleA_S1_R1_001.fastq")
            write_fastq(raw / "SampleA_S1_R2_001.fastq")
            write_gzip_fastq(raw / "SampleB_S2_L001_R1_001.fastq.gz")
            write_gzip_fastq(raw / "SampleB_S2_L001_R2_001.fastq.gz")
            write_fastq(raw / "Undetermined_S0_R1_001.fastq")
            write_fastq(raw / "Undetermined_S0_R2_001.fastq")

            manifest = root / "fastq_manifest.tsv"
            sample_names = root / "sample.names"
            report = root / "fastq_validation_report.txt"
            ok = root / "fastq_validation.ok"
            code = SamplePrep.main(
                [
                    "validate",
                    "--raw-dir",
                    str(raw),
                    "--manifest",
                    str(manifest),
                    "--sample-names",
                    str(sample_names),
                    "--report",
                    str(report),
                    "--ok",
                    str(ok),
                ]
            )

            self.assertEqual(code, 0)
            self.assertEqual(sample_names.read_text().splitlines(), ["SampleA", "SampleB"])
            rows = read_tsv(manifest)
            self.assertEqual([row["SampleID"] for row in rows], ["SampleA", "SampleB"])
            self.assertEqual(rows[0]["R1Compressed"], "false")
            self.assertEqual(rows[1]["R1Compressed"], "true")
            self.assertIn("Ignored Undetermined FASTQ", report.read_text())
            self.assertTrue(ok.exists())

    def test_validate_reports_missing_duplicate_and_malformed_fastqs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw"
            raw.mkdir()
            write_fastq(raw / "MissingR2_S1_R1_001.fastq")
            write_fastq(raw / "DupSample_S1_L001_R1_001.fastq")
            write_fastq(raw / "DupSample_S1_L001_R2_001.fastq")
            write_fastq(raw / "DupSample_S1_L002_R1_001.fastq")
            write_fastq(raw / "DupSample_S1_L002_R2_001.fastq")
            write_fastq(raw / "not_a_low_bm_fastq_name.fastq")

            report = root / "report.txt"
            code = SamplePrep.main(
                [
                    "validate",
                    "--raw-dir",
                    str(raw),
                    "--manifest",
                    str(root / "manifest.tsv"),
                    "--sample-names",
                    str(root / "sample.names"),
                    "--report",
                    str(report),
                    "--ok",
                    str(root / "ok"),
                ]
            )

            text = report.read_text()
            self.assertEqual(code, 1)
            self.assertIn("MissingR2: missing R2 FASTQ", text)
            self.assertIn("DupSample: multiple R1 FASTQs", text)
            self.assertIn("DupSample: multiple R2 FASTQs", text)
            self.assertIn("Malformed FASTQ filename", text)

    def test_validate_rejects_mismatched_r1_r2_source_prefixes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw"
            raw.mkdir()
            write_fastq(raw / "SampleA_S1_R1_001.fastq")
            write_fastq(raw / "SampleA_S2_R2_001.fastq")

            result = SamplePrep.build_fastq_manifest(raw)

            self.assertEqual(result.rows, [])
            self.assertTrue(
                any("R1/R2 raw prefixes differ" in error for error in result.errors),
                result.errors,
            )

    def test_norm_fastq_links_plain_fastq_and_decompresses_gzip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw"
            raw.mkdir()
            plain = raw / "Plain_R1_001.fastq"
            zipped = raw / "Plain_R2_001.fastq.gz"
            write_fastq(plain)
            write_gzip_fastq(zipped)
            manifest = root / "manifest.tsv"
            SamplePrep.write_manifest(
                [
                    {
                        "SampleID": "Plain",
                        "R1Path": str(plain),
                        "R2Path": str(zipped),
                        "R1Compressed": "false",
                        "R2Compressed": "true",
                        "R1SourceName": plain.name,
                        "R2SourceName": zipped.name,
                    }
                ],
                manifest,
            )

            out_dir = root / "norm"
            code = SamplePrep.main(
                [
                    "norm-fastq",
                    "--manifest",
                    str(manifest),
                    "--out-dir",
                    str(out_dir),
                    "--threads",
                    "1",
                    "--done",
                    str(out_dir / ".norm_fastq.done"),
                ]
            )

            self.assertEqual(code, 0)
            self.assertTrue((out_dir / "Plain_R1_001.fastq").is_symlink())
            self.assertEqual((out_dir / "Plain_R2_001.fastq").read_text(), FASTQ_TEXT)
            self.assertTrue((out_dir / ".norm_fastq.done").exists())

    def test_umi_dedup_invokes_ampumi_through_python_module(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.tsv"
            selected_dir = root / "selected"
            selected_dir.mkdir()
            out_dir = root / "deduped"
            log_dir = root / "logs"
            SamplePrep.write_manifest(
                [
                    {
                        "SampleID": "SampleA",
                        "R1Path": str(root / "raw_R1.fastq"),
                        "R2Path": str(root / "raw_R2.fastq"),
                        "R1Compressed": "false",
                        "R2Compressed": "false",
                        "R1SourceName": "raw_R1.fastq",
                        "R2SourceName": "raw_R2.fastq",
                    }
                ],
                manifest,
            )
            write_fastq(selected_dir / "Selected.SampleA.UMI_R1.fastq")
            commands: list[list[str]] = []

            def fake_run(command: list[str], **_kwargs: object) -> SimpleNamespace:
                commands.append(command)
                return SimpleNamespace(returncode=0)

            with (
                patch.object(SamplePrep, "ampumi_module_available", return_value=True),
                patch.object(SamplePrep.subprocess, "run", side_effect=fake_run),
            ):
                code = SamplePrep.main(
                    [
                        "umi-dedup",
                        "--manifest",
                        str(manifest),
                        "--selected-dir",
                        str(selected_dir),
                        "--out-dir",
                        str(out_dir),
                        "--sample-log-dir",
                        str(log_dir),
                        "--umi-regex",
                        "^IIII",
                        "--threads",
                        "1",
                        "--done",
                        str(out_dir / ".umi_dedup.done"),
                    ]
                )

            self.assertEqual(code, 0)
            self.assertTrue((out_dir / ".umi_dedup.done").exists())
            self.assertEqual(commands[0][:4], [SamplePrep.sys.executable, "-m", "AmpUMI.AmpUMI", "Process"])
            self.assertNotEqual(commands[0][0], "AmpUMI")

    def test_umi_dedup_reports_missing_ampumi_environment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.tsv"
            selected_dir = root / "selected"
            selected_dir.mkdir()
            out_dir = root / "deduped"
            log_dir = root / "logs"
            SamplePrep.write_manifest(
                [
                    {
                        "SampleID": "SampleA",
                        "R1Path": str(root / "raw_R1.fastq"),
                        "R2Path": str(root / "raw_R2.fastq"),
                        "R1Compressed": "false",
                        "R2Compressed": "false",
                        "R1SourceName": "raw_R1.fastq",
                        "R2SourceName": "raw_R2.fastq",
                    }
                ],
                manifest,
            )
            write_fastq(selected_dir / "Selected.SampleA.UMI_R1.fastq")

            with patch.object(SamplePrep, "ampumi_module_available", return_value=False):
                code = SamplePrep.main(
                    [
                        "umi-dedup",
                        "--manifest",
                        str(manifest),
                        "--selected-dir",
                        str(selected_dir),
                        "--out-dir",
                        str(out_dir),
                        "--sample-log-dir",
                        str(log_dir),
                        "--umi-regex",
                        "^IIII",
                        "--threads",
                        "1",
                        "--done",
                        str(out_dir / ".umi_dedup.done"),
                    ]
                )

            self.assertEqual(code, 1)
            text = (log_dir / "02_umi_dedup.SampleA.log").read_text()
            self.assertIn(f"Python executable: {SamplePrep.sys.executable}", text)
            self.assertIn("AmpUMI is not importable from this Python environment", text)


def write_fastq(path: Path) -> None:
    path.write_text(FASTQ_TEXT)


def write_gzip_fastq(path: Path) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(FASTQ_TEXT)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


FASTQ_TEXT = "@read1\nACGT\n+\n!!!!\n"


if __name__ == "__main__":
    unittest.main()
