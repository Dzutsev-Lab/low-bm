#!/usr/bin/env python3
"""Batch-level FASTQ validation and lightweight sample-prep helpers."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


RAW_FASTQ_RE = re.compile(r"^(?P<prefix>.+)_R(?P<read>[12])_001\.fastq(?:\.gz)?$")
MANIFEST_FIELDS = [
    "SampleID",
    "R1Path",
    "R2Path",
    "R1Compressed",
    "R2Compressed",
    "R1SourceName",
    "R2SourceName",
]


@dataclass(frozen=True)
class FastqManifestResult:
    rows: list[dict[str, str]]
    errors: list[str]
    warnings: list[str]
    ignored: list[str]
    malformed: list[str]


@dataclass(frozen=True)
class StageResult:
    sample_id: str
    ok: bool
    detail: str


def parse_raw_fastq_name(path: str | Path) -> tuple[str, int, str] | None:
    name = Path(path).name
    match = RAW_FASTQ_RE.match(name)
    if not match:
        return None

    source_prefix = match.group("prefix")
    sample = source_prefix
    while True:
        normalized = re.sub(r"_(?:S\d+|L\d{3})$", "", sample)
        if normalized == sample:
            break
        sample = normalized
    return sample, int(match.group("read")), source_prefix


def iter_fastq_files(raw_dir: str | Path) -> list[Path]:
    raw_path = Path(raw_dir)
    candidates = set(raw_path.glob("*.fastq"))
    candidates.update(raw_path.glob("*.fastq.gz"))
    return sorted(candidates)


def build_fastq_manifest(raw_dir: str | Path) -> FastqManifestResult:
    grouped: dict[str, dict[int, list[tuple[Path, str]]]] = {}
    ignored: list[str] = []
    malformed: list[str] = []

    for fastq in iter_fastq_files(raw_dir):
        parsed = parse_raw_fastq_name(fastq)
        if parsed is None:
            malformed.append(str(fastq))
            continue

        sample_id, read, source_prefix = parsed
        if "Undetermined" in sample_id or "Undetermined" in fastq.name:
            ignored.append(str(fastq))
            continue

        grouped.setdefault(sample_id, {1: [], 2: []})[read].append((fastq, source_prefix))

    rows: list[dict[str, str]] = []
    errors: list[str] = []

    for sample_id in sorted(grouped):
        reads = grouped[sample_id]
        r1s = sorted(reads[1], key=lambda item: str(item[0]))
        r2s = sorted(reads[2], key=lambda item: str(item[0]))

        if len(r1s) != 1 or len(r2s) != 1:
            if not r1s:
                errors.append(f"{sample_id}: missing R1 FASTQ.")
            if not r2s:
                errors.append(f"{sample_id}: missing R2 FASTQ.")
            if len(r1s) > 1:
                errors.append(
                    f"{sample_id}: multiple R1 FASTQs after sample-name normalization: "
                    + ", ".join(str(path) for path, _prefix in r1s)
                )
            if len(r2s) > 1:
                errors.append(
                    f"{sample_id}: multiple R2 FASTQs after sample-name normalization: "
                    + ", ".join(str(path) for path, _prefix in r2s)
                )
            continue

        r1_path, r1_prefix = r1s[0]
        r2_path, r2_prefix = r2s[0]
        if r1_prefix != r2_prefix:
            errors.append(
                f"{sample_id}: R1/R2 raw prefixes differ after sample-name normalization: "
                f"{r1_path.name} vs {r2_path.name}"
            )
            continue

        rows.append(
            {
                "SampleID": sample_id,
                "R1Path": str(r1_path.resolve()),
                "R2Path": str(r2_path.resolve()),
                "R1Compressed": str(str(r1_path).endswith(".gz")).lower(),
                "R2Compressed": str(str(r2_path).endswith(".gz")).lower(),
                "R1SourceName": r1_path.name,
                "R2SourceName": r2_path.name,
            }
        )

    if malformed:
        errors.append(
            "Malformed FASTQ filename(s): "
            + ", ".join(malformed)
            + ". Expected *_R1_001.fastq(.gz) and *_R2_001.fastq(.gz)."
        )
    if not rows:
        errors.append(f"No valid paired samples found in {Path(raw_dir)}.")

    warnings = [
        "Ignored Undetermined FASTQ(s): " + ", ".join(ignored),
    ] if ignored else []
    return FastqManifestResult(rows=rows, errors=errors, warnings=warnings, ignored=ignored, malformed=malformed)


def write_manifest(rows: list[dict[str, str]], path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_manifest(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    missing = [field for field in MANIFEST_FIELDS if field not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"Manifest is missing required column(s): {', '.join(missing)}")
    return rows


def write_sample_names(rows: list[dict[str, str]], path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(row["SampleID"] for row in rows) + ("\n" if rows else ""))


def write_validation_report(result: FastqManifestResult, raw_dir: str | Path, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "FASTQ validation report",
        f"Raw FASTQ directory: {Path(raw_dir)}",
        f"Valid paired samples: {len(result.rows)}",
        "",
    ]
    if result.rows:
        lines.append("Samples:")
        lines.extend(f"  - {row['SampleID']}" for row in result.rows)
        lines.append("")
    if result.warnings:
        lines.append("Warnings:")
        lines.extend(f"  - {warning}" for warning in result.warnings)
        lines.append("")
    if result.errors:
        lines.append("Errors:")
        lines.extend(f"  - {error}" for error in result.errors)
    else:
        lines.append("No validation errors detected.")
    destination.write_text("\n".join(lines) + "\n")


def cmd_validate(args: argparse.Namespace) -> int:
    result = build_fastq_manifest(args.raw_dir)
    write_manifest(result.rows, args.manifest)
    write_sample_names(result.rows, args.sample_names)
    write_validation_report(result, args.raw_dir, args.report)

    ok_path = Path(args.ok)
    if result.errors:
        if ok_path.exists():
            ok_path.unlink()
        print(f"FASTQ validation failed. See {args.report}.", file=sys.stderr)
        for error in result.errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    ok_path.parent.mkdir(parents=True, exist_ok=True)
    ok_path.write_text("ok\n")
    print(f"Validated {len(result.rows)} paired FASTQ sample(s).")
    return 0


def is_gzip_path(path: str | Path) -> bool:
    return str(path).endswith(".gz")


def replace_symlink(source: Path, destination: Path) -> None:
    if destination.exists() or destination.is_symlink():
        if destination.is_dir() and not destination.is_symlink():
            raise IsADirectoryError(destination)
        destination.unlink()
    destination.symlink_to(source.resolve())


def normalize_one_fastq(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if is_gzip_path(source):
        tmp = destination.with_suffix(destination.suffix + ".tmp")
        with gzip.open(source, "rb") as src, tmp.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        tmp.replace(destination)
    else:
        replace_symlink(source, destination)


def normalized_fastq_path(norm_dir: str | Path, sample_id: str, read: int) -> Path:
    return Path(norm_dir) / f"{sample_id}_R{read}_001.fastq"


def selected_umi_path(selected_dir: str | Path, sample_id: str) -> Path:
    return Path(selected_dir) / f"Selected.{sample_id}.UMI_R1.fastq"


def count_summary_path(selected_dir: str | Path, sample_id: str) -> Path:
    return Path(selected_dir) / f"CountSummary.{sample_id}.tsv"


def deduped_fastq_path(dedup_dir: str | Path, sample_id: str) -> Path:
    return Path(dedup_dir) / f"Deduped.{sample_id}.fastq"


def run_stage(
    rows: Iterable[dict[str, str]],
    threads: int,
    worker: Callable[[dict[str, str]], StageResult],
    done_path: str | Path,
    stage_name: str,
) -> int:
    row_list = list(rows)
    workers = max(1, int(threads or 1))
    failures: list[StageResult] = []

    def guarded(row: dict[str, str]) -> StageResult:
        try:
            return worker(row)
        except Exception as exc:  # noqa: BLE001 - summarize any sample-level failure.
            return StageResult(row.get("SampleID", "<unknown>"), False, str(exc))

    if workers == 1:
        results = [guarded(row) for row in row_list]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(guarded, row) for row in row_list]
            results = [future.result() for future in as_completed(futures)]

    for result in sorted(results, key=lambda item: item.sample_id):
        marker = "ok" if result.ok else "fail"
        print(f"[{marker}] {stage_name}: {result.sample_id}: {result.detail}")
        if not result.ok:
            failures.append(result)

    if failures:
        print(f"{stage_name} failed for {len(failures)} sample(s):", file=sys.stderr)
        for failure in sorted(failures, key=lambda item: item.sample_id):
            print(f"- {failure.sample_id}: {failure.detail}", file=sys.stderr)
        return 1

    done = Path(done_path)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("ok\n")
    return 0


def cmd_norm_fastq(args: argparse.Namespace) -> int:
    rows = read_manifest(args.manifest)

    def worker(row: dict[str, str]) -> StageResult:
        sample_id = row["SampleID"]
        normalize_one_fastq(Path(row["R1Path"]), normalized_fastq_path(args.out_dir, sample_id, 1))
        normalize_one_fastq(Path(row["R2Path"]), normalized_fastq_path(args.out_dir, sample_id, 2))
        return StageResult(sample_id, True, "normalized R1/R2")

    return run_stage(rows, args.threads, worker, args.done, "norm_fastq")


def cmd_umi_selection(args: argparse.Namespace) -> int:
    rows = read_manifest(args.manifest)
    log_dir = Path(args.sample_log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    def worker(row: dict[str, str]) -> StageResult:
        sample_id = row["SampleID"]
        out_umi = selected_umi_path(args.out_dir, sample_id)
        out_summary = count_summary_path(args.out_dir, sample_id)
        out_umi.parent.mkdir(parents=True, exist_ok=True)
        out_summary.parent.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / f"01_umi_select.{sample_id}.log"
        command = [
            sys.executable,
            args.umi_selection_script,
            "--sample-name",
            sample_id,
            "--r1",
            str(normalized_fastq_path(args.norm_dir, sample_id, 1)),
            "--r2",
            str(normalized_fastq_path(args.norm_dir, sample_id, 2)),
            "--r2-primer-motif",
            args.r2_primer_motif,
            "--poly-G-threshold",
            str(args.poly_g_threshold),
            "--umi-len",
            str(args.umi_len),
            "--max-offset",
            str(args.max_offset),
            "--out-count-summary",
            str(out_summary),
            "--out-umi-r1",
            str(out_umi),
        ]
        if args.r2_primer_skip:
            command.extend(["--r2-primer-skip"])
        with log_path.open("w") as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        if completed.returncode != 0:
            return StageResult(sample_id, False, f"see {log_path}")
        return StageResult(sample_id, True, f"see {log_path}")

    return run_stage(rows, args.threads, worker, args.done, "umi_selection")


def count_lines(path: str | Path) -> int:
    with Path(path).open("rb") as handle:
        return sum(1 for _ in handle)


def cmd_umi_dedup(args: argparse.Namespace) -> int:
    rows = read_manifest(args.manifest)
    log_dir = Path(args.sample_log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    def worker(row: dict[str, str]) -> StageResult:
        sample_id = row["SampleID"]
        selected = selected_umi_path(args.selected_dir, sample_id)
        output = deduped_fastq_path(args.out_dir, sample_id)
        output.parent.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / f"02_umi_dedup.{sample_id}.log"
        n_lines = count_lines(selected)
        with log_path.open("w") as log:
            if n_lines == 0:
                print(f"Input FASTQ ({selected}) is empty - skipping AmpUMI", file=log)
                output.write_text("")
                return StageResult(sample_id, True, f"empty input; see {log_path}")

            print(f"Running AmpUMI on {n_lines} lines from {selected}", file=log)
            completed = subprocess.run(
                [
                    "AmpUMI",
                    "Process",
                    "--fastq",
                    str(selected),
                    "--fastq_out",
                    str(output),
                    "--umi_regex",
                    args.umi_regex,
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode != 0:
            return StageResult(sample_id, False, f"see {log_path}")
        return StageResult(sample_id, True, f"see {log_path}")

    return run_stage(rows, args.threads, worker, args.done, "umi_dedup")


def cmd_no_umi_count_summary(args: argparse.Namespace) -> int:
    rows = read_manifest(args.manifest)

    def worker(row: dict[str, str]) -> StageResult:
        sample_id = row["SampleID"]
        r1 = normalized_fastq_path(args.norm_dir, sample_id, 1)
        n_lines = count_lines(r1)
        if n_lines % 4 != 0:
            return StageResult(sample_id, False, f"{r1} has {n_lines} lines, not a multiple of 4")
        out_summary = count_summary_path(args.out_dir, sample_id)
        out_summary.parent.mkdir(parents=True, exist_ok=True)
        reads = n_lines // 4
        out_summary.write_text(f"SampleID\tRaw_reads\tSelected_reads\n{sample_id}\t{reads}\t{reads}\n")
        return StageResult(sample_id, True, f"{reads} reads")

    return run_stage(rows, args.threads, worker, args.done, "no_umi_count_summary")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate")
    validate.add_argument("--raw-dir", required=True)
    validate.add_argument("--manifest", required=True)
    validate.add_argument("--sample-names", required=True)
    validate.add_argument("--report", required=True)
    validate.add_argument("--ok", required=True)
    validate.set_defaults(func=cmd_validate)

    norm = subparsers.add_parser("norm-fastq")
    norm.add_argument("--manifest", required=True)
    norm.add_argument("--out-dir", required=True)
    norm.add_argument("--threads", type=int, default=1)
    norm.add_argument("--done", required=True)
    norm.set_defaults(func=cmd_norm_fastq)

    umi_select = subparsers.add_parser("umi-selection")
    umi_select.add_argument("--manifest", required=True)
    umi_select.add_argument("--norm-dir", required=True)
    umi_select.add_argument("--out-dir", required=True)
    umi_select.add_argument("--sample-log-dir", required=True)
    umi_select.add_argument("--umi-selection-script", required=True)
    umi_select.add_argument("--r2-primer-motif", required=True)
    umi_select.add_argument("--r2-primer-skip", action="store_true")
    umi_select.add_argument("--poly-G-threshold", dest="poly_g_threshold", type=float, required=True)
    umi_select.add_argument("--umi-len", type=int, required=True)
    umi_select.add_argument("--max-offset", type=int, required=True)
    umi_select.add_argument("--threads", type=int, default=1)
    umi_select.add_argument("--done", required=True)
    umi_select.set_defaults(func=cmd_umi_selection)

    umi_dedup = subparsers.add_parser("umi-dedup")
    umi_dedup.add_argument("--manifest", required=True)
    umi_dedup.add_argument("--selected-dir", required=True)
    umi_dedup.add_argument("--out-dir", required=True)
    umi_dedup.add_argument("--sample-log-dir", required=True)
    umi_dedup.add_argument("--umi-regex", required=True)
    umi_dedup.add_argument("--threads", type=int, default=1)
    umi_dedup.add_argument("--done", required=True)
    umi_dedup.set_defaults(func=cmd_umi_dedup)

    no_umi = subparsers.add_parser("no-umi-count-summary")
    no_umi.add_argument("--manifest", required=True)
    no_umi.add_argument("--norm-dir", required=True)
    no_umi.add_argument("--out-dir", required=True)
    no_umi.add_argument("--threads", type=int, default=1)
    no_umi.add_argument("--done", required=True)
    no_umi.set_defaults(func=cmd_no_umi_count_summary)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
