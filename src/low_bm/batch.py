"""Batch-table parsing and per-row config generation for low-bm runs."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


CANONICAL_COLUMNS = [
    "trialID",
    "trial_descript",
    "exp_dir",
    "metadata",
    "host",
]

OPTIONAL_CONFIG_COLUMNS = [
    "process_umis",
]

VALID_HOSTS = ("human", "mouse")


@dataclass(frozen=True)
class BatchRow:
    """One normalized row from the experiment batch table."""

    row_number: int
    trialID: str
    trial_descript: str
    exp_dir: str
    metadata: str
    host: str
    process_umis: str | None = None

    @property
    def trial_name(self) -> str:
        return f"{self.trialID}_{self.trial_descript}"

    @property
    def config_stem(self) -> str:
        return f"{self.trialID}_runconfig"

    def config_items(self) -> list[tuple[str, str]]:
        """Return the fields that should be written to a per-row run config."""
        items = [
            ("trialID", self.trialID),
            ("trial_descript", self.trial_descript),
            ("exp_dir", self.exp_dir),
            ("metadata", self.metadata),
            ("host", self.host),
        ]
        if self.process_umis not in (None, ""):
            items.append(("process_umis", self.process_umis or ""))
        return items


def read_batch_table(path: str | Path) -> list[BatchRow]:
    """Read the canonical headered batch TSV."""
    table_path = Path(path)
    with table_path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))

    if not rows:
        raise ValueError(f"Batch table is empty: {table_path}")

    return _read_headered_table(table_path, rows)


def write_run_config(row: BatchRow, output_dir: str | Path) -> Path:
    """Write a Snakemake config override for one batch-table row."""
    config_dir = Path(output_dir)
    config_dir.mkdir(parents=True, exist_ok=True)
    output = config_dir / f"{row.config_stem}.yaml"
    with output.open("w", newline="\n") as handle:
        for key, value in row.config_items():
            handle.write(f"{key}: {_yaml_quote(value)}\n")
    return output


def load_simple_run_config(path: str | Path) -> dict[str, str]:
    """Read the simple scalar YAML emitted by write_run_config."""
    values: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        values[key.strip()] = _yaml_unquote(value.strip())
    return values


def _read_headered_table(table_path: Path, raw_rows: list[list[str]]) -> list[BatchRow]:
    with table_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        records = list(reader)

    fieldnames = reader.fieldnames or []
    missing_cols = [key for key in CANONICAL_COLUMNS if key not in fieldnames]
    if missing_cols:
        suspicious = [
            field
            for field in fieldnames
            if any(key in field for key in CANONICAL_COLUMNS + OPTIONAL_CONFIG_COLUMNS)
            and field not in CANONICAL_COLUMNS
            and field not in OPTIONAL_CONFIG_COLUMNS
        ]
        hint = ""
        if suspicious:
            hint = (
                " Suspicious header value(s): "
                + ", ".join(repr(field) for field in suspicious)
                + ". Check that columns are separated by tabs, not spaces."
            )
        raise ValueError(f"Missing required batch table column(s): {', '.join(missing_cols)}.{hint}")

    parsed: list[BatchRow] = []
    for index, record in enumerate(records, start=1):
        raw_row = raw_rows[index]
        if len(raw_row) != len(fieldnames):
            raise ValueError(
                f"Batch row {index} has {len(raw_row)} tab-separated field(s), "
                f"but the header has {len(fieldnames)} column(s)."
            )
        row = {key: (record.get(key, "") or "").strip() for key in CANONICAL_COLUMNS}
        for required in CANONICAL_COLUMNS:
            if not row.get(required):
                raise ValueError(f"Missing required batch table value: {required} in row {index}")
        row["host"] = normalize_host(row["host"], context=f"batch row {index}")
        optional = {
            key: (record.get(key, "") or "").strip()
            for key in OPTIONAL_CONFIG_COLUMNS
            if key in fieldnames
        }
        parsed.append(BatchRow(row_number=index, **row, **optional))
    return parsed


def normalize_host(value: str, context: str = "host") -> str:
    """Return a normalized host value or raise a clear batch-config error."""
    host = str(value).strip().lower()
    if host not in VALID_HOSTS:
        raise ValueError(
            f"Invalid host value for {context}: {value!r}. "
            f"Expected one of: {', '.join(VALID_HOSTS)}."
        )
    return host


def _yaml_quote(value: object) -> str:
    text = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{text}"'


def _yaml_unquote(value: str) -> str:
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]
    return value.replace('\\"', '"').replace("\\\\", "\\")
