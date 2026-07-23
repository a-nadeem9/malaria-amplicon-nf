"""Expected KELT inline-barcode map parsing and sample-sheet enrichment."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DNA_RE = re.compile(r"^[ACGT]+$")
SAMPLE_ALIASES = {"sample", "sampleid", "samplename", "id"}
LIBRARY_ALIASES = {"library", "libraryid", "lib"}
FORWARD_ALIASES = {"forward", "forwardbarcode", "fwd", "fwdbarcode", "barcodef"}
REVERSE_ALIASES = {"reverse", "reversebarcode", "rev", "revbarcode", "barcoder"}


@dataclass(frozen=True)
class KeltBarcode:
    sample_id: str
    forward: str
    reverse: str
    library: str = ""
    row: int = 0


def _header_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value or "").strip().lower())


def _value_key(value: str) -> str:
    return str(value or "").strip().casefold()


def _find_column(fieldnames: Iterable[str], aliases: set[str]) -> str:
    for field in fieldnames:
        if _header_key(field) in aliases:
            return field
    return ""


def _delimiter(path: Path) -> str:
    if path.suffix.lower() == ".tsv":
        return "\t"
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        sample = handle.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except csv.Error:
        return ","


def load_kelt_barcode_map(path: Path | str) -> list[KeltBarcode]:
    source = Path(path).expanduser().resolve()
    if source.suffix.lower() not in {".csv", ".tsv"}:
        raise ValueError("KELT barcode map must be a CSV or TSV file")
    if not source.exists():
        raise ValueError(f"KELT barcode map not found: {source}")

    with source.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=_delimiter(source))
        fields = list(reader.fieldnames or [])
        sample_col = _find_column(fields, SAMPLE_ALIASES)
        forward_col = _find_column(fields, FORWARD_ALIASES)
        reverse_col = _find_column(fields, REVERSE_ALIASES)
        library_col = _find_column(fields, LIBRARY_ALIASES)
        missing = [
            label
            for label, column in (
                ("sample_id", sample_col),
                ("forward_barcode", forward_col),
                ("reverse_barcode", reverse_col),
            )
            if not column
        ]
        if missing:
            raise ValueError("KELT barcode map is missing columns: " + ", ".join(missing))

        records: list[KeltBarcode] = []
        seen: dict[tuple[str, str], KeltBarcode] = {}
        for row_number, row in enumerate(reader, start=2):
            sample_id = str(row.get(sample_col, "") or "").strip()
            forward = re.sub(r"\s+", "", str(row.get(forward_col, "") or "")).upper()
            reverse = re.sub(r"\s+", "", str(row.get(reverse_col, "") or "")).upper()
            library = str(row.get(library_col, "") or "").strip() if library_col else ""
            if not any((sample_id, forward, reverse, library)):
                continue
            if not sample_id or not forward or not reverse:
                raise ValueError(f"KELT barcode map row {row_number} has an incomplete barcode pair")
            if not DNA_RE.fullmatch(forward) or not DNA_RE.fullmatch(reverse):
                raise ValueError(f"KELT barcode map row {row_number} contains a non-ACGT barcode")
            record = KeltBarcode(sample_id, forward, reverse, library, row_number)
            key = (_value_key(sample_id), _value_key(library))
            previous = seen.get(key)
            if previous and (previous.forward, previous.reverse) != (forward, reverse):
                suffix = f" in {library}" if library else ""
                raise ValueError(f"KELT barcode map has conflicting pairs for {sample_id}{suffix}")
            if not previous:
                seen[key] = record
                records.append(record)
    if not records:
        raise ValueError("KELT barcode map contains no barcode pairs")
    return records


def enrich_rows_with_kelt_barcodes(
    rows: list[dict[str, str]],
    barcode_map_path: Path | str,
) -> list[dict[str, str]]:
    records = load_kelt_barcode_map(barcode_map_path)
    by_sample_library = {
        (_value_key(record.sample_id), _value_key(record.library)): record
        for record in records
        if record.library
    }
    by_sample: dict[str, list[KeltBarcode]] = {}
    for record in records:
        by_sample.setdefault(_value_key(record.sample_id), []).append(record)

    enriched: list[dict[str, str]] = []
    missing: list[str] = []
    ambiguous: list[str] = []
    for source_row in rows:
        row = dict(source_row)
        library = str(row.get("library", "") or "").strip()
        identities = [
            str(row.get("sample_id", "") or "").strip(),
            str(row.get("biological_sample_id", "") or "").strip(),
        ]
        identities = list(dict.fromkeys(identity for identity in identities if identity))
        match: KeltBarcode | None = None
        for identity in identities:
            match = by_sample_library.get((_value_key(identity), _value_key(library)))
            if match:
                break
        if match is None:
            candidates = {
                candidate
                for identity in identities
                for candidate in by_sample.get(_value_key(identity), [])
            }
            if len(candidates) == 1:
                match = next(iter(candidates))
            elif len(candidates) > 1:
                ambiguous.append(str(row.get("sample_id", "")))
        if match is None:
            if str(row.get("sample_id", "")) not in ambiguous:
                missing.append(str(row.get("sample_id", "")))
        else:
            row["expected_fwd_barcode"] = match.forward
            row["expected_rev_barcode"] = match.reverse
            row["kelt_barcode_source"] = f"row {match.row}"
        enriched.append(row)

    if ambiguous:
        raise ValueError(
            "KELT barcode map has multiple possible rows for: " + ", ".join(ambiguous[:12])
            + ("..." if len(ambiguous) > 12 else "")
        )
    if missing:
        raise ValueError(
            f"KELT barcode map is missing {len(missing)} selected sample(s): "
            + ", ".join(missing[:12])
            + ("..." if len(missing) > 12 else "")
        )
    return enriched


def inspect_kelt_barcode_map(path: Path | str) -> dict[str, object]:
    records = load_kelt_barcode_map(path)
    return {
        "path": str(Path(path).expanduser().resolve()),
        "pairs": len(records),
        "libraries": sorted({record.library for record in records if record.library}),
        "forward_barcodes": len({record.forward for record in records}),
        "reverse_barcodes": len({record.reverse for record in records}),
    }
