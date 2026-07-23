"""FASTQ discovery and SIMPLseq sample sheet writing."""

from __future__ import annotations

import csv
import hashlib
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from typing import Any
from pathlib import Path

from .metadata import METADATA_FIELDS, enrich_rows_with_metadata
from .kelt import enrich_rows_with_kelt_barcodes
from .pathutils import user_path


MONTH_ALIASES = {
    "jan": "01",
    "january": "01",
    "feb": "02",
    "february": "02",
    "mar": "03",
    "march": "03",
    "apr": "04",
    "april": "04",
    "may": "05",
    "jun": "06",
    "june": "06",
    "jul": "07",
    "july": "07",
    "aug": "08",
    "august": "08",
    "sep": "09",
    "sept": "09",
    "september": "09",
    "oct": "10",
    "october": "10",
    "nov": "11",
    "november": "11",
    "dec": "12",
    "december": "12",
}
DEFAULT_COLLECTION_DAY = "27"
MONTH_PATTERN = (
    r"January|February|March|April|June|July|August|September|October|November|December|"
    r"Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sept|Sep|Oct|Nov|Dec"
)
NON_PARTICIPANT_TOKEN_RE = re.compile(r"^(run|lane|pool|amplicon|l)[0-9A-Za-z]*$", re.IGNORECASE)
NEGATIVE_SAMPLE_RE = re.compile(
    r"(^|[^a-z0-9])(ctrl|control|ntc|negative|neg|blank|no[-_ ]?template)([^a-z0-9]|$)",
    re.IGNORECASE,
)
POSITIVE_SAMPLE_RE = re.compile(
    r"(^|[^a-z0-9])(positive|pos)([^a-z0-9]|$)",
    re.IGNORECASE,
)

SAMPLE_FIELDS = [
    "sample_id",
    "biological_sample_id",
    "fastq_1",
    "fastq_2",
    "sample_type",
    "library",
    "participant_id",
    "collection_date",
    "replicate",
    "collection_date_source",
    "inferred_year",
    "inferred_day",
    "date_note",
]
KELT_FIELDS = ["expected_fwd_barcode", "expected_rev_barcode", "kelt_barcode_source"]

READ_SUFFIXES = [
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    ("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
    ("_R1.fq.gz", "_R2.fq.gz"),
    ("_R1_001.fq.gz", "_R2_001.fq.gz"),
]


@dataclass(frozen=True)
class FastqPair:
    sample_id: str
    fastq_1: Path
    fastq_2: Path
    sample_type: str
    biological_sample_id: str = ""
    library: str = ""
    participant_id: str = ""
    collection_date: str = ""
    collection_date_inferred: bool = False
    replicate: str = ""
    collection_date_source: str = ""
    inferred_year: bool = False
    inferred_day: bool = False
    date_note: str = ""


@dataclass(frozen=True)
class FastqScan:
    fastq_dir: Path
    pairs: list[FastqPair]
    missing_r2: list[str]
    orphan_r2: list[str]
    md5_files: int
    total_fastq_bytes: int
    duplicate_sample_ids: list[str]
    auto_disambiguated_sample_ids: int = 0


def split_read_suffix(name: str) -> tuple[str, str, str] | None:
    for r1_suffix, r2_suffix in READ_SUFFIXES:
        if name.endswith(r1_suffix):
            return name[: -len(r1_suffix)], "R1", r1_suffix
        if name.endswith(r2_suffix):
            return name[: -len(r2_suffix)], "R2", r2_suffix
    return None


def _month_number(value: str) -> str:
    return MONTH_ALIASES.get(value.lower(), "")


def _format_date(year: str, month: str, day: str = "") -> str:
    month_int = int(month)
    if day:
        return f"{int(year):04d}-{month_int:02d}-{int(day):02d}"
    return f"{int(year):04d}-{month_int:02d}"


def _normalize_fallback_collection_year(value: str | int | None) -> str:
    raw = str(value or "").strip()
    if not raw:
        return ""
    if not re.fullmatch(r"(19|20)[0-9]{2}", raw):
        raise ValueError("fallback_collection_year must be a four-digit year between 1900 and 2099")
    return raw


def _normalize_fallback_collection_day(value: str | int | None) -> str:
    raw = str(value or DEFAULT_COLLECTION_DAY).strip()
    if not re.fullmatch(r"[0-9]{1,2}", raw):
        raise ValueError("fallback_collection_day must be an integer from 1 to 31")
    day = int(raw)
    if day < 1 or day > 31:
        raise ValueError("fallback_collection_day must be an integer from 1 to 31")
    return f"{day:02d}"


def _format_date_with_default_day(year: str, month: str, fallback_collection_day: str = DEFAULT_COLLECTION_DAY) -> str:
    return _format_date(year, month, fallback_collection_day)


def _set_collection_date(
    parsed: dict[str, Any],
    collection_date: str,
    source: str,
    *,
    inferred_year: bool = False,
    inferred_day: bool = False,
    note: str = "",
) -> None:
    parsed["collection_date"] = collection_date
    parsed["collection_date_source"] = source
    parsed["inferred_year"] = inferred_year
    parsed["inferred_day"] = inferred_day
    parsed["collection_date_inferred"] = inferred_year or inferred_day
    parsed["date_note"] = note


def _mark_month_missing_year(parsed: dict[str, Any]) -> None:
    parsed["collection_date_source"] = "filename_month_missing_year"
    parsed["date_note"] = "Filename provided month only; add metadata or set a fallback year to write collection_date."


def _looks_like_date_token(token: str) -> bool:
    return bool(
        re.fullmatch(rf"({MONTH_PATTERN})[0-9]{{4}}", token, re.IGNORECASE)
        or re.fullmatch(rf"[0-9]{{4}}({MONTH_PATTERN})", token, re.IGNORECASE)
        or re.fullmatch(r"[0-9]{6,8}", token)
    )


def infer_sample_type(sample_id: str) -> str:
    if POSITIVE_SAMPLE_RE.search(sample_id):
        return "positive"
    if NEGATIVE_SAMPLE_RE.search(sample_id):
        return "negative"
    return "sample"


def derive_biological_sample_id(label: str) -> str:
    """Remove a technical-replicate token while preserving the sample label."""
    replicate = re.search(r"Rep(?:licate)?[-_ .]*[0-9]+[A-Za-z]?", label, re.IGNORECASE)
    if not replicate:
        return label
    before = label[:replicate.start()].rstrip("-_ .")
    after = label[replicate.end():].lstrip("-_ .")
    biological_id = "-".join(part for part in (before, after) if part)
    return biological_id or label


def parse_label_metadata(
    label: str,
    *,
    fallback_collection_year: str | int | None = "",
    fallback_collection_day: str | int | None = DEFAULT_COLLECTION_DAY,
) -> dict[str, Any]:
    fallback_collection_year = _normalize_fallback_collection_year(fallback_collection_year)
    fallback_collection_day = _normalize_fallback_collection_day(fallback_collection_day)
    parsed = {
        "participant_id": "",
        "collection_date": "",
        "collection_date_inferred": False,
        "replicate": "",
        "collection_date_source": "",
        "inferred_year": False,
        "inferred_day": False,
        "date_note": "",
    }

    compact = re.match(
        rf"^(?P<participant>[A-Za-z]+[0-9]+)(?P<month>{MONTH_PATTERN})(?P<year>[0-9]{{4}})(?P<replicate>Rep[0-9A-Za-z]+)$",
        label,
        re.IGNORECASE,
    )
    if compact:
        parsed["participant_id"] = compact.group("participant")
        _set_collection_date(
            parsed,
            _format_date_with_default_day(
                compact.group("year"),
                _month_number(compact.group("month")),
                fallback_collection_day,
            ),
            f"filename_month_year_inferred_day_{int(fallback_collection_day)}",
            inferred_day=True,
            note=f"Filename provided year and month; inferred day {int(fallback_collection_day)}.",
        )
        parsed["replicate"] = compact.group("replicate")
        return parsed

    replicate = re.search(r"(Rep(?:licate)?[-_ .]*[0-9]+[A-Za-z]?)", label, re.IGNORECASE)
    if replicate:
        parsed["replicate"] = re.sub(
            r"^Replicate",
            "Rep",
            re.sub(r"[-_ .]+", "", replicate.group(1)),
            flags=re.IGNORECASE,
        )

    month_year = re.search(
        rf"(?P<month>{MONTH_PATTERN})[-_ .]*(?P<year>[0-9]{{4}})",
        label,
        re.IGNORECASE,
    )
    year_month = re.search(
        rf"(?P<year>[0-9]{{4}})[-_ .]*(?P<month>{MONTH_PATTERN})",
        label,
        re.IGNORECASE,
    )
    iso_date = re.search(
        r"(?P<year>20[0-9]{2}|19[0-9]{2})[-_ .](?P<month>[0-9]{1,2})(?:[-_ .](?P<day>[0-9]{1,2}))?",
        label,
    )
    compact_date = re.search(r"(?P<year>20[0-9]{2}|19[0-9]{2})(?P<month>[0-9]{2})(?P<day>[0-9]{2})", label)
    month_only = re.search(rf"(^|[^A-Za-z0-9])(?P<month>{MONTH_PATTERN})([^A-Za-z0-9]|$)", label, re.IGNORECASE)

    if month_year:
        _set_collection_date(
            parsed,
            _format_date_with_default_day(
                month_year.group("year"),
                _month_number(month_year.group("month")),
                fallback_collection_day,
            ),
            f"filename_month_year_inferred_day_{int(fallback_collection_day)}",
            inferred_day=True,
            note=f"Filename provided year and month; inferred day {int(fallback_collection_day)}.",
        )
    elif year_month:
        _set_collection_date(
            parsed,
            _format_date_with_default_day(
                year_month.group("year"),
                _month_number(year_month.group("month")),
                fallback_collection_day,
            ),
            f"filename_year_month_inferred_day_{int(fallback_collection_day)}",
            inferred_day=True,
            note=f"Filename provided year and month; inferred day {int(fallback_collection_day)}.",
        )
    elif iso_date:
        if iso_date.group("day"):
            _set_collection_date(
                parsed,
                _format_date(iso_date.group("year"), iso_date.group("month"), iso_date.group("day")),
                "filename_exact_date",
            )
        else:
            _set_collection_date(
                parsed,
                _format_date_with_default_day(
                    iso_date.group("year"),
                    iso_date.group("month"),
                    fallback_collection_day,
                ),
                f"filename_year_month_inferred_day_{int(fallback_collection_day)}",
                inferred_day=True,
                note=f"Filename provided year and month; inferred day {int(fallback_collection_day)}.",
            )
    elif compact_date:
        _set_collection_date(
            parsed,
            _format_date(compact_date.group("year"), compact_date.group("month"), compact_date.group("day")),
            "filename_exact_date",
        )
    elif month_only:
        if fallback_collection_year:
            _set_collection_date(
                parsed,
                _format_date_with_default_day(
                    fallback_collection_year,
                    _month_number(month_only.group("month")),
                    fallback_collection_day,
                ),
                f"filename_month_inferred_year_{fallback_collection_year}_day_{int(fallback_collection_day)}",
                inferred_year=True,
                inferred_day=True,
                note=(
                    "Filename provided month only; inferred year "
                    f"{fallback_collection_year} and day {int(fallback_collection_day)}."
                ),
            )
        else:
            _mark_month_missing_year(parsed)

    tokens = [token for token in re.split(r"[^A-Za-z0-9]+", label) if token]
    for token in tokens:
        if re.fullmatch(rf"{MONTH_PATTERN}", token, re.IGNORECASE):
            continue
        if re.fullmatch(r"Rep(?:licate)?[0-9A-Za-z]*", token, re.IGNORECASE):
            continue
        if re.fullmatch(r"(20[0-9]{2}|19[0-9]{2}|[0-9]{1,2})", token):
            continue
        if _looks_like_date_token(token):
            continue
        if NON_PARTICIPANT_TOKEN_RE.fullmatch(token):
            continue
        if re.search(r"[A-Za-z]", token) and re.search(r"[0-9]", token):
            parsed["participant_id"] = token
            break

    return parsed


def _mpg_sample_wrapper(stripped: str) -> dict[str, str]:
    patterns = (
        (
            re.compile(
                r"^mpg_(?P<run>[^_]+)_Lib(?P<group>[0-9]+)-Pool-"
                r"(?P<label>.+?)(?:_S(?P<index>[0-9]+))?$",
                re.IGNORECASE,
            ),
            "Lib",
        ),
        (
            re.compile(
                r"^mpg_(?P<run>[^_]+)_Amplicon-Pool-(?P<group>[0-9]+)-"
                r"(?P<label>.+?)(?:_S(?P<index>[0-9]+))?$",
                re.IGNORECASE,
            ),
            "Pool",
        ),
        (
            re.compile(
                r"^mpg_(?P<run>[^_]+)_LIB-(?P<group>[0-9]+)-"
                r"(?P<label>.+?)(?:_S(?P<index>[0-9]+))?$",
                re.IGNORECASE,
            ),
            "Lib",
        ),
    )
    for pattern, group_prefix in patterns:
        match = pattern.match(stripped)
        if match:
            return {
                "label": match.group("label"),
                "source_group": f"{group_prefix}{match.group('group')}",
                "library": f"LIB-{match.group('group')}" if group_prefix == "Lib" else f"Pool-{match.group('group')}",
                "source_run": match.group("run"),
                "sequencing_index": match.group("index") or "",
            }
    return {"label": stripped, "source_group": "", "library": "", "source_run": "", "sequencing_index": ""}


def parse_fastq_name(
    name: str,
    include_pool_in_sample_id: bool = False,
    *,
    fallback_collection_year: str | int | None = "",
    fallback_collection_day: str | int | None = DEFAULT_COLLECTION_DAY,
) -> dict[str, Any]:
    base = os.path.basename(name)
    read_parts = split_read_suffix(base)
    stripped = read_parts[0] if read_parts else re.sub(r"_R[12](?:_001)?\.f(?:ast)?q\.gz$", "", base)
    parsed = {
        "sample_id": stripped,
        "biological_sample_id": stripped,
        "participant_id": "",
        "collection_date": "",
        "collection_date_inferred": False,
        "replicate": "",
        "collection_date_source": "",
        "inferred_year": False,
        "inferred_day": False,
        "date_note": "",
        "source_group": "",
        "library": "",
        "source_run": "",
        "sequencing_index": "",
    }
    wrapper = _mpg_sample_wrapper(stripped)
    label = wrapper["label"]
    parsed.update(wrapper)
    parsed["biological_sample_id"] = derive_biological_sample_id(label)
    parsed["sample_id"] = (
        f"{label}_{wrapper['source_group']}"
        if include_pool_in_sample_id and wrapper["source_group"]
        else label
    )

    parsed.update(
        parse_label_metadata(
            label,
            fallback_collection_year=fallback_collection_year,
            fallback_collection_day=fallback_collection_day,
        )
    )
    return parsed


def _duplicate_ids(pairs: list[FastqPair]) -> list[str]:
    return sorted(
        sample_id for sample_id, count in Counter(pair.sample_id for pair in pairs).items() if count > 1
    )


def _disambiguate_sample_ids(pairs: list[FastqPair]) -> tuple[list[FastqPair], int]:
    duplicate_bases = set(_duplicate_ids(pairs))
    if not duplicate_bases:
        return pairs, 0

    parsed_pairs = [
        (pair, parse_fastq_name(pair.fastq_1.name))
        for pair in pairs
    ]
    tag_levels = (
        ("source_group",),
        ("source_group", "source_run"),
        ("source_group", "source_run", "sequencing_index"),
    )
    candidate_pairs = pairs
    for fields in tag_levels:
        candidate_pairs = []
        for pair, parsed in parsed_pairs:
            base = pair.sample_id
            if base not in duplicate_bases:
                candidate_pairs.append(pair)
                continue
            tags = [str(parsed.get(field) or "").strip() for field in fields]
            tags = [tag for tag in tags if tag]
            candidate = "_".join([base, *tags]) if tags else base
            candidate_pairs.append(replace(pair, sample_id=candidate, sample_type=infer_sample_type(candidate)))
        if not _duplicate_ids(candidate_pairs):
            changed = sum(left.sample_id != right.sample_id for left, right in zip(pairs, candidate_pairs))
            return candidate_pairs, changed

    # Unusual wrappers can still collide. A short deterministic file-derived tag
    # preserves every library without asking the user to resolve sequencing names.
    final_pairs: list[FastqPair] = []
    for index, (pair, parsed) in enumerate(parsed_pairs):
        base = pair.sample_id
        if base not in duplicate_bases:
            final_pairs.append(pair)
            continue
        digest = hashlib.sha1(pair.fastq_1.name.encode("utf-8")).hexdigest()[:8]
        candidate = f"{candidate_pairs[index].sample_id}_{digest}"
        final_pairs.append(replace(pair, sample_id=candidate, sample_type=infer_sample_type(candidate)))
    changed = sum(left.sample_id != right.sample_id for left, right in zip(pairs, final_pairs))
    return final_pairs, changed


def _safe_exists(path: Path) -> bool:
    try:
        return path.exists()
    except OSError:
        return False


def _safe_file_records(root: Path) -> list[tuple[Path, int]]:
    files: list[tuple[Path, int]] = []
    try:
        with os.scandir(root) as entries:
            for entry in entries:
                try:
                    if entry.is_file(follow_symlinks=False):
                        files.append((Path(entry.path), entry.stat(follow_symlinks=False).st_size))
                except OSError:
                    continue
    except OSError:
        return files
    return files


def scan_fastqs(
    fastq_dir: Path | str,
    *,
    include_pool_in_sample_id: bool = False,
    fallback_collection_year: str | int | None = "",
    fallback_collection_day: str | int | None = DEFAULT_COLLECTION_DAY,
    libraries: list[str] | tuple[str, ...] | set[str] | None = None,
) -> FastqScan:
    fallback_collection_year = _normalize_fallback_collection_year(fallback_collection_year)
    fallback_collection_day = _normalize_fallback_collection_day(fallback_collection_day)
    try:
        root = user_path(fastq_dir).resolve()
    except OSError:
        root = user_path(fastq_dir).expanduser()
    if not _safe_exists(root):
        return FastqScan(root, [], [], [], 0, 0, [])
    file_records = _safe_file_records(root)
    files = [path for path, _size in file_records]
    r1_candidates: dict[str, list[tuple[Path, str]]] = defaultdict(list)
    r2_candidates: dict[str, list[tuple[Path, str]]] = defaultdict(list)
    total_bytes = 0
    for path, size in file_records:
        read_parts = split_read_suffix(path.name)
        if not read_parts:
            continue
        total_bytes += size
        prefix, read, suffix = read_parts
        if read == "R1":
            r1_candidates[prefix].append((path, suffix))
        else:
            r2_candidates[prefix].append((path, suffix))

    duplicate_reads: list[str] = []
    for read, candidates in (("R1", r1_candidates), ("R2", r2_candidates)):
        for prefix, matches in sorted(candidates.items()):
            if len(matches) > 1:
                filenames = ", ".join(sorted(path.name for path, _suffix in matches))
                duplicate_reads.append(f"{prefix} {read}: {filenames}")
    if duplicate_reads:
        raise ValueError(
            "Multiple FASTQ files resolve to the same sample/read. "
            "Remove or rename duplicate naming variants: " + "; ".join(duplicate_reads)
        )

    r1 = {prefix: matches[0] for prefix, matches in r1_candidates.items()}
    r2 = {prefix: matches[0] for prefix, matches in r2_candidates.items()}
    pairs: list[FastqPair] = []
    missing_r2: list[str] = []
    for prefix, (f1, _suffix) in sorted(r1.items()):
        if prefix not in r2:
            missing_r2.append(f1.name)
            continue
        f2 = r2[prefix][0]
        parsed = parse_fastq_name(
            f1.name,
            include_pool_in_sample_id=include_pool_in_sample_id,
            fallback_collection_year=fallback_collection_year,
            fallback_collection_day=fallback_collection_day,
        )
        sample_id = parsed["sample_id"]
        pairs.append(
            FastqPair(
                sample_id=sample_id,
                fastq_1=f1,
                fastq_2=f2,
                sample_type=infer_sample_type(sample_id),
                biological_sample_id=str(parsed.get("biological_sample_id", "")) or sample_id,
                library=str(parsed.get("library", "")),
                participant_id=parsed["participant_id"],
                collection_date=parsed["collection_date"],
                collection_date_inferred=bool(parsed["collection_date_inferred"]),
                replicate=parsed["replicate"],
                collection_date_source=parsed["collection_date_source"],
                inferred_year=bool(parsed["inferred_year"]),
                inferred_day=bool(parsed["inferred_day"]),
                date_note=parsed["date_note"],
            )
        )
    orphan_r2 = sorted(path.name for prefix, (path, _suffix) in r2.items() if prefix not in r1)
    selected_libraries = {str(value).strip() for value in (libraries or []) if str(value).strip()}
    if selected_libraries:
        pairs = [pair for pair in pairs if pair.library in selected_libraries]
    pairs, auto_disambiguated = _disambiguate_sample_ids(pairs)
    duplicate_ids = _duplicate_ids(pairs)
    md5_files = sum(1 for p in files if p.name.endswith(".md5"))
    return FastqScan(
        root,
        pairs,
        missing_r2,
        orphan_r2,
        md5_files,
        total_bytes,
        duplicate_ids,
        auto_disambiguated,
    )


def pair_to_row(pair: FastqPair, output_root: Path, absolute: bool) -> dict[str, str]:
    def output_path(path: Path) -> str:
        if absolute:
            return str(path)
        try:
            return os.path.relpath(path, output_root).replace(os.sep, "/")
        except ValueError:
            return str(path)

    if absolute:
        fq1 = output_path(pair.fastq_1)
        fq2 = output_path(pair.fastq_2)
    else:
        fq1 = output_path(pair.fastq_1)
        fq2 = output_path(pair.fastq_2)
    return {
        "sample_id": pair.sample_id,
        "biological_sample_id": pair.biological_sample_id or pair.sample_id,
        "fastq_1": fq1,
        "fastq_2": fq2,
        "sample_type": pair.sample_type,
        "library": pair.library,
        "participant_id": pair.participant_id,
        "collection_date": pair.collection_date,
        "replicate": pair.replicate,
        "collection_date_source": pair.collection_date_source,
        "inferred_year": str(pair.inferred_year).lower(),
        "inferred_day": str(pair.inferred_day).lower(),
        "date_note": pair.date_note,
    }


def write_samples_csv(
    fastq_dir: Path | str,
    output_csv: Path | str,
    *,
    include_pool_in_sample_id: bool = False,
    absolute: bool = False,
    fallback_collection_year: str | int | None = "",
    fallback_collection_day: str | int | None = DEFAULT_COLLECTION_DAY,
    metadata_path: Path | str | None = None,
    metadata_sheet: str = "",
    metadata_date_order: str = "auto",
    metadata_columns: dict[str, str] | None = None,
    kelt_barcode_map: Path | str | None = None,
    libraries: list[str] | tuple[str, ...] | set[str] | None = None,
) -> tuple[int, list[str]]:
    output = user_path(output_csv).resolve()
    scan = scan_fastqs(
        fastq_dir,
        include_pool_in_sample_id=include_pool_in_sample_id,
        fallback_collection_year=fallback_collection_year,
        fallback_collection_day=fallback_collection_day,
        libraries=libraries,
    )
    if scan.duplicate_sample_ids:
        return 0, scan.duplicate_sample_ids
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = [pair_to_row(pair, output.parent, absolute) for pair in scan.pairs]
    rows = enrich_rows_with_metadata(
        rows,
        metadata_path,
        metadata_sheet=metadata_sheet,
        date_order=metadata_date_order,
        column_overrides=metadata_columns,
    )
    if kelt_barcode_map:
        rows = enrich_rows_with_kelt_barcodes(rows, kelt_barcode_map)
    fieldnames = list(SAMPLE_FIELDS)
    if metadata_path:
        fieldnames.extend(field for field in METADATA_FIELDS if field not in fieldnames)
    if kelt_barcode_map:
        fieldnames.extend(field for field in KELT_FIELDS if field not in fieldnames)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return len(scan.pairs), []
