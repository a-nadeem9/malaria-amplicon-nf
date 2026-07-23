"""Build reproducible longitudinal inputs from sequencing calls and metadata.

The sequencing pipeline and the visit calendar answer different questions:
sequencing determines which alleles were retained, while metadata determines
which participant visits occurred.  This module joins those sources without
turning an unobserved genotype into a biological absence.
"""

from __future__ import annotations

import calendar as month_calendar
import hashlib
import json
import re
from collections import Counter
from dataclasses import asdict
from datetime import date
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .metadata import (
    NEGATIVE_DETECTION_TOKENS,
    POSITIVE_DETECTION_TOKENS,
    MetadataCatalog,
    MetadataRecord,
    normalize_detection_token,
    normalize_key,
    normalize_metadata_contract,
    normalized_identifier_collisions,
)


MODEL_VISIT_CLASSES = {
    "genotyped",
    "pcr_negative",
    "pcr_positive_genotype_missing",
}
MISSING_DETECTION_VALUES = {"", "absent", "missing", "na", "nan", "nat", "none", "unknown"}


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def normalized_token(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", clean_text(value).lower())


def first_nonempty(values: Iterable[object]) -> str:
    for value in values:
        cleaned = clean_text(value)
        if cleaned:
            return cleaned
    return ""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_detection(
    values: Iterable[object],
    value_map: dict[str, str] | None = None,
) -> tuple[set[str], list[str]]:
    configured = {
        normalize_detection_token(raw): str(state).strip().lower()
        for raw, state in (value_map or {}).items()
        if clean_text(raw)
    }
    states: set[str] = set()
    raw: list[str] = []
    for value in values:
        display = clean_text(value)
        token = normalize_detection_token(display)
        if not display:
            continue
        configured_state = configured.get(token, "")
        if not configured_state and token in MISSING_DETECTION_VALUES:
            continue
        raw.append(display)
        if configured_state == "ignore":
            continue
        if configured_state in {"positive", "negative"}:
            states.add(configured_state)
        elif configured_state == "review":
            states.add("unknown")
        elif token in POSITIVE_DETECTION_TOKENS:
            states.add("positive")
        elif token in NEGATIVE_DETECTION_TOKENS:
            states.add("negative")
        else:
            states.add("unknown")
    return states, sorted(set(raw))


def _fallback_date(record: MetadataRecord, year: str, day: str) -> tuple[str, str]:
    if record.collection_date:
        return record.collection_date, "metadata_exact_date"
    if not year or not record.month:
        return "", "metadata_date_missing"
    try:
        numeric_year = int(year)
        numeric_month = int(record.month)
        numeric_day = int(day)
        last_day = month_calendar.monthrange(numeric_year, numeric_month)[1]
        if numeric_day < 1 or numeric_day > last_day:
            return "", "metadata_fallback_date_invalid"
        return date(numeric_year, numeric_month, numeric_day).isoformat(), "metadata_month_fallback_date"
    except (TypeError, ValueError):
        return "", "metadata_fallback_date_invalid"


def _metadata_visit_rows(
    catalog: MetadataCatalog | None,
    participant_names: dict[str, str],
    fallback_year: str,
    fallback_day: str,
    excluded_status_values: Iterable[str] = (),
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    columns = [
        "participant_id", "collection_date", "metadata_pcr_values", "metadata_status_values",
        "metadata_season", "metadata_rows", "metadata_source_rows", "metadata_date_source",
    ]
    if catalog is None:
        return pd.DataFrame(columns=columns), []

    rows: list[dict[str, object]] = []
    excluded: list[dict[str, object]] = []
    excluded_status_tokens = {
        normalized_token(value) for value in excluded_status_values if clean_text(value)
    }
    for record in catalog.records:
        participant_key = normalize_key(record.participant_id)
        if participant_key not in participant_names:
            continue
        status = record.values.get("metadata_status", "")
        if normalized_token(status) in excluded_status_tokens:
            excluded.append({
                "participant_id": participant_names[participant_key],
                "collection_date": record.collection_date,
                "reason": "administrative_row_not_a_specimen",
                "metadata_row": record.row,
            })
            continue
        collection_date, date_source = _fallback_date(record, fallback_year, fallback_day)
        if not collection_date:
            excluded.append({
                "participant_id": participant_names[participant_key],
                "collection_date": "",
                "reason": date_source,
                "metadata_row": record.row,
            })
            continue
        rows.append({
            "participant_id": participant_names[participant_key],
            "collection_date": collection_date,
            "metadata_pcr": record.values.get("metadata_pcr", ""),
            "metadata_status": status,
            "metadata_season": record.values.get("metadata_season", ""),
            "metadata_row": record.row,
            "metadata_date_source": date_source,
        })
    if not rows:
        return pd.DataFrame(columns=columns), excluded

    source = pd.DataFrame(rows)
    grouped_rows: list[dict[str, object]] = []
    for (participant, collection_date), group in source.groupby(
        ["participant_id", "collection_date"], sort=True, dropna=False
    ):
        grouped_rows.append({
            "participant_id": participant,
            "collection_date": collection_date,
            "metadata_pcr_values": [value for value in group["metadata_pcr"] if clean_text(value)],
            "metadata_status_values": [value for value in group["metadata_status"] if clean_text(value)],
            "metadata_season": first_nonempty(group["metadata_season"]),
            "metadata_rows": int(len(group)),
            "metadata_source_rows": "|".join(str(value) for value in sorted(group["metadata_row"].unique())),
            "metadata_date_source": "|".join(sorted(set(group["metadata_date_source"]))),
        })
    return pd.DataFrame(grouped_rows, columns=columns), excluded


def _biological_samples(samples: pd.DataFrame) -> pd.DataFrame:
    frame = samples.copy()
    frame.columns = [clean_text(column).lower() for column in frame.columns]
    for column in ("sample_id", "participant_id", "collection_date"):
        if column not in frame.columns:
            frame[column] = ""
        frame[column] = frame[column].map(clean_text)
    if "sample_type" not in frame.columns:
        frame["sample_type"] = "sample"
    sample_type = frame["sample_type"].map(clean_text).str.lower()
    return frame.loc[
        sample_type.eq("sample") & frame["participant_id"].ne("") & frame["collection_date"].ne("")
    ].copy()


def build_visit_calendar(
    calls: pd.DataFrame,
    samples: pd.DataFrame,
    catalog: MetadataCatalog | None = None,
    *,
    fallback_year: str = "",
    fallback_day: str = "27",
    metadata_contract: dict[str, Any] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return the full calendar, its review rows, and excluded metadata rows."""
    biological = _biological_samples(samples)
    participant_collisions = normalized_identifier_collisions(
        biological["participant_id"].drop_duplicates()
    )
    if participant_collisions:
        displays = next(iter(participant_collisions.values()))
        raise ValueError(
            "Participant identifiers normalize to the same matching key: "
            f"{', '.join(displays)}. Make these identifiers unambiguous before continuing."
        )
    participant_names = {
        normalize_key(value): value for value in biological["participant_id"].drop_duplicates()
    }
    contract = normalize_metadata_contract(metadata_contract)
    metadata_visits, excluded_metadata = _metadata_visit_rows(
        catalog,
        participant_names,
        fallback_year,
        fallback_day,
        contract["excluded_status_values"],
    )

    for column in ("metadata_pcr", "metadata_status", "metadata_season"):
        if column not in biological.columns:
            biological[column] = ""
    sequencing_rows: list[dict[str, object]] = []
    for (participant, collection_date), group in biological.groupby(
        ["participant_id", "collection_date"], sort=True
    ):
        sequencing_rows.append({
            "participant_id": participant,
            "collection_date": collection_date,
            "sequencing_pcr_values": [value for value in group["metadata_pcr"] if clean_text(value)],
            "sequencing_status_values": [value for value in group["metadata_status"] if clean_text(value)],
            "sequencing_season": first_nonempty(group["metadata_season"]),
            "sequencing_rows": int(len(group)),
        })
    sequencing_visits = pd.DataFrame(sequencing_rows)
    if sequencing_visits.empty:
        sequencing_visits = pd.DataFrame(columns=[
            "participant_id", "collection_date", "sequencing_pcr_values",
            "sequencing_status_values", "sequencing_season", "sequencing_rows",
        ])

    keys = pd.concat(
        [
            metadata_visits[["participant_id", "collection_date"]],
            sequencing_visits[["participant_id", "collection_date"]],
        ],
        ignore_index=True,
    ).drop_duplicates()
    calendar = keys.merge(
        metadata_visits, on=["participant_id", "collection_date"], how="left"
    ).merge(
        sequencing_visits, on=["participant_id", "collection_date"], how="left"
    )

    genotype_visits = calls[["participant_id", "collection_date"]].drop_duplicates().copy()
    genotype_visits["has_genotype"] = True
    calendar = calendar.merge(
        genotype_visits, on=["participant_id", "collection_date"], how="left"
    )
    calendar["has_genotype"] = calendar["has_genotype"].fillna(False).astype(bool)

    visit_classes: list[str] = []
    pcr_states: list[str] = []
    raw_pcr_values: list[str] = []
    notes: list[str] = []
    review_required: list[bool] = []
    for row in calendar.itertuples(index=False):
        metadata_pcr = row.metadata_pcr_values if isinstance(row.metadata_pcr_values, list) else []
        sequencing_pcr = row.sequencing_pcr_values if isinstance(row.sequencing_pcr_values, list) else []
        states, raw = normalize_detection(
            [*metadata_pcr, *sequencing_pcr],
            contract["detection_value_map"],
        )
        if row.has_genotype:
            if "negative" in states or ("positive" in states and "negative" in states):
                visit_class = "conflicting_genotype_detection"
                note = (
                    "A retained genotype conflicts with a negative or mixed detection result; "
                    "this visit is audited and excluded rather than guessed."
                )
                needs_review = True
            elif "unknown" in states:
                visit_class = "genotyped_detection_review"
                note = (
                    "A retained genotype is available, but its detection value is mapped to review; "
                    "confirm the metadata value before analysis."
                )
                needs_review = True
            else:
                visit_class = "genotyped"
                note = "A retained allele genotype is available."
                needs_review = False
        elif "positive" in states and "negative" in states:
            visit_class = "conflicting_pcr"
            note = "Positive and negative detection labels conflict; this visit is excluded rather than guessed."
            needs_review = True
        elif states == {"positive"}:
            visit_class = "pcr_positive_genotype_missing"
            note = "Infection was detected but no retained genotype is available; DINEMITES may impute this visit."
            needs_review = False
        elif states == {"negative"}:
            visit_class = "pcr_negative"
            note = "Infection was not detected; alleles are treated as observed absent."
            needs_review = False
        else:
            visit_class = "unknown_no_genotype"
            note = "No retained genotype or unambiguous detection state is available; this visit is excluded."
            needs_review = True
        visit_classes.append(visit_class)
        pcr_states.append("|".join(sorted(states)))
        raw_pcr_values.append("|".join(raw))
        notes.append(note)
        review_required.append(needs_review)

    calendar["visit_class"] = visit_classes
    calendar["pcr_state"] = pcr_states
    calendar["raw_pcr_values"] = raw_pcr_values
    calendar["audit_note"] = notes
    calendar["review_required"] = review_required
    calendar["source"] = "unknown"
    has_metadata = calendar["metadata_rows"].notna()
    has_sequencing = calendar["sequencing_rows"].notna()
    calendar.loc[has_metadata & has_sequencing, "source"] = "metadata+sequencing"
    calendar.loc[has_metadata & ~has_sequencing, "source"] = "metadata-only"
    calendar.loc[~has_metadata & has_sequencing, "source"] = "sequencing-only"
    calendar["season"] = calendar["metadata_season"].fillna("")
    missing_season = calendar["season"].eq("")
    calendar.loc[missing_season, "season"] = calendar.loc[missing_season, "sequencing_season"].fillna("")
    calendar["status"] = calendar["metadata_status_values"].map(
        lambda value: "|".join(sorted(set(value))) if isinstance(value, list) else ""
    )
    missing_status = calendar["status"].eq("")
    calendar.loc[missing_status, "status"] = calendar.loc[missing_status, "sequencing_status_values"].map(
        lambda value: "|".join(sorted(set(value))) if isinstance(value, list) else ""
    )
    calendar["included_in_dinemites"] = calendar["visit_class"].isin(MODEL_VISIT_CLASSES)

    participants_with_genotypes = set(
        calendar.loc[calendar["has_genotype"], "participant_id"].astype(str)
    )
    calendar["participant_has_genotype"] = calendar["participant_id"].isin(participants_with_genotypes)
    calendar["included_in_dinemites"] &= calendar["participant_has_genotype"]
    included = calendar.loc[calendar["included_in_dinemites"]].copy()
    if not included.empty:
        study_start = pd.to_datetime(included["collection_date"]).min()
        included["time"] = (
            pd.to_datetime(included["collection_date"]) - study_start
        ).dt.days.astype(int)
        calendar = calendar.merge(
            included[["participant_id", "collection_date", "time"]],
            on=["participant_id", "collection_date"], how="left",
        )
    else:
        calendar["time"] = pd.NA
    calendar = calendar.sort_values(["participant_id", "collection_date"]).reset_index(drop=True)
    audit = calendar.loc[
        calendar["review_required"] | ~calendar["participant_has_genotype"]
    ].copy()
    return calendar, audit, pd.DataFrame(excluded_metadata)


def build_dinemites_inputs(
    calls: pd.DataFrame, calendar: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    included = calendar.loc[calendar["included_in_dinemites"] & calendar["time"].notna()].copy()
    observed = calls.merge(
        included[[
            "participant_id", "collection_date", "time", "season", "status", "visit_class",
        ]],
        on=["participant_id", "collection_date"], how="inner", validate="many_to_one",
    ).rename(columns={"participant_id": "subject", "collection_date": "date_full"})
    observed["date_label"] = pd.to_datetime(observed["date_full"]).dt.strftime("%d %b %Y")
    observed = observed[[
        "allele", "time", "subject", "locus", "date_full", "date_label", "season",
        "status", "visit_class",
    ]]

    empty_visits = included.loc[
        included["visit_class"].isin({"pcr_negative", "pcr_positive_genotype_missing"})
    ].copy()
    empty_rows = pd.DataFrame({
        "allele": pd.NA,
        "time": empty_visits["time"].astype(int),
        "subject": empty_visits["participant_id"],
        "locus": first_nonempty(calls["locus"]) if not calls.empty else "locus",
        "date_full": empty_visits["collection_date"],
        "date_label": pd.to_datetime(empty_visits["collection_date"]).dt.strftime("%d %b %Y"),
        "season": empty_visits["season"],
        "status": empty_visits["status"],
        "visit_class": empty_visits["visit_class"],
    })
    result = pd.concat([observed, empty_rows], ignore_index=True)
    result = result.sort_values(["subject", "time", "allele"], na_position="last").reset_index(drop=True)
    result["date_full"] = pd.to_datetime(result["date_full"]).dt.strftime("%Y-%m-%d")
    season_text = result["season"].fillna("").astype(str).str.strip().str.lower()
    result["covariate_season"] = season_text.str.contains("wet|rain", regex=True).astype(int)
    result["covariate_season_missing"] = season_text.eq("").astype(int)

    qpcr_only = included.loc[
        included["visit_class"].eq("pcr_positive_genotype_missing"),
        ["participant_id", "time", "collection_date"],
    ].rename(columns={"participant_id": "subject", "collection_date": "date_full"})
    qpcr_only["time"] = qpcr_only["time"].astype(int)
    qpcr_only["date_full"] = pd.to_datetime(qpcr_only["date_full"]).dt.strftime("%Y-%m-%d")
    return result, qpcr_only


def build_observed_dynamics(calls: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Calculate transparent, non-modelled changes between genotyped visits."""
    visits = calls[["participant_id", "collection_date"]].drop_duplicates()
    host = visits.groupby("participant_id", as_index=False).agg(
        first_genotyped_visit=("collection_date", "min"),
        last_genotyped_visit=("collection_date", "max"),
        genotyped_visits=("collection_date", "nunique"),
    )
    host["genotyped_intervals"] = host["genotyped_visits"] - 1
    lifetimes = calls.groupby(["participant_id", "locus", "allele"], as_index=False).agg(
        first_observed_visit=("collection_date", "min"),
        last_observed_visit=("collection_date", "max"),
        observed_visits=("collection_date", "nunique"),
        max_abundance_pct=("max_abundance_pct", "max"),
    ).merge(host, on="participant_id", how="left", validate="many_to_one")
    lifetimes["acquired_after_first_visit"] = (
        lifetimes["first_observed_visit"] > lifetimes["first_genotyped_visit"]
    )
    lifetimes["cleared_before_last_visit"] = (
        lifetimes["last_observed_visit"] < lifetimes["last_genotyped_visit"]
    )
    turnover = lifetimes.groupby("participant_id", as_index=False).agg(
        acquisitions=("acquired_after_first_visit", "sum"),
        clearances=("cleared_before_last_visit", "sum"),
        distinct_alleles=("allele", "nunique"),
        genotyped_visits=("genotyped_visits", "first"),
        genotyped_intervals=("genotyped_intervals", "first"),
    )
    denominator = turnover["genotyped_intervals"].clip(lower=1)
    turnover["acquisitions_per_interval"] = turnover["acquisitions"] / denominator
    turnover["turnover_per_interval"] = (
        turnover["acquisitions"] + turnover["clearances"]
    ) / denominator
    turnover["eligible_longitudinal"] = turnover["genotyped_intervals"] > 0

    transition_rows: list[dict[str, object]] = []
    for participant, group in calls.groupby("participant_id", sort=True):
        visit_sets = {
            collection_date: set(zip(visit["locus"], visit["allele"]))
            for collection_date, visit in group.groupby("collection_date", sort=True)
        }
        dates = sorted(visit_sets)
        for before_date, after_date in zip(dates, dates[1:]):
            before = visit_sets[before_date]
            after = visit_sets[after_date]
            transition_rows.append({
                "participant_id": participant,
                "from_date": before_date,
                "to_date": after_date,
                "gap_days": (pd.Timestamp(after_date) - pd.Timestamp(before_date)).days,
                "alleles_before": len(before),
                "alleles_after": len(after),
                "retained": len(before & after),
                "acquired": len(after - before),
                "cleared": len(before - after),
            })
    transitions = pd.DataFrame(transition_rows)
    return {
        "host_time_summary": host,
        "allele_lifetimes": lifetimes,
        "turnover_metrics": turnover,
        "adjacent_visit_transitions": transitions,
    }


def longitudinal_summary(calendar: pd.DataFrame, calls: pd.DataFrame) -> dict[str, object]:
    class_counts = Counter(calendar["visit_class"])
    modeled = calendar.loc[calendar["included_in_dinemites"]]
    participant_visits = calls[["participant_id", "collection_date"]].drop_duplicates()
    repeated = participant_visits.groupby("participant_id")["collection_date"].nunique()
    return {
        "schema_version": 1,
        "sequenced_cohort_participants": int(calendar["participant_id"].nunique()),
        "calendar_visits": int(len(calendar)),
        "modeled_visits": int(len(modeled)),
        "genotyped_visits": int(class_counts.get("genotyped", 0)),
        "pcr_negative_visits": int(class_counts.get("pcr_negative", 0)),
        "pcr_positive_genotype_missing_visits": int(
            class_counts.get("pcr_positive_genotype_missing", 0)
        ),
        "excluded_conflicting_visits": int(class_counts.get("conflicting_pcr", 0)),
        "excluded_genotype_detection_conflicts": int(
            class_counts.get("conflicting_genotype_detection", 0)
        ),
        "excluded_genotype_detection_review": int(
            class_counts.get("genotyped_detection_review", 0)
        ),
        "excluded_unknown_visits": int(class_counts.get("unknown_no_genotype", 0)),
        "participants_with_retained_genotypes": int(calls["participant_id"].nunique()),
        "participants_with_2plus_genotyped_visits": int((repeated >= 2).sum()),
        "visit_class_counts": dict(sorted(class_counts.items())),
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def catalog_profile(catalog: MetadataCatalog | None, metadata_path: Path | None) -> dict[str, Any]:
    if catalog is None or metadata_path is None:
        return {"available": False, "source": "", "sha256": "", "mapping": {}}
    return {
        "available": True,
        "source": str(metadata_path),
        "filename": metadata_path.name,
        "sha256": sha256_file(metadata_path),
        "sheet": catalog.sheet,
        "header_row": catalog.header_row,
        "mapping": catalog.columns,
        "records": len(catalog.records),
        "issues": [asdict(issue) for issue in catalog.issues],
    }
