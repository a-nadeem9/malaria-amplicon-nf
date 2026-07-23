#!/usr/bin/env python3
"""Build one filtered, replicate-merged allele table for downstream analyses."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import pandas as pd


RUNTIME_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RUNTIME_ROOT / "src"))

from simplseq.longitudinal import (  # noqa: E402
    build_dinemites_inputs,
    build_observed_dynamics,
    build_visit_calendar,
    catalog_profile,
    longitudinal_summary,
    write_json,
)
from simplseq.metadata import inspect_metadata, normalize_metadata_contract  # noqa: E402
from simplseq.panel import panel_profile  # noqa: E402
from simplseq.samplesheet import derive_biological_sample_id  # noqa: E402


REQUIRED_SAMPLE_COLUMNS = {"sample_id", "participant_id", "collection_date", "replicate"}


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def first_nonempty(values: pd.Series) -> object:
    for value in values:
        if clean_text(value):
            return value
    return ""


def biological_metadata(samples: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    samples = samples.copy()
    samples.columns = [clean_text(column).lower() for column in samples.columns]
    missing = sorted(REQUIRED_SAMPLE_COLUMNS - set(samples.columns))
    if missing:
        raise ValueError(f"samples.csv is missing required columns: {', '.join(missing)}")

    for column in REQUIRED_SAMPLE_COLUMNS:
        samples[column] = samples[column].map(clean_text)
    if samples["sample_id"].eq("").any() or samples["sample_id"].duplicated().any():
        raise ValueError("samples.csv sample_id values must be non-empty and unique.")
    if "biological_sample_id" not in samples.columns:
        samples["biological_sample_id"] = samples["sample_id"].map(derive_biological_sample_id)
        repaired_biological_ids = int(
            samples["biological_sample_id"].ne(samples["sample_id"]).sum()
        )
    else:
        samples["biological_sample_id"] = samples["biological_sample_id"].map(clean_text)
        # Older app releases copied sample_id into biological_sample_id. Repair
        # only those placeholder values; preserve explicit specimen IDs supplied
        # by users even when several specimens share a participant/date.
        placeholder_biological = (
            samples["biological_sample_id"].eq("")
            | samples["biological_sample_id"].eq(samples["sample_id"])
        )
        repaired = samples.loc[placeholder_biological, "sample_id"].map(
            derive_biological_sample_id
        )
        repaired_biological_ids = int(
            repaired.ne(samples.loc[placeholder_biological, "biological_sample_id"]).sum()
        )
        samples.loc[placeholder_biological, "biological_sample_id"] = repaired
    if "library" not in samples.columns:
        samples["library"] = ""
    samples["library"] = samples["library"].map(clean_text)
    if "sample_type" not in samples.columns:
        samples["sample_type"] = "sample"
    sample_type = samples["sample_type"].map(clean_text).str.lower()
    controls = ~sample_type.eq("sample")
    missing_identity = samples["participant_id"].eq("") | samples["collection_date"].eq("")
    missing_biological = samples["biological_sample_id"].eq("")
    eligible = ~controls & ~missing_identity & ~missing_biological & samples["sample_id"].ne("")
    audit = {
        "sample_sheet_rows": int(len(samples)),
        "controls_excluded": int(controls.sum()),
        "unassigned_or_undated_excluded": int((~controls & missing_identity).sum()),
        "missing_biological_sample_id_excluded": int((~controls & ~missing_identity & missing_biological).sum()),
        "biological_sequencing_rows": int(eligible.sum()),
        "legacy_biological_ids_repaired": repaired_biological_ids,
    }
    return samples.loc[eligible].copy(), audit


def parse_allele_columns(columns: list[str]) -> pd.DataFrame:
    rows = []
    for column in columns:
        locus, separator, allele = str(column).partition(",")
        if not separator or not locus or not allele:
            raise ValueError(
                f"Allele column {column!r} is not encoded as LOCUS,ALLELE."
            )
        rows.append({"allele_column": column, "locus": locus, "allele": allele})
    return pd.DataFrame(rows)


def build_analysis_table(args: argparse.Namespace) -> dict[str, object]:
    source = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    if source.empty or len(source.columns) < 2:
        raise ValueError(f"Allele count table is empty: {args.input}")
    sample_column = str(source.columns[0])
    source = source.rename(columns={sample_column: "sample_id"})
    source["sample_id"] = source["sample_id"].map(clean_text)
    if source["sample_id"].eq("").any() or source["sample_id"].duplicated().any():
        raise ValueError("Allele count table sample IDs must be non-empty and unique.")
    allele_columns = [str(column) for column in source.columns if column != "sample_id"]
    parsed = parse_allele_columns(allele_columns)
    for column in allele_columns:
        numeric = pd.to_numeric(source[column], errors="coerce")
        invalid = numeric.isna() | ~numeric.map(math.isfinite) | (numeric < 0) | (numeric % 1 != 0)
        if invalid.any():
            bad_samples = ", ".join(source.loc[invalid, "sample_id"].head(5))
            raise ValueError(
                f"Allele column {column!r} contains non-integer, negative, or missing counts"
                + (f" for: {bad_samples}" if bad_samples else "")
            )
        source[column] = numeric.astype("int64")

    raw_samples = pd.read_csv(args.samples, dtype=str, keep_default_na=False)
    samples, exclusion_audit = biological_metadata(raw_samples)
    known_sample_ids = set(raw_samples["sample_id"].map(clean_text))
    unknown_source_ids = sorted(set(source["sample_id"]) - known_sample_ids)
    if unknown_source_ids:
        raise ValueError(
            "Allele count table contains IDs absent from samples.csv: "
            + ", ".join(unknown_source_ids[:20])
        )
    matched = samples[samples["sample_id"].isin(set(source["sample_id"]))].copy()
    if matched.empty:
        raise ValueError("No biological sample IDs matched between samples.csv and the allele count table.")
    source = source[source["sample_id"].isin(set(matched["sample_id"]))].copy()
    source = source.merge(
        matched[
            [
                "sample_id", "biological_sample_id", "participant_id", "collection_date",
                "replicate", "library",
            ]
        ],
        on="sample_id",
        how="inner",
        validate="one_to_one",
    )

    identity_counts = source.groupby("biological_sample_id").agg(
        participants=("participant_id", "nunique"),
        dates=("collection_date", "nunique"),
    )
    conflicting_identity = identity_counts[
        (identity_counts["participants"] != 1) | (identity_counts["dates"] != 1)
    ]
    if not conflicting_identity.empty:
        raise ValueError(
            "A biological_sample_id must map to exactly one participant and collection date: "
            + ", ".join(conflicting_identity.index.astype(str).tolist()[:20])
        )

    source["replicate_key"] = source["replicate"].map(clean_text)
    source.loc[source["replicate_key"].eq(""), "replicate_key"] = "unlabelled"

    observation_long = source.melt(
        id_vars=[
            "sample_id", "biological_sample_id", "participant_id", "collection_date",
            "replicate_key", "library",
        ],
        value_vars=allele_columns,
        var_name="allele_column",
        value_name="reads",
    ).merge(parsed, on="allele_column", how="left", validate="many_to_one")

    # One technical replicate can be sequenced in more than one library/lane. Sum
    # those observations within each declared replicate before testing replicate
    # agreement and merging the confirmed counts into one biological sample.
    long = (
        observation_long.groupby(
            [
                "biological_sample_id", "participant_id", "collection_date", "replicate_key",
                "locus", "allele", "allele_column",
            ],
            as_index=False,
        )
        .agg(
            reads=("reads", "sum"),
            sequencing_observations=("sample_id", "nunique"),
            libraries=("library", lambda values: ";".join(sorted({clean_text(value) for value in values if clean_text(value)}))),
        )
    )

    biological_replicates = (
        long.groupby(["biological_sample_id", "participant_id", "collection_date"], as_index=False)
        .agg(n_replicates=("replicate_key", "nunique"))
    )

    # Apply the abundance rule within every technical replicate. An allele must pass
    # in every available replicate before its counts can enter the biological call.
    # The 100-read depth requirement belongs to the merged biological sample/locus.
    if args.denominator == "locus":
        long["replicate_filter_total_reads"] = long.groupby(
            ["biological_sample_id", "replicate_key", "locus"]
        )["reads"].transform("sum")
    else:
        long["replicate_filter_total_reads"] = long.groupby(
            ["biological_sample_id", "replicate_key"]
        )["reads"].transform("sum")
    replicate_required = long["replicate_filter_total_reads"] * (args.min_abundance_pct / 100.0)
    long["replicate_min_reads_required"] = replicate_required.map(
        lambda value: max(1, math.ceil(value))
    )
    long["replicate_abundance_pct"] = 0.0
    positive_replicate_total = long["replicate_filter_total_reads"] > 0
    long.loc[positive_replicate_total, "replicate_abundance_pct"] = (
        100.0
        * long.loc[positive_replicate_total, "reads"]
        / long.loc[positive_replicate_total, "replicate_filter_total_reads"]
    )
    long["replicate_pass"] = (
        (long["reads"] > 0)
        & (long["reads"] >= long["replicate_min_reads_required"])
    )

    filter_summary = (
        long.groupby(
            ["biological_sample_id", "participant_id", "collection_date", "replicate_key"],
            as_index=False,
        )
        .agg(total_nonzero_alleles=("reads", lambda values: int((values > 0).sum())),
             passing_filter=("replicate_pass", "sum"),
             sequencing_observations=("sequencing_observations", "max"))
    )
    filter_summary["removed_by_threshold"] = (
        filter_summary["total_nonzero_alleles"] - filter_summary["passing_filter"]
    )

    strict_agreement = (
        long.groupby(
            ["biological_sample_id", "participant_id", "collection_date", "locus", "allele", "allele_column"],
            as_index=False,
        )
        .agg(
            strict_n_replicates=("replicate_key", "nunique"),
            strict_n_present=("replicate_pass", "sum"),
        )
    )

    merged = (
        long.groupby(
            ["biological_sample_id", "participant_id", "collection_date", "locus", "allele", "allele_column"],
            as_index=False,
        )
        .agg(
            n_replicates=("replicate_key", "nunique"),
            n_detected_replicates=("reads", lambda values: int((values > 0).sum())),
            read_count=("reads", "sum"),
        )
    )
    merged = merged.merge(
        strict_agreement,
        on=["biological_sample_id", "participant_id", "collection_date", "locus", "allele", "allele_column"],
        how="left",
        validate="one_to_one",
    )
    merged["locus_total_reads"] = merged.groupby(
        ["biological_sample_id", "locus"]
    )["read_count"].transform("sum")
    if args.denominator == "locus":
        merged["filter_total_reads"] = merged["locus_total_reads"]
    else:
        merged["filter_total_reads"] = merged.groupby("biological_sample_id")[
            "read_count"
        ].transform("sum")
    merged["min_reads_required"] = (
        merged["filter_total_reads"] * (args.min_abundance_pct / 100.0)
    ).map(lambda value: max(1, math.ceil(value)))
    merged["abundance_pct"] = 0.0
    positive_total = merged["filter_total_reads"] > 0
    merged.loc[positive_total, "abundance_pct"] = (
        100.0
        * merged.loc[positive_total, "read_count"]
        / merged.loc[positive_total, "filter_total_reads"]
    )
    merged["passes_depth"] = merged["locus_total_reads"] >= args.min_locus_reads
    merged["passes_abundance"] = (
        merged["passes_depth"]
        & (merged["read_count"] > 0)
        & (merged["read_count"] >= merged["min_reads_required"])
    )
    merged["replicate_confirmed"] = (
        (merged["strict_n_present"] == merged["strict_n_replicates"])
        & (merged["strict_n_present"] > 0)
    )
    merged["passes_analysis_filter"] = (
        merged["passes_abundance"] & merged["replicate_confirmed"]
    )

    major_reads = merged.groupby(["biological_sample_id", "locus"])["read_count"].transform("max")
    merged["major_in_sample"] = (
        merged["passes_depth"]
        & (merged["read_count"] > 0)
        & merged["read_count"].eq(major_reads)
    )
    candidate_calls = merged[merged["passes_analysis_filter"]].copy()
    allele_support = (
        candidate_calls.groupby(["locus", "allele", "allele_column"], as_index=False)
        .agg(
            biological_sample_occurrences=("biological_sample_id", "nunique"),
            major_in_any_sample=("major_in_sample", "max"),
        )
    )
    merged = merged.merge(
        allele_support,
        on=["locus", "allele", "allele_column"],
        how="left",
        validate="many_to_one",
    )
    merged["biological_sample_occurrences"] = (
        merged["biological_sample_occurrences"].fillna(0).astype(int)
    )
    merged["major_in_any_sample"] = merged["major_in_any_sample"].fillna(False).astype(bool)
    merged["passes_recurrence"] = (
        (merged["biological_sample_occurrences"] >= args.min_biological_samples)
        | merged["major_in_any_sample"]
    )
    merged["retained"] = merged["passes_analysis_filter"] & merged["passes_recurrence"]
    merged["retention_reason"] = "not retained"
    merged.loc[
        merged["retained"] & (merged["biological_sample_occurrences"] >= args.min_biological_samples),
        "retention_reason",
    ] = "observed in multiple biological samples"
    merged.loc[
        merged["retained"] & (merged["biological_sample_occurrences"] < args.min_biological_samples),
        "retention_reason",
    ] = "major haplotype in at least one biological sample"

    retained = merged[merged["retained"]].copy()
    canonical = retained.rename(columns={"biological_sample_id": "sample_id"})[
        [
            "sample_id", "participant_id", "collection_date", "locus", "allele", "read_count",
            "abundance_pct", "locus_total_reads", "n_replicates", "n_detected_replicates",
            "replicate_confirmed", "biological_sample_occurrences", "major_in_any_sample",
            "min_reads_required", "retention_reason",
        ]
    ].sort_values(["sample_id", "locus", "allele"])
    # Longitudinal summaries use this stable name because an allele can be
    # represented by several sequencing observations before biological merge.
    canonical["max_abundance_pct"] = canonical["abundance_pct"]

    retained_columns = sorted(retained["allele_column"].unique())
    if retained_columns:
        wide = retained.pivot(index="biological_sample_id", columns="allele_column", values="read_count").fillna(0)
        wide = wide.reindex(columns=retained_columns).astype(int).reset_index().rename(columns={"biological_sample_id": "sample"})
    else:
        wide = pd.DataFrame({"sample": sorted(source["biological_sample_id"].unique())})

    merged_samples = matched.copy()
    aggregation = {
        column: first_nonempty
        for column in merged_samples.columns
        if column not in {"biological_sample_id", "sample_id", "replicate", "library"}
    }
    merged_samples = merged_samples.groupby("biological_sample_id", as_index=False).agg(aggregation)
    merged_samples["sample_id"] = merged_samples["biological_sample_id"]
    merged_samples["replicate"] = "merged"
    merged_samples["library"] = "merged"
    merged_samples["sample_type"] = "sample"

    replicate_summary = merged[
        [
            "biological_sample_id", "participant_id", "collection_date", "locus", "allele",
            "n_replicates", "n_detected_replicates", "replicate_confirmed", "passes_depth",
            "passes_abundance", "passes_recurrence", "retained", "retention_reason",
        ]
    ].rename(columns={"biological_sample_id": "sample_id"})

    for path in [args.output_wide, args.output_long, args.output_samples, args.filter_summary, args.replicate_summary, args.summary]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(args.output_wide, sep="\t", index=False)
    canonical.to_csv(args.output_long, sep="\t", index=False)
    merged_samples.to_csv(args.output_samples, index=False)
    filter_summary.to_csv(args.filter_summary, sep="\t", index=False)
    replicate_summary.to_csv(args.replicate_summary, sep="\t", index=False)

    summary: dict[str, object] = {
        "schema_version": 3,
        "mode": args.mode,
        "replicate_policy": "strict",
        "source_table": str(Path(args.input).resolve()),
        "biological_call_profile": {
            "minimum_sample_locus_reads": args.min_locus_reads,
            "minimum_abundance_percent": args.min_abundance_pct,
            "abundance_denominator": args.denominator,
            "recurrence_rule": (
                f"observed in at least {args.min_biological_samples} biological samples "
                "or major haplotype in at least one biological sample"
            ),
        },
        "abundance_filter": {
            "threshold_percent": args.min_abundance_pct,
            "denominator": args.denominator,
            "minimum_sample_locus_reads": args.min_locus_reads,
            "positive_calls_before": int((merged["read_count"] > 0).sum()),
            "positive_calls_after": int(merged["passes_analysis_filter"].sum()),
        },
        "recurrence_filter": {
            "minimum_biological_samples": args.min_biological_samples,
            "major_haplotype_rescue": True,
            "calls_before": int(merged["passes_analysis_filter"].sum()),
            "calls_after": int(merged["retained"].sum()),
            "alleles_retained_by_major_rescue": int(
                merged.loc[
                    merged["retained"]
                    & (merged["biological_sample_occurrences"] < args.min_biological_samples),
                    ["locus", "allele"],
                ].drop_duplicates().shape[0]
            ),
        },
        "exclusions": exclusion_audit,
        "replicate_merge": {
            "observation_rule": "sum repeated library/lane observations within each declared technical replicate",
            "rule": "require the abundance-filtered allele in every technical replicate",
            "count_rule": "sum read counts across technical replicates",
            "biological_samples": int(len(biological_replicates)),
            "participant_visits": int(
                biological_replicates[["participant_id", "collection_date"]].drop_duplicates().shape[0]
            ),
            "visits_with_multiple_replicates": int((biological_replicates["n_replicates"] > 1).sum()),
            "maximum_replicates_per_visit": int(biological_replicates["n_replicates"].max()),
        },
        "output": {
            "participants": int(canonical["participant_id"].nunique()) if not canonical.empty else 0,
            "participant_visits": int(
                canonical[["participant_id", "collection_date"]].drop_duplicates().shape[0]
            ) if not canonical.empty else 0,
            "biological_samples": int(canonical["sample_id"].nunique()) if not canonical.empty else 0,
            "loci": int(canonical["locus"].nunique()) if not canonical.empty else 0,
            "alleles": int(canonical[["locus", "allele"]].drop_duplicates().shape[0]) if not canonical.empty else 0,
            "allele_calls": int(len(canonical)),
        },
        "panel": panel_profile(canonical["locus"] if not canonical.empty else []),
        "files": {
            "wide": str(Path(args.output_wide).resolve()),
            "long": str(Path(args.output_long).resolve()),
            "samples": str(Path(args.output_samples).resolve()),
            "filter_summary": str(Path(args.filter_summary).resolve()),
            "replicate_summary": str(Path(args.replicate_summary).resolve()),
        },
    }
    longitudinal_dir_raw = clean_text(getattr(args, "longitudinal_dir", ""))
    if longitudinal_dir_raw:
        longitudinal_dir = Path(longitudinal_dir_raw).resolve()
        longitudinal_dir.mkdir(parents=True, exist_ok=True)
        metadata_profile: dict[str, object] = {}
        metadata_profile_raw = clean_text(getattr(args, "metadata_profile", ""))
        if metadata_profile_raw:
            metadata_profile_path = Path(metadata_profile_raw).resolve()
            if metadata_profile_path.is_file():
                metadata_profile = json.loads(metadata_profile_path.read_text(encoding="utf-8"))
        metadata_contract = normalize_metadata_contract(metadata_profile.get("contract", {}))
        metadata_raw = clean_text(getattr(args, "metadata", ""))
        metadata_path = Path(metadata_raw).resolve() if metadata_raw else None
        catalog = None
        if metadata_path is not None:
            catalog = inspect_metadata(
                metadata_path,
                clean_text(getattr(args, "metadata_sheet", "")),
                date_order=clean_text(getattr(args, "metadata_date_order", "auto")) or "auto",
                column_overrides=metadata_contract["columns"],
            )
            fatal = [issue.message for issue in catalog.issues if issue.severity == "error"]
            if fatal:
                raise ValueError(" ".join(fatal))

        calendar, visit_audit, excluded_metadata = build_visit_calendar(
            canonical,
            raw_samples,
            catalog,
            fallback_year=clean_text(getattr(args, "fallback_year", "")),
            fallback_day=clean_text(getattr(args, "fallback_day", "27")) or "27",
            metadata_contract=metadata_contract,
        )
        dinemites_input, qpcr_only = build_dinemites_inputs(canonical, calendar)
        dynamics = build_observed_dynamics(canonical) if not canonical.empty else {
            "host_time_summary": pd.DataFrame(),
            "allele_lifetimes": pd.DataFrame(),
            "turnover_metrics": pd.DataFrame(),
            "adjacent_visit_transitions": pd.DataFrame(),
        }

        longitudinal_files = {
            "visit_calendar": longitudinal_dir / "visit_calendar.tsv",
            "visit_audit": longitudinal_dir / "visit_audit.tsv",
            "excluded_metadata_rows": longitudinal_dir / "excluded_metadata_rows.tsv",
            "dinemites_input": longitudinal_dir / "dinemites_observed_input.tsv",
            "qpcr_only": longitudinal_dir / "qpcr_positive_genotype_missing.tsv",
        }
        calendar.to_csv(longitudinal_files["visit_calendar"], sep="\t", index=False)
        visit_audit.to_csv(longitudinal_files["visit_audit"], sep="\t", index=False)
        excluded_metadata.to_csv(longitudinal_files["excluded_metadata_rows"], sep="\t", index=False)
        dinemites_input.to_csv(longitudinal_files["dinemites_input"], sep="\t", index=False)
        qpcr_only.to_csv(longitudinal_files["qpcr_only"], sep="\t", index=False)
        for name, table in dynamics.items():
            target = longitudinal_dir / f"{name}.tsv"
            table.to_csv(target, sep="\t", index=False)
            longitudinal_files[name] = target

        longitudinal_payload = longitudinal_summary(calendar, canonical)
        longitudinal_payload["panel"] = summary["panel"]
        longitudinal_payload["metadata"] = catalog_profile(catalog, metadata_path)
        longitudinal_payload["metadata"]["contract"] = metadata_contract
        longitudinal_payload["files"] = {
            name: str(path) for name, path in longitudinal_files.items()
        }
        longitudinal_payload["method"] = {
            "negative_visit_rule": "PCR-negative visits are observed allele absences.",
            "missing_genotype_rule": (
                "PCR-positive visits without a retained genotype are marked for DINEMITES imputation."
            ),
            "unknown_visit_rule": "Ambiguous visits are audited and excluded rather than guessed.",
            "observed_dynamics_rule": (
                "Acquisition and clearance are calculated from first and last observed allele visits; "
                "rates use the number of intervals between genotyped visits."
            ),
        }
        write_json(longitudinal_dir / "longitudinal_summary.json", longitudinal_payload)
        analysis_profile = {
            "schema_version": 2,
            "analysis_mode": args.mode,
            "allele_definition": "cdhit_cluster_sum" if args.mode == "cdhit98" else "exact_cigar",
            "abundance_threshold_percent": args.min_abundance_pct,
            "abundance_denominator": args.denominator,
            "minimum_sample_locus_reads": args.min_locus_reads,
            "recurrence_rule": "two_biological_samples_or_major_in_one",
            "replicate_rule": "strict",
            "count_rule": "sum_read_counts_across_technical_replicates",
            "metadata": longitudinal_payload["metadata"],
            "longitudinal": longitudinal_payload,
            "panel": summary["panel"],
        }
        profile_path = Path(
            clean_text(getattr(args, "analysis_profile", ""))
            or str(longitudinal_dir / "analysis_profile.json")
        ).resolve()
        write_json(profile_path, analysis_profile)
        summary["longitudinal"] = longitudinal_payload
        summary["files"]["analysis_profile"] = str(profile_path)

    Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--input", required=True)
    result.add_argument("--samples", required=True)
    result.add_argument("--output-wide", required=True)
    result.add_argument("--output-long", required=True)
    result.add_argument("--output-samples", required=True)
    result.add_argument("--filter-summary", required=True)
    result.add_argument("--replicate-summary", required=True)
    result.add_argument("--summary", required=True)
    result.add_argument("--mode", choices=("primary", "cdhit98"), default="primary")
    result.add_argument("--min-abundance-pct", type=float, default=1.0)
    result.add_argument("--min-locus-reads", type=int, default=100)
    result.add_argument("--min-biological-samples", type=int, default=2)
    result.add_argument("--denominator", choices=("locus", "sample"), default="locus")
    result.add_argument("--metadata", default="")
    result.add_argument("--metadata-sheet", default="")
    result.add_argument("--metadata-date-order", choices=("auto", "mdy", "dmy"), default="auto")
    result.add_argument("--metadata-profile", default="")
    result.add_argument("--fallback-year", default="")
    result.add_argument("--fallback-day", default="27")
    result.add_argument("--longitudinal-dir", default="")
    result.add_argument("--analysis-profile", default="")
    return result


def main() -> int:
    args = parser().parse_args()
    if not 0 <= args.min_abundance_pct <= 100:
        raise SystemExit("--min-abundance-pct must be between 0 and 100")
    if args.min_locus_reads < 0:
        raise SystemExit("--min-locus-reads must be non-negative")
    if args.min_biological_samples < 1:
        raise SystemExit("--min-biological-samples must be at least 1")
    summary = build_analysis_table(args)
    print(json.dumps(summary["output"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
