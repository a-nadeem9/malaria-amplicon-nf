#!/usr/bin/env python3
"""Cluster filter-pass nucleotide ASVs with CD-HIT-EST and aggregate counts."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


CLUSTER_HEADER = re.compile(r"^>Cluster\s+(\d+)\s*$")
MEMBER_ID = re.compile(r">(.+?)\.\.\.")
MEMBER_LENGTH = re.compile(r"\s(\d+)nt,")
MEMBER_IDENTITY = re.compile(r"([0-9]+(?:\.[0-9]+)?)%")


def truthy(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "t", "yes", "y"}


def numeric(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_asv_sequences(path: Path) -> dict[str, str]:
    """Read a FASTA catalog without allowing duplicate feature identifiers."""
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"ASV FASTA is empty: {path}")
    identifiers = [str(record.id).strip() for record in records]
    duplicate_ids = sorted(identifier for identifier, count in Counter(identifiers).items() if count > 1)
    if duplicate_ids:
        raise ValueError(
            "ASV FASTA has duplicate feature IDs: " + ", ".join(duplicate_ids[:5])
        )
    if any(not identifier for identifier in identifiers):
        raise ValueError(f"ASV FASTA has an empty feature ID: {path}")
    return {
        identifier: str(record.seq).upper()
        for identifier, record in zip(identifiers, records)
    }


def read_count_table(path: Path, eligible_samples: set[str] | None = None) -> pd.DataFrame:
    """Read a count table while preserving and validating its original headers."""
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    if len(header) < 2:
        raise ValueError(f"ASV count table has no feature columns: {path}")
    blank_columns = [index + 1 for index, label in enumerate(header) if not str(label).strip()]
    if blank_columns:
        raise ValueError(
            "ASV count table has blank column names at position(s): "
            + ", ".join(str(index) for index in blank_columns[:5])
        )
    duplicate_columns = sorted(label for label, count in Counter(header).items() if count > 1)
    if duplicate_columns:
        raise ValueError(
            "ASV count table has duplicate feature columns: "
            + ", ".join(duplicate_columns[:5])
        )

    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if table.empty:
        raise ValueError(f"ASV count table is empty: {path}")
    table = table.rename(columns={str(table.columns[0]): "sample"})
    table["sample"] = table["sample"].astype(str)
    if eligible_samples is not None:
        table = table[table["sample"].isin(eligible_samples)].copy()
    if table.empty:
        raise ValueError("No biological sample rows remain for CD-HIT clustering.")

    feature_columns = list(table.columns[1:])
    numeric_counts = table[feature_columns].apply(pd.to_numeric, errors="coerce")
    invalid = numeric_counts.isna() | ~np.isfinite(numeric_counts)
    negative = numeric_counts < 0
    fractional = numeric_counts != numeric_counts.round()
    if invalid.to_numpy().any() or negative.to_numpy().any() or fractional.to_numpy().any():
        bad = invalid | negative | fractional
        row_index, column_index = np.argwhere(bad.to_numpy())[0]
        sample = table.iloc[int(row_index)]["sample"]
        column = feature_columns[int(column_index)]
        value = table.iloc[int(row_index)][column]
        raise ValueError(
            "ASV count table contains a non-finite, negative, or non-integer count "
            f"for sample {sample!r}, feature {column!r}: {value!r}"
        )
    table[feature_columns] = numeric_counts.astype("int64")
    return table


def resolve_feature_columns(
    feature_columns: list[str],
    asv_sequences: dict[str, str],
) -> dict[str, str]:
    """Resolve every table feature to one, and only one, FASTA feature."""
    sequence_to_ids: dict[str, list[str]] = defaultdict(list)
    for asv_id, sequence in asv_sequences.items():
        sequence_to_ids[str(sequence).upper()].append(str(asv_id))

    resolved: dict[str, str] = {}
    unresolved: list[str] = []
    ambiguous: dict[str, list[str]] = {}
    for column in feature_columns:
        label = str(column)
        candidates: set[str] = set()
        if label in asv_sequences:
            candidates.add(label)
        candidates.update(sequence_to_ids.get(label.upper(), []))
        if not candidates:
            unresolved.append(label)
        elif len(candidates) > 1:
            ambiguous[label] = sorted(candidates)
        else:
            resolved[label] = next(iter(candidates))

    columns_by_asv: dict[str, list[str]] = defaultdict(list)
    for column, asv_id in resolved.items():
        columns_by_asv[asv_id].append(column)
    multiply_resolved = {
        asv_id: columns
        for asv_id, columns in columns_by_asv.items()
        if len(columns) > 1
    }

    problems: list[str] = []
    if unresolved:
        problems.append("unresolved: " + ", ".join(sorted(unresolved)[:5]))
    if ambiguous:
        preview = "; ".join(
            f"{column} -> {', '.join(ids)}"
            for column, ids in sorted(ambiguous.items())[:5]
        )
        problems.append("ambiguous: " + preview)
    if multiply_resolved:
        preview = "; ".join(
            f"{asv_id} <- {', '.join(columns)}"
            for asv_id, columns in sorted(multiply_resolved.items())[:5]
        )
        problems.append("multiply resolved: " + preview)
    if problems:
        raise ValueError(
            "Every ASV count-table feature must resolve exactly once against the FASTA ("
            + " | ".join(problems)
            + ")"
        )
    return resolved


def biological_sample_ids(path: Path) -> tuple[set[str], dict[str, int]]:
    samples = pd.read_csv(path, dtype=str).fillna("")
    samples.columns = [str(column).strip().lower() for column in samples.columns]
    required = {"sample_id", "participant_id", "collection_date"}
    missing = sorted(required - set(samples.columns))
    if missing:
        raise ValueError(f"Sample sheet is missing: {', '.join(missing)}")
    sample_type = (
        samples["sample_type"].astype(str).str.strip().str.lower()
        if "sample_type" in samples.columns
        else pd.Series("sample", index=samples.index)
    )
    participant = samples["participant_id"].astype(str).str.strip()
    collection_date = samples["collection_date"].astype(str).str.strip()
    sample_id = samples["sample_id"].astype(str).str.strip()
    biological_id = (
        samples["biological_sample_id"].astype(str).str.strip()
        if "biological_sample_id" in samples.columns
        else sample_id
    )
    controls = ~sample_type.eq("sample")
    missing_identity = participant.eq("") | collection_date.eq("") | biological_id.eq("")
    unassigned = ~controls & missing_identity
    eligible = ~controls & ~missing_identity & sample_id.ne("")
    return set(sample_id[eligible]), {
        "controls_excluded": int(controls.sum()),
        "unassigned_excluded": int(unassigned.sum()),
        "biological_rows": int(eligible.sum()),
        "biological_samples": int(biological_id[eligible].nunique()),
    }


def mapped_asv_rows(
    path: Path,
    *,
    min_total_reads: int,
    min_samples: int,
    exclude_bimeras: bool,
    biological_metrics: dict[str, dict[str, int]] | None = None,
) -> tuple[list[dict[str, str]], int]:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        ref_fields = [field for field in fields if field.startswith("refid_")]
        rows = list(reader)

    if "hapid" not in fields:
        raise ValueError(f"Mapped ASV table has no hapid column: {path}")
    mapped_ids = [str(row.get("hapid", "")).strip() for row in rows]
    duplicate_ids = sorted(asv_id for asv_id, count in Counter(mapped_ids).items() if asv_id and count > 1)
    if duplicate_ids:
        raise ValueError(
            "Mapped ASV table has duplicate feature IDs: " + ", ".join(duplicate_ids[:5])
        )
    if any(not asv_id for asv_id in mapped_ids):
        raise ValueError(f"Mapped ASV table has an empty hapid value: {path}")

    retained: list[dict[str, str]] = []
    mapped_count = 0
    for row in rows:
        locus = next(
            (
                str(row.get(field, "")).strip()
                for field in ref_fields
                if str(row.get(field, "")).strip().upper() not in {"", "NA", "NAN"}
            ),
            "",
        )
        if ref_fields and not locus:
            continue
        mapped_count += 1
        asv_id = str(row.get("hapid", ""))
        metrics = (biological_metrics or {}).get(asv_id, {})
        total_reads = numeric(metrics.get("total_reads", row.get("total_reads")))
        total_samples = numeric(metrics.get("total_samples", row.get("total_samples")))
        if total_reads < min_total_reads:
            continue
        if total_samples < min_samples:
            continue
        if str(row.get("snv_filter", "")).strip().upper() == "FAIL":
            continue
        if str(row.get("indel_filter", "")).strip().upper() == "FAIL":
            continue
        if exclude_bimeras and truthy(row.get("bimera")):
            continue
        row = dict(row)
        row["_locus"] = locus
        row["total_reads"] = str(int(total_reads))
        row["total_samples"] = str(int(total_samples))
        retained.append(row)
    return retained, mapped_count


def write_filtered_fasta(source: Path, retained_ids: set[str], output: Path) -> dict[str, str]:
    catalog = read_asv_sequences(source)
    records = [record for record in SeqIO.parse(source, "fasta") if record.id in retained_ids]
    found = {record.id for record in records}
    missing = sorted(retained_ids - found)
    if missing:
        raise ValueError(f"Filter-pass ASV IDs missing from FASTA: {', '.join(missing[:5])}")
    if not records:
        raise ValueError("No ASV sequences remain after pipeline QC.")
    SeqIO.write(records, output, "fasta")
    return {asv_id: catalog[asv_id] for asv_id in retained_ids}


def parse_clstr(path: Path) -> list[dict[str, object]]:
    members: list[dict[str, object]] = []
    cluster_number: int | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            header = CLUSTER_HEADER.match(line)
            if header:
                cluster_number = int(header.group(1))
                continue
            if not line or cluster_number is None:
                continue
            member_match = MEMBER_ID.search(line)
            length_match = MEMBER_LENGTH.search(line)
            if not member_match or not length_match:
                raise ValueError(f"Could not parse CD-HIT cluster line: {line}")
            is_representative = line.endswith("*")
            identity_match = MEMBER_IDENTITY.search(line)
            members.append(
                {
                    "cluster_number": cluster_number,
                    "cluster_id": f"CDHIT_{cluster_number + 1:04d}",
                    "member_asv": member_match.group(1),
                    "length_nt": int(length_match.group(1)),
                    "identity_pct": 100.0 if is_representative else numeric(identity_match.group(1) if identity_match else ""),
                    "is_representative": is_representative,
                }
            )
    if not members:
        raise ValueError(f"CD-HIT produced no cluster members: {path}")
    return members


def read_cigar_annotations(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {
            str(row.get("ASV", "")): {
                "amplicon": str(row.get("Amplicon", "")),
                "cigar": str(row.get("CIGAR", "")),
            }
            for row in reader
            if row.get("ASV")
        }


def membership_table(
    members: list[dict[str, object]],
    mapped_rows: list[dict[str, str]],
    cigar_annotations: dict[str, dict[str, str]],
) -> pd.DataFrame:
    mapped = {str(row["hapid"]): row for row in mapped_rows}
    representatives = {
        int(member["cluster_number"]): str(member["member_asv"])
        for member in members
        if bool(member["is_representative"])
    }
    rows: list[dict[str, object]] = []
    for member in members:
        asv_id = str(member["member_asv"])
        metadata = mapped.get(asv_id, {})
        annotation = cigar_annotations.get(asv_id, {})
        cluster_number = int(member["cluster_number"])
        rows.append(
            {
                "cluster_id": member["cluster_id"],
                "representative_asv": representatives.get(cluster_number, ""),
                "member_asv": asv_id,
                "identity_pct": member["identity_pct"],
                "is_representative": bool(member["is_representative"]),
                "length_nt": member["length_nt"],
                "amplicon": annotation.get("amplicon") or metadata.get("_locus", ""),
                "cigar": annotation.get("cigar", ""),
                "total_reads": int(numeric(metadata.get("total_reads"))),
                "total_samples": int(numeric(metadata.get("total_samples"))),
            }
        )
    return pd.DataFrame(rows)


def validate_cluster_members(
    members: list[dict[str, object]],
    expected_asv_ids: set[str],
) -> None:
    """Require CD-HIT to assign every retained ASV to exactly one cluster."""
    observed = [str(member.get("member_asv", "")).strip() for member in members]
    duplicate_members = sorted(
        member for member, count in Counter(observed).items() if member and count > 1
    )
    missing = sorted(expected_asv_ids - set(observed))
    unexpected = sorted(set(observed) - expected_asv_ids)
    empty = sum(not member for member in observed)

    representatives: dict[int, int] = Counter(
        int(member["cluster_number"])
        for member in members
        if bool(member.get("is_representative"))
    )
    cluster_numbers = {int(member["cluster_number"]) for member in members}
    invalid_representatives = sorted(
        cluster for cluster in cluster_numbers if representatives.get(cluster, 0) != 1
    )

    problems: list[str] = []
    if duplicate_members:
        problems.append("multiply assigned: " + ", ".join(duplicate_members[:5]))
    if missing:
        problems.append("unresolved retained ASVs: " + ", ".join(missing[:5]))
    if unexpected:
        problems.append("unknown CD-HIT members: " + ", ".join(unexpected[:5]))
    if empty:
        problems.append(f"empty member IDs: {empty}")
    if invalid_representatives:
        problems.append(
            "clusters without exactly one representative: "
            + ", ".join(str(cluster) for cluster in invalid_representatives[:5])
        )
    if problems:
        raise ValueError(
            "CD-HIT membership must resolve every retained ASV exactly once ("
            + " | ".join(problems)
            + ")"
        )


def validate_read_count_conservation(
    samples: pd.Series,
    input_totals: pd.Series,
    output_totals: pd.Series,
) -> tuple[int, int]:
    """Require clustered counts to preserve reads for every row and globally."""
    input_values = input_totals.astype("int64").reset_index(drop=True)
    output_values = output_totals.astype("int64").reset_index(drop=True)
    if len(input_values) != len(output_values):
        raise ValueError(
            "CD-HIT read-count conservation failed: input and output sample counts differ"
        )
    mismatches = input_values != output_values
    if mismatches.any():
        position = int(np.flatnonzero(mismatches.to_numpy())[0])
        sample = str(samples.reset_index(drop=True).iloc[position])
        raise ValueError(
            "CD-HIT per-sample read-count conservation failed for "
            f"{sample!r}: {int(input_values.iloc[position])} input reads != "
            f"{int(output_values.iloc[position])} clustered reads"
        )
    input_global = int(input_values.sum())
    output_global = int(output_values.sum())
    if input_global != output_global:
        raise ValueError(
            "CD-HIT global read-count conservation failed: "
            f"{input_global} input reads != {output_global} clustered reads"
        )
    return input_global, output_global


def aggregate_cluster_counts(
    seqtab_path: Path,
    asv_sequences: dict[str, str],
    membership: pd.DataFrame,
    eligible_samples: set[str] | None = None,
) -> pd.DataFrame:
    seqtab = read_count_table(seqtab_path, eligible_samples)
    resolved = resolve_feature_columns(list(seqtab.columns[1:]), asv_sequences)
    required_membership_columns = {"member_asv", "cluster_id"}
    missing_columns = sorted(required_membership_columns - set(membership.columns))
    if missing_columns:
        raise ValueError(
            "CD-HIT membership table is missing columns: " + ", ".join(missing_columns)
        )
    if membership.empty:
        raise ValueError("CD-HIT membership table is empty")
    member_ids = membership["member_asv"].astype(str).str.strip()
    duplicate_members = sorted(
        member for member, count in Counter(member_ids).items() if member and count > 1
    )
    if duplicate_members:
        raise ValueError(
            "CD-HIT membership multiply resolves ASVs: " + ", ".join(duplicate_members[:5])
        )
    unknown_members = sorted(set(member_ids) - set(asv_sequences))
    if unknown_members:
        raise ValueError(
            "CD-HIT membership contains ASVs absent from the FASTA: "
            + ", ".join(unknown_members[:5])
        )
    columns_by_asv = {asv_id: column for column, asv_id in resolved.items()}
    unresolved_members = sorted(set(member_ids) - set(columns_by_asv))
    if unresolved_members:
        raise ValueError(
            "CD-HIT membership ASVs have no exactly resolved count-table column: "
            + ", ".join(unresolved_members[:5])
        )

    cluster_for_asv = dict(zip(member_ids, membership["cluster_id"].astype(str)))
    cluster_columns: dict[str, list[str]] = {}
    for asv_id, cluster_id in cluster_for_asv.items():
        if not str(cluster_id).strip():
            raise ValueError(f"CD-HIT membership has an empty cluster ID for {asv_id}")
        cluster_columns.setdefault(str(cluster_id), []).append(columns_by_asv[asv_id])

    output = pd.DataFrame({"sample": seqtab["sample"].astype(str)})
    for cluster_id in sorted(cluster_columns):
        output[cluster_id] = seqtab[cluster_columns[cluster_id]].sum(axis=1).astype("int64")

    input_columns = [columns_by_asv[asv_id] for asv_id in member_ids]
    input_totals = seqtab[input_columns].sum(axis=1).astype("int64")
    output_totals = output.iloc[:, 1:].sum(axis=1).astype("int64")
    input_global, output_global = validate_read_count_conservation(
        output["sample"], input_totals, output_totals
    )
    output.attrs["input_cluster_reads"] = input_global
    output.attrs["total_cluster_reads"] = output_global
    return output


def run_clustering(args: argparse.Namespace) -> dict[str, object]:
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    eligible_samples, sample_audit = biological_sample_ids(Path(args.samples))
    all_asv_sequences = read_asv_sequences(Path(args.fasta))
    mapped_rows, mapped_count = mapped_asv_rows(
        Path(args.mapped),
        min_total_reads=args.min_total_reads,
        min_samples=args.min_samples,
        exclude_bimeras=args.exclude_bimeras,
    )
    retained_ids = {str(row["hapid"]) for row in mapped_rows}
    cigar_annotations = read_cigar_annotations(Path(args.asv_to_cigar))
    cigar_ids = set(cigar_annotations)
    if retained_ids != cigar_ids:
        missing_from_cdhit = sorted(cigar_ids - retained_ids)
        extra_in_cdhit = sorted(retained_ids - cigar_ids)
        details = []
        if missing_from_cdhit:
            details.append("missing from CD-HIT input: " + ", ".join(missing_from_cdhit[:10]))
        if extra_in_cdhit:
            details.append("extra in CD-HIT input: " + ", ".join(extra_in_cdhit[:10]))
        raise ValueError(
            "CD-HIT must consume exactly the ASV set represented by the primary CIGAR mapping ("
            + " | ".join(details)
            + ")"
        )
    filtered_fasta = outdir / "filtered_ASVSeqs.fasta"
    asv_sequences = write_filtered_fasta(Path(args.fasta), retained_ids, filtered_fasta)

    representatives = outdir / "cdhit_representatives.fasta"
    raw_clusters = Path(str(representatives) + ".clstr")
    records_by_id = {record.id: record for record in SeqIO.parse(filtered_fasta, "fasta")}
    locus_for_asv = {str(row["hapid"]): str(row.get("_locus", "")).strip() for row in mapped_rows}
    loci = sorted({locus for locus in locus_for_asv.values() if locus})
    if not loci:
        raise ValueError("No locus assignments are available for filter-pass ASVs.")

    members: list[dict[str, object]] = []
    representative_records = []
    raw_cluster_lines: list[str] = []
    commands: list[list[str]] = []
    log_parts: list[str] = []
    global_cluster_number = 0

    with tempfile.TemporaryDirectory(prefix="cdhit_by_locus_", dir=outdir) as temporary:
        temporary_dir = Path(temporary)
        for locus_index, locus in enumerate(loci, start=1):
            locus_ids = sorted(asv_id for asv_id, value in locus_for_asv.items() if value == locus)
            if not locus_ids:
                continue
            locus_slug = re.sub(r"[^A-Za-z0-9]+", "_", locus).strip("_") or f"LOCUS_{locus_index}"
            locus_fasta = temporary_dir / f"{locus_slug}.fasta"
            locus_representatives = temporary_dir / f"{locus_slug}_representatives.fasta"
            SeqIO.write([records_by_id[asv_id] for asv_id in locus_ids], locus_fasta, "fasta")
            command = [
                args.executable,
                "-i",
                str(locus_fasta),
                "-o",
                str(locus_representatives),
                "-c",
                f"{args.identity:.4f}",
                "-n",
                str(args.word_size),
                "-G",
                "1",
                "-g",
                "0",
                "-d",
                "0",
                "-T",
                str(max(1, args.threads)),
                "-M",
                "0",
            ]
            commands.append(command)
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            log_parts.append(
                f"## Locus: {locus}\n$ {shlex.join(command)}\n"
                + (completed.stdout or "")
                + (completed.stderr or "")
            )
            if completed.returncode != 0:
                (outdir / "cdhit.log").write_text("\n".join(log_parts), encoding="utf-8")
                raise RuntimeError(
                    f"CD-HIT-EST failed for locus {locus} with exit code "
                    f"{completed.returncode}. See cdhit.log."
                )

            local_members = parse_clstr(Path(str(locus_representatives) + ".clstr"))
            local_cluster_numbers = sorted({int(member["cluster_number"]) for member in local_members})
            local_to_global: dict[int, int] = {}
            for local_number in local_cluster_numbers:
                local_to_global[local_number] = global_cluster_number
                global_cluster_number += 1

            for member in local_members:
                local_number = int(member["cluster_number"])
                global_number = local_to_global[local_number]
                member["cluster_number"] = global_number
                member["cluster_id"] = f"CDHIT_{locus_slug}_{global_number + 1:04d}"
                members.append(member)

            representative_records.extend(SeqIO.parse(locus_representatives, "fasta"))
            for local_number in local_cluster_numbers:
                global_number = local_to_global[local_number]
                raw_cluster_lines.append(f">Cluster {global_number}")
                cluster_members = [
                    member for member in local_members
                    if int(member["cluster_number"]) == global_number
                ]
                for member_index, member in enumerate(cluster_members):
                    suffix = "*" if bool(member["is_representative"]) else f"at +/{float(member['identity_pct']):.2f}%"
                    raw_cluster_lines.append(
                        f"{member_index}\t{int(member['length_nt'])}nt, "
                        f">{member['member_asv']}... {suffix}"
                    )

    (outdir / "cdhit.log").write_text("\n".join(log_parts), encoding="utf-8")
    validate_cluster_members(members, retained_ids)
    SeqIO.write(representative_records, representatives, "fasta")
    raw_clusters.write_text("\n".join(raw_cluster_lines) + "\n", encoding="utf-8")
    membership = membership_table(members, mapped_rows, cigar_annotations)
    membership_path = outdir / "cdhit_cluster_membership.tsv"
    membership.to_csv(membership_path, sep="\t", index=False)

    counts = aggregate_cluster_counts(
        Path(args.seqtab), all_asv_sequences, membership, eligible_samples=eligible_samples
    )
    counts_path = outdir / "cdhit_cluster_counts.tsv"
    counts.to_csv(counts_path, sep="\t", index=False)

    sizes = membership.groupby("cluster_id", sort=True).size()
    locus_counts = membership.groupby("cluster_id")["amplicon"].nunique()
    summary: dict[str, object] = {
        "method": "CD-HIT-EST",
        "identity_threshold": args.identity,
        "identity_definition": "global identity over the full length of the shorter sequence",
        "word_size": args.word_size,
        "assignment_mode": "first qualifying cluster",
        "clustering_scope": "within_locus",
        "loci_clustered": len(loci),
        "mapped_asvs": mapped_count,
        "input_asvs": len(membership),
        "clusters": int(sizes.size),
        "merged_asvs": int(len(membership) - sizes.size),
        "singleton_clusters": int((sizes == 1).sum()),
        "largest_cluster_size": int(sizes.max()),
        "mixed_amplicon_clusters": int((locus_counts > 1).sum()),
        "samples": int(len(counts)),
        "input_cluster_reads": int(counts.attrs["input_cluster_reads"]),
        "total_cluster_reads": int(counts.attrs["total_cluster_reads"]),
        "read_count_conserved": True,
        "unresolved_asv_columns": [],
        "unresolved_asv_column_count": 0,
        "sample_exclusions": sample_audit,
        "pipeline_qc": {
            "min_total_reads": args.min_total_reads,
            "min_samples": args.min_samples,
            "exclude_bimeras": args.exclude_bimeras,
            "exclude_snv_or_indel_failures": True,
        },
        "commands": commands,
    }
    (outdir / "cdhit_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--mapped", required=True)
    parser.add_argument("--seqtab", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--asv-to-cigar", required=True)
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--identity", type=float, default=0.989)
    parser.add_argument("--word-size", type=int, default=10)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-total-reads", type=int, default=100)
    parser.add_argument("--min-samples", type=int, default=2)
    parser.add_argument("--exclude-bimeras", action="store_true")
    parser.add_argument("--executable", default="cd-hit-est")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not 0.8 <= args.identity <= 1.0:
        raise SystemExit("--identity must be between 0.8 and 1.0")
    if not 4 <= args.word_size <= 12:
        raise SystemExit("--word-size must be between 4 and 12")
    summary = run_clustering(args)
    print(
        f"CD-HIT-EST retained {summary['input_asvs']} ASVs in "
        f"{summary['clusters']} clusters at {float(summary['identity_threshold']):.1%} identity."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
