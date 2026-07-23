#!/usr/bin/env python3
"""CLI helpers used by the first-pass SIMPLseq Nextflow port."""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
from collections import Counter
from pathlib import Path

import pandas as pd


MONTH_PATTERN = (
    r"Jan(?:uary)?|Feb(?:ruary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|"
    r"Jul(?:y)?|Aug(?:ust)?|Sep(?:t(?:ember)?)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?"
)
TECHNICAL_SAMPLE_TOKEN_RE = re.compile(
    r"^(?:mpg|lib|library|lane|pool|amplicon|run|plate|barcode|bc|index|l|s|r)[0-9A-Za-z-]*$",
    re.IGNORECASE,
)


def is_month_token(token: str) -> bool:
    return bool(re.fullmatch(rf"{MONTH_PATTERN}", token or "", re.IGNORECASE))


def is_replicate_token(token: str) -> bool:
    return bool(re.fullmatch(r"Rep(?:licate)?[0-9A-Za-z]*", token or "", re.IGNORECASE))


def is_date_like_token(token: str) -> bool:
    token = token or ""
    return bool(
        re.fullmatch(r"(?:19|20)[0-9]{2}", token)
        or re.fullmatch(r"[0-9]{1,2}", token)
        or re.fullmatch(r"[0-9]{6,8}", token)
        or re.fullmatch(rf"(?:{MONTH_PATTERN})[0-9]{{4}}", token, re.IGNORECASE)
        or re.fullmatch(rf"[0-9]{{4}}(?:{MONTH_PATTERN})", token, re.IGNORECASE)
    )


def strip_sample_label_noise(label: str) -> str:
    """Remove read, pool, replicate, and sequencer suffixes that are not participant IDs."""
    cleaned = os.path.basename(str(label or "")).strip()
    cleaned = re.sub(r"\.(?:f(?:ast)?q)(?:\.gz)?$", "", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"[-_ .]*R[12](?:[-_ .]*001)?$", "", cleaned, flags=re.IGNORECASE)
    for _ in range(4):
        updated = re.sub(r"[-_ .]*Pool[0-9A-Za-z]+$", "", cleaned, flags=re.IGNORECASE)
        updated = re.sub(r"[-_ .]*S[0-9]+[A-Za-z]*$", "", updated, flags=re.IGNORECASE)
        updated = re.sub(r"[-_ .]*Rep(?:licate)?[-_ .]*[0-9A-Za-z]+$", "", updated, flags=re.IGNORECASE)
        updated = updated.strip("-_ .")
        if updated == cleaned:
            break
        cleaned = updated
    return cleaned


def participant_token_from_label(label: str) -> str:
    """Return the strongest participant-like token embedded in a mixed sample label."""
    tokens = [token for token in re.split(r"[^A-Za-z0-9]+", label or "") if token]
    for token in tokens:
        if is_month_token(token) or is_replicate_token(token) or is_date_like_token(token):
            continue
        if TECHNICAL_SAMPLE_TOKEN_RE.fullmatch(token):
            continue
        if re.search(r"[A-Za-z]", token) and re.search(r"[0-9]", token):
            return token
    return ""


def clean_patient_candidate(candidate: str) -> str:
    cleaned = strip_sample_label_noise(candidate)
    cleaned = re.sub(
        rf"[-_ .]*(?:{MONTH_PATTERN})[-_ .]*(?:19|20)?[0-9]{{0,2}}$",
        "",
        cleaned,
        flags=re.IGNORECASE,
    ).strip("-_ .")

    token_candidate = participant_token_from_label(cleaned)
    if token_candidate:
        return token_candidate

    pieces = [piece for piece in re.split(r"[^A-Za-z0-9]+", cleaned) if piece]
    while pieces and (TECHNICAL_SAMPLE_TOKEN_RE.fullmatch(pieces[0]) or is_date_like_token(pieces[0])):
        pieces.pop(0)
    while pieces and (
        TECHNICAL_SAMPLE_TOKEN_RE.fullmatch(pieces[-1])
        or is_month_token(pieces[-1])
        or is_replicate_token(pieces[-1])
        or is_date_like_token(pieces[-1])
    ):
        pieces.pop()
    if pieces:
        return "-".join(pieces)
    return cleaned


def read_fasta_lengths(path: str) -> list[tuple[str, int]]:
    records = []
    name = None
    seq = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, len("".join(seq))))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        records.append((name, len("".join(seq))))
    return records


def resolve_manifest_path(value: str, manifest_path: str, manifest_root: str | None = None) -> str:
    if not value:
        return value
    path = Path(value)
    if path.is_absolute():
        return str(path)
    base = Path(manifest_root).resolve() if manifest_root else Path(manifest_path).resolve().parent
    return str((base / path).resolve())


def cmd_preflight(args: argparse.Namespace) -> int:
    read_length = args.read_length
    required = ["sample_id", "fastq_1", "fastq_2"]
    rows = []

    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    with open(args.samples, newline="") as handle:
        reader = csv.DictReader(handle)
        missing_cols = [col for col in required if col not in (reader.fieldnames or [])]
        if missing_cols:
            rows.append(("ERROR", "samples_csv", "Missing required columns: " + ",".join(missing_cols)))
            sample_rows = []
        else:
            sample_rows = list(reader)

    if not sample_rows:
        rows.append(("ERROR", "samples_csv", "No samples found. Use `simplseq scan --fastq-dir data --out samples.csv`, or create samples.csv manually."))

    ids = [row.get("sample_id", "") for row in sample_rows]
    for sample_id, count in Counter(ids).items():
        if sample_id and count > 1:
            rows.append(("ERROR", sample_id, "Duplicate sample_id"))

    for row in sample_rows:
        sample_id = row.get("sample_id", "")
        for col in ["fastq_1", "fastq_2"]:
            value = row.get(col, "")
            resolved = resolve_manifest_path(value, args.samples, args.samples_root)
            if not value:
                rows.append(("ERROR", sample_id, f"Missing {col}"))
            elif not os.path.exists(resolved):
                rows.append(("ERROR", sample_id, f"{col} not found: {value}"))
            elif not os.path.exists(resolved + ".md5"):
                rows.append(("INFO", sample_id, f"{col} has no optional .md5 sidecar; SIMPLseq will compute input MD5: {value}.md5"))
        if row.get("sample_type", "").lower() in {"negative", "ntc", "control_negative"}:
            rows.append(("INFO", sample_id, "Negative control marked in manifest"))

    barcode_rows = [
        row for row in sample_rows
        if row.get("expected_fwd_barcode", "") or row.get("expected_rev_barcode", "")
    ]
    if args.inline_barcodes_enabled.lower() == "false":
        pass
    elif barcode_rows:
        incomplete = [
            row.get("sample_id", "") for row in sample_rows
            if bool(row.get("expected_fwd_barcode", "")) != bool(row.get("expected_rev_barcode", ""))
        ]
        if incomplete:
            rows.append(("ERROR", "inline_barcodes", "Incomplete barcode pairs: " + ",".join(incomplete[:20])))
        else:
            rows.append(("INFO", "inline_barcodes", f"Barcode columns populated for {len(barcode_rows)} samples"))
    else:
        rows.append(("WARN", "inline_barcodes", "No expected barcode pairs found; contamination module cannot run yet"))

    with open(args.geometry, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["amplicon", "length_bp", "pe150_raw_overlap_bp", "recommended_handling"])
        for name, length in read_fasta_lengths(args.amplicons):
            overlap = (2 * read_length) - length
            if overlap >= 20:
                handling = "merge"
            elif overlap >= 12:
                handling = "borderline_validate"
            else:
                handling = "concatenate"
            writer.writerow([name, length, overlap, handling])

    with open(args.barcode, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["status", "metric", "value"])
        writer.writerow(["INFO", "samples_with_expected_barcode_pairs", len(barcode_rows)])
        writer.writerow(["INFO", "sentinel_locus", args.sentinel_locus])

    if not any(status == "ERROR" for status, _, _ in rows):
        rows.insert(0, ("OK", "preflight", f"{len(sample_rows)} samples validated"))

    with open(args.report, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["status", "scope", "message"])
        writer.writerows(rows)

    if any(status == "ERROR" for status, _, _ in rows):
        raise SystemExit("Preflight failed. See " + args.report)
    return 0


def cmd_write_meta(args: argparse.Namespace) -> int:
    with open(args.samples, newline="") as handle, open(args.out, "w", newline="") as out:
        reader = csv.DictReader(handle)
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for row in reader:
            writer.writerow([
                row["sample_id"],
                resolve_manifest_path(row["fastq_1"], args.samples, args.samples_root),
                resolve_manifest_path(row["fastq_2"], args.samples, args.samples_root),
            ])
    return 0


def cmd_write_pipeline_json(args: argparse.Namespace) -> int:
    out = {
        "path_to_meta": args.meta,
        "Class": args.pipeline_class,
        "maxEE": args.max_ee,
        "trimRight": args.trim_right,
        "minLen": args.min_len,
        "truncQ": args.trunc_q,
        "max_consist": args.max_consist,
        "omegaA": args.omega_a,
        "justConcatenate": args.just_concatenate,
        "saveRdata": args.save_rdata,
        "dada2_randomize": args.dada2_randomize,
        "dada2_multithread": args.dada2_multithread,
        "dada2_seed": args.dada2_seed,
        "pr1": args.primers_fwd,
        "pr2": args.primers_rev,
        "overlap_pr1": args.overlap_primers_fwd,
        "overlap_pr2": args.overlap_primers_rev,
    }
    with open(args.out, "w") as handle:
        json.dump(out, handle, indent=2)
        handle.write("\n")
    return 0


def read_seqtab(path: str) -> pd.DataFrame:
    with open(path, newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    duplicate_columns = sorted(name for name, count in Counter(header[1:]).items() if count > 1)
    if duplicate_columns:
        preview = ", ".join(duplicate_columns[:5])
        raise ValueError(f"ASV count table has duplicate sequence columns: {preview}")

    df = pd.read_csv(path, sep="\t")
    if df.empty:
        raise SystemExit(f"Empty seqtab: {path}")
    df = df.rename(columns={df.columns[0]: "sample"}).set_index("sample")
    return df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)


def read_bimera(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        return pd.DataFrame(columns=["sequence", "bimera"])
    frame = pd.read_csv(path, sep="\t")
    missing = sorted({"sequence", "bimera"} - set(frame.columns))
    if missing:
        raise ValueError(f"Chimera table {path} is missing columns: {', '.join(missing)}")
    return frame.astype({"sequence": str})


def parse_bimera_status(value: object, *, source: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "1", "yes", "y"}:
        return True
    if normalized in {"false", "f", "0", "no", "n"}:
        return False
    raise ValueError(f"Invalid chimera status {value!r} for {source}; expected TRUE or FALSE")


def validate_cdhit_summary(summary: dict[str, object]) -> dict[str, int]:
    """Validate CD-HIT output invariants owned by report preparation."""
    if not summary:
        return {"input_asvs": 0, "clusters": 0, "merged_asvs": 0}

    try:
        input_asvs = int(summary.get("input_asvs", 0) or 0)
        clusters = int(summary.get("clusters", 0) or 0)
        merged_asvs = int(summary.get("merged_asvs", input_asvs - clusters) or 0)
    except (TypeError, ValueError) as exc:
        raise ValueError("CD-HIT summary contains non-integer ASV counts") from exc

    if min(input_asvs, clusters, merged_asvs) < 0:
        raise ValueError("CD-HIT summary contains negative ASV counts")
    if input_asvs != clusters + merged_asvs:
        raise ValueError(
            "CD-HIT ASV-count conservation failed: "
            f"{input_asvs} input ASVs != {clusters} output clusters + "
            f"{merged_asvs} ASVs consolidated"
        )

    unresolved = summary.get("unresolved_asv_columns", [])
    if isinstance(unresolved, str):
        unresolved = [item.strip() for item in unresolved.split(",") if item.strip()]
    unresolved_count = int(summary.get("unresolved_asv_column_count", len(unresolved or [])) or 0)
    if unresolved_count:
        preview = ", ".join(str(item) for item in list(unresolved or [])[:5])
        detail = f": {preview}" if preview else ""
        raise ValueError(
            f"CD-HIT output has {unresolved_count} unresolved ASV-column mapping(s){detail}"
        )

    if "input_cluster_reads" in summary and "total_cluster_reads" in summary:
        input_reads = int(summary.get("input_cluster_reads", 0) or 0)
        output_reads = int(summary.get("total_cluster_reads", 0) or 0)
        if input_reads != output_reads:
            raise ValueError(
                "CD-HIT read-count conservation failed: "
                f"{input_reads} input reads != {output_reads} clustered output reads"
            )

    return {
        "input_asvs": input_asvs,
        "clusters": clusters,
        "merged_asvs": merged_asvs,
    }


def cmd_prepare_stage2(args: argparse.Namespace) -> int:
    seqtab = read_seqtab(args.seqtab)
    valid_cols = [
        col for col in seqtab.columns
        if set(str(col)).issubset(set("ACGTN")) and len(str(col)) >= args.strict_min_asv_length
    ]
    strict = seqtab.loc[:, valid_cols]
    input_retained_reads = int(seqtab.loc[:, valid_cols].to_numpy().sum())
    strict.to_csv(args.strict_seqtab, sep="\t", index_label="sample")
    written_strict = read_seqtab(args.strict_seqtab)
    output_retained_reads = int(written_strict.to_numpy().sum())
    if (
        input_retained_reads != output_retained_reads
        or list(strict.columns) != list(written_strict.columns)
        or list(strict.index.astype(str)) != list(written_strict.index.astype(str))
        or strict.to_numpy().tolist() != written_strict.to_numpy().tolist()
    ):
        raise ValueError(
            "Stage-2 ASV count conservation failed after writing the strict table: "
            f"{input_retained_reads} retained input reads != {output_retained_reads} output reads"
        )

    op_bimera = read_bimera(args.op_bimera)
    nop_bimera = read_bimera(args.nop_bimera)
    corrected = pd.read_csv(args.corrected_asv, sep="\t", dtype=str, keep_default_na=False)
    corrected_missing = sorted({"ASV", "correctedASV"} - set(corrected.columns))
    if corrected_missing:
        raise ValueError(
            f"Corrected-ASV table {args.corrected_asv} is missing columns: "
            + ", ".join(corrected_missing)
        )

    mapped_rows = []
    if not corrected.empty and not nop_bimera.empty:
        bmap: dict[str, list[bool]] = {}
        for _, row in nop_bimera.iterrows():
            sequence = str(row["sequence"])
            bmap.setdefault(sequence, []).append(
                parse_bimera_status(row["bimera"], source=f"non-overlap ASV {sequence}")
            )
        for _, row in corrected.iterrows():
            corrected_asv = str(row.get("correctedASV", "")).strip()
            original_asv = str(row.get("ASV", ""))
            if corrected_asv.upper() in {"", "NA", "NAN"} or corrected_asv not in strict.columns:
                continue
            if original_asv not in bmap:
                raise ValueError(
                    "Unable to resolve corrected ASV-column mapping: "
                    f"source ASV {original_asv!r} (corrected to {corrected_asv!r}) "
                    "is absent from the non-overlap chimera table"
                )
            mapped_rows.append({"sequence": corrected_asv, "bimera": any(bmap[original_asv])})

    op_rows = []
    for _, row in op_bimera.iterrows():
        sequence = str(row["sequence"])
        if sequence in strict.columns:
            op_rows.append({
                "sequence": sequence,
                "bimera": parse_bimera_status(row["bimera"], source=f"overlap ASV {sequence}"),
            })
    bimera = pd.DataFrame(mapped_rows + op_rows, columns=["sequence", "bimera"])
    if not bimera.empty:
        bimera = bimera.groupby("sequence", sort=False, as_index=False)["bimera"].any()
        bimera = bimera[bimera["sequence"].isin(strict.columns)]

    status_by_sequence = {
        str(row["sequence"]): bool(row["bimera"])
        for _, row in bimera.iterrows()
    }
    unresolved_columns = [str(column) for column in strict.columns if str(column) not in status_by_sequence]
    if unresolved_columns:
        preview = ", ".join(unresolved_columns[:5])
        raise ValueError(
            f"Unable to resolve chimera status for {len(unresolved_columns)} retained ASV column(s): {preview}"
        )
    if len(status_by_sequence) != len(strict.columns):
        raise ValueError(
            "Stage-2 ASV feature-count conservation failed: "
            f"{len(strict.columns)} retained ASV columns != {len(status_by_sequence)} chimera decisions"
        )

    with open(args.strict_bimera, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sequence", "bimera"])
        for sequence in strict.columns:
            writer.writerow([sequence, "TRUE" if status_by_sequence[str(sequence)] else "FALSE"])
    return 0


def cmd_prepare_stage3(args: argparse.Namespace) -> int:
    mapped = pd.read_csv(args.mapped, sep="\t", dtype=str, keep_default_na=False)
    if "hapid" not in mapped.columns:
        raise ValueError(f"Mapped ASV table has no hapid column: {args.mapped}")
    mapped_ids = mapped["hapid"].astype(str)
    duplicate_ids = sorted(mapped_ids[mapped_ids.duplicated(keep=False)].unique())
    if duplicate_ids:
        raise ValueError(f"Mapped ASV table has duplicate hapid values: {', '.join(duplicate_ids[:5])}")

    seqtab = pd.read_csv(args.seqtab, sep="\t")
    seqtab = seqtab.rename(columns={seqtab.columns[0]: "sample"})
    count_cols = [col for col in seqtab.columns if col != "sample"]
    expected_ids = {f"ASV{index}" for index in range(1, len(count_cols) + 1)}
    unresolved_ids = sorted(expected_ids - set(mapped_ids))
    unexpected_ids = sorted(set(mapped_ids) - expected_ids)
    if unresolved_ids or unexpected_ids:
        details = []
        if unresolved_ids:
            details.append("missing " + ", ".join(unresolved_ids[:5]))
        if unexpected_ids:
            details.append("unexpected " + ", ".join(unexpected_ids[:5]))
        raise ValueError(
            "Unable to resolve ASV-column mappings for CD-HIT preparation: " + "; ".join(details)
        )

    if "mapping_status" in mapped.columns:
        filtered = mapped[mapped["mapping_status"].astype(str).str.upper().eq("MAPPED")]
    elif "refid_3D7" in mapped.columns:
        filtered = mapped[mapped["refid_3D7"].astype(str) != "NA"]
    else:
        filtered = mapped
    filtered.to_csv(args.filtered_mapped, sep="\t", index=False)

    seqtab[count_cols] = seqtab[count_cols].fillna(0).astype(int)
    seqtab.to_csv(args.fixed_seqtab, sep="\t", index=False)
    written = pd.read_csv(args.fixed_seqtab, sep="\t")
    written_count_cols = [col for col in written.columns if col != "sample"]
    if (
        count_cols != written_count_cols
        or seqtab["sample"].astype(str).tolist() != written["sample"].astype(str).tolist()
        or seqtab[count_cols].to_numpy().tolist() != written[written_count_cols].to_numpy().tolist()
    ):
        input_reads = int(seqtab[count_cols].to_numpy().sum())
        output_reads = int(written[written_count_cols].to_numpy().sum())
        raise ValueError(
            "CD-HIT preparation count conservation failed: "
            f"{input_reads} input reads != {output_reads} written output reads"
        )
    return 0


def row_value(row: list[str], header: list[str], name: str, fallback: int) -> str:
    if name in header:
        return row[header.index(name)]
    return row[fallback]


def passes_filters(row: list[str], header: list[str], min_reads: int, min_samples: int, include_failed: bool, exclude_bimeras: bool) -> bool:
    total_reads = int(row_value(row, header, "total_reads", 2))
    total_samples = int(row_value(row, header, "total_samples", 3))
    if total_reads < min_reads or total_samples < min_samples:
        return False
    if not include_failed:
        snv_filter = row_value(row, header, "snv_filter", -3)
        indel_filter = row_value(row, header, "indel_filter", -2)
        if snv_filter == "FAIL" or indel_filter == "FAIL":
            return False
    if exclude_bimeras:
        bimera = row_value(row, header, "bimera", -1).upper()
        if bimera == "TRUE":
            return False
    return True


def cmd_check_cigar_inputs(args: argparse.Namespace) -> int:
    counts: dict[str, dict[str, int]] = {}
    with open(args.table, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        for row in reader:
            if not row:
                continue
            amplicon = row_value(row, header, "refid_3D7", 5)
            if amplicon == "NA":
                continue
            if amplicon not in counts:
                counts[amplicon] = {"input_asvs": 0, "passing_asvs": 0}
            counts[amplicon]["input_asvs"] += 1
            if passes_filters(row, header, args.min_reads, args.min_samples, args.include_failed, args.exclude_bimeras):
                counts[amplicon]["passing_asvs"] += 1

    with open(args.summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["amplicon", "input_asvs", "passing_asvs"])
        for amplicon in sorted(counts):
            writer.writerow([amplicon, counts[amplicon]["input_asvs"], counts[amplicon]["passing_asvs"]])

    passing_total = sum(item["passing_asvs"] for item in counts.values())
    if passing_total == 0:
        raise SystemExit("No ASVs pass CIGAR filters.")
    return 0


def read_tsv(path: str, limit: int = 12) -> tuple[list[str], list[list[str]]]:
    if not os.path.exists(path):
        return [], []
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        return [], []
    return rows[0], rows[1:limit + 1]


def format_number(value: object) -> str:
    try:
        if isinstance(value, float) and not value.is_integer():
            return f"{value:,.1f}"
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value)


def html_table(path: str, title: str, limit: int = 12, collapsed: bool = False) -> str:
    header, rows = read_tsv(path, limit=limit)
    if not header:
        return f"<h2>{html.escape(title)}</h2><p>Not available.</p>"
    head = "".join(f"<th>{html.escape(str(col))}</th>" for col in header)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(cell))}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    table = f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    if collapsed:
        return f"<details><summary>{html.escape(title)}</summary>{table}</details>"
    return f"<h2>{html.escape(title)}</h2>{table}"


def variant_amplicon(name: str) -> str:
    text = str(name)
    for sep in [":", ","]:
        if sep in text:
            return text.split(sep, 1)[0]
    return text


def metric_card(label: str, value: object, note: str = "") -> str:
    note_html = f"<span>{html.escape(note)}</span>" if note else ""
    return f"<div class=\"card\"><b>{html.escape(format_number(value))}</b><em>{html.escape(label)}</em>{note_html}</div>"


def metric_value(value: object) -> str:
    return html.escape(format_number(value))


def bar_rows(rows: list[dict[str, object]], value_key: str, label_key: str, columns: list[tuple[str, str]]) -> str:
    if not rows:
        return "<p class=\"note\">Not available.</p>"
    max_value = max(float(row.get(value_key, 0) or 0) for row in rows) or 1.0
    header = "".join(f"<th>{html.escape(title)}</th>" for title, _ in columns)
    body = []
    for row in rows:
        width = 100.0 * float(row.get(value_key, 0) or 0) / max_value
        cells = []
        for title, key in columns:
            value = row.get(key, "")
            if key == label_key:
                cells.append(
                    f"<td><div class=\"bar-label\">{html.escape(str(value))}</div>"
                    f"<div class=\"bar\"><span style=\"width:{width:.1f}%\"></span></div></td>"
                )
            else:
                cells.append(f"<td>{metric_value(value)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def horizontal_chart(rows: list[dict[str, object]], label_key: str, value_key: str, title: str, max_items: int = 12) -> str:
    if not rows:
        return ""
    data = rows[:max_items]
    max_value = max(float(row.get(value_key, 0) or 0) for row in data) or 1.0
    chart_rows = []
    for row in data:
        label = html.escape(str(row.get(label_key, "")))
        value = float(row.get(value_key, 0) or 0)
        width = 100.0 * value / max_value
        chart_rows.append(
            f'<div class="chart-row"><span class="chart-label">{label}</span>'
            f'<span class="chart-track"><span style="width:{width:.1f}%"></span></span>'
            f'<span class="chart-value">{html.escape(format_number(value))}</span></div>'
        )
    return f'<div class="chart"><h3>{html.escape(title)}</h3>{"".join(chart_rows)}</div>'


def locus_recovery_chart(sample_rows: list[dict[str, object]], total_loci: int = 6) -> str:
    if not sample_rows:
        return ""
    total_loci = max(0, int(total_loci or 0))
    counts = Counter(int(row.get("loci_detected", 0) or 0) for row in sample_rows)
    rows = [{"loci": f"{i}/{total_loci}", "samples": counts.get(i, 0)} for i in range(total_loci, -1, -1)]
    return horizontal_chart(rows, "loci", "samples", "Sample locus recovery", max_items=total_loci + 1)


def infer_patient_id_from_sample(sample: object) -> str:
    """Infer a patient/participant ID from a longitudinal sample label."""
    label = str(sample or "").strip()
    if not label:
        return ""
    label = strip_sample_label_noise(label)

    token_candidate = participant_token_from_label(label)
    if token_candidate:
        return token_candidate

    compact = re.match(
        rf"^(?P<participant>.+?)(?P<month>{MONTH_PATTERN})(?P<year>[0-9]{{4}})"
        r"(?:[-_ .]*Rep(?:licate)?[-_ .]*[0-9A-Za-z]+)?$",
        label,
        re.IGNORECASE,
    )
    if compact:
        return clean_patient_candidate(compact.group("participant")) or label

    date_match = re.search(
        rf"(?P<date>(?:{MONTH_PATTERN})[-_ .]*[0-9]{{4}}|"
        rf"[0-9]{{4}}[-_ .]*(?:{MONTH_PATTERN}|[0-9]{{1,2}})(?:[-_ .]*[0-9]{{1,2}})?|"
        r"[0-9]{8})",
        label,
        re.IGNORECASE,
    )
    if date_match:
        candidates = [label[: date_match.start()], label[date_match.end() :]]
        for candidate in candidates:
            cleaned = clean_patient_candidate(candidate)
            if cleaned:
                return cleaned

    return clean_patient_candidate(label) or label


def summarize_cigar_matrix(
    path: str,
    samples_path: str = "",
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    if not os.path.exists(path):
        return {}, [], []
    df = pd.read_csv(path, sep="\t")
    if df.empty or len(df.columns) < 2:
        return {}, [], []
    sample_col = df.columns[0]
    count_cols = [col for col in df.columns if col != sample_col]
    counts = df[count_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    amp_by_col = {col: variant_amplicon(col) for col in count_cols}

    sample_rows = []
    for idx, row in counts.iterrows():
        positive = [col for col in count_cols if int(row[col]) > 0]
        sample_rows.append({
            "sample": df.loc[idx, sample_col],
            "total_reads": int(row.sum()),
            "loci_detected": len({amp_by_col[col] for col in positive}),
            "haplotypes_detected": len(positive),
        })
    sample_rows.sort(key=lambda item: int(item["total_reads"]), reverse=True)

    amp_rows = []
    for amplicon in sorted(set(amp_by_col.values())):
        cols = [col for col in count_cols if amp_by_col[col] == amplicon]
        sub = counts[cols]
        per_sample = sub.sum(axis=1)
        amp_rows.append({
            "amplicon": amplicon,
            "final_reads": int(per_sample.sum()),
            "samples_detected": int((per_sample >= 100).sum()),
            "haplotypes": int((sub.sum(axis=0) > 0).sum()),
            "median_reads_detected_samples": int(per_sample[per_sample > 0].median()) if (per_sample > 0).any() else 0,
        })
    amp_rows.sort(key=lambda item: int(item["final_reads"]), reverse=True)

    if samples_path and os.path.exists(samples_path):
        sample_sheet = pd.read_csv(samples_path, dtype=str).fillna("")
        participant_by_sample = {
            str(row["sample_id"]).strip(): str(row["participant_id"]).strip()
            for _, row in sample_sheet.iterrows()
            if str(row.get("sample_id", "")).strip()
        }
        patients = len({
            participant_by_sample.get(str(value).strip(), "")
            for value in df[sample_col]
            if participant_by_sample.get(str(value).strip(), "")
        })
    else:
        patients = len({
            patient
            for patient in (infer_patient_id_from_sample(value) for value in df[sample_col])
            if patient
        })

    overview = {
        "samples": int(len(df)),
        "patients": patients,
        "total_final_reads": int(counts.to_numpy().sum()),
        "samples_with_reads": sum(1 for row in sample_rows if int(row["total_reads"]) > 0),
        "median_loci_detected": float(pd.Series([row["loci_detected"] for row in sample_rows]).median()) if sample_rows else 0,
        "total_haplotypes": len(count_cols),
    }
    return overview, amp_rows, sample_rows


def summarize_mapped_asvs(path: str) -> tuple[dict[str, object], list[dict[str, object]]]:
    if not os.path.exists(path):
        return {}, []
    df = pd.read_csv(path, sep="\t")
    if df.empty or "refid_3D7" not in df.columns:
        return {}, []
    df = df[df["refid_3D7"].astype(str) != "NA"].copy()
    if df.empty:
        return {}, []
    df["total_reads"] = pd.to_numeric(df.get("total_reads", 0), errors="coerce").fillna(0).astype(int)
    df["total_samples"] = pd.to_numeric(df.get("total_samples", 0), errors="coerce").fillna(0).astype(int)
    rows = []
    for amplicon, group in df.groupby("refid_3D7"):
        pass_mask = pd.Series([True] * len(group), index=group.index)
        if "snv_filter" in group:
            pass_mask &= group["snv_filter"].astype(str).str.upper().eq("PASS")
        if "indel_filter" in group:
            pass_mask &= group["indel_filter"].astype(str).str.upper().eq("PASS")
        if "bimera" in group:
            pass_mask &= ~group["bimera"].astype(str).str.upper().eq("TRUE")
        rows.append({
            "amplicon": amplicon,
            "mapped_asvs": int(len(group)),
            "filter_pass_asvs": int(pass_mask.sum()),
            "mapped_reads": int(group["total_reads"].sum()),
            "max_samples_per_asv": int(group["total_samples"].max()),
        })
    rows.sort(key=lambda item: int(item["mapped_reads"]), reverse=True)
    overview = {
        "mapped_asvs": int(len(df)),
        "mapped_reads": int(df["total_reads"].sum()),
        "amplicons_with_asvs": int(df["refid_3D7"].nunique()),
    }
    return overview, rows


def summarize_asv_filter_steps(
    path: str,
    *,
    min_total_reads: int = 100,
    min_samples: int = 2,
    exclude_bimeras: bool = True,
) -> list[dict[str, object]]:
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return []

    ref_cols = [col for col in df.columns if str(col).startswith("refid_")]
    mapped_mask = pd.Series([True] * len(df), index=df.index)
    if "mapping_status" in df.columns:
        mapped_mask = df["mapping_status"].astype(str).str.upper().eq("MAPPED")
    elif ref_cols:
        mapped_mask = pd.Series([False] * len(df), index=df.index)
        for col in ref_cols:
            values = df[col].astype(str).str.strip().str.upper()
            mapped_mask |= ~values.isin({"", "NA", "NAN"})

    total_reads = pd.to_numeric(df.get("total_reads", 0), errors="coerce").fillna(0)
    total_samples = pd.to_numeric(df.get("total_samples", 0), errors="coerce").fillna(0)
    steps: list[dict[str, object]] = [
        {
            "step": "DADA2/post-processing ASVs",
            "requirement": "All ASV sequence features entering reference/locus mapping.",
            "count": int(len(df)),
        },
        {
            "step": "Mapped ASVs",
            "requirement": "ASV has an assigned amplicon/locus reference.",
            "count": int(mapped_mask.sum()),
        },
    ]
    current = mapped_mask.copy()

    feature_filters = []
    if min_total_reads > 0:
        feature_filters.append((
            "Cohort-wide read filter",
            f"Keep ASVs with total_reads >= {min_total_reads:,} across the run.",
            total_reads >= min_total_reads,
        ))
    if min_samples > 0:
        feature_filters.append((
            "Blanket recurrence filter",
            f"Keep ASVs observed in at least {min_samples:,} distinct biological samples.",
            total_samples >= min_samples,
        ))
    for label, requirement, mask in feature_filters:
        current &= mask
        steps.append({"step": label, "requirement": requirement, "count": int(current.sum())})

    if "snv_filter" in df.columns:
        current &= df["snv_filter"].astype(str).str.upper().eq("PASS")
        steps.append({
            "step": "SNV-distance filter",
            "requirement": "Keep ASVs marked PASS by the configured SNV/edit-distance rules.",
            "count": int(current.sum()),
        })

    if "indel_filter" in df.columns:
        current &= df["indel_filter"].astype(str).str.upper().eq("PASS")
        steps.append({
            "step": "INDEL-distance filter",
            "requirement": "Keep ASVs within the configured locus-specific inserted/deleted-base limit.",
            "count": int(current.sum()),
        })

    if exclude_bimeras and "bimera" in df.columns:
        current &= ~df["bimera"].astype(str).str.upper().eq("TRUE")
        steps.append({
            "step": "Chimera filter",
            "requirement": "Remove ASVs marked as bimeras/chimeras.",
            "count": int(current.sum()),
        })

    return steps


def summarize_asv_to_cigar(path: str) -> list[dict[str, object]]:
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path, sep="\t")
    if df.empty or "Amplicon" not in df.columns:
        return []
    rows = []
    for amplicon, group in df.groupby("Amplicon"):
        rows.append({
            "amplicon": amplicon,
            "asvs_with_cigar": int(len(group)),
            "unique_cigars": int(group["CIGAR"].nunique()) if "CIGAR" in group else "",
        })
    rows.sort(key=lambda item: int(item["asvs_with_cigar"]), reverse=True)
    return rows


def write_asv_filtering_summary(
    path: str,
    *,
    project_name: str,
    mapped_path: str,
    asv_to_cigar_path: str,
    cdhit_summary: dict[str, object],
    min_total_reads: int = 100,
    min_samples: int = 2,
    exclude_bimeras: bool = True,
) -> None:
    validated_cdhit = validate_cdhit_summary(cdhit_summary)
    steps = summarize_asv_filter_steps(
        mapped_path,
        min_total_reads=min_total_reads,
        min_samples=min_samples,
        exclude_bimeras=exclude_bimeras,
    )
    asv_cigar_rows = summarize_asv_to_cigar(asv_to_cigar_path)
    cigar_asvs = sum(int(row.get("asvs_with_cigar", 0) or 0) for row in asv_cigar_rows)
    unique_cigars = sum(int(row.get("unique_cigars", 0) or 0) for row in asv_cigar_rows)

    lines = [
        f"{project_name} ASV filtering summary",
        "=" * (len(project_name) + 22),
        "",
        "Purpose",
        "-------",
        "This file tracks how many ASV sequence features remain after each feature-level QC step.",
        "CD-HIT is reported separately because it does not alter the primary exact-CIGAR analysis.",
        "",
        "Pipeline ASV feature filters",
        "----------------------------",
    ]

    previous: int | None = None
    baseline: int | None = None
    for item in steps:
        count = int(item["count"])
        if previous is None:
            change = "baseline"
            baseline = count
        else:
            removed = max(0, previous - count)
            percent = (100.0 * count / previous) if previous else 0.0
            baseline_percent = (100.0 * count / baseline) if baseline else 0.0
            change = (
                f"removed {removed:,}; retained {percent:.1f}% of previous step; "
                f"{baseline_percent:.1f}% of starting ASVs remain"
            )
        lines.extend([
            f"- {item['step']}: {count:,} ASVs ({change})",
            f"  Requirement: {item['requirement']}",
        ])
        previous = count

    if steps:
        filter_pass = int(steps[-1]["count"])
        lines.extend([
            "",
            "CIGAR conversion",
            "----------------",
            f"- Filter-pass ASVs submitted to CIGAR conversion: {filter_pass:,}",
            f"- ASVs represented in ASV-to-CIGAR map: {cigar_asvs:,}",
            f"- Unique CIGAR haplotypes after conversion: {unique_cigars:,}",
        ])

    if cdhit_summary:
        input_asvs = validated_cdhit["input_asvs"]
        clusters = validated_cdhit["clusters"]
        merged = validated_cdhit["merged_asvs"]
        identity = float(cdhit_summary.get("identity_threshold", 0.989) or 0.989)
        identity_label = f"{identity:.1%}"
        lines.extend([
            "",
            f"CD-HIT {identity_label} sensitivity clustering",
            "------------------------------------",
            f"- Input to CD-HIT: {input_asvs:,} filter-pass nucleotide ASVs.",
            f"- CD-HIT-EST threshold: {identity_label} global sequence identity.",
            f"- Output clusters: {clusters:,}.",
            f"- Sequence features grouped away by clustering: {merged:,}.",
            "",
            "Interpretation: CD-HIT groups similar, QC-passed ASV sequences into locus-specific clusters.",
            "It does not overwrite the primary exact-CIGAR table. Optional sensitivity analyses use the summed cluster counts separately.",
        ])

    Path(path).write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def read_tsv_dicts(path: str) -> list[dict[str, str]]:
    if not os.path.exists(path):
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compact_number(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if abs(number) >= 1_000_000:
        return f"{number / 1_000_000:.1f}M"
    if abs(number) >= 1_000:
        return f"{number:,.0f}"
    if number.is_integer():
        return f"{int(number):,}"
    return f"{number:,.1f}"


def compact_table_cell(value: object) -> str:
    text = str(value)
    numeric_text = text.replace(",", "")
    if numeric_text.replace(".", "", 1).lstrip("-").isdigit():
        return compact_number(numeric_text)
    return text


def is_optional_md5_notice(row: dict[str, str]) -> bool:
    return ".md5 sidecar" in row.get("message", "").lower()


def report_metric_card(label: str, value: object, note: str = "", tone: str = "") -> str:
    tone_class = f" {html.escape(tone)}" if tone else ""
    note_html = f"<span>{html.escape(note)}</span>" if note else ""
    return (
        f'<article class="metric{tone_class}">'
        "<i></i>"
        f"<span>{html.escape(label)}</span>"
        f"<b>{html.escape(str(value))}</b>"
        f"{note_html}"
        "</article>"
    )


def svg_bar_chart(
    rows: list[dict[str, object]],
    label_key: str,
    value_key: str,
    title: str,
    subtitle: str = "",
    color: str = "#20988c",
    max_items: int = 12,
    wide: bool = False,
) -> str:
    data = rows[:max_items]
    if not data:
        return '<div class="chart"><h3>No data available</h3></div>'
    max_value = max(float(row.get(value_key, 0) or 0) for row in data) or 1.0
    chart_rows = []
    for row in data:
        label = str(row.get(label_key, ""))
        value = float(row.get(value_key, 0) or 0)
        width_pct = 100 * value / max_value
        chart_rows.append(
            '<div class="chart-row">'
            f'<span class="chart-label" title="{html.escape(label)}">{html.escape(label)}</span>'
            '<span class="chart-track">'
            f'<span style="width:{width_pct:.1f}%"></span>'
            "</span>"
            f'<span class="chart-value">{html.escape(compact_number(value))}</span>'
            "</div>"
        )
    classes = "chart chart-wide" if wide else "chart"
    subtitle_html = f'<p class="note">{html.escape(subtitle)}</p>' if subtitle else ""
    return (
        f'<div class="{classes}">'
        f"<h3>{html.escape(title)}</h3>"
        f"{subtitle_html}"
        f'{"".join(chart_rows)}'
        "</div>"
    )


def svg_funnel(items: list[tuple[str, int, str]]) -> str:
    max_value = max((value for _, value, _ in items), default=1) or 1
    steps = []
    previous_value = None
    previous_label = ""
    for label, value, _color in items:
        width = max(2, 100 * value / max_value)
        if previous_value is None:
            note = "baseline"
        elif previous_value:
            note = f"{value / previous_value:.0%} of {previous_label}"
        else:
            note = f"0% of {previous_label}"
        previous_value = value
        previous_label = label
        steps.append(
            '<div class="chart-row funnel-row">'
            f'<span class="chart-label" title="{html.escape(label)}">{html.escape(label)}</span>'
            '<span class="chart-track">'
            f'<span style="width:{width:.1f}%"></span>'
            "</span>"
            f'<span class="chart-value">{html.escape(compact_number(value))}<small>{html.escape(note)}</small></span>'
            "</div>"
        )
    return (
        "<h2>2. ASV Processing Funnel</h2>"
        '<p class="note">How many sequence features survive each major step.</p>'
        '<div class="chart funnel-card">'
        f'{"".join(steps)}'
        "</div>"
    )


def report_table(
    rows: list[dict[str, object]],
    columns: list[tuple[str, str]],
    title: str,
    subtitle: str = "",
    bar_key: str = "",
) -> str:
    if not rows:
        return f"<h2>{html.escape(title)}</h2><p class=\"note\">No data available.</p>"
    header = "".join(f"<th>{html.escape(label)}</th>" for label, _ in columns)
    max_value = max((float(row.get(bar_key, 0) or 0) for row in rows), default=0) or 1.0
    body = []
    for row in rows:
        cells = []
        for index, (_, key) in enumerate(columns):
            value = row.get(key, "")
            if index == 0 and bar_key:
                bar_value = float(row.get(bar_key, 0) or 0)
                width = 100 * bar_value / max_value
                cell = (
                    f'<div class="bar-label">{html.escape(str(value))}</div>'
                    '<div class="bar">'
                    f'<span style="width:{width:.1f}%"></span>'
                    "</div>"
                )
            else:
                if isinstance(value, (int, float)):
                    value = compact_number(value)
                cell = html.escape(str(value))
            cells.append(f"<td>{cell}</td>")
        body.append(f"<tr>{''.join(cells)}</tr>")
    subtitle_html = f'<p class="note">{html.escape(subtitle)}</p>' if subtitle else ""
    return (
        f"<h2>{html.escape(title)}</h2>"
        f"{subtitle_html}"
        f'<table><thead><tr>{header}</tr></thead><tbody>{"".join(body)}</tbody></table>'
    )


def report_details(path: str, title: str, limit: int = 12, note: str = "") -> str:
    header, rows = read_tsv(path, limit=limit)
    if not header:
        return ""
    head = "".join(f"<th>{html.escape(str(col))}</th>" for col in header)
    status_idx = header.index("status") if "status" in header else -1
    message_idx = header.index("message") if "message" in header else -1
    body_rows = []
    for row in rows:
        cells = list(row)
        if status_idx >= 0 and message_idx >= 0 and ".md5 sidecar" in str(cells[message_idx]).lower():
            cells[status_idx] = "INFO"
        body_rows.append(
            "<tr>" + "".join(f"<td>{html.escape(compact_table_cell(cell))}</td>" for cell in cells) + "</tr>"
        )
    body = "".join(body_rows)
    note_html = f'<p class="detail-note">{note}</p>' if note else ""
    return (
        f"<details><summary>{html.escape(title)}</summary>"
        f"{note_html}"
        f'<div class="technical"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'
        "</details>"
    )


def cmd_make_report(args: argparse.Namespace) -> int:
    cigar_overview, cigar_amplicons, sample_depth = summarize_cigar_matrix(args.cigar, args.samples)
    mapped_overview, mapped_amplicons = summarize_mapped_asvs(args.mapped)
    asv_cigar_rows = summarize_asv_to_cigar(args.asv_to_cigar)
    try:
        cdhit_summary = json.loads(Path(args.cdhit_summary).read_text(encoding="utf-8"))
    except (OSError, TypeError, json.JSONDecodeError):
        cdhit_summary = {}
    validate_cdhit_summary(cdhit_summary)
    preflight_rows = read_tsv_dicts(args.preflight)
    status_counts = Counter(row.get("status", "") for row in preflight_rows)
    optional_md5_notices = sum(1 for row in preflight_rows if is_optional_md5_notice(row))
    filter_pass_total = sum(int(row.get("filter_pass_asvs", 0) or 0) for row in mapped_amplicons)
    locus_counts = Counter(int(row.get("loci_detected", 0) or 0) for row in sample_depth)
    total_loci = len(cigar_amplicons)
    locus_recovery_rows = [
        {"bucket": f"{count}/{total_loci}", "samples": locus_counts.get(count, 0)}
        for count in range(total_loci, -1, -1)
    ]
    samples = cigar_overview.get("samples", "NA")
    patients = cigar_overview.get("patients", "NA")
    median_loci = cigar_overview.get("median_loci_detected", "NA")
    total_final_reads = cigar_overview.get("total_final_reads", "NA")
    samples_with_reads = cigar_overview.get("samples_with_reads", "NA")
    mapped_asvs = mapped_overview.get("mapped_asvs", "NA")
    amplicons_with_asvs = mapped_overview.get("amplicons_with_asvs", "NA")
    total_haplotypes = cigar_overview.get("total_haplotypes", "NA")
    unique_cigars = sum(int(row.get("unique_cigars", 0) or 0) for row in asv_cigar_rows)
    error_count = status_counts.get("ERROR", 0)
    warning_count = sum(
        1
        for row in preflight_rows
        if row.get("status", "") == "WARN" and not is_optional_md5_notice(row)
    )
    if error_count:
        run_status = "Needs review"
    elif warning_count:
        run_status = "Passed with warnings"
    else:
        run_status = "Passed"

    def safe_int(value: object) -> int:
        try:
            return int(value)
        except (TypeError, ValueError):
            return 0

    cdhit_identity = float(cdhit_summary.get("identity_threshold", 0.989) or 0.989)
    cdhit_panel = ""
    if cdhit_summary:
        cdhit_panel = f"""
    <section class="report-section">
      <div class="section-heading">
        <div>
          <h2>{cdhit_identity:.1%} ASV sequence clusters</h2>
          <p>CD-HIT-EST clustering of ASVs that passed the same feature-level QC used before CIGAR conversion. Clustering is a separate descriptive output and does not replace CIGAR alleles.</p>
        </div>
      </div>
      <div class="cards">
        <div class="card"><b>{safe_int(cdhit_summary.get('input_asvs')):,}</b><em>ASVs clustered</em><span>filter-pass nucleotide sequences</span></div>
        <div class="card"><b>{safe_int(cdhit_summary.get('clusters')):,}</b><em>Clusters</em><span>{cdhit_identity:.1%} global identity</span></div>
        <div class="card"><b>{safe_int(cdhit_summary.get('singleton_clusters')):,}</b><em>Singleton clusters</em><span>one ASV only</span></div>
        <div class="card"><b>{safe_int(cdhit_summary.get('largest_cluster_size')):,}</b><em>Largest cluster</em><span>member ASVs</span></div>
      </div>
    </section>"""

    cards = [
        report_metric_card("Patients", compact_number(patients), "unique patient IDs", "teal"),
        report_metric_card("Samples", compact_number(samples), "final CIGAR table", "teal"),
        report_metric_card("Final reads", compact_number(total_final_reads), "post-CIGAR reads", "green"),
        report_metric_card("Median loci detected", f"{format_number(median_loci)}/{total_loci}", "per sample", "gold"),
        report_metric_card("ASVs mapped to loci", compact_number(mapped_asvs), "amplicon-mapped variants", "teal"),
        report_metric_card("CIGAR haplotypes", compact_number(total_haplotypes), f"{compact_number(unique_cigars)} unique CIGAR strings", "pink"),
    ]
    summary_path = Path(args.out).with_name("asv_filtering_summary.txt")
    write_asv_filtering_summary(
        str(summary_path),
        project_name=args.project_name,
        mapped_path=args.mapped,
        asv_to_cigar_path=args.asv_to_cigar,
        cdhit_summary=cdhit_summary,
        min_total_reads=args.min_total_reads,
        min_samples=args.min_samples,
        exclude_bimeras=args.exclude_bimeras,
    )
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Run Summary</title>
  <style>
    @import url("https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700;800&display=swap");
    * {{ box-sizing: border-box; }}
    :root {{
      --ink: #17383d;
      --muted: #356f75;
      --subtle: #5c8589;
      --line: rgba(0, 92, 104, 0.18);
      --line-dark: rgba(0, 92, 104, 0.34);
      --paper: #ffffff;
      --page: #ffffff;
      --table-head: #eef8f7;
      --panel-soft: #fbfdfd;
      --accent: #20988c;
      --accent-dark: #005c68;
      --accent-light: #dff4f1;
    }}
    body {{ font-family: "Source Sans 3", "Segoe UI", Arial, sans-serif; margin: 0; color: var(--ink); background: var(--page); font-size: 14px; line-height: 1.5; }}
    .report-shell {{ max-width: 1180px; margin: 0 auto; padding: 22px 24px 36px; }}
    h1, h2, h3, p {{ margin-top: 0; }}
    h1 {{ margin: 0; font-size: 26px; line-height: 1.15; letter-spacing: 0; font-weight: 700; }}
    h2 {{ margin: 0 0 12px; font-size: 21px; line-height: 1.25; letter-spacing: 0; font-weight: 700; }}
    h3 {{ margin: 0 0 10px; font-size: 15px; line-height: 1.3; letter-spacing: 0; font-weight: 700; }}
    .report-header {{ display: grid; grid-template-columns: minmax(0, 1fr) minmax(220px, 280px); gap: 18px; align-items: start; background: transparent; border: 0; border-bottom: 1px solid var(--line); border-radius: 0; padding: 4px 0 16px; }}
    .eyebrow {{ margin: 0 0 6px; color: var(--muted); font-size: 11px; font-weight: 700; text-transform: uppercase; }}
    .lede {{ max-width: 720px; margin: 7px 0 0; color: var(--muted); font-size: 14px; }}
    .status-card {{ border: 1px solid var(--line); background: var(--paper); border-radius: 8px; padding: 10px 12px; }}
    .status-card span, .metric span {{ display: block; color: var(--muted); font-size: 11px; font-weight: 700; text-transform: uppercase; }}
    .status-card strong {{ display: block; margin-top: 3px; font-size: 15px; line-height: 1.2; }}
    .status-card dl {{ display: grid; grid-template-columns: repeat(2, 1fr); gap: 8px 12px; margin: 9px 0 0; padding-top: 9px; border-top: 1px solid var(--line); }}
    .status-card dt {{ color: var(--subtle); font-size: 11px; text-transform: uppercase; }}
    .status-card dd {{ margin: 2px 0 0; font-weight: 700; font-variant-numeric: tabular-nums; }}
    .metrics {{ display: grid; grid-template-columns: repeat(6, minmax(0, 1fr)); gap: 10px; margin-top: 18px; }}
    .metric {{ min-height: 112px; border: 1px solid var(--line); background: var(--paper); border-radius: 8px; padding: 13px; }}
    .metric i {{ display: none; }}
    .metric b {{ display: block; font-size: 25px; line-height: 1; margin-top: 7px; font-variant-numeric: tabular-nums; }}
    .metric span + b {{ margin-top: 7px; }}
    .metric b + span {{ color: var(--muted); font-size: 12px; font-weight: 400; text-transform: none; margin-top: 8px; }}
    .report-section {{ border: 1px solid var(--line); background: var(--paper); border-radius: 8px; padding: 20px; margin-top: 18px; }}
    .note, .detail-note, .details-note, .footer-note {{ color: var(--muted); }}
    .section-heading {{ display: flex; align-items: start; justify-content: space-between; gap: 18px; margin-bottom: 16px; }}
    .section-heading p {{ margin: 6px 0 0; color: var(--muted); max-width: 760px; }}
    table {{ border-collapse: separate; border-spacing: 0; width: 100%; margin-top: 10px; font-size: 13px; background: var(--paper); border: 1px solid var(--line); border-radius: 8px; overflow: hidden; }}
    th, td {{ border: 1px solid var(--line); padding: 7px 8px; text-align: left; vertical-align: top; }}
    th {{ background: var(--table-head); color: var(--ink); font-size: 11px; text-transform: uppercase; }}
    tbody tr:nth-child(even) {{ background: var(--panel-soft); }}
    td:not(:first-child), th:not(:first-child) {{ text-align: right; font-variant-numeric: tabular-nums; }}
    .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 10px; margin-top: 16px; }}
    .card {{ border: 1px solid var(--line); background: var(--panel-soft); border-radius: 8px; padding: 12px; }}
    .card b {{ display: block; font-size: 22px; line-height: 1; color: var(--ink); font-variant-numeric: tabular-nums; }}
    .card em {{ display: block; font-style: normal; font-weight: 700; margin-top: 6px; }}
    .card span {{ display: block; color: var(--muted); font-size: 12px; margin-top: 5px; }}
    .bar {{ height: 6px; background: rgba(0, 92, 104, 0.12); border-radius: 999px; margin-top: 6px; overflow: hidden; }}
    .bar span {{ display: block; height: 100%; background: var(--accent); }}
    .bar-label {{ font-weight: 700; color: var(--ink); }}
    .viz-grid {{ display: grid; grid-template-columns: repeat(2, minmax(320px, 1fr)); gap: 14px; }}
    .chart {{ border: 1px solid var(--line); border-radius: 8px; padding: 14px; background: var(--panel-soft); }}
    .chart-wide {{ grid-column: 1 / -1; }}
    .chart-row {{ display: grid; grid-template-columns: 122px 1fr 82px; gap: 10px; align-items: center; margin: 8px 0; }}
    .chart-label {{ font-weight: 700; color: var(--ink); overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
    .chart-track {{ display: block; height: 9px; background: rgba(0, 92, 104, 0.12); border-radius: 999px; overflow: hidden; }}
    .chart-track span {{ display: block; height: 100%; background: var(--accent); }}
    .chart-value {{ text-align: right; font-variant-numeric: tabular-nums; color: var(--ink); font-weight: 700; }}
    .chart-value small {{ display: block; color: var(--subtle); font-size: 10px; font-weight: 400; }}
    .funnel-row {{ grid-template-columns: 142px 1fr 94px; }}
    .funnel-row .chart-label {{ white-space: normal; line-height: 1.15; }}
    .funnel-row .chart-track span {{ background: var(--accent); }}
    .funnel-card {{ margin-top: 0; }}
    .funnel-card .chart-row {{ grid-template-columns: 190px 1fr 150px; max-width: 1000px; }}
    .table-section {{ overflow-x: auto; }}
    .table-section table, .technical table {{ min-width: 720px; }}
    .table-section h2, .technical-section h2 {{ margin-top: 24px; padding-top: 24px; border-top: 1px solid var(--line); }}
    .table-section h2:first-child, .technical-section h2:first-child {{ margin-top: 0; padding-top: 0; border-top: 0; }}
    details {{ border: 1px solid var(--line); background: var(--panel-soft); border-radius: 8px; margin-top: 10px; padding: 0; overflow: hidden; }}
    summary {{ cursor: pointer; font-weight: 700; color: var(--ink); font-size: 15px; padding: 12px 14px; }}
    details[open] summary {{ border-bottom: 1px solid var(--line); background: var(--table-head); }}
    .technical {{ overflow: auto; max-height: 380px; padding: 0 14px 14px; }}
    .detail-note {{ margin: 14px 16px 0; }}
    .details-note {{ margin: 0 0 14px; }}
    .footer-note {{ margin: 22px 0 0; font-size: 12px; text-align: center; }}
    @media (max-width: 980px) {{ .report-header {{ grid-template-columns: 1fr; }} .metrics {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }} .viz-grid {{ grid-template-columns: 1fr; }} .section-heading {{ display: block; }} }}
    @media (max-width: 720px) {{ .report-shell {{ padding: 18px; }} .report-header {{ padding-bottom: 14px; }} .report-section {{ padding: 16px; }} .metrics {{ grid-template-columns: 1fr; }} .chart-row, .funnel-row, .funnel-card .chart-row {{ grid-template-columns: 1fr; gap: 6px; }} .chart-value {{ text-align: left; }} }}
    @media print {{
      body {{ background: #fff; color: #111; margin: 0.35in; }}
      .report-shell {{ max-width: none; padding: 0; }}
      .report-header {{ color: #111; background: #fff; }}
      .chart, .card, table, details {{ break-inside: avoid; page-break-inside: avoid; }}
    }}
  </style>
</head>
<body>
  <main class="report-shell">
    <header class="report-header">
      <div>
        <p class="eyebrow">SIMPLseq workflow report</p>
        <h1>Run Summary</h1>
        <p class="lede">Quality-control summary of read depth, locus recovery, ASV filtering, and CIGAR conversion for this SIMPLseq run.</p>
      </div>
      <aside class="status-card" aria-label="Run status">
        <span>Run status</span>
        <strong>{html.escape(run_status)}</strong>
        <dl>
          <div><dt>Samples</dt><dd>{html.escape(compact_number(samples))}</dd></div>
          <div><dt>Patients</dt><dd>{html.escape(compact_number(patients))}</dd></div>
          <div><dt>Errors</dt><dd>{error_count}</dd></div>
          <div><dt>Warnings</dt><dd>{warning_count}</dd></div>
        </dl>
      </aside>
    </header>
    <section class="metrics" aria-label="Run metrics">{''.join(cards)}</section>
    <section class="report-section">
      <div class="section-heading">
        <div>
          <h2>1. Visual Summary</h2>
          <p>Read-depth and recovery summaries for rapid review before inspecting the tabular outputs.</p>
        </div>
      </div>
      <div class="viz-grid">
        {svg_bar_chart(cigar_amplicons, "amplicon", "final_reads", "Final reads by amplicon")}
        {svg_bar_chart(mapped_amplicons, "amplicon", "filter_pass_asvs", "Filter-pass ASVs by amplicon")}
        {svg_bar_chart(locus_recovery_rows, "bucket", "samples", "Sample locus recovery")}
        {svg_bar_chart(asv_cigar_rows, "amplicon", "unique_cigars", "Unique CIGAR strings")}
      </div>
    </section>
    <section class="report-section">
      {svg_funnel([("Mapped ASVs", safe_int(mapped_asvs), "#343832"), ("Filter-pass ASVs", filter_pass_total, "#343832"), ("CIGAR haplotypes", safe_int(total_haplotypes), "#343832")])}
    </section>
    {cdhit_panel}
    <section class="report-section table-section">
      {report_table(cigar_amplicons, [("Amplicon", "amplicon"), ("Final reads", "final_reads"), ("Samples >=100 reads", "samples_detected"), ("Haplotypes", "haplotypes"), ("Median reads in detected samples", "median_reads_detected_samples")], "3. Final Read Depth by Amplicon", "Final CIGAR/haplotype reads grouped by SIMPLseq target. Samples detected uses a 100-read threshold, matching the paper's locus-detection summaries.", "final_reads")}
    </section>
    <section class="report-section table-section">
      {report_table(mapped_amplicons, [("Amplicon", "amplicon"), ("Mapped reads", "mapped_reads"), ("Mapped ASVs", "mapped_asvs"), ("Filter-pass ASVs", "filter_pass_asvs"), ("Max samples per ASV", "max_samples_per_asv")], "4. Mapped ASVs by Amplicon", "", "mapped_reads")}
    </section>
    <section class="report-section table-section">
      {report_table(sample_depth, [("Sample", "sample"), ("Final reads", "total_reads"), ("Loci detected", "loci_detected"), ("Haplotypes detected", "haplotypes_detected")], "5. Sample Read Depth and Locus Recovery", "", "total_reads")}
    </section>
    <section class="report-section table-section">
      {report_table(asv_cigar_rows, [("Amplicon", "amplicon"), ("ASVs with CIGAR", "asvs_with_cigar"), ("Unique CIGAR strings", "unique_cigars")], "6. ASV to CIGAR Conversion", "", "asvs_with_cigar")}
    </section>
    <section class="report-section technical-section">
      <h2>7. Technical Previews</h2>
      <p class="details-note">Reference and output table previews are collapsed by default. Full TSV files are saved in the run folder.</p>
      {report_details(args.geometry, "Amplicon reference check")}
      {report_details(args.cigar_summary, "CIGAR input check")}
      {report_details(args.mapped, "Mapped ASVs", limit=20, note='Preview only. The full table is saved in the run folder.')}
      {report_details(args.cigar, "CIGAR count matrix", limit=6, note='Preview only. The full table is saved in the run folder.')}
    </section>
    <p class="footer-note">Generated by the Nextflow SIMPLseq workflow. {error_count} errors, {warning_count} warnings, {optional_md5_notices} optional checksum notices.</p>
  </main>
</body>
</html>
"""
    with open(args.out, "w") as handle:
        handle.write(html_text)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("preflight")
    p.add_argument("--samples", required=True)
    p.add_argument("--amplicons", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--geometry", required=True)
    p.add_argument("--barcode", required=True)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--inline-barcodes-enabled", default="false")
    p.add_argument("--sentinel-locus", default="KELT")
    p.add_argument("--samples-root", default="")
    p.set_defaults(func=cmd_preflight)

    p = sub.add_parser("write-meta")
    p.add_argument("--samples", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--samples-root", default="")
    p.set_defaults(func=cmd_write_meta)

    p = sub.add_parser("write-pipeline-json")
    p.add_argument("--meta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--pipeline-class", default="parasite")
    p.add_argument("--max-ee", required=True)
    p.add_argument("--trim-right", required=True)
    p.add_argument("--min-len", required=True)
    p.add_argument("--trunc-q", required=True)
    p.add_argument("--max-consist", required=True)
    p.add_argument("--omega-a", required=True)
    p.add_argument("--just-concatenate", required=True)
    p.add_argument("--save-rdata", default="")
    p.add_argument("--dada2-randomize", default="")
    p.add_argument("--dada2-multithread", default="")
    p.add_argument("--dada2-seed", default="")
    p.add_argument("--primers-fwd", required=True)
    p.add_argument("--primers-rev", required=True)
    p.add_argument("--overlap-primers-fwd", required=True)
    p.add_argument("--overlap-primers-rev", required=True)
    p.set_defaults(func=cmd_write_pipeline_json)

    p = sub.add_parser("prepare-stage2")
    p.add_argument("--seqtab", required=True)
    p.add_argument("--op-bimera", required=True)
    p.add_argument("--nop-bimera", required=True)
    p.add_argument("--corrected-asv", required=True)
    p.add_argument("--strict-seqtab", required=True)
    p.add_argument("--strict-bimera", required=True)
    p.add_argument("--strict-min-asv-length", type=int, default=100)
    p.set_defaults(func=cmd_prepare_stage2)

    p = sub.add_parser("prepare-stage3")
    p.add_argument("--mapped", required=True)
    p.add_argument("--seqtab", required=True)
    p.add_argument("--filtered-mapped", required=True)
    p.add_argument("--fixed-seqtab", required=True)
    p.set_defaults(func=cmd_prepare_stage3)

    p = sub.add_parser("check-cigar-inputs")
    p.add_argument("--table", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--min-reads", type=int, default=100)
    p.add_argument("--min-samples", type=int, default=2)
    p.add_argument("--include-failed", action="store_true")
    p.add_argument("--exclude-bimeras", action="store_true")
    p.set_defaults(func=cmd_check_cigar_inputs)

    p = sub.add_parser("make-report")
    p.add_argument("--project-name", default="SIMPLseq")
    p.add_argument("--preflight", required=True)
    p.add_argument("--geometry", required=True)
    p.add_argument("--cigar-summary", required=True)
    p.add_argument("--mapped", required=True)
    p.add_argument("--asv-to-cigar", required=True)
    p.add_argument("--cigar", required=True)
    p.add_argument("--cdhit-summary", required=True)
    p.add_argument("--min-total-reads", type=int, default=100)
    p.add_argument("--min-samples", type=int, default=2)
    p.add_argument("--exclude-bimeras", action="store_true")
    p.add_argument("--samples", default="")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_make_report)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
