#!/usr/bin/env python
"""Detect unexpected well-specific KELT inline-barcode pairs."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass
class SampleResult:
    sample_id: str
    library: str
    expected_forward: str
    expected_reverse: str
    merged_reads: int
    kelt_pair_reads: int
    expected_pair_reads: int
    unexpected_pair_reads: int
    partial_barcode_reads: int
    ambiguous_barcode_reads: int
    unexpected_fraction: float
    qc_status: str
    signal_band: str


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def hamming_within(left: str, right: str, maximum: int) -> bool:
    return len(left) == len(right) and sum(a != b for a, b in zip(left, right)) <= maximum


def matching_prefixes(sequence: str, barcodes: list[str], mismatches: int) -> list[str]:
    return [barcode for barcode in barcodes if hamming_within(sequence[: len(barcode)], barcode, mismatches)]


def matching_suffixes(sequence: str, reverse_barcodes: list[str], mismatches: int) -> list[str]:
    return [
        barcode
        for barcode in reverse_barcodes
        if hamming_within(sequence[-len(barcode) :], reverse_complement(barcode), mismatches)
    ]


def iter_fastq_sequences(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline().strip().upper()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"Truncated FASTQ record in {path}")
            yield sequence


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {"sample_id", "expected_fwd_barcode", "expected_rev_barcode"}
    missing = sorted(required.difference(rows[0].keys() if rows else set()))
    if missing:
        raise ValueError("KELT-enabled sample sheet is missing columns: " + ", ".join(missing))
    incomplete = [
        str(row.get("sample_id", ""))
        for row in rows
        if not str(row.get("expected_fwd_barcode", "")).strip()
        or not str(row.get("expected_rev_barcode", "")).strip()
    ]
    if incomplete:
        raise ValueError("Expected KELT barcode pair is missing for: " + ", ".join(incomplete[:20]))
    return rows


def count_sample(
    row: dict[str, str],
    fastq_dir: Path,
    temp_dir: Path,
    forward_barcodes: list[str],
    reverse_barcodes: list[str],
    pair_owners: dict[tuple[str, str], list[str]],
    mismatches: int,
    trace_max_reads: int,
) -> tuple[SampleResult, list[dict[str, object]]]:
    sample_id = str(row["sample_id"]).strip()
    expected_forward = str(row["expected_fwd_barcode"]).strip().upper()
    expected_reverse = str(row["expected_rev_barcode"]).strip().upper()
    read_1 = fastq_dir / f"{sample_id}_val_1.fq.gz"
    read_2 = fastq_dir / f"{sample_id}_val_2.fq.gz"
    if not read_1.exists() or not read_2.exists():
        raise ValueError(f"Adapter-trimmed FASTQs were not found for {sample_id}")
    merged = temp_dir / f"{sample_id}.merged.fastq.gz"
    command = [
        "bbmerge.sh",
        f"in1={read_1}",
        f"in2={read_2}",
        f"out={merged}",
        "threads=1",
        "overwrite=t",
    ]
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        detail = " ".join(completed.stdout.splitlines()[-8:])
        raise RuntimeError(f"BBMerge failed for {sample_id}: {detail}")

    observed_pairs: Counter[tuple[str, str]] = Counter()
    partial = 0
    ambiguous = 0
    merged_reads = 0
    try:
        for sequence in iter_fastq_sequences(merged):
            merged_reads += 1
            forward_matches = matching_prefixes(sequence, forward_barcodes, mismatches)
            reverse_matches = matching_suffixes(sequence, reverse_barcodes, mismatches)
            if len(forward_matches) == 1 and len(reverse_matches) == 1:
                observed_pairs[(forward_matches[0], reverse_matches[0])] += 1
            elif len(forward_matches) > 1 or len(reverse_matches) > 1:
                ambiguous += 1
            elif forward_matches or reverse_matches:
                partial += 1
    finally:
        merged.unlink(missing_ok=True)

    expected_reads = observed_pairs.get((expected_forward, expected_reverse), 0)
    unexpected_reads = sum(
        count for pair, count in observed_pairs.items() if pair != (expected_forward, expected_reverse)
    )
    kelt_pair_reads = sum(observed_pairs.values())
    if kelt_pair_reads == 0:
        qc_status = "not_evaluable"
        signal_band = "no_kelt_pairs"
    elif unexpected_reads == 0:
        qc_status = "clear"
        signal_band = "none"
    elif unexpected_reads <= trace_max_reads:
        qc_status = "review"
        signal_band = "trace"
    else:
        qc_status = "review"
        signal_band = "detected"
    result = SampleResult(
        sample_id=sample_id,
        library=str(row.get("library", "") or ""),
        expected_forward=expected_forward,
        expected_reverse=expected_reverse,
        merged_reads=merged_reads,
        kelt_pair_reads=kelt_pair_reads,
        expected_pair_reads=expected_reads,
        unexpected_pair_reads=unexpected_reads,
        partial_barcode_reads=partial,
        ambiguous_barcode_reads=ambiguous,
        unexpected_fraction=(unexpected_reads / kelt_pair_reads) if kelt_pair_reads else 0.0,
        qc_status=qc_status,
        signal_band=signal_band,
    )
    counts = []
    for (forward, reverse), count in observed_pairs.most_common():
        expected = forward == expected_forward and reverse == expected_reverse
        owners = pair_owners.get((forward, reverse), [])
        counts.append(
            {
                "sample_id": sample_id,
                "observed_forward": forward,
                "observed_reverse": reverse,
                "read_pairs": count,
                "pair_status": "expected" if expected else "unexpected",
                "mapped_samples": ";".join(owners),
            }
        )
    return result, counts


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, summary: dict[str, object], results: list[SampleResult]) -> None:
    review_rows = [result for result in results if result.qc_status == "review"]
    table_rows = "".join(
        "<tr>"
        f"<td>{html.escape(result.sample_id)}</td>"
        f"<td>{html.escape(result.library or '--')}</td>"
        f"<td>{result.expected_pair_reads:,}</td>"
        f"<td>{result.unexpected_pair_reads:,}</td>"
        f"<td>{result.unexpected_fraction:.2%}</td>"
        f"<td><span class='status {html.escape(result.signal_band)}'>{html.escape(result.signal_band.replace('_', ' '))}</span></td>"
        "</tr>"
        for result in (review_rows or results)
    )
    path.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<title>KELT contamination QC</title><style>"
        ":root{font-family:Arial,sans-serif;color:#12383d;background:#fff}body{margin:0;padding:32px}"
        "main{max-width:1060px;margin:auto}h1{font-size:28px;font-weight:500;margin:0 0 8px}p{color:#526c70;line-height:1.55}"
        ".metrics{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));border-block:1px solid #d7e3e4;margin:26px 0}"
        ".metric{padding:18px 12px}.metric span{display:block;color:#61777a;font-size:11px;text-transform:uppercase}.metric strong{font-size:24px;font-weight:500}"
        "table{width:100%;border-collapse:collapse;font-size:13px}th,td{padding:11px 12px;border-bottom:1px solid #e4ecec;text-align:left}th{color:#61777a;font-size:11px;font-weight:500;text-transform:uppercase}"
        ".status{font-size:11px;text-transform:capitalize}.trace{color:#9b6904}.detected{color:#d12b4f}.none{color:#087c62}.no_kelt_pairs{color:#61777a}"
        ".note{border-left:3px solid #159d91;padding:2px 0 2px 14px;margin:24px 0}.note strong{font-weight:500}"
        "@media(max-width:700px){body{padding:20px}.metrics{grid-template-columns:repeat(2,1fr)}table{display:block;overflow:auto}}"
        "</style></head><body><main>"
        "<h1>KELT contamination QC</h1>"
        "<p>Unexpected well-specific KELT inline-barcode pairs detected after adapter trimming and paired-read merging.</p>"
        "<div class='metrics'>"
        f"<div class='metric'><span>Samples checked</span><strong>{summary['samples_checked']}</strong></div>"
        f"<div class='metric'><span>Clear</span><strong>{summary['samples_clear']}</strong></div>"
        f"<div class='metric'><span>Trace signal</span><strong>{summary['samples_trace']}</strong></div>"
        f"<div class='metric'><span>Detected signal</span><strong>{summary['samples_detected']}</strong></div>"
        "</div>"
        "<div class='note'><strong>Interpretation</strong><p>Any unexpected barcode pair is retained for review. One or two reads are labelled trace signal; three or more are labelled detected signal. KELT is a PCR1 sentinel and does not exclude contamination at unbarcoded loci or contamination introduced before barcoding.</p></div>"
        f"<h2>{'Samples requiring review' if review_rows else 'All samples'}</h2>"
        "<table><thead><tr><th>Sample</th><th>Library</th><th>Expected pairs</th><th>Unexpected pairs</th><th>Unexpected fraction</th><th>Signal</th></tr></thead>"
        f"<tbody>{table_rows}</tbody></table>"
        "</main></body></html>",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq-dir", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--mismatches", type=int, default=0)
    parser.add_argument("--trace-max-reads", type=int, default=2)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.mismatches < 0 or args.mismatches > 2:
        parser.error("--mismatches must be 0, 1, or 2")
    if not shutil.which("bbmerge.sh"):
        parser.error("bbmerge.sh is required for KELT contamination QC")

    samples = read_samples(Path(args.samples))
    fastq_dir = Path(args.fastq_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    forward_barcodes = sorted({str(row["expected_fwd_barcode"]).strip().upper() for row in samples}, key=len, reverse=True)
    reverse_barcodes = sorted({str(row["expected_rev_barcode"]).strip().upper() for row in samples}, key=len, reverse=True)
    pair_owners: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in samples:
        pair = (str(row["expected_fwd_barcode"]).strip().upper(), str(row["expected_rev_barcode"]).strip().upper())
        pair_owners[pair].append(str(row["sample_id"]))

    results: list[SampleResult] = []
    count_rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="kelt-", dir=outdir) as temporary:
        with ThreadPoolExecutor(max_workers=max(1, min(args.workers, len(samples)))) as executor:
            futures = [
                executor.submit(
                    count_sample,
                    row,
                    fastq_dir,
                    Path(temporary),
                    forward_barcodes,
                    reverse_barcodes,
                    pair_owners,
                    args.mismatches,
                    args.trace_max_reads,
                )
                for row in samples
            ]
            for future in as_completed(futures):
                result, counts = future.result()
                results.append(result)
                count_rows.extend(counts)
    results.sort(key=lambda item: item.sample_id)
    count_rows.sort(key=lambda item: (str(item["sample_id"]), -int(item["read_pairs"])))

    call_rows = [asdict(result) for result in results]
    write_tsv(outdir / "kelt_contamination_calls.tsv", call_rows, list(SampleResult.__dataclass_fields__))
    write_tsv(
        outdir / "kelt_barcode_counts.tsv",
        count_rows,
        ["sample_id", "observed_forward", "observed_reverse", "read_pairs", "pair_status", "mapped_samples"],
    )
    summary = {
        "samples_checked": len(results),
        "samples_clear": sum(result.qc_status == "clear" for result in results),
        "samples_review": sum(result.qc_status == "review" for result in results),
        "samples_trace": sum(result.signal_band == "trace" for result in results),
        "samples_detected": sum(result.signal_band == "detected" for result in results),
        "samples_not_evaluable": sum(result.qc_status == "not_evaluable" for result in results),
        "unexpected_barcode_pair_reads": sum(result.unexpected_pair_reads for result in results),
        "matching_mode": "exact" if args.mismatches == 0 else f"up_to_{args.mismatches}_mismatches",
        "trace_max_reads": args.trace_max_reads,
    }
    (outdir / "kelt_contamination_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_report(outdir / "kelt_contamination_report.html", summary, results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
