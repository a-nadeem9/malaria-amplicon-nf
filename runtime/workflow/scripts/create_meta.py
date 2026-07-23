#!/usr/bin/env python3
"""Create the three-column FASTQ manifest consumed by AmpliconPipeline."""

from __future__ import annotations

import argparse
import csv
import fnmatch
from pathlib import Path


def pattern_suffix(pattern: str) -> str:
    """Return the literal suffix from the simple glob patterns used here."""
    if not pattern.startswith("*") or any(char in pattern[1:] for char in "*?["):
        raise ValueError(f"FASTQ pattern must be a single leading '*' plus a literal suffix: {pattern}")
    suffix = pattern[1:]
    if not suffix:
        raise ValueError("FASTQ pattern suffix cannot be empty.")
    return suffix


def discover_pairs(path: Path, forward_pattern: str, reverse_pattern: str) -> list[tuple[str, Path, Path]]:
    if not path.is_dir():
        raise ValueError(f"FASTQ directory does not exist: {path}")
    forward_suffix = pattern_suffix(forward_pattern)
    reverse_suffix = pattern_suffix(reverse_pattern)
    pairs: list[tuple[str, Path, Path]] = []
    missing_reverse: list[str] = []
    seen_ids: set[str] = set()

    for forward in sorted(path.iterdir(), key=lambda item: item.name.lower()):
        if not forward.is_file() or not fnmatch.fnmatch(forward.name, forward_pattern):
            continue
        sample_id = forward.name[: -len(forward_suffix)]
        if not sample_id:
            raise ValueError(f"Could not derive a sample ID from {forward.name}")
        if sample_id in seen_ids:
            raise ValueError(f"Duplicate sample ID generated from FASTQ names: {sample_id}")
        reverse = path / f"{sample_id}{reverse_suffix}"
        if not reverse.is_file():
            missing_reverse.append(reverse.name)
            continue
        seen_ids.add(sample_id)
        pairs.append((sample_id, forward.resolve(), reverse.resolve()))

    if missing_reverse:
        preview = ", ".join(missing_reverse[:5])
        remainder = f" and {len(missing_reverse) - 5} more" if len(missing_reverse) > 5 else ""
        raise ValueError(f"Missing reverse FASTQ file(s): {preview}{remainder}")
    if not pairs:
        raise ValueError(f"No paired FASTQ files matched {forward_pattern} and {reverse_pattern} in {path}")
    return pairs


def write_manifest(pairs: list[tuple[str, Path, Path]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows((sample_id, str(forward), str(reverse)) for sample_id, forward, reverse in pairs)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path_to_fq", required=True, help="Directory containing FASTQ files")
    parser.add_argument("--output_file", required=True, help="Output three-column manifest")
    parser.add_argument("--pattern_fw", required=True, help="Forward FASTQ glob pattern")
    parser.add_argument("--pattern_rv", required=True, help="Reverse FASTQ glob pattern")
    args = parser.parse_args()

    try:
        pairs = discover_pairs(Path(args.path_to_fq), args.pattern_fw, args.pattern_rv)
        output = Path(args.output_file)
        write_manifest(pairs, output)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"Meta file generated at {output} ({len(pairs)} paired samples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
