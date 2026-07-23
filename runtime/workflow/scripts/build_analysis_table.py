#!/usr/bin/env python3
"""Build the shared allele-count table used by downstream analyses."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing input table: {path}")
    table = pd.read_csv(path, sep="\t")
    if table.empty:
        raise ValueError(f"Input table is empty: {path}")
    return table


def build_off_mode(cigar: Path, output: Path) -> dict[str, object]:
    if not cigar.exists():
        raise FileNotFoundError(f"Missing CIGAR count table: {cigar}")
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(cigar, output)
    table = read_table(output)
    return {
        "mode": "off",
        "source": str(cigar),
        "output": str(output),
        "samples": int(len(table)),
        "alleles": int(max(table.shape[1] - 1, 0)),
        "description": "Exact CIGAR allele table copied for downstream analysis.",
    }


def cluster_locus_map(membership: pd.DataFrame) -> dict[str, str]:
    required = {"cluster_id", "amplicon"}
    missing = required - set(membership.columns)
    if missing:
        raise ValueError(f"CD-HIT membership is missing columns: {', '.join(sorted(missing))}")
    cleaned = membership.assign(
        cluster_id=membership["cluster_id"].astype(str).str.strip(),
        amplicon=membership["amplicon"].astype(str).str.strip(),
    )
    bad = cleaned.loc[(cleaned["cluster_id"] == "") | (cleaned["amplicon"] == "")]
    if not bad.empty:
        raise ValueError("CD-HIT membership has clusters with missing cluster_id or amplicon.")
    locus_counts = cleaned.groupby("cluster_id")["amplicon"].nunique()
    mixed = locus_counts[locus_counts > 1]
    if not mixed.empty:
        preview = ", ".join(mixed.index.astype(str).tolist()[:8])
        raise ValueError(
            "CD-HIT clusters cross loci, so a downstream locus cannot be assigned safely: "
            f"{preview}"
        )
    return cleaned.groupby("cluster_id")["amplicon"].first().to_dict()


def build_summed_mode(counts_path: Path, membership_path: Path, output: Path) -> dict[str, object]:
    counts = read_table(counts_path)
    membership = read_table(membership_path)
    sample_column = str(counts.columns[0])
    locus_for_cluster = cluster_locus_map(membership)
    cluster_columns = [str(column) for column in counts.columns[1:]]
    missing = [cluster for cluster in cluster_columns if cluster not in locus_for_cluster]
    if missing:
        raise ValueError(
            "CD-HIT count table has clusters missing from membership: "
            + ", ".join(missing[:8])
        )

    output_table = pd.DataFrame({"sample": counts[sample_column].astype(str)})
    for cluster_id in sorted(cluster_columns):
        locus = locus_for_cluster[cluster_id]
        output_table[f"{locus},{cluster_id}"] = pd.to_numeric(
            counts[cluster_id], errors="coerce"
        ).fillna(0).astype(int)

    output.parent.mkdir(parents=True, exist_ok=True)
    output_table.to_csv(output, sep="\t", index=False)
    return {
        "mode": "cdhit_summed",
        "source_counts": str(counts_path),
        "source_membership": str(membership_path),
        "output": str(output),
        "samples": int(len(output_table)),
        "alleles": int(max(output_table.shape[1] - 1, 0)),
        "clusters": int(len(cluster_columns)),
        "description": "CD-HIT clusters represented as summed cluster-count alleles.",
    }


def build(args: argparse.Namespace) -> dict[str, object]:
    mode = str(args.mode).strip().lower()
    output = Path(args.out).resolve()
    if mode == "off":
        summary = build_off_mode(Path(args.cigar).resolve(), output)
    elif mode == "summed":
        summary = build_summed_mode(
            Path(args.cdhit_counts).resolve(),
            Path(args.cdhit_membership).resolve(),
            output,
        )
    else:
        raise ValueError(f"Unsupported analysis table mode: {args.mode}")

    if args.summary:
        summary_path = Path(args.summary).resolve()
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=["off", "summed"], required=True)
    p.add_argument("--cigar", required=True)
    p.add_argument("--cdhit-counts", default="")
    p.add_argument("--cdhit-membership", default="")
    p.add_argument("--out", required=True)
    p.add_argument("--summary", default="")
    return p


def main() -> int:
    summary = build(parser().parse_args())
    print(
        f"Built analysis table ({summary['mode']}): "
        f"{summary['samples']} samples x {summary['alleles']} alleles"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
