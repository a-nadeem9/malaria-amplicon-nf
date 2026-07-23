"""Command line interface for malaria-amplicon-nf."""

from __future__ import annotations

import argparse
import datetime as dt
import importlib.util
import json
import os
import re
import socket
import sys
import time
from pathlib import Path

from .job_state import read_json, write_json
from .pathutils import user_path
from .resources import human_bytes
from .runner import check_environment, progress_summary, project_root, results_manifest, run_nextflow
from .samplesheet import DEFAULT_COLLECTION_DAY, scan_fastqs, write_samples_csv
from . import __version__


APP_NAME = "malaria-amplicon-nf"
APP_VERSION = os.environ.get(
    "SIMPLSEQ_VERSION",
    "v1.0" if __version__ == "1.0.0" else f"v{__version__}",
)
APP_SUBTITLE = "Linux / WSL / macOS browser workflow"


def help_description() -> str:
    line = "=" * 54
    return "\n".join(
        [
            line,
            f"  >_ {APP_NAME} {APP_VERSION}",
            f"     {APP_SUBTITLE}",
            "     Nextflow + Conda/Mamba runtime",
            line,
        ]
    )


def use_color() -> bool:
    return os.environ.get("NO_COLOR") is None and (sys.stdout.isatty() or os.environ.get("FORCE_COLOR"))


def color(text: str, code: str) -> str:
    if not use_color():
        return text
    return f"\033[{code}m{text}\033[0m"


def tag(label: str, code: str) -> str:
    return color(f"[{label}]", code)


def print_banner(title: str, subtitle: str = "") -> None:
    line = "=" * 54
    print(color(line, "90"))
    print(f"  >_ {color(f'{APP_NAME} {APP_VERSION}', '1;37')}")
    print(f"     {color(APP_SUBTITLE, '36')}")
    if title != APP_NAME:
        print(f"     {color(title, '1;37')}")
    if subtitle:
        print(f"     {color(subtitle, '36')}")
    print(color(line, "90"))
    print()


def print_check_rows(rows: list[dict[str, str]]) -> int:
    failed = 0
    for row in rows:
        status = row["status"]
        if status == "ok":
            marker = tag("OK", "32") + "  "
        elif status == "warn":
            marker = tag("INFO", "34")
        else:
            marker = tag("MISS", "31")
        if status not in {"ok", "warn"}:
            failed += 1
        name = color(row["name"], "1")
        print(f"{marker} {name}: {row['detail']}")
    return failed


def collection_year_arg(value: str) -> str:
    raw = str(value or "").strip()
    if raw and not re.fullmatch(r"(19|20)[0-9]{2}", raw):
        raise argparse.ArgumentTypeError("use a four-digit year from 1900 to 2099")
    return raw


def collection_day_arg(value: str) -> str:
    raw = str(value or "").strip()
    if not re.fullmatch(r"[0-9]{1,2}", raw):
        raise argparse.ArgumentTypeError("use a day from 1 to 31")
    day = int(raw)
    if day < 1 or day > 31:
        raise argparse.ArgumentTypeError("use a day from 1 to 31")
    return f"{day:02d}"


def cmd_scan(args: argparse.Namespace) -> int:
    print_banner("FASTQ pairing", "Sample sheet preparation")
    scan = scan_fastqs(
        args.fastq_dir,
        include_pool_in_sample_id=args.include_pool_in_sample_id,
        fallback_collection_year=args.fallback_collection_year,
        fallback_collection_day=args.fallback_collection_day,
    )
    count, duplicates = write_samples_csv(
        args.fastq_dir,
        args.out,
        include_pool_in_sample_id=args.include_pool_in_sample_id,
        absolute=args.absolute,
        fallback_collection_year=args.fallback_collection_year,
        fallback_collection_day=args.fallback_collection_day,
        metadata_path=args.metadata or None,
        metadata_sheet=args.metadata_sheet,
    )
    if duplicates:
        print(f"{tag('ERROR', '31')} Duplicate sample IDs found:")
        for item in duplicates[:20]:
            print(f"  {item}")
        print(f"{tag('INFO', '34')} Re-run with --include-pool-in-sample-id or edit the sample names.")
        return 2
    print(f"{tag('OK', '32')} Wrote {count} sample rows to {Path(args.out).resolve()}")
    print(f"{tag('INFO', '34')} FASTQ pairs:       {len(scan.pairs)}")
    print(f"{tag('INFO', '34')} Missing R2 mates:  {len(scan.missing_r2)}")
    print(f"{tag('INFO', '34')} Orphan R2 files:   {len(scan.orphan_r2)}")
    print(f"{tag('INFO', '34')} Duplicate IDs:     {len(scan.duplicate_sample_ids)}")
    print(f"{tag('INFO', '34')} MD5 files:         {scan.md5_files}")
    print(f"{tag('INFO', '34')} Total FASTQ size:  {human_bytes(scan.total_fastq_bytes)}")
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    print_banner("Runtime checks", "Python / R / DADA2 / Nextflow")
    root = project_root()
    samples = user_path(args.samples).resolve() if args.samples else None
    outdir = user_path(args.out).resolve()
    rows = check_environment(root, samples, outdir=outdir)
    failed = print_check_rows(rows)
    if failed:
        print(f"\n{tag('ERROR', '31')} {failed} checks need attention before a full run.")
        return 1
    print(f"\n{tag('OK', '32')} malaria-amplicon-nf environment looks ready.")
    return 0


def cmd_run_direct(args: argparse.Namespace) -> int:
    outdir = user_path(args.out).resolve()
    worker_token = str(getattr(args, "worker_token", "") or "").strip()
    worker_state_raw = str(getattr(args, "worker_state", "") or "").strip()
    worker_state = user_path(worker_state_raw).resolve() if worker_state_raw else None
    if bool(worker_token) != bool(worker_state):
        raise ValueError("Internal run worker token and state path must be provided together")
    if worker_state is not None:
        try:
            worker_state.relative_to(outdir)
        except ValueError as exc:
            raise ValueError("Internal run worker state must stay inside the run output folder") from exc
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            parent_state = read_json(worker_state)
            if parent_state.get("worker_pid") == os.getpid() and parent_state.get("worker_token") == worker_token:
                break
            time.sleep(0.02)
        started_at = read_json(worker_state).get("started_at") or dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        write_json(worker_state, {
            "schema_version": 1,
            "status": "running",
            "started_at": started_at,
            "worker_pid": os.getpid(),
            "worker_token": worker_token,
            "outdir": str(outdir),
        })

    try:
        result = run_nextflow(
            user_path(args.samples),
            outdir,
            profile=args.profile,
            resume=not args.no_resume,
            work_dir=user_path(args.work_dir) if args.work_dir else None,
            dry_run=args.dry_run,
            cpus=args.cpus,
            memory=args.memory,
            kelt_enabled=args.kelt,
            kelt_barcode_map=user_path(args.kelt_barcode_map) if args.kelt_barcode_map else None,
        )
    except Exception as exc:
        completed_at = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        if worker_state is not None:
            write_json(worker_state, {
                "schema_version": 1,
                "status": "failed",
                "started_at": read_json(worker_state).get("started_at"),
                "completed_at": completed_at,
                "worker_token": worker_token,
                "outdir": str(outdir),
                "detail": str(exc),
            })
        run_state = read_json(outdir / "run_state.json")
        write_json(outdir / "run_state.json", {
            **run_state,
            "status": "failed",
            "completed_at": completed_at,
            "outdir": str(outdir),
            "detail": str(exc),
        })
        print(f"{tag('ERROR', '31')} {exc}", file=sys.stderr)
        return 1

    if worker_state is not None:
        write_json(worker_state, {
            "schema_version": 1,
            "status": "complete" if result.returncode == 0 else "failed",
            "started_at": read_json(worker_state).get("started_at"),
            "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "returncode": result.returncode,
            "worker_token": worker_token,
            "outdir": str(outdir),
        })
    print(f"{tag('INFO', '34')} Technical log: {result.technical_log}")
    return result.returncode


def cmd_status(args: argparse.Namespace) -> int:
    outdir = user_path(args.out).resolve()
    state_file = outdir / "run_state.json"
    if state_file.exists():
        print(state_file.read_text(encoding="utf-8", errors="replace"))
    summary = progress_summary(outdir)
    print(f"{tag('INFO', '34')} Progress: {summary['completed_stages']}/{summary['total_stages']} stages")
    print(f"{tag('INFO', '34')} Current:  {summary['current_stage']}")
    return 0


def cmd_results(args: argparse.Namespace) -> int:
    manifest = results_manifest(user_path(args.out))
    if args.json:
        print(json.dumps(manifest, indent=2))
        return 0
    print_banner("Run outputs", "Reports and final tables")
    print(f"{tag('INFO', '34')} Output folder: {manifest['outdir']}")
    state = manifest.get("state") or {}
    if state:
        print(f"{tag('INFO', '34')} Run status:    {state.get('status', 'unknown')}")
    print()
    missing = 0
    for row in manifest["files"]:
        marker = tag("OK", "32") + "  " if row["exists"] else tag("MISS", "31")
        if not row["exists"]:
            missing += 1
        print(f"{marker} {row['label']}: {row['path']}")
    return 1 if missing else 0


def cmd_app(args: argparse.Namespace) -> int:
    return cmd_flask_app(args)


def cmd_flask_app(args: argparse.Namespace) -> int:
    if args.host not in {"127.0.0.1", "localhost", "::1"}:
        print(f"{tag('ERROR', '31')} malaria-amplicon-nf only serves the browser UI on a loopback host.", file=sys.stderr)
        return 1
    root = project_root()
    app = root / "gui" / "flask_app.py"
    if not app.exists():
        print(f"Flask app not found: {app}", file=sys.stderr)
        return 1
    port = find_free_port(args.port)
    if port is None:
        print(
            f"{tag('ERROR', '31')} No free port found from {args.port} to {args.port + 49}.",
            file=sys.stderr,
        )
        print(f"{tag('INFO', '34')} Stop another app or run: simplseq run --port 8600", file=sys.stderr)
        return 1
    if port != args.port:
        print(f"{tag('INFO', '34')} Port {args.port} is busy; using {port} instead.")
        print(f"{tag('INFO', '34')} Opening malaria-amplicon-nf at http://{args.host}:{port}")

    try:
        spec = importlib.util.spec_from_file_location("simplseq_flask_app", app)
        if spec is None or spec.loader is None:
            raise ImportError(f"could not load {app}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except ImportError as exc:
        print(f"{tag('ERROR', '31')} Flask GUI dependencies are not available: {exc}", file=sys.stderr)
        print(f"{tag('INFO', '34')} Re-run the installer or install the managed runtime.", file=sys.stderr)
        return 1
    return int(module.run_server(root=root, host=args.host, port=port, open_browser=not args.no_browser))


def find_free_port(start: int, attempts: int = 50) -> int | None:
    for port in range(start, start + attempts):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.settimeout(0.2)
            try:
                sock.bind(("127.0.0.1", port))
            except OSError:
                continue
            return port
    return None


def add_direct_run_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--samples", required=True)
    parser.add_argument("--out", "--outdir", dest="out", default="results")
    parser.add_argument("--profile", choices=["local", "reproducible"], default="local")
    parser.add_argument("--work-dir", default="")
    parser.add_argument("--cpus", type=int, default=0, help="CPUs for heavy local stages")
    parser.add_argument("--memory", default="", help="Memory for heavy local stages, e.g. '12 GB'")
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--kelt", action="store_true", help="Run KELT inline-barcode contamination QC")
    parser.add_argument("--kelt-barcode-map", default="", help="Source barcode map recorded in run provenance")
    parser.add_argument("--worker-token", default="", help=argparse.SUPPRESS)
    parser.add_argument("--worker-state", default="", help=argparse.SUPPRESS)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="simplseq",
        description=help_description(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("scan", help="Pair FASTQs and write samples.csv")
    p.add_argument("--fastq-dir", default="data")
    p.add_argument("--out", default="samples.csv")
    p.add_argument("--absolute", action="store_true")
    p.add_argument("--include-pool-in-sample-id", action="store_true")
    p.add_argument("--metadata", default="", help="Optional CSV/TSV/XLSX metadata file to enrich samples.csv")
    p.add_argument("--metadata-sheet", default="", help="Optional worksheet name for XLSX metadata")
    p.add_argument(
        "--fallback-collection-year",
        type=collection_year_arg,
        default="",
        help="Year to use when filenames contain a month but no year",
    )
    p.add_argument(
        "--fallback-collection-day",
        type=collection_day_arg,
        default=DEFAULT_COLLECTION_DAY,
        help=f"Day to use when collection dates have year/month but no day [default {DEFAULT_COLLECTION_DAY}]",
    )
    p.set_defaults(func=cmd_scan)

    p = sub.add_parser("check", help="Check local runtime and optional inputs")
    p.add_argument("--samples", default=None)
    p.add_argument("--out", "--outdir", dest="out", default="results")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("run", help="Open the malaria-amplicon-nf browser interface")
    p.add_argument("--port", type=int, default=8501, help="Preferred local browser port")
    p.add_argument("--host", default="127.0.0.1", help=argparse.SUPPRESS)
    p.add_argument("--no-browser", action="store_true", help=argparse.SUPPRESS)
    p.set_defaults(func=cmd_app)

    p = sub.add_parser("run-headless", help="Run the workflow without the browser GUI")
    add_direct_run_args(p)
    p.set_defaults(func=cmd_run_direct)

    p = sub.add_parser("status", help="Show local run status")
    p.add_argument("--out", "--outdir", dest="out", default="results")
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("results", help="List expected output files")
    p.add_argument("--out", "--outdir", dest="out", default="results")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_results)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    return args.func(args)
