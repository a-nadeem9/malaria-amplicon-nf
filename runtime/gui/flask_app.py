"""Flask browser app for malaria-amplicon-nf."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid
import webbrowser
import zipfile
from collections import Counter
from io import BytesIO
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

from flask import Flask, abort, jsonify, render_template, request, send_file, send_from_directory, url_for
from waitress import serve

from simplseq import __version__
from simplseq.analysis_jobs import launch_job, reconcile_state
from simplseq.job_state import process_matches, read_json, terminate_process_group, write_json
from simplseq.pathutils import user_path
from simplseq.progress import read_events
from simplseq.resources import human_bytes
from simplseq.runner import (
    STAGES,
    check_environment,
    local_runtime_env,
    progress_summary,
    project_root,
    results_manifest,
)
from simplseq.metadata import (
    detection_value_options,
    discover_metadata_file,
    enrich_rows_with_metadata,
    infer_metadata_year,
    inspect_metadata,
    normalize_metadata_contract,
)
from simplseq.kelt import inspect_kelt_barcode_map
from simplseq.panel import panel_profile
from simplseq.samplesheet import (
    DEFAULT_COLLECTION_DAY,
    SAMPLE_FIELDS,
    FastqPair,
    FastqScan,
    pair_to_row,
    scan_fastqs,
    write_samples_csv,
)


RUN_PROCESSES: dict[str, subprocess.Popen[str]] = {}
RUN_LOCK = threading.Lock()
DOWNLOAD_SLUG_RE = re.compile(r"[^a-z0-9]+")
RUN_NAME_RE = re.compile(r"[^A-Za-z0-9._-]+")
ANSI_RE = re.compile(r"\x1B\][^\x07]*(?:\x07|\x1B\\)|\x1B\[[0-?]*[ -/]*[@-~]|\x1B[@-Z\\-_]")
ASV_FILTERING_SUMMARY_LABEL = "ASV filtering summary"
CORE_RESULT_LABELS = {
    ASV_FILTERING_SUMMARY_LABEL,
    "ASV count table",
    "Mapped ASV table",
    "ASV to CIGAR map",
    "CIGAR count table",
}
KELT_REPORT_LABEL = "KELT contamination report"
KELT_TABLE_LABELS = {"KELT contamination calls", "KELT barcode counts"}
CDHIT_TABLE_LABELS = {"CD-HIT cluster membership", "CD-HIT cluster count table"}
CDHIT_RESULT_LABELS = CDHIT_TABLE_LABELS | {
    "CD-HIT representative FASTA",
    "CD-HIT filtered ASV FASTA",
    "CD-HIT raw clusters",
    "CD-HIT cluster summary",
}
PREFERENCE_KEYS = {
    "metadataDateOrder", "fallbackCollectionYear", "fallbackCollectionDay",
    "resumeRun", "dryRun", "cpus", "memory",
    "analysisMinAbundancePct", "analysisAbundanceDenominator", "analysisCdhitMode",
    "dinemitesEnabled", "dinemitesModel", "dinemitesNLags", "dinemitesTLag",
    "dinemitesMinAbundancePct", "dinemitesAbundanceDenominator", "dinemitesNoDayCutoff",
    "dinemitesImputations",
    "dinemitesSeed", "dinemitesRefresh", "dinemitesBayesianLagDays", "dinemitesBayesianChains",
    "dinemitesBayesianParallelChains", "dinemitesBayesianWarmup", "dinemitesBayesianSampling",
    "dinemitesBayesianAdaptDelta", "dinemitesUseSeasonCovariate", "dinemitesUseAgeCovariate",
    "dinemitesUseGenderCovariate", "dinemitesCustomCovariates", "dinemitesBayesianDropOut",
    "dciferEnabled", "dciferMinAbundancePct", "dciferAbundanceDenominator", "dciferCoiLrank",
    "dciferIbdGridNr", "dciferAlpha",
}


def analysis_runtime_env(root: Path) -> dict[str, str]:
    """Return the managed analysis environment used by R subprocesses."""
    env = os.environ.copy()
    env["SIMPLSEQ_PROJECT_ROOT"] = str(root)
    env_dir_raw = env.get("SIMPLSEQ_ENV_DIR", "").strip()
    if env_dir_raw:
        env_dir = Path(env_dir_raw)
        env.setdefault("CONDA_PREFIX", str(env_dir))
        cmdstan = env_dir / "bin" / "cmdstan"
        if cmdstan.is_dir():
            env.setdefault("CMDSTAN", str(cmdstan))
        env_bin = str(env_dir / "bin")
        path_parts = env.get("PATH", "").split(os.pathsep)
        if env_bin not in path_parts:
            env["PATH"] = os.pathsep.join([env_bin, *path_parts])
    return env
REPORT_LABEL = "Run summary"
VIEWABLE_REPORT_LABELS = {REPORT_LABEL, KELT_REPORT_LABEL, ASV_FILTERING_SUMMARY_LABEL}
TABLE_RESULT_LABELS = (CORE_RESULT_LABELS - {ASV_FILTERING_SUMMARY_LABEL}) | KELT_TABLE_LABELS | CDHIT_TABLE_LABELS
BUNDLE_RESULT_LABELS = CORE_RESULT_LABELS | KELT_TABLE_LABELS | CDHIT_RESULT_LABELS | {
    REPORT_LABEL,
    ASV_FILTERING_SUMMARY_LABEL,
    KELT_REPORT_LABEL,
    "KELT contamination summary",
    "Input FASTQ MD5s",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def result_panel_profile(outdir: Path) -> dict[str, Any]:
    """Describe recovered loci without assuming every run is SIMPLseq."""
    cigar_path = outdir / "run_dada2" / "seqtab_cigar.tsv"
    loci: set[str] = set()
    if safe_is_file(cigar_path):
        try:
            with cigar_path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
                header = next(csv.reader(handle, delimiter="\t"), [])
            loci = {
                column.split(",", 1)[0].strip()
                for column in header[1:]
                if "," in column and column.split(",", 1)[0].strip()
            }
        except OSError:
            loci = set()
    return panel_profile(loci)


def freeze_run_configuration(
    outdir: Path,
    samples: Path,
    *,
    metadata_path: Path | None,
    metadata_sheet: str,
    metadata_date_order: str,
    fallback_year: str,
    fallback_day: str,
    metadata_contract: dict[str, Any],
    selected_libraries: list[str],
    request_payload: dict[str, Any],
) -> dict[str, Any]:
    """Freeze run inputs so downstream analyses can be recreated later."""
    configuration_dir = outdir / "configuration"
    configuration_dir.mkdir(parents=True, exist_ok=True)
    metadata_profile: dict[str, Any] = {
        "available": False,
        "source_path": "",
        "frozen_path": "",
        "sha256": "",
        "sheet": metadata_sheet,
        "date_order": metadata_date_order,
        "fallback_year": fallback_year,
        "fallback_day": fallback_day,
        "mapping": {},
        "contract": normalize_metadata_contract(metadata_contract),
        "issues": [],
    }
    if metadata_path is not None:
        suffix = metadata_path.suffix.lower() or ".dat"
        frozen_metadata = configuration_dir / f"metadata_source{suffix}"
        shutil.copy2(metadata_path, frozen_metadata)
        normalized_contract = normalize_metadata_contract(metadata_contract)
        catalog = inspect_metadata(
            frozen_metadata,
            metadata_sheet,
            date_order=metadata_date_order,
            column_overrides=normalized_contract["columns"],
        )
        catalog_summary = catalog.summary()
        resolved_detection_map = {
            option["value"]: option["state"]
            for option in detection_value_options(
                catalog.value_counts.get("metadata_pcr", {}),
                normalized_contract["detection_value_map"],
            )
        }
        normalized_contract["columns"] = dict(catalog.columns)
        normalized_contract["detection_value_map"] = resolved_detection_map
        metadata_profile.update({
            "available": True,
            "source_path": str(metadata_path),
            "frozen_path": str(frozen_metadata),
            "sha256": file_sha256(frozen_metadata),
            "sheet": catalog.sheet,
            "header_row": catalog.header_row,
            "mapping": catalog.columns,
            "contract": normalized_contract,
            "records": len(catalog.records),
            "participants": catalog_summary.get("participants", 0),
            "issues": catalog_summary.get("issues", []),
            "issue_count_total": catalog_summary.get("issue_count_total", 0),
        })
    write_json(configuration_dir / "metadata_profile.json", metadata_profile)

    run_manifest = {
        "schema_version": 1,
        "app_version": __version__,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "samples_path": str(samples),
        "samples_sha256": file_sha256(samples) if samples.is_file() else "",
        "metadata": metadata_profile,
        "selected_libraries": selected_libraries or ["all"],
        "runtime": {
            "cpus": request_payload.get("cpus", ""),
            "memory": request_payload.get("memory", ""),
            "resume": payload_resume_enabled(request_payload),
            "dry_run": bool_payload(request_payload, "dry_run", False),
        },
    }
    write_json(configuration_dir / "run_manifest.json", run_manifest)
    return run_manifest


def frozen_metadata_settings(outdir: Path) -> dict[str, Any]:
    profile_path = outdir / "configuration" / "metadata_profile.json"
    profile = read_json(profile_path)
    frozen_path = str(profile.get("frozen_path", "") or "")
    return {
        "path": frozen_path if frozen_path and Path(frozen_path).is_file() else "",
        "sheet": str(profile.get("sheet", "") or ""),
        "date_order": str(profile.get("date_order", "auto") or "auto"),
        "fallback_year": str(profile.get("fallback_year", "") or ""),
        "fallback_day": str(profile.get("fallback_day", "27") or "27"),
        "contract": normalize_metadata_contract(profile.get("contract", {})),
        "profile_path": str(profile_path) if profile_path.is_file() else "",
    }


def parse_asv_filtering_summary_text(content: str) -> dict[str, Any]:
    """Convert the saved human-readable ASV audit into structured UI data."""
    step_pattern = re.compile(r"^- (.+?): ([0-9,]+) ASVs \((.+)\)$")
    requirement_pattern = re.compile(r"^\s+Requirement: (.+)$")
    value_pattern = re.compile(r"^- (.+?): ([0-9,]+)(?:\.|$)")
    steps: list[dict[str, Any]] = []
    sections: dict[str, dict[str, Any]] = {"cigar": {}, "cdhit": {}}
    current_section = ""
    for line in content.splitlines():
        stripped = line.rstrip()
        if stripped == "CIGAR conversion":
            current_section = "cigar"
            continue
        if stripped.startswith("CD-HIT ") and "sensitivity clustering" in stripped:
            current_section = "cdhit"
            continue
        match = step_pattern.match(stripped)
        if match:
            count = int(match.group(2).replace(",", ""))
            previous = steps[-1]["count"] if steps else count
            baseline = steps[0]["count"] if steps else count
            steps.append({
                "step": match.group(1),
                "count": count,
                "removed": max(0, previous - count),
                "retained_previous_pct": round(100.0 * count / previous, 1) if previous else 0.0,
                "retained_start_pct": round(100.0 * count / baseline, 1) if baseline else 0.0,
                "requirement": "",
            })
            current_section = "pipeline"
            continue
        requirement = requirement_pattern.match(stripped)
        if requirement and steps and current_section == "pipeline":
            steps[-1]["requirement"] = requirement.group(1)
            continue
        value = value_pattern.match(stripped)
        if value and current_section in sections:
            key = re.sub(r"[^a-z0-9]+", "_", value.group(1).lower()).strip("_")
            sections[current_section][key] = int(value.group(2).replace(",", ""))
    return {"steps": steps, **sections}


def resolve_app_path(root: Path, value: str | os.PathLike[str] | None, default: str | Path) -> Path:
    raw = str(value if value not in {None, ""} else default).strip()
    path = user_path(raw)
    if not path.is_absolute():
        path = root / path
    return path.expanduser().resolve()


def resolve_run_output_parent(value: str | os.PathLike[str] | None) -> Path:
    raw = str(value or "").strip()
    if not raw:
        raise ValueError("Choose an output folder before starting the run.")
    path = user_path(raw).expanduser()
    if not path.is_absolute():
        raise ValueError("Output folder must be a full absolute path. Use Choose folder to select it.")
    return path.resolve()


def resolve_fastq_folder(value: str | os.PathLike[str] | None) -> Path:
    raw = str(value or "").strip()
    if not raw:
        raise ValueError("Choose a FASTQ folder before scanning.")
    path = user_path(raw).expanduser()
    if not path.is_absolute():
        raise ValueError("FASTQ folder must be a full absolute path. Use Choose folder to select it.")
    return path.resolve()


def preferences_path() -> Path:
    configured = os.environ.get("SIMPLSEQ_CONFIG_DIR", "").strip()
    if configured:
        base = Path(configured).expanduser()
    elif is_windows():
        base = Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming")) / "malaria-amplicon-nf"
    elif is_macos():
        base = Path.home() / "Library" / "Application Support" / "malaria-amplicon-nf"
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")) / "malaria-amplicon-nf"
    return base / "preferences.json"


def sanitized_preferences(value: Any) -> dict[str, Any]:
    if not isinstance(value, dict):
        return {}
    clean: dict[str, Any] = {}
    for key in PREFERENCE_KEYS:
        item = value.get(key)
        if isinstance(item, bool) or item is None:
            clean[key] = item
        elif isinstance(item, (str, int, float)):
            clean[key] = str(item)[:4096] if isinstance(item, str) else item
    return clean


def rel_or_abs(root: Path, path: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def bool_payload(data: dict[str, Any], key: str, default: bool = False) -> bool:
    value = data.get(key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in {"1", "true", "yes", "on"}
    return bool(value)


def int_payload(data: dict[str, Any], key: str, default: int = 0) -> int:
    try:
        return int(data.get(key, default) or default)
    except (TypeError, ValueError):
        return default


def collection_date_defaults(data: dict[str, Any]) -> tuple[str, str]:
    year = str(data.get("fallback_collection_year", "") or "").strip()
    if year and not re.fullmatch(r"(19|20)[0-9]{2}", year):
        raise ValueError("Fallback year must be a four-digit year from 1900 to 2099.")
    raw_day = str(data.get("fallback_collection_day", DEFAULT_COLLECTION_DAY) or DEFAULT_COLLECTION_DAY).strip()
    if not re.fullmatch(r"[0-9]{1,2}", raw_day):
        raise ValueError("Fallback day must be from 1 to 31.")
    day = int(raw_day)
    if day < 1 or day > 31:
        raise ValueError("Fallback day must be from 1 to 31.")
    return year, f"{day:02d}"


def dinemites_model_settings(data: dict[str, Any]) -> dict[str, int | float | str | bool]:
    def require_int(key: str, default: int, minimum: int) -> int:
        try:
            value = int(data.get(key, default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be at least {minimum}.") from exc
        if value < minimum:
            raise ValueError(f"{key} must be at least {minimum}.")
        return value

    def require_float_range(key: str, default: float, minimum: float, maximum: float) -> float:
        try:
            value = float(data.get(key, default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be between {minimum} and {maximum}.") from exc
        if value < minimum or value > maximum:
            raise ValueError(f"{key} must be between {minimum} and {maximum}.")
        return value

    try:
        n_lags = int(data.get("n_lags", 3))
    except (TypeError, ValueError) as exc:
        raise ValueError("n_lags must be at least 1.") from exc
    if n_lags < 1:
        raise ValueError("n_lags must be at least 1.")

    try:
        min_abundance_pct = float(data.get("min_abundance_pct", 1.0))
    except (TypeError, ValueError) as exc:
        raise ValueError("min_abundance_pct must be between 0 and 100.") from exc
    if min_abundance_pct < 0 or min_abundance_pct > 100:
        raise ValueError("min_abundance_pct must be between 0 and 100.")
    abundance_denominator = str(data.get("abundance_denominator", "locus")).strip().lower()
    if abundance_denominator not in {"locus", "sample"}:
        raise ValueError("abundance_denominator must be either locus or sample.")

    seed = require_int("seed", 1, 1)
    refresh = require_int("refresh", 100, 0)
    n_imputations = require_int("n_imputations", 10, 1)
    if n_imputations > 100:
        raise ValueError("n_imputations cannot exceed 100.")
    bayesian_lag_days = require_int("bayesian_lag_days", 30, 1)
    bayesian_chains = require_int("bayesian_chains", 4, 1)
    bayesian_parallel_chains = require_int("bayesian_parallel_chains", 2, 1)
    if bayesian_parallel_chains > bayesian_chains:
        raise ValueError("bayesian_parallel_chains cannot exceed bayesian_chains.")
    bayesian_iter_warmup = require_int("bayesian_iter_warmup", 500, 1)
    bayesian_iter_sampling = require_int("bayesian_iter_sampling", 500, 1)
    bayesian_adapt_delta = require_float_range("bayesian_adapt_delta", 0.99, 0.000001, 0.999999)
    bayesian_drop_out = bool_payload(data, "bayesian_drop_out", False)
    infection_general_covariates = str(data.get("infection_general_covariates", "none") or "none").strip()
    if not infection_general_covariates:
        infection_general_covariates = "none"
    if not re.fullmatch(r"[A-Za-z0-9_, .-]+", infection_general_covariates):
        raise ValueError("infection_general_covariates must be none, auto for season, or comma-separated column names.")

    common_settings: dict[str, int | float | str | bool] = {
        "min_abundance_pct": min_abundance_pct,
        "abundance_denominator": abundance_denominator,
        "seed": seed,
        "refresh": refresh,
        "n_imputations": n_imputations,
        "bayesian_lag_days": bayesian_lag_days,
        "bayesian_chains": bayesian_chains,
        "bayesian_parallel_chains": bayesian_parallel_chains,
        "bayesian_iter_warmup": bayesian_iter_warmup,
        "bayesian_iter_sampling": bayesian_iter_sampling,
        "bayesian_adapt_delta": bayesian_adapt_delta,
        "bayesian_drop_out": bayesian_drop_out,
        "infection_general_covariates": infection_general_covariates,
    }

    no_day_cutoff = bool_payload(data, "no_day_cutoff", False)
    raw_t_lag = str(data.get("t_lag", "Inf")).strip()
    if no_day_cutoff or raw_t_lag.lower() in {"", "inf", "infinity", "none"}:
        return {
            "n_lags": n_lags,
            "t_lag": "Inf",
            **common_settings,
        }

    try:
        numeric_t_lag = float(raw_t_lag)
    except ValueError as exc:
        raise ValueError("t_lag must be a non-negative number or Inf.") from exc
    if numeric_t_lag < 0:
        raise ValueError("t_lag must be a non-negative number or Inf.")

    return {
        "n_lags": n_lags,
        "t_lag": str(int(numeric_t_lag) if numeric_t_lag.is_integer() else numeric_t_lag),
        **common_settings,
    }


def dcifer_settings(data: dict[str, Any]) -> dict[str, int | float | str]:
    def require_int(key: str, default: int, minimum: int) -> int:
        try:
            value = int(data.get(key, default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be at least {minimum}.") from exc
        if value < minimum:
            raise ValueError(f"{key} must be at least {minimum}.")
        return value

    try:
        min_abundance_pct = float(data.get("min_abundance_pct", 1.0))
    except (TypeError, ValueError) as exc:
        raise ValueError("min_abundance_pct must be between 0 and 100.") from exc
    if min_abundance_pct < 0 or min_abundance_pct > 100:
        raise ValueError("min_abundance_pct must be between 0 and 100.")

    abundance_denominator = str(data.get("abundance_denominator", "locus")).strip().lower()
    if abundance_denominator not in {"locus", "sample"}:
        raise ValueError("abundance_denominator must be either locus or sample.")

    try:
        alpha = float(data.get("alpha", 0.05))
    except (TypeError, ValueError) as exc:
        raise ValueError("alpha must be greater than 0 and less than 1.") from exc
    if alpha <= 0 or alpha >= 1:
        raise ValueError("alpha must be greater than 0 and less than 1.")

    afreq_mode = str(data.get("afreq_mode", "current_run")).strip().lower()
    if afreq_mode != "current_run":
        raise ValueError("afreq_mode currently supports only current_run.")

    return {
        "min_abundance_pct": min_abundance_pct,
        "abundance_denominator": abundance_denominator,
        "coi_lrank": require_int("coi_lrank", 2, 1),
        "ibd_grid_nr": require_int("ibd_grid_nr", 1000, 1),
        "alpha": alpha,
        "afreq_mode": afreq_mode,
    }


def slugify(label: str) -> str:
    slug = DOWNLOAD_SLUG_RE.sub("-", label.lower()).strip("-")
    return slug or "file"


def default_run_name(now: dt.datetime | None = None) -> str:
    stamp = (now or dt.datetime.now()).strftime("%Y-%m-%d_%H%M%S")
    return f"SIMPLseq_{stamp}"


def safe_run_name(value: str | None, *, now: dt.datetime | None = None) -> str:
    raw = str(value or "").strip() or default_run_name(now)
    name = RUN_NAME_RE.sub("_", raw).strip("._-")
    return name or default_run_name(now)


def allocate_run_outdir(
    parent: Path,
    run_name: str | None,
    *,
    now: dt.datetime | None = None,
    reuse_existing: bool = False,
) -> Path:
    parent.mkdir(parents=True, exist_ok=True)
    name = safe_run_name(run_name, now=now)
    candidate = parent / name
    if reuse_existing:
        candidate.mkdir(parents=True, exist_ok=True)
        return candidate.resolve()
    for index in range(1, 1000):
        unique = candidate if index == 1 else parent / f"{name}_{index:02d}"
        try:
            unique.mkdir(parents=True, exist_ok=False)
            return unique.resolve()
        except FileExistsError:
            continue
    raise RuntimeError(f"Could not allocate a unique run folder under {parent}")


def json_error(message: str, status: int = 400, **extra: Any):
    payload = {"ok": False, "error": message}
    payload.update(extra)
    return jsonify(payload), status


def log_failure_detail(log_path: Path, fallback: str) -> str:
    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return fallback
    for line in reversed(lines):
        clean = ANSI_RE.sub("", line).strip()
        if not clean:
            continue
        if "Execution halted" in clean:
            continue
        if "ERROR:" in clean:
            return clean.split("ERROR:", 1)[1].strip()
        if clean.startswith("Error"):
            return clean
    return fallback


def sample_pair_json(pair: FastqPair, root: Path) -> dict[str, str]:
    return {
        "sample_id": pair.sample_id,
        "biological_sample_id": pair.biological_sample_id or pair.sample_id,
        "library": pair.library,
        "participant_id": pair.participant_id,
        "collection_date": pair.collection_date,
        "collection_date_inferred": pair.collection_date_inferred,
        "collection_date_source": pair.collection_date_source,
        "inferred_year": pair.inferred_year,
        "inferred_day": pair.inferred_day,
        "date_note": pair.date_note,
        "replicate": pair.replicate,
        "sample_type": pair.sample_type,
        "fastq_1": rel_or_abs(root, pair.fastq_1),
        "fastq_2": rel_or_abs(root, pair.fastq_2),
    }


def scan_json(
    scan: FastqScan,
    root: Path,
    *,
    preview_limit: int = 100,
    fallback_collection_year: str = "",
    fallback_collection_day: str = DEFAULT_COLLECTION_DAY,
) -> dict[str, Any]:
    missing_pairs = len(scan.missing_r2) + len(scan.orphan_r2)
    return {
        "fastq_dir": str(scan.fastq_dir),
        "pair_count": len(scan.pairs),
        "md5_files": scan.md5_files,
        "total_fastq_bytes": scan.total_fastq_bytes,
        "total_fastq_size": human_bytes(scan.total_fastq_bytes),
        "missing_pairs": missing_pairs,
        "missing_r2": scan.missing_r2[:100],
        "orphan_r2": scan.orphan_r2[:100],
        "duplicate_sample_ids": scan.duplicate_sample_ids,
        "auto_disambiguated_sample_ids": scan.auto_disambiguated_sample_ids,
        "collection_month_without_year": sum(1 for pair in scan.pairs if pair.inferred_year),
        "missing_collection_year_count": sum(
            1 for pair in scan.pairs if pair.collection_date_source == "filename_month_missing_year"
        ),
        "inferred_collection_dates": sum(1 for pair in scan.pairs if pair.collection_date_inferred),
        "inferred_year_count": sum(1 for pair in scan.pairs if pair.inferred_year),
        "inferred_day_count": sum(1 for pair in scan.pairs if pair.inferred_day),
        "fallback_collection_year": fallback_collection_year,
        "fallback_collection_day": fallback_collection_day,
        "preview": [sample_pair_json(pair, root) for pair in scan.pairs[:preview_limit]],
    }


def read_samples_preview(path: Path, limit: int = 100) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        rows = []
        for index, row in enumerate(csv.DictReader(handle)):
            if index >= limit:
                break
            rows.append({field: row.get(field, "") for field in SAMPLE_FIELDS})
        return rows


def sample_metadata_match_counts(path: Path) -> dict[str, int]:
    if not path.exists():
        return {}
    counts: dict[str, int] = {}
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle)
            if "metadata_match_status" not in (reader.fieldnames or []):
                return {}
            for row in reader:
                status = str(row.get("metadata_match_status", "") or "unknown").strip() or "unknown"
                counts[status] = counts.get(status, 0) + 1
    except OSError:
        return {}
    return counts


def sample_metadata_match_counts_from_rows(rows: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        status = str(row.get("metadata_match_status", "") or "").strip()
        if not status:
            continue
        counts[status] = counts.get(status, 0) + 1
    return counts


def analysis_readiness(outdir: Path, mode: str = "primary") -> dict[str, Any]:
    """Summarize prerequisites from the exact prepared table selected for analysis."""
    normalized_mode = "cdhit98" if str(mode).strip().lower() in {
        "cdhit98", "cdhit_summed", "summed"
    } else "primary"
    stem = "cdhit98" if normalized_mode == "cdhit98" else "primary"
    prepared_samples = outdir / "analysis_input" / f"{stem}_samples.csv"
    prepared_table = outdir / "analysis_input" / f"{stem}_seqtab.tsv"
    samples_path = prepared_samples if safe_is_file(prepared_samples) else outdir / "samples.csv"
    cigar_path = prepared_table if safe_is_file(prepared_table) else outdir / "run_dada2" / "seqtab_cigar.tsv"

    table_sample_ids: set[str] = set()
    loci: set[str] = set()
    if safe_is_file(cigar_path):
        try:
            with cigar_path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, [])
                loci = {
                    column.split(",", 1)[0].strip()
                    for column in header[1:]
                    if "," in column and column.split(",", 1)[0].strip()
                }
                table_sample_ids = {
                    str(row[0]).strip() for row in reader if row and str(row[0]).strip()
                }
        except OSError:
            table_sample_ids = set()
            loci = set()

    samples: list[dict[str, str]] = []
    if safe_is_file(samples_path):
        try:
            with samples_path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
                for row in csv.DictReader(handle):
                    sample_id = str(row.get("sample_id", "") or "").strip()
                    participant = str(row.get("participant_id", "") or "").strip()
                    sample_type = str(row.get("sample_type", "sample") or "sample").strip().lower()
                    if sample_type != "sample" or not participant:
                        continue
                    if table_sample_ids and sample_id not in table_sample_ids:
                        continue
                    samples.append(row)
        except OSError:
            samples = []

    participant_visits: dict[str, set[str]] = {}
    missing_dates = 0
    for row in samples:
        participant = str(row.get("participant_id", "") or "").strip()
        collection_date = str(row.get("collection_date", "") or "").strip()
        if not collection_date:
            missing_dates += 1
        if participant:
            participant_visits.setdefault(participant, set()).add(collection_date)
    repeated_subjects = sum(1 for visits in participant_visits.values() if len({value for value in visits if value}) > 1)

    env_dir = Path(os.environ.get("SIMPLSEQ_ENV_DIR", "")) if os.environ.get("SIMPLSEQ_ENV_DIR") else None
    cmdstan_available = bool(
        (os.environ.get("CMDSTAN") and Path(os.environ["CMDSTAN"]).is_dir())
        or (env_dir and (env_dir / "bin" / "cmdstan").is_dir())
    )
    return {
        "samples": len(samples),
        "participants": len(participant_visits),
        "repeated_subjects": repeated_subjects,
        "missing_dates": missing_dates,
        "loci": len(loci),
        "cmdstan_available": cmdstan_available,
        "analysis_mode": normalized_mode,
        "dinemites_ready": bool(samples and repeated_subjects and missing_dates == 0 and safe_is_file(cigar_path)),
        "dcifer_ready": bool(len(samples) >= 2 and len(loci) >= 2 and safe_is_file(cigar_path)),
        "dcifer_single_locus": len(loci) == 1,
    }


def sample_rows_preview(
    scan: FastqScan,
    output_root: Path,
    *,
    absolute: bool,
    metadata_path: Path | None,
    metadata_sheet: str,
    metadata_date_order: str = "auto",
    metadata_columns: dict[str, str] | None = None,
    limit: int = 100,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, int], dict[str, Any]]:
    rows = [pair_to_row(pair, output_root, absolute) for pair in scan.pairs]
    rows = enrich_rows_with_metadata(
        rows,
        metadata_path,
        metadata_sheet=metadata_sheet,
        date_order=metadata_date_order,
        column_overrides=metadata_columns,
    )
    preview_fields = [*SAMPLE_FIELDS, "metadata_match_status", "metadata_status"]
    preview = [{field: row.get(field, "") for field in preview_fields} for row in rows[:limit]]

    def review_status(row: dict[str, str]) -> str:
        sample_type = str(row.get("sample_type", "sample")).strip().lower()
        if sample_type in {"negative", "positive", "control"}:
            return "excluded"
        inferred = any(str(row.get(field, "")).lower() == "true" for field in ("inferred_year", "inferred_day"))
        ambiguous = str(row.get("metadata_match_status", "")).lower() == "ambiguous"
        missing_date = not str(row.get("collection_date", "")).strip()
        return "review" if inferred or ambiguous or missing_date else "ready"

    def review_reason(row: dict[str, str]) -> str:
        status = review_status(row)
        sample_type = str(row.get("sample_type", "sample")).strip().lower()
        if status == "excluded":
            return f"{sample_type.title()} control" if sample_type in {"negative", "positive"} else "Control"
        if not str(row.get("collection_date", "")).strip():
            return "Collection date missing"
        if str(row.get("metadata_match_status", "")).lower() == "ambiguous":
            return "Multiple metadata rows"
        inferred_year = str(row.get("inferred_year", "")).lower() == "true"
        inferred_day = str(row.get("inferred_day", "")).lower() == "true"
        if inferred_year and inferred_day:
            return "Year + day filled automatically"
        if inferred_year:
            return "Year filled automatically"
        if inferred_day:
            return "Day filled automatically"
        return ""

    scan_fields = [
        "sample_id",
        "biological_sample_id",
        "participant_id",
        "library",
        "collection_date",
        "replicate",
        "sample_type",
        "inferred_year",
        "inferred_day",
        "metadata_match_status",
    ]
    scan_rows = []
    for row in rows:
        concise = {field: row.get(field, "") for field in scan_fields}
        concise["review_status"] = review_status(row)
        concise["review_reason"] = review_reason(row)
        scan_rows.append(concise)

    sample_status_counts = {
        status: sum(review_status(row) == status for row in rows)
        for status in ("ready", "review", "excluded")
    }
    date_counts = {
        "unresolved_collection_date_count": sum(not str(row.get("collection_date", "")).strip() for row in rows),
        "inferred_year_count": sum(str(row.get("inferred_year", "")).lower() == "true" for row in rows),
        "inferred_day_count": sum(str(row.get("inferred_day", "")).lower() == "true" for row in rows),
        "participant_count": len({
            str(row.get("participant_id", "")).strip()
            for row in rows
            if str(row.get("participant_id", "")).strip()
            and str(row.get("sample_type", "sample")).strip().lower() == "sample"
        }),
        "sample_status_counts": sample_status_counts,
        "library_counts": dict(sorted(Counter(
            str(row.get("library", "")).strip()
            for row in rows
            if str(row.get("library", "")).strip()
        ).items())),
    }
    date_counts["inferred_collection_dates"] = sum(
        str(row.get("inferred_year", "")).lower() == "true"
        or str(row.get("inferred_day", "")).lower() == "true"
        for row in rows
    )
    return preview, scan_rows, sample_metadata_match_counts_from_rows(rows), date_counts


def file_tail(path: Path, max_bytes: int) -> tuple[str, bool]:
    if not path.exists():
        return "", False
    size = path.stat().st_size
    with path.open("rb") as handle:
        if size > max_bytes:
            handle.seek(size - max_bytes)
            data = handle.read()
            return data.decode("utf-8", errors="replace"), True
        return handle.read().decode("utf-8", errors="replace"), False


def clean_log_text(text: str) -> str:
    cleaned = ANSI_RE.sub("", text)
    return "".join(ch for ch in cleaned if ch in {"\n", "\r", "\t"} or ord(ch) >= 32).replace("\r", "\n")


def run_log_tail(outdir: Path, max_bytes: int) -> tuple[str, bool]:
    text, truncated = file_tail(outdir / "technical_log.txt", max_bytes)
    if text.strip():
        return text, truncated

    parts = []
    for path in (outdir / "logs" / "flask-run.stdout.log", outdir / "logs" / "flask-run.stderr.log"):
        log_text, log_truncated = file_tail(path, max_bytes)
        if log_text.strip():
            parts.append(log_text.rstrip())
        truncated = truncated or log_truncated
    return "\n".join(parts), truncated


def safe_exists(path: Path) -> bool:
    try:
        return path.exists()
    except OSError:
        return False


def safe_is_dir(path: Path) -> bool:
    try:
        return path.is_dir()
    except OSError:
        return False


def safe_is_file(path: Path) -> bool:
    try:
        return path.is_file()
    except OSError:
        return False


def safe_resolve(path: Path) -> Path | None:
    try:
        return path.expanduser().resolve()
    except (OSError, RuntimeError):
        return None


def nearest_existing_directory(path: Path, fallback: Path | None = None) -> Path:
    """Return the closest existing directory without leaving users at a dead path."""
    candidate = safe_resolve(path) or path.expanduser()
    while True:
        if safe_is_dir(candidate):
            return candidate
        parent = candidate.parent
        if parent == candidate:
            break
        candidate = parent

    fallback_path = safe_resolve(fallback) if fallback is not None else None
    if fallback_path is not None and safe_is_dir(fallback_path):
        return fallback_path
    return Path.cwd().resolve()


def is_wsl() -> bool:
    try:
        proc_version = Path("/proc/version").read_text(encoding="utf-8", errors="ignore").lower()
    except OSError:
        return False
    return "microsoft" in proc_version or "wsl" in proc_version


def path_style() -> str:
    if os.name == "nt":
        return "windows"
    if is_wsl():
        return "wsl"
    return "posix"


def is_macos() -> bool:
    return sys.platform == "darwin"


def is_windows() -> bool:
    return os.name == "nt"


def wsl_to_windows_path(path: Path) -> str:
    text = str(path)
    match = re.match(r"^/mnt/([a-zA-Z])(?:/(.*))?$", text)
    if match:
        drive = match.group(1).upper()
        rest = (match.group(2) or "").replace("/", "\\")
        return f"{drive}:\\{rest}" if rest else f"{drive}:\\"
    try:
        completed = subprocess.run(
            ["wslpath", "-w", text],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=3,
        )
    except (OSError, subprocess.TimeoutExpired):
        return text
    return completed.stdout.strip() or text


def windows_to_wsl_path(value: str) -> str:
    text = value.strip().replace("\\", "/")
    match = re.match(r"^([A-Za-z]):/?(.*)$", text)
    if not match:
        return value.strip()
    drive = match.group(1).lower()
    rest = match.group(2).strip("/")
    return f"/mnt/{drive}/{rest}" if rest else f"/mnt/{drive}"


def is_windows_network_path(value: str) -> bool:
    text = value.strip().replace("/", "\\")
    return text.startswith("\\\\")


def unsupported_network_path_result() -> dict[str, Any]:
    return {
        "ok": False,
        "selected": False,
        "error": (
            "Windows network-share folders are not directly available inside WSL. "
            "Copy the data to a local Windows drive, or mount the share inside WSL "
            "and paste its Linux mount path manually."
        ),
    }


def select_tkinter_folder_dialog(
    initial: Path | None = None,
    *,
    prompt: str = "Select a folder",
    allow_new_folder: bool = False,
) -> dict[str, Any]:
    """Open a foreground native picker for a Windows-hosted Python runtime."""
    try:
        import tkinter  # noqa: F401
    except (ImportError, RuntimeError) as exc:
        return {"ok": False, "error": f"Windows folder picker is unavailable: {exc}"}

    start_dir = nearest_existing_directory(initial or Path.cwd(), Path.cwd())
    # Flask handles requests on worker threads. Tk is much more reliable when its
    # event loop owns the main thread, so the chooser runs in a short-lived child.
    script = r"""
import sys
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
root.update_idletasks()
selected = filedialog.askdirectory(
    parent=root,
    title=sys.argv[2],
    initialdir=sys.argv[1],
    mustexist=sys.argv[3] == "1",
)
print(selected or "")
root.destroy()
"""
    try:
        completed = subprocess.run(
            [sys.executable, "-c", script, str(start_dir), prompt, "0" if allow_new_folder else "1"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=115,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "error": f"Windows folder picker could not open: {exc}"}
    selected = completed.stdout.strip().splitlines()[-1].strip() if completed.stdout.strip() else ""
    if completed.returncode != 0 and not selected:
        detail = " ".join(completed.stderr.split())[:220] if completed.stderr else "Tk folder picker failed"
        return {"ok": False, "error": detail}

    if not selected:
        return {"ok": True, "selected": False}
    return {"ok": True, "selected": True, "path": str(Path(selected).resolve())}


def select_windows_folder_dialog(
    initial: Path | None = None,
    *,
    prompt: str = "Select a folder",
    allow_new_folder: bool = False,
    convert_to_wsl: bool = True,
) -> dict[str, Any]:
    env = os.environ.copy()
    if initial:
        env["SIMPLSEQ_PICKER_INITIAL"] = wsl_to_windows_path(initial) if convert_to_wsl else str(initial)
    env["SIMPLSEQ_PICKER_PROMPT"] = prompt
    env["SIMPLSEQ_PICKER_ALLOW_CREATE"] = "1" if allow_new_folder else "0"
    script = r"""
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;

[ComImport, Guid("DC1C5A9C-E88A-4DDE-A5A1-60F82A20AEF7")]
internal class FileOpenDialogCom { }

[ComImport, Guid("43826D1E-E718-42EE-BC55-A1E261C37BFE"), InterfaceType(ComInterfaceType.InterfaceIsIUnknown)]
internal interface IShellItem {
    void BindToHandler(IntPtr pbc, ref Guid bhid, ref Guid riid, out IntPtr ppv);
    void GetParent(out IShellItem parent);
    void GetDisplayName(uint sigdnName, out IntPtr name);
    void GetAttributes(uint mask, out uint attributes);
    void Compare(IShellItem item, uint hint, out int order);
}

[ComImport, Guid("42F85136-DB7E-439C-85F1-E4075D135FC8"), InterfaceType(ComInterfaceType.InterfaceIsIUnknown)]
internal interface IFileDialog {
    [PreserveSig] int Show(IntPtr parent);
    void SetFileTypes(uint count, IntPtr filterSpec);
    void SetFileTypeIndex(uint index);
    void GetFileTypeIndex(out uint index);
    void Advise(IntPtr events, out uint cookie);
    void Unadvise(uint cookie);
    void SetOptions(uint options);
    void GetOptions(out uint options);
    void SetDefaultFolder(IShellItem folder);
    void SetFolder(IShellItem folder);
    void GetFolder(out IShellItem folder);
    void GetCurrentSelection(out IShellItem item);
    void SetFileName([MarshalAs(UnmanagedType.LPWStr)] string name);
    void GetFileName([MarshalAs(UnmanagedType.LPWStr)] out string name);
    void SetTitle([MarshalAs(UnmanagedType.LPWStr)] string title);
    void SetOkButtonLabel([MarshalAs(UnmanagedType.LPWStr)] string label);
    void SetFileNameLabel([MarshalAs(UnmanagedType.LPWStr)] string label);
    void GetResult(out IShellItem item);
    void AddPlace(IShellItem item, int alignment);
    void SetDefaultExtension([MarshalAs(UnmanagedType.LPWStr)] string extension);
    void Close(int result);
    void SetClientGuid(ref Guid guid);
    void ClearClientData();
    void SetFilter(IntPtr filter);
}

public static class SimplseqFolderPicker {
    private const uint PickFolders = 0x20;
    private const uint ForceFileSystem = 0x40;
    private const uint PathMustExist = 0x800;
    private const uint FileSystemPath = 0x80058000;
    private const int Cancelled = unchecked((int)0x800704C7);

    [DllImport("shell32.dll", CharSet = CharSet.Unicode, PreserveSig = false)]
    private static extern void SHCreateItemFromParsingName(
        [MarshalAs(UnmanagedType.LPWStr)] string path,
        IntPtr bindContext,
        ref Guid interfaceId,
        [MarshalAs(UnmanagedType.Interface)] out IShellItem item);

    [DllImport("user32.dll")]
    private static extern IntPtr GetForegroundWindow();

    public static string Pick(string initialPath, string title) {
        IFileDialog dialog = (IFileDialog)new FileOpenDialogCom();
        try {
            uint options;
            dialog.GetOptions(out options);
            dialog.SetOptions(options | PickFolders | ForceFileSystem | PathMustExist);
            dialog.SetTitle(String.IsNullOrWhiteSpace(title) ? "Select a folder" : title);
            dialog.SetOkButtonLabel("Select folder");
            if (!String.IsNullOrWhiteSpace(initialPath) && System.IO.Directory.Exists(initialPath)) {
                Guid shellItemId = new Guid("43826D1E-E718-42EE-BC55-A1E261C37BFE");
                IShellItem initialItem;
                SHCreateItemFromParsingName(initialPath, IntPtr.Zero, ref shellItemId, out initialItem);
                dialog.SetFolder(initialItem);
                Marshal.ReleaseComObject(initialItem);
            }
            int result = dialog.Show(GetForegroundWindow());
            if (result == Cancelled) return null;
            Marshal.ThrowExceptionForHR(result);
            IShellItem selectedItem;
            dialog.GetResult(out selectedItem);
            try {
                IntPtr pathPointer;
                selectedItem.GetDisplayName(FileSystemPath, out pathPointer);
                try { return Marshal.PtrToStringUni(pathPointer); }
                finally { Marshal.FreeCoTaskMem(pathPointer); }
            } finally {
                Marshal.ReleaseComObject(selectedItem);
            }
        } finally {
            Marshal.ReleaseComObject(dialog);
        }
    }
}
'@

$title = if ($env:SIMPLSEQ_PICKER_PROMPT) { $env:SIMPLSEQ_PICKER_PROMPT } else { "Select a folder" }
$selected = [SimplseqFolderPicker]::Pick($env:SIMPLSEQ_PICKER_INITIAL, $title)
if ($selected) { Write-Output $selected }
"""
    try:
        completed = subprocess.run(
            ["powershell.exe", "-NoLogo", "-NoProfile", "-Sta", "-ExecutionPolicy", "Bypass", "-Command", script],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=25,
            env=env,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "error": f"Folder picker could not open: {exc}"}
    selected = completed.stdout.strip().splitlines()[-1].strip() if completed.stdout.strip() else ""
    if completed.returncode != 0 and not selected:
        detail = " ".join(completed.stderr.split())[:220] if completed.stderr else "PowerShell folder picker failed"
        return {"ok": False, "error": detail}
    if not selected:
        return {"ok": True, "selected": False}
    if convert_to_wsl and is_windows_network_path(selected):
        return unsupported_network_path_result()
    path = windows_to_wsl_path(selected) if convert_to_wsl else selected
    return {"ok": True, "selected": True, "path": path, "windows_path": selected}


def select_windows_file_dialog(
    initial: Path | None = None,
    *,
    prompt: str = "Select metadata file",
    kind: str = "metadata",
    convert_to_wsl: bool = True,
) -> dict[str, Any]:
    env = os.environ.copy()
    if initial:
        env["SIMPLSEQ_PICKER_INITIAL"] = wsl_to_windows_path(initial) if convert_to_wsl else str(initial)
    env["SIMPLSEQ_PICKER_PROMPT"] = prompt
    env["SIMPLSEQ_PICKER_FILTER"] = (
        "KELT barcode maps (*.csv;*.tsv)|*.csv;*.tsv|All files (*.*)|*.*"
        if kind == "kelt"
        else "Metadata files (*.csv;*.tsv;*.xlsx;*.xlsm)|*.csv;*.tsv;*.xlsx;*.xlsm|All files (*.*)|*.*"
    )
    script = r"""
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
[System.Windows.Forms.Application]::EnableVisualStyles()

$owner = New-Object System.Windows.Forms.Form
$owner.Text = "malaria-amplicon-nf file picker"
$owner.Size = New-Object System.Drawing.Size(1, 1)
$owner.StartPosition = "CenterScreen"
$owner.ShowInTaskbar = $false
$owner.TopMost = $true
$owner.Opacity = 0
$owner.Show()
$owner.Activate()

$dialog = New-Object System.Windows.Forms.OpenFileDialog
$dialog.Title = if ($env:SIMPLSEQ_PICKER_PROMPT) { $env:SIMPLSEQ_PICKER_PROMPT } else { "Select metadata file" }
$dialog.Filter = $env:SIMPLSEQ_PICKER_FILTER
$dialog.CheckFileExists = $true
$dialog.Multiselect = $false
if ($env:SIMPLSEQ_PICKER_INITIAL) {
  if (Test-Path -LiteralPath $env:SIMPLSEQ_PICKER_INITIAL -PathType Leaf) {
    $dialog.InitialDirectory = Split-Path -Parent $env:SIMPLSEQ_PICKER_INITIAL
    $dialog.FileName = Split-Path -Leaf $env:SIMPLSEQ_PICKER_INITIAL
  } elseif (Test-Path -LiteralPath $env:SIMPLSEQ_PICKER_INITIAL -PathType Container) {
    $dialog.InitialDirectory = $env:SIMPLSEQ_PICKER_INITIAL
  }
}
$result = $dialog.ShowDialog($owner)
$owner.Dispose()
if ($result -eq [System.Windows.Forms.DialogResult]::OK) {
  Write-Output $dialog.FileName
}
"""
    try:
        completed = subprocess.run(
            ["powershell.exe", "-NoLogo", "-NoProfile", "-Sta", "-ExecutionPolicy", "Bypass", "-Command", script],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=25,
            env=env,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "error": f"File picker could not open: {exc}"}
    selected = completed.stdout.strip().splitlines()[-1].strip() if completed.stdout.strip() else ""
    if completed.returncode != 0 and not selected:
        detail = " ".join(completed.stderr.split())[:220] if completed.stderr else "PowerShell file picker failed"
        return {"ok": False, "error": detail}
    if not selected:
        return {"ok": True, "selected": False}
    if convert_to_wsl and is_windows_network_path(selected):
        return unsupported_network_path_result()
    path = windows_to_wsl_path(selected) if convert_to_wsl else selected
    return {"ok": True, "selected": True, "path": path, "windows_path": selected}


def select_macos_folder_dialog(initial: Path | None = None, *, prompt: str = "Select a folder") -> dict[str, Any]:
    initial_path = initial.as_posix() if initial and safe_exists(initial) else ""
    script = r"""
on run argv
  set promptText to item 1 of argv
  if (count of argv) > 1 and item 2 of argv is not "" then
    try
      set initialAlias to POSIX file (item 2 of argv) as alias
      set chosenFolder to choose folder with prompt promptText default location initialAlias
    on error
      set chosenFolder to choose folder with prompt promptText
    end try
  else
    set chosenFolder to choose folder with prompt promptText
  end if
  return POSIX path of chosenFolder
end run
"""
    try:
        completed = subprocess.run(
            ["osascript", "-e", script, "--", prompt, initial_path],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=25,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "error": f"Folder picker could not open: {exc}"}
    selected = completed.stdout.strip().splitlines()[-1].strip() if completed.stdout.strip() else ""
    if completed.returncode != 0 and not selected:
        detail = " ".join(completed.stderr.split())
        if "User canceled" in detail or "(-128)" in detail:
            return {"ok": True, "selected": False}
        return {"ok": False, "error": detail[:220] or "macOS folder picker failed"}
    if not selected:
        return {"ok": True, "selected": False}
    path = selected.rstrip("/") or "/"
    return {"ok": True, "selected": True, "path": path}


def select_macos_file_dialog(initial: Path | None = None, *, prompt: str = "Select metadata file") -> dict[str, Any]:
    initial_path = initial.as_posix() if initial and safe_exists(initial) else ""
    script = r"""
on run argv
  set promptText to item 1 of argv
  if (count of argv) > 1 and item 2 of argv is not "" then
    try
      set initialAlias to POSIX file (item 2 of argv) as alias
      set chosenFile to choose file with prompt promptText default location initialAlias
    on error
      set chosenFile to choose file with prompt promptText
    end try
  else
    set chosenFile to choose file with prompt promptText
  end if
  return POSIX path of chosenFile
end run
"""
    try:
        completed = subprocess.run(
            ["osascript", "-e", script, "--", prompt, initial_path],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=25,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "error": f"File picker could not open: {exc}"}
    selected = completed.stdout.strip().splitlines()[-1].strip() if completed.stdout.strip() else ""
    if completed.returncode != 0 and not selected:
        detail = " ".join(completed.stderr.split())
        if "User canceled" in detail or "(-128)" in detail:
            return {"ok": True, "selected": False}
        return {"ok": False, "error": detail[:220] or "macOS file picker failed"}
    if not selected:
        return {"ok": True, "selected": False}
    return {"ok": True, "selected": True, "path": selected.rstrip("/") or "/"}


def select_linux_folder_dialog(initial: Path | None = None, *, prompt: str = "Select a folder") -> dict[str, Any]:
    initial_path = str(initial) if initial and safe_exists(initial) else str(Path.home())
    if shutil.which("zenity"):
        command = ["zenity", "--file-selection", "--directory", f"--title={prompt}", f"--filename={initial_path.rstrip('/')}/"]
    elif shutil.which("kdialog"):
        command = ["kdialog", "--getexistingdirectory", initial_path, "--title", prompt]
    else:
        return {
            "ok": False,
            "error": "No native Linux folder chooser is available. Install zenity or kdialog, then try again.",
        }
    try:
        completed = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=115,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "error": f"Folder picker could not open: {exc}"}
    selected = completed.stdout.strip().splitlines()[-1].strip() if completed.stdout.strip() else ""
    if completed.returncode in {1, 130} and not selected:
        return {"ok": True, "selected": False}
    if completed.returncode != 0 and not selected:
        detail = " ".join(completed.stderr.split())[:220]
        return {"ok": False, "error": detail or "Linux folder picker failed"}
    return {"ok": True, "selected": bool(selected), "path": selected.rstrip("/") or "/"}


def select_linux_file_dialog(
    initial: Path | None = None,
    *,
    prompt: str = "Select metadata file",
    kind: str = "metadata",
) -> dict[str, Any]:
    initial_path = str(initial) if initial and safe_exists(initial) else str(Path.home())
    filter_label = "KELT barcode maps" if kind == "kelt" else "Metadata files"
    filter_glob = "*.csv *.tsv" if kind == "kelt" else "*.csv *.tsv *.xlsx *.xlsm"
    if shutil.which("zenity"):
        command = [
            "zenity", "--file-selection", f"--title={prompt}", f"--filename={initial_path}",
            f"--file-filter={filter_label} | {filter_glob}",
        ]
    elif shutil.which("kdialog"):
        command = [
            "kdialog", "--getopenfilename", initial_path,
            f"{filter_label} ({filter_glob})", "--title", prompt,
        ]
    else:
        return {
            "ok": False,
            "error": "No native Linux file chooser is available. Install zenity or kdialog, then try again.",
        }
    try:
        completed = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=115,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"ok": False, "error": f"File picker could not open: {exc}"}
    selected = completed.stdout.strip().splitlines()[-1].strip() if completed.stdout.strip() else ""
    if completed.returncode in {1, 130} and not selected:
        return {"ok": True, "selected": False}
    if completed.returncode != 0 and not selected:
        detail = " ".join(completed.stderr.split())[:220]
        return {"ok": False, "error": detail or "Linux file picker failed"}
    return {"ok": True, "selected": bool(selected), "path": selected}


def select_desktop_picker(
    picker_type: str,
    *,
    initial: Path | None = None,
    prompt: str,
    allow_new_folder: bool = False,
    kind: str = "metadata",
) -> dict[str, Any] | None:
    bridge_raw = os.environ.get("SIMPLSEQ_PICKER_BRIDGE_DIR", "").strip()
    if not bridge_raw:
        return None
    bridge = user_path(bridge_raw)
    if not bridge.is_dir():
        return {"ok": False, "error": "Desktop picker bridge is unavailable. Reopen the desktop app and try again."}

    request_id = uuid.uuid4().hex
    request_path = bridge / f"{request_id}.request.json"
    response_path = bridge / f"{request_id}.response.json"
    temporary_path = bridge / f"{request_id}.request.json.tmp"
    initial_path = ""
    if initial:
        initial_path = wsl_to_windows_path(initial) if is_wsl() else str(initial)
    payload = {
        "picker_type": picker_type,
        "initial": initial_path,
        "prompt": prompt,
        "allow_new_folder": allow_new_folder,
        "kind": kind,
    }
    try:
        temporary_path.write_text(json.dumps(payload), encoding="utf-8")
        os.replace(temporary_path, request_path)
    except OSError as exc:
        return {"ok": False, "error": f"Desktop picker request could not be created: {exc}"}

    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        if response_path.exists():
            response = read_json(response_path)
            response_path.unlink(missing_ok=True)
            request_path.unlink(missing_ok=True)
            if not isinstance(response, dict):
                return {"ok": False, "error": "Desktop picker returned an invalid response."}
            selected_path = str(response.get("path", "") or "").strip()
            if selected_path and is_wsl():
                if is_windows_network_path(selected_path):
                    return unsupported_network_path_result()
                selected_path = windows_to_wsl_path(selected_path)
            return {
                "ok": bool(response.get("ok", False)),
                "selected": bool(response.get("selected", False)),
                "path": selected_path,
                "error": str(response.get("error", "") or ""),
            }
        time.sleep(0.1)

    request_path.unlink(missing_ok=True)
    temporary_path.unlink(missing_ok=True)
    return {"ok": False, "error": "Desktop picker timed out. Close any hidden dialog and try again."}


def select_folder_dialog(
    initial: Path | None = None,
    *,
    prompt: str = "Select a folder",
    allow_new_folder: bool = False,
) -> dict[str, Any]:
    bridged = select_desktop_picker(
        "folder",
        initial=initial,
        prompt=prompt,
        allow_new_folder=allow_new_folder,
    )
    if bridged is not None:
        return bridged
    if is_windows():
        native_result = select_tkinter_folder_dialog(
            initial,
            prompt=prompt,
            allow_new_folder=allow_new_folder,
        )
        if native_result.get("ok"):
            return native_result
        return select_windows_folder_dialog(initial, prompt=prompt, allow_new_folder=allow_new_folder, convert_to_wsl=False)
    if is_macos():
        return select_macos_folder_dialog(initial, prompt=prompt)
    if is_wsl():
        return select_windows_folder_dialog(
            initial,
            prompt=prompt,
            allow_new_folder=allow_new_folder,
            convert_to_wsl=True,
        )
    return select_linux_folder_dialog(initial, prompt=prompt)


def select_file_dialog(
    initial: Path | None = None,
    *,
    prompt: str = "Select metadata file",
    kind: str = "metadata",
) -> dict[str, Any]:
    bridged = select_desktop_picker("file", initial=initial, prompt=prompt, kind=kind)
    if bridged is not None:
        return bridged
    if is_windows():
        return select_windows_file_dialog(initial, prompt=prompt, kind=kind, convert_to_wsl=False)
    if is_macos():
        return select_macos_file_dialog(initial, prompt=prompt)
    if is_wsl():
        return select_windows_file_dialog(initial, prompt=prompt, kind=kind, convert_to_wsl=True)
    return select_linux_file_dialog(initial, prompt=prompt, kind=kind)


def common_paths(workspace_root: Path, app_root: Path) -> list[dict[str, str]]:
    paths: list[tuple[str, Path]] = [
        ("Current folder", workspace_root),
        ("Data in current folder", workspace_root / "data"),
        ("Home", Path.home()),
    ]
    if safe_exists(workspace_root / "test-data"):
        paths.append(("Test data", workspace_root / "test-data"))
    desktop = Path.home() / "Desktop"
    if safe_exists(desktop):
        paths.append(("Desktop", desktop))
    windows_home = Path("/mnt/c/Users") / Path.home().name
    for label, path in [
        ("Windows home", windows_home),
        ("Windows Desktop", windows_home / "Desktop"),
        ("Windows Downloads", windows_home / "Downloads"),
        ("Windows Documents", windows_home / "Documents"),
    ]:
        if safe_exists(path):
            paths.append((label, path))
    for mount in [Path("/mnt/c"), Path("/mnt/d")]:
        if safe_exists(mount):
            paths.append((str(mount), mount))
    windows_users = Path("/mnt/c/Users")
    if safe_exists(windows_users):
        try:
            candidates = sorted(windows_users.iterdir(), key=lambda item: item.name.lower())
        except OSError:
            candidates = []
        for candidate in candidates:
            desktop_dir = candidate / "Desktop"
            if safe_exists(desktop_dir):
                paths.append((f"{candidate.name} Desktop", desktop_dir))
    seen: set[str] = set()
    result = []
    for label, path in paths:
        resolved_path = safe_resolve(path)
        if resolved_path is None:
            continue
        resolved = str(resolved_path)
        if resolved in seen:
            continue
        seen.add(resolved)
        result.append({"label": label, "path": resolved})
    return result


def fastq_count(path: Path) -> int:
    count = 0
    try:
        for item in path.iterdir():
            if item.name.endswith((".fastq.gz", ".fq.gz")) and safe_is_file(item):
                count += 1
    except OSError:
        return 0
    return count


def browse_payload(path: Path) -> dict[str, Any]:
    if not safe_exists(path):
        return {"path": str(path), "exists": False, "directories": [], "parent": str(path.parent)}
    if not safe_is_dir(path):
        return {"path": str(path), "exists": True, "is_dir": False, "directories": [], "parent": str(path.parent)}
    directories = []
    try:
        children = sorted((item for item in path.iterdir() if safe_is_dir(item)), key=lambda item: item.name.lower())
    except OSError:
        children = []
    for child in children[:250]:
        directories.append(
            {
                "name": child.name,
                "path": str(child),
                "fastq_files": fastq_count(child),
            }
        )
    current_fastq_files = fastq_count(path)
    scan = scan_fastqs(path) if current_fastq_files else FastqScan(path, [], [], [], 0, 0, [])
    return {
        "path": str(path),
        "exists": True,
        "is_dir": True,
        "parent": str(path.parent),
        "directories": directories,
        "fastq_files": current_fastq_files,
        "pair_count": len(scan.pairs),
        "missing_pairs": len(scan.missing_r2) + len(scan.orphan_r2),
    }


def result_files_with_downloads(root: Path, outdir: Path) -> dict[str, Any]:
    manifest = results_manifest(outdir)
    files = []
    for item in manifest["files"]:
        label = str(item["label"])
        slug = slugify(label)
        exists = bool(item["exists"])
        size_bytes = int(item["size_bytes"])
        files.append(
            {
                "label": label,
                "slug": slug,
                "path": item["path"],
                "relative_path": rel_or_abs(root, Path(str(item["path"]))),
                "exists": exists,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "status": "ready" if exists and size_bytes else "missing",
                "download_url": url_for("download_result", file_key=slug, out=str(outdir)) if exists else "",
                "view_url": url_for("download_result", file_key=slug, out=str(outdir), inline=1)
                if exists and label in VIEWABLE_REPORT_LABELS
                else "",
                "table_url": url_for("api_result_table", file_key=slug, out=str(outdir))
                if exists and label in TABLE_RESULT_LABELS
                else "",
            }
        )
    manifest["files"] = files
    manifest["report"] = next((item for item in files if item["label"] == REPORT_LABEL), None)
    manifest["core_files"] = [item for item in files if item["label"] in CORE_RESULT_LABELS]
    manifest["kelt_report"] = next((item for item in files if item["label"] == KELT_REPORT_LABEL), None)
    manifest["kelt_files"] = [item for item in files if item["label"] in KELT_TABLE_LABELS]
    manifest["cdhit_files"] = [item for item in files if item["label"] in CDHIT_TABLE_LABELS]
    manifest["cdhit_summary"] = read_json(outdir / "cdhit" / "cdhit_summary.json")
    manifest["panel"] = result_panel_profile(outdir)
    manifest["support_files"] = [
        item for item in files if item["label"] != REPORT_LABEL and item["label"] not in CORE_RESULT_LABELS
    ]
    manifest["ready_counts"] = {
        "core": sum(1 for item in manifest["core_files"] if item["status"] == "ready"),
        "support": sum(1 for item in manifest["support_files"] if item["status"] == "ready"),
    }
    manifest["bundle_ready"] = any(
        item["label"] in BUNDLE_RESULT_LABELS and item["status"] == "ready" for item in files
    )
    manifest["bundle_url"] = url_for("download_bundle", out=str(outdir)) if manifest["bundle_ready"] else ""
    manifest["outdir"] = str(outdir)
    manifest["run_name"] = outdir.name or str(outdir)
    return manifest


def bundle_result_paths(outdir: Path) -> list[tuple[Path, str]]:
    outdir = outdir.resolve()
    bundled: list[tuple[Path, str]] = []
    for item in results_manifest(outdir)["files"]:
        label = str(item["label"])
        if label not in BUNDLE_RESULT_LABELS:
            continue
        path = Path(str(item["path"])).resolve()
        try:
            arcname = path.relative_to(outdir).as_posix()
        except ValueError:
            continue
        if path.exists() and path.is_file() and path.stat().st_size > 0:
            bundled.append((path, arcname))
    return bundled


def downloads_dir() -> Path:
    path = Path.home() / "Downloads"
    path.mkdir(parents=True, exist_ok=True)
    return path


def unique_download_path(filename: str) -> Path:
    safe_name = Path(filename).name or "download"
    target = downloads_dir() / safe_name
    if not target.exists():
        return target
    stem = target.stem
    suffix = target.suffix
    for index in range(1, 1000):
        candidate = target.with_name(f"{stem}_{index:02d}{suffix}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"Could not create a unique download filename for {safe_name}")


def result_file_path(outdir: Path, file_key: str) -> tuple[Path, str]:
    outdir = outdir.resolve()
    manifest = results_manifest(outdir)
    for item in manifest["files"]:
        label = str(item["label"])
        if slugify(label) != file_key:
            continue
        path = Path(str(item["path"])).resolve()
        try:
            path.relative_to(outdir)
        except ValueError:
            abort(404)
        if not path.exists() or not path.is_file():
            abort(404)
        return path, label
    abort(404)


def open_directory(path: Path) -> None:
    if is_windows():
        os.startfile(str(path))  # type: ignore[attr-defined]
        return
    if is_macos():
        subprocess.Popen(["open", str(path)])
        return
    if is_wsl():
        subprocess.Popen(["explorer.exe", wsl_to_windows_path(path)])
        return
    subprocess.Popen(["xdg-open", str(path)])


def active_state(outdir: Path) -> bool:
    state = read_json(outdir / "run_state.json")
    return state.get("status") in {"starting", "running"}


def active_run_registry_path() -> Path:
    return preferences_path().with_name("active_run.json")


def run_worker_state_path(outdir: Path) -> Path:
    return outdir / ".run-worker.json"


def run_worker_is_active(outdir: Path, state: dict[str, Any] | None = None) -> bool:
    worker_state = state if state is not None else read_json(run_worker_state_path(outdir))
    if worker_state.get("status") not in {"starting", "running"}:
        return False
    try:
        pid = int(worker_state.get("worker_pid", 0))
    except (TypeError, ValueError):
        return False
    token = str(worker_state.get("worker_token", ""))
    state_path = str(run_worker_state_path(outdir).resolve())
    return process_matches(pid, "simplseq", "run-headless", token, state_path)


def reconcile_run_worker(outdir: Path) -> tuple[dict[str, Any], bool]:
    state_path = run_worker_state_path(outdir)
    worker_state = read_json(state_path)
    active = run_worker_is_active(outdir, worker_state)
    if worker_state.get("status") in {"starting", "running"} and not active:
        completed_at = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        worker_state = {
            **worker_state,
            "status": "failed",
            "completed_at": completed_at,
            "detail": "The workflow worker stopped before the run completed.",
        }
        worker_state.pop("worker_pid", None)
        write_json(state_path, worker_state)
        run_state_path = outdir / "run_state.json"
        run_state = read_json(run_state_path)
        if run_state.get("status") in {"starting", "running"}:
            write_json(run_state_path, {
                **run_state,
                "status": "failed",
                "completed_at": completed_at,
                "detail": worker_state["detail"],
            })
    return worker_state, active


def process_active(outdir: Path) -> bool:
    _state, active = reconcile_run_worker(outdir.resolve())
    return active


def active_process_outdir() -> Path | None:
    registry = read_json(active_run_registry_path())
    raw_outdir = str(registry.get("outdir", "")).strip()
    if raw_outdir:
        outdir = Path(raw_outdir).expanduser().resolve()
        if process_active(outdir):
            return outdir
    with RUN_LOCK:
        for key, process in list(RUN_PROCESSES.items()):
            if process.poll() is None:
                return Path(key)
            RUN_PROCESSES.pop(key, None)
    return None


def any_process_active() -> bool:
    return active_process_outdir() is not None


def claim_run_slot(outdir: Path) -> str | None:
    """Atomically reserve the single local workflow slot across GUI instances."""
    registry_path = active_run_registry_path()
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    payload = {
        "schema_version": 1,
        "status": "reserving",
        "reservation_token": token,
        "outdir": str(outdir.resolve()),
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    }
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    for _attempt in range(3):
        try:
            descriptor = os.open(registry_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        except FileExistsError:
            existing = read_json(registry_path)
            existing_outdir = str(existing.get("outdir", "")).strip()
            if existing_outdir and process_active(Path(existing_outdir)):
                return None
            try:
                age_seconds = max(0.0, time.time() - registry_path.stat().st_mtime)
            except OSError:
                age_seconds = 0.0
            if existing.get("status") == "reserving" and age_seconds < 120:
                return None
            registry_path.unlink(missing_ok=True)
            continue
        try:
            os.write(descriptor, encoded)
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        return token
    return None


def release_run_slot(token: str) -> None:
    registry_path = active_run_registry_path()
    registry = read_json(registry_path)
    if token and token in {registry.get("reservation_token"), registry.get("worker_token")}:
        registry_path.unlink(missing_ok=True)


def stop_tracked_process(outdir: Path) -> bool:
    outdir = outdir.resolve()
    worker_state_path = run_worker_state_path(outdir)
    worker_state, active = reconcile_run_worker(outdir)
    if not active:
        return False
    pid = int(worker_state["worker_pid"])
    token = str(worker_state.get("worker_token", ""))
    if not terminate_process_group(pid):
        return False

    with RUN_LOCK:
        RUN_PROCESSES.pop(str(outdir), None)

    completed_at = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    write_json(worker_state_path, {
        **worker_state,
        "status": "stopped",
        "completed_at": completed_at,
        "detail": "Run stopped by the user.",
    })
    release_run_slot(token)

    state_file = outdir / "run_state.json"
    state = read_json(state_file)
    state.update(
        {
            "status": "stopped",
            "completed_at": completed_at,
            "outdir": str(outdir),
            "detail": "Run stopped by user from the malaria-amplicon-nf browser app.",
        }
    )
    write_json(state_file, state)
    return True



def payload_resume_enabled(data: dict[str, Any]) -> bool:
    if "resume" in data:
        return bool_payload(data, "resume", False)
    if "clean" in data:
        return not bool_payload(data, "clean", True)
    return False


def headless_run_command(
    samples: Path,
    outdir: Path,
    data: dict[str, Any],
    *,
    worker_token: str = "",
    worker_state: Path | None = None,
) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "simplseq",
        "run-headless",
        "--samples",
        str(samples),
        "--out",
        str(outdir),
    ]
    cpus = int_payload(data, "cpus", 0)
    memory = str(data.get("memory", "")).strip()
    if cpus:
        command.extend(["--cpus", str(cpus)])
    if memory:
        command.extend(["--memory", memory])
    if not payload_resume_enabled(data):
        command.append("--no-resume")
    if bool_payload(data, "dry_run", False):
        command.append("--dry-run")
    kelt_barcode_map = str(data.get("kelt_barcode_map", "") or "").strip()
    if kelt_barcode_map:
        command.extend(["--kelt", "--kelt-barcode-map", kelt_barcode_map])
    if worker_token and worker_state is not None:
        command.extend(["--worker-token", worker_token, "--worker-state", str(worker_state)])
    return command


def start_run_process(
    root: Path,
    samples: Path,
    outdir: Path,
    data: dict[str, Any],
    *,
    worker_token: str,
) -> subprocess.Popen[str]:
    outdir.mkdir(parents=True, exist_ok=True)
    logs_dir = outdir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    worker_state = run_worker_state_path(outdir)
    command = headless_run_command(
        samples,
        outdir,
        data,
        worker_token=worker_token,
        worker_state=worker_state,
    )
    env = local_runtime_env(root)
    env["SIMPLSEQ_PROJECT_ROOT"] = str(root)
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    stdout_handle = (logs_dir / "flask-run.stdout.log").open("w", encoding="utf-8")
    stderr_handle = (logs_dir / "flask-run.stderr.log").open("w", encoding="utf-8")
    kwargs: dict[str, Any] = {}
    if os.name == "nt":
        kwargs["creationflags"] = (
            getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
            | getattr(subprocess, "CREATE_NO_WINDOW", 0)
        )
    else:
        kwargs["start_new_session"] = True
    started_at = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    write_json(worker_state, {
        "schema_version": 1,
        "status": "starting",
        "started_at": started_at,
        "worker_token": worker_token,
        "outdir": str(outdir.resolve()),
    })
    process = subprocess.Popen(
        command,
        cwd=root,
        env=env,
        text=True,
        stdout=stdout_handle,
        stderr=stderr_handle,
        **kwargs,
    )
    stdout_handle.close()
    stderr_handle.close()
    write_json(worker_state, {
        "schema_version": 1,
        "status": "starting",
        "started_at": started_at,
        "worker_pid": process.pid,
        "worker_token": worker_token,
        "outdir": str(outdir.resolve()),
    })
    write_json(active_run_registry_path(), {
        "schema_version": 1,
        "status": "running",
        "worker_pid": process.pid,
        "worker_token": worker_token,
        "worker_state": str(worker_state),
        "outdir": str(outdir.resolve()),
        "started_at": started_at,
    })
    with RUN_LOCK:
        RUN_PROCESSES[str(outdir)] = process
    return process


def create_app(root: Path | None = None, workspace_root: Path | None = None) -> Flask:
    app_root = (root or project_root()).resolve()
    workspace = (workspace_root or Path.cwd()).expanduser().resolve()
    gui_root = Path(__file__).resolve().parent
    app = Flask(
        __name__,
        static_folder=str(gui_root / "static"),
        template_folder=str(gui_root / "templates"),
    )
    app.config["SIMPLSEQ_PROJECT_ROOT"] = app_root
    app.config["SIMPLSEQ_WORKSPACE_ROOT"] = workspace
    app.config["TEMPLATES_AUTO_RELOAD"] = True

    @app.after_request
    def disable_local_ui_cache(response):
        if request.path == "/" or request.path.startswith("/static/"):
            response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
            response.headers["Pragma"] = "no-cache"
            response.headers["Expires"] = "0"
        return response

    @app.before_request
    def protect_loopback_app():
        try:
            hostname = (urlsplit(f"//{request.host}").hostname or "").lower()
        except ValueError:
            abort(400)
        if hostname not in {"127.0.0.1", "localhost", "::1"}:
            abort(400)
        if request.method in {"GET", "HEAD", "OPTIONS"} or not request.path.startswith("/api/"):
            return None
        if request.headers.get("Sec-Fetch-Site", "").lower() == "cross-site":
            abort(403)
        origin = request.headers.get("Origin", "").strip()
        if origin and urlsplit(origin).netloc.lower() != request.host.lower():
            abort(403)
        return None

    @app.get("/")
    def index():
        static_paths = [gui_root / "static" / "css" / "app.css", gui_root / "static" / "js" / "app.js"]
        asset_version = str(max((path.stat().st_mtime_ns for path in static_paths if path.exists()), default=0))
        return render_template(
            "index.html",
            workspace_root=str(workspace),
            asset_version=asset_version,
        )

    @app.get("/assets/<path:filename>")
    def assets(filename: str):
        assets_dir = (app_root / "assets").resolve()
        requested = (assets_dir / filename).resolve()
        try:
            requested.relative_to(assets_dir)
        except ValueError:
            abort(404)
        return send_from_directory(assets_dir, filename)

    @app.get("/api/health")
    def api_health():
        return jsonify(
            {
                "ok": True,
                "app": "malaria-amplicon-nf",
                "version": f"v{__version__}-dev" if __version__ == "0.1.0" else f"v{__version__}",
                "app_root": str(app_root),
                "workspace_root": str(workspace),
                "path_style": path_style(),
                "common_paths": common_paths(workspace, app_root),
            }
        )

    @app.get("/api/preferences")
    def api_preferences_get():
        path = preferences_path()
        return jsonify({"ok": True, "settings": sanitized_preferences(read_json(path))})

    @app.post("/api/preferences")
    def api_preferences_save():
        data = request.get_json(silent=True) or {}
        settings = sanitized_preferences(data.get("settings", data))
        write_json(preferences_path(), settings)
        return jsonify({"ok": True})

    @app.get("/api/browse")
    def api_browse():
        requested = resolve_app_path(workspace, request.args.get("path"), workspace)
        path = nearest_existing_directory(requested, workspace)
        return jsonify({"ok": True, **browse_payload(path)})

    @app.post("/api/select-folder")
    def api_select_folder():
        data = request.get_json(silent=True) or {}
        initial = resolve_app_path(workspace, data.get("initial"), workspace)
        prompt = str(data.get("prompt") or "Select a folder")
        allow_new_folder = bool_payload(data, "allow_new_folder", False)
        return jsonify(select_folder_dialog(initial, prompt=prompt, allow_new_folder=allow_new_folder))

    @app.post("/api/select-file")
    def api_select_file():
        data = request.get_json(silent=True) or {}
        initial_raw = str(data.get("initial") or "").strip()
        initial = resolve_app_path(workspace, initial_raw, workspace) if initial_raw else workspace
        prompt = str(data.get("prompt") or "Select metadata file")
        kind = str(data.get("kind") or "metadata").strip().lower()
        if kind not in {"metadata", "kelt"}:
            return json_error("Unsupported file picker type.", 400)
        return jsonify(select_file_dialog(initial, prompt=prompt, kind=kind))

    @app.post("/api/kelt/inspect")
    def api_kelt_inspect():
        data = request.get_json(silent=True) or {}
        path_raw = str(data.get("path", "") or "").strip()
        if not path_raw:
            return json_error("Choose a KELT barcode map first.", 400)
        path = resolve_app_path(workspace, path_raw, path_raw)
        try:
            return jsonify({"ok": True, **inspect_kelt_barcode_map(path)})
        except (OSError, ValueError) as exc:
            return json_error(str(exc), 400)

    @app.post("/api/metadata/inspect")
    def api_metadata_inspect():
        data = request.get_json(silent=True) or {}
        metadata_contract = normalize_metadata_contract(data.get("metadata_contract", {}))
        path_raw = str(data.get("path", "") or "").strip()
        if not path_raw:
            return json_error("Choose a metadata file first.", 400)
        path = resolve_app_path(workspace, path_raw, path_raw)
        try:
            catalog = inspect_metadata(
                path,
                str(data.get("sheet", "") or "").strip(),
                date_order=str(data.get("date_order", "auto") or "auto").lower(),
                column_overrides=metadata_contract["columns"],
            )
        except (OSError, RuntimeError, ValueError) as exc:
            return json_error(str(exc), 400)
        summary = catalog.summary()
        summary["detection_values"] = detection_value_options(
            catalog.value_counts.get("metadata_pcr", {}),
            metadata_contract["detection_value_map"],
        )
        metadata_contract["columns"] = dict(catalog.columns)
        metadata_contract["detection_value_map"] = {
            item["value"]: item["state"] for item in summary["detection_values"]
        }
        return jsonify({
            "ok": not any(issue["severity"] == "error" for issue in summary["issues"]),
            "metadata_contract": metadata_contract,
            **summary,
        })

    @app.get("/api/metadata/template")
    def api_metadata_template():
        header = (
            "participant_id,visit_date,visit_month,pcr_result,season,"
            "visit_status,age,sex\r\n"
        )
        return send_file(
            BytesIO(header.encode("utf-8")),
            mimetype="text/csv",
            as_attachment=True,
            download_name="simplseq_metadata_template.csv",
        )

    @app.post("/api/scan")
    def api_scan():
        data = request.get_json(silent=True) or {}
        try:
            fastq_dir = resolve_fastq_folder(data.get("fastq_dir"))
        except ValueError as exc:
            return json_error(str(exc), 400)
        metadata_raw = str(data.get("metadata_path", "") or "").strip()
        metadata_autodetected = not metadata_raw
        metadata_path = (
            resolve_app_path(workspace, metadata_raw, metadata_raw)
            if metadata_raw
            else discover_metadata_file(fastq_dir)
        )
        metadata_sheet = str(data.get("metadata_sheet", "") or "").strip()
        metadata_date_order = str(data.get("metadata_date_order", "auto") or "auto").lower()
        metadata_contract = normalize_metadata_contract(data.get("metadata_contract", {}))
        include_pool = bool_payload(data, "include_pool_in_sample_id", False)
        absolute = bool_payload(data, "absolute_paths", True)
        write_samples = bool_payload(data, "write_samples", True)
        samples_out = (
            resolve_app_path(workspace, data.get("samples_out"), "samples.csv")
            if write_samples or str(data.get("samples_out") or "").strip()
            else None
        )
        try:
            fallback_year, fallback_day = collection_date_defaults(data)
        except ValueError as exc:
            return json_error(str(exc), 400)
        if metadata_path is not None and not metadata_path.exists():
            return json_error(f"Metadata file not found: {metadata_path}", 400)
        fallback_year_autodetected = False
        if not fallback_year and metadata_path is not None:
            try:
                fallback_year = infer_metadata_year(
                    metadata_path,
                    metadata_sheet,
                    date_order=metadata_date_order,
                    column_overrides=metadata_contract["columns"],
                )
            except (OSError, RuntimeError, ValueError) as exc:
                return json_error(f"Metadata could not be used: {exc}", 400)
            fallback_year_autodetected = bool(fallback_year)

        scan = scan_fastqs(
            fastq_dir,
            include_pool_in_sample_id=include_pool,
            fallback_collection_year=fallback_year,
            fallback_collection_day=fallback_day,
        )
        written = False
        duplicates: list[str] = scan.duplicate_sample_ids
        count = 0
        try:
            preview, scan_rows, metadata_counts, date_counts = sample_rows_preview(
                scan,
                samples_out.parent if samples_out is not None else workspace,
                absolute=absolute,
                metadata_path=metadata_path,
                metadata_sheet=metadata_sheet,
                metadata_date_order=metadata_date_order,
                metadata_columns=metadata_contract["columns"],
            )
        except (OSError, RuntimeError, ValueError) as exc:
            return json_error(f"Metadata could not be used: {exc}", 400)
        if write_samples and samples_out is not None and not duplicates:
            count, duplicates = write_samples_csv(
                fastq_dir,
                samples_out,
                include_pool_in_sample_id=include_pool,
                absolute=absolute,
                fallback_collection_year=fallback_year,
                fallback_collection_day=fallback_day,
                metadata_path=metadata_path,
                metadata_sheet=metadata_sheet,
                metadata_date_order=metadata_date_order,
                metadata_columns=metadata_contract["columns"],
            )
            written = not duplicates
        response = scan_json(
            scan,
            workspace,
            fallback_collection_year=fallback_year,
            fallback_collection_day=fallback_day,
        )
        response.update(
            {
                "ok": True,
                "samples_out": str(samples_out) if samples_out is not None else "",
                "samples_relative": rel_or_abs(workspace, samples_out) if samples_out is not None else "samples.csv",
                "samples_written": written,
                "sample_rows_written": count,
                "sample_rows_previewed": len(scan.pairs),
                "metadata_path": str(metadata_path) if metadata_path is not None else "",
                "metadata_autodetected": bool(metadata_path is not None and metadata_autodetected),
                "fallback_year_autodetected": fallback_year_autodetected,
                "metadata_sheet": metadata_sheet,
                "metadata_date_order": metadata_date_order,
                "metadata_contract": metadata_contract,
                "metadata_match_counts": metadata_counts,
                "sample_preview": preview,
                "scan_rows": scan_rows,
                **date_counts,
            }
        )
        return jsonify(response)

    @app.post("/api/check")
    def api_check():
        data = request.get_json(silent=True) or {}
        samples = data.get("samples")
        samples_path = resolve_app_path(workspace, samples, "samples.csv") if samples else None
        outdir = resolve_app_path(workspace, data.get("outdir"), "results")
        rows = check_environment(app_root, samples_path, outdir=outdir)
        failed = sum(1 for row in rows if row.get("status") not in {"ok", "warn"})
        return jsonify({"ok": failed == 0, "failed": failed, "checks": rows})

    @app.post("/api/run")
    def api_run():
        data = request.get_json(silent=True) or {}
        metadata_path: Path | None = None
        metadata_sheet = ""
        metadata_date_order = "auto"
        metadata_contract = normalize_metadata_contract({})
        fallback_year = ""
        fallback_day = DEFAULT_COLLECTION_DAY
        selected_libraries: list[str] = []
        try:
            output_parent = resolve_run_output_parent(data.get("outdir"))
        except ValueError as exc:
            return json_error(str(exc), 400)
        resume = payload_resume_enabled(data)
        if any_process_active():
            return json_error("A malaria-amplicon-nf run is already active.", 409, outdir=str(active_process_outdir() or output_parent))
        try:
            outdir = allocate_run_outdir(output_parent, str(data.get("run_name", "")), reuse_existing=resume)
        except OSError as exc:
            return json_error(f"Output folder is not writable: {exc}", 400, outdir=str(output_parent))
        dry_run = bool_payload(data, "dry_run", False)
        fastq_raw = str(data.get("fastq_dir", "") or "").strip()
        samples = outdir / "samples.csv"
        if fastq_raw:
            try:
                fastq_dir = resolve_fastq_folder(fastq_raw)
            except ValueError as exc:
                return json_error(str(exc), 400, outdir=str(outdir))
            metadata_raw = str(data.get("metadata_path", "") or "").strip()
            metadata_path = (
                resolve_app_path(workspace, metadata_raw, metadata_raw)
                if metadata_raw
                else discover_metadata_file(fastq_dir)
            )
            metadata_sheet = str(data.get("metadata_sheet", "") or "").strip()
            metadata_date_order = str(data.get("metadata_date_order", "auto") or "auto").lower()
            metadata_contract = normalize_metadata_contract(data.get("metadata_contract", {}))
            kelt_barcode_raw = str(data.get("kelt_barcode_map", "") or "").strip()
            kelt_barcode_map = (
                resolve_app_path(workspace, kelt_barcode_raw, kelt_barcode_raw)
                if kelt_barcode_raw
                else None
            )
            include_pool = bool_payload(data, "include_pool_in_sample_id", False)
            absolute = bool_payload(data, "absolute_paths", True)
            try:
                fallback_year, fallback_day = collection_date_defaults(data)
            except ValueError as exc:
                return json_error(str(exc), 400)
            if metadata_path is not None and not metadata_path.exists():
                return json_error(f"Metadata file not found: {metadata_path}", 400, outdir=str(outdir))
            if kelt_barcode_map is not None and not kelt_barcode_map.exists():
                return json_error(f"KELT barcode map not found: {kelt_barcode_map}", 400, outdir=str(outdir))
            if not fallback_year and metadata_path is not None:
                try:
                    fallback_year = infer_metadata_year(
                        metadata_path,
                        metadata_sheet,
                        date_order=metadata_date_order,
                        column_overrides=metadata_contract["columns"],
                    )
                except (OSError, RuntimeError, ValueError) as exc:
                    return json_error(f"Metadata could not be used: {exc}", 400, outdir=str(outdir))
            raw_libraries = data.get("libraries", [])
            if isinstance(raw_libraries, str):
                raw_libraries = [raw_libraries]
            selected_libraries = [
                str(value).strip()
                for value in raw_libraries
                if str(value).strip() and str(value).strip().lower() != "all"
            ]
            try:
                count, duplicates = write_samples_csv(
                    fastq_dir,
                    samples,
                    include_pool_in_sample_id=include_pool,
                    absolute=absolute,
                    fallback_collection_year=fallback_year,
                    fallback_collection_day=fallback_day,
                    metadata_path=metadata_path,
                    metadata_sheet=metadata_sheet,
                    metadata_date_order=metadata_date_order,
                    metadata_columns=metadata_contract["columns"],
                    kelt_barcode_map=kelt_barcode_map,
                    libraries=selected_libraries,
                )
            except (OSError, RuntimeError, ValueError) as exc:
                return json_error(str(exc), 400, outdir=str(outdir), samples=str(samples))
            if duplicates:
                return json_error(
                    f"Duplicate sample IDs: {', '.join(duplicates[:8])}",
                    400,
                    outdir=str(outdir),
                    samples=str(samples),
                )
            if count == 0 and not dry_run:
                return json_error(f"No paired FASTQ files found in {fastq_dir}", 400, outdir=str(outdir))
        elif data.get("samples"):
            samples = resolve_app_path(workspace, data.get("samples"), "samples.csv")
        if not samples.exists() and not dry_run:
            return json_error(f"Sample sheet not found: {samples}", 400)
        if samples.exists():
            try:
                freeze_run_configuration(
                    outdir,
                    samples,
                    metadata_path=metadata_path,
                    metadata_sheet=metadata_sheet,
                    metadata_date_order=metadata_date_order,
                    fallback_year=fallback_year,
                    fallback_day=fallback_day,
                    metadata_contract=metadata_contract,
                    selected_libraries=selected_libraries,
                    request_payload=data,
                )
            except (OSError, RuntimeError, ValueError) as exc:
                return json_error(
                    f"Run inputs could not be frozen reproducibly: {exc}",
                    400,
                    outdir=str(outdir),
                )
        worker_token = claim_run_slot(outdir)
        if worker_token is None:
            return json_error(
                "A malaria-amplicon-nf run is already active.",
                409,
                outdir=str(active_process_outdir() or outdir),
            )
        try:
            process = start_run_process(
                app_root,
                samples,
                outdir,
                data,
                worker_token=worker_token,
            )
        except Exception as exc:  # pragma: no cover - reported to browser
            release_run_slot(worker_token)
            return json_error(str(exc), 500)
        return jsonify(
            {
                "ok": True,
                "pid": process.pid,
                "outdir": str(outdir),
                "output_parent": str(output_parent),
                "samples": str(samples),
                "dry_run": dry_run,
                "status_url": url_for("api_status", out=str(outdir)),
            }
        )

    @app.post("/api/stop-run")
    def api_stop_run():
        data = request.get_json(silent=True) or {}
        outdir = resolve_app_path(workspace, data.get("outdir") or data.get("out"), "results")
        stopped = stop_tracked_process(outdir)
        if not stopped:
            return json_error("No active malaria-amplicon-nf run was found.", 409, outdir=str(outdir))
        state = read_json(outdir / "run_state.json")
        summary = progress_summary(outdir)
        summary = dict(summary)
        summary["status"] = "stopped"
        return jsonify(
            {
                "ok": True,
                "stopped": True,
                "active": False,
                "outdir": str(outdir),
                "state": state,
                "summary": summary,
            }
        )

    @app.get("/api/active-run")
    def api_active_run():
        outdir = active_process_outdir()
        if outdir is None:
            return jsonify({"ok": True, "active": False})
        return jsonify(
            {
                "ok": True,
                "active": True,
                "outdir": str(outdir),
                "status_url": url_for("api_status", out=str(outdir)),
            }
        )

    @app.get("/api/status")
    def api_status():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        state = read_json(outdir / "run_state.json")
        summary = progress_summary(outdir)
        tracked_active = process_active(outdir)
        if tracked_active:
            state = dict(state)
            summary = dict(summary)
            status = str(state.get("status") or summary.get("status") or "pending")
            if status == "pending":
                status = "starting"
            state.setdefault("status", status)
            summary["status"] = status
            if str(summary.get("current_stage", "")) not in STAGES:
                summary["current_stage"] = "pending"
        if active_state(outdir) and not tracked_active:
            state = read_json(outdir / "run_state.json")
            summary = dict(summary)
            summary["status"] = str(state.get("status") or "failed")
        return jsonify(
            {
                "ok": True,
                "outdir": str(outdir),
                "state": state,
                "summary": summary,
                "active": tracked_active,
            }
        )

    @app.get("/api/progress")
    def api_progress():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        events = read_events(outdir / "progress.jsonl")
        return jsonify(
            {
                "ok": True,
                "outdir": str(outdir),
                "events": events,
                "summary": progress_summary(outdir),
            }
        )

    @app.get("/api/results")
    def api_results():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        manifest = result_files_with_downloads(workspace, outdir)
        return jsonify({"ok": True, "outdir_exists": safe_is_dir(outdir), **manifest})

    @app.get("/api/result-table/<file_key>")
    def api_result_table(file_key: str):
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        path, label = result_file_path(outdir, file_key)
        suffix = path.suffix.lower()
        if suffix not in {".tsv", ".csv"}:
            abort(404)
        delimiter = "\t" if suffix == ".tsv" else ","
        limit = max(1, min(int_payload(request.args, "limit", 500), 1000))
        rows: list[dict[str, str]] = []
        truncated = False
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            columns = list(reader.fieldnames or [])
            for index, row in enumerate(reader):
                if index >= limit:
                    truncated = True
                    break
                rows.append({column: row.get(column, "") for column in columns})
        return jsonify(
            {
                "ok": True,
                "label": label,
                "filename": path.name,
                "size": human_bytes(path.stat().st_size),
                "columns": columns,
                "rows": rows,
                "shown_rows": len(rows),
                "truncated": truncated,
            }
        )

    @app.post("/api/open-folder")
    def api_open_folder():
        data = request.get_json(silent=True) or {}
        path = resolve_app_path(workspace, data.get("path"), "results").resolve()
        if not path.exists() or not path.is_dir():
            abort(404)
        try:
            open_directory(path)
        except OSError as exc:
            return jsonify({"ok": False, "error": f"Could not open folder: {exc}"}), 500
        return jsonify({"ok": True, "path": str(path)})

    @app.get("/api/logs")
    def api_logs():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        max_bytes = max(1000, min(int_payload(request.args, "max_bytes", 50000), 250000))
        log_text, truncated = run_log_tail(outdir, max_bytes)
        return jsonify(
            {
                "ok": True,
                "outdir": str(outdir),
                "path": str(outdir / "technical_log.txt"),
                "text": clean_log_text(log_text),
                "truncated": truncated,
            }
        )

    @app.get("/download/<file_key>")
    def download_result(file_key: str):
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        path, label = result_file_path(outdir, file_key)
        inline = request.args.get("inline") == "1" and label in VIEWABLE_REPORT_LABELS
        return send_file(path, as_attachment=not inline, download_name=path.name)

    @app.post("/api/save-result/<file_key>")
    def api_save_result(file_key: str):
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        path, _label = result_file_path(outdir, file_key)
        target = unique_download_path(path.name)
        shutil.copy2(path, target)
        return jsonify({"ok": True, "path": str(target), "filename": target.name})

    @app.get("/download-bundle")
    def download_bundle():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        bundle_paths = bundle_result_paths(outdir)
        if not bundle_paths:
            abort(404)
        archive = BytesIO()
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
            for path, arcname in bundle_paths:
                zip_handle.write(path, arcname)
        archive.seek(0)
        bundle_name = f"{outdir.name or 'simplseq-results'}-output-bundle.zip"
        return send_file(archive, as_attachment=True, download_name=bundle_name, mimetype="application/zip")

    @app.post("/api/save-bundle")
    def api_save_bundle():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        bundle_paths = bundle_result_paths(outdir)
        if not bundle_paths:
            abort(404)
        bundle_name = f"{outdir.name or 'malaria-amplicon-nf'}-output-bundle.zip"
        target = unique_download_path(bundle_name)
        with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
            for path, arcname in bundle_paths:
                zip_handle.write(path, arcname)
        return jsonify({"ok": True, "path": str(target), "filename": target.name})

    # ------------------------------------------------------------------
    # Shared downstream analysis inputs
    # ------------------------------------------------------------------

    def downstream_analysis_mode(value: Any) -> str:
        normalized = str(value or "primary").strip().lower()
        if normalized in {"cdhit98", "cdhit989", "cdhit_summed", "summed"}:
            return "cdhit98"
        return "primary"

    def analysis_stem(mode: str) -> str:
        normalized = downstream_analysis_mode(mode)
        return "cdhit98" if normalized == "cdhit98" else "primary"

    def analysis_input_path(outdir: Path, mode: str) -> Path:
        return outdir / "analysis_input" / f"{analysis_stem(mode)}_seqtab.tsv"

    def analysis_source_path(outdir: Path, mode: str) -> Path:
        if downstream_analysis_mode(mode) == "cdhit98":
            return outdir / "analysis_input" / "cdhit98_unfiltered_seqtab.tsv"
        return outdir / "run_dada2" / "seqtab_cigar.tsv"

    def analysis_long_path(outdir: Path, mode: str) -> Path:
        return outdir / "analysis_input" / f"{analysis_stem(mode)}_alleles.tsv"

    def analysis_samples_path(outdir: Path, mode: str) -> Path:
        return outdir / "analysis_input" / f"{analysis_stem(mode)}_samples.csv"

    def analysis_summary_path(outdir: Path, mode: str) -> Path:
        return outdir / "analysis_input" / f"{analysis_stem(mode)}_summary.json"

    def table_dimensions(path: Path) -> tuple[int, int]:
        if not safe_is_file(path):
            return 0, 0
        try:
            with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, [])
                rows = sum(1 for _ in reader)
            return rows, max(len(header) - 1, 0)
        except OSError:
            return 0, 0

    def analysis_input_payload(outdir: Path) -> dict[str, Any]:
        primary = analysis_input_path(outdir, "primary")
        sensitivity = analysis_input_path(outdir, "cdhit98")
        primary_samples, primary_alleles = table_dimensions(primary)
        sensitivity_samples, sensitivity_alleles = table_dimensions(sensitivity)
        primary_summary = read_json(analysis_summary_path(outdir, "primary"))
        summary = read_json(analysis_summary_path(outdir, "cdhit98"))
        return {
            "primary": {
                "available": safe_is_file(primary),
                "path": str(primary),
                "samples": primary_samples,
                "alleles": primary_alleles,
                "label": "Exact CIGAR alleles",
                "summary": primary_summary,
                "table_url": url_for("api_analysis_table_preview", out=str(outdir), mode="primary")
                if safe_is_file(primary) else "",
                "download_url": url_for("download_analysis_table", out=str(outdir), mode="primary")
                if safe_is_file(primary) else "",
            },
            "cdhit98": {
                "available": safe_is_file(sensitivity),
                "path": str(sensitivity),
                "samples": sensitivity_samples,
                "alleles": sensitivity_alleles,
                "label": "CD-HIT 98.9% sensitivity",
                "summary": summary,
                "table_url": url_for("api_analysis_table_preview", out=str(outdir), mode="cdhit98")
                if safe_is_file(sensitivity) else "",
                "download_url": url_for("download_analysis_table", out=str(outdir), mode="cdhit98")
                if safe_is_file(sensitivity) else "",
            },
        }

    @app.get("/api/analysis-table/status")
    def api_analysis_table_status():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        return jsonify({"ok": True, "outdir": str(outdir), "inputs": analysis_input_payload(outdir)})

    @app.get("/api/asv-filtering-audit")
    def api_asv_filtering_audit():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        path = outdir / "reports" / "asv_filtering_summary.txt"
        if not safe_is_file(path):
            abort(404)
        payload = parse_asv_filtering_summary_text(path.read_text(encoding="utf-8", errors="replace"))
        payload.update({
            "ok": True,
            "filename": path.name,
            "size": human_bytes(path.stat().st_size),
        })
        return jsonify(payload)

    @app.post("/api/analysis-table/build")
    def api_analysis_table_build():
        data = request.get_json(silent=True) or {}
        outdir = resolve_app_path(workspace, data.get("outdir"), "results")
        mode = downstream_analysis_mode(data.get("mode"))
        try:
            min_abundance_pct = float(data.get("min_abundance_pct", 1.0))
        except (TypeError, ValueError):
            return json_error("Minimum allele abundance must be between 0 and 100.", 400)
        denominator = str(data.get("abundance_denominator", "locus")).strip().lower()
        if not 0 <= min_abundance_pct <= 100:
            return json_error("Minimum allele abundance must be between 0 and 100.", 400)
        if denominator not in {"locus", "sample"}:
            return json_error("Abundance denominator must be locus or sample.", 400)

        analysis_dir = outdir / "analysis_input"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        samples = outdir / "samples.csv"
        if not safe_is_file(samples):
            return json_error("The completed run has no samples.csv file.", 400)

        if mode == "cdhit98":
            counts = outdir / "cdhit" / "cdhit_cluster_counts.tsv"
            membership = outdir / "cdhit" / "cdhit_cluster_membership.tsv"
            if not safe_is_file(counts) or not safe_is_file(membership):
                return json_error(
                    "CD-HIT QC outputs are not available for this run. Complete the main workflow first.",
                    400,
                )
            source = analysis_source_path(outdir, mode)
            cluster_command = [
                sys.executable,
                str(app_root / "workflow" / "scripts" / "build_analysis_table.py"),
                "--mode", "summed",
                "--cigar", str(outdir / "run_dada2" / "seqtab_cigar.tsv"),
                "--cdhit-counts", str(counts),
                "--cdhit-membership", str(membership),
                "--out", str(source),
                "--summary", str(analysis_dir / "cdhit98_cluster_summary.json"),
            ]
            cluster_completed = subprocess.run(
                cluster_command, cwd=app_root, env=analysis_runtime_env(app_root),
                text=True, capture_output=True, check=False, timeout=600,
            )
            if cluster_completed.returncode != 0:
                detail = (cluster_completed.stderr or cluster_completed.stdout or "Could not build the CD-HIT table.").strip()
                return json_error(detail, 400)
        else:
            source = analysis_source_path(outdir, mode)
            if not safe_is_file(source):
                return json_error("The completed run has no exact-CIGAR count table.", 400)

        stem = analysis_stem(mode)
        summary_path = analysis_summary_path(outdir, mode)
        metadata_settings = frozen_metadata_settings(outdir)
        longitudinal_dir = analysis_dir / f"{stem}_longitudinal"
        command = [
            sys.executable,
            str(app_root / "workflow" / "scripts" / "prepare_analysis_table.py"),
            "--input", str(source),
            "--samples", str(samples),
            "--output-wide", str(analysis_input_path(outdir, mode)),
            "--output-long", str(analysis_long_path(outdir, mode)),
            "--output-samples", str(analysis_samples_path(outdir, mode)),
            "--filter-summary", str(analysis_dir / f"{stem}_abundance_filter.tsv"),
            "--replicate-summary", str(analysis_dir / f"{stem}_replicate_merge.tsv"),
            "--summary", str(summary_path),
            "--mode", mode,
            "--min-abundance-pct", str(min_abundance_pct),
            "--min-locus-reads", "100",
            "--min-biological-samples", "2",
            "--denominator", denominator,
            "--metadata-sheet", metadata_settings["sheet"],
            "--metadata-date-order", metadata_settings["date_order"],
            "--metadata-profile", metadata_settings["profile_path"],
            "--fallback-year", metadata_settings["fallback_year"],
            "--fallback-day", metadata_settings["fallback_day"],
            "--longitudinal-dir", str(longitudinal_dir),
            "--analysis-profile", str(outdir / "configuration" / f"analysis_profile_{stem}.json"),
        ]
        if metadata_settings["path"]:
            command.extend(["--metadata", metadata_settings["path"]])
        try:
            completed = subprocess.run(
                command,
                cwd=app_root,
                env=analysis_runtime_env(app_root),
                text=True,
                capture_output=True,
                check=False,
                timeout=600,
            )
        except subprocess.TimeoutExpired:
            return json_error(
                "Building the CD-HIT sensitivity input exceeded 10 minutes. "
                "The primary exact-allele results were not changed.",
                504,
            )
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "Could not build the CD-HIT sensitivity table.").strip()
            return json_error(detail, 400)
        summary = read_json(summary_path)
        if mode == "cdhit98":
            membership_copy = analysis_dir / "cdhit98_membership.tsv"
            shutil.copy2(outdir / "cdhit" / "cdhit_cluster_membership.tsv", membership_copy)
            cdhit_summary = read_json(outdir / "cdhit" / "cdhit_summary.json")
            summary.update({
                "identity_threshold": cdhit_summary.get("identity_threshold", 0.989),
                "clustering_scope": cdhit_summary.get("clustering_scope", "within_locus"),
                "input_asvs": cdhit_summary.get("input_asvs", 0),
                "clusters": cdhit_summary.get("clusters", summary.get("output", {}).get("alleles", 0)),
                "membership": str(membership_copy),
                "primary_input": str(analysis_input_path(outdir, "primary")),
            })
        write_json(summary_path, summary)
        return jsonify({
            "ok": True,
            "mode": mode,
            "summary": summary,
            "inputs": analysis_input_payload(outdir),
        })

    @app.get("/api/analysis-table/table")
    def api_analysis_table_preview():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        mode = downstream_analysis_mode(request.args.get("mode"))
        path = analysis_long_path(outdir, mode)
        if not safe_is_file(path):
            abort(404)
        limit = max(1, min(int_payload(request.args, "limit", 200), 1000))
        with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = list(reader.fieldnames or [])
            rows = []
            truncated = False
            for index, row_item in enumerate(reader):
                if index >= limit:
                    truncated = True
                    break
                rows.append({column: row_item.get(column, "") for column in columns})
        return jsonify({
            "ok": True,
            "label": "CD-HIT 98.9% sensitivity input" if mode == "cdhit98" else "Primary exact-CIGAR analysis input",
            "filename": path.name,
            "size": human_bytes(path.stat().st_size),
            "columns": columns,
            "rows": rows,
            "shown_rows": len(rows),
            "truncated": truncated,
        })

    @app.get("/download/analysis-table")
    def download_analysis_table():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        mode = downstream_analysis_mode(request.args.get("mode"))
        path = analysis_long_path(outdir, mode)
        if not safe_is_file(path):
            abort(404)
        return send_file(path, as_attachment=True, download_name=path.name)

    # ------------------------------------------------------------------
    # DINEMITES analysis endpoints
    # ------------------------------------------------------------------

    DINEMITES_FILE_KEYS = {
        "input": "dinemites_input.tsv",
        "allele_probabilities": "dinemites_allele_probabilities.tsv",
        "allele_key": "dinemites_allele_key.tsv",
        "molfoi": "dinemites_molfoi.tsv",
        "new_infections": "dinemites_new_infections.tsv",
    }
    DINEMITES_FILE_LABELS = {
        "input": "Model input",
        "allele_probabilities": "Allele probabilities",
        "allele_key": "Allele key",
        "molfoi": "Molecular force of infection",
        "new_infections": "New infection events",
    }

    def dinemites_outdir(outdir: Path, mode: str = "primary") -> Path:
        normalized = downstream_analysis_mode(mode)
        return outdir / ("dinemites_cdhit98" if normalized == "cdhit98" else "dinemites")

    def dinemites_plots_dir(outdir: Path, mode: str = "primary") -> Path:
        return dinemites_outdir(outdir, mode) / "dinemites_plots"

    def dinemites_state_path(outdir: Path, mode: str = "primary") -> Path:
        return dinemites_outdir(outdir, mode) / "dinemites_state.json"

    def dinemites_row_value(row: dict[str, str], *keys: str) -> str:
        for key in keys:
            value = row.get(key, "")
            if value not in {"", None}:
                return value
        return ""

    def dinemites_tsv_rows(path: Path, limit: int = 200) -> list[dict[str, str]]:
        if not safe_is_file(path):
            return []
        rows: list[dict[str, str]] = []
        try:
            with path.open(newline="", encoding="utf-8-sig", errors="replace") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for index, row in enumerate(reader):
                    if index >= limit:
                        break
                    rows.append({
                        str(key or "").lstrip("\ufeff"): str(value or "")
                        for key, value in row.items()
                    })
        except Exception:
            return []
        return rows

    def dinemites_input_summary(path: Path) -> dict[str, int]:
        if not safe_is_file(path):
            return {}
        subjects: set[str] = set()
        visits: set[tuple[str, str]] = set()
        loci: set[str] = set()
        alleles: set[str] = set()
        total_rows = 0
        genotype_rows = 0
        try:
            with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
                for row_item in csv.DictReader(handle, delimiter="\t"):
                    total_rows += 1
                    subject = dinemites_row_value(row_item, "subject", "participant_id")
                    collection_date = dinemites_row_value(row_item, "collection_date", "date_full")
                    locus = dinemites_row_value(row_item, "locus")
                    allele = dinemites_row_value(row_item, "allele")
                    if subject:
                        subjects.add(subject)
                        if collection_date:
                            visits.add((subject, collection_date))
                    if locus:
                        loci.add(locus)
                    if allele:
                        alleles.add(allele)
                        genotype_rows += 1
        except OSError:
            return {}
        return {
            "subjects": len(subjects),
            "participant_visits": len(visits),
            "loci": len(loci),
            "unique_alleles": len(alleles),
            "genotype_rows": genotype_rows,
            "empty_visit_rows": total_rows - genotype_rows,
            "total_rows": total_rows,
        }

    def dinemites_plot_subject(filename: str) -> str:
        stem = Path(filename).stem
        return stem[len("subject_") :] if stem.startswith("subject_") else stem

    def png_dimensions(path: Path) -> tuple[int, int]:
        try:
            with path.open("rb") as handle:
                header = handle.read(24)
            if header[:8] == b"\x89PNG\r\n\x1a\n" and header[12:16] == b"IHDR":
                return int.from_bytes(header[16:20], "big"), int.from_bytes(header[20:24], "big")
        except OSError:
            pass
        return 0, 0

    def dinemites_plot_entries(outdir: Path, mode: str = "primary") -> list[dict[str, Any]]:
        plots_dir = dinemites_plots_dir(outdir, mode)
        if not safe_is_dir(plots_dir):
            return []
        entries: list[dict[str, Any]] = []
        try:
            paths = sorted(plots_dir.glob("*.png"), key=lambda item: item.name.lower())
        except OSError:
            return []
        for path in paths:
            if not safe_is_file(path):
                continue
            size_bytes = path.stat().st_size
            filename = path.name
            width, height = png_dimensions(path)
            entries.append({
                "filename": filename,
                "subject": dinemites_plot_subject(filename),
                "exists": True,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "width": width,
                "height": height,
                "view_url": url_for(
                    "download_dinemites_plot",
                    filename=filename,
                    out=str(outdir),
                    mode=downstream_analysis_mode(mode),
                    inline=1,
                ),
                "download_url": url_for(
                    "download_dinemites_plot",
                    filename=filename,
                    out=str(outdir),
                    mode=downstream_analysis_mode(mode),
                ),
            })
        return entries

    def start_dinemites_job(
        root: Path,
        outdir: Path,
        samples: Path,
        model_type: str,
        model_settings: dict[str, int | float | str | bool],
        analysis_mode: str = "primary",
    ) -> int:
        analysis_mode = downstream_analysis_mode(analysis_mode)
        dm_dir = dinemites_outdir(outdir, analysis_mode)
        dm_dir.mkdir(parents=True, exist_ok=True)
        state_file = dinemites_state_path(outdir, analysis_mode)
        cigar_path = analysis_input_path(outdir, analysis_mode)
        prepared_samples = analysis_samples_path(outdir, analysis_mode)
        dm_input = dm_dir / "dinemites_input.tsv"
        qpcr_times = dm_dir / "qpcr_positive_genotype_missing.tsv"
        dm_preparation = dm_dir / "dinemites_preparation.json"
        log_path = dm_dir / "dinemites.log"
        stem = analysis_stem(analysis_mode)
        longitudinal_dir = outdir / "analysis_input" / f"{stem}_longitudinal"
        prepared_longitudinal = longitudinal_dir / "dinemites_observed_input.tsv"
        prepared_qpcr = longitudinal_dir / "qpcr_positive_genotype_missing.tsv"
        steps: list[dict[str, Any]] = []
        cmd_convert = [
            "Rscript",
            str(root / "workflow" / "scripts" / "simplseq_to_dinemites.R"),
            "--cigar", str(cigar_path),
            "--samples", str(prepared_samples),
            "--out", str(dm_input),
            "--preparation_summary", str(dm_preparation),
            "--min_abundance_pct", "0",
            "--abundance_denominator", "locus",
        ]
        if safe_is_file(prepared_longitudinal):
            shutil.copy2(prepared_longitudinal, dm_input)
            if safe_is_file(prepared_qpcr):
                shutil.copy2(prepared_qpcr, qpcr_times)
            longitudinal_summary_path = longitudinal_dir / "longitudinal_summary.json"
            if safe_is_file(longitudinal_summary_path):
                shutil.copy2(longitudinal_summary_path, dm_preparation)
        else:
            steps.append({
                "label": "Preparing DINEMITES input",
                "command": cmd_convert,
                "timeout_seconds": 600,
                "failure_detail": "Could not prepare the DINEMITES input. Check dinemites.log.",
            })
        cmd_run = [
            "Rscript",
            str(root / "workflow" / "scripts" / "run_dinemites.R"),
            "--input", str(dm_input),
            "--model", model_type,
            "--outdir", str(dm_dir),
            "--n_lags", str(model_settings["n_lags"]),
            "--t_lag", str(model_settings["t_lag"]),
            "--seed", str(model_settings["seed"]),
            "--refresh", str(model_settings["refresh"]),
            "--n_imputations", str(model_settings["n_imputations"]),
            "--bayesian_lag_days", str(model_settings["bayesian_lag_days"]),
            "--bayesian_chains", str(model_settings["bayesian_chains"]),
            "--bayesian_parallel_chains", str(model_settings["bayesian_parallel_chains"]),
            "--bayesian_iter_warmup", str(model_settings["bayesian_iter_warmup"]),
            "--bayesian_iter_sampling", str(model_settings["bayesian_iter_sampling"]),
            "--bayesian_adapt_delta", str(model_settings["bayesian_adapt_delta"]),
            "--bayesian_drop_out", str(model_settings["bayesian_drop_out"]).lower(),
            "--infection_general_covariates", str(model_settings["infection_general_covariates"]),
        ]
        if safe_is_file(qpcr_times):
            cmd_run.extend(["--qpcr_times", str(qpcr_times)])
        steps.append({
            "label": f"Running the {model_type} model",
            "command": cmd_run,
            "timeout_seconds": 7200,
            "failure_detail": "DINEMITES did not complete. Check dinemites.log.",
        })
        return launch_job(
            job_type="dinemites",
            root=root,
            output_dir=dm_dir,
            state_file=state_file,
            log_file=log_path,
            state_payload={
                "model": model_type,
                "analysis_mode": analysis_mode,
                **model_settings,
                "outdir": str(dm_dir),
            },
            steps=steps,
            env=analysis_runtime_env(root),
        )

    @app.post("/api/dinemites/run")
    def api_dinemites_run():
        data = request.get_json(silent=True) or {}
        model_type = str(data.get("model_type", "simple")).strip()
        if model_type not in {"simple", "clustering", "bayesian"}:
            return json_error(f"Invalid DINEMITES model type: {model_type}", 400)
        try:
            model_settings = dinemites_model_settings(data)
        except ValueError as exc:
            return json_error(str(exc), 400)

        outdir = resolve_app_path(workspace, data.get("outdir"), "results")
        samples = resolve_app_path(workspace, data.get("samples"), "samples.csv")
        analysis_mode = downstream_analysis_mode(data.get("analysis_mode", data.get("analysis_cdhit_mode")))

        cigar_path = analysis_input_path(outdir, analysis_mode)
        if not cigar_path.exists():
            if analysis_mode == "cdhit98":
                return json_error(
                    "Build the CD-HIT 98.9% sensitivity input in Quality control before running this analysis.",
                    400,
                )
            return json_error(
                "Pipeline results not found. Run the main SIMPLseq pipeline first "
                "(seqtab_cigar.tsv is required).",
                400,
            )
        prepared_samples = analysis_samples_path(outdir, analysis_mode)
        if not prepared_samples.exists():
            return json_error("Build the shared analysis input in Quality control before running DINEMITES.", 400)

        _, active = reconcile_state(dinemites_state_path(outdir, analysis_mode))
        if active:
            return json_error("A DINEMITES analysis is already running for this output folder.", 409)

        worker_pid = start_dinemites_job(
            app_root, outdir, samples, model_type, model_settings, analysis_mode
        )

        return jsonify({
            "ok": True,
            "status": "running",
            "worker_pid": worker_pid,
            "model": model_type,
            "analysis_mode": analysis_mode,
            **model_settings,
            "outdir": str(outdir),
        })

    @app.get("/api/dinemites/status")
    def api_dinemites_status():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        analysis_mode = downstream_analysis_mode(request.args.get("mode"))
        state, active = reconcile_state(dinemites_state_path(outdir, analysis_mode))
        status = state.get("status", "idle")
        if active and status not in {"running"}:
            status = "running"
        return jsonify({
            "ok": True,
            "status": status,
            "active": active,
            "state": state,
            "outdir": str(outdir),
            "analysis_mode": analysis_mode,
        })

    @app.get("/api/dinemites/results")
    def api_dinemites_results():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        analysis_mode = downstream_analysis_mode(request.args.get("mode"))
        dm_dir = dinemites_outdir(outdir, analysis_mode)
        state, _ = reconcile_state(dinemites_state_path(outdir, analysis_mode))

        files: dict[str, dict[str, Any]] = {}
        for key, filename in DINEMITES_FILE_KEYS.items():
            path = dm_dir / filename
            exists = safe_is_file(path)
            size_bytes = path.stat().st_size if exists else 0
            files[key] = {
                "exists": exists,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "download_url": url_for(
                    "download_dinemites", file_key=key, out=str(outdir), mode=analysis_mode
                ) if exists else "",
                "table_url": url_for(
                    "api_dinemites_table", file_key=key, out=str(outdir), mode=analysis_mode
                )
                if exists
                else "",
            }

        plots = dinemites_plot_entries(outdir, analysis_mode)
        model_file_keys = set(DINEMITES_FILE_KEYS) - {"input"}
        has_outputs = any(files[key]["exists"] for key in model_file_keys) or bool(plots)
        if state.get("status") in {None, "", "idle"} and has_outputs:
            state = {
                **state,
                "status": "complete",
                "detail": "Existing DINEMITES outputs found.",
                "outdir": str(dm_dir),
            }

        # Read summary from results files if available
        summary: dict[str, Any] = {"new_infections": "--", "molfoi": "--", "subjects": "--"}
        model_summary = read_json(dm_dir / "dinemites_summary.json")
        preparation = read_json(dm_dir / "dinemites_preparation.json")
        parameters = read_json(outdir / "parameters.json")
        analysis_parameters = parameters.get("analysis_parameters", {})
        preparation["pipeline_qc"] = {
            "min_total_reads": analysis_parameters.get(
                "cigar_min_total_reads", parameters.get("cigar_min_total_reads", 100)
            ),
            "min_sequencing_samples": analysis_parameters.get(
                "cigar_min_samples", parameters.get("cigar_min_samples", 2)
            ),
            "exclude_bimeras": analysis_parameters.get(
                "cigar_exclude_bimeras", parameters.get("cigar_exclude_bimeras", True)
            ),
        }
        preparation.setdefault("abundance_filter", {
            "threshold_percent": state.get("min_abundance_pct", 1.0),
            "denominator": state.get("abundance_denominator", "locus"),
        })
        preparation.setdefault("replicate_merge", {"rule": "intersection"})
        preparation.setdefault(
            "submitted",
            dinemites_input_summary(dm_dir / DINEMITES_FILE_KEYS["input"]),
        )
        new_inf_path = dm_dir / "dinemites_new_infections.tsv"
        molfoi_path = dm_dir / "dinemites_molfoi.tsv"
        subjects_data: list[dict[str, str]] = []
        molfoi_by_subject: dict[str, str] = {}
        allele_key_rows: list[dict[str, str]] = []

        for row_item in dinemites_tsv_rows(dm_dir / "dinemites_allele_key.tsv", limit=500):
            allele_key_rows.append({
                "short_allele_id": dinemites_row_value(row_item, "short_allele_id", "short_id", "allele_id"),
                "locus": dinemites_row_value(row_item, "locus"),
                "allele": dinemites_row_value(row_item, "allele"),
            })

        if safe_is_file(molfoi_path):
            try:
                with molfoi_path.open(newline="", encoding="utf-8-sig", errors="replace") as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    rows = list(reader)
                    vals = []
                    for r in rows:
                        subject = dinemites_row_value(r, "subject", "participant_id")
                        molfoi_value = dinemites_row_value(r, "molFOI", "molfoi", "mol_foi")
                        if subject:
                            molfoi_by_subject[subject] = molfoi_value
                        try:
                            vals.append(float(molfoi_value))
                        except (ValueError, TypeError):
                            pass
                    if vals:
                        summary["molfoi"] = round(sum(vals) / len(vals), 3)
                    if rows:
                        summary["subjects"] = len({dinemites_row_value(r, "subject", "participant_id") for r in rows})
            except Exception:
                pass

        if safe_is_file(new_inf_path):
            try:
                with new_inf_path.open(newline="", encoding="utf-8-sig", errors="replace") as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    rows = list(reader)
                    new_infection_values = []
                    summary["subjects"] = len({dinemites_row_value(r, "subject", "participant_id") for r in rows})
                    for r in rows:
                        new_value = dinemites_row_value(r, "new_infections", "n_new")
                        try:
                            new_infection_values.append(float(new_value))
                        except (ValueError, TypeError):
                            pass
                    for r in rows[:200]:
                        subject = dinemites_row_value(r, "subject", "participant_id")
                        new_value = dinemites_row_value(r, "new_infections", "n_new")
                        subjects_data.append({
                            "subject": subject,
                            "new_infections": new_value,
                            "molfoi": molfoi_by_subject.get(subject, dinemites_row_value(r, "molFOI", "molfoi", "mol_foi")),
                            "time_points": dinemites_row_value(r, "time_points", "n_timepoints"),
                        })
                    if new_infection_values:
                        total_new = sum(new_infection_values)
                        summary["new_infections"] = int(total_new) if total_new.is_integer() else round(total_new, 3)
            except Exception:
                pass

        if not subjects_data and molfoi_by_subject:
            for subject, molfoi_value in list(molfoi_by_subject.items())[:200]:
                subjects_data.append({
                    "subject": subject,
                    "new_infections": "",
                    "molfoi": molfoi_value,
                    "time_points": "",
                })

        return jsonify({
            "ok": True,
            "state": state,
            "files": files,
            "plots": plots,
            "summary": summary,
            "model_summary": model_summary,
            "preparation": preparation,
            "readiness": analysis_readiness(outdir, analysis_mode),
            "subjects": subjects_data,
            "allele_key": allele_key_rows,
            "outdir": str(outdir),
            "analysis_mode": analysis_mode,
            "output_dir": str(dm_dir),
        })

    @app.get("/api/dinemites/table/<file_key>")
    def api_dinemites_table(file_key: str):
        if file_key not in DINEMITES_FILE_KEYS:
            abort(404)
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        analysis_mode = downstream_analysis_mode(request.args.get("mode"))
        dm_dir = dinemites_outdir(outdir, analysis_mode).resolve()
        path = (dm_dir / DINEMITES_FILE_KEYS[file_key]).resolve()
        try:
            path.relative_to(dm_dir)
        except ValueError:
            abort(404)
        if not safe_is_file(path):
            abort(404)
        limit = max(1, min(int_payload(request.args, "limit", 500), 1000))
        rows: list[dict[str, str]] = []
        truncated = False
        with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = list(reader.fieldnames or [])
            for index, row_item in enumerate(reader):
                if index >= limit:
                    truncated = True
                    break
                rows.append({column: row_item.get(column, "") for column in columns})
        return jsonify({
            "ok": True,
            "label": DINEMITES_FILE_LABELS.get(file_key, "DINEMITES result table"),
            "filename": path.name,
            "size": human_bytes(path.stat().st_size),
            "columns": columns,
            "rows": rows,
            "shown_rows": len(rows),
            "truncated": truncated,
        })

    @app.get("/download/dinemites/<file_key>")
    def download_dinemites(file_key: str):
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        analysis_mode = downstream_analysis_mode(request.args.get("mode"))
        filename = DINEMITES_FILE_KEYS.get(file_key)
        if not filename:
            abort(404)
        dm_dir = dinemites_outdir(outdir, analysis_mode).resolve()
        path = (dm_dir / filename).resolve()
        try:
            path.relative_to(dm_dir)
        except ValueError:
            abort(404)
        if not path.exists() or not path.is_file():
            abort(404)
        return send_file(path, as_attachment=True, download_name=path.name)

    @app.get("/download/dinemites-plot/<path:filename>")
    def download_dinemites_plot(filename: str):
        if Path(filename).name != filename or "\\" in filename:
            abort(404)
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        analysis_mode = downstream_analysis_mode(request.args.get("mode"))
        plots_dir = dinemites_plots_dir(outdir, analysis_mode).resolve()
        path = (plots_dir / filename).resolve()
        try:
            path.relative_to(plots_dir)
        except ValueError:
            abort(404)
        if path.suffix.lower() != ".png" or not path.exists() or not path.is_file():
            abort(404)
        inline = request.args.get("inline") == "1"
        return send_file(path, as_attachment=not inline, download_name=path.name, mimetype="image/png")

    # ------------------------------------------------------------------
    # dcifer analysis endpoints
    # ------------------------------------------------------------------

    DCIFER_FILE_KEYS = {
        "input": "dcifer_input_long.tsv",
        "filter_summary": "dcifer_filter_summary.tsv",
        "replicate_summary": "dcifer_replicate_summary.tsv",
        "coi": "dcifer_coi.tsv",
        "allele_frequencies": "dcifer_allele_frequencies.tsv",
        "pairwise_relatedness": "dcifer_pairwise_relatedness.tsv",
        "relatedness_matrix": "dcifer_relatedness_matrix.tsv",
        "pvalue_matrix": "dcifer_pvalue_matrix.tsv",
        "summary": "dcifer_summary.json",
    }

    def dcifer_outdir(outdir: Path, mode: str = "primary") -> Path:
        normalized = downstream_analysis_mode(mode)
        return outdir / ("dcifer_cdhit98" if normalized == "cdhit98" else "dcifer")

    def dcifer_plots_dir(outdir: Path, mode: str = "primary") -> Path:
        return dcifer_outdir(outdir, mode) / "dcifer_plots"

    def dcifer_state_path(outdir: Path, mode: str = "primary") -> Path:
        return dcifer_outdir(outdir, mode) / "dcifer_state.json"

    def clear_dcifer_outputs(dc_dir: Path) -> None:
        for filename in DCIFER_FILE_KEYS.values():
            path = dc_dir / filename
            try:
                if path.is_file():
                    path.unlink()
            except OSError:
                pass
        plots_dir = dc_dir / "dcifer_plots"
        try:
            plot_paths = list(plots_dir.glob("*.png")) if plots_dir.is_dir() else []
        except OSError:
            plot_paths = []
        for path in plot_paths:
            try:
                path.unlink()
            except OSError:
                pass

    def dcifer_plot_entries(outdir: Path, mode: str = "primary") -> list[dict[str, Any]]:
        plots_dir = dcifer_plots_dir(outdir, mode)
        if not safe_is_dir(plots_dir):
            return []
        entries: list[dict[str, Any]] = []
        try:
            paths = sorted(plots_dir.glob("*.png"), key=lambda item: item.name.lower())
        except OSError:
            return []
        for path in paths:
            if not safe_is_file(path):
                continue
            size_bytes = path.stat().st_size
            filename = path.name
            entries.append({
                "filename": filename,
                "title": Path(filename).stem.replace("_", " "),
                "exists": True,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "view_url": url_for(
                    "download_dcifer_plot",
                    filename=filename,
                    out=str(outdir),
                    mode=downstream_analysis_mode(mode),
                    inline=1,
                ),
                "download_url": url_for(
                    "download_dcifer_plot",
                    filename=filename,
                    out=str(outdir),
                    mode=downstream_analysis_mode(mode),
                ),
            })
        return entries

    def dcifer_matrix_payload(path: Path, value_label: str, limit: int = 80) -> dict[str, Any]:
        if not safe_is_file(path):
            return {"labels": [], "rows": [], "value_label": value_label, "truncated": False}
        try:
            with path.open(newline="", encoding="utf-8-sig", errors="replace") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                fieldnames = reader.fieldnames or []
                row_label_field = str(fieldnames[0] or "sample_id") if fieldnames else "sample_id"
                all_labels = fieldnames[1:]
                labels = all_labels[:limit]
                rows: list[dict[str, Any]] = []
                row_count_seen = 0
                truncated = len(all_labels) > limit
                for index, row_item in enumerate(reader):
                    row_count_seen = index + 1
                    if index >= limit:
                        truncated = True
                        break
                    values: list[float | None] = []
                    for label in labels:
                        raw_value = (row_item.get(label) or "").strip()
                        if raw_value.upper() in {"", "NA", "NAN"}:
                            values.append(None)
                            continue
                        try:
                            values.append(float(raw_value))
                        except ValueError:
                            values.append(None)
                    row_label = dinemites_row_value(row_item, "sample_id", row_label_field)
                    rows.append({"sample_id": row_label, "values": values})
                return {
                    "labels": labels,
                    "rows": rows,
                    "value_label": value_label,
                    "total_columns": len(all_labels),
                    "total_rows": row_count_seen,
                    "truncated": truncated,
                }
        except Exception:
            return {"labels": [], "rows": [], "value_label": value_label, "truncated": False}

    def start_dcifer_job(
        root: Path,
        outdir: Path,
        samples: Path,
        settings: dict[str, int | float | str],
        analysis_mode: str = "primary",
    ) -> int:
        analysis_mode = downstream_analysis_mode(analysis_mode)
        dc_dir = dcifer_outdir(outdir, analysis_mode)
        dc_dir.mkdir(parents=True, exist_ok=True)
        clear_dcifer_outputs(dc_dir)
        state_file = dcifer_state_path(outdir, analysis_mode)
        cigar_path = analysis_input_path(outdir, analysis_mode)
        prepared_samples = analysis_samples_path(outdir, analysis_mode)
        dc_input = dc_dir / "dcifer_input_long.tsv"
        log_path = dc_dir / "dcifer.log"
        cmd_convert = [
            "Rscript",
            str(root / "workflow" / "scripts" / "simplseq_to_dcifer.R"),
            "--cigar", str(cigar_path),
            "--samples", str(prepared_samples),
            "--out", str(dc_input),
            "--filter_summary", str(dc_dir / "dcifer_filter_summary.tsv"),
            "--replicate_summary", str(dc_dir / "dcifer_replicate_summary.tsv"),
            "--min_abundance_pct", "0",
            "--abundance_denominator", "locus",
        ]
        cmd_run = [
            "Rscript",
            str(root / "workflow" / "scripts" / "run_dcifer.R"),
            "--input", str(dc_input),
            "--outdir", str(dc_dir),
            "--coi_lrank", str(settings["coi_lrank"]),
            "--ibd_grid_nr", str(settings["ibd_grid_nr"]),
            "--alpha", str(settings["alpha"]),
            "--afreq_mode", str(settings["afreq_mode"]),
        ]
        return launch_job(
            job_type="dcifer",
            root=root,
            output_dir=dc_dir,
            state_file=state_file,
            log_file=log_path,
            state_payload={**settings, "analysis_mode": analysis_mode, "outdir": str(dc_dir)},
            steps=[
                {
                    "label": "Preparing Dcifer input",
                    "command": cmd_convert,
                    "timeout_seconds": 600,
                    "failure_detail": "Could not prepare the Dcifer input. Check dcifer.log.",
                },
                {
                    "label": "Estimating relatedness",
                    "command": cmd_run,
                    "timeout_seconds": 7200,
                    "failure_detail": "Dcifer did not complete. Check dcifer.log.",
                },
            ],
            env=analysis_runtime_env(root),
        )

    @app.post("/api/dcifer/run")
    def api_dcifer_run():
        data = request.get_json(silent=True) or {}
        try:
            settings = dcifer_settings(data)
        except ValueError as exc:
            return json_error(str(exc), 400)

        outdir = resolve_app_path(workspace, data.get("outdir"), "results")
        samples = resolve_app_path(workspace, data.get("samples"), "samples.csv")
        analysis_mode = downstream_analysis_mode(data.get("analysis_mode", data.get("analysis_cdhit_mode")))
        cigar_path = analysis_input_path(outdir, analysis_mode)
        if not cigar_path.exists():
            if analysis_mode == "cdhit98":
                return json_error(
                    "Build the CD-HIT 98.9% sensitivity input in Quality control before running this analysis.",
                    400,
                )
            return json_error(
                "Pipeline results not found. Run the main SIMPLseq pipeline first "
                "(run_dada2/seqtab_cigar.tsv is required).",
                400,
            )
        prepared_samples = analysis_samples_path(outdir, analysis_mode)
        if not prepared_samples.exists():
            return json_error("Build the shared analysis input in Quality control before running Dcifer.", 400)
        _, active = reconcile_state(dcifer_state_path(outdir, analysis_mode))
        if active:
            return json_error("A dcifer analysis is already running for this output folder.", 409)

        worker_pid = start_dcifer_job(app_root, outdir, samples, settings, analysis_mode)

        return jsonify({
            "ok": True,
            "status": "running",
            "worker_pid": worker_pid,
            "analysis_mode": analysis_mode,
            **settings,
            "outdir": str(outdir),
        })

    @app.get("/api/dcifer/status")
    def api_dcifer_status():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        analysis_mode = downstream_analysis_mode(request.args.get("mode"))
        state, active = reconcile_state(dcifer_state_path(outdir, analysis_mode))
        status = state.get("status", "idle")
        if active and status != "running":
            status = "running"
        return jsonify({
            "ok": True,
            "status": status,
            "active": active,
            "state": state,
            "outdir": str(outdir),
            "analysis_mode": analysis_mode,
        })

    @app.get("/api/dcifer/results")
    def api_dcifer_results():
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        analysis_mode = downstream_analysis_mode(request.args.get("mode"))
        dc_dir = dcifer_outdir(outdir, analysis_mode)
        state, _ = reconcile_state(dcifer_state_path(outdir, analysis_mode))

        files: dict[str, dict[str, Any]] = {}
        for key, filename in DCIFER_FILE_KEYS.items():
            path = dc_dir / filename
            exists = safe_is_file(path)
            size_bytes = path.stat().st_size if exists else 0
            files[key] = {
                "exists": exists,
                "size_bytes": size_bytes,
                "size": human_bytes(size_bytes),
                "download_url": url_for(
                    "download_dcifer", file_key=key, out=str(outdir), mode=analysis_mode
                ) if exists else "",
            }

        plots = dcifer_plot_entries(outdir, analysis_mode)
        has_outputs = any(item["exists"] for item in files.values()) or bool(plots)
        if state.get("status") in {None, "", "idle"} and has_outputs:
            state = {
                **state,
                "status": "complete",
                "detail": "Existing dcifer outputs found.",
                "outdir": str(dc_dir),
            }

        summary: dict[str, Any] = {
            "samples": "--",
            "pairs": "--",
            "max_relatedness": "--",
            "raw_p_le_alpha": "--",
            "q_le_alpha": "--",
        }
        summary_json = read_json(dc_dir / "dcifer_summary.json")
        for key in ("samples", "loci", "pairs", "max_relatedness", "raw_p_le_alpha", "q_le_alpha", "adequacy", "caveat"):
            if key in summary_json:
                summary[key] = summary_json[key]

        coi_rows = dinemites_tsv_rows(dc_dir / "dcifer_coi.tsv", limit=10000)
        if coi_rows and "samples" not in summary_json:
            summary["samples"] = len({
                dinemites_row_value(row_item, "sample_id")
                for row_item in coi_rows
                if dinemites_row_value(row_item, "sample_id")
            })

        pair_rows = dinemites_tsv_rows(dc_dir / "dcifer_pairwise_relatedness.tsv", limit=10000)
        pairs: list[dict[str, str]] = []
        estimates: list[float] = []
        raw_p_count = 0
        alpha = float(state.get("alpha", summary_json.get("alpha", 0.05)) or 0.05)
        for index, row_item in enumerate(pair_rows):
            if index < 200:
                pairs.append({
                    "sample_a": dinemites_row_value(row_item, "sample_a"),
                    "sample_b": dinemites_row_value(row_item, "sample_b"),
                    "estimate": dinemites_row_value(row_item, "estimate"),
                    "p_value": dinemites_row_value(row_item, "p_value"),
                    "q_value": dinemites_row_value(row_item, "q_value"),
                    "ci_lower": dinemites_row_value(row_item, "ci_lower", "CI_lower"),
                    "ci_upper": dinemites_row_value(row_item, "ci_upper", "CI_upper"),
                    "comparison_type": dinemites_row_value(row_item, "comparison_type"),
                })
            try:
                estimates.append(float(dinemites_row_value(row_item, "estimate")))
            except (TypeError, ValueError):
                pass
            raw_flag = dinemites_row_value(row_item, "raw_p_le_alpha").lower()
            if raw_flag in {"true", "t", "1", "yes"}:
                raw_p_count += 1
            elif not raw_flag:
                try:
                    if float(dinemites_row_value(row_item, "p_value")) <= alpha:
                        raw_p_count += 1
                except (TypeError, ValueError):
                    pass
        if pair_rows and "pairs" not in summary_json:
            summary["pairs"] = len(pair_rows)
        if pair_rows and "raw_p_le_alpha" not in summary_json:
            summary["raw_p_le_alpha"] = raw_p_count
        if estimates and "max_relatedness" not in summary_json:
            summary["max_relatedness"] = round(max(estimates), 6)

        return jsonify({
            "ok": True,
            "state": state,
            "files": files,
            "plots": plots,
            "matrices": {
                "relatedness": dcifer_matrix_payload(dc_dir / "dcifer_relatedness_matrix.tsv", "Relatedness"),
                "pvalue": dcifer_matrix_payload(dc_dir / "dcifer_pvalue_matrix.tsv", "p-value"),
            },
            "summary": summary,
            "readiness": analysis_readiness(outdir, analysis_mode),
            "pairs": pairs,
            "outdir": str(outdir),
            "analysis_mode": analysis_mode,
            "output_dir": str(dc_dir),
        })

    @app.get("/download/dcifer/<file_key>")
    def download_dcifer(file_key: str):
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        analysis_mode = downstream_analysis_mode(request.args.get("mode"))
        filename = DCIFER_FILE_KEYS.get(file_key)
        if not filename:
            abort(404)
        dc_dir = dcifer_outdir(outdir, analysis_mode).resolve()
        path = (dc_dir / filename).resolve()
        try:
            path.relative_to(dc_dir)
        except ValueError:
            abort(404)
        if not path.exists() or not path.is_file():
            abort(404)
        return send_file(path, as_attachment=True, download_name=path.name)

    @app.get("/download/dcifer-plot/<path:filename>")
    def download_dcifer_plot(filename: str):
        if Path(filename).name != filename or "\\" in filename:
            abort(404)
        outdir = resolve_app_path(workspace, request.args.get("out"), "results")
        analysis_mode = downstream_analysis_mode(request.args.get("mode"))
        plots_dir = dcifer_plots_dir(outdir, analysis_mode).resolve()
        path = (plots_dir / filename).resolve()
        try:
            path.relative_to(plots_dir)
        except ValueError:
            abort(404)
        if path.suffix.lower() != ".png" or not path.exists() or not path.is_file():
            abort(404)
        inline = request.args.get("inline") == "1"
        return send_file(path, as_attachment=not inline, download_name=path.name, mimetype="image/png")

    return app


def open_browser_later(url: str) -> None:
    time.sleep(1.0)
    if is_wsl():
        try:
            subprocess.Popen(["cmd.exe", "/c", "start", "", url], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return
        except Exception:
            pass
    try:
        webbrowser.open(url)
    except Exception:
        pass


def run_server(root: Path | None = None, host: str = "127.0.0.1", port: int = 8501, open_browser: bool = True) -> int:
    if host not in {"127.0.0.1", "localhost", "::1"}:
        raise ValueError("malaria-amplicon-nf only serves the browser UI on a loopback host")
    app = create_app(root, Path.cwd())
    url = f"http://{host}:{port}"
    if open_browser:
        threading.Thread(target=open_browser_later, args=(url,), daemon=True).start()
    serve(app, host=host, port=port, threads=6, clear_untrusted_proxy_headers=True)
    return 0


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Run the malaria-amplicon-nf Flask GUI.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8501)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()
    return run_server(host=args.host, port=args.port, open_browser=not args.no_browser)


if __name__ == "__main__":
    raise SystemExit(main())
