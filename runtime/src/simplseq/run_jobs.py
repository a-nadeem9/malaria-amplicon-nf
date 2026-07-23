"""Durable detached workers for primary Nextflow pipeline runs."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import signal
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Any

from .job_state import read_json, write_json


SCHEMA_VERSION = 1
ACTIVE_STATUSES = {"starting", "running", "stopping"}
FINAL_STATUSES = {"complete", "failed", "stopped", "dry_run"}
ALLOWED_PROFILES = {"local", "reproducible"}
WORKER_STATE_NAME = "run_job_state.json"
WORKER_SPEC_NAME = ".run-job.json"
MEMORY_RE = re.compile(r"^[1-9][0-9]*(?:\.[0-9]+)?\s*(?:KB|MB|GB|TB)$", re.IGNORECASE)
REQUIRED_SAMPLE_COLUMNS = {"sample_id", "fastq_1", "fastq_2"}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def _resolved_path(value: object, label: str) -> Path:
    text = str(value or "").strip()
    if not text:
        raise ValueError(f"{label} is required.")
    return Path(text).expanduser().resolve()


def _validate_samples(samples: Path) -> None:
    if not samples.is_file():
        raise ValueError(f"Sample sheet is unavailable: {samples}")
    try:
        with samples.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            columns = {str(column or "").strip() for column in (reader.fieldnames or [])}
            missing = sorted(REQUIRED_SAMPLE_COLUMNS - columns)
            if missing:
                raise ValueError(f"Sample sheet is missing required columns: {', '.join(missing)}")
            sample_ids: set[str] = set()
            row_count = 0
            for row_number, row in enumerate(reader, start=2):
                row_count += 1
                sample_id = str(row.get("sample_id", "") or "").strip()
                if not sample_id:
                    raise ValueError(f"Sample sheet row {row_number} has no sample_id.")
                if sample_id in sample_ids:
                    raise ValueError(f"Sample sheet contains duplicate sample_id: {sample_id}")
                sample_ids.add(sample_id)
                for field in ("fastq_1", "fastq_2"):
                    if not str(row.get(field, "") or "").strip():
                        raise ValueError(f"Sample sheet row {row_number} has no {field} path.")
            if row_count == 0:
                raise ValueError("Sample sheet contains no samples.")
    except UnicodeError as exc:
        raise ValueError("Sample sheet is not valid UTF-8 CSV.") from exc


def _validate_options(raw: object, outdir: Path) -> dict[str, Any]:
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ValueError("Run options must be a JSON object.")
    allowed = {
        "resume",
        "work_dir",
        "dry_run",
        "cpus",
        "memory",
        "kelt_enabled",
        "kelt_barcode_map",
    }
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError(f"Unsupported run options: {', '.join(unknown)}")

    resume = raw.get("resume", True)
    dry_run = raw.get("dry_run", False)
    kelt_enabled = raw.get("kelt_enabled", False)
    if not isinstance(resume, bool) or not isinstance(dry_run, bool) or not isinstance(kelt_enabled, bool):
        raise ValueError("resume, dry_run, and kelt_enabled must be true or false.")

    cpus = raw.get("cpus")
    if cpus is not None:
        if isinstance(cpus, bool):
            raise ValueError("cpus must be a positive integer.")
        try:
            cpus = int(cpus)
        except (TypeError, ValueError) as exc:
            raise ValueError("cpus must be a positive integer.") from exc
        if not 1 <= cpus <= 1024:
            raise ValueError("cpus must be between 1 and 1024.")

    memory = raw.get("memory")
    if memory is not None:
        memory = str(memory).strip().upper()
        if not MEMORY_RE.fullmatch(memory):
            raise ValueError("memory must use a value such as '12 GB'.")

    work_dir_value = str(raw.get("work_dir", "") or "").strip()
    work_dir = Path(work_dir_value).expanduser().resolve() if work_dir_value else (outdir / ".nextflow_work").resolve()
    if work_dir.exists() and not work_dir.is_dir():
        raise ValueError("Nextflow work_dir must be a folder.")

    barcode_value = str(raw.get("kelt_barcode_map", "") or "").strip()
    barcode_map = Path(barcode_value).expanduser().resolve() if barcode_value else None
    if barcode_map is not None and not barcode_map.is_file():
        raise ValueError(f"KELT barcode map is unavailable: {barcode_map}")

    return {
        "resume": resume,
        "work_dir": str(work_dir),
        "dry_run": dry_run,
        "cpus": cpus,
        "memory": memory,
        "kelt_enabled": kelt_enabled,
        "kelt_barcode_map": str(barcode_map) if barcode_map else "",
    }


def validate_job_spec(raw: dict[str, Any]) -> dict[str, Any]:
    if raw.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported pipeline run job specification version.")

    root = _resolved_path(raw.get("root"), "Runtime root")
    samples = _resolved_path(raw.get("samples"), "Sample sheet")
    outdir = _resolved_path(raw.get("outdir"), "Output folder")
    state_file = _resolved_path(raw.get("state_file"), "Worker state file")
    spec_path = _resolved_path(raw.get("spec_path"), "Worker specification file")
    profile = str(raw.get("profile", "")).strip().lower()
    token = str(raw.get("worker_token", "")).strip()

    if not (root / "main.nf").is_file():
        raise ValueError("Runtime root does not contain main.nf.")
    _validate_samples(samples)
    if outdir.exists() and not outdir.is_dir():
        raise ValueError("Output path must be a folder.")
    if profile not in ALLOWED_PROFILES:
        raise ValueError("Pipeline profile must be local or reproducible.")
    if not re.fullmatch(r"[0-9a-f]{32}", token):
        raise ValueError("Pipeline worker token is invalid.")
    if not _is_within(state_file, outdir) or not _is_within(spec_path, outdir):
        raise ValueError("Pipeline worker state and specification must stay inside the output folder.")
    if state_file.name.casefold() == "run_state.json":
        raise ValueError("Worker status must not overwrite the pipeline run_state.json file.")

    options = _validate_options(raw.get("options"), outdir)
    return {
        "schema_version": SCHEMA_VERSION,
        "root": str(root),
        "samples": str(samples),
        "outdir": str(outdir),
        "profile": profile,
        "options": options,
        "state_file": str(state_file),
        "spec_path": str(spec_path),
        "worker_token": token,
    }


def _state(spec: dict[str, Any], status: str, **extra: Any) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": status,
        "samples": spec["samples"],
        "outdir": spec["outdir"],
        "profile": spec["profile"],
        "job_spec": spec["spec_path"],
        "worker_token": spec["worker_token"],
    }
    payload.update(extra)
    return payload


def _wait_for_parent_state(state_file: Path, pid: int, token: str, timeout: float = 10.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        state = read_json(state_file)
        if state.get("worker_pid") == pid and state.get("worker_token") == token:
            return
        time.sleep(0.02)
    raise RuntimeError("Pipeline worker could not verify its persisted PID and token.")


def execute_job(spec_path: Path, expected_token: str) -> int:
    spec_path = spec_path.expanduser().resolve()
    raw = read_json(spec_path)
    if str(raw.get("worker_token", "")) != expected_token:
        raise ValueError("Pipeline worker token did not match its job specification.")
    spec = validate_job_spec({**raw, "spec_path": str(spec_path)})
    state_file = Path(spec["state_file"])
    pid = os.getpid()
    _wait_for_parent_state(state_file, pid, expected_token)
    started_at = read_json(state_file).get("started_at") or utc_now()

    write_json(
        state_file,
        _state(spec, "running", started_at=started_at, worker_pid=pid),
    )
    try:
        from .runner import run_nextflow

        options = spec["options"]
        result = run_nextflow(
            Path(spec["samples"]),
            Path(spec["outdir"]),
            profile=spec["profile"],
            resume=options["resume"],
            work_dir=Path(options["work_dir"]),
            root=Path(spec["root"]),
            dry_run=options["dry_run"],
            cpus=options["cpus"],
            memory=options["memory"],
            kelt_enabled=options["kelt_enabled"],
            kelt_barcode_map=Path(options["kelt_barcode_map"]) if options["kelt_barcode_map"] else None,
        )
        returncode = int(result.returncode)
        status = "dry_run" if options["dry_run"] and returncode == 0 else ("complete" if returncode == 0 else "failed")
        write_json(
            state_file,
            _state(
                spec,
                status,
                started_at=started_at,
                completed_at=utc_now(),
                returncode=returncode,
                command=result.command,
                technical_log=str(result.log_path),
            ),
        )
        return returncode
    except Exception as exc:
        write_json(
            state_file,
            _state(
                spec,
                "failed",
                started_at=started_at,
                completed_at=utc_now(),
                returncode=1,
                detail=str(exc),
            ),
        )
        return 1


def launch_job(
    *,
    samples: Path,
    outdir: Path,
    root: Path,
    profile: str = "local",
    resume: bool = True,
    work_dir: Path | None = None,
    dry_run: bool = False,
    cpus: int | None = None,
    memory: str | None = None,
    kelt_enabled: bool = False,
    kelt_barcode_map: Path | None = None,
    env: dict[str, str] | None = None,
) -> int:
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    state_file = outdir / WORKER_STATE_NAME
    previous_state, active = reconcile_state(state_file)
    if active:
        raise RuntimeError(f"A pipeline worker is already active with PID {previous_state.get('worker_pid')}.")

    token = uuid.uuid4().hex
    spec_path = outdir / WORKER_SPEC_NAME
    spec = validate_job_spec(
        {
            "schema_version": SCHEMA_VERSION,
            "root": str(root.expanduser().resolve()),
            "samples": str(samples.expanduser().resolve()),
            "outdir": str(outdir),
            "profile": profile,
            "options": {
                "resume": resume,
                "work_dir": str(work_dir.expanduser().resolve()) if work_dir else "",
                "dry_run": dry_run,
                "cpus": cpus,
                "memory": memory,
                "kelt_enabled": kelt_enabled,
                "kelt_barcode_map": str(kelt_barcode_map.expanduser().resolve()) if kelt_barcode_map else "",
            },
            "state_file": str(state_file),
            "spec_path": str(spec_path),
            "worker_token": token,
        }
    )
    write_json(spec_path, spec)
    started_at = utc_now()
    write_json(state_file, _state(spec, "starting", started_at=started_at))

    command = [sys.executable, "-m", "simplseq.run_jobs", "--spec", str(spec_path), "--token", token]
    kwargs: dict[str, Any] = {
        "cwd": Path(spec["root"]),
        "env": env or os.environ.copy(),
        "stdin": subprocess.DEVNULL,
        "stdout": subprocess.DEVNULL,
        "stderr": subprocess.DEVNULL,
        "text": True,
        "close_fds": True,
    }
    if os.name == "nt":
        kwargs["creationflags"] = (
            getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
            | getattr(subprocess, "DETACHED_PROCESS", 0)
            | getattr(subprocess, "CREATE_NO_WINDOW", 0)
        )
    else:
        kwargs["start_new_session"] = True

    try:
        process = subprocess.Popen(command, **kwargs)
    except Exception as exc:
        write_json(
            state_file,
            _state(spec, "failed", started_at=started_at, completed_at=utc_now(), detail=str(exc)),
        )
        raise
    write_json(
        state_file,
        _state(spec, "starting", started_at=started_at, worker_pid=process.pid),
    )
    return process.pid


def _process_command(pid: int) -> str:
    proc_cmdline = Path(f"/proc/{pid}/cmdline")
    if proc_cmdline.is_file():
        try:
            return proc_cmdline.read_bytes().replace(b"\0", b" ").decode(errors="replace")
        except OSError:
            return ""
    if os.name != "nt":
        try:
            result = subprocess.run(
                ["ps", "-p", str(pid), "-o", "command="],
                capture_output=True,
                text=True,
                check=False,
                timeout=2,
            )
            return result.stdout.strip() if result.returncode == 0 else ""
        except (OSError, subprocess.TimeoutExpired):
            return ""
    try:
        result = subprocess.run(
            [
                "powershell.exe",
                "-NoLogo",
                "-NoProfile",
                "-Command",
                f"(Get-CimInstance Win32_Process -Filter 'ProcessId = {pid}').CommandLine",
            ],
            capture_output=True,
            text=True,
            check=False,
            timeout=3,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        return result.stdout.strip() if result.returncode == 0 else ""
    except (OSError, subprocess.TimeoutExpired):
        return ""


def _normalized_command(value: object) -> str:
    return str(value or "").replace("\\", "/").casefold()


def job_is_active(state: dict[str, Any]) -> bool:
    if state.get("status") not in ACTIVE_STATUSES:
        return False
    try:
        pid = int(state.get("worker_pid", 0))
    except (TypeError, ValueError):
        return False
    if pid <= 0 or pid == os.getpid():
        return False
    command = _process_command(pid)
    normalized = _normalized_command(command)
    spec_path = _normalized_command(state.get("job_spec"))
    token = str(state.get("worker_token", ""))
    return bool(
        command
        and "simplseq.run_jobs" in normalized
        and spec_path
        and spec_path in normalized
        and re.fullmatch(r"[0-9a-f]{32}", token)
        and token in command
    )


def reconcile_state(state_file: Path) -> tuple[dict[str, Any], bool]:
    state = read_json(state_file)
    active = job_is_active(state)
    if state.get("status") in ACTIVE_STATUSES and not active:
        state = {
            **state,
            "status": "failed",
            "completed_at": utc_now(),
            "returncode": 1,
            "detail": "The pipeline worker stopped before completion. Pipeline run_state.json was left unchanged for diagnosis.",
        }
        state.pop("worker_pid", None)
        write_json(state_file, state)
    return state, active


def _pid_exists(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


def _terminate_process_tree(state: dict[str, Any], timeout_seconds: float = 5.0) -> None:
    if not job_is_active(state):
        raise RuntimeError("Refusing to stop a process that is not the verified pipeline worker.")
    pid = int(state["worker_pid"])
    if os.name == "nt":
        completed = subprocess.run(
            ["taskkill", "/PID", str(pid), "/T", "/F"],
            capture_output=True,
            text=True,
            check=False,
            timeout=max(1.0, timeout_seconds),
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        if completed.returncode != 0 and _process_command(pid):
            detail = (completed.stderr or completed.stdout).strip()
            raise RuntimeError(detail or "Windows could not stop the pipeline worker tree.")
        return

    try:
        process_group = os.getpgid(pid)
    except ProcessLookupError:
        return
    if process_group != pid:
        raise RuntimeError("Refusing to stop a pipeline worker outside its detached process group.")
    os.killpg(process_group, signal.SIGTERM)
    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        if not _pid_exists(pid):
            return
        time.sleep(0.05)
    if _pid_exists(pid):
        os.killpg(process_group, signal.SIGKILL)


def stop_job(state_file: Path, timeout_seconds: float = 5.0) -> bool:
    state_file = state_file.expanduser().resolve()
    state = read_json(state_file)
    if state.get("status") not in ACTIVE_STATUSES:
        return False
    if not job_is_active(state):
        reconcile_state(state_file)
        return False

    stopping = {**state, "status": "stopping", "stop_requested_at": utc_now()}
    write_json(state_file, stopping)
    _terminate_process_tree(stopping, timeout_seconds=timeout_seconds)
    stopped = {
        **stopping,
        "status": "stopped",
        "completed_at": utc_now(),
        "detail": "Pipeline run stopped by the user.",
    }
    stopped.pop("worker_pid", None)
    write_json(state_file, stopped)
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run one detached malaria-amplicon-nf pipeline job.")
    parser.add_argument("--spec", required=True, type=Path)
    parser.add_argument("--token", required=True)
    args = parser.parse_args(argv)
    return execute_job(args.spec, args.token)


if __name__ == "__main__":
    raise SystemExit(main())
