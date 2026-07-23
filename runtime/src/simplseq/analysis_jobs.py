"""Durable subprocess jobs for downstream scientific analyses."""

from __future__ import annotations

import argparse
import datetime as dt
import os
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Any

from .job_state import process_command, read_json, write_json


SCHEMA_VERSION = 1
ACTIVE_STATUSES = {"starting", "running"}
ALLOWED_SCRIPTS = {
    "dinemites": {"simplseq_to_dinemites.R", "run_dinemites.R"},
    "dcifer": {"simplseq_to_dcifer.R", "run_dcifer.R"},
}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def validate_job_spec(raw: dict[str, Any]) -> dict[str, Any]:
    if raw.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported analysis job specification version.")

    job_type = str(raw.get("job_type", "")).strip().lower()
    if job_type not in ALLOWED_SCRIPTS:
        raise ValueError("Unsupported analysis job type.")

    root = Path(str(raw.get("root", ""))).expanduser().resolve()
    output_dir = Path(str(raw.get("output_dir", ""))).expanduser().resolve()
    state_file = Path(str(raw.get("state_file", ""))).expanduser().resolve()
    log_file = Path(str(raw.get("log_file", ""))).expanduser().resolve()
    scripts_dir = (root / "workflow" / "scripts").resolve()

    if not scripts_dir.is_dir():
        raise ValueError("Workflow scripts directory is unavailable.")
    if not _is_within(state_file, output_dir) or not _is_within(log_file, output_dir):
        raise ValueError("Analysis state and log files must stay inside the analysis output folder.")

    steps = raw.get("steps")
    if not isinstance(steps, list) or not 1 <= len(steps) <= 4:
        raise ValueError("Analysis job must contain between one and four steps.")

    validated_steps: list[dict[str, Any]] = []
    for step in steps:
        if not isinstance(step, dict):
            raise ValueError("Invalid analysis job step.")
        command = step.get("command")
        if not isinstance(command, list) or len(command) < 2 or not all(isinstance(item, str) for item in command):
            raise ValueError("Analysis commands must be argument lists.")
        if Path(command[0]).name.lower() not in {"rscript", "rscript.exe"}:
            raise ValueError("Only managed Rscript analysis steps are allowed.")
        script = Path(command[1]).expanduser().resolve()
        if script.parent != scripts_dir or script.name not in ALLOWED_SCRIPTS[job_type]:
            raise ValueError(f"Unapproved {job_type} analysis script: {script.name}")
        timeout = int(step.get("timeout_seconds", 0))
        if timeout < 1 or timeout > 86_400:
            raise ValueError("Analysis step timeout must be between 1 second and 24 hours.")
        validated_steps.append({
            "label": str(step.get("label", script.stem)).strip() or script.stem,
            "command": command,
            "timeout_seconds": timeout,
            "failure_detail": str(step.get("failure_detail", "Analysis step failed.")).strip(),
        })

    state_payload = raw.get("state_payload", {})
    if not isinstance(state_payload, dict):
        raise ValueError("Invalid analysis state payload.")

    return {
        **raw,
        "job_type": job_type,
        "root": str(root),
        "output_dir": str(output_dir),
        "state_file": str(state_file),
        "log_file": str(log_file),
        "steps": validated_steps,
        "state_payload": state_payload,
    }


def _state(spec: dict[str, Any], status: str, **extra: Any) -> dict[str, Any]:
    payload = {
        **spec["state_payload"],
        "status": status,
        "job_type": spec["job_type"],
        "output_dir": spec["output_dir"],
        "job_spec": spec["spec_path"],
        "worker_token": spec["worker_token"],
    }
    payload.update(extra)
    return payload


def _log_tail(path: Path, max_chars: int = 1800) -> str:
    try:
        text = path.read_text(encoding="utf-8", errors="replace").strip()
    except OSError:
        return ""
    return text[-max_chars:]


def _wait_for_parent_state(state_file: Path, pid: int, timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if read_json(state_file).get("worker_pid") == pid:
            return
        time.sleep(0.02)


def execute_job(spec_path: Path, expected_token: str) -> int:
    raw = read_json(spec_path)
    if str(raw.get("worker_token", "")) != expected_token:
        raise ValueError("Analysis worker token did not match its job specification.")
    spec = validate_job_spec({**raw, "spec_path": str(spec_path.resolve())})
    state_file = Path(spec["state_file"])
    log_file = Path(spec["log_file"])
    root = Path(spec["root"])
    pid = os.getpid()
    _wait_for_parent_state(state_file, pid)

    started_at = read_json(state_file).get("started_at") or utc_now()
    try:
        for index, step in enumerate(spec["steps"], start=1):
            write_json(state_file, _state(
                spec,
                "running",
                started_at=started_at,
                worker_pid=pid,
                stage_index=index,
                stage_count=len(spec["steps"]),
                stage=step["label"],
            ))
            mode = "w" if index == 1 else "a"
            with log_file.open(mode, encoding="utf-8") as log_handle:
                log_handle.write(f"\n== {step['label']} ==\n")
                log_handle.flush()
                completed = subprocess.run(
                    step["command"],
                    cwd=root,
                    env=os.environ.copy(),
                    text=True,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=step["timeout_seconds"],
                )
            if completed.returncode != 0:
                detail = step["failure_detail"]
                tail = _log_tail(log_file)
                if tail:
                    detail = f"{detail}\n\n{tail}"
                write_json(state_file, _state(
                    spec,
                    "failed",
                    started_at=started_at,
                    completed_at=utc_now(),
                    detail=detail,
                    failed_stage=step["label"],
                    returncode=completed.returncode,
                ))
                return completed.returncode or 1

        write_json(state_file, _state(
            spec,
            "complete",
            started_at=started_at,
            completed_at=utc_now(),
        ))
        return 0
    except subprocess.TimeoutExpired as exc:
        write_json(state_file, _state(
            spec,
            "failed",
            started_at=started_at,
            completed_at=utc_now(),
            detail=f"{spec['job_type'].upper()} exceeded its {int(exc.timeout)} second step limit.",
        ))
        return 124
    except Exception as exc:
        write_json(state_file, _state(
            spec,
            "failed",
            started_at=started_at,
            completed_at=utc_now(),
            detail=str(exc),
        ))
        return 1


def launch_job(
    *,
    job_type: str,
    root: Path,
    output_dir: Path,
    state_file: Path,
    log_file: Path,
    state_payload: dict[str, Any],
    steps: list[dict[str, Any]],
    env: dict[str, str] | None = None,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    spec_path = output_dir / f".{job_type}-job.json"
    spec = validate_job_spec({
        "schema_version": SCHEMA_VERSION,
        "job_type": job_type,
        "root": str(root.resolve()),
        "output_dir": str(output_dir.resolve()),
        "state_file": str(state_file.resolve()),
        "log_file": str(log_file.resolve()),
        "state_payload": state_payload,
        "steps": steps,
        "worker_token": token,
        "spec_path": str(spec_path.resolve()),
    })
    write_json(spec_path, spec)
    started_at = utc_now()
    write_json(state_file, _state(spec, "starting", started_at=started_at))

    command = [sys.executable, "-m", "simplseq.analysis_jobs", "--spec", str(spec_path), "--token", token]
    kwargs: dict[str, Any] = {
        "cwd": root,
        "env": env or os.environ.copy(),
        "stdin": subprocess.DEVNULL,
        "stdout": subprocess.DEVNULL,
        "stderr": subprocess.DEVNULL,
        "text": True,
    }
    if os.name == "nt":
        kwargs["creationflags"] = (
            getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
            | getattr(subprocess, "CREATE_NO_WINDOW", 0)
        )
    else:
        kwargs["start_new_session"] = True
    process = subprocess.Popen(command, **kwargs)
    write_json(state_file, _state(
        spec,
        "starting",
        started_at=started_at,
        worker_pid=process.pid,
    ))
    return process.pid


def job_is_active(state: dict[str, Any]) -> bool:
    if state.get("status") not in ACTIVE_STATUSES:
        return False
    try:
        pid = int(state.get("worker_pid", 0))
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False
    command = process_command(pid)
    if not command or "simplseq.analysis_jobs" not in command:
        return False
    spec_path = str(state.get("job_spec", ""))
    token = str(state.get("worker_token", ""))
    return bool(spec_path and token and spec_path in command and token in command)


def reconcile_state(state_file: Path) -> tuple[dict[str, Any], bool]:
    state = read_json(state_file)
    active = job_is_active(state)
    if state.get("status") in ACTIVE_STATUSES and not active:
        state = {
            **state,
            "status": "failed",
            "completed_at": utc_now(),
            "detail": "The analysis worker stopped before completion. Run the analysis again; existing pipeline outputs were not changed.",
        }
        state.pop("worker_pid", None)
        write_json(state_file, state)
    return state, active


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run one managed downstream analysis job.")
    parser.add_argument("--spec", required=True, type=Path)
    parser.add_argument("--token", required=True)
    args = parser.parse_args(argv)
    return execute_job(args.spec.expanduser().resolve(), args.token)


if __name__ == "__main__":
    raise SystemExit(main())
