"""600s public-only jobs; terminate only own process groups and retain progress."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
from queue import SimpleQueue
import signal
import subprocess
import sys
from time import monotonic, sleep

from common import ROOT, read, save, atomic, sha, utc


def group_rss(group):
    """Sample sum of process-group RSS, including external flatter children."""
    total = 0
    for path in Path("/proc").iterdir():
        if not path.name.isdigit():
            continue
        try:
            stat = (path / "stat").read_text()
            fields = stat[stat.rfind(")") + 2:].split()
            if int(fields[2]) == group:
                total += int(fields[21]) * os.sysconf("SC_PAGE_SIZE")
        except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError):
            pass
    return total


def stop_group(process):
    """Caller owns this start_new_session process; never target broad users."""
    group = process.pid
    try:
        os.killpg(group, signal.SIGTERM)
    except ProcessLookupError:
        return
    until = monotonic() + 2
    while monotonic() < until:
        try:
            os.killpg(group, 0)
        except ProcessLookupError:
            break
        sleep(.05)
    try:
        os.killpg(group, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def final_status(interruption, returncode, progress, output, timeout, campaign_deadline):
    # A polling delay cannot turn an on-time completed output into a timeout.
    if returncode == 0 and output and progress.get("complete"):
        if progress["elapsed_seconds"] <= timeout and progress["monotonic"] <= campaign_deadline:
            return "COMPLETE"
        return interruption or "LATE_COMPLETE"
    return interruption or "ERROR"


def execute(record, run, cpu, manifest, campaign_deadline):
    job = run / "jobs" / record["id"]
    job.mkdir(exist_ok=False)
    # Only the mathematical public instance is passed to the solver.
    save(job / "public.json", record["public"])
    if monotonic() >= campaign_deadline:
        result = {"id": record["id"], "status": "NOT_RUN_BUDGET", "metadata": record["metadata"]}
        save(job / "result.json", result)
        return result
    command = ["taskset", "-c", str(cpu), sys.executable, str(ROOT / "worker.py"),
               "--public", str(job / "public.json"), "--out", str(job),
               "--backend", manifest["backend"], "--address-space", str(manifest["address_space_limit_bytes"]),
               "--implementation", record["metadata"].get("implementation", "ours")]
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    launch_time = monotonic()
    peak, progress, status = 0, {}, None
    with (job / "console.log").open("x") as console:
        process = subprocess.Popen(command, stdout=console, stderr=subprocess.STDOUT, env=env, start_new_session=True)
        save(job / "process.json", {"pid": process.pid, "process_group": process.pid, "cpu": cpu,
                                   "command": command, "created_utc": utc()})
        while True:
            if (job / "progress.json").exists():
                progress = read(job / "progress.json")
            peak = max(peak, group_rss(process.pid))
            now = monotonic()
            source_start = progress.get("source_start_monotonic")
            if now >= campaign_deadline:
                status = "CAMPAIGN_BUDGET_REACHED"
            elif source_start is not None and now - source_start >= manifest["timeout_seconds"]:
                status = "TIMEOUT"
            elif source_start is None and now - launch_time >= 120:
                status = "STARTUP_TIMEOUT"
            elif peak > manifest["peak_rss_limit_bytes"]:
                status = "MEMORY_LIMIT"
            code = process.poll()
            if code is not None:
                break
            if status is not None:
                stop_group(process)
                break
            sleep(.2)
    if (job / "progress.json").exists():
        progress = read(job / "progress.json")
    output_path = job / "solver_result.json"
    solver_result = read(output_path) if output_path.exists() else None
    status = final_status(status, process.returncode, progress, solver_result,
                          manifest["timeout_seconds"], campaign_deadline)
    result = {"id": record["id"], "metadata": record["metadata"], "status": status,
              "returncode": process.returncode, "last_phase": progress.get("phase", "startup"),
              "source_elapsed_seconds": monotonic() - progress["source_start_monotonic"] if "source_start_monotonic" in progress else None,
              "last_checkpoint_elapsed_seconds": progress.get("elapsed_seconds"),
              "launch_elapsed_seconds": monotonic() - launch_time,
              "peak_sampled_process_group_rss_bytes": peak,
              "memory_measurement": "sampled summed RSS, shared pages may be double-counted",
              "progress_counters_are_lower_bounds": status != "COMPLETE",
              "last_checkpoint": progress, "solver_result": solver_result,
              "created_utc": utc()}
    save(job / "result.json", result)
    print(json.dumps({"id": record["id"], **record["metadata"], "status": status,
                      "last_phase": result["last_phase"]}), flush=True)
    return result


def run_campaign(run_path, cpus):
    from solver import FLATTER, FLATTER_LIBRARY
    run = Path(run_path).resolve()
    manifest = read(run / "manifest.json")
    assert manifest["public_sha256"] == sha(run / "public.jsonl")
    assert all(sha(ROOT / name) == value for name, value in manifest["source_sha256"].items())
    assert all(sha(path) == value for path, value in manifest["dependency_sha256"].items())
    assert sha(FLATTER) == manifest["flatter_executable_sha256"]
    assert sha(FLATTER_LIBRARY) == manifest["flatter_library_sha256"]
    assert sha(ROOT.parent / "HNP_Near_Limit_Experiment_Prompt_KO.md") == manifest["prompt_sha256"]
    assert sha(ROOT.parent / "STOC2027_HNPs (19).pdf") == manifest["paper_sha256"]
    assert len(cpus) == len(set(cpus)) <= manifest["maximum_workers"]
    assert set(cpus) <= os.sched_getaffinity(0)
    records = [json.loads(line) for line in (run / "public.jsonl").read_text().splitlines()]
    assert len(records) == manifest["planned_trials"]
    started = monotonic()
    clock_path = run.parent / "campaign_clock.json"
    boot_id = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    if not clock_path.exists():
        save(clock_path, {"created_utc": utc(), "boot_id": boot_id,
                          "campaign_start_monotonic": started,
                          "campaign_deadline_monotonic": started + manifest["campaign_budget_seconds"],
                          "budget_seconds": manifest["campaign_budget_seconds"]})
    clock = read(clock_path)
    assert clock["boot_id"] == boot_id, "Host reboot: campaign clock needs explicit reconciliation"
    assert clock["budget_seconds"] == manifest["campaign_budget_seconds"]
    deadline = clock["campaign_deadline_monotonic"]
    save(run / "execution.json", {"created_utc": utc(), "cpus": cpus, "pid": os.getpid(),
         "campaign_start_monotonic": clock["campaign_start_monotonic"], "campaign_deadline_monotonic": deadline,
         "stage_start_monotonic": started, "shared_campaign_clock_sha256": sha(clock_path),
         "timeout_seconds": manifest["timeout_seconds"], "budget_seconds": manifest["campaign_budget_seconds"]})
    slots = SimpleQueue()
    for cpu in cpus:
        slots.put(cpu)
    def one(record):
        cpu = slots.get()
        try:
            return execute(record, run, cpu, manifest, deadline)
        finally:
            slots.put(cpu)
    with ThreadPoolExecutor(max_workers=len(cpus)) as executor:
        futures = [executor.submit(one, record) for record in records]
        results = [future.result() for future in as_completed(futures)]
    save(run / "finished.json", {"created_utc": utc(), "jobs": len(results),
                                "elapsed_seconds": monotonic() - started,
                                "statuses": {s: sum(r["status"] == s for r in results) for s in {r["status"] for r in results}}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--cpus", default="14,15")
    args = parser.parse_args()
    run_campaign(args.run, [int(v) for v in args.cpus.split(",")])
