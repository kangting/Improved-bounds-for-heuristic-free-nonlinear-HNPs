"""Versioned controller: no per-job timer; adopt workers without restarting them.

The original committed sources/configuration are not modified. An explicit
user-authorized amendment supersedes only the per-job/startup time limits.
The shared 24-hour campaign clock and memory/CPU limits remain in force.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from time import monotonic, sleep

from common import ROOT, read, save, sha, utc
from runner import group_rss


def process_state(pid):
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
        fields = stat[stat.rfind(")") + 2:].split()
        return {"state": fields[0], "group": int(fields[2]), "start_ticks": int(fields[19])}
    except (FileNotFoundError, ProcessLookupError):
        return None


def alive(pid, start_ticks):
    state = process_state(pid)
    return state is not None and state["start_ticks"] == start_ticks and state["state"] != "Z"


def terminate_owned_group(pid, start_ticks):
    if not alive(pid, start_ticks):
        return
    assert process_state(pid)["group"] == pid
    os.killpg(pid, signal.SIGTERM)
    until = monotonic() + 2
    while alive(pid, start_ticks) and monotonic() < until:
        sleep(.05)
    if alive(pid, start_ticks):
        os.killpg(pid, signal.SIGKILL)


def valid_complete(progress, output, deadline):
    return bool(progress.get("complete") and output and
                progress.get("monotonic", float("inf")) <= deadline)


def take_over(run):
    """Stop the exact old supervisor only; its new-session workers continue."""
    old = read(run / "execution.json")
    pid = old["pid"]
    state = process_state(pid)
    adopted = []
    for path in (run / "jobs").glob("*/process.json"):
        if (path.parent / "result.json").exists():
            continue
        info = read(path)
        current = process_state(info["pid"])
        if current and current["state"] != "Z":
            cmd = Path(f"/proc/{info['pid']}/cmdline").read_bytes().split(b"\0")
            assert str(ROOT / "worker.py").encode() in cmd
            assert str(path.parent / "public.json").encode() in cmd
            assert current["group"] == info["pid"]
            adopted.append({"id": path.parent.name, "pid": info["pid"], **current})
    if state and state["state"] != "Z":
        cmd = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
        assert any(x.endswith(b"/runner.py") for x in cmd)
        assert any(b"pilot-v1" in x for x in cmd)
        # Never send a process-group signal here: preserve the live computations.
        os.kill(pid, signal.SIGTERM)
        until = monotonic() + 5
        while alive(pid, state["start_ticks"]) and monotonic() < until:
            sleep(.05)
        assert not alive(pid, state["start_ticks"]), "Old supervisor still alive"
    return adopted


def verify(run, manifest):
    assert sha(run / "public.jsonl") == manifest["public_sha256"]
    assert all(sha(ROOT / key) == value for key, value in manifest["source_sha256"].items())
    assert all(sha(path) == value for path, value in manifest["dependency_sha256"].items())
    from solver import FLATTER, FLATTER_LIBRARY
    assert sha(FLATTER) == manifest["flatter_executable_sha256"]
    assert sha(FLATTER_LIBRARY) == manifest["flatter_library_sha256"]
    assert sha(ROOT.parent / "HNP_Near_Limit_Experiment_Prompt_KO.md") == manifest["prompt_sha256"]
    assert sha(ROOT.parent / "STOC2027_HNPs (19).pdf") == manifest["paper_sha256"]


def execute(record, run, manifest, deadline, cpu):
    job = run / "jobs" / record["id"]
    if (job / "result.json").exists():
        return read(job / "result.json")
    adopted = job.exists()
    console = process = None
    if adopted:
        assert read(job / "public.json") == record["public"]
        info = read(job / "process.json")
        pid = info["pid"]
        state = process_state(pid)
        if state and state["state"] != "Z":
            cmd = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
            assert str(ROOT / "worker.py").encode() in cmd
            assert str(job / "public.json").encode() in cmd
            assert state["group"] == pid
        start_ticks = state["start_ticks"] if state else None
        cpu = info["cpu"]
    else:
        job.mkdir(exist_ok=False)
        save(job / "public.json", record["public"])
        if monotonic() >= deadline:
            result = {"id": record["id"], "metadata": record["metadata"], "status": "NOT_RUN_BUDGET"}
            save(job / "result.json", result)
            return result
        command = ["taskset", "-c", str(cpu), sys.executable, str(ROOT / "worker.py"),
                   "--public", str(job / "public.json"), "--out", str(job),
                   "--backend", manifest["backend"], "--address-space", str(manifest["address_space_limit_bytes"]),
                   "--implementation", record["metadata"].get("implementation", "ours")]
        env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        console = (job / "console.log").open("x")
        process = subprocess.Popen(command, stdout=console, stderr=subprocess.STDOUT,
                                   env=env, start_new_session=True)
        pid = process.pid
        state = process_state(pid)
        start_ticks = state["start_ticks"]
        save(job / "process.json", {"pid": pid, "process_group": pid, "cpu": cpu,
             "command": command, "created_utc": utc(), "start_ticks": start_ticks,
             "controller": "runner_unbounded.py"})
    peak, interruption, progress = 0, None, {}
    while True:
        if (job / "progress.json").exists():
            progress = read(job / "progress.json")
        output = read(job / "solver_result.json") if (job / "solver_result.json").exists() else None
        peak = max(peak, group_rss(pid))
        # Final output timestamp, not interpreter teardown, defines completion.
        if valid_complete(progress, output, deadline):
            break
        if not alive(pid, start_ticks):
            break
        if monotonic() >= deadline:
            interruption = "CAMPAIGN_BUDGET_REACHED"
        elif peak > manifest["peak_rss_limit_bytes"]:
            interruption = "MEMORY_LIMIT"
        if interruption:
            terminate_owned_group(pid, start_ticks)
            break
        sleep(.2)
    if process is not None:
        if process.poll() is not None:
            process.wait()
        if console:
            console.close()
    if (job / "progress.json").exists():
        progress = read(job / "progress.json")
    output = read(job / "solver_result.json") if (job / "solver_result.json").exists() else None
    # A recorded memory-limit interruption is never silently changed to success.
    complete = not interruption and valid_complete(progress, output, deadline)
    status = "COMPLETE" if complete else interruption or "ERROR"
    result = {"id": record["id"], "metadata": record["metadata"], "status": status,
              "returncode": process.poll() if process is not None else None,
              "returncode_unknown_for_adopted_worker": adopted,
              "last_phase": progress.get("phase", "startup"),
              "solver_output_elapsed_seconds": progress.get("elapsed_seconds") if complete else None,
              "monitored_source_elapsed_seconds": monotonic() - progress["source_start_monotonic"] if "source_start_monotonic" in progress else None,
              "last_checkpoint_elapsed_seconds": progress.get("elapsed_seconds"),
              "peak_sampled_process_group_rss_bytes": peak,
              "peak_memory_is_lower_bound": adopted,
              "memory_measurement": "sampled summed RSS; adopted-worker peak excludes pre-handoff samples",
              "progress_counters_are_lower_bounds": not complete,
              "last_checkpoint": progress, "solver_result": output,
              "effective_per_job_timeout_seconds": None,
              "policy_amendment": "timeout-amendment-v1.json", "created_utc": utc()}
    save(job / "result.json", result)
    print(json.dumps({"id": record["id"], **record["metadata"], "status": status,
                      "last_phase": result["last_phase"]}), flush=True)
    return result


def main(run_path, cpus, adopt):
    run = Path(run_path).resolve()
    manifest = read(run / "manifest.json")
    verify(run, manifest)
    clock = read(run.parent / "campaign_clock.json")
    assert clock["boot_id"] == Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    assert len(cpus) == len(set(cpus)) <= manifest["maximum_workers"]
    assert set(cpus) <= os.sched_getaffinity(0)
    adopted = take_over(run) if adopt else []
    save(run / "timeout-amendment-v1.json", {"created_utc": utc(), "user_authorized": True,
         "request": "Remove per-job time limit where possible; check flatter path overhead",
         "original_manifest_sha256": sha(run / "manifest.json"),
         "controller_sha256": sha(__file__), "per_job_timeout_seconds": None,
         "startup_timeout_seconds": None, "campaign_clock_preserved": clock,
         "memory_and_cpu_limits_unchanged": True, "adopted_workers": adopted,
         "historical_peak_rss_unavailable_for_adopted_workers": True})
    save(run / "execution-unbounded-v1.json", {"created_utc": utc(), "pid": os.getpid(),
         "cpus": cpus, "per_job_timeout_seconds": None, "shared_campaign_clock": clock})
    records = [json.loads(line) for line in (run / "public.jsonl").read_text().splitlines()]
    assert len(records) == manifest["planned_trials"]
    # Partition by original CPU so an adopted worker and a new job never share a lease.
    partitions = {cpu: [] for cpu in cpus}
    for index, record in enumerate(records):
        path = run / "jobs" / record["id"] / "process.json"
        cpu = read(path)["cpu"] if path.exists() else cpus[index % len(cpus)]
        assert cpu in partitions
        partitions[cpu].append(record)
    def lane(cpu, rows):
        return [execute(row, run, manifest, clock["campaign_deadline_monotonic"], cpu) for row in rows]
    with ThreadPoolExecutor(max_workers=len(cpus)) as pool:
        pending = [pool.submit(lane, cpu, rows) for cpu, rows in partitions.items()]
        results = [r for future in as_completed(pending) for r in future.result()]
    save(run / "finished-unbounded-v1.json", {"created_utc": utc(), "jobs": len(results),
         "statuses": {s: sum(r["status"] == s for r in results) for s in {r["status"] for r in results}}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--cpus", default="14,15")
    parser.add_argument("--adopt", action="store_true")
    args = parser.parse_args()
    main(args.run, [int(x) for x in args.cpus.split(",")], args.adopt)
