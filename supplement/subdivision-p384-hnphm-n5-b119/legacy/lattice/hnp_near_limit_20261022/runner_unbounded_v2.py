"""No job timers; preserve completed output and supervise entire worker groups."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from time import monotonic, sleep

from common import ROOT, read, save, sha, utc, digest
from runner_unbounded import process_state, alive, verify, valid_complete


def group_snapshot(group):
    members, rss = [], 0
    for path in Path("/proc").iterdir():
        if not path.name.isdigit():
            continue
        try:
            stat = (path / "stat").read_text()
            fields = stat[stat.rfind(")") + 2:].split()
            if int(fields[2]) == group and fields[0] != "Z":
                members.append(int(path.name))
                rss += int(fields[21]) * os.sysconf("SC_PAGE_SIZE")
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return members, rss


def stop_owned(group):
    """The caller has validated this job's recorded, dedicated session group."""
    if not group_snapshot(group)[0]:
        return
    try:
        os.killpg(group, signal.SIGTERM)
    except ProcessLookupError:
        return
    until = monotonic() + 2
    while monotonic() < until and group_snapshot(group)[0]:
        sleep(.05)
    if group_snapshot(group)[0]:
        try:
            os.killpg(group, signal.SIGKILL)
        except ProcessLookupError:
            pass


def status_of(interruption, progress, output, deadline):
    if interruption == "MEMORY_LIMIT":
        return interruption
    if valid_complete(progress, output, deadline):
        return "COMPLETE"
    return interruption or ("LATE_COMPLETE" if progress.get("complete") and output else "ERROR")


def takeover(run):
    execution = read(run / "execution-unbounded-v1.json")
    pid = execution["pid"]
    state = process_state(pid)
    receipt = {"old_controller_pid": pid, "numerical_workers_restarted": False}
    if state and state["state"] != "Z":
        argv = [x.decode() for x in Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0") if x]
        assert any(Path(x).resolve() == ROOT / "runner_unbounded.py" for x in argv[1:] if not x.startswith("--"))
        assert "--run" in argv and Path(argv[argv.index("--run") + 1]).resolve() == run
        try:
            os.kill(pid, signal.SIGTERM)  # This PID only; never the worker groups.
        except ProcessLookupError:
            pass
        until = monotonic() + 5
        while alive(pid, state["start_ticks"]) and monotonic() < until:
            sleep(.05)
        assert not alive(pid, state["start_ticks"])
    return receipt


def execute(record, run, manifest, deadline, cpu):
    job = run / "jobs" / record["id"]
    if (job / "result.json").exists():
        return read(job / "result.json")
    adopted = job.exists()
    process, console = None, None
    if adopted:
        assert digest(read(job / "public.json")) == record["public_sha256"]
        info = read(job / "process.json")
        pid = info["pid"]
        state = process_state(pid)
        if state and state["state"] != "Z":
            argv = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
            assert str(ROOT / "worker.py").encode() in argv
            assert str(job / "public.json").encode() in argv
            assert state["group"] == pid
            assert set(os.sched_getaffinity(pid)) == {cpu}
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
        save(job / "process.json", {"pid": pid, "process_group": pid, "cpu": cpu,
             "command": command, "start_ticks": state["start_ticks"], "created_utc": utc(),
             "controller": "runner_unbounded_v2.py"})
    peak, interruption, cleanup_reason, progress, output = 0, None, None, {}, None
    while True:
        if (job / "progress.json").exists():
            progress = read(job / "progress.json")
        if (job / "solver_result.json").exists():
            output = read(job / "solver_result.json")
            assert output["public_sha256"] == record["public_sha256"]
        members, rss = group_snapshot(pid)
        peak = max(peak, rss)
        if rss > manifest["peak_rss_limit_bytes"]:
            interruption = "MEMORY_LIMIT"
            stop_owned(pid)
        elif monotonic() >= deadline and members:
            if valid_complete(progress, output, deadline):
                cleanup_reason = "CAMPAIGN_DEADLINE_DURING_POST_OUTPUT_TEARDOWN"
            else:
                interruption = "CAMPAIGN_BUDGET_REACHED"
            stop_owned(pid)
        if not members or interruption or cleanup_reason:
            # Reap newly owned children; do not release the CPU while descendants live.
            while group_snapshot(pid)[0]:
                sleep(.05)
            break
        sleep(.2)
    if process is not None:
        process.wait()
        console.close()
    if (job / "progress.json").exists():
        progress = read(job / "progress.json")
    output = read(job / "solver_result.json") if (job / "solver_result.json").exists() else None
    status = status_of(interruption, progress, output, deadline)
    result = {"id": record["id"], "metadata": record["metadata"], "status": status,
              "returncode": process.returncode if process else None,
              "returncode_unknown_for_adopted_worker": adopted,
              "last_phase": progress.get("phase", "startup"),
              "solver_output_elapsed_seconds": progress.get("elapsed_seconds") if status == "COMPLETE" else None,
              "monitored_source_elapsed_seconds": monotonic() - progress["source_start_monotonic"] if "source_start_monotonic" in progress else None,
              "last_checkpoint_elapsed_seconds": progress.get("elapsed_seconds"),
              "peak_sampled_process_group_rss_bytes": peak, "peak_memory_is_lower_bound": adopted,
              "memory_measurement": "sampled summed RSS; adopted peak excludes pre-handoff samples",
              "post_output_cleanup_interruption": cleanup_reason,
              "progress_counters_are_lower_bounds": status != "COMPLETE",
              "last_checkpoint": progress, "solver_result": output,
              "effective_per_job_timeout_seconds": None, "effective_cpu": cpu,
              "policy_amendment": "timeout-amendment-v2.json", "created_utc": utc()}
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
    receipt = takeover(run) if adopt else None
    save(run / "timeout-amendment-v2.json", {"created_utc": utc(), "user_authorized": True,
         "per_job_timeout_seconds": None, "startup_timeout_seconds": None,
         "original_manifest_sha256": sha(run / "manifest.json"), "controller_sha256": sha(__file__),
         "controller_dependency_sha256": {"runner_unbounded.py": sha(ROOT / "runner_unbounded.py")},
         "campaign_clock_preserved": clock, "memory_limit_unchanged": True,
         "distinct_physical_core_cpus": cpus, "takeover": receipt,
         "historical_peak_rss_unavailable_for_adopted_workers": True})
    save(run / "execution-unbounded-v2.json", {"created_utc": utc(), "pid": os.getpid(),
         "cpus": cpus, "per_job_timeout_seconds": None, "shared_campaign_clock": clock})
    records = [json.loads(line) for line in (run / "public.jsonl").read_text().splitlines()]
    assert len(records) == manifest["planned_trials"]
    partitions = {cpu: [] for cpu in cpus}
    for index, record in enumerate(records):
        path = run / "jobs" / record["id"] / "process.json"
        cpu = read(path)["cpu"] if path.exists() else cpus[index % len(cpus)]
        if cpu == 15 and (run / "affinity-amendment-v1.json").exists():
            cpu = 16
        assert cpu in partitions
        partitions[cpu].append(record)
    def lane(cpu, rows):
        return [execute(row, run, manifest, clock["campaign_deadline_monotonic"], cpu) for row in rows]
    with ThreadPoolExecutor(max_workers=len(cpus)) as pool:
        pending = [pool.submit(lane, cpu, rows) for cpu, rows in partitions.items()]
        results = [r for future in as_completed(pending) for r in future.result()]
    save(run / "finished-unbounded-v2.json", {"created_utc": utc(), "jobs": len(results),
         "statuses": {s: sum(r["status"] == s for r in results) for s in {r["status"] for r in results}}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--cpus", default="14,16")
    parser.add_argument("--adopt", action="store_true")
    args = parser.parse_args()
    main(args.run, [int(x) for x in args.cpus.split(",")], args.adopt)
