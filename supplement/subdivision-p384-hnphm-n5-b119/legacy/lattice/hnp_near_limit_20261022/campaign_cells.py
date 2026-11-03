"""Paired adaptive-cell pilots/main/auxiliaries, shared leases and24h deadline."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import os
from pathlib import Path
from queue import SimpleQueue
import subprocess
import sys
from time import monotonic, sleep

from common import ROOT, read, save, atomic, sha, utc, digest
from runner_unbounded import verify
from runner_unbounded_v2 import group_snapshot, stop_owned, valid_complete


def write_skip(run, record, status, reason, evidence=None):
    job = run / "jobs" / record["id"]
    if (job / "result.json").exists():
        return read(job / "result.json")
    job.mkdir(exist_ok=False)
    save(job / "public.json", record["public"])
    result = {"id": record["id"], "metadata": record["metadata"], "status": status,
              "reason": reason, "preflight_evidence": evidence, "solver_result": None,
              "source_construction_started": False, "created_utc": utc(),
              "complete_search": False, "effective_per_job_timeout_seconds": None}
    save(job / "result.json", result)
    return result


def public_prep_reuse(record):
    reference = ROOT / "runs" / record["metadata"]["reference_run"] / "jobs" / record["metadata"]["reference_record_id"]
    progress = reference / "progress.json"
    if progress.exists() and (reference / "certificate.json.gz").exists():
        data = read(progress)
        if "preparation" in data:
            assert digest(read(reference / "public.json")) == record["public_sha256"]
            return reference
    return None


def execute(run, record, manifest, cpu, clock):
    job = run / "jobs" / record["id"]
    if (job / "result.json").exists():
        return read(job / "result.json")
    deadline = clock["campaign_deadline_monotonic"]
    lease_dir = ROOT / "runs/cpu-leases"
    lease_dir.mkdir(exist_ok=True)
    with (lease_dir / f"cpu{cpu}.lock").open("a+") as lease:
        while True:
            try:
                fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if monotonic() >= deadline:
                    return write_skip(run, record, "NOT_RUN_BUDGET", "Shared24h budget expired while waiting for CPU lease")
                sleep(.2)
        if monotonic() >= deadline:
            return write_skip(run, record, "NOT_RUN_BUDGET", "Shared24h budget expired before launch")
        job.mkdir(exist_ok=False)
        save(job / "public.json", record["public"])
        command = ["taskset", "-c", str(cpu), sys.executable, str(ROOT / "worker_cells.py"),
                   "--public", str(job / "public.json"), "--out", str(job),
                   "--backend", manifest["backend"], "--address-space", str(manifest["address_space_limit_bytes"]),
                   "--maximum-nodes", str(manifest["maximum_nodes"]),
                   "--maximum-targets", str(manifest["maximum_target_fallback"])]
        reuse = public_prep_reuse(record)
        if reuse:
            command.extend(["--reuse-preparation", str(reuse)])
        env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        peak, interruption, progress, output = 0, None, {}, None
        launched = monotonic()
        with (job / "console.log").open("x") as console:
            process = subprocess.Popen(command, stdout=console, stderr=subprocess.STDOUT,
                                       env=env, start_new_session=True)
            save(job / "process.json", {"pid": process.pid, "process_group": process.pid,
                 "command": command, "cpu": cpu, "created_utc": utc(),
                 "worker_sha256": sha(ROOT / "worker_cells.py"),
                 "solver_sha256": sha(ROOT / "solver_cells.py"),
                 "lease": str(lease_dir / f"cpu{cpu}.lock")})
            while True:
                if (job / "progress.json").exists():
                    progress = read(job / "progress.json")
                if (job / "solver_result.json").exists():
                    output = read(job / "solver_result.json")
                    assert output["public_sha256"] == record["public_sha256"]
                members, rss = group_snapshot(process.pid)
                peak = max(peak, rss)
                if rss > manifest["peak_rss_limit_bytes"]:
                    interruption = "MEMORY_LIMIT"
                elif monotonic() >= deadline and members:
                    interruption = "CAMPAIGN_BUDGET_REACHED"
                if interruption:
                    stop_owned(process.pid)
                if not members or interruption:
                    while group_snapshot(process.pid)[0]:
                        sleep(.05)
                    break
                sleep(.3)
            process.wait()
        if (job / "progress.json").exists():
            progress = read(job / "progress.json")
        output = read(job / "solver_result.json") if (job / "solver_result.json").exists() else None
        if interruption == "MEMORY_LIMIT":
            status = interruption
        elif valid_complete(progress, output, deadline):
            status = "COMPLETE"
        elif output and progress.get("terminal") and progress.get("monotonic", float("inf")) <= deadline:
            status = output["status"]
        else:
            status = interruption or "ERROR"
        result = {"id": record["id"], "metadata": record["metadata"], "status": status,
                  "returncode": process.returncode, "last_phase": progress.get("phase"),
                  "last_checkpoint": progress, "solver_result": output,
                  "solver_output_elapsed_seconds": progress.get("elapsed_seconds") if status == "COMPLETE" else None,
                  "monitored_source_elapsed_seconds": monotonic() - progress["source_start_monotonic"] if "source_start_monotonic" in progress else None,
                  "new_worker_monitored_seconds": monotonic() - launched,
                  "peak_sampled_process_group_rss_bytes": peak,
                  "peak_memory_is_lower_bound": bool(reuse),
                  "progress_counters_are_lower_bounds": status != "COMPLETE",
                  "effective_per_job_timeout_seconds": None, "created_utc": utc()}
        save(job / "result.json", result)
        print(json.dumps({"run": run.name, "id": record["id"], "problem": record["metadata"]["problem"],
                          "unknown_bits": record["metadata"]["unknown_bits"], "status": status,
                          "visited_nodes": (output or {}).get("visited_nodes")}), flush=True)
        return result


def main(runs_root, cpus):
    root = Path(runs_root).resolve()
    assert root == ROOT / "runs"
    assert len(set(cpus)) == len(cpus) == 2
    assert set(cpus) <= os.sched_getaffinity(0)
    # Include the preserved MI preparation in the global three-core bound.
    topology = [(Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/physical_package_id").read_text().strip(),
                 Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/core_id").read_text().strip()) for cpu in [14] + cpus]
    assert len(set(topology)) == 3
    clock = read(root / "campaign_clock.json")
    assert clock["boot_id"] == Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    stages = {}
    for name in ["pilot", "controls", "main", "auxiliary"]:
        run = root / f"cells-{name}-v1"
        manifest = read(run / "manifest.json")
        verify(run, manifest)
        assert manifest["algorithm_arm"] == "adaptive_cell_tree"
        assert manifest["effective_per_job_timeout_seconds"] is None
        rows = [json.loads(x) for x in (run / "public.jsonl").read_text().splitlines()]
        assert len(rows) == manifest["planned_trials"]
        save(run / "execution-cells-v1.json", {"created_utc": utc(), "pid": os.getpid(),
             "cpus": cpus, "shared_campaign_clock": clock, "per_job_timeout_seconds": None,
             "maximum_nodes": manifest["maximum_nodes"], "maximum_target_fallback": manifest["maximum_target_fallback"],
             "global_slots_including_preserved_preparation": 3, "source_worker_receives_no_private_opening": True})
        stages[name] = (run, manifest, rows)
    slots = SimpleQueue()
    for cpu in cpus:
        slots.put(cpu)
    def batch(stage, rows):
        run, manifest, _ = stages[stage]
        def one(record):
            cpu = slots.get()
            try:
                return execute(run, record, manifest, cpu, clock)
            finally:
                slots.put(cpu)
        with ThreadPoolExecutor(max_workers=len(cpus)) as pool:
            return list(pool.map(one, rows))
    def family_rows(stage, family):
        return [r for r in stages[stage][2] if r["metadata"]["problem"] == family]
    for family in ["hnphm", "echnp", "mihnp"]:
        pilots = family_rows("pilot", family)
        if family == "mihnp":
            # Smaller source dimensions remain useful if the primary preparation is slow.
            batch("auxiliary", family_rows("auxiliary", family))
            old_job = ROOT / "runs/pilot-v1/jobs/fbdf8aa9a5989914b119e546"
            old_progress = read(old_job / "progress.json")
            elapsed = monotonic() - old_progress["source_start_monotonic"]
            if "preparation" not in old_progress and elapsed >= stages["pilot"][1]["maximum_new_preparation_cost_preflight_seconds"]:
                for stage in ["pilot", "controls", "main"]:
                    for record in family_rows(stage, family):
                        write_skip(stages[stage][0], record, "NOT_RUN_RESOURCE_PREFLIGHT",
                            "Same-dimensional primary preparation is still running beyond3600s; no additional heavy preparations launched",
                            {"reference_job": str(old_job), "observed_preparation_elapsed_lower_bound": elapsed})
                continue
        pilot_results = batch("pilot", pilots)
        controls = family_rows("controls", family)
        control_results = batch("controls", controls)
        lookup = {record["id"]: result for record, result in zip(pilots + controls, pilot_results + control_results)}
        main_rows = family_rows("main", family)
        for bits in sorted({r["metadata"]["unknown_bits"] for r in main_rows}, reverse=True):
            prerequisites = [r for r in pilots + controls if r["metadata"]["unknown_bits"] == bits]
            assert len(prerequisites) == 3
            evidence = [{"run": "pilot" if r in pilots else "controls", "id": r["id"],
                         "status": lookup[r["id"]]["status"]} for r in prerequisites]
            candidates = [r for r in main_rows if r["metadata"]["unknown_bits"] == bits]
            assert len(candidates) == 30
            if all(item["status"] == "COMPLETE" for item in evidence):
                batch("main", candidates)
            else:
                for record in candidates:
                    write_skip(stages["main"][0], record, "NOT_RUN_RESOURCE_PREFLIGHT",
                        "Corresponding paired pilots did not all complete within committed resource work caps; not a mathematical failure",
                        evidence)
        if family != "mihnp":
            batch("auxiliary", family_rows("auxiliary", family))
    for stage, (run, manifest, records) in stages.items():
        statuses = {}
        for record in records:
            result = read(run / "jobs" / record["id"] / "result.json")
            statuses[result["status"]] = statuses.get(result["status"], 0) + 1
        save(run / "finished-cells-v1.json", {"created_utc": utc(), "planned_trials": len(records), "statuses": statuses})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-root", default=str(ROOT / "runs"))
    parser.add_argument("--cpus", default="16,18")
    args = parser.parse_args()
    main(args.runs_root, [int(x) for x in args.cpus.split(",")])
