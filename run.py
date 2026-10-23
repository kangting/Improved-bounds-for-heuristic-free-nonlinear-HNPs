"""Run the stored MIHNP, EC-HNP, and HNPHM comparison inputs."""

import argparse
import concurrent.futures
import json
import os
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from hnp_native import key


def load_jsonl(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def postrun_quotients(record, coefficients, truth, modulus):
    errors = [truth["errors"][i] for i in coefficients["error_indices"]]
    quotients = []
    for polynomial in coefficients["relevant"]:
        value = polynomial["constant"]
        value += sum(c * x for c, x in zip(polynomial["linear"], errors))
        value += sum(c * errors[i] * errors[j]
                     for c, (i, j) in zip(polynomial["quadratic"], coefficients["edges"]))
        if value % modulus:
            raise ArithmeticError("True error does not satisfy a reduced congruence")
        quotients.append(value // modulus)
    readout = sum(coefficients["alpha"][f["index"]] * q
                  for f, q in zip(coefficients["relevant"], quotients))
    if readout != errors[0]:
        raise ArithmeticError("Post-run target identity failed")
    record["true_relevant_nonzero_quotients"] = sum(q != 0 for q in quotients)


def run_one(case, truth, timeout):
    record = {name: value for name, value in case.items() if name != "public"}
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="hnp_") as temporary:
        folder = Path(temporary)
        request = folder / "request.json"
        request.write_text(json.dumps({"method": case["method"], "public": case["public"]}))
        process = subprocess.Popen(
            [sys.executable, "-B", str(ROOT / "worker.py"), str(request), str(folder)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=os.name == "posix")
        timed_out = False
        try:
            _, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
            _, stderr = process.communicate()
        result_path = folder / "result.json"
        progress_path = folder / "progress.json"
        if progress_path.exists():
            record.update(json.loads(progress_path.read_text()))
        result = json.loads(result_path.read_text()) if result_path.exists() else {}
        coefficients = result.pop("_postrun", None)
        if timed_out:
            record.update(status="TIMEOUT", secrets=None)
        elif process.returncode:
            record.update({**result, "status": "ERROR", "secrets": None})
            record.setdefault("error", stderr.decode(errors="replace").strip()[-2000:])
        elif result:
            record.update(result)
        else:
            record.update(status="ERROR", secrets=None, error="Worker returned no result")

        record["parent_elapsed_seconds"] = time.perf_counter() - started
        record["contains_true_secret"] = None
        record["singleton_true"] = None
        record["output_size"] = None
        if record.get("secrets") is not None:
            secrets = record["secrets"]
            record["contains_true_secret"] = key(truth["secret"]) in {key(y) for y in secrets}
            record["singleton_true"] = len(secrets) == 1 and record["contains_true_secret"]
            record["output_size"] = len(secrets)
            if record["status"] == "COMPLETE" and coefficients is not None:
                postrun_quotients(record, coefficients, truth, case["public"]["p"])
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=("all", "ours", "prior"), default="all")
    parser.add_argument("--family", choices=("all", "mihnp", "echnp", "hnphm"), default="all")
    parser.add_argument("--model", choices=("all", "paired", "independent"), default="all")
    parser.add_argument("--unknown-bits", type=int)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--out", type=Path, default=ROOT / "results" / "new.jsonl")
    args = parser.parse_args()
    if args.workers < 1 or args.timeout <= 0 or (args.limit is not None and args.limit < 0):
        parser.error("workers and timeout must be positive, and limit must be nonnegative")
    if args.out.exists():
        parser.error(f"Output already exists; choose a new --out path: {args.out}")
    cases = [case for case in load_jsonl(ROOT / "config" / "cases.jsonl")
             if (args.method == "all" or case["method"] == args.method)
             and (args.family == "all" or case["family"] == args.family)
             and (args.model == "all" or case["observation_model"] == args.model)
             and (args.unknown_bits is None or case["unknown_bits"] == args.unknown_bits)]
    if args.limit is not None:
        cases = cases[:args.limit]
    truths = {row["id"]: row for row in load_jsonl(ROOT / "config" / "truth.jsonl")}
    missing = [case["id"] for case in cases if case["id"] not in truths]
    if missing:
        parser.error(f"Missing truth for {len(missing)} cases")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", buffering=1) as output:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_one, case, truths[case["id"]], args.timeout): case
                       for case in cases}
            for future in concurrent.futures.as_completed(futures):
                record = future.result()
                output.write(json.dumps(record, separators=(",", ":")) + "\n")
                print(json.dumps({name: record.get(name) for name in
                                  ("id", "method", "family", "unknown_bits", "status",
                                   "parent_elapsed_seconds", "singleton_true")}), flush=True)


if __name__ == "__main__":
    main()
