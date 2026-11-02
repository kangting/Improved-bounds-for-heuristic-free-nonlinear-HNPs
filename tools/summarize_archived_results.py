#!/usr/bin/env python3
"""Validate and summarize the archived 1,050-run release; run no experiments.

Only the Python standard library is used. Inputs are read from --repo-root;
only configurations.csv, endpoints.csv and provenance.json under --out are
written. Neither secrets nor absolute input paths are written to the reports.
"""

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics


SOURCE_FILES = {
    "trials": "results/trials.jsonl",
    "cases": "config/cases.jsonl",
    "truth": "config/truth.jsonl",
}
FINISHED = {"COMPLETE", "RETURNED", "NO_CANDIDATE"}
INTERRUPTED = {"CENSORED_WORK"}
METADATA = ("family", "N", "unknown_bits", "observation_model",
            "secret_index", "trial", "method")
MEASUREMENTS = ("parent_elapsed_seconds", "worker_elapsed_seconds",
                "nodes_visited", "decoder_calls")
IMPLEMENTATIONS = {
    ("mihnp", "prior"): "lssw12",
    ("echnp", "prior"): "shani17",
    ("hnphm", "prior"): "hgns03",
}
ENDPOINTS = (
    ("mihnp", "ours", "independent", 173),
    ("echnp", "ours", "paired", 124),
    ("echnp", "ours", "independent", 124),
    ("hnphm", "ours", "independent", 123),
)
P384 = 2 ** 384 - 2 ** 128 - 2 ** 96 + 2 ** 32 - 1


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def read_records(path):
    payload = path.read_bytes()
    records = {}
    for line_number, line in enumerate(payload.decode("utf-8").splitlines(), 1):
        if not line.strip():
            continue
        row = json.loads(line)
        require(isinstance(row, dict), f"{path.name}:{line_number}: expected object")
        identifier = row.get("id")
        require(isinstance(identifier, str), f"{path.name}:{line_number}: missing id")
        require(identifier not in records, f"{path.name}: duplicate id {identifier}")
        records[identifier] = row
    return records, hashlib.sha256(payload).hexdigest()


def group_key(row):
    return (row["family"], row["method"], row["observation_model"],
            row["unknown_bits"])


def constant(rows, field):
    values = {canonical(row.get(field)) for row in rows}
    require(len(values) == 1, f"Nonconstant {field} in {group_key(rows[0])}")
    return rows[0].get(field)


def validate(trials, cases, truths):
    require(set(trials) == set(cases) == set(truths), "Input ID sets differ")
    require(len(trials) == 1050, "This reporter expects exactly 1,050 archived runs")
    groups = collections.defaultdict(list)
    output_checks = 0
    for identifier, row in trials.items():
        case = cases[identifier]
        for field in METADATA:
            require(row.get(field) == case.get(field),
                    f"Metadata mismatch: {identifier}, {field}")
        require(row["method"] in {"ours", "prior"}, f"Unknown method: {identifier}")
        require(row["family"] in {"mihnp", "echnp", "hnphm"},
                f"Unknown family: {identifier}")
        require(row["status"] in FINISHED | INTERRUPTED,
                f"Unexpected archived status: {identifier}")
        expected_statuses = ({"COMPLETE", "CENSORED_WORK"}
                             if row["method"] == "ours"
                             else {"RETURNED", "NO_CANDIDATE"})
        require(row["status"] in expected_statuses,
                f"Method/status mismatch: {identifier}")
        require(case["public"]["p"] == P384, f"Unexpected field: {identifier}")
        require(case["public"]["N"] == row["N"], f"Public N mismatch: {identifier}")
        require("secret" in truths[identifier], f"Missing truth: {identifier}")
        output = row.get("secrets")
        flag_fields = ("output_size", "contains_true_secret", "singleton_true")
        if output is None:
            require(row["status"] in INTERRUPTED,
                    f"Finished run has no saved output: {identifier}")
            require(all(row.get(field) is None for field in flag_fields),
                    f"Unavailable output has non-null flags: {identifier}")
        else:
            require(row["status"] in FINISHED and isinstance(output, list),
                    f"Invalid saved output: {identifier}")
            truth = canonical(truths[identifier]["secret"])
            contains = any(canonical(candidate) == truth for candidate in output)
            singleton = len(output) == 1 and contains
            require(type(row.get("output_size")) is int
                    and row["output_size"] == len(output),
                    f"Output size mismatch: {identifier}")
            require(row.get("contains_true_secret") is contains,
                    f"True-secret inclusion mismatch: {identifier}")
            require(row.get("singleton_true") is singleton,
                    f"Singleton mismatch: {identifier}")
            output_checks += 3
        for field in MEASUREMENTS[:2]:
            value = row.get(field)
            require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
                    f"Invalid recorded time: {identifier}, {field}")
        if row["method"] == "ours":
            for field in MEASUREMENTS[2:]:
                require(type(row.get(field)) is int and row[field] >= 0,
                        f"Invalid recorded counter: {identifier}, {field}")
        groups[group_key(row)].append(row)

    require(len(groups) == 35, "This reporter expects exactly 35 configurations")
    for key, rows in groups.items():
        require(len(rows) == 30, f"Expected 30 runs: {key}")
        by_secret = collections.defaultdict(list)
        for row in rows:
            by_secret[row["secret_index"]].append(row)
        require(set(by_secret) == {0, 1, 2}, f"Expected three secret indices: {key}")
        for index, cohort in by_secret.items():
            require(len(cohort) == 10 and {r["trial"] for r in cohort} == set(range(10)),
                    f"Expected ten inputs per secret: {key}, {index}")
            require(len({canonical(truths[r["id"]]["secret"]) for r in cohort}) == 1,
                    f"Truth changes within secret index: {key}, {index}")

    status_counts = dict(sorted(collections.Counter(r["status"] for r in trials.values()).items()))
    require(status_counts == {"CENSORED_WORK": 150, "COMPLETE": 240,
                              "NO_CANDIDATE": 353, "RETURNED": 307},
            "Archived status totals differ from the expected release")
    completed = [r for r in trials.values() if r["status"] == "COMPLETE"]
    censored = [r for r in trials.values() if r["status"] == "CENSORED_WORK"]
    require(all(r["nodes_visited"] == 1 and r["decoder_calls"] == 1
                and r["singleton_true"] is True for r in completed),
            "Completed ours runs do not all have the reported root-node behavior")
    require(all(r["nodes_visited"] == 4096 and r["secrets"] is None for r in censored),
            "Censored runs do not all match the observed 4,096-node stop")
    return groups, {
        "matching_unique_ids": len(trials),
        "configurations": len(groups),
        "runs_per_configuration": 30,
        "secrets_per_configuration": 3,
        "inputs_per_secret": 10,
        "output_flag_equality_checks": output_checks,
        "singleton_true_recounted": sum(r["singleton_true"] is True for r in trials.values()),
        "status_counts": status_counts,
        "ours_completed_with_nodes_1_and_decoder_calls_1": len(completed),
        "censored_with_nodes_4096_and_null_output": len(censored),
    }


def add_statistics(summary, label, rows):
    for field in MEASUREMENTS:
        values = [row[field] for row in rows if row.get(field) is not None]
        summary[f"{label}.{field}.n_recorded"] = len(values)
        summary[f"{label}.{field}.median"] = statistics.median(values) if values else "NA"
        summary[f"{label}.{field}.max"] = max(values) if values else "NA"


def summarize(rows):
    first = rows[0]
    family, method, model, bits = group_key(first)
    statuses = collections.Counter(row["status"] for row in rows)
    finished = [r for r in rows if r["status"] in FINISHED]
    singleton = [r for r in finished if r["singleton_true"] is True]
    interrupted = [r for r in rows if r["status"] in INTERRUPTED]
    full_count = constant(rows, "N")
    summary = {
        "family": family,
        "method": method,
        "implementation": "recover_pruned" if method == "ours" else IMPLEMENTATIONS[(family, method)],
        "query_model": model,
        "field_bits": 384,
        "unknown_bits": bits,
        "revealed_fraction": (384 - bits) / 384,
        "full_input_count": full_count,
        "full_input_unit": "pairs" if family == "hnphm" else "observations",
        "lattice_input_count": constant(rows, "source_observations_used"),
        "full_scalar_report_count": 2 * full_count + 1 if family == "hnphm" else full_count,
        "source_rank": constant(rows, "source_rank"),
        "ambient_dimension": constant(rows, "ambient_dimension"),
        "records": len(rows),
        "finished_count": len(finished),
        "singleton_true_count": len(singleton),
        "other_finished_count": len(finished) - len(singleton),
        "interrupted_count": len(interrupted),
        "complete_count": statuses["COMPLETE"],
        "returned_count": statuses["RETURNED"],
        "no_candidate_count": statuses["NO_CANDIDATE"],
        "censored_work_count": statuses["CENSORED_WORK"],
        "output_unavailable_count": sum(r["secrets"] is None for r in rows),
        "root_node_complete_count": sum(r["status"] == "COMPLETE" and r.get("nodes_visited") == 1
                                        for r in rows),
    }
    for label, subset in (("finished", finished), ("singleton", singleton),
                          ("interrupted", interrupted)):
        add_statistics(summary, label, subset)
    return summary


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True,
                        help="repository containing the three archived JSONL inputs")
    parser.add_argument("--out", type=Path, required=True,
                        help="directory for the three derived reporting files")
    args = parser.parse_args()
    loaded = {}
    sources = {}
    for role, relative_path in SOURCE_FILES.items():
        records, digest = read_records(args.repo_root / relative_path)
        loaded[role] = records
        sources[relative_path] = {"sha256": digest, "records": len(records)}
    groups, checks = validate(loaded["trials"], loaded["cases"], loaded["truth"])
    summaries = {key: summarize(rows) for key, rows in sorted(groups.items())}
    endpoints = [summaries[key] for key in ENDPOINTS]
    require(all(r["singleton_true_count"] == 30 and r["finished_count"] == 30 for r in endpoints),
            "Expected four endpoints with 30/30 singleton-true recoveries")
    output_names = ("configurations.csv", "endpoints.csv", "provenance.json")
    input_paths = {(args.repo_root / relative_path).resolve() for relative_path in SOURCE_FILES.values()}
    require(all((args.out / name).resolve() not in input_paths for name in output_names),
            "A reporting output would overwrite an input")
    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / output_names[0], list(summaries.values()))
    write_csv(args.out / output_names[1], endpoints)
    provenance = {
        "schema_version": 1,
        "scope": "archived_1050_run_release_only",
        "source_files": sources,
        "validation": checks,
        "units": {"elapsed_time": "seconds", "revealed_fraction": "fraction of 384 bits"},
        "missing_value": "NA means absent or inapplicable; it is not zero",
        "statistics": "n_recorded, median and max; finished, singleton and interrupted are separate",
        "verification_scope": "Saved-output/truth equality only; recorded mathematical verification flags are not recomputed",
        "execution": "No attack, lattice reduction, decoder, input generation or experiment script is executed",
    }
    (args.out / output_names[2]).write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n",
                                          encoding="utf-8")
    print(json.dumps(checks, sort_keys=True))


if __name__ == "__main__":
    main()
