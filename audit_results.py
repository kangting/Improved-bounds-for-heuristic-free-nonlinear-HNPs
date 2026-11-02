"""Portable, post-hoc aggregation of stored records; never run an attack.

Only the Python standard library is imported. Stored outputs are compared
with stored evaluation secrets by equality, not by a new recovery or a
mathematical check of the original constraints.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent
RELEASE_BASE_COMMIT = "4a54acf8b3692c59091e3f4270f388077dd616a1"
FINISHED = {"COMPLETE", "RETURNED", "NO_CANDIDATE"}
GROUP = ("family", "method", "observation_model", "unknown_bits", "N", "source_observations_used")
IMPLEMENTATION = {"mihnp": "lssw12", "echnp": "shani17", "hnphm": "hgns03"}
RAW_FIELDS = (
    "id", "family", "method", "observation_model", "N", "unknown_bits", "secret_index", "trial",
    "status", "parent_elapsed_seconds", "worker_elapsed_seconds", "source_rank", "ambient_dimension",
    "source_observations_used", "nodes_visited", "decoder_calls", "output_size", "contains_true_secret",
    "singleton_true", "true_relevant_nonzero_quotients", "true_relevant_quotients_all_zero",
    "postrun_true_readout_verified", "outputs_original_constraints_valid",
)
COUNTERS = ("nodes_visited", "decoder_calls")
TIMES = ("parent_elapsed_seconds", "worker_elapsed_seconds")


def sha_bytes(value):
    return hashlib.sha256(value).hexdigest()


def load_lines(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def keyed(rows, label):
    result = {row["id"]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate IDs in {label}")
    return result


def csv_value(value):
    if value is None:
        return "NA"
    if isinstance(value, (list, dict, bool)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def csv_bytes(rows):
    fields = sorted(set().union(*(row.keys() for row in rows))) if rows else []
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({name: csv_value(row.get(name)) for name in fields})
    return stream.getvalue().encode("utf-8")


def recount(outputs, secret):
    if outputs is None:
        return None, None, None
    contains = any(candidate == secret for candidate in outputs)
    return len(outputs), contains, contains and len(outputs) == 1


def stats(rows, field, prefix, out):
    values = [row[field] for row in rows if isinstance(row.get(field), (int, float))
              and not isinstance(row.get(field), bool)]
    out[prefix + "_n_recorded"] = len(values)
    out[prefix + "_median"] = median(values) if values else None
    out[prefix + "_max"] = max(values) if values else None


def make_rows(cases, originals, truths):
    case_map = keyed(cases, "cases")
    result_map = keyed(originals, "results")
    truth_map = keyed(truths, "truth")
    if set(case_map) != set(result_map) or set(case_map) != set(truth_map):
        raise ValueError("This stored-study audit requires identical case/result/truth ID sets")
    output = []
    for original in sorted(originals, key=lambda row: row["id"]):
        case, truth = case_map[original["id"]], truth_map[original["id"]]
        for field in ("family", "method", "observation_model", "unknown_bits", "N", "secret_index", "trial"):
            if original[field] != case[field]:
                raise ValueError(f"Case/result {field} mismatch: {original['id']}")
        public = case["public"]
        if public["family"] != original["family"] or public["N"] != original["N"]:
            raise ValueError(f"Public/metadata mismatch: {original['id']}")
        if public["X"] != 1 << original["unknown_bits"]:
            raise ValueError(f"Width mismatch: {original['id']}")
        size, contains, singleton = recount(original.get("secrets"), truth["secret"])
        for name, value in [("output_size", size), ("contains_true_secret", contains), ("singleton_true", singleton)]:
            if original.get(name) is not None and original[name] != value:
                raise ValueError(f"Stored {name} disagrees with output/truth equality: {original['id']}")
        row = {field: original.get(field) for field in RAW_FIELDS}
        used = original.get("source_observations_used")
        ours = original["method"] == "ours"
        bits = public["p"].bit_length()
        row.update({
            "audit.field_bits": bits,
            "audit.revealed_fraction": (bits - original["unknown_bits"]) / bits,
            "audit.scalar_report_count": len(public["reports"]),
            "audit.input_unit": "multiplier_product_pairs" if original["family"] == "hnphm" else "observations",
            "audit.implementation_from_release": "recover_pruned" if ours else IMPLEMENTATION[original["family"]],
            "audit.source_indices0_from_release_code": list(range(used)) if used is not None else None,
            "audit.source_indices_evidence": "release_prefix_code_not_per_trial_indices_log",
            "audit.seed_recorded": truth.get("observation_seed", original.get("seed")),
            "audit.output_size_recounted": size,
            "audit.contains_true_recounted": contains,
            "audit.singleton_true_recounted": singleton,
            "audit.algorithm_finished": original["status"] in FINISHED,
            "audit.completion_contract": "complete_subset_list_then_full_input_filter" if ours else "single_answer_algorithm",
            "audit.root_node_completed": (original["status"] == "COMPLETE" and original["nodes_visited"] == 1)
                if original.get("nodes_visited") is not None else None,
            "audit.first_grid_cell_completion_recorded": None,
            "audit.distinct_tested_targets_recorded": None,
            "audit.distinct_targets_if_release_code_matches": original.get("decoder_calls") if ours else None,
            "audit.timeout_seconds_recorded": None,
            "audit.runner_default_timeout_seconds_from_release": 600,
            "audit.max_nodes_from_release": 4096 if ours else None,
        })
        for phase in ("source", "reduction", "search", "verification"):
            row[f"audit.{phase}_seconds_recorded"] = None
        output.append(row)
    return output


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[name] for name in GROUP)].append(row)
    output = []
    for key, group in sorted(groups.items()):
        done = [row for row in group if row["audit.algorithm_finished"]]
        singleton = [row for row in done if row["audit.singleton_true_recounted"] is True]
        interrupted = [row for row in group if not row["audit.algorithm_finished"]]
        counts = Counter(row["status"] for row in group)
        row = dict(zip(GROUP, key))
        row.update({
            "runs": len(group), "algorithm_finished": len(done), "singleton_true": len(singleton),
            "other_finished": len(done) - len(singleton), "interrupted": len(interrupted),
            "status_counts": dict(counts), "revealed_fraction": group[0]["audit.revealed_fraction"],
            "source_rank_values": sorted({r["source_rank"] for r in group if r["source_rank"] is not None}),
            "ambient_dimension_values": sorted({r["ambient_dimension"] for r in group if r["ambient_dimension"] is not None}),
            "first_grid_cell_completed_count": None,
            "distinct_tested_targets_recorded_median": None,
            "distinct_tested_targets_recorded_max": None,
        })
        for status in ("COMPLETE", "RETURNED", "NO_CANDIDATE", "CENSORED_WORK", "TIMEOUT", "ERROR"):
            row[status + "_count"] = counts[status]
        row["other_status_count"] = sum(count for name, count in counts.items()
            if name not in {"COMPLETE", "RETURNED", "NO_CANDIDATE", "CENSORED_WORK", "TIMEOUT", "ERROR"})
        roots = [r["audit.root_node_completed"] for r in done if r["audit.root_node_completed"] is not None]
        row["root_node_completed_n_recorded"] = len(roots)
        row["root_node_completed_count"] = sum(roots) if roots else None
        for label, subset in (("finished", done), ("singleton", singleton), ("interrupted", interrupted)):
            for field in TIMES + COUNTERS:
                stats(subset, field, label + "_" + field, row)
        output.append(row)
    return output


def cohorts(rows, truths):
    truth_map = keyed(truths, "truth")
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[name] for name in GROUP)].append(row)
    output = []
    for key, group in sorted(groups.items()):
        ids = defaultdict(list)
        for row in group:
            ids[row["secret_index"]].append(row)
        counts = {str(index): len(values) for index, values in sorted(ids.items())}
        trial_indices = {str(index): sorted(row["trial"] for row in values) for index, values in sorted(ids.items())}
        fixed = all(len({json.dumps(truth_map[row["id"]]["secret"], sort_keys=True) for row in values}) == 1
                    for values in ids.values())
        distinct = len({json.dumps(truth_map[values[0]["id"]]["secret"], sort_keys=True) for values in ids.values()})
        out = dict(zip(GROUP, key))
        out.update(secret_index_counts=counts, trial_indices=trial_indices,
            secret_fixed_within_index=fixed, distinct_secret_count=distinct,
            is_three_fixed_secrets_ten_trials=(len(ids) == 3 and distinct == 3 and fixed
                and all(indices == list(range(10)) for indices in trial_indices.values())))
        output.append(out)
    return output


def artifacts(root):
    inputs = ["config/cases.jsonl", "config/truth.jsonl", "results/trials.jsonl"]
    cases, truths, originals = [load_lines(root / name) for name in inputs]
    rows = make_rows(cases, originals, truths)
    summary, cohort_rows = summarize(rows), cohorts(rows, truths)
    coverage = [{"original_field": field, "present": sum(field in row for row in originals),
                 "nonnull": sum(row.get(field) is not None for row in originals),
                 "missing": sum(field not in row for row in originals)} for field in RAW_FIELDS]
    source_paths = ["run.py", "worker.py", "summarize.py", "plot.py", "requirements.txt"]
    source_paths += [path.relative_to(root).as_posix() for path in sorted((root / "src").glob("*.py"))]
    provenance = {
        "scope": "post-hoc stored-primary-study aggregation; no attack execution",
        "publication_base_commit": RELEASE_BASE_COMMIT,
        "publication_base_commit_is_not_a_claim_of_historical_execution_provenance": True,
        "input_sha256": {name: sha_bytes((root / name).read_bytes()) for name in inputs},
        "release_source_sha256_not_measured_code_certificates": {
            name: sha_bytes((root / name).read_bytes()) for name in source_paths},
        "aggregation_script_sha256": sha_bytes((root / "audit_results.py").read_bytes()),
        "runs": len(rows), "configurations": len(summary),
        "status_counts": dict(sorted(Counter(row["status"] for row in rows).items())),
        "singleton_true_by_stored_output_truth_equality": sum(row["audit.singleton_true_recounted"] is True for row in rows),
        "all_cohorts_three_fixed_secrets_ten_trials": all(row["is_three_fixed_secrets_ten_trials"] for row in cohort_rows),
        "original_constraint_or_quotient_math_reverified": False,
        "original_inputs_results_or_attacks_modified": False,
    }
    return {
        "runs.csv": csv_bytes(rows), "summary.csv": csv_bytes(summary),
        "cohorts.csv": csv_bytes(cohort_rows), "field_coverage.csv": csv_bytes(coverage),
        "provenance.json": (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--out", type=Path, help="New output directory; never overwritten")
    parser.add_argument("--check", action="store_true", help="Compare with results/audit without writing files")
    args = parser.parse_args()
    if args.check and args.out is not None:
        parser.error("Choose --check or --out")
    root = args.root.resolve()
    output = artifacts(root)
    if args.check:
        folder = root / "results/audit"
        failed = [name for name, data in output.items()
                  if not (folder / name).exists() or (folder / name).read_bytes() != data]
        if failed:
            parser.exit(1, "Stored audit differs: " + ", ".join(failed) + "\n")
        print("Stored audit matches all five regenerated files; no attack was run.")
    else:
        folder = args.out if args.out is not None else root / "work/audit"
        folder.mkdir(parents=True, exist_ok=False)
        for name, data in output.items():
            with (folder / name).open("xb") as handle:
                handle.write(data)
        print(f"Wrote five post-hoc audit files to {folder}; no attack was run.")


if __name__ == "__main__":
    main()
