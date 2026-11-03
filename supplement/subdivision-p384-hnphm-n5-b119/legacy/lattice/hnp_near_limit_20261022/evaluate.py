"""Private, post-solver evaluation; never imported by a recovery worker."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path

from common import read, save, sha, digest, utc
from domain import PARAMETERS, validate
from solver import core, sources19


def valid_output(public, secret):
    """Original canonical clipped leading intervals, including hidden multipliers."""
    p, width, reports = public["p"], public["X"], public["reports"]
    family = public["family"]
    def inside(z, report):
        return report <= z <= min(report + width - 1, p - 1)
    try:
        if family == "mihnp":
            if not isinstance(secret, int) or not 0 <= secret < p:
                return False
            for shift, report in zip(public["shifts"], reports):
                denominator = (secret + shift) % p
                if not denominator or not inside(pow(denominator, -1, p), report):
                    return False
            return True
        if family == "hnphm":
            if not isinstance(secret, int) or not 0 < secret < p or not inside(secret, reports[0]):
                return False
            n = public["N"]
            for i in range(n):
                lo, hi = max(1, reports[i + 1]), min(p - 1, reports[i + 1] + width - 1)
                low = reports[n + i + 1]
                high = min(p - 1, low + width - 1)
                if core.modular_interval_count(secret, p, lo, hi, low, high) == 0:
                    return False
            return True
        ec = core.curve(public)
        candidate = core.point(ec, secret)
        if not (public["M"] * candidate).is_zero():
            return False
        for encoded, report in zip(public["queries"], reports):
            shifted = candidate + core.point(ec, encoded)
            if report is None:
                if not shifted.is_zero():
                    return False
            elif shifted.is_zero() or not inside(int(shifted[0]), report):
                return False
        return True
    except (ValueError, TypeError, ArithmeticError):
        return False


def regenerate(public, opening):
    rebuilt, private = core.generate(public["family"], public["p"], public["N"],
        public["X"].bit_length() - 1, opening["observation_seed"], opening["secret"],
        {key: public[key] for key in PARAMETERS} if public["family"] == "echnp" else None,
        public.get("observation_model", "independent"))
    assert rebuilt == public
    for key in ["errors", "multipliers", "query_scalars", "rejected_tuples"]:
        if key in opening:
            assert private[key] == opening[key]
    assert valid_output(public, opening["secret"])
    return True


def quotient_diagnostic(public, opening, result, certificate):
    """Forbidden until the complete-search result has actually been committed."""
    assert result["status"] == "COMPLETE"
    output = result["solver_result"]
    if output.get("preparation", {}).get("search_mode") == "infinity":
        return {"performed": False, "reason": "infinity shortcut, no reduced source"}
    with gzip.open(certificate, "rt") as handle:
        proof = json.load(handle)
    source = sources19(public)
    transform, alpha = proof["transform"], proof["alpha"]
    size = len(source.polynomials)
    assert len(transform) == len(alpha) == size
    assert all(len(row) == size for row in transform)
    assert [sum(alpha[i] * transform[i][j] for i in range(size)) for j in range(size)] == [0] * (size - 1) + [1]
    errors = opening["errors"]
    assert all(lo <= value <= hi for value, (lo, hi) in zip(errors, source.box))
    values = [polynomial.evaluate(errors, source.edges) for polynomial in source.polynomials]
    assert all(value % public["p"] == 0 for value in values)
    relevant = output["preparation"]["relevant_rows"]
    assert relevant == [i for i, a in enumerate(alpha) if a]
    quotients = {}
    for i in relevant:
        numerator = sum(c * value for c, value in zip(transform[i], values))
        assert numerator % public["p"] == 0
        quotients[str(i)] = numerator // public["p"]
    assert sum(alpha[i] * quotients[str(i)] for i in relevant) == errors[0]
    return {"performed": True, "performed_after_complete_search": True,
            "certificate_sha256": sha(certificate), "relevant_count": len(relevant),
            "nonzero_relevant_quotient_count": sum(value != 0 for value in quotients.values()),
            "has_nonzero_relevant_quotient": any(quotients.values()),
            "relevant_quotients": quotients, "target_readout_matches_true_error": True}


def evaluate_one(record, opening, job, verify_generation=True):
    public, metadata = record["public"], record["metadata"]
    assert record["id"] == opening["id"]
    assert digest(public) == record["public_sha256"]
    generated = regenerate(public, opening) if verify_generation else None
    row = {"id": record["id"], **metadata, "p": public["p"],
           "public_sha256": record["public_sha256"],
           "seed": opening["observation_seed"], "generation_verified": generated,
           "infinity_observations": public["reports"].count(None),
           "rejected_query_tuples": opening.get("rejected_tuples", 0),
           "contains_true_secret": None, "singleton_true": None,
           "complete_search": False, "output_size": None,
           "has_partial_valid_candidate": False,
           "has_partial_true_candidate": False,
           "q_diagnostic_performed": False,
           "has_nonzero_relevant_quotient": None}
    if public["family"] == "echnp":
        row.update(curve_parameters={key: public[key] for key in PARAMETERS})
    result_path = job / "result.json"
    if not result_path.exists():
        row["status"] = "RUNNING_OR_QUEUED" if job.exists() else "NOT_STARTED"
        return row, None
    result = read(result_path)
    row.update(status=result["status"], phase=result.get("last_phase"),
               result_sha256=sha(result_path),
               peak_sampled_process_group_rss_bytes=result.get("peak_sampled_process_group_rss_bytes"),
               peak_memory_is_lower_bound=result.get("peak_memory_is_lower_bound", False))
    completed = result["status"] == "COMPLETE"
    solver_output = result.get("solver_result") or {}
    checkpoint = result.get("last_checkpoint") or {}
    prep = solver_output.get("preparation") or checkpoint.get("preparation") or {}
    row.update({key: prep[key] for key in ["source_rank", "ambient_dimension", "planned_iterations",
         "planned_cells", "planned_target_proposals", "search_mode", "subdivisions", "W_over_p", "phase_seconds",
         "logical_leaf_count", "tree_worst_case_nodes", "planned_iterations_kind",
         "reference_literal_search_mode", "reference_literal_planned_iterations",
         "maximum_nodes", "original_interval_pruning", "preparation_reuse"] if key in prep})
    counts = solver_output if completed else checkpoint
    row.update({key: counts.get(key) for key in ["visited_outer_iterations", "visited_cells", "tested_target_proposals",
          "visited_nodes", "logical_cells_covered", "pending_frontier_nodes", "singleton_quotient_subtrees",
          "singleton_target_interval_subtrees", "empty_target_interval_subtrees", "termination_reason",
          "target_switch_budget", "effective_node_budget", "target_fallback_visited",
          "target_fallback_skipped_already_tested", "rejected_subtrees", "exhaustive_target_scan_used", "search_completion_certificate"]})
    row["progress_counters_are_lower_bounds"] = not completed
    row["total_output_seconds"] = checkpoint.get("elapsed_seconds") if completed else None
    row["monitored_source_seconds"] = result.get("monitored_source_elapsed_seconds", result.get("source_elapsed_seconds"))
    row["effective_per_job_timeout_seconds"] = result.get("effective_per_job_timeout_seconds")
    # Only actual on-time committed complete output counts as a complete-search metric.
    if completed:
        outputs = solver_output["outputs"]
        assert solver_output["public_sha256"] == record["public_sha256"]
        validities = [valid_output(public, candidate) for candidate in outputs]
        contains = core.key(opening["secret"]) in {core.key(candidate) for candidate in outputs}
        row.update(complete_search=True, output_size=len(outputs),
                   all_outputs_original_constraint_valid=all(validities),
                   contains_true_secret=contains, singleton_true=contains and len(outputs) == 1,
                   outcome="SINGLETON_TRUE" if contains and len(outputs) == 1 else
                           "CONTAINS_TRUE_LIST" if contains else "NO_ANSWER" if not outputs else "WRONG_CANDIDATE")
        diagnostic = quotient_diagnostic(public, opening, result, job / "certificate.json.gz")
        row["q_diagnostic_performed"] = diagnostic["performed"]
        row["has_nonzero_relevant_quotient"] = diagnostic.get("has_nonzero_relevant_quotient")
        row["nonzero_relevant_quotient_count"] = diagnostic.get("nonzero_relevant_quotient_count")
        return row, diagnostic
    # Conservatively include only checkpointed candidates at/before campaign cutoff.
    clock = read(job.parent.parent.parent / "campaign_clock.json")
    events = job / "events.jsonl"
    partial = {}
    if events.exists():
        for line in events.read_text().splitlines():
            event = json.loads(line)
            if event["monotonic"] <= clock["campaign_deadline_monotonic"]:
                for candidate in event.get("partial_outputs", []):
                    partial[core.key(candidate)] = candidate
    validities = [valid_output(public, candidate) for candidate in partial.values()]
    row.update(partial_output_size=len(partial),
               has_partial_valid_candidate=any(validities),
               has_partial_true_candidate=core.key(opening["secret"]) in partial,
               all_partial_outputs_original_constraint_valid=all(validities),
               outcome=result["status"])
    return row, None


def main(run_path, out_path, verify_generation=True):
    run, out = Path(run_path).resolve(), Path(out_path).resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / "private").mkdir()
    validate(prove=False)
    manifest = read(run / "manifest.json")
    assert sha(run / "public.jsonl") == manifest["public_sha256"]
    assert sha(run / "private/truth.jsonl") == manifest["truth_sha256"]
    records = [json.loads(line) for line in (run / "public.jsonl").read_text().splitlines()]
    openings = {row["id"]: row for row in map(json.loads, (run / "private/truth.jsonl").read_text().splitlines())}
    assert len(records) == len(openings) == manifest["planned_trials"]
    rows = []
    for record in records:
        row, diagnostic = evaluate_one(record, openings[record["id"]], run / "jobs" / record["id"], verify_generation)
        rows.append(row)
        if diagnostic is not None:
            save(out / "private" / (record["id"] + "-quotients.json"), diagnostic)
    with (out / "trials.jsonl").open("x") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with (out / "trials.csv").open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else value
                         for key, value in row.items()} for row in rows)
    summary = {"created_utc": utc(), "run": str(run), "evaluator_sha256": sha(__file__),
               "raw_sha256": sha(out / "trials.jsonl"), "planned_trials": len(rows),
               "statuses": dict(Counter(row["status"] for row in rows)),
               "completed_searches": sum(row["complete_search"] for row in rows),
               "complete_contains_true": sum(row["contains_true_secret"] is True for row in rows),
               "singleton_true": sum(row["singleton_true"] is True for row in rows),
               "private_diagnostics_after_complete_only": True,
               "all_inputs_regenerated": verify_generation}
    save(out / "summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--skip-regeneration", action="store_true")
    args = parser.parse_args()
    main(args.run, args.out, not args.skip_regeneration)
