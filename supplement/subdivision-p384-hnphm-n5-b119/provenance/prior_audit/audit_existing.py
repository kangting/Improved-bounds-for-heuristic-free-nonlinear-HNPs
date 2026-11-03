"""Audit existing logs/certificates only; no LLL, decoder, traversal, or new inputs."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from math import prod
from pathlib import Path
from statistics import median

WORKSPACE = Path(__file__).resolve().parents[2]
P384 = 2**384 - 2**128 - 2**96 + 2**32 - 1
PRIMARY = Path("Improved-bounds-for-heuristic-free-nonlinear-HNPs")
AUX = Path("lattice/hnp_near_limit_20261022/runs/cells-auxiliary-v1")
LEGACY = Path("lattice/hnp_threshold_experiments")


def read(path):
    return json.loads(path.read_text())


def records(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def determinant(rows):
    """Fraction-free determinant, used only for already stored transformations."""
    a = [list(row) for row in rows]
    n, sign, previous = len(a), 1, 1
    assert all(len(row) == n for row in a)
    for k in range(n - 1):
        if a[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if a[i][k]), None)
            if swap is None:
                return 0
            a[k], a[swap] = a[swap], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                value, remainder = divmod(numerator, previous)
                assert remainder == 0
                a[i][j] = value
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def hn_sources(public, target_last):
    """Reconstruct published HN relations, not a reduced basis or an attack."""
    p, n, reports = public["p"], public["N"], public["reports"]
    k = 2 * n + 1
    rows = []
    for i in range(n):
        row = [0] * (1 + n + k)  # constant, N cross monomials, 2N+1 linear terms
        row[0] = reports[0] * reports[i + 1] - reports[n + i + 1]
        row[1 + i] = 1
        row[1 + n] = reports[i + 1]
        row[1 + n + i + 1] = reports[0]
        row[1 + n + n + i + 1] = -1
        rows.append(row)
    order = list(range(1, n + 1)) + [0] if target_last else list(range(n + 1))
    for i in order:
        row = [0] * (1 + n + k)
        row[1 + n + i] = p
        rows.append(row)
    return rows


def polynomial_interval(row, n, box):
    low = high = row[0]
    for coefficient, (a, b) in zip(row[1 + n:], box):
        low += min(coefficient * a, coefficient * b)
        high += max(coefficient * a, coefficient * b)
    for i, coefficient in enumerate(row[1:1 + n]):
        corners = [coefficient * a * b for a in box[0] for b in box[i + 1]]
        low += min(corners)
        high += max(corners)
    return low, high


def weighted_target(alpha, intervals, first_coordinate):
    low = sum(min(alpha[j] * a, alpha[j] * b) for j, (a, b) in intervals.items())
    high = sum(max(alpha[j] * a, alpha[j] * b) for j, (a, b) in intervals.items())
    return max(low, first_coordinate[0]), min(high, first_coordinate[1])


def certificate_root(public, proof, target_last, clip_original):
    n, p, width = public["N"], public["p"], public["X"]
    source = hn_sources(public, target_last)
    u, alpha = proof.get("transform", proof.get("U")), proof["alpha"]
    m = 2 * n + 1
    assert len(u) == len(alpha) == m
    assert abs(determinant(u)) == 1
    target_index = m - 1 if target_last else n
    assert [sum(alpha[i] * u[i][j] for i in range(m)) for j in range(m)] == [int(j == target_index) for j in range(m)]
    transformed = [[sum(u[i][r] * source[r][c] for r in range(m))
                    for c in range(len(source[0]))] for i in range(m)]
    scaled = [[value * (width**2 if c < n else width) for c, value in enumerate(row[1:])]
              for row in transformed]
    if "R" in proof:
        assert scaled == proof["R"]
    relevant = [j for j, coefficient in enumerate(alpha) if coefficient]
    w = max(sum(map(abs, scaled[j])) for j in relevant)
    # Historical source19 declares signed boxes but intersects them with public
    # original leading-bit intervals. Legacy source9 uses its leading box.
    box = [(max(0, 1 - report), min(width - 1, p - 1 - report)) for report in public["reports"]] if clip_original else [(0, width - 1)] * m
    intervals = {}
    for j in relevant:
        low, high = polynomial_interval(transformed[j], n, box)
        interval = (-((-low) // p), high // p)
        assert interval[0] <= interval[1]
        intervals[j] = interval
    first = (max(box[0][0], 0), min(box[0][1], p - 1 - public["reports"][0]))
    low, high = weighted_target(alpha, intervals, first)
    return {"root_target_lower": low, "root_target_upper": high,
        "root_target_cardinality": max(0, high - low + 1),
        "ambiguous_relevant_quotient_rows": [j for j, (a, b) in intervals.items() if a != b],
        "source_rank": m, "ambient_dimension": 3 * n + 1,
        "relevant_rows": relevant, "W_rebuilt": w,
        "unimodular_and_target_identity_verified": True,
        "stored_R_identity_verified": "R" in proof}


def main_groups(workspace):
    original = records(workspace / PRIMARY / "results/trials.jsonl")
    groups = defaultdict(list)
    for row in original:
        if row["method"] == "ours":
            groups[(row["family"], row["observation_model"], row["unknown_bits"])].append(row)
    result = []
    for (family, model, bits), group in sorted(groups.items()):
        completed = [row for row in group if row["status"] == "COMPLETE"]
        capped = [row for row in group if row["status"] == "CENSORED_WORK"]
        times = [row["parent_elapsed_seconds"] for row in completed]
        result.append({"family": family, "query_model": model, "b": bits,
            "revealed_fraction": (384 - bits) / 384, "runs": len(group),
            "complete": len(completed), "singleton_true_recorded": sum(row.get("singleton_true") is True for row in group),
            "root_node_completions": sum(row["nodes_visited"] == 1 for row in completed),
            "nontrivial_tree_completions": sum(row["nodes_visited"] > 1 for row in completed),
            "work_censored": len(capped), "censored_nodes_values": sorted({row["nodes_visited"] for row in capped}),
            "completed_parent_seconds_median": median(times) if times else None,
            "completed_parent_seconds_max": max(times) if times else None,
            "other_statuses": dict(Counter(row["status"] for row in group if row["status"] not in {"COMPLETE", "CENSORED_WORK"}))})
    return result


def p384_auxiliary(workspace):
    run = workspace / AUX
    manifest = read(run / "manifest.json")
    for name, checksum in manifest["source_sha256"].items():
        assert sha(run.parent.parent / name) == checksum
    for name, checksum in manifest["dependency_sha256"].items():
        assert sha(Path(name)) == checksum
    truths = {row["id"]: row for row in records(run / "private/truth.jsonl")}
    evaluated = {row["id"]: row for row in records(workspace / "lattice/hnp_near_limit_20261022/evaluations/cells-auxiliary-checkpoint-v1/trials.jsonl")}
    output = []
    for record in records(run / "public.jsonl"):
        public = record["public"]
        if public["family"] != "hnphm" or public["N"] != 5:
            continue
        assert public["p"] == P384 and public["X"] == 2**119
        job = run / "jobs" / record["id"]
        result = read(job / "result.json")
        solver = result["solver_result"]
        with gzip.open(job / "certificate.json.gz", "rt") as handle:
            proof = json.load(handle)
        audit = certificate_root(public, proof, target_last=True, clip_original=True)
        prep = solver["preparation"]
        assert audit["W_rebuilt"] == prep["W"]
        assert audit["relevant_rows"] == prep["relevant_rows"]
        subdivisions = -((-12 * prep["W"]) // public["p"])
        assert subdivisions == prep["subdivisions"]
        assert prod(min(subdivisions, 2 * public["X"] + 1) for _ in public["reports"]) == prep["planned_cells"]
        assert result["status"] == "COMPLETE" and solver["outputs"] == [truths[record["id"]]["secret"]]
        assert evaluated[record["id"]]["singleton_true"] is True
        assert evaluated[record["id"]]["all_outputs_original_constraint_valid"] is True
        output.append({"id": record["id"], "family": "hnphm", "N": 5, "b": 119,
            "revealed_fraction": 265/384, "backend_recorded": "flatter-ntl", **audit,
            "subdivisions": subdivisions, "logical_cells": prep["planned_cells"],
            "status_recorded": result["status"], "visited_nodes_recorded": solver["visited_nodes"],
            "empty_target_subtrees_recorded": solver["empty_target_interval_subtrees"],
            "decoder_target_calls_recorded": solver["tested_target_proposals"],
            "pending_frontier_recorded": solver["pending_frontier_nodes"],
            "target_fallback_recorded": solver["exhaustive_target_scan_used"],
            "singleton_matches_stored_truth": True,
            "original_constraints_valid_by_existing_evaluator": True,
            "result_file": str(AUX / "jobs" / record["id"] / "result.json"),
            "certificate_sha256": sha(job / "certificate.json.gz")})
    assert len(output) == 3
    return output


def legacy_cohort(workspace):
    base = workspace / LEGACY
    public = {row["id"]: row["public"] for row in records(base / "runs/leakage_inputs_v2/public.jsonl")}
    selected = {row["id"] for row in read(base / "runs/decoded/leakage_v1/verification.json")["nontrivial_cells_beyond_singleton"]}
    replay = {row["id"]: row for row in read(base / "runs/replays/leakage_v3.json")["rows"]}
    result = []
    for native in records(base / "runs/decoded/leakage_v1/results.jsonl"):
        p = public[native["id"]]
        if native["route"] != "new" or p["family"] != "hnphm" or p["N"] != 1 or p["X"] != 2**15:
            continue
        assert p["p"] == 2**61 - 1
        certificate = base / "runs/decoded/leakage_v1/certificates" / (native["id"] + ".json")
        audit = certificate_root(p, read(certificate), target_last=False, clip_original=False)
        assert audit["W_rebuilt"] == native["preparation"]["W"]
        result.append({"id": native["id"], "N": 1, "b": 15, "p_bits": 61,
            **audit, "native_status": native["status"], "native_mode": native["preparation"]["search_mode"],
            "native_cells_or_targets_visited": native["visited"], "native_decoder_calls": native["decoder_calls"],
            "native_output_size": len(native["secrets"]),
            "selected_cell_completion_beyond_root": native["id"] in selected,
            "replay_status": replay[native["id"]]["status"], "replay_shortcut_used": replay[native["id"]]["shortcut_used"],
            "certificate_sha256": sha(certificate),
            "certificate_file": str(LEGACY / "runs/decoded/leakage_v1/certificates" / certificate.name)})
    assert len(result) == 20
    observed = [row for row in result if row["selected_cell_completion_beyond_root"]]
    assert len(observed) == 14
    assert all(row["root_target_cardinality"] > 1 and row["native_status"] == "COMPLETE"
               and row["native_mode"] == "cells" and row["native_decoder_calls"] == 1
               and row["replay_shortcut_used"] is False for row in observed)
    return result


def collect(workspace):
    main = main_groups(workspace)
    aux = p384_auxiliary(workspace)
    legacy = legacy_cohort(workspace)
    branching = [row for row in aux if row["visited_nodes_recorded"] > 1]
    assert len(branching) == 1 and branching[0]["root_target_cardinality"] > 1
    assert sum(row["complete"] for row in main) == 240
    assert sum(row["root_node_completions"] for row in main) == 240
    assert sum(row["work_censored"] for row in main) == 150
    return {"main_mechanism.csv": main,
        "ec_models.csv": [row for row in main if row["family"] == "echnp" and row["b"] in {124, 125, 126}],
        "p384_auxiliary.csv": aux, "legacy_n1.csv": legacy}


def write_csv(path, rows):
    keys = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else json.dumps(row[key], sort_keys=True)
                if isinstance(row[key], (dict, list, bool)) else row[key] for key in keys})


def main(workspace, out):
    tables = collect(workspace)
    out.mkdir(parents=True, exist_ok=False)
    for name, rows in tables.items():
        write_csv(out / name, rows)
    paths = [PRIMARY / "results/trials.jsonl", AUX / "public.jsonl", AUX / "manifest.json",
        LEGACY / "runs/decoded/leakage_v1/results.jsonl", LEGACY / "runs/decoded/leakage_v1/verification.json",
        LEGACY / "runs/replays/leakage_v3.json"]
    receipt = {"audit_script_sha256": sha(Path(__file__)),
        "original_input_record_sha256": {str(path): sha(workspace / path) for path in paths},
        "current_main_complete_root_only": 240, "current_main_work_censored": 150,
        "separate_p384_auxiliary_inputs": 3, "separate_p384_nontrivial_tree_completion": 1,
        "legacy_cohort_inputs": 20, "legacy_cell_completions_with_ambiguous_root_target": 14,
        "no_new_inputs_lll_decoder_or_traversal_executed": True,
        "no_new_timing_or_causal_speedup_measurement": True,
        "public_github_repository_not_modified": True}
    with (out / "verification.json").open("x") as handle:
        json.dump(receipt, handle, indent=2)
    print(json.dumps({"out": str(out), "tables": {name: len(rows) for name, rows in tables.items()},
                      "new_attacks_executed": False}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=WORKSPACE)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.workspace.resolve(), args.out.resolve())
