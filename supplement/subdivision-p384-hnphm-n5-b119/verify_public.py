"""Exact checks and bounded verification replay of three stored public cases.

No truth file, legacy solver, lattice reducer, decoder, clock, or subprocess
is imported or used. The replay checks the historical splitting rule; it is
not a historical per-node trace or a new timed attack experiment.
"""
import csv
import gzip
import hashlib
import json
from math import prod
from pathlib import Path

IDS = ("7129544b15d539fb71e5e983", "210cec052c05c4c5e7ac8845", "4b2bf95e230d9cfa2a6c1929")
EXPECTED_NODES = dict(zip(IDS, (1, 5, 1)))
P = 2**384 - 2**128 - 2**96 + 2**32 - 1
X = 2**119


class VerificationError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def load(path):
    return json.loads(Path(path).read_text())


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def determinant(rows):
    a = [list(row) for row in rows]
    n, sign, previous = len(a), 1, 1
    require(n > 0 and all(len(row) == n for row in a), "Non-square determinant")
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
                value, remainder = divmod(a[i][j] * pivot - a[i][k] * a[k][j], previous)
                require(remainder == 0, "Non-exact Bareiss division")
                a[i][j] = value
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def multiply(a, b):
    return [[sum(row[k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for row in a]


def build_model(public, certificate):
    """Only public input and stored integer certificate enter interval maths."""
    require(public["family"] == "hnphm" and public["N"] == 5, "Wrong auxiliary setting")
    require(public["p"] == P and public["X"] == X and public["report_model"] == "leading_bits", "Wrong field/width/model")
    reports = public["reports"]
    require(len(reports) == 11 and all(type(v) is int and 0 <= v < P and v % X == 0 for v in reports), "Invalid reports")
    n, variables, rank, ambient = 5, 11, 11, 16
    # Coefficient layout: constant, five x0*x(i+1), eleven linear terms.
    source = []
    for i in range(n):
        row = [0] * (1 + ambient)
        row[0] = reports[0] * reports[i + 1] - reports[n + i + 1]
        row[1 + i] = 1
        row[1 + n] = reports[i + 1]
        row[1 + n + i + 1] = reports[0]
        row[1 + n + n + i + 1] = -1
        source.append(row)
    for i in list(range(1, n + 1)) + [0]:
        row = [0] * (1 + ambient)
        row[1 + n + i] = P
        source.append(row)
    b0 = [[value * (X**2 if j < n else X) for j, value in enumerate(row[1:])] for row in source]
    require(determinant([row[:rank] for row in b0]) != 0, "Source rank is not11")
    u, alpha = certificate["transform"], certificate["alpha"]
    require(len(u) == len(alpha) == rank and all(len(row) == rank for row in u), "Wrong certificate shape")
    require(all(type(v) is int for row in u for v in row) and all(type(v) is int for v in alpha), "Noninteger certificate")
    det_u = determinant(u)
    require(abs(det_u) == 1, "U is not unimodular")
    unit = [0] * (rank - 1) + [1]
    require(multiply([alpha], u)[0] == unit, "alpha U != e_target")
    h = multiply(u, source)  # Constants retained; never reduce these coefficients modp.
    r = multiply(u, b0)
    scaled_h = [[value * (X**2 if j < n else X) for j, value in enumerate(row[1:])] for row in h]
    require(r == scaled_h, "R != U B0 in full ambient space")
    require(multiply([alpha], h)[0] == source[-1], "Full target polynomial identity failed")
    relevant = [j for j, a in enumerate(alpha) if a]
    w = max(sum(map(abs, r[j])) for j in relevant)
    subdivisions = max(1, -((-12 * w) // P))
    declared = [(-X, X)] * variables
    original = [(max(0, 1 - v), min(X - 1, P - 1 - v)) for v in reports]
    return {"public": public, "source": source, "B0": b0, "U": u, "alpha": alpha,
        "R": r, "H": h, "det_U": det_u, "relevant": relevant, "W": w,
        "subdivisions": subdivisions, "declared_box": declared, "original_box": original,
        "parts": [min(subdivisions, 2 * X + 1)] * variables}


def polynomial_interval(polynomial, box):
    low = high = polynomial[0]
    for c, (a, b) in zip(polynomial[6:], box):
        low += min(c * a, c * b)
        high += max(c * a, c * b)
    for i, c in enumerate(polynomial[1:6]):
        endpoints = [c * a * b for a in box[0] for b in box[i + 1]]
        low += min(endpoints)
        high += max(endpoints)
    return low, high


def partition_union(interval, count, first, last):
    lo, hi = interval
    length = hi - lo + 1
    require(0 <= first <= last < count <= length, "Invalid balanced grid range")
    return lo + first * length // count, lo + (last + 1) * length // count - 1


def choose_split(model, node, box, ambiguous):
    available = [i for i, (a, b) in enumerate(node) if a < b]
    require(bool(available), "Unresolved original grid leaf")
    scores = {i: 0 for i in available}
    for j in ambiguous:
        polynomial = model["H"][j]
        for i in available:
            scores[i] += abs(polynomial[6 + i]) * (box[i][1] - box[i][0])
        for i, c in enumerate(polynomial[1:6]):
            a, b = 0, i + 1
            if a in scores:
                scores[a] += abs(c) * max(map(abs, box[b])) * (box[a][1] - box[a][0])
            if b in scores:
                scores[b] += abs(c) * max(map(abs, box[a])) * (box[b][1] - box[b][0])
    coordinate = max(available, key=lambda i: (scores[i], node[i][1] - node[i][0], -i))
    return coordinate, scores


def replay(model, maximum_recorded_nodes):
    """Public-only deterministic verification; never call reconstruction.

    At most the recorded1/5nodes are evaluated. Divergence is an error,
    not permission to launch a larger search or a target scan.
    """
    require(maximum_recorded_nodes in (1, 5), "Unexpected replay budget")
    root = tuple((0, count - 1) for count in model["parts"])
    stack, trace, targets = [("root", root)], [], set()
    covered = leaves = rejected = empty_target = singleton_q = singletons = repeats = 0
    while stack:
        require(len(trace) < maximum_recorded_nodes, "Verification traversal exceeds stored node count")
        name, node = stack.pop()
        weight = prod(b - a + 1 for a, b in node)
        leaf = all(a == b for a, b in node)
        leaves += int(leaf)
        raw = [partition_union(interval, count, a, b) for interval, count, (a, b)
               in zip(model["declared_box"], model["parts"], node)]
        box = [(max(a, lo), min(b, hi)) for (a, b), (lo, hi) in zip(raw, model["original_box"])]
        row = {"path": name, "grid_index_ranges": [list(v) for v in node],
            "logical_cell_weight": weight, "original_grid_leaf": leaf,
            "declared_union_box": [list(v) for v in raw], "public_intersected_box": [list(v) for v in box],
            "relevant_intervals": {}}
        if any(a > b for a, b in box):
            row["decision"] = "EMPTY_PUBLIC_BOX"
            covered += weight
            rejected += 1
            trace.append(row)
            continue
        intervals, ambiguous = {}, []
        for j in model["relevant"]:
            low, high = polynomial_interval(model["H"][j], box)
            qlo, qhi = -((-low) // P), high // P
            row["relevant_intervals"][str(j)] = {"polynomial_bounds": [low, high], "quotient_interval": [qlo, qhi]}
            if qlo > qhi:
                row["decision"] = "EMPTY_QUOTIENT_INTERVAL"
                break
            intervals[j] = (qlo, qhi)
            if qlo != qhi:
                ambiguous.append(j)
        if len(intervals) != len(model["relevant"]):
            covered += weight
            rejected += 1
        else:
            require(not (leaf and ambiguous), "Original grid leaf has ambiguous quotient")
            low = sum(min(model["alpha"][j] * a, model["alpha"][j] * b) for j, (a, b) in intervals.items())
            high = sum(max(model["alpha"][j] * a, model["alpha"][j] * b) for j, (a, b) in intervals.items())
            final = (max(low, box[0][0]), min(high, box[0][1]))
            row.update(weighted_target_before_intersection=[low, high], final_target_interval=list(final),
                       final_target_cardinality=max(0, final[1] - final[0] + 1), ambiguous_quotient_rows=ambiguous)
            if final[0] > final[1]:
                row["decision"] = "EMPTY_TARGET_INTERVAL"
                covered += weight
                rejected += 1
                empty_target += 1
            elif final[0] == final[1]:
                row["decision"] = "SINGLETON_TARGET"
                covered += weight
                singleton_q += int(not ambiguous)
                singletons += 1
                repeats += int(final[0] in targets)
                targets.add(final[0])
            else:
                coordinate, scores = choose_split(model, node, box, ambiguous)
                a, b = node[coordinate]
                middle = (a + b) // 2
                left, right = list(node), list(node)
                left[coordinate], right[coordinate] = (a, middle), (middle + 1, b)
                row.update(decision="SPLIT", split_coordinate0=coordinate, split_grid_midpoint=middle,
                           sensitivity_scores={str(i): score for i, score in scores.items()})
                stack.extend([(name + "/R", tuple(right)), (name + "/L", tuple(left))])
        trace.append(row)
    logical_cells = prod(model["parts"])
    require(covered == logical_cells, "Verification frontier did not cover original grid")
    return {"artifact_kind": "new_verification_replay_not_historical_trace", "decoder_executed": False,
        "new_timing_measurements": None, "nodes": trace, "visited_nodes": len(trace),
        "visited_grid_leaves": leaves, "logical_cells_covered": covered,
        "rejected_subtrees": rejected, "empty_target_interval_subtrees": empty_target,
        "singleton_quotient_subtrees": singleton_q, "singleton_target_interval_subtrees": singletons,
        "distinct_singleton_targets": sorted(targets), "duplicate_target_requests_avoided": repeats,
        "pending_frontier_nodes": 0, "exhaustive_target_scan_executed": False}


def floor_sum(n, m, a, b):
    answer = 0
    require(n >= 0 and m > 0, "Invalid floor-sum range")
    while True:
        qa, a = divmod(a, m)
        qb, b = divmod(b, m)
        answer += qa * n * (n - 1) // 2 + qb * n
        ymax = a * n + b
        if ymax < m:
            return answer
        n, b, m, a = ymax // m, ymax % m, a, m


def valid_stored_output(public, secret):
    """Check a stored candidate's constraints, not generate new candidates."""
    reports = public["reports"]
    if type(secret) is not int or not 0 < secret < P or not reports[0] <= secret <= min(P - 1, reports[0] + X - 1):
        return False
    for i in range(5):
        lo, hi = max(1, reports[i + 1]), min(P - 1, reports[i + 1] + X - 1)
        lower, upper = reports[6 + i], min(P - 1, reports[6 + i] + X - 1)
        count = floor_sum(hi - lo + 1, P, secret, secret * lo - lower) - floor_sum(hi - lo + 1, P, secret, secret * lo - upper - 1)
        if count <= 0:
            return False
    return True


def verify_trial(root, identifier, compare_published_replay=True):
    require(identifier in IDS, "Trial outside this three-input supplement")
    job = root / "records/trials" / identifier
    public, solver, result = [load(job / name) for name in ("public.json", "solver_result.json", "result.json")]
    with gzip.open(job / "certificate.json.gz", "rt") as handle:
        certificate = json.load(handle)
    model = build_model(public, certificate)
    replayed = replay(model, EXPECTED_NODES[identifier])
    prep = solver["preparation"]
    require(solver["public_sha256"] == digest(public), "Saved public digest mismatch")
    require(result["status"] == solver["status"] == "COMPLETE" and result["solver_result"] == solver, "Saved result inconsistency")
    require(result["id"] == identifier, "Wrong saved result ID")
    selected = {row["id"]: row for row in [json.loads(line) for line in (root / "records/public.selected.jsonl").read_text().splitlines()]}
    require(set(selected) == set(IDS), "Public selection differs from requested three inputs")
    require(selected[identifier]["public"] == public and selected[identifier]["public_sha256"] == digest(public), "Original public-row binding differs")
    require(selected[identifier]["metadata"] == result["metadata"], "Historical selection metadata differs")
    require(prep["source_rank"] == 11 and prep["ambient_dimension"] == 16 and prep["W"] == model["W"], "Saved lattice summary mismatch")
    require(prep["relevant_rows"] == model["relevant"] and prep["subdivisions"] == model["subdivisions"], "Relevant rows/subdivisions mismatch")
    require(prep["planned_cells"] == prod(model["parts"]), "Saved logical cell count mismatch")
    require(prep["preparation_reuse"] is None and solver["original_interval_pruning"] is True, "Unexpected historical preparation policy")
    require(solver["maximum_nodes"] == 10000 and result["effective_per_job_timeout_seconds"] is None, "Unexpected resource policy")
    require(solver["first_answer_early_stop"] is False and solver["exhaustive_target_scan_used"] is False, "Unexpected first-hit/fallback route")
    require(solver["complete_search"] is True and solver["termination_reason"] == "FRONTIER_EXHAUSTED"
            and solver["search_completion_certificate"] == "certified_cell_frontier_exhausted", "Saved completeness declaration differs")
    require(solver["target_switch_budget"] == prep["planned_target_proposals"] == 2 * X + 1
            and solver["effective_node_budget"] == 10000, "Historical target/node budget mismatch")
    for old, new in [("visited_nodes", "visited_nodes"), ("visited_cells", "visited_grid_leaves"),
                     ("logical_cells_covered", "logical_cells_covered"), ("rejected_subtrees", "rejected_subtrees"),
                     ("empty_target_interval_subtrees", "empty_target_interval_subtrees"),
                     ("singleton_quotient_subtrees", "singleton_quotient_subtrees"),
                     ("singleton_target_interval_subtrees", "singleton_target_interval_subtrees"),
                     ("repeated_target_proposals_avoided", "duplicate_target_requests_avoided"),
                     ("pending_frontier_nodes", "pending_frontier_nodes")]:
        require(solver[old] == replayed[new], "Historical aggregate differs: " + old)
    targets = replayed["distinct_singleton_targets"]
    require(targets == solver["tested_target_values"] and len(targets) == solver["tested_target_proposals"] == 1, "Historical target-call mismatch")
    require(solver["visited_outer_iterations"] == replayed["visited_nodes"], "Historical outer counter differs")
    require(solver["outputs"] == [public["reports"][0] + targets[0]], "Stored output differs from target readout")
    require(all(valid_stored_output(public, value) for value in solver["outputs"]), "Stored output violates original intervals")
    progress, process = load(job / "progress.json"), load(job / "process.json")
    require(progress == result["last_checkpoint"] and progress["complete"] is True, "Historical final checkpoint differs")
    events = [json.loads(line) for line in (job / "events.jsonl").read_text().splitlines()]
    require(events[-1] == progress and any(row.get("phase") == "adaptive_cell_recovery" for row in events), "Historical event/checkpoint mismatch")
    require(all("nodes" not in row and "box" not in row for row in events), "Unexpected historical per-node trace schema")
    manifest = load(root / "records/campaign/manifest.json")
    require(process["solver_sha256"] == manifest["source_sha256"]["solver_cells.py"] and process["worker_sha256"] == manifest["source_sha256"]["worker_cells.py"], "Historical process source binding differs")
    artifact = {"trial_id": identifier, "public_sha256": sha(job / "public.json"),
        "original_certificate_sha256": sha(job / "certificate.json.gz"),
        "basis_evidence": {"certificate_stores_U_and_alpha_not_R": True,
            "B0_reconstructed_now": model["B0"], "R_reconstructed_now": model["R"],
            "full_H_reconstructed_now": model["H"], "det_U": model["det_U"],
            "alpha": model["alpha"], "relevant_rows": model["relevant"], "W": model["W"],
            "subdivisions": model["subdivisions"], "declared_source_box": [list(v) for v in model["declared_box"]]},
        "verification_replay": replayed, "stored_outputs_publicly_validated": True}
    if compare_published_replay:
        require((root / "verification_replay" / (identifier + ".json")).read_bytes() == json_bytes(artifact), "Published verification replay differs")
    return artifact


def verify_integrity(root, evaluation_only=False):
    manifest = load(root / "manifest.json")
    count = 0
    for entry in manifest["files"]:
        name = entry["path"]
        is_evaluation = name.startswith("evaluation/")
        if is_evaluation != evaluation_only:
            continue
        relative = Path(name)
        require(not relative.is_absolute() and ".." not in relative.parts, "Unsafe manifest path")
        path = root / relative
        require(path.is_file() and sha(path) == entry["published_sha256"], "Public file hash mismatch: " + name)
        count += 1
    return count


def verify_legacy_bindings(root):
    manifest = load(root / "records/campaign/manifest.json")
    count = 0
    for name, checksum in manifest["source_sha256"].items():
        require(Path(name).name == name, "Unsafe historical source filename")
        require(sha(root / "legacy/lattice/hnp_near_limit_20261022" / name) == checksum, "Historical source hash mismatch: " + name)
        count += 1
    for name, checksum in manifest["dependency_sha256"].items():
        suffix = name[name.index("lattice/"):]
        require(not Path(suffix).is_absolute() and ".." not in Path(suffix).parts, "Unsafe historical dependency path")
        require(sha(root / "legacy" / suffix) == checksum, "Historical dependency hash mismatch: " + suffix)
        count += 1
    require(count == 15, "Expected15historical source bindings")
    execution = load(root / "records/campaign/execution-cells-v1.json")
    require(execution["maximum_nodes"] == 10000 and execution["per_job_timeout_seconds"] is None, "Historical execution policy mismatch")
    return count


def check_public(root):
    verify_integrity(root)
    verify_legacy_bindings(root)
    outputs = {identifier: verify_trial(root, identifier) for identifier in IDS}
    with (root / "provenance/prior_audit/p384_auxiliary.csv").open(newline="") as handle:
        previous = {row["id"]: row for row in csv.DictReader(handle)}
    for identifier, artifact in outputs.items():
        require(previous[identifier]["certificate_sha256"] == artifact["original_certificate_sha256"], "Recovered certificate differs from earlier audit")
        root_node = artifact["verification_replay"]["nodes"][0]
        require(int(previous[identifier]["root_target_cardinality"]) == root_node["final_target_cardinality"], "Earlier root-cardinality audit differs")
    return outputs
