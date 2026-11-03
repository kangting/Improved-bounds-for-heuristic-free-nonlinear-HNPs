"""Complete public adaptive traversal of the original balanced cell family.

Unlike literal Algorithm1, certified nodes can cover entire subtrees of C.
No quotient Cartesian product, planted pruning or first-hit stop is used.
C is tried first; only a publicly small target fallback may finish a capped
tree. A nonempty frontier without that fallback is NOT completion.
"""
import gzip
import json
from math import prod
from pathlib import Path
import shutil
from time import monotonic

from sage.all import ZZ, matrix, vector
from common import read, digest, sha
from solver import core, sources19, prepare, solve as literal_solve


def public_error_box(public):
    assert public["report_model"] == "leading_bits"
    nonzero = public["family"] in {"mihnp", "hnphm"}
    return tuple((max(0, int(nonzero) - report),
                  min(public["X"] - 1, public["p"] - 1 - report))
                 for report in public["reports"])


def partition_union(interval, count, first, last):
    lo, hi = interval
    length = hi - lo + 1
    assert 0 <= first <= last < count <= length
    return lo + first * length // count, lo + (last + 1) * length // count - 1


def prepare_cells(public, backend, event, certificate_path, reuse=None):
    if reuse is None:
        prep = prepare(public, backend, event, certificate_path)
        reused = None
    else:
        reuse = Path(reuse)
        assert digest(read(reuse / "public.json")) == digest(public)
        progress = read(reuse / "progress.json")
        old_summary = progress["preparation"]
        event(phase="restore_public_preparation")
        start = monotonic()
        with gzip.open(reuse / "certificate.json.gz", "rt") as handle:
            proof = json.load(handle)
        source = sources19(public)
        transform, alpha = matrix(ZZ, proof["transform"]), proof["alpha"]
        m = len(source.polynomials)
        assert transform.nrows() == transform.ncols() == len(alpha) == m
        assert abs(transform.det()) == 1
        assert vector(ZZ, alpha) * transform == vector(ZZ, [0] * (m - 1) + [1])
        relevant = [i for i, a in enumerate(alpha) if a]
        assert relevant == old_summary["relevant_rows"]
        reduced = transform * source.scaled_matrix()
        assert old_summary["W"] == max(sum(abs(int(v)) for v in reduced.row(i)) for i in relevant)
        coefficients = transform * source.coefficient_matrix()
        constants = transform * vector(ZZ, [f.constant for f in source.polynomials])
        q = len(source.edges)
        polys = tuple(core.Polynomial(int(constants[i]), tuple(map(int, coefficients.row(i)[q:])),
                           tuple(map(int, coefficients.row(i)[:q]))) for i in range(m))
        summary = dict(old_summary)
        summary["phase_seconds"] = {k: v for k, v in old_summary["phase_seconds"].items() if k != "recovery"}
        summary["phase_seconds"]["restore_public_preparation"] = monotonic() - start
        shutil.copyfile(reuse / "certificate.json.gz", certificate_path)
        reused = {"directory": str(reuse.resolve()), "certificate_sha256": sha(reuse / "certificate.json.gz"),
                  "original_source_start_monotonic": progress["source_start_monotonic"],
                  "original_source_phase_seconds": old_summary["phase_seconds"],
                  "same_mathematical_public_input": True}
        prep = {"sources": source, "transform": transform, "alpha": alpha,
                "reduced": reduced, "polynomials": polys, "summary": summary}
    summary = prep["summary"]
    assert summary["relevant_rows"] == [i for i, a in enumerate(prep["alpha"]) if a]
    summary["reference_literal_search_mode"] = summary["search_mode"]
    summary["reference_literal_planned_iterations"] = summary["planned_iterations"]
    summary["search_mode"] = "adaptive_cell_tree"
    summary["logical_leaf_count"] = summary["planned_cells"]
    summary["tree_worst_case_nodes"] = 2 * summary["planned_cells"] - 1
    summary["planned_iterations"] = summary["tree_worst_case_nodes"]
    summary["planned_iterations_kind"] = "worst_case_tree_nodes_before_work_cap"
    summary["declared_cell_family_unchanged"] = True
    summary["preparation_reuse"] = reused
    return prep


def choose_split(public, prep, node, box, ambiguous):
    """Public coefficient sensitivity; ties are deterministic, no true witness."""
    available = [i for i, (a, b) in enumerate(node) if a < b]
    if not available:
        raise ArithmeticError("An original leaf cell failed quotient isolation")
    scores = {i: 0 for i in available}
    for j in ambiguous:
        polynomial = prep["polynomials"][j]
        for i in available:
            scores[i] += abs(polynomial.linear[i]) * (box[i][1] - box[i][0])
        for c, (a, b) in zip(polynomial.quadratic, prep["sources"].edges):
            if a in scores:
                scores[a] += abs(c) * max(map(abs, box[b])) * (box[a][1] - box[a][0])
            if b in scores:
                scores[b] += abs(c) * max(map(abs, box[a])) * (box[b][1] - box[b][0])
    return max(available, key=lambda i: (scores[i], node[i][1] - node[i][0], -i))


def target_interval(alpha, quotient_intervals, coordinate_interval):
    """Sign-aware interval image of the exact target readout, no tuples."""
    lower = upper = 0
    assert set(quotient_intervals) == {i for i, a in enumerate(alpha) if a}
    for i, (lo, hi) in quotient_intervals.items():
        assert lo <= hi
        a = alpha[i]
        lower += min(a * lo, a * hi)
        upper += max(a * lo, a * hi)
    return max(lower, coordinate_interval[0]), min(upper, coordinate_interval[1])


def recover_cells(public, prep, event, maximum_nodes=10000, tighten_original=True,
                  target_switch_budget=None):
    assert maximum_nodes >= 0
    source, summary = prep["sources"], prep["summary"]
    relevant = summary["relevant_rows"]
    assert relevant == [i for i, a in enumerate(prep["alpha"]) if a]
    parts = [min(summary["subdivisions"], hi - lo + 1) for lo, hi in source.box]
    assert prod(parts) == summary["planned_cells"]
    root = tuple((0, count - 1) for count in parts)
    stack, answers, decoded = [root], {}, set()
    original = public_error_box(public) if tighten_original else source.box
    visited = leaves = covered = rejected = certified = calls = repeats = 0
    target_empty = target_singleton = 0
    node_budget = min(maximum_nodes, target_switch_budget) if target_switch_budget is not None else maximum_nodes
    last = monotonic()

    def checkpoint(force=False):
        nonlocal last
        now = monotonic()
        if force or now - last >= 1:
            event(phase="adaptive_cell_recovery", visited_outer_iterations=visited,
                  visited_nodes=visited, visited_cells=leaves,
                  logical_cells_covered=covered, rejected_subtrees=rejected,
                  singleton_quotient_subtrees=certified, tested_target_proposals=calls,
                  empty_target_interval_subtrees=target_empty,
                  singleton_target_interval_subtrees=target_singleton,
                  repeated_target_proposals_avoided=repeats, pending_frontier_nodes=len(stack),
                  partial_outputs=list(answers.values()))
            last = now

    while stack and visited < node_budget:
        node = stack.pop()
        visited += 1
        is_leaf = all(lo == hi for lo, hi in node)
        leaves += is_leaf
        weight = prod(hi - lo + 1 for lo, hi in node)
        raw_box = [partition_union(interval, count, lo, hi)
                   for interval, count, (lo, hi) in zip(source.box, parts, node)]
        box = [(max(lo, bound[0]), min(hi, bound[1])) for (lo, hi), bound in zip(raw_box, original)]
        if any(lo > hi for lo, hi in box):
            covered += weight
            rejected += 1
            checkpoint()
            continue
        intervals, ambiguous, impossible = {}, [], False
        for j in relevant:
            lo, hi = prep["polynomials"][j].interval(box, source.edges)
            low, high = core.ceildiv(lo, public["p"]), hi // public["p"]
            if low > high:
                impossible = True
                break
            intervals[j] = (low, high)
            if low != high:
                ambiguous.append(j)
        if impossible:
            covered += weight
            rejected += 1
        else:
            if is_leaf and ambiguous:
                raise ArithmeticError("An original cell contains multiple possible relevant quotients")
            target_lo, target_hi = target_interval(prep["alpha"], intervals, box[0])
            if target_lo > target_hi:
                covered += weight
                rejected += 1
                target_empty += 1
            elif target_lo == target_hi:
                covered += weight
                certified += not ambiguous
                target_singleton += 1
                target = target_lo
                if target not in decoded:
                    decoded.add(target)
                    calls += 1
                    for secret in core.decode(public, target):
                        answers[core.key(secret)] = secret
                else:
                    repeats += 1
            else:
                coordinate = choose_split(public, prep, node, box, ambiguous)
                lo, hi = node[coordinate]
                middle = (lo + hi) // 2
                left, right = list(node), list(node)
                left[coordinate], right[coordinate] = (lo, middle), (middle + 1, hi)
                stack.extend([tuple(right), tuple(left)])
        checkpoint()
    checkpoint(force=True)
    complete = not stack
    if complete:
        assert covered == summary["planned_cells"]
    return {"status": "COMPLETE" if complete else "RESOURCE_LIMITED_PARTIAL",
            "complete_search": complete, "outputs": core.unique(answers.values()),
            "visited_outer_iterations": visited, "visited_nodes": visited,
            "visited_cells": leaves, "logical_cells_covered": covered,
            "rejected_subtrees": rejected, "singleton_quotient_subtrees": certified,
            "empty_target_interval_subtrees": target_empty,
            "singleton_target_interval_subtrees": target_singleton,
            "tested_target_proposals": calls, "pending_frontier_nodes": len(stack),
            "tested_target_values": sorted(decoded),
            "repeated_target_proposals_avoided": repeats,
            "maximum_nodes": maximum_nodes, "original_interval_pruning": tighten_original,
            "effective_node_budget": node_budget, "target_switch_budget": target_switch_budget,
            "termination_reason": "FRONTIER_EXHAUSTED" if complete else
                "TARGET_BUDGET_REACHED" if target_switch_budget is not None and visited >= target_switch_budget else "NODE_RESOURCE_CAP_REACHED",
            "exhaustive_target_scan_used": False, "first_answer_early_stop": False}


def solve_cells(public, backend, event, certificate_path, maximum_nodes=10000,
                tighten_original=True, reuse=None, maximum_targets=10000):
    if public["family"] == "echnp" and None in public["reports"]:
        direct = literal_solve(public, backend, event, certificate_path)
        direct.update(complete_search=True, algorithm_arm="adaptive_cell_tree",
                      infinity_exception=True, visited_nodes=0,
                      exhaustive_target_scan_used=False, first_answer_early_stop=False)
        direct["preparation"].update(algorithm_arm="adaptive_cell_tree", literal_algorithm1=False)
        return direct
    prep = prepare_cells(public, backend, event, certificate_path, reuse)
    prep["summary"].update(maximum_nodes=maximum_nodes, original_interval_pruning=tighten_original,
                           algorithm_arm="adaptive_cell_tree", literal_algorithm1=False)
    targets = prep["summary"]["planned_target_proposals"]
    prep["summary"]["target_switch_budget"] = targets
    prep["summary"]["planned_iterations"] = min(prep["summary"]["tree_worst_case_nodes"], maximum_nodes, targets)
    prep["summary"]["planned_iterations_kind"] = "committed_node_budget_not_expected_completion_cost"
    event(phase="preparation_finished", preparation=prep["summary"])
    start = monotonic()
    result = recover_cells(public, prep, event, maximum_nodes, tighten_original, targets)
    result["target_fallback_limit"] = maximum_targets
    if result["termination_reason"] == "TARGET_BUDGET_REACHED" and targets <= maximum_targets:
        # Only after trying C; never turn a huge target space into a new long job.
        event(phase="bounded_target_fallback", target_fallback_planned=prep["summary"]["planned_target_proposals"])
        answers = {core.key(secret): secret for secret in result["outputs"]}
        visited = skipped = positions = 0
        decoded = set(result["tested_target_values"])
        for target in range(prep["sources"].box[0][0], prep["sources"].box[0][1] + 1):
            positions += 1
            if target in decoded:
                skipped += 1
            else:
                visited += 1
                decoded.add(target)
                for secret in core.decode(public, target):
                    answers[core.key(secret)] = secret
            if positions % 256 == 0:
                event(phase="bounded_target_fallback", target_fallback_visited=visited,
                      partial_outputs=list(answers.values()))
        result.update(status="COMPLETE", complete_search=True, outputs=core.unique(answers.values()),
                      cell_tree_termination_reason=result["termination_reason"],
                      termination_reason="BOUNDED_TARGET_FALLBACK_EXHAUSTED",
                      exhaustive_target_scan_used=True, target_fallback_visited=visited,
                      target_fallback_range_positions=positions,
                      target_fallback_skipped_already_tested=skipped,
                      tested_target_values=sorted(decoded),
                      tested_target_proposals=result["tested_target_proposals"] + visited,
                      search_completion_certificate="exhaustive_bounded_target_fallback")
        prep["summary"]["search_mode"] = "adaptive_cells_then_bounded_target_fallback"
    else:
        result["search_completion_certificate"] = "certified_cell_frontier_exhausted" if result["complete_search"] else None
        result["target_fallback_not_run_reason"] = (
            "cell_search_completed" if result["complete_search"] else
            "node_resource_cap_reached_before_target_route_budget" if result["termination_reason"] == "NODE_RESOURCE_CAP_REACHED" else
            "target_count_exceeds_public_work_cap")
    prep["summary"]["phase_seconds"]["recovery"] = monotonic() - start
    return {**result, "preparation": prep["summary"], "infinity_observations": 0}
