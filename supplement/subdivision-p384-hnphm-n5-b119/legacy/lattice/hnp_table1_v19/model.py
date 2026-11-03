"""Literal source-(19) ordering, declared boxes, and uncapped Algorithm 1.

The historical source-(9) arithmetic is reused without modifying it. No
generator opening, planted secret, or observation seed enters this module.
"""
from __future__ import annotations

import importlib.util
import json
from itertools import product
from math import prod
from pathlib import Path
import sys
import subprocess
from time import perf_counter

from sage.all import ZZ, matrix

CORE_PATH = Path(__file__).resolve().parents[1] / "hnp_threshold_experiments/hnp.py"
MI_BINARY = Path(__file__).resolve().parent / "bin/mi_targets"
spec = importlib.util.spec_from_file_location("hnp_source9_arithmetic", CORE_PATH)
core = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = core
spec.loader.exec_module(core)


def sources19(public):
    old = core.construct_sources(public)
    order = [i for i in range(len(old.polynomials)) if i != old.target_row]
    order.append(old.target_row)
    box = old.box
    if public["family"] == "hnphm":
        box = tuple((-public["X"], public["X"]) for _ in old.box)
    return core.Sources(tuple(old.polynomials[i] for i in order), old.edges,
                        len(order) - 1, box, old.scale)


def prepare19(public, event=lambda **kw: None):
    started = perf_counter()
    sources = sources19(public)
    basis = sources.scaled_matrix()
    timings = {"source": perf_counter() - started}
    event(stage="lll_started", rank=basis.nrows(), ambient=basis.ncols())
    start = perf_counter()
    reduced, transform = core.tracked_lll(basis)
    timings["lll"] = perf_counter() - start
    event(stage="lll_finished", seconds=timings["lll"])
    start = perf_counter()
    alpha = list(map(int, transform.inverse().change_ring(ZZ).row(sources.target_row)))
    coefficient = transform * sources.coefficient_matrix()
    constants = transform * matrix(ZZ, [[f.constant] for f in sources.polynomials])
    qcount = len(sources.edges)
    polynomials = tuple(core.Polynomial(
        int(constants[i, 0]), tuple(map(int, coefficient.row(i)[qcount:])),
        tuple(map(int, coefficient.row(i)[:qcount]))) for i in range(basis.nrows()))
    relevant = [i for i, a in enumerate(alpha) if a]
    width_norm = max(sum(abs(int(c)) for c in reduced.row(i)) for i in relevant)
    subdivisions = max(1, core.ceildiv(12 * width_norm, public["p"]))
    ncells = prod(min(subdivisions, b - a + 1) for a, b in sources.box)
    count_target = sources.box[0][1] - sources.box[0][0] + 1
    timings["readout"] = perf_counter() - start
    event(stage="audit_started", search_count=min(ncells, count_target))
    start = perf_counter()
    assert transform * basis == reduced
    assert abs(transform.det()) == 1
    assert matrix(ZZ, [alpha]) * transform == matrix(
        ZZ, [[int(j == sources.target_row) for j in range(basis.nrows())]])
    certification = core.exact_lll_audit(reduced)
    timings["exact_audit"] = perf_counter() - start
    return {"sources": sources, "basis": basis, "reduced": reduced,
            "transform": transform, "alpha": alpha, "relevant": relevant,
            "polynomials": polynomials, "W": width_norm, "b": subdivisions,
            "cell_count": ncells, "target_count": count_target,
            "search_count": min(ncells, count_target),
            "search_mode": "target" if count_target <= ncells else "cells",
            "timing": timings, "certification": certification}


def recover19(public, prep, event=lambda **kw: None):
    """All cells or all declared targets, without time/work limits or early stop."""
    if public["family"] == "echnp" and None in public["reports"]:
        return core.recover(public)
    start = perf_counter()
    answers, seen = [], set()
    visited = calls = rejected = 0
    lower_target, upper_target = prep["sources"].box[0]
    if public["family"] == "mihnp" and prep["search_mode"] == "target":
        return mi_targets(public, event)
    if prep["search_mode"] == "target":
        for target in range(lower_target, upper_target + 1):
            answers.extend(core.decode(public, target))
            visited += 1
            calls += 1
            if visited % 100000 == 0:
                event(stage="search", visited=visited, planned=prep["search_count"],
                      seconds=perf_counter() - start)
    else:
        parts = [core.balanced_parts(a, b, prep["b"]) for a, b in prep["sources"].box]
        for cell in product(*parts):
            visited += 1
            target, possible = 0, True
            for j in prep["relevant"]:
                lo, hi = prep["polynomials"][j].interval(cell, prep["sources"].edges)
                lower, upper = core.ceildiv(lo, public["p"]), hi // public["p"]
                if lower > upper:
                    possible = False
                    rejected += 1
                    break
                if lower != upper:
                    raise ArithmeticError("Cell did not isolate its relevant quotient")
                target += prep["alpha"][j] * lower
            if possible and lower_target <= target <= upper_target:
                # Caching identical public reconstruction calls preserves the list.
                if target not in seen:
                    seen.add(target)
                    answers.extend(core.decode(public, target))
                    calls += 1
            if visited % 10000 == 0:
                event(stage="search", visited=visited, planned=prep["search_count"],
                      seconds=perf_counter() - start)
    return {"status": "COMPLETE", "secrets": core.unique(answers),
            "visited": visited, "decoder_calls": calls, "cells_rejected": rejected,
            "search_seconds": perf_counter() - start}


def mi_targets(public, event=lambda **kw: None):
    """Exact batch inversion for the full declared MI target loop, no early stop."""
    data = " ".join(map(str, [public["p"], public["N"], public["X"],
                              *public["reports"], *public["shifts"]])) + "\n"
    process = subprocess.Popen([str(MI_BINARY)], stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    process.stdin.write(data)
    process.stdin.close()
    result = None
    for line in process.stdout:
        row = json.loads(line)
        if row["stage"] == "native_result":
            result = row
        else:
            event(**row)
    code = process.wait()  # No timeout.
    error = process.stderr.read()
    process.stdout.close()
    process.stderr.close()
    if code or result is None:
        raise RuntimeError(f"MI batch reconstruction failed: {code}: {error}")
    assert result["visited"] == result["decoder_calls"] == public["X"]
    for secret in result["secrets"]:
        target = pow((secret + public["shifts"][0]) % public["p"], -1, public["p"]) - public["reports"][0]
        assert secret in core.decode(public, target)
    result.pop("stage")
    result["secrets"] = core.unique(result["secrets"])
    result["decoder_backend"] = "exact batched modular inversions; every target visited"
    return result
