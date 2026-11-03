"""Public-only complete Algorithm1 with flatter and exact tracked NTL cleanup."""
from __future__ import annotations

import gzip
import json
from math import prod
import os
from pathlib import Path
import re
import subprocess
import sys
from time import monotonic

from sage.all import QQ, ZZ, matrix, vector

ADAPTER_DIR = Path(__file__).resolve().parent.parent / "hnp_table1_v19"
sys.path.insert(0, str(ADAPTER_DIR))
from model import core, sources19

FLATTER = Path("/home/pqsnc/.local/bin/flatter")
FLATTER_LIBRARY = Path("/home/pqsnc/.local/lib/libflatter.so")


def matrix_text(basis):
    return "[" + "\n".join("[" + " ".join(map(str, row)) + "]" for row in basis.rows()) + "]\n"


def parse_matrix(text, rows, columns):
    parsed = [list(map(int, row.split())) for row in re.findall(r"\[([^\[\]]*)\]", text)]
    if len(parsed) != rows or any(len(row) != columns for row in parsed):
        raise ValueError("Flatter returned an unexpected matrix shape")
    return matrix(ZZ, parsed)


def recover_transform(original, reduced):
    """All sources have independent first m columns, including rectangular HN."""
    m = original.nrows()
    pivot = original.matrix_from_columns(range(m))
    target = reduced.matrix_from_columns(range(m))
    transform = (target.change_ring(QQ) * pivot.change_ring(QQ).inverse()).change_ring(ZZ)
    if transform * original != reduced or abs(transform.det()) != 1:
        raise ArithmeticError("Reducer did not preserve the exact row lattice")
    return transform


def reduce_basis(basis, backend, event):
    timing = {}
    if backend == "flatter-ntl":
        event(phase="flatter")
        env = {k: v for k, v in os.environ.items() if not k.startswith("FLATTER_")}
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        start = monotonic()
        process = subprocess.run([str(FLATTER)], input=matrix_text(basis),
                                 text=True, capture_output=True, check=True, env=env)
        reduced = parse_matrix(process.stdout, basis.nrows(), basis.ncols())
        timing["flatter"] = monotonic() - start
        event(phase="transformation_recovery")
        start = monotonic()
        first_transform = recover_transform(basis, reduced)
        timing["transformation_recovery"] = monotonic() - start
    elif backend == "ntl-exact":
        reduced = basis
        first_transform = None
        timing.update(flatter=0.0, transformation_recovery=0.0)
    else:
        raise ValueError("Unknown reduction backend")
    event(phase="exact_ntl_lll")
    start = monotonic()
    reduced, cleanup = reduced.LLL(delta=QQ(3) / 4, algorithm="NTL:LLL", fp=None,
                                   transformation=True)
    timing["exact_ntl_lll"] = monotonic() - start
    start = monotonic()
    transform = cleanup if first_transform is None else cleanup * first_transform
    assert transform * basis == reduced
    assert abs(transform.det()) == 1
    timing["transformation_check"] = monotonic() - start
    return reduced, transform, timing


def prepare(public, backend, event, certificate_path=None):
    event(phase="source_construction")
    start = monotonic()
    sources = sources19(public)
    basis = sources.scaled_matrix()
    n, family = public["N"], public["family"]
    expected_rank = 2 * n + 1 if family == "hnphm" else n * (n + 1) // 2
    expected_ambient = 3 * n + 1 if family == "hnphm" else expected_rank
    assert (basis.nrows(), basis.ncols()) == (expected_rank, expected_ambient)
    assert sources.target_row == expected_rank - 1
    target = sources.polynomials[-1]
    assert target.constant == 0 and not any(target.quadratic)
    assert target.linear == (public["p"],) + (0,) * (len(target.linear) - 1)
    timing = {"source_construction": monotonic() - start}
    event(phase="source_ready", source_rank=expected_rank, ambient_dimension=expected_ambient)
    reduced, transform, reduction_timing = reduce_basis(basis, backend, event)
    timing.update(reduction_timing)
    event(phase="readout")
    start = monotonic()
    unit = vector(ZZ, [0] * (expected_rank - 1) + [1])
    alpha = list(map(int, transform.transpose().solve_right(unit).change_ring(ZZ)))
    assert vector(ZZ, alpha) * transform == unit
    coefficient = transform * sources.coefficient_matrix()
    constants = transform * vector(ZZ, [f.constant for f in sources.polynomials])
    quadratic_count = len(sources.edges)
    polynomials = tuple(core.Polynomial(
        int(constants[i]), tuple(map(int, coefficient.row(i)[quadratic_count:])),
        tuple(map(int, coefficient.row(i)[:quadratic_count]))) for i in range(expected_rank))
    relevant = [i for i, a in enumerate(alpha) if a]
    width_norm = max(sum(abs(int(c)) for c in reduced.row(i)) for i in relevant)
    subdivisions = max(1, core.ceildiv(12 * width_norm, public["p"]))
    planned_cells = prod(min(subdivisions, upper - lower + 1) for lower, upper in sources.box)
    targets = sources.box[0][1] - sources.box[0][0] + 1
    timing["readout"] = monotonic() - start
    if certificate_path is not None:
        event(phase="certificate_write")
        start = monotonic()
        with gzip.open(certificate_path, "xt") as handle:
            json.dump({"transform": [list(map(int, row)) for row in transform.rows()],
                       "alpha": alpha}, handle, separators=(",", ":"))
        timing["certificate_write"] = monotonic() - start
    summary = {"source_rank": expected_rank, "ambient_dimension": expected_ambient,
               "relevant_rows": relevant, "W": width_norm, "W_over_p": str(QQ(width_norm) / public["p"]),
               "subdivisions": subdivisions, "error_interval_lengths": [hi - lo + 1 for lo, hi in sources.box],
               "planned_cells": planned_cells, "planned_target_proposals": targets,
               "planned_iterations": min(planned_cells, targets),
               "search_mode": "target_scan" if targets <= planned_cells else "cell_search",
               "phase_seconds": timing, "backend": backend,
               "basis_change_exact": True, "target_readout_exact": True,
               "lll_contract": "NTL exact LLL, delta3/4, size1/2; trusted library contract",
               "independent_large_gso_audit": False}
    return {"sources": sources, "reduced": reduced, "transform": transform,
            "alpha": alpha, "polynomials": polynomials, "summary": summary}


def recover(public, prep, event):
    summary, sources = prep["summary"], prep["sources"]
    lower_target, upper_target = sources.box[0]
    answers = {}
    visited = calls = rejected = 0
    last_checkpoint = monotonic()

    def checkpoint(force=False):
        nonlocal last_checkpoint
        now = monotonic()
        if force or now - last_checkpoint >= 1:
            event(phase="recovery", visited_outer_iterations=visited,
                  visited_cells=visited if summary["search_mode"] == "cell_search" else 0,
                  tested_target_proposals=calls, rejected_cells=rejected,
                  partial_outputs=list(answers.values()))
            last_checkpoint = now

    def reconstruct(target):
        nonlocal calls
        calls += 1
        for secret in core.decode(public, target):
            answers[core.key(secret)] = secret

    if summary["search_mode"] == "target_scan":
        for target in range(lower_target, upper_target + 1):
            visited += 1
            reconstruct(target)
            if visited % 256 == 0:
                checkpoint()
    else:
        from itertools import product
        parts = [core.balanced_parts(lo, hi, summary["subdivisions"]) for lo, hi in sources.box]
        for cell in product(*parts):
            visited += 1
            target, possible = 0, True
            for j in summary["relevant_rows"]:
                lo, hi = prep["polynomials"][j].interval(cell, sources.edges)
                low, high = core.ceildiv(lo, public["p"]), hi // public["p"]
                if low > high:
                    possible = False
                    rejected += 1
                    break
                if low != high:
                    raise ArithmeticError("Cell interval did not isolate a quotient")
                target += prep["alpha"][j] * low
            if possible and lower_target <= target <= upper_target:
                reconstruct(target)
            checkpoint()
    checkpoint(force=True)
    assert visited == summary["planned_iterations"]
    return {"status": "COMPLETE", "outputs": core.unique(answers.values()),
            "visited_outer_iterations": visited,
            "visited_cells": visited if summary["search_mode"] == "cell_search" else 0,
            "tested_target_proposals": calls, "rejected_cells": rejected}


def solve(public, backend, event, certificate_path=None):
    if public["family"] == "echnp" and None in public["reports"]:
        event(phase="infinity_reconstruction")
        direct = core.recover(public)
        return {"status": "COMPLETE", "outputs": direct["secrets"], "infinity_observations": public["reports"].count(None),
                "visited_outer_iterations": 0, "visited_cells": 0, "tested_target_proposals": 0,
                "preparation": {"search_mode": "infinity", "planned_iterations": 0}}
    prep = prepare(public, backend, event, certificate_path)
    event(phase="preparation_finished", preparation=prep["summary"])
    start = monotonic()
    result = recover(public, prep, event)
    prep["summary"]["phase_seconds"]["recovery"] = monotonic() - start
    return {**result, "infinity_observations": 0, "preparation": prep["summary"]}
