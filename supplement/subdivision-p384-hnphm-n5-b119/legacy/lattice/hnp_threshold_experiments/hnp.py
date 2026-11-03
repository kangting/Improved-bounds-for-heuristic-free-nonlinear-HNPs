"""Exact sources, target decoders and Algorithm 1 of the supplied HNP PDF.

Run with Sage's Python. Generation and recovery have separate interfaces;
the recovery routines receive no planted secret or generator seed.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from math import prod
import random
from time import perf_counter

from sage.all import EllipticCurve, GF, QQ, ZZ, matrix, set_random_seed


def ceildiv(a: int, b: int) -> int:
    return -(-a // b)


def floor_sum(n: int, m: int, a: int, b: int) -> int:
    """Sum floor((a*i+b)/m), 0<=i<n, including signed a and b.

    Euclidean recursion, O(log(max(m,abs(a)))) arithmetic steps; no scan of
    the hidden-multiplier interval. Appendix D of the supplied paper.
    """
    if n < 0 or m <= 0:
        raise ValueError("Require n>=0 and m>0")
    ans = 0
    while True:
        qa, a = divmod(a, m)
        qb, b = divmod(b, m)
        ans += qa * n * (n - 1) // 2 + qb * n
        ymax = a * n + b
        if ymax < m:
            return ans
        n, b, m, a = ymax // m, ymax % m, a, m


def modular_interval_count(y: int, p: int, lo: int, hi: int,
                           lower: int, upper: int) -> int:
    """Number of t in [lo,hi] with lower <= y*t mod p <= upper."""
    if lo > hi or lower > upper:
        return 0
    n = hi - lo + 1
    return (floor_sum(n, p, y, y * lo - lower)
            - floor_sum(n, p, y, y * lo - upper - 1))


@dataclass(frozen=True)
class Polynomial:
    constant: int
    linear: tuple[int, ...]
    quadratic: tuple[int, ...]

    def evaluate(self, x, edges) -> int:
        return (self.constant + sum(c * t for c, t in zip(self.linear, x))
                + sum(c * x[i] * x[j]
                      for c, (i, j) in zip(self.quadratic, edges)))

    def interval(self, box, edges) -> tuple[int, int]:
        lo = hi = self.constant
        for c, (a, b) in zip(self.linear, box):
            lo += min(c * a, c * b)
            hi += max(c * a, c * b)
        for c, (i, j) in zip(self.quadratic, edges):
            if c:
                values = [c * a * b for a in box[i] for b in box[j]]
                lo += min(values)
                hi += max(values)
        return lo, hi


@dataclass(frozen=True)
class Sources:
    polynomials: tuple[Polynomial, ...]
    edges: tuple[tuple[int, int], ...]
    target_row: int
    box: tuple[tuple[int, int], ...]
    scale: int

    def coefficient_matrix(self):
        return matrix(ZZ, [list(f.quadratic) + list(f.linear)
                           for f in self.polynomials])

    def scaled_matrix(self):
        return matrix(ZZ, [[self.scale ** 2 * c for c in f.quadratic]
                           + [self.scale * c for c in f.linear]
                           for f in self.polynomials])


def curve(public):
    return EllipticCurve(GF(public["p"]), public["curve"])


def point(ec, encoded):
    return ec(0) if encoded is None else ec(encoded)


def encode_point(q):
    return None if q.is_zero() else [int(q[0]), int(q[1])]


def construct_sources(public) -> Sources:
    p, xbound = int(public["p"]), int(public["X"])
    n, family = public["N"], public["family"]
    reports = public["reports"]
    if family == "hnphm":
        k = 2 * n + 1
        edges = tuple((0, i + 1) for i in range(n))
    else:
        k = n
        edges = tuple(combinations(range(n), 2))
    edge_index = {e: j for j, e in enumerate(edges)}
    polynomials = []

    def add(constant=0, linear=None, quadratic=None):
        polynomials.append(Polynomial(int(constant), tuple(linear or [0] * k),
                                      tuple(quadratic or [0] * len(edges))))

    if family == "mihnp":
        shifts = public["shifts"]
        for i, j in edges:
            a = pow((shifts[i] - shifts[j]) % p, -1, p)
            if 2 * a > p:
                a -= p
            v, w = reports[i], reports[j]
            lin, quad = [0] * k, [0] * len(edges)
            lin[i], lin[j] = w + a, v - a
            quad[edge_index[i, j]] = 1
            add(v * w + a * (v - w), lin, quad)
    elif family == "echnp":
        ec = curve(public)
        queries = [point(ec, q) for q in public["queries"]]
        if any(q.is_zero() for q in queries):
            raise ValueError("EC queries must be nonzero")
        coords, aedge, bedge = [int(q[0]) for q in queries], {}, {}
        for i, j in edges:
            d = queries[j] - queries[i]
            if (2 * d).is_zero():
                raise ValueError("EC query difference lies in E[2]")
            aedge[i, j], bedge[i, j] = int(d[0]), int(d[1])
        star = n - 1
        for i, j in combinations(range(star), 2):
            weights = {(i, j): 1,
                       (i, star): (-bedge[i, j]
                                   * pow(bedge[i, star], -1, p)) % p,
                       (j, star): (bedge[i, j]
                                   * pow(bedge[j, star], -1, p)) % p}
            lin, quad = [0] * k, [0] * len(edges)
            kappa = sum(c * (coords[u] * coords[v]
                             - aedge[u, v] * (coords[u] + coords[v]))
                        for (u, v), c in weights.items()) % p
            constant = -kappa
            for (u, v), c in weights.items():
                a = aedge[u, v]
                quad[edge_index[u, v]] = c
                lin[u] += c * (reports[v] - a)
                lin[v] += c * (reports[u] - a)
                constant += c * (reports[u] * reports[v]
                                 - a * (reports[u] + reports[v]))
            add(constant, lin, quad)
        for i in range(star):
            quad = [0] * len(edges)
            quad[edge_index[i, star]] = p
            add(quadratic=quad)
    elif family == "hnphm":
        for i in range(n):
            lin, quad = [0] * k, [0] * len(edges)
            lin[0], lin[i + 1], lin[n + i + 1] = reports[i + 1], reports[0], -1
            quad[i] = 1
            add(reports[0] * reports[i + 1] - reports[n + i + 1], lin, quad)
    else:
        raise ValueError(f"Unknown family: {family}")
    target_row = len(polynomials)
    for i in range(n + 1 if family == "hnphm" else n):
        lin = [0] * k
        lin[i] = p
        add(linear=lin)
    box = tuple((0, xbound - 1) for _ in range(k))
    if public.get("report_model") == "radius":
        if family != "hnphm":
            raise ValueError("Radius reports are implemented only for HNPHM")
        box = tuple((-xbound, xbound) for _ in range(k))
    return Sources(tuple(polynomials), edges, target_row, box, xbound)


def target_interval(public):
    if public.get("report_model") == "radius":
        return (max(-public["X"], 1 - public["reports"][0]),
                min(public["X"], public["p"] - 1 - public["reports"][0]))
    return (0, min(public["X"] - 1, public["p"] - 1 - public["reports"][0]))


def decode(public, target):
    """Exact original-input verifier, no planted information."""
    lo, hi = target_interval(public)
    if not lo <= target <= hi:
        return []
    p, xbound, n = public["p"], public["X"], public["N"]
    reports, family = public["reports"], public["family"]
    value = reports[0] + target
    if family == "mihnp":
        if not 0 < value < p:
            return []
        secret = (pow(value, -1, p) - public["shifts"][0]) % p
        for shift, report in zip(public["shifts"], reports):
            if (secret + shift) % p == 0:
                return []
            completed = pow((secret + shift) % p, -1, p)
            if not report <= completed < min(p, report + xbound):
                return []
        return [secret]
    if family == "hnphm":
        if not 0 < value < p:
            return []
        for i in range(n):
            radius = public.get("report_model") == "radius"
            lower_error, upper_error = (-xbound, xbound) if radius else (0, xbound - 1)
            left = max(1, reports[i + 1] + lower_error)
            right = min(p - 1, reports[i + 1] + upper_error)
            lower = max(0, reports[n + i + 1] + lower_error)
            upper = min(p - 1, reports[n + i + 1] + upper_error)
            if modular_interval_count(value, p, left, right, lower, upper) <= 0:
                return []
        return [value]
    ec = curve(public)
    if not 0 <= value < p:
        return []
    queries = [point(ec, q) for q in public["queries"]]
    out = []
    for q in ec.lift_x(GF(p)(value), all=True):
        secret = q - queries[0]
        good = True
        for query, report in zip(queries, reports):
            shifted = secret + query
            if report is None:
                good = shifted.is_zero()
            else:
                good = (not shifted.is_zero()
                        and report <= int(shifted[0]) < min(p, report + xbound))
            if not good:
                break
        if good:
            out.append(encode_point(secret))
    return out


def key(value):
    return tuple(value) if isinstance(value, list) else value


def unique(values):
    return sorted({key(v): v for v in values}.values(), key=str)


def exact_lll_audit(reduced):
    """Exact rational standard (3/4,1/2)-LLL inequalities."""
    gs, mu = reduced.change_ring(QQ).gram_schmidt()
    norms = [row.dot_product(row) for row in gs.rows()]
    for i in range(reduced.nrows()):
        if norms[i] <= 0:
            raise ArithmeticError("Rank-deficient basis")
        if any(abs(mu[i, j]) > QQ(1) / 2 for j in range(i)):
            raise ArithmeticError("Basis is not exactly size reduced")
        if i and norms[i] < (QQ(3) / 4 - mu[i, i - 1] ** 2) * norms[i - 1]:
            raise ArithmeticError("Exact Lovasz condition failed")
    return {"size_reduced_half": True, "lovasz_three_quarters": True}


def tracked_lll(basis):
    """Robust floating reducer followed by exact size reduction to 1/2.

    fpLLL can fail at the exact floating boundary eta=0.5. A stronger Lovasz
    request plus rational row operations gives the paper's standard LLL
    contract; prepare() independently certifies it over QQ.
    """
    reduced, transform = basis.LLL(delta=0.99, eta=0.501, transformation=True)
    gs, mus = reduced.change_ring(QQ).gram_schmidt()
    if all(abs(mus[i, j]) <= QQ(1) / 2 for i in range(reduced.nrows()) for j in range(i)):
        return reduced, transform
    norms = [row.dot_product(row) for row in gs.rows()]
    for i in range(reduced.nrows()):
        for j in reversed(range(i)):
            c = reduced.row(i).dot_product(gs.row(j)) / norms[j]
            if abs(c) > QQ(1) / 2:
                nearest = (c + QQ(1) / 2).floor()
                reduced.add_multiple_of_row(i, j, -nearest)
                transform.add_multiple_of_row(i, j, -nearest)
    return reduced, transform


def prepare(public, audit=True):
    start = perf_counter()
    sources = construct_sources(public)
    basis = sources.scaled_matrix()
    source_seconds = perf_counter() - start
    start = perf_counter()
    reduced, transform = tracked_lll(basis)
    lll_seconds = perf_counter() - start
    start = perf_counter()
    inverse = transform.inverse().change_ring(ZZ)
    alpha = list(map(int, inverse.row(sources.target_row)))
    coefficient = transform * sources.coefficient_matrix()
    constants = transform * matrix(ZZ, [[f.constant] for f in sources.polynomials])
    qcount = len(sources.edges)
    polynomials = tuple(Polynomial(int(constants[i, 0]),
                                   tuple(map(int, coefficient.row(i)[qcount:])),
                                   tuple(map(int, coefficient.row(i)[:qcount])))
                        for i in range(basis.nrows()))
    relevant = [i for i, a in enumerate(alpha) if a]
    width_norm = max(sum(abs(int(c)) for c in reduced.row(i)) for i in relevant)
    subdivisions = max(1, ceildiv(12 * width_norm, public["p"]))
    ncells = prod(min(subdivisions, b - a + 1) for a, b in sources.box)
    target_lo, target_hi = target_interval(public)
    count_target = max(0, target_hi - target_lo + 1)
    readout_seconds = perf_counter() - start
    start = perf_counter()
    if audit:
        assert transform * basis == reduced
        assert abs(transform.det()) == 1
        assert matrix(ZZ, [alpha]) * transform == matrix(
            ZZ, [[int(j == sources.target_row) for j in range(basis.nrows())]])
        certification = exact_lll_audit(reduced)
    else:
        certification = None
    return {"sources": sources, "basis": basis, "reduced": reduced,
            "transform": transform, "alpha": alpha, "relevant": relevant,
            "polynomials": polynomials, "W": width_norm, "b": subdivisions,
            "cell_count": ncells, "target_count": count_target,
            "search_count": min(ncells, count_target),
            "search_mode": "target" if count_target <= ncells else "cells",
            "timing": {"source": source_seconds, "lll": lll_seconds,
                       "readout": readout_seconds,
                       "exact_audit": perf_counter() - start},
            "certification": certification}


def balanced_parts(lo, hi, b):
    n = hi - lo + 1
    b = min(n, b)
    return [(lo + i * n // b, lo + (i + 1) * n // b - 1) for i in range(b)]


def recover(public, prepared=None, max_work=None):
    """Algorithm 1. A cap returns CENSORED without an incomplete secret list."""
    if public["family"] == "echnp" and None in public["reports"]:
        idx = public["reports"].index(None)
        ec = curve(public)
        secret = -point(ec, public["queries"][idx])
        for q, v in zip(public["queries"], public["reports"]):
            shifted = secret + point(ec, q)
            if ((v is None and not shifted.is_zero())
                    or (v is not None and (shifted.is_zero()
                        or not v <= int(shifted[0]) < min(public["p"], v + public["X"])))):
                return {"status": "COMPLETE", "secrets": [], "visited": 0,
                        "decoder_calls": 0, "search_seconds": 0.0}
        return {"status": "COMPLETE", "secrets": [encode_point(secret)],
                "visited": 0, "decoder_calls": 0, "search_seconds": 0.0}
    prep = prepared if prepared is not None else prepare(public)
    if max_work is not None and prep["search_count"] > max_work:
        return {"status": "CENSORED_WORK", "secrets": None, "visited": 0,
                "decoder_calls": 0, "search_seconds": 0.0}
    start = perf_counter()
    answers, seen = [], set()
    visited = calls = rejected = 0
    if prep["search_mode"] == "target":
        proposals = range(*[target_interval(public)[0], target_interval(public)[1] + 1])
        for target in proposals:
            answers.extend(decode(public, target))
            visited += 1
            calls += 1
    else:
        parts = [balanced_parts(a, b, prep["b"]) for a, b in prep["sources"].box]
        for cell in product(*parts):
            visited += 1
            target, possible = 0, True
            for j in prep["relevant"]:
                lo, hi = prep["polynomials"][j].interval(cell, prep["sources"].edges)
                lower, upper = ceildiv(lo, public["p"]), hi // public["p"]
                if lower > upper:
                    possible = False
                    rejected += 1
                    break
                if lower != upper:
                    raise ArithmeticError("Cell did not isolate the relevant quotient")
                target += prep["alpha"][j] * lower
            if possible and target not in seen:
                seen.add(target)
                answers.extend(decode(public, target))
                calls += 1
    return {"status": "COMPLETE", "secrets": unique(answers), "visited": visited,
            "decoder_calls": calls, "cells_rejected": rejected,
            "search_seconds": perf_counter() - start}


def singleton_ablation(public, prepared):
    """Same-source ablation, NOT a published prior algorithm."""
    target = 0
    for j in prepared["relevant"]:
        lo, hi = prepared["polynomials"][j].interval(
            prepared["sources"].box, prepared["sources"].edges)
        lower, upper = ceildiv(lo, public["p"]), hi // public["p"]
        if lower != upper:
            return {"status": "NOT_SINGLETON", "secrets": []}
        target += prepared["alpha"][j] * lower
    return {"status": "SINGLETON", "secrets": decode(public, target)}


def recover_singleton_first(public, prepared=None, max_work=None):
    """Exact implementation shortcut, then unchanged paper Algorithm 1.

    All witnesses have the same relevant quotients when the whole-box
    intervals are singletons, hence the same target. Its exact decoder
    therefore already returns the complete list. The added polynomial check
    preserves the average-time theorem; fallback wall time can gain overhead.
    """
    if public["family"] == "echnp" and None in public["reports"]:
        return {**recover(public), "shortcut_used": False,
                "effective_search_mode": "infinity"}
    prep = prepared if prepared is not None else prepare(public)
    start = perf_counter()
    whole = singleton_ablation(public, prep)
    if whole["status"] == "SINGLETON":
        return {"status": "COMPLETE", "secrets": unique(whole["secrets"]),
                "shortcut_used": True, "effective_search_mode": "whole_box_singleton",
                "effective_search_work": 1, "visited": 0, "decoder_calls": 1,
                "search_seconds": perf_counter() - start}
    result = recover(public, prep, max_work=max_work)
    result["search_seconds"] = perf_counter() - start
    return {**result, "shortcut_used": False,
            "effective_search_mode": prep["search_mode"],
            "effective_search_work": prep["search_count"]}


def finite_exponent(family: str, n: int) -> Fraction:
    if family == "mihnp":
        k = 0
        while (k + 1) * (k + 2) <= n:
            k += 1
        return Fraction(n + k * k - 1, 2 * (n + k * (k + 1)))
    if family == "echnp":
        return Fraction(n - 2, 5 * n - 6)
    return Fraction(n, 3 * n + 1)


def generate(family, p, n, unknown_bits, observation_seed, secret,
             ec_parameters=None, observation_model="independent"):
    """Public bucket reports and a separate private audit opening.

    The true secret is fixed by the caller across observation seeds. EC
    query tuples are sampled iid then rejected as whole tuples, preserving
    exactly the conditioned distribution of Definition 5.1.
    """
    rng = random.Random(observation_seed)
    xbound = 1 << unknown_bits
    if not 1 <= xbound < p:
        raise ValueError("Bucket width must be less than p")
    public = {"family": family, "p": int(p), "N": n, "X": xbound,
              "report_model": "leading_bits", "observation_model": observation_model}
    private = {"secret": secret, "observation_seed": observation_seed}
    if family == "mihnp":
        shifts, used = [], set()
        while len(shifts) < n:
            t = rng.randrange(p)
            if t not in used and (secret + t) % p:
                shifts.append(t)
                used.add(t)
        values = [pow((secret + t) % p, -1, p) for t in shifts]
        public["shifts"] = shifts
    elif family == "hnphm":
        multipliers = [rng.randrange(1, p) for _ in range(n)]
        values = [secret] + multipliers + [secret * t % p for t in multipliers]
        private["multipliers"] = multipliers
    elif family == "echnp":
        if ec_parameters is None:
            raise ValueError("Need a known-order subgroup")
        public.update(ec_parameters)
        ec = curve(public)
        generator = point(ec, public["generator"])
        if generator.order() != public["M"]:
            raise ValueError("Incorrect subgroup order")
        fixed_secret = point(ec, secret)
        rejected = 0
        while True:
            if observation_model == "paired":
                if n % 2 != 1:
                    raise ValueError("Paired EC observations require odd N")
                # A fixed nonzero base translates the published 0,+t,-t
                # oracle model while keeping the new source constructible.
                base_scalar = 3
                samples = [rng.randrange(public["M"]) for _ in range((n - 1) // 2)]
                scalars = [base_scalar] + [s for t in samples
                                          for s in ((base_scalar + t) % public["M"],
                                                    (base_scalar - t) % public["M"])]
            else:
                scalars = [rng.randrange(public["M"]) for _ in range(n)]
            queries = [s * generator for s in scalars]
            if (all(not q.is_zero() for q in queries)
                    and all(not (2 * (b - a)).is_zero()
                            for a, b in combinations(queries, 2))):
                break
            rejected += 1
        completed = [fixed_secret + q for q in queries]
        values = [None if q.is_zero() else int(q[0]) for q in completed]
        public["queries"] = [encode_point(q) for q in queries]
        private.update({"query_scalars": scalars, "rejected_tuples": rejected})
    else:
        raise ValueError(family)
    public["reports"] = [None if z is None else (z // xbound) * xbound for z in values]
    private["errors"] = [None if z is None else z % xbound for z in values]
    return public, private


def supersingular_parameters(p, seed=137):
    """Cyclic E: w²=u³+u for p=3 mod 4, order p+1.

    With exactly one rational nontrivial 2-torsion point, the supersingular
    group's invariant factors (first divides gcd(p-1,p+1)) force cyclicity.
    Verify the generator's exact order independently using Sage.
    Parameter construction is outside per-input recovery timings.
    """
    if p % 4 != 3:
        raise ValueError("Require p=3 mod 4")
    set_random_seed(seed)
    ec = EllipticCurve(GF(p), [1, 0])
    ec.set_order(p + 1)
    while True:
        generator = ec.random_point()
        if generator.order() == p + 1:
            return {"curve": [1, 0], "generator": encode_point(generator), "M": p + 1}


def exhaustive_oracle(public):
    """Independent small-field enumeration over the secret domain.

    HNPHM deliberately scans multipliers, independent of floor-sum code.
    EC enumerates curve points, independent of target lift_x proposals.
    """
    p, family, reports = public["p"], public["family"], public["reports"]
    width = public["X"]
    if family == "echnp":
        ec = curve(public)
        queries = [point(ec, q) for q in public["queries"]]
        answers = []
        for secret in ec.points():
            good = True
            for q, v in zip(queries, reports):
                z = secret + q
                if v is None:
                    good = z.is_zero()
                else:
                    good = not z.is_zero() and v <= int(z[0]) < min(p, v + width)
                if not good:
                    break
            if good:
                answers.append(encode_point(secret))
        return unique(answers)
    answers = []
    for secret in range(1 if family == "hnphm" else 0, p):
        if family == "mihnp":
            good = True
            for t, v in zip(public["shifts"], reports):
                if not (secret + t) % p:
                    good = False
                    break
                z = pow((secret + t) % p, -1, p)
                if not v <= z < min(p, v + width):
                    good = False
                    break
        else:
            radius = public.get("report_model") == "radius"
            low_error, high_error = (-width, width) if radius else (0, width - 1)
            good = max(1, reports[0] + low_error) <= secret <= min(p - 1, reports[0] + high_error)
            for i in range(public["N"]):
                if not good:
                    break
                good = any(max(0, reports[public["N"] + i + 1] + low_error)
                           <= secret * t % p
                           <= min(p - 1, reports[public["N"] + i + 1] + high_error)
                           for t in range(max(1, reports[i + 1] + low_error),
                                          min(p - 1, reports[i + 1] + high_error) + 1))
        if good:
            answers.append(secret)
    return unique(answers)
