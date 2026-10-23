"""Native ports of the attached prior LLL/Babai instantiations.

Preserves the lattice formulas in hnp2/src/baselines.py. These are polynomial
LLL/Babai implementations, not exact SVP/CVP or a reproduction of every
constant in the original recovery guarantees. Only public input is accepted.

Primary descriptions:
LSSW12, Section 2.1, dimension 2N:
https://users.monash.edu.au/~rste/HNP-Invers-JSC.pdf
Shani17, Sections 4.1.2/5.1, dimension 3n+3 (paired observations):
https://eprint.iacr.org/2016/1189
HGNS03, Lemma 2 and Theorem 1:
https://doi.org/10.1090/S0025-5718-03-01495-9
"""
from fractions import Fraction
from math import exp, log, sqrt
from time import perf_counter

from flint import fmpz_mat, fmpq

from ec_native import curve, encode_point, point
from hnp_native import decode, unique
from lattice_native import tracked_lll as original_tracked_lll, LAST_LLL as SHARED_LLL


LAST_PRIOR_LLL = {}

def tracked_lll_std75(basis):
    """Standard (3/4,1/2)-LLL for all published prior instantiations."""
    try:
        reduced,transform=basis.lll(transform=True,delta=.75,eta=.5,gram='approx')
        if transform*basis != reduced:raise ArithmeticError('Prior LLL transformation identity failed')
    except Exception:
        from exact_lll_fallback import exact_lll
        reduced,transform=exact_lll(basis)
    LAST_PRIOR_LLL.clear();LAST_PRIOR_LLL.update({'backend':'FLINT','gram':'approx',
                                  'delta':.75,'eta':.5,'exact_transform_identity':True})
    return reduced,transform


def tracked_lll(basis):
    """Same flatter + standard LLL backend as the compared quotient search."""
    reduced,transform=original_tracked_lll(basis)
    LAST_PRIOR_LLL.clear()
    LAST_PRIOR_LLL.update(dict(SHARED_LLL))
    LAST_PRIOR_LLL['comparison_backend']='shared_flatter_finishing_lll'
    return reduced,transform


def midpoints(public):
    p, width = public["p"], public["X"]
    centers = [(v + min(p - 1, v + width - 1)) // 2 for v in public["reports"]]
    # Same Delta across observations, including a final clipped bucket.
    return centers, 1 + width // 2


def shortest_approximation(basis):
    """First vector of the shared standard LLL implementation."""
    reduced, _ = tracked_lll(basis)
    return [int(reduced[0, j]) for j in range(reduced.ncols())], reduced


def lssw12(public):
    start = perf_counter()
    if public["family"] != "mihnp":
        raise ValueError("LSSW is a MIHNP algorithm")
    p, n = public["p"], public["N"]
    shifts = public["shifts"]
    centers, delta = midpoints(public)
    d = 2 * n
    basis = fmpz_mat(d, d)
    basis[0, 0] = delta ** 2
    for i in range(n):
        basis[1 + i, 1 + i] = delta
    for j in range(1, n):
        col = n + j
        c = (shifts[0] - shifts[j]) % p
        inv_c = pow(c, -1, p)
        aa = c * centers[0] * centers[j] + centers[0] - centers[j]
        b1, bj = -c * centers[j] - 1, -c * centers[0] + 1
        basis[0, col] = -inv_c * aa % p
        basis[1, col] = -inv_c * b1 % p
        basis[1 + j, col] = -inv_c * bj % p
        basis[col, col] = p
    construction_seconds = perf_counter() - start
    start = perf_counter()
    shortest, _ = shortest_approximation(basis)
    lll_seconds = perf_counter() - start
    start = perf_counter()
    answers = []
    if shortest[0]:
        for i in range(n):
            correction = Fraction(delta * shortest[1 + i], shortest[0])
            if correction.denominator != 1:
                continue
            completed = (centers[i] - int(correction)) % p
            if not completed:
                continue
            secret = (pow(completed, -1, p) - shifts[i]) % p
            target_completed = (secret + shifts[0]) % p
            if target_completed:
                answers.extend(decode(public, pow(target_completed, -1, p)
                                       - public["reports"][0]))
                if answers:
                    break
    return {"status": "RETURNED" if answers else "NO_CANDIDATE",
            "secrets": unique(answers), "dimension": d,
            "Delta": delta, "svp_gamma_log2": (d - 1) / 2,
            "timing": {"source": construction_seconds, "lll": lll_seconds,
                       "readout_verify": perf_counter() - start}}


def exact_gram_schmidt_rows(basis, audit_lll=False):
    """Exact rational row Gram--Schmidt; intended for the prior's small CVP basis."""
    rows = [[fmpq(basis[i, j]) for j in range(basis.ncols())]
            for i in range(basis.nrows())]
    gs, norms = [], []
    for index,raw in enumerate(rows):
        row = list(raw)
        last_mu=None
        for j,(previous, norm) in enumerate(zip(gs, norms)):
            coefficient = sum(a * b for a, b in zip(raw, previous)) / norm
            if audit_lll and abs(coefficient)>fmpq(1,2):
                raise ArithmeticError('Prior exact LLL size condition failed')
            if j==index-1:last_mu=coefficient
            if coefficient:
                row = [a - coefficient * b for a, b in zip(row, previous)]
        norm = sum(c * c for c in row)
        if not norm:
            raise ValueError("Rank-deficient nearest-plane basis")
        if audit_lll and index and norm < (fmpq(3,4)-last_mu**2)*norms[-1]:
            raise ArithmeticError('Prior exact LLL Lovasz condition failed')
        gs.append(row)
        norms.append(norm)
    return rows, gs, norms


def babai(basis, target, audit_lll=True):
    """Exact rational nearest plane after LLL, same half-up tie rule as original."""
    reduced, _ = tracked_lll(basis)
    rows, gs, norms = exact_gram_schmidt_rows(reduced, audit_lll=audit_lll)
    residual = [fmpq(c) for c in target]
    for i in reversed(range(reduced.nrows())):
        coefficient = sum(a * b for a, b in zip(residual, gs[i])) / norms[i]
        nearest = int((coefficient + fmpq(1, 2)).floor())
        if nearest:
            residual = [a - nearest * b for a, b in zip(residual, rows[i])]
    result = [fmpq(a) - b for a, b in zip(target, residual)]
    if any(x.q != 1 for x in result):
        raise ArithmeticError("Nearest-plane result not integral")
    return [int(x) for x in result]


def small_multiplier(a, p, delta):
    """Last rational-CF denominator strictly below sqrt(p/Delta), HGNS03 Lemma 2."""
    modulus = p
    previous = 1
    qminus2, qminus1 = 1, 0
    while p:
        term, remainder = divmod(a, p)
        denominator = term * qminus1 + qminus2
        if denominator * denominator * delta >= modulus:
            break
        previous = denominator
        a, p = p, remainder
        qminus2, qminus1 = qminus1, denominator
    return previous


def hgns03(public):
    start = perf_counter()
    if public["family"] != "hnphm":
        raise ValueError("HGNS is a HNPHM algorithm")
    p, n = public["p"], public["N"]
    if public.get("report_model") == "radius":
        centers, delta = public["reports"], public["X"]
    else:
        centers, delta = midpoints(public)
    multiplier = small_multiplier(centers[0], p, delta)
    v = [multiplier * t % p for t in centers[1:n + 1]]
    u = [multiplier * z % p for z in centers[n + 1:]]
    basis = fmpz_mat(n + 1, n + 1)
    for i in range(n):
        basis[i, i] = p ** 2
        basis[n, i] = p * v[i]
    basis[n, n] = 1
    construction_seconds = perf_counter() - start
    start = perf_counter()
    closest = babai(basis, [p * z for z in u] + [0])
    cvp_seconds = perf_counter() - start
    start = perf_counter()
    secret = int(closest[-1]) % p
    answers = decode(public, secret - public["reports"][0])
    return {"status": "RETURNED" if answers else "NO_CANDIDATE",
            "secrets": answers, "dimension": n + 1, "Delta": delta,
            "lambda": multiplier, "cvp_factor_log2": (n + 1) / 2,
            "timing": {"source": construction_seconds, "lll_babai": cvp_seconds,
                       "readout_verify": perf_counter() - start}}


def _infinity_echnp(public):
    ec = curve(public)
    idx = public["reports"].index(None)
    secret = -point(ec, public["queries"][idx])
    for query, report in zip(public["queries"], public["reports"]):
        shifted = secret + point(ec, query)
        if report is None:
            if not shifted.is_zero():
                return []
        elif shifted.is_zero() or not report <= int(shifted[0]) < min(public["p"], report + public["X"]):
            return []
    return [encode_point(secret)]


def shani17(public):
    start = perf_counter()
    if public.get("observation_model") != "paired":
        return {"status": "NOT_APPLICABLE_MODEL", "secrets": None, "timing": {}}
    if None in public["reports"]:
        return {"status": "RETURNED", "secrets": _infinity_echnp(public),
                "exception": "infinity_determines_secret", "dimension": 0,
                "Delta": 1, "svp_gamma_log2": 0,
                "timing": {"exception_decode": perf_counter() - start}}
    p, count = public["p"], public["N"]
    if count % 2 != 1:
        raise ValueError("Paired EC input needs an odd observation count")
    n = (count - 1) // 2
    centers, delta = midpoints(public)
    h0 = centers[0]
    ec = curve(public)
    base = point(ec, public["queries"][0])
    queries = [point(ec, q) for q in public["queries"]]
    d = 3 * n + 3
    basis = fmpz_mat(d, d)
    basis[0, 0] = delta ** 3
    for i in range(n + 1):
        basis[1 + i, 1 + i] = delta ** 2
        basis[n + 2 + i, n + 2 + i] = delta
    aa_curve, b_curve = public["curve"]
    for i in range(n):
        plus, minus = queries[1 + 2 * i], queries[2 + 2 * i]
        if plus + minus != 2 * base:
            raise ValueError("Not a paired EC query tuple")
        q = int((plus - base)[0])
        h = centers[1 + 2 * i] + centers[2 + 2 * i]
        aa = h - 2 * q
        ap = 2 * (h0 - q)
        bb = 2 * (h * (h0 - q) - 2 * h0 * q - aa_curve - q ** 2)
        bp = (h0 - q) ** 2
        cc = h * (h0 - q) ** 2 - 2 * ((h0 ** 2 + aa_curve) * q
                                      + (aa_curve + q ** 2) * h0 + 2 * b_curve)
        col = 2 * n + 3 + i
        basis[0, col] = -cc % p
        basis[1, col] = -bb % p
        basis[2 + i, col] = -bp % p
        basis[n + 2, col] = -aa % p
        basis[n + 3 + i, col] = -ap % p
        basis[col, col] = p
    construction_seconds = perf_counter() - start
    start = perf_counter()
    shortest, _ = shortest_approximation(basis)
    lll_seconds = perf_counter() - start
    start = perf_counter()
    answers = []
    if shortest[0]:
        correction = Fraction(delta * shortest[1], shortest[0])
        if correction.denominator == 1:
            completed = (h0 + int(correction)) % p
            answers = decode(public, completed - public["reports"][0])
    return {"status": "RETURNED" if answers else "NO_CANDIDATE",
            "secrets": unique(answers), "dimension": d, "Delta": delta,
            "svp_gamma_log2": (d - 1) / 2,
            "timing": {"source": construction_seconds, "lll": lll_seconds,
                       "readout_verify": perf_counter() - start}}


def prior(public):
    return {"mihnp": lssw12, "echnp": shani17, "hnphm": hgns03}[public["family"]](public)

def prior_finite_bound(public, result):
    """Original conservative probability upper bound, logarithmic to avoid overflow.

    A positive log bound is vacuous, not a guarantee of failure. These include
    the LLL approximation constants, rather than pretending to use exact SVP.
    """
    if result["status"] == "NOT_APPLICABLE_MODEL":
        return None
    p, n, delta = public["p"], public["N"], result["Delta"]
    family = public["family"]
    lp = log(p)
    if family == "mihnp":
        h = 2 ** result["svp_gamma_log2"] * sqrt(2 * n)
        a = n * log(2) + (3 * n - 2) * log(4 * h * delta + 1) - (n - 1) * log(p - 1)
        b = log(4) + 3 * log(4 * h * delta + 1) - log(p - 1)
        terms = [a, b]
    elif family == "echnp":
        n = (n - 1) // 2
        eta = 2 * 2 ** result["svp_gamma_log2"] * sqrt(3 * n + 1)
        denominator = p - 2 * sqrt(p) - 2
        if denominator <= 0:
            return {"log_failure_upper_bound": None, "certifies_nontrivial_probability": False}
        ld = log(denominator)
        scale = log(6 * eta * delta + 1)
        terms = [n * log(8) + (6 * n + 3) * scale - n * ld,
                 log(16) + 6 * scale - ld,
                 log(2 * n + 3) - log(p - 2 * sqrt(p)),
                 n * log(16 * delta) - n * ld]
    else:
        # Conservative adaptation to the explicitly advertised dimension-N+1
        # Babai factor, rather than the paper's sharper CVP constants.
        terms = [lp + n * (log(33) + log(n) / 2 + (n + 1) * log(2) / 2
                          + lp / 2 + 2.5 * log(delta) - log(p - 1))]
    m = max(terms)
    from math import exp
    bound = m + log(sum(exp(t - m) for t in terms))
    return {"log_failure_upper_bound": bound,
            "failure_upper_bound_clipped": exp(min(0.0, bound)),
            "reference_bound_nontrivial": bound < 0,
            "reference_only": True,
            "model_note": "Primary unconditioned probability expression; not a certificate for the accepted sampling law."}
