"""Exact standard LLL fallback with integer row-operation tracking.

This intentionally simple implementation is reserved for a failure of the
fast lattice backend. It accepts a full-row-rank rectangular integer basis,
uses rational delta=3/4 and exact size reduction to eta=1/2, and returns
(reduced, U) with reduced=U*basis and U unimodular.

Gram--Schmidt data are recomputed after each nontrivial row operation. This
is slower than optimized LLL but avoids floating arithmetic and external LLL
backend behavior. Standard LLL termination and polynomial bit bounds apply;
the implementation places no iteration or coefficient cutoff on recovery.
"""
from fractions import Fraction
from flint import fmpz_mat


_HALF = Fraction(1, 2)
_DELTA = Fraction(3, 4)


def _gram_schmidt(rows):
    """Return exact projection coefficients and squared orthogonal norms."""
    orthogonal, norms, coefficients = [], [], []
    for row in rows:
        current = list(map(Fraction, row))
        mu = []
        for previous, norm in zip(orthogonal, norms):
            coefficient = sum(a*b for a,b in zip(row,previous)) / norm
            mu.append(coefficient)
            if coefficient:
                current = [a-coefficient*b for a,b in zip(current,previous)]
        norm = sum(value*value for value in current)
        if norm <= 0:
            raise ValueError('Exact LLL requires a full-row-rank basis')
        orthogonal.append(current)
        norms.append(norm)
        coefficients.append(mu)
    return coefficients, norms


def exact_lll(basis):
    """Reduce a flint integer row basis exactly and track a unimodular U.

    Half-integer ties are left in place since their magnitudes already equal
    eta=1/2. For a reduction, nearest-integer rounding uses floor(mu+1/2).
    """
    n, dimension = basis.nrows(), basis.ncols()
    if n == 0 or n > dimension:
        raise ValueError('Require a nonempty full-row-rank integer basis')
    rows = [[int(basis[i,j]) for j in range(dimension)] for i in range(n)]
    transform = [[int(i==j) for j in range(n)] for i in range(n)]
    mu, norms = _gram_schmidt(rows)  # Also reject dependent input before swaps.
    k = 1
    while k < n:
        for j in reversed(range(k)):
            if abs(mu[k][j]) > _HALF:
                nearest_fraction = mu[k][j] + _HALF
                nearest = nearest_fraction.numerator // nearest_fraction.denominator
                rows[k] = [a-nearest*b for a,b in zip(rows[k],rows[j])]
                transform[k] = [a-nearest*b for a,b in zip(transform[k],transform[j])]
                mu, norms = _gram_schmidt(rows)
        if norms[k] >= (_DELTA-mu[k][k-1]**2)*norms[k-1]:
            k += 1
        else:
            rows[k],rows[k-1] = rows[k-1],rows[k]
            transform[k],transform[k-1] = transform[k-1],transform[k]
            mu, norms = _gram_schmidt(rows)
            k = max(1,k-1)
    reduced = fmpz_mat(rows)
    unimodular = fmpz_mat(transform)
    # Keep the exact identity check in the fallback, too.
    if unimodular*basis != reduced:
        raise ArithmeticError('Exact LLL row tracking failed')
    return reduced,unimodular
