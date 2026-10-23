"""Exact Babai nearest plane using fraction-free Gram elimination.

This is an arithmetic optimization of baselines_native.babai. It uses the
same LLL basis, rational projections and half-up integer rounding. It does
not change the lattice, nearest-plane order or approximation factor.

For G=B B^T positive definite, fraction-free LU has no pivot permutation.
If Delta_j=U[j,j] and Delta_{-1}=1, Gram--Schmidt gives
  mu[i,j]=U[j,i]/Delta_j,  ||b_j^*||^2=Delta_j/Delta_{j-1}.
Thus the projections of a target and the reverse nearest-plane updates can
be computed in O(n^2) rational operations, without constructing n ambient
rational Gram--Schmidt vectors.
"""
from time import perf_counter
from flint import fmpz_mat, fmpq
from lattice_native import tracked_lll


def nearest_plane_reduced(basis, target, *, return_info=False):
    """Nearest plane in an already reduced, full-row-rank integer basis.

    The arithmetic remains exact for rectangular bases. ``return_info`` adds
    timing and structural diagnostics; it never receives evaluation truth.
    """
    n, d = basis.nrows(), basis.ncols()
    if len(target) != d:
        raise ValueError('Target dimension does not match the basis')
    if n == 0 or n > d:
        raise ValueError('Require a nonempty, full-row-rank lattice basis')
    start = perf_counter()
    gram = basis * basis.transpose()
    gram_time = perf_counter() - start
    start = perf_counter()
    permutation, lower, diagonal, upper = gram.fflu()
    fflu_time = perf_counter() - start
    if any(int(permutation[i,j]) != int(i==j) for i in range(n) for j in range(n)):
        raise ArithmeticError('Positive-definite Gram matrix required an unexpected pivot')
    pivots = [upper[i,i] for i in range(n)]
    if any(pivot <= 0 for pivot in pivots):
        raise ValueError('Basis Gram matrix is not positive definite')
    start = perf_counter()
    mu = [[fmpq(upper[j,i],pivots[j]) for j in range(i)] for i in range(n)]
    dots = fmpz_mat(1,d,list(map(int,target))) * basis.transpose()
    unnormalized = []
    coefficients = []
    previous_pivot = 1
    for i in range(n):
        projection = fmpq(dots[0,i])
        for j in range(i):
            projection -= mu[i][j] * unnormalized[j]
        unnormalized.append(projection)
        coefficients.append(projection * fmpq(previous_pivot,pivots[i]))
        previous_pivot = pivots[i]
    projection_time = perf_counter() - start
    start = perf_counter()
    integers = [0]*n
    for i in reversed(range(n)):
        nearest = int((coefficients[i]+fmpq(1,2)).floor())
        integers[i] = nearest
        if nearest:
            for j in range(i):
                coefficients[j] -= nearest * mu[i][j]
    closest = fmpz_mat(1,n,integers) * basis
    result = [int(closest[0,j]) for j in range(d)]
    rounding_time = perf_counter() - start
    if return_info:
        return result, {'gram_seconds':gram_time,'fflu_seconds':fflu_time,
                        'projection_seconds':projection_time,
                        'rounding_seconds':rounding_time,
                        'rank':n,'ambient_dimension':d,
                        'max_gram_minor_bits':max(int(p).bit_length() for p in pivots)}
    return result


def babai(basis,target):
    """Drop-in replacement for baselines_native.babai, including shared LLL."""
    reduced,_=tracked_lll(basis)
    return nearest_plane_reduced(reduced,target)
