"""flatter preconditioning followed by tracked LLL and exact identities.

flatter's approximate quality parameter is not used as an LLL guarantee.
The full ambient basis is finished by the existing FLINT LLL routine.
"""
from fractions import Fraction
from math import gcd
from time import perf_counter
from fpylll import IntegerMatrix, LLL
from flint import fmpz_mat


def rows(a):
    return [[int(a[i,j]) for j in range(a.ncols())] for i in range(a.nrows())]


LAST_LLL = {}

def _finishing_lll(basis):
    # FLINT eta=.5 avoids fpLLL's Babai-boundary loop. Some FLINT 0.9
    # exact-Gram transformations fail their identity; reject them and retry.
    failures=[]
    for gram,delta in [('exact',.99),('approx',.99),('exact',.75)]:
        try:
            reduced,transform=basis.lll(transform=True,delta=delta,eta=.5,gram=gram)
        except Exception as exc:
            failures.append({'gram':gram,'delta':delta,'reason':repr(exc)})
            continue
        if transform*basis==reduced:
            return reduced,transform,{'backend':'FLINT','gram':gram,
                         'delta':delta,'eta':.5,'rejected_transformations':failures}
        failures.append({'gram':gram,'delta':delta,'reason':'B != U*B0'})
    from exact_lll_fallback import exact_lll
    reduced,transform=exact_lll(basis)
    return reduced,transform,{'backend':'Fraction exact fallback',
                             'delta':.75,'eta':.5,'rejected_transformations':failures}


def tracked_lll(basis, phase=None):
    from flatter_backend import precondition
    LAST_LLL.clear()
    total_start = perf_counter()
    start = perf_counter()
    source = rows(basis)
    scalar = 0
    largest = 0
    for row in source:
        for entry in row:
            if entry:
                scalar = gcd(scalar, abs(entry))
                largest = max(largest, abs(entry).bit_length())
    if not scalar:
        raise ArithmeticError('Zero source basis')
    # A uniform positive scalar preserves both the lattice geometry up to
    # scale and every row operation.  Unlike selecting a coordinate minor,
    # this does not alter the ambient Euclidean metric relative to itself.
    normalized = (fmpz_mat([[entry // scalar for entry in row] for row in source])
                  if scalar != 1 else basis)
    normalization_seconds = perf_counter() - start
    initial, first, stats = precondition(normalized, phase=phase)
    if phase:
        phase('finishing_lll_start', {'flatter_seconds':stats['flatter_seconds'],
                                     'rank':initial.nrows(),
                                     'ambient_dimension':initial.ncols()})
    start = perf_counter()
    normalized_reduced, final, finishing = _finishing_lll(initial)
    finishing_seconds = perf_counter() - start
    if phase:
        phase('finishing_lll_done', {'finishing_lll_seconds':finishing_seconds,
                                    'finishing_lll':finishing})
    start = perf_counter()
    if abs(int(final.det())) != 1:
        raise ArithmeticError('Finishing LLL transform is not unimodular')
    transform = final * first
    reduced = normalized_reduced * scalar
    if transform * basis != reduced:
        raise ArithmeticError('Composed flatter/LLL identity failed')
    LAST_LLL.update({**stats, 'backend':'flatter + ' + finishing['backend'],
                    'delta':finishing['delta'], 'eta':finishing['eta'],
                    'finishing_lll':finishing,
                    'uniform_scale_factor':scalar,
                    'uniform_normalization_seconds':normalization_seconds,
                    'input_basis_max_bits':largest,
                    'normalized_input_max_bits':max(abs(int(normalized[i,j])).bit_length()
                        for i in range(normalized.nrows()) for j in range(normalized.ncols())),
                    'flatter_transform_identity_domain':'uniformly normalized full ambient basis',
                    'exact_restored_basis_identity':True,
                    'finishing_lll_seconds':finishing_seconds,
                    'composed_transform_verification_seconds':perf_counter()-start,
                    'exact_composed_transform_identity':True,
                    'exact_finishing_unimodular_transform':True,
                    'total_reduction_seconds':perf_counter()-total_start})
    if phase:
        phase('transform_done', {'total_reduction_seconds':LAST_LLL['total_reduction_seconds'],
                                 'exact_composed_transform_identity':True,
                                 'exact_unimodular_transform':True})
    return reduced,transform


def exact_lll_audit(reduced):
    """Independent exact (3/4,1/2)-LLL audit, intended for small cases."""
    r = rows(reduced)
    gs, norms = [], []
    for i, b in enumerate(r):
        g = list(map(Fraction,b))
        for j, previous in enumerate(gs):
            mu = sum(x*y for x,y in zip(b,previous)) / norms[j]
            if abs(mu)>Fraction(1,2):
                raise ArithmeticError('Exact size reduction failed')
            if j==i-1:
                last_mu=mu
            g = [x-mu*y for x,y in zip(g,previous)]
        norm=sum(x*x for x in g)
        if norm<=0:
            raise ArithmeticError('Rank deficient')
        if i and norm < (Fraction(3,4)-last_mu**2)*norms[-1]:
            raise ArithmeticError('Exact Lovasz condition failed')
        gs.append(g);norms.append(norm)
    return {'independent_exact_lll': True, 'size_reduced_half': True,
            'lovasz_three_quarters': True}
