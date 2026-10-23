"""flatter CLI preconditioning with exact recovery of its row operations.

The CLI returns only the reduced basis.  We recover and certify an integer
unimodular U with R = U B, including for full-row-rank rectangular B.  No
reduction-quality claim is made here: lattice_native finishes with LLL.
"""
from heapq import heapify, heappop, heappush
from pathlib import Path
from time import perf_counter
import os
import re
import shutil
import subprocess

from flint import fmpq_mat, fmpz_mat


def _rows(matrix):
    return [[int(matrix[i, j]) for j in range(matrix.ncols())]
            for i in range(matrix.nrows())]


def _peeling_pivots(source):
    """Find columns with one as-yet-unassigned source row.

    The HNP sources have private quadratic coordinates followed by p-pivots.
    Peeling their supports avoids a large rational inverse and also handles
    the EC star columns, which are not initially private coordinates.
    """
    m, d = len(source), len(source[0])
    supports = [{i for i in range(m) if source[i][j]} for j in range(d)]
    row_columns = [[j for j in range(d) if source[i][j]] for i in range(m)]
    ready = [j for j, support in enumerate(supports) if len(support) == 1]
    heapify(ready)
    assigned = set()
    order = []
    while ready:
        column = heappop(ready)
        if len(supports[column]) != 1:
            continue
        row = next(iter(supports[column]))
        if row in assigned:
            continue
        assigned.add(row)
        order.append((row, column))
        for j in row_columns[row]:
            supports[j].discard(row)
            if len(supports[j]) == 1:
                heappush(ready, j)
    return order if len(order) == m else None


def recover_transform(basis, reduced):
    """Return U and stats, certifying U integral, unimodular, and UB=R."""
    m, d = basis.nrows(), basis.ncols()
    if not m or m > d:
        raise ValueError('Expected a nonempty full-row-rank lattice basis')
    if (reduced.nrows(), reduced.ncols()) != (m, d):
        raise ArithmeticError('flatter changed the basis shape or rank')
    start = perf_counter()
    source, target = _rows(basis), _rows(reduced)
    order = _peeling_pivots(source)
    if order is not None:
        coefficients = [[0] * m for _ in range(m)]
        for row, column in order:
            previous = [(j, source[j][column]) for j in range(m)
                        if j != row and source[j][column]]
            pivot = source[row][column]
            for i in range(m):
                numerator = target[i][column] - sum(
                    coefficients[i][j] * entry for j, entry in previous)
                quotient, remainder = divmod(numerator, pivot)
                if remainder:
                    raise ArithmeticError('Nonintegral flatter row transformation')
                coefficients[i][row] = quotient
        transform = fmpz_mat(coefficients)
        method = 'sparse triangular pivots'
        pivots = [column for _, column in order]
    else:
        # General full-row-rank fallback.  A coordinate minor recovers row
        # operations only; flatter itself still sees the full ambient metric.
        echelon, _, rank = basis.rref()
        if rank != m:
            raise ArithmeticError('Rank-deficient source basis')
        pivots = [next(j for j in range(d) if echelon[i, j])
                  for i in range(m)]
        minor = fmpz_mat([[source[i][j] for j in pivots] for i in range(m)])
        image = fmpz_mat([[target[i][j] for j in pivots] for i in range(m)])
        # FLINT 0.9's integer=True inverse has a sign issue for det=-1.
        rational = fmpq_mat(image) * minor.inv()
        if any(rational[i, j].denominator != 1
               for i in range(m) for j in range(m)):
            raise ArithmeticError('Nonintegral flatter row transformation')
        transform = fmpz_mat([[int(rational[i, j]) for j in range(m)]
                              for i in range(m)])
        method = 'RREF pivot minor and rational inverse'
    recovery_seconds = perf_counter() - start
    start = perf_counter()
    if transform * basis != reduced:
        raise ArithmeticError('flatter full ambient identity R = U B failed')
    determinant = int(transform.det())
    if abs(determinant) != 1:
        raise ArithmeticError('flatter transformation is not unimodular')
    return transform, {
        'transform_recovery_method': method,
        'transform_recovery_seconds': recovery_seconds,
        'transform_verification_seconds': perf_counter() - start,
        'exact_transform_identity': True,
        'transform_determinant': determinant,
        'exact_unimodular_transform': True,
        'pivot_columns': pivots,
    }


def _binary():
    override = os.environ.get('HNP_FLATTER_BIN')
    if override:
        selected = shutil.which(override) or override
        path = Path(selected).expanduser().resolve()
        if not path.is_file() or not os.access(path, os.X_OK):
            raise FileNotFoundError(f'HNP_FLATTER_BIN is not executable: {path}')
        return path
    installed = shutil.which('flatter')
    if installed:
        return Path(installed).resolve()
    raise FileNotFoundError('flatter not installed; set HNP_FLATTER_BIN')


def _parse_basis(text, m, d):
    text = text.strip()
    if not (text.startswith('[') and text.endswith(']')):
        raise ValueError('Missing FPLLL outer brackets in flatter output')
    inner = text[1:-1]
    result, position = [], 0
    for match in re.finditer(r'\[([^\[\]]*)\]', inner):
        if inner[position:match.start()].strip():
            raise ValueError('Unexpected text in flatter output')
        tokens = match.group(1).split()
        if len(tokens) != d or any(not re.fullmatch(r'-?[0-9]+', s) for s in tokens):
            raise ValueError('Invalid row in flatter output')
        result.append([int(s) for s in tokens])
        position = match.end()
    if inner[position:].strip() or len(result) != m:
        raise ValueError('Unexpected shape or trailing text in flatter output')
    return fmpz_mat(result)


def precondition(basis, phase=None):
    """Execute the CLI in this worker's process group and recover exact U."""
    executable = _binary()
    args = [str(executable), '-delta', '0.99', '-v']
    environment = os.environ.copy()
    threads = {name: '1' for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS']}
    environment.update(threads)
    library_directories = [executable.parent.parent / 'lib']
    library_paths = [str(path) for path in library_directories if path.is_dir()]
    if environment.get('LD_LIBRARY_PATH'):
        library_paths.append(environment['LD_LIBRARY_PATH'])
    if library_paths:
        environment['LD_LIBRARY_PATH'] = os.pathsep.join(library_paths)
    matrix = '[\n' + '\n'.join('[' + ' '.join(map(str, row)) + ']'
                               for row in _rows(basis)) + '\n]\n'
    if phase:
        phase('flatter_start', {'flatter_binary': str(executable),
                               'flatter_args': args, 'thread_environment': threads})
    start = perf_counter()
    process = subprocess.run(args, input=matrix, capture_output=True, text=True,
                             env=environment)
    flatter_seconds = perf_counter() - start
    if phase:
        phase('flatter_done', {'flatter_seconds': flatter_seconds,
                              'returncode': process.returncode})
    if process.returncode:
        raise RuntimeError(f'flatter exit {process.returncode}: {process.stderr[-8192:]}')
    reduced = _parse_basis(process.stdout, basis.nrows(), basis.ncols())
    transform, stats = recover_transform(basis, reduced)
    if phase:
        phase('flatter_transform_done', stats)
    stats.update({
        'flatter_binary': str(executable),
        'flatter_args': args,
        'flatter_thread_environment': threads, 'flatter_seconds': flatter_seconds,
        'flatter_ld_library_path': environment.get('LD_LIBRARY_PATH'),
        'flatter_stderr': process.stderr[-65536:],
        'flatter_requested_quality_is_approximate': True,
        'input_rank': basis.nrows(), 'ambient_dimension': basis.ncols(),
    })
    return reduced, transform, stats
