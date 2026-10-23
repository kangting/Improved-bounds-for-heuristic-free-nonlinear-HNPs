"""Exact extrema for a multilinear star polynomial on an integer box.

This optional range evaluator uses public polynomial coefficients only. It
changes neither the lattice nor the final subdivision grid. Unlike generic
termwise intervals, it combines each leaf's linear and bilinear coefficient.
"""
from __future__ import annotations


def star_interval(f, box, edges):
    """Return exact min/max of f on box when every quadratic edge meets x_0.

    f = c + l_0*x_0 + sum_{i>0} (l_i + q_i*x_0)*x_i.
    For fixed x_0, choose each leaf endpoint according to its coefficient.
    The resulting minimum is concave in x_0 and the maximum is convex, so
    both global extrema occur at an endpoint of the central interval.

    The implementation also handles reversed/duplicate star edges. Squares
    and nonstar quadratic terms are rejected, rather than silently bounded.
    Complexity is O(len(box)+len(edges)) integer operations.
    """
    k = len(box)
    if not k or len(f.linear) != k or len(f.quadratic) != len(edges):
        raise ValueError("Polynomial/box dimensions do not agree")
    if any(len(bounds) != 2 or bounds[0] > bounds[1] for bounds in box):
        raise ValueError("Require nonempty closed coordinate intervals")
    q = [0] * k
    for c, (i, j) in zip(f.quadratic, edges):
        if not c:
            continue
        if not 0 <= i < k or not 0 <= j < k:
            raise ValueError("Quadratic coordinate outside box")
        if i == j or (i != 0 and j != 0):
            raise ValueError("Require cross terms on a star centered at x_0")
        q[j if i == 0 else i] += c
    low = high = None
    for center in set(box[0]):
        current_low = current_high = f.constant + f.linear[0] * center
        for i in range(1, k):
            c = f.linear[i] + q[i] * center
            a, b = box[i]
            if c >= 0:
                current_low += c * a
                current_high += c * b
            else:
                current_low += c * b
                current_high += c * a
        low = current_low if low is None else min(low, current_low)
        high = current_high if high is None else max(high, current_high)
    return low, high
