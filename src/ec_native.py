"""Small exact prime-field EC backend for HNP experiments.

Affine arithmetic is deliberately simple and uses Python arbitrary precision
integers.  This is an experiment backend, not a constant-time crypto library.
P-384 constants: NIST SP 800-186, Section 3.2.1.4 (February 2023).
https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-186.pdf
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache


def sqrt_mod(value: int, p: int) -> int | None:
    """One square root modulo an odd prime, or None for a nonresidue."""
    value %= p
    if value == 0:
        return 0
    if pow(value, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        root = pow(value, (p + 1) // 4, p)
        return min(root, p - root)
    # Tonelli--Shanks, including primes p = 1 mod 4 used in small tests.
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    from secrets import randbelow
    z = 2 + randbelow(p-2)
    while pow(z, (p - 1) // 2, p) != p - 1:
        z = 2 + randbelow(p-2)
    c = pow(z, q, p)
    root = pow(value, (q + 1) // 2, p)
    t, m = pow(value, q, p), s
    while t != 1:
        i, power = 0, t
        while power != 1 and i < m:
            power = power * power % p
            i += 1
        if i == m:
            raise ArithmeticError("Tonelli--Shanks invariant failed")
        b = pow(c, 1 << (m - i - 1), p)
        root = root * b % p
        c = b * b % p
        t = t * c % p
        m = i
    return min(root, p - root)


@dataclass(frozen=True)
class EC:
    p: int
    a: int
    b: int

    def __post_init__(self):
        p = int(self.p)
        if p <= 3:
            raise ValueError("Require an odd prime p > 3")
        object.__setattr__(self, "p", p)
        object.__setattr__(self, "a", int(self.a) % p)
        object.__setattr__(self, "b", int(self.b) % p)
        if (4 * self.a ** 3 + 27 * self.b ** 2) % p == 0:
            raise ValueError("Singular elliptic curve")

    def __call__(self, coordinates=0):
        if coordinates is None or coordinates == 0:
            return Point(self, None, None)
        if isinstance(coordinates, Point):
            if coordinates.ec != self:
                raise ValueError("Different curves")
            return coordinates
        return Point(self, int(coordinates[0]), int(coordinates[1]))

    def lift_x(self, x: int, all: bool = False):
        x = int(x) % self.p
        root = sqrt_mod((x ** 3 + self.a * x + self.b) % self.p, self.p)
        roots = [] if root is None else [root] if root == 0 else [root, self.p - root]
        points = [Point(self, x, y) for y in roots]
        if all:
            return points
        if not points:
            raise ValueError("No curve point with this x-coordinate")
        return points[0]

    def points(self):
        """Independent tiny-field oracle only; never used by large solvers."""
        if self.p > 100_000:
            raise ValueError("Full point enumeration is restricted to tiny tests")
        yield self(0)
        for x in range(self.p):
            yield from self.lift_x(x, all=True)


@dataclass(frozen=True)
class Point:
    ec: EC
    x: int | None
    y: int | None

    def __post_init__(self):
        if self.x is None or self.y is None:
            if self.x is not None or self.y is not None:
                raise ValueError("Infinity needs both coordinates None")
            return
        p = self.ec.p
        x, y = int(self.x) % p, int(self.y) % p
        if (y * y - x ** 3 - self.ec.a * x - self.ec.b) % p:
            raise ValueError("Point is not on the curve")
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "y", y)

    def is_zero(self):
        return self.x is None

    def __getitem__(self, index):
        if self.is_zero():
            raise ValueError("Infinity has no affine coordinates")
        return (self.x, self.y)[index]

    def xy(self):
        if self.is_zero():
            raise ValueError("Infinity has no affine coordinates")
        return self.x, self.y

    def __neg__(self):
        return self if self.is_zero() else Point(self.ec, self.x, -self.y)

    def __add__(self, other):
        if not isinstance(other, Point):
            return NotImplemented
        if self.ec != other.ec:
            raise ValueError("Points belong to different curves")
        if self.is_zero():
            return other
        if other.is_zero():
            return self
        p = self.ec.p
        if self.x == other.x:
            if (self.y + other.y) % p == 0:
                return self.ec(0)
            slope = (3 * self.x ** 2 + self.ec.a) * pow(2 * self.y, -1, p) % p
        else:
            slope = (other.y - self.y) * pow(other.x - self.x, -1, p) % p
        x = (slope * slope - self.x - other.x) % p
        y = (slope * (self.x - x) - self.y) % p
        return Point(self.ec, x, y)

    def __sub__(self, other):
        return self + (-other)

    def __rmul__(self, scalar):
        scalar = int(scalar)
        if scalar < 0:
            return (-scalar) * (-self)
        answer, term = self.ec(0), self
        while scalar:
            if scalar & 1:
                answer = answer + term
            scalar >>= 1
            if scalar:
                term = term + term
        return answer

    __mul__ = __rmul__


@lru_cache(maxsize=32)
def _curve(p, a, b):
    return EC(int(p), int(a), int(b))


def curve(public):
    return _curve(int(public["p"]), *map(int, public["curve"]))


def point(ec, encoded):
    return ec(encoded)


def encode_point(q):
    return None if q.is_zero() else [q.x, q.y]


P384_P = 2 ** 384 - 2 ** 128 - 2 ** 96 + 2 ** 32 - 1
P384_PARAMS = {
    "curve": [-3, int("B3312FA7E23EE7E4988E056BE3F82D19181D9C6EFE814112"
                      "0314088F5013875AC656398D8A2ED19D2A85C8EDD3EC2AEF", 16)],
    "generator": [
        int("AA87CA22BE8B05378EB1C71EF320AD746E1D3B628BA79B98"
            "59F741E082542A385502F25DBF55296C3A545E3872760AB7", 16),
        int("3617DE4A96262C6F5D9E98BF9292DC29F8F41DBD289A147C"
            "E9DA3113B5F0B8C00A60B1CE1D7E819D7A431D7C90EA0E5F", 16)],
    "M": int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
             "C7634D81F4372DDF581A0DB248B0A77AECEC196ACCC52973", 16),
    "curve_name": "NIST P-384",
    "cofactor": 1,
    "parameter_source": "NIST SP 800-186, Section 3.2.1.4",
}


@lru_cache(maxsize=1)
def validate_p384():
    """Check prime p,M, nonsingularity, membership, and exact generator order."""
    from flint import fmpz
    ec = curve({"p": P384_P, **P384_PARAMS})
    g = point(ec, P384_PARAMS["generator"])
    m = P384_PARAMS["M"]
    checks = {
        "p_prime": bool(fmpz(P384_P).is_prime()),
        "M_prime": bool(fmpz(m).is_prime()),
        "generator_nonzero": not g.is_zero(),
        "generator_on_curve": (g.y * g.y - g.x ** 3 - ec.a * g.x - ec.b) % ec.p == 0,
        "M_times_generator_infinity": (m * g).is_zero(),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Invalid P-384 parameters: {checks}")
    # Since M is prime and G != O, [M]G=O proves that G has exact order M.
    return {**checks, "generator_exact_order": m, "source": P384_PARAMS["parameter_source"]}


def small_parameters(p=127, seed=None):
    """Deterministically find a full-order generator for y^2=x^3+x, p=127.

    Enumeration and factorization are confined to this tiny test helper.
    ``seed`` is accepted to ease migration from the old Sage helper.
    """
    if p != 127:
        raise ValueError("This fixture is defined only for p=127")
    ec = _curve(p, 1, 0)
    points = list(ec.points())
    if len(points) != p + 1:
        raise ArithmeticError("Unexpected test curve order")
    # p+1=128 is a power of two: [128]G=O and [64]G!=O prove order 128.
    g = next(q for q in points if not q.is_zero() and not (64 * q).is_zero())
    if not (128 * g).is_zero():
        raise ArithmeticError("Invalid tiny generator")
    return {"curve": [1, 0], "generator": encode_point(g), "M": 128}


# Former artifact helper name; only the tiny fixture is supported here.
supersingular_parameters = small_parameters
