"""Official P-384 domain, NIST SP800-186 Section3.2.1.4 (pp11–12)."""
P = 2**384 - 2**128 - 2**96 + 2**32 - 1
A = -3
B = int("b3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef", 16)
GX = int("aa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7", 16)
GY = int("3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f", 16)
M = int("ffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc52973", 16)
SOURCE_URL = "https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-186.pdf"
PARAMETERS = {"curve": [A, B], "generator": [GX, GY], "M": M}


def validate(prove=True):
    from sage.all import ZZ
    from solver import core
    assert P.bit_length() == M.bit_length() == 384
    if prove:
        assert ZZ(P).is_prime(proof=True) and ZZ(M).is_prime(proof=True)
    ec = core.curve({"p": P, **PARAMETERS})
    generator = core.point(ec, PARAMETERS["generator"])
    assert not generator.is_zero() and (M * generator).is_zero()
    assert (4 * A**3 + 27 * B**2) % P != 0
    # M divides #E. Hasse's interval contains only this multiple of prime M.
    radius = int(ZZ(4 * P).isqrt())
    assert P + 1 - radius <= M <= P + 1 + radius
    assert 2 * M > P + 1 + radius
    ec.set_order(M)
    generator.set_order(M)
    assert generator.order() == M
    return {"p": P, **PARAMETERS, "cofactor": 1, "p_prime_proven": prove,
            "M_prime_proven": prove, "generator_on_curve": True,
            "generator_exact_order": True, "curve_order_by_hasse": M,
            "ec_polylog_subgroup_condition_c": 1, "source_url": SOURCE_URL}
