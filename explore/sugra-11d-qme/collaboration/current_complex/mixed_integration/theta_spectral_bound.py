"""A rational operator-norm bound on the four real, positive-resistance Wick flows."""

from fractions import Fraction
from math import isqrt


def root_upper(square):
    scale = 4096
    scaled = square * scale**2
    lo = isqrt(scaled)
    return Fraction(lo + (lo * lo != scaled), scale)


def pairing(u, v):
    return sum(x * y for x, y in zip(u, v))


def wick_bound(momenta, polarizations):
    assert all(sum(axis) == 0 for axis in zip(*momenta))
    subsets = [(0,) * 5]
    for wave in momenta:
        subsets.extend([tuple(x + y for x, y in zip(v, wave)) for v in subsets])
    radius = root_upper(max(pairing(v, v) for v in subsets))
    total = Fraction(0)
    for weight, waves in polarizations:
        transfers = Fraction(1)
        for u, p in zip(waves[:9], momenta[:9]):
            assert pairing(u, p) == 0
            trace_twice = sum(abs(pairing(u, q)) for q in momenta)
            rank_one = pairing(p, p) * pairing(u, u)
            transfers *= (root_upper(trace_twice**2 + rank_one) + root_upper(rank_one)) / 2
        eta, q = waves[9], momenta[9]
        wedge = root_upper(pairing(eta, eta) * pairing(q, q) - pairing(eta, q)**2)
        total += abs(weight) * transfers * wedge
    return 2 * (8 + root_upper(48)) * radius**3 * total
