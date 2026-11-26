"""Universal component denominator from resistance and conductance degrees."""

from fractions import Fraction
from functools import cache
from itertools import permutations, combinations
from math import comb, factorial, gcd, lcm, prod


SHAPES = (
    ((1, 6, 6), 158400, 4),
    ((2, 5, 6), 285120, 2),
    ((3, 4, 6), 332640, 2),
    ((3, 5, 5), 798336, 4),
    ((4, 4, 5), 997920, 4),
)


def compositions(total, length):
    if length == 0:
        if total == 0:
            yield ()
    elif length == 1:
        yield (total,)
    else:
        for i in range(total + 1):
            for tail in compositions(total - i, length - 1):
                yield (i, *tail)


@cache
def gap_factor(count, degree):
    # The initial prefix vanishes. Each other prefix is attached to one gap.
    numer = gcd(
        *(prod(factorial(t) for t in parts) for parts in compositions(degree, count))
    )
    return Fraction(factorial(count) * numer, factorial(count + degree))


@cache
def moment(shape, b, c):
    return Fraction(
        prod(13 - 2 * shape[1] + 2 * j for j in range(b))
        * prod(13 - 2 * shape[2] + 2 * j for j in range(c)),
        prod(13 + 2 * j for j in range(b + c)),
    )


@cache
def angular(shape, bcounts, ccounts):
    # A first current on B contributes w_B-delta_{B,P} on mean P.
    # Likewise for C. Both first currents have no other constant terms.
    b0 = sum(bcounts) - bcounts[1]
    c0 = sum(ccounts) - ccounts[2]
    return sum(
        (-1) ** (bcounts[1] + ccounts[2] - i - j)
        * comb(bcounts[1], i)
        * comb(ccounts[2], j)
        * moment(shape, b0 + i, c0 + j)
        for i in range(bcounts[1] + 1)
        for j in range(ccounts[2] + 1)
    )


def add_one(counts, index):
    return tuple(n + int(i == index) for i, n in enumerate(counts))


def coordinate_bound(shape, raw, group):
    normalized = physical = 1
    m = tuple(n - 1 for n in shape)
    for degree in range(9):
        # Each degree-eight term has a two-by-two current minor, whose
        # conductance quadratic cancels over the integers.
        for six in compositions(degree - 2 if degree == 8 else degree, 6):
            b, c = six[:3], six[3:]
            if shape == (1, 6, 6) and max(sum(b), sum(c)) + (degree == 8) > 4:
                continue
            pairs = combinations(range(3), 2) if degree == 8 else [None]
            for pair in pairs:
                if pair is None:
                    bb, cc = b, c
                    angular_value = angular(shape, b, c)
                else:
                    i, j = pair
                    bb, cc = add_one(b, i), add_one(c, j)
                    angular_value = angular(shape, bb, cc) - angular(
                        shape, add_one(b, j), add_one(c, i)
                    )
                d = tuple(i + j for i, j in zip(bb, cc))
                if any(di and not mi for mi, di in zip(m, d)):
                    continue
                value = angular_value * prod(gap_factor(mi, di) for mi, di in zip(m, d))
                normalized = lcm(normalized, value.denominator)
                physical = lcm(physical, (value / raw / group).denominator)
    return normalized, physical


def component_denominator(shape, raw, group):
    bounds = [coordinate_bound(s, raw, group) for s in set(permutations(shape))]
    return tuple(gcd(*values) for values in zip(*bounds))


if __name__ == "__main__":
    expected = (
        293772917487329280000,
        489621529145548800000,
        244810764572774400000,
        979243058291097600000,
        489621529145548800000,
    )
    for (shape, raw, group), want in zip(SHAPES, expected):
        norm, den = component_denominator(shape, raw, group)
        assert den == want and norm % 11 and den % 121
        print(shape, "normalized", norm, "coefficient with symmetry", den, flush=True)
