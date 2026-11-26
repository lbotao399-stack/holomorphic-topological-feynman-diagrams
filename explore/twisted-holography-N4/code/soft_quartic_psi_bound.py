"""Upper control for the mixed two-soft fourfold vertex via compact psi intersections."""

from fractions import Fraction
from functools import cache
from itertools import combinations, product
from math import comb, factorial, prod


def odd_factorial(k):
    return prod(range(k, 0, -2))


def intersection(genus, degrees, last=False):
    return virasoro(genus, tuple(sorted(degrees)), last)


@cache
def virasoro(genus, degrees, last):
    n = len(degrees)
    if (
        genus < 0
        or 2 * genus - 2 + n <= 0
        or min(degrees, default=0) < 0
        or sum(degrees) != 3 * genus - 3 + n
    ):
        return Fraction(0)
    if genus == 0:
        return Fraction(factorial(n - 3), prod(map(factorial, degrees)))
    if genus == 1 and degrees == (1,):
        return Fraction(1, 24)

    def integral(g, ds):
        return intersection(g, ds, last)

    k, *other = degrees
    if k == 0:
        return sum(
            integral(genus, other[:i] + [d - 1] + other[i + 1 :])
            for i, d in enumerate(other)
            if d
        )
    if k == 1:
        return (2 * genus - 2 + len(other)) * integral(genus, other)
    if last:
        *other, k = degrees
    result = sum(
        Fraction(odd_factorial(2 * k + 2 * d - 1), odd_factorial(2 * d - 1))
        * integral(genus, other[:i] + [k + d - 1] + other[i + 1 :])
        for i, d in enumerate(other)
    )
    for a in range(k - 1):
        b = k - 2 - a
        boundary = integral(genus - 1, [a, b, *other])
        for mask in range(1 << len(other)):
            left = [a] + [d for i, d in enumerate(other) if mask & (1 << i)]
            right = [b] + [d for i, d in enumerate(other) if not mask & (1 << i)]
            boundary += sum(
                integral(h, left) * integral(genus - h, right) for h in range(genus + 1)
            )
        result += (
            Fraction(odd_factorial(2 * a + 1) * odd_factorial(2 * b + 1), 2) * boundary
        )
    return result / odd_factorial(2 * k + 1)


def cut_intersection(genus, degrees, first, second, last):
    """The two named marks choose the sides of each separating node uniquely."""
    n = len(degrees)
    free = sorted(set(range(n)) - {first, second})
    value = Fraction(0)
    for size in range(len(free) + 1):
        for subset in combinations(free, size):
            left = {first, *subset}
            right = set(range(n)) - left
            for h in range(genus + 1):
                value += intersection(
                    h, [0, *[degrees[i] for i in left]], last
                ) * intersection(genus - h, [0, *[degrees[i] for i in right]], last)
    return value


def upper_bound(links, last):
    # The first link remains an actual Betti form in the last intersection.
    first, second, multiple = links[0]
    powers = [(first, second, multiple - 1), *links[1:]]
    total = Fraction(0)
    correction = Fraction(0)
    for choices in product(*(range(power + 1) for _, _, power in powers)):
        degrees = [0] * 5
        weight = 1
        for k, (i, j, power) in zip(choices, powers):
            degrees[i] += k
            degrees[j] += power - k
            weight *= comb(power, k)
        for mark in (first, second):
            degrees[mark] += 1
            total += weight * intersection(4, degrees, last)
            degrees[mark] -= 1
        correction += weight * cut_intersection(4, degrees, first, second, last)
    # Two soft integrations give 2!; the fourteen single-coordinate cotangent averages each give 1/2.
    divisor = 2**14 * prod(factorial(p) for _, _, p in links)
    return 2 * total / divisor, 2 * (total - correction) / divisor


def check():
    links = [(0, 1, 4), (1, 4, 2), (2, 4, 4), (3, 4, 4)]
    for last in (False, True):
        for g in range(1, 5):
            assert intersection(g, [3 * g - 2], last) == Fraction(
                1, 24**g * factorial(g)
            )
        # Mbar_0,4 has two boundary points separating a specified pair.
        assert cut_intersection(0, [0] * 4, 0, 1, last) == 2
        all_psi, one_betti = upper_bound(links, last)
        assert all_psi == Fraction(1535473758853, 157803970402713600)
        assert one_betti == Fraction(379101134941, 39450992600678400)
    print("Mixed fourfold coefficient upper bound:", one_betti)


if __name__ == "__main__":
    check()
