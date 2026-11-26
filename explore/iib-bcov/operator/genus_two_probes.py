"""Certify unisolvent eight-primary transverse tests of the genus-two BCOV potential."""
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, prod
from collections import Counter

from array import array


def det(a):
    rows = [list(row) for row in a]
    denominator = sign = 1
    for column in range(len(rows) - 1):
        pivot = next((i for i in range(column, len(rows)) if rows[i][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            sign = -sign
        pivot_value = rows[column][column]
        for i in range(column + 1, len(rows)):
            for j in range(column + 1, len(rows)):
                rows[i][j] = (pivot_value * rows[i][j] - rows[i][column] * rows[column][j]) // denominator
        denominator = pivot_value
    return sign * rows[-1][-1]


def partitions(n, least=1):
    if not n:
        yield []
    for size in range(least, n + 1):
        for rest in partitions(n - size, size):
            yield [size, *rest]


def character_multiplicities(partition):
    trivial = permutation = Fraction()
    for cycle in partitions(8):
        degree = max(partition) + len(partition) - 1
        h = [1, -1] + [0] * (degree - 1)
        for length in cycle:
            for j in range(length, degree + 1):
                h[j] += h[j - length]
        value = det([[h[n - i + j] if n - i + j >= 0 else 0
                      for j in range(len(partition))] for i, n in enumerate(partition)])
        factor = prod(n**m * factorial(m) for n, m in Counter(cycle).items())
        trivial += Fraction(value, factor)
        permutation += Fraction(cycle.count(1) * value, factor)
    return trivial, permutation


def momentum_pair(a, b):
    rows = [tuple(int(i == j) for i in range(5)) for j in range(5)]
    rows.extend([a, b])
    rows.append(tuple(-sum(col) for col in zip(*rows)))
    return tuple(rows)


def orbit_data():
    subsets = list(combinations(range(8), 5))
    subset_number = {sum(1 << i for i in s): j for j, s in enumerate(subsets)}
    signed = {}
    for subset in subsets:
        at = 1 + subset_number[sum(1 << i for i in subset)]
        for ordered in permutations(subset):
            odd = sum(ordered[i] > ordered[j] for i,j in combinations(range(5), 2)) % 2
            signed[ordered] = -at if odd else at
    order = list(permutations(range(8)))
    image = [array('b', (signed[tuple(g[i] for i in subset)] for g in order)) for subset in subsets]
    return subsets, image, [g[0]*8+g[1] for g in order]


def potential_samples(p, data, monomials):
    # Average a product of four five-row minors times p_1 in slot 0 over all 8! relabelings.
    subsets, image, pair = data
    minors = [0, *(det([p[i] for i in subset]) for subset in subsets)]
    signed = minors + [-x for x in reversed(minors[1:])]
    d = [array('i', (signed[x] for x in col)) for col in image]
    out = []
    for slots in monomials:
        weights = [0] * 64
        for label, a, b, c, e in zip(pair, *(d[i] for i in slots)):
            weights[label] += a*b*c*e
        out.append([sum(weights[8*i+j]*p[j][k] for j in range(8))
                    for i in range(8) for k in range(5)])
    return out


# All tests use the top four-vector at leg 7 and coordinate vectors at legs 0 and j.
MOMENTA = momentum_pair((1, -1, 1, 0, -1), (0, 1, 1, -1, -1))
TEST_SLOTS = (
    (1, 1, 0), (1, 1, 2), (1, 1, 4), (1, 2, 0), (1, 4, 0),
    (2, 1, 0), (2, 1, 1), (2, 1, 4),
    (3, 1, 0), (3, 1, 1), (3, 1, 2), (3, 1, 4),
    (4, 1, 0), (4, 1, 1))
ORBITS = (
    (1, 7, 40, 43), (24, 33, 42, 52), (24, 37, 42, 55), (8, 13, 14, 29),
    (5, 15, 27, 31), (1, 24, 42, 55), (19, 25, 37, 46), (8, 22, 41, 41),
    (19, 38, 40, 54), (14, 22, 39, 39), (11, 34, 46, 50), (6, 6, 37, 46),
    (0, 5, 26, 40), (14, 31, 34, 48))


def evaluate_potential(potential):
    """Fourteen rational linear functionals of an 8 by 5 polynomial potential at MOMENTA."""
    assert len(potential) == 8 and all(len(row) == 5 for row in potential)
    result = []
    for j, first, second in TEST_SLOTS:
        assert MOMENTA[0][first] == MOMENTA[j][second] == 0
        result.append((potential[0][first] - potential[7][first]) * MOMENTA[7][second]
                      - (potential[j][second] - potential[7][second]) * MOMENTA[7][first])
    return tuple(result)


SECOND_MOMENTA = momentum_pair((1, 2, -1, 0, -1), (-2, 0, 1, -1, 2))
STRONG_ROWS = (0, 1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 32, 33)
EXTRA_ORBITS = ((6, 19, 25, 46), (2, 35, 39, 45), (22, 27, 37, 46),
                (4, 8, 34, 46), (24, 24, 38, 41), (11, 12, 15, 36))


def component_tests(p, potential):
    total = [sum(col) for col in zip(*potential)]
    values = []
    for leg, k in zip(potential, p):
        pivot = next(i for i, x in enumerate(k) if x)
        for i in range(5):
            if pivot != i:
                values.append(k[pivot] * (10 * leg[i] - total[i])
                              - k[i] * (10 * leg[pivot] - total[pivot]))
    return values


def evaluate_strong(first_potential, second_potential):
    """Test 10 E_i minus the sum of all eight E rows against transverse waves."""
    rows = component_tests(MOMENTA, first_potential) + component_tests(SECOND_MOMENTA, second_potential)
    return tuple(rows[i] for i in STRONG_ROWS)


def check():
    assert character_multiplicities((5, 4, 4, 4, 4)) == (6, 32)
    assert character_multiplicities((4, 4, 4, 4, 4)) == (4, 12)
    data = orbit_data()
    samples = potential_samples(MOMENTA, data, ORBITS + EXTRA_ORBITS)
    tests = [evaluate_potential([col[i:i+5] for i in range(0, 40, 5)]) for col in samples[:14]]
    residue = det(tests) % 1009
    assert residue == 96
    second = potential_samples(SECOND_MOMENTA, data, ORBITS + EXTRA_ORBITS)
    strong = [evaluate_strong([a[i:i+5] for i in range(0, 40, 5)],
                              [b[i:i+5] for i in range(0, 40, 5)])
              for a, b in zip(samples, second)]
    strong_residue = det(strong) % 1009
    assert strong_residue == 927
    print('Twenty-probe determinant modulo 1009:', strong_residue)
    print('Potential, common and longitudinal multiplicities: 32, 6, 12.')
    print('Fourteen-probe determinant modulo 1009:', residue)


if __name__ == '__main__':
    check()
