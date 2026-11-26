"""Direct Cartan-weight Koszul matrices for the sl2 zero-holomorphic-jet sector."""

from itertools import combinations
from quadratic import exact_rank

# Coordinates in the order E, H, F with [H,E]=2E, [H,F]=-2F, [E,F]=H.
BRACKET = (((1, 0, 2), (0, 1, -2)),
           ((0, 2, 1), (2, 0, -1)),
           ((2, 1, 2), (1, 2, -2)))


def triples(total):
    for a in range(total + 1):
        for b in range(total - a + 1):
            yield a, b, total - a - b


def basis(p, q, degree, charge):
    if min(p, q) < degree:
        return []
    return [(x, y, b) for x in triples(p - degree) for y in triples(q - degree)
            for b in combinations(range(3), degree)
            if x[0] - x[2] + y[0] - y[2] + sum(1 - i for i in b) == charge]


def differential(monomial):
    x, y, b = monomial
    row = {}
    for pos, color in enumerate(b):
        for i, j, coefficient in BRACKET[color]:
            xx = tuple(v + int(k == i) for k, v in enumerate(x))
            yy = tuple(v + int(k == j) for k, v in enumerate(y))
            target = (xx, yy, b[:pos] + b[pos + 1:])
            row[target] = row.get(target, 0) + (-1) ** pos * coefficient
    return row


def homology(p, q, charge):
    spaces = [basis(p, q, k, charge) for k in range(4)]
    ranks = [exact_rank(map(differential, space)) for space in spaces] + [0]
    return tuple(len(space) - ranks[k] - ranks[k + 1] for k, space in enumerate(spaces))


if __name__ == "__main__":
    total = 0
    for p in range(7):
        for q in range(7):
            zero, raising = homology(p, q, 0), homology(p, q, 1)
            invariant = tuple(a - b for a, b in zip(zero, raising))
            expected = (int((p + q) % 2 == 0), int(min(p, q) >= 1 and (p + q) % 2 == 1), 0, 0)
            assert invariant == expected, ((p, q), invariant, expected)
            total += sum(invariant)
    print(f"49 sl2 bidegrees, {total} relative classes; all four b-degrees checked")
