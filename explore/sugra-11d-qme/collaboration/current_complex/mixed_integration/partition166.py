"""Export the separated path factors, order signs and prefixes for exact contraction."""

from collections import Counter, defaultdict
from itertools import islice, permutations
from math import prod
from exact_scalar import MASK, N, collect, epsilon, fixture, graph_factors, multiply


def determinant(momenta, order):
    def det(group):
        return sum(
            epsilon(p) * prod(momenta[j][a] for j, a in zip(group, p))
            for p in permutations(range(N))
        )

    group = sorted(order[:5]) + sorted(order[5:])
    return epsilon(group) * det(group[:5]) * det(group[5:])


def export(
    data,
    order,
    permute=True,
    limit=120,
    lie=True,
    polarizations=None,
    modulus=0,
    det=None,
):
    print("D", determinant(data[0], order) if det is None else det)
    for h in (1, 2):
        group = order[:5] if h == 1 else order[5:]
        orders = islice(permutations(group), limit) if permute else [group]
        for n, path in enumerate(orders):
            labels = path + order[5:] if h == 1 else order[:5] + path
            gamma_sign = (-1) ** (9 - labels.index(9)) if 9 in group else 1
            print("R", h, n, *(v for i in path for v in data[0][i]))
            sums = Counter()
            inputs = polarizations[h - 1] if polarizations else [(1, data[1])]
            for weight, polar in inputs:
                for flow, (_, factors, _) in enumerate(
                    graph_factors((data[0], polar), labels, (h,), False, lie)
                ):
                    for i, row in enumerate(factors[h]):
                        for j, polynomial in enumerate(row):
                            for key, coeff in polynomial.items():
                                if coeff:
                                    sums[
                                        flow, i, j, key & MASK if h == 1 else key >> 20
                                    ] += weight * gamma_sign * coeff
            for key, value in sums.items():
                c = value % modulus if modulus else value
                if c:
                    print("P", h, *key, n, c)
    for flow, (wick, _, (left, right)) in enumerate(graph_factors(data, order, ())):
        ends = defaultdict(dict)
        for (i, j, k), x in left.items():
            for (h, s, t), y in right.items():
                if h == i:
                    collect(ends[j, s, k, t], multiply(x, y), wick)
        for (j, s, k, t), poly in ends.items():
            for key, c in poly.items():
                if c:
                    print("E", flow, j, s, k, t, key & MASK, key >> 20, c)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--order", default="0146723589")
    parser.add_argument("--single", action="store_true")
    parser.add_argument("--limit", type=int, default=120)
    args = parser.parse_args()
    order = tuple(map(int, args.order))
    assert sorted(order) == list(range(10))
    export(fixture(), order, not args.single, args.limit)
