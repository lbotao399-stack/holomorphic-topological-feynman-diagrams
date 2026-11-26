"""Exact separation of the mixed detector's polarization weights across a path cut."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

WEIGHTS = Path(__file__).resolve().parent.parent / "mixed_graph_probe.tsv"
SINGLE_WEIGHTS = WEIGHTS.with_name("mixed_single_weights.tsv")


def key(group, bits, eta):
    return (
        bits & sum(1 << (j - 4) for j in group if 4 <= j <= 8),
        eta if 9 in group else -1,
    )


def integral(v):
    denominator = lcm(*(x.denominator for x in v.values()))
    values = {k: int(x * denominator) for k, x in v.items() if x}
    content = reduce(gcd, map(abs, values.values()))
    return {k: c // content for k, c in values.items()}, Fraction(content, denominator)


def factors(order, prime=None, weights_file=WEIGHTS):
    tables = defaultdict(Counter)
    for line in Path(weights_file).read_text().splitlines():
        values = tuple(map(int, line.split()))
        if len(values) == 3:
            weight, bits, eta = values
            function = 0
        else:
            weight, bits, function, eta = values
        tables[function][key(order[:5], bits, eta), key(order[5:], bits, eta)] += weight
    output = {}
    for function, table in tables.items():
        residual = {
            k: (Fraction(v) if prime is None else v % prime) for k, v in table.items()
        }
        residual = {k: v for k, v in residual.items() if v}
        product = []
        while residual:
            # Favor a sparse row to keep the path components small.
            counts = Counter(a for a, b in residual)
            i, j = min(
                residual,
                key=lambda t: (counts[t[0]], abs(residual[t].numerator).bit_length()),
            )
            pivot = residual[i, j]
            left = {a: v for (a, b), v in residual.items() if b == j}
            right = {b: v for (a, b), v in residual.items() if a == i}
            if prime is None:
                left, scale_l = integral(left)
                right, scale_r = integral(right)
                scale = scale_l * scale_r / pivot
            else:
                scale = pow(pivot, -1, prime)
            for a, u in left.items():
                for b, v in right.items():
                    k = a, b
                    residual[k] = residual.get(k, 0) - scale * u * v
                    if prime is not None:
                        residual[k] %= prime
                    if not residual[k]:
                        del residual[k]
            product.append((scale, left, right))
        rebuilt = Counter()
        for scale, left, right in product:
            for a, u in left.items():
                for b, v in right.items():
                    rebuilt[a, b] += scale * u * v
        difference = (rebuilt[k] - table[k] for k in set(rebuilt) | set(table))
        assert all((x % prime if prime is not None else x) == 0 for x in difference)
        output[function] = product
    return output


if __name__ == "__main__":
    hist = Counter()
    for b in combinations(range(9), 5):
        order = b + tuple(x for x in range(10) if x not in b)
        r = factors(order)
        hist[tuple(len(r[f]) for f in sorted(r))] += 1
    for ranks, count in sorted(hist.items()):
        print("ranks for function selectors 0,1,3", ranks, "partitions", count)
