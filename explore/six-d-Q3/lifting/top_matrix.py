"""Alternating color traces for homogeneous matrices over an exterior algebra."""

from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "s2code"))
from salg import add, const, gen, mul, scale


def word(entries):
    result = const(1)
    for entry in entries:
        result = mul(result, entry)
    return result


def trace_word(matrices):
    n = len(matrices[0])
    total = {}
    for indices in product(range(n), repeat=len(matrices)):
        entries = [matrix[indices[j]][indices[(j + 1) % len(indices)]]
                   for j, matrix in enumerate(matrices)]
        total = add(total, word(entries))
    return total


def cycles(perm):
    seen = set()
    result = []
    for i in range(len(perm)):
        if i in seen:
            continue
        orbit = []
        while i not in seen:
            seen.add(i)
            orbit.append(i)
            i = perm[i]
        result.append(orbit)
    return result


def alternating_trace(matrices, odd, trace_weight=1):
    result = {}
    traces = {}
    for perm in permutations(range(len(matrices))):
        parts = cycles(perm)
        sign = (-1) ** (len(perm) - len(parts))
        slots = [j for part in parts for j in part if odd[j]]
        sign *= (-1) ** sum(a > b for i, a in enumerate(slots) for b in slots[i + 1:])
        factors = []
        for part in parts:
            key = tuple(part)
            if key not in traces:
                traces[key] = scale(trace_word([matrices[j] for j in part]), trace_weight)
            factors.append(traces[key])
        result = add(result, word(factors), sign)
    return result


def example(n, k, diagonal):
    odd = [i % 3 != 0 for i in range(k)]
    matrices = []
    for r in range(k):
        matrix = []
        for i in range(n):
            row = []
            for j in range(n):
                z = 2 + i + 3 * j + r * (j + 1)
                entry = scale(gen((1, r, i, j, 0)), z) if odd[r] else const(z)
                row.append({} if diagonal and i != j else entry)
            matrix.append(row)
        matrices.append(matrix)
    actual = alternating_trace(matrices, odd)
    if diagonal:
        expected = {}
        for labels in permutations(range(n), k):
            expected = add(expected, word([M[i][i] for M, i in zip(matrices, labels)]))
        assert actual == expected
        assert bool(actual) == (k <= n)
    else:
        assert k > n and not actual


def symplectic(rank):
    # M = J^{-1}K for homogeneous skew matrices K, using
    # J = [[0,I],[-I,0]]. The coefficient algebra may contain odd generators.
    n = 2 * rank
    odd = [i % 2 == 0 for i in range(rank + 1)]
    matrices = []
    for r, parity in enumerate(odd):
        skew = [[{} for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                value = 1 + i + r * (1 + j)
                entry = scale(gen((1, r, i, j, 0)), value) if parity else const(value)
                skew[i][j], skew[j][i] = entry, scale(entry, -1)
        matrices.append([[scale(skew[(i + rank) % n][j], -1 if i < rank else 1)
                          for j in range(n)] for i in range(n)])
    assert not alternating_trace(matrices, odd, Fraction(1, 2))


if __name__ == "__main__":
    for n in (2, 3):
        for k in range(1, n + 2):
            example(n, k, True)
        example(n, n + 1, False)
    for rank in (1, 2, 3):
        symplectic(rank)
    print("matrix trace relations at N=2,3 and symplectic half-trace identities at ranks 1,2,3 passed")
