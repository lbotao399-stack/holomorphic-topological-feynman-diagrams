"""Exact SL(2)-invariant CE ranks in arities three and four on traceless matrices."""
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import prod
from pathlib import Path
from random import Random
import sys

from sympy import Matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'deform' / 'code'))
import ce


@cache
def words(k):
    return tuple(t for t in ce.all_taus(k) if not ce.has_singleton(t))


def metric(first, second, k, n):
    """N^k times the Hermitian Gram form with every matrix projected to sl_N."""
    successor = {i: cycle[(j + 1) % len(cycle)]
                 for cycle in second for j, i in enumerate(cycle)}
    value = 0
    for mask in range(1 << k):
        removed = [i for i in range(k) if mask >> i & 1]
        cycles = [[i for i in cycle if i not in removed] for cycle in first]
        empty = cycles.count([])
        predecessor = {i: cycle[j - 1] for cycle in cycles for j, i in enumerate(cycle)}
        predecessor.update((i, i) for i in removed)
        todo = set(range(k))
        count = 0
        while todo:
            i = next(iter(todo))
            while i in todo:
                todo.remove(i)
                i = predecessor[successor[i]]
            count += 1
        value += (-1) ** len(removed) * n ** (k - len(removed) + empty + count)
    return value


@cache
def projection(k, n):
    if n >= k:
        return Matrix.eye(len(words(k)))
    gram = Matrix([[metric(t, u, k, n) for t in words(k)] for u in words(k)])
    return gram[list(gram.rref()[1]), :]


def bracket(p, q):
    return p[0] * q[1] - p[1] * q[0]


def seeds(k, w, central=0):
    representatives = {}
    for t in words(k - central):
        representatives.setdefault(tuple(sorted(map(len, t))), t)
    pairs = list(combinations(range(k - 1), 2))
    return [(t, m) for t in representatives.values()
            for m in combinations_with_replacement(pairs, w)]


@cache
def actions(k, t, central):
    traceless = k - central
    indices = {u: i for i, u in enumerate(words(traceless))}
    orders = (a + b for a, b in product(permutations(range(traceless)),
                                         permutations(range(traceless, k))))
    return [(s, ce.perm_sign(s), indices[ce.relabel(t, [(j,) for j in s])])
            for s in orders]


def evaluate(k, generators, p, central=0):
    brackets = {(i, j): bracket(p[i], p[j]) for i in range(k) for j in range(k)}
    vals = [[0] * len(generators) for t in words(k-central)]
    for col, (t, monomial) in enumerate(generators):
        for s, sign, row in actions(k, t, central):
            vals[row][col] += sign * prod(brackets[s[i], s[j]] for i, j in monomial)
    return Matrix(vals)


def differential(k, generators, p, central=0):
    target = {t: i for i, t in enumerate(words(k+1-central))}
    result = [[0] * len(generators) for t in target]
    for i, j in combinations(range(k+1-central), 2):
        rest = [a for a in range(k+1) if a not in (i, j)]
        merged = tuple(x+y for x, y in zip(p[i], p[j]))
        ev = evaluate(k, generators, [merged] + [p[a] for a in rest], central).tolist()
        for order, sign in (((i, j), 1), ((j, i), -1)):
            labels = [order] + [(a,) for a in rest]
            for t, values in zip(words(k-central), ev):
                dest = result[target[ce.relabel(t, labels)]]
                for col, x in enumerate(values):
                    dest[col] += (-1) ** (i+j) * sign * x
    return Matrix(result)


def balance(p):
    return p + [tuple(-sum(a) for a in zip(*p))]


def independent_indices(matrix):
    return matrix.to_DM().rref()[1]


def dimension(k, w, n, central=0):
    if central > k or not words(k-central):
        return [], []
    generators = seeds(k, w, central)
    samples = [balance([(1, 0), (0, 1)])] if k == 3 else [
        balance([(1, 0), (0, 1), (-z, y)])
        for y in range(w+1) for z in range(w+1-y)]
    evaluation = Matrix.vstack(*(projection(k-central, n)*evaluate(k, generators, p, central) for p in samples))
    return [generators[i] for i in independent_indices(evaluation)], samples


def certify(w, n):
    three, _ = dimension(3, w, n)
    four, samples = dimension(4, w, n)
    # Homogeneity in [01],[02],[12] makes the triangular grid unisolvent.
    r3 = len(independent_indices(Matrix.vstack(*(projection(4, n)*differential(3, three, p)
                                                  for p in samples)))) if three else 0
    # d^2=0 bounds rank d4 from above. At w=2, N>=3 the closed double trace
    # evaluates nontrivially on two orthogonal Cartan elements, unlike all boundaries.
    predicted_h4 = int(w == 2 and n >= 3)
    upper = len(four) - r3 - predicted_h4
    rng = Random(1949)
    rows = Matrix.zeros(0, len(four))
    r4 = 0
    for _ in range(30):
        if r4 == upper:
            break
        p = balance([(rng.randrange(-5, 6), rng.randrange(-5, 6)) for _ in range(4)])
        rows = rows.col_join(projection(5, n)*differential(4, four, p))
        r4 = len(independent_indices(rows))
    assert r4 == upper, (w, n, len(four), r3, r4, upper)
    assert len(three) == r3
    return len(three), len(four), r3, r4, len(four)-r3-r4


if __name__ == '__main__':
    tests = [tuple(map(int, pair.split(','))) for pair in sys.argv[1:]]
    if not tests:
        tests = [(w, n) for w in range(2, 7) for n in (2, 3, 4, 5)]
    for w, n in tests:
        print(w, n, *certify(w, n), flush=True)
