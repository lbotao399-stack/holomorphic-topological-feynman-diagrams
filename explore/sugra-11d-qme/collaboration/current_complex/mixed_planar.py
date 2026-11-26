"""Exact planar marked CE and induced divergence for one triangular perturbation."""

from collections import Counter, defaultdict
from functools import cache
from itertools import combinations
from fractions import Fraction as F
from math import gcd, lcm
from functools import reduce


def weight(v):
    return 2 * v[2] - v[1] - 2


def basis(g, wh):
    states = defaultdict(list)

    def build(left, normal, minslot):
        for a in range(left + 3):
            marked = (1, left, a)
            if sum(map(weight, normal)) + weight(marked) == wh:
                states[len(normal) + 1].append(tuple(normal + [marked]))
        for h in range(1, left + 2):
            for i in range(h + 3):
                nxt = (0, h, i)
                if nxt > minslot:
                    build(left - h, normal + [nxt], nxt)

    build(g, [], (-1, 0, 0))
    return states


def diff(word):
    out = Counter()
    for i, j in combinations(range(len(word)), 2):
        a, h, u = word[i]
        b, t, v = word[j]
        c = u * (t + 2) - v * (h + 2)
        if not c:
            continue
        merged = (a + b, h + t, u + v - 1)
        if h + t < -1:
            continue
        rest = [word[k] for k in range(len(word)) if k not in (i, j)]
        if merged in rest:
            continue
        sign = (-1) ** (i + j + sum(v < merged for v in rest))
        out[tuple(sorted([merged, *rest]))] += sign * c
    return out


def rank(words):
    piv = {}
    for word in words:
        row = {k: F(v) for k, v in diff(word).items() if v}
        while row:
            key = min(row)
            v = row[key]
            if key not in piv:
                piv[key] = {k: t / v for k, t in row.items()}
                break
            for k, b in piv[key].items():
                row[k] = row.get(k, F(0)) - v * b
                if not row[k]:
                    del row[k]
    return len(piv)


def result(g, hw):
    out = []
    for wh in (hw, hw + 2):
        states = basis(g, wh)
        dims = Counter({q: len(co) for q, co in states.items()})
        ranks = Counter({q: rank(co) for q, co in states.items()})
        out.append((dims, ranks))
    rows = []
    for q in range(1, g + 2):
        d = out[0][0][q] - out[1][0][q]
        r = out[0][1][q] - out[1][1][q]
        incoming = out[0][1][q + 1] - out[1][1][q + 1]
        if d:
            rows.append((q, d, r, d - r - incoming))
    return rows


def normalize(t):
    if len(set(t)) != len(t):
        return None, 0
    return tuple(sorted(t)), (-1) ** sum(
        a > b for i, a in enumerate(t) for b in t[i + 1 :]
    )


def linear(word, kind):
    out = Counter()
    for j, (m, h, i) in enumerate(word):
        if kind == "raise":
            new = (m, h, i + 1)
            factor = h + 2 - i
        elif kind == "dx":
            new = (m, h - 1, i - 1)
            factor = i
        else:
            new = (m, h - 1, i)
            factor = h + 2 - i
        if not factor or (kind != "raise" and new[1] < (1 if not m else -1)):
            continue
        key, sign = normalize((*word[:j], new, *word[j + 1 :]))
        if sign:
            out[key] += sign * factor
    return {k: v for k, v in out.items() if v}


def subtract(v, m, a):
    for k, val in m.items():
        v[k] = v.get(k, F(0)) - a * val
        if not v[k]:
            del v[k]


def lead(v, p):
    while v:
        key = min(v)
        c = v[key]
        if key not in p:
            return key, c
        subtract(v, p[key], c)
    return None, None


def invariant_cycles(states):
    piv = {}
    cycles = []
    for word in states:
        image = {(0, k): F(v) for k, v in diff(word).items() if v}
        image.update({(1, k): F(v) for k, v in linear(word, "raise").items()})
        chain = {word: F(1)}
        while image:
            k = min(image)
            a = image[k]
            if k not in piv:
                piv[k] = (
                    {u: v / a for u, v in image.items()},
                    {u: v / a for u, v in chain.items()},
                )
                break
            p, c = piv[k]
            subtract(image, p, a)
            subtract(chain, c, a)
        else:
            cycles.append(chain)
    return cycles


def graph_test():
    # Three unmarked planar fields of Taylor orders 2, 2, 3 and marked order 3.
    states = (
        word
        for word in basis(6, 0)[4]
        if tuple(slot[1] for slot in word) == (1, 1, 2, 2)
    )
    cycles = invariant_cycles(states)
    boundaries = {}
    for word in basis(6, 0)[5]:
        column = {k: F(v) for k, v in diff(word).items() if v}
        k, a = lead(column, boundaries)
        if k is not None:
            boundaries[k] = {j: c / a for j, c in column.items()}
    classes = [
        cycle for cycle in cycles if lead(cycle.copy(), boundaries)[0] is not None
    ]
    assert len(classes) == 1
    cycle = classes[0]
    scale = lcm(*(c.denominator for c in cycle.values()))
    coefficients = {w: int(c * scale) for w, c in cycle.items()}
    divisor = reduce(gcd, coefficients.values())
    if coefficients[min(coefficients)] < 0:
        divisor = -divisor
    result = {w: c // divisor for w, c in coefficients.items()}
    assert len(result) == 24
    return result


def test():
    q = 4
    cycles = invariant_cycles(basis(6, 0)[q])
    # Current chains at each coefficient weight. Closure after projection is tested separately.
    out_piv = {}
    for typ, weight in [("dx", -1), ("dy", 1)]:
        for w in basis(5, weight)[q + 1]:
            v = {(typ, k): F(c) for k, c in diff(w).items() if c}
            k, a = lead(v, out_piv)
            if k is not None:
                out_piv[k] = {j: x / a for j, x in v.items()}
    rank = 0
    for c in cycles:
        v = Counter()
        boundary = Counter()
        for w, a in c.items():
            for typ in ["dx", "dy"]:
                for k, x in linear(w, typ).items():
                    v[typ, k] += a * x
                    for kk, xx in diff(k).items():
                        boundary[typ, kk] += a * x * xx
        assert not any(boundary.values())
        v = {k: a for k, a in v.items() if a}
        k, a = lead(v, out_piv)
        if k is not None:
            rank += 1
            out_piv[k] = {j: x / a for j, x in v.items()}
    assert (len(cycles), rank) == (9, 1)
    print(
        "invariant scalar four-cycles",
        len(cycles),
        "dual divergence rank on homology",
        rank,
    )


if __name__ == "__main__":
    scalar = result(6, 0)
    current = result(5, 1)
    assert scalar == [
        (2, 1, 0, 0),
        (3, 6, 1, 0),
        (4, 14, 5, 1),
        (5, 9, 8, 0),
        (6, 1, 1, 0),
    ]
    assert current == [(2, 2, 0, 0), (3, 11, 2, 1), (4, 15, 8, 1), (5, 6, 6, 0)]
    print(
        "planar marked scalar/current (arity, dimension, boundary rank, homology)",
        scalar,
        current,
    )
    test()
    cycle = graph_test()
    print(
        "planar graph probe: terms",
        len(cycle),
        "largest coefficient",
        max(map(abs, cycle.values())),
    )
