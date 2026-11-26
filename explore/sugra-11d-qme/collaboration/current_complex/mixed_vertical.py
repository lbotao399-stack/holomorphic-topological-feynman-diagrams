"""Three-dimensional relative CE characters for one mixed mode."""

from itertools import permutations, combinations
from collections import Counter
from fractions import Fraction as F
from functools import cache


def comps(d):
    for i in range(d + 1):
        for j in range(d - i + 1):
            yield (i, j, d - i - j)


def shift(a, i, j):
    b = list(a)
    b[i] -= 1
    b[j] += 1
    return tuple(b)


def basis(d):
    states = []
    for a in comps(d):
        for i in range(3):
            if i == 2 and a[i] > 0:
                continue
            # Div of z^a e_i is a_i z^(a-e_i); cancel with the linked e2 monomial.
            # These are ordinary monomials; division by a_2+1 cancels the divergence.
            pol = {(a, i): F(1)}
            if a[i]:
                pol[shift(a, i, 2), 2] = -F(a[i], a[2] + 1)
            states.append((a, i, pol))
    return states


from collections import defaultdict


@cache
def jets(h):
    return basis(h + 1)


@cache
def general_bracket(m, n):
    h, i = m
    t, j = n
    a = jets(h)[i][2]
    b = jets(t)[j][2]
    slots = {(p, l): s for s, (p, l, _) in enumerate(jets(h + t))}
    out = Counter()
    for left, right, sgn in [(a, b, 1), (b, a, -1)]:
        for (u, axis), x in left.items():
            for (v, k), y in right.items():
                if v[axis]:
                    w = tuple(
                        r + s - int(j == axis) for j, (r, s) in enumerate(zip(u, v))
                    )
                    out[w, k] += sgn * x * y * v[axis]
    return {(h + t, slots[s]): c for s, c in out.items() if c and s in slots}


@cache
def ws(h, i):
    a, j, _ = jets(h)[i]
    return tuple(int(j == k) - a[k] for k in range(3))


def parts(g, n, m=1):
    if not n:
        if not g:
            yield ()
        return
    for i in range(m, g // n + 1):
        for t in parts(g - i, n - 1, i):
            yield (i, *t)


@cache
def groups(h, n):
    d = defaultdict(list)
    for comb in combinations(range(len(jets(h))), n):
        weight = tuple(sum(ws(h, i)[k] for i in comb) for k in range(3))
        d[weight].append(tuple((h, i) for i in comb))
    return d


def slots(g, q, w):
    for p in parts(g, q):
        groups_list = [groups(h, n) for h, n in Counter(p).items()]

        def expand(rest, weight, word):
            if len(rest) == 1:
                for v in rest[0].get(weight, ()):
                    yield word + v
            else:
                for u, words in rest[0].items():
                    t = tuple(a - b for a, b in zip(weight, u))
                    for v in words:
                        yield from expand(rest[1:], t, word + v)

        yield from expand(groups_list, w, ())


@cache
def ranks(g, q, w, prime):
    import numpy as np

    rows = {}
    cols = []
    for word in slots(g, q, w):
        col = Counter()
        for i, j in combinations(range(q), 2):
            rest = word[:i] + word[i + 1 : j] + word[j + 1 :]
            for k, v in general_bracket(word[i], word[j]).items():
                if k in rest:
                    continue
                ind = rows.setdefault(tuple(sorted((k, *rest))), len(rows))
                col[ind] += (-1) ** (i + j + sum(a < k for a in rest)) * v
        cols.append(col)
    mat = np.zeros((len(rows), len(cols)), dtype=np.int32)
    for j, col in enumerate(cols):
        for i, v in col.items():
            mat[i, j] = v.numerator * pow(v.denominator, -1, prime) % prime
    rank = 0
    for j in range(len(cols)):
        nz = np.flatnonzero(mat[rank:, j])
        if not len(nz):
            continue
        i = rank + int(nz[0])
        mat[[rank, i]] = mat[[i, rank]]
        a = pow(int(mat[rank, j]), -1, prime)
        mat[rank] = (mat[rank] * a) % prime
        below = np.flatnonzero(mat[rank + 1 :, j]) + rank + 1
        mat[below, j:] = (
            mat[below, j:] - mat[below, j, None] * mat[rank, None, j:]
        ) % prime
        rank += 1
    return len(cols), rank


def hdim(g, q, hw, prime):
    a = b = c = 0
    rho = (2, 1, 0)
    for p in permutations(rho):
        w = tuple(x + y - z for x, y, z in zip(hw, rho, p))
        # Weyl signs are measured with respect to rho in decreasing order.
        face = (-1) ** sum(i < j for k, i in enumerate(p) for j in p[k + 1 :])
        t, r = ranks(g, q, w, prime)
        _, s = ranks(g, q + 1, w, prime) if q + 1 <= g else (0, 0)
        a += face * t
        b += face * r
        c += face * s
    return a, b, c, a - b - c


if __name__ == "__main__":
    expected = {
        (4, 4, (2, -3, -3)): (1, 1, 0, 0),
        (4, 4, (1, -2, -3)): (4, 4, 0, 0),
        (5, 4, (1, -3, -3)): (10, 6, 3, 1),
        (5, 4, (0, -2, -3)): (16, 9, 6, 1),
        (5, 5, (1, -3, -3)): (3, 3, 0, 0),
        (5, 5, (0, -2, -3)): (6, 6, 0, 0),
        (6, 4, (0, -3, -3)): (22, 8, 14, 0),
        (6, 4, (-1, -2, -3)): (27, 8, 19, 0),
        (6, 6, (0, -3, -3)): (7, 7, 0, 0),
        (6, 5, (0, -3, -3)): (21, 14, 7, 0),
        (6, 5, (-1, -2, -3)): (24, 19, 5, 0),
        (7, 6, (-1, -3, -3)): (28, 22, 6, 0),
        (7, 5, (-1, -3, -3)): (51, 29, 22, 0),
        (7, 3, (-2, -2, -3)): (4, 0, 4, 0),
        (7, 4, (-2, -2, -3)): (21, 4, 17, 0),
        (7, 5, (-2, -2, -3)): (33, 17, 15, 1),
        (7, 6, (-2, -2, -3)): (18, 15, 3, 0),
        (7, 7, (-2, -2, -3)): (3, 3, 0, 0),
    }
    for prime in (251, 257):
        for key, value in expected.items():
            actual = hdim(*key, prime)
            assert actual == value, (key, actual)
        print(
            "prime",
            prime,
            "(dimension, outgoing rank, incoming rank, homology)",
            list(expected.values()),
        )
