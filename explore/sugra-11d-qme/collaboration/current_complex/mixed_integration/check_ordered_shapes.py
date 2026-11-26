"""Independent cubature replay of each ordered mixed fixture over F_1009."""

from collections import defaultdict
from functools import cache
from fractions import Fraction
from itertools import product
from math import factorial, prod

import numpy as np

from exact_scalar import MASK, components, fixture
from ordered_shapes import numerator, paths
from theta_denominators import SHAPES, compositions

P = 1009


def field(q):
    return q.numerator * pow(q.denominator, -1, P) % P


@cache
def simplex(m, s):
    if not m:
        return {(1,): 1}
    out = defaultdict(int)
    for i in range(s + 1):
        ell = s - i
        h = m + 1 + 2 * ell
        q = Fraction(
            (-1) ** i * factorial(m) * h ** (2 * s + 1),
            4**s * factorial(i) * factorial(m + 1 + 2 * s - i),
        )
        for a in compositions(ell, m + 1):
            t = tuple((1 + 2 * v) * pow(h, -1, P) % P for v in a)
            out[t] = (out[t] + field(q)) % P
    for d in range(2 * s + 2):
        for e in compositions(d, m + 1):
            if tuple(sorted(e)) != e:
                continue
            test = (
                sum(
                    w * prod(pow(x, n, P) for x, n in zip(t, e)) for t, w in out.items()
                )
                % P
            )
            exact = Fraction(factorial(m) * prod(map(factorial, e)), factorial(m + d))
            assert test == field(exact)
    return out


def rule(shape):
    # Tensor products of consecutive rule differences. Level four covers
    # every total-degree-eight power, as the one-path rule is exact through 2s+1.
    out = defaultdict(int)
    for level in (2, 3, 4):
        for orders in compositions(level, 3):
            tables = [simplex(n - 1, s).items() for n, s in zip(shape, orders)]
            for selected in product(*tables):
                pts, ws = zip(*selected)
                out[pts] = (out[pts] + (1 if level != 3 else -2) * prod(ws)) % P
    return {k: v for k, v in out.items() if v}


def angular(shape):
    # Build weights on the integer triangular grid by solving all monomials.
    grid = [pair for d in range(8) for pair in compositions(d, 2)]
    mat = []
    for b, c in grid:
        rhs = Fraction(
            prod(13 - 2 * shape[1] + 2 * j for j in range(b))
            * prod(13 - 2 * shape[2] + 2 * j for j in range(c)),
            prod(13 + 2 * j for j in range(b + c)),
        )
        mat.append([pow(i, b, P) * pow(j, c, P) % P for i, j in grid] + [field(rhs)])
    for j in range(len(grid)):
        k = next(h for h in range(j, len(grid)) if mat[h][j])
        mat[j], mat[k] = mat[k], mat[j]
        mat[j] = [x * pow(mat[j][j], -1, P) % P for x in mat[j]]
        for i in range(len(grid)):
            if i != j:
                coeff = mat[i][j]
                mat[i] = [(x - coeff * y) % P for x, y in zip(mat[i], mat[j])]
    return [(r[0], r[1], a[-1]) for r, a in zip(grid, mat) if a[-1]]


def cubature(poly, data, shape, order):
    resist = rule(shape)
    ps = paths(shape, order)
    means = np.empty((3, len(resist), 5), np.int64)
    weights = np.array(list(resist.values()), dtype=np.int64)
    for ix, points in enumerate(resist):
        for h, (pt, ls) in enumerate(zip(points, ps)):
            q = np.zeros(5, dtype=np.int64)
            val = q.copy()
            for t, i in zip(pt[1:], ls):
                q += data[i]
                val += t * q
            means[h, ix, :] = val % P
    msum = means.sum(axis=0) % P
    left = sorted({k & MASK for k in poly})
    right = sorted({k >> 20 for k in poly})
    l_id = {k: i for i, k in enumerate(left)}
    r_id = {k: i for i, k in enumerate(right)}
    coeff = np.zeros((len(left), len(right)), dtype=np.int64)
    for key, v in poly.items():
        coeff[l_id[key & MASK], r_id[key >> 20]] = v % P

    def monomials(v, keys):
        ans = np.ones((len(keys), v.shape[0]), dtype=np.int64)
        for j in range(5):
            axis = np.array([components(k)[j] for k in keys])
            powers = [np.ones(v.shape[0], dtype=np.int64)]
            for _ in range(int(axis.max())):
                powers.append((powers[-1] * v[:, j]) % P)
            ans *= np.array(powers)[axis, :]
            ans %= P
        return ans

    total = 0
    # Every matrix inner product fits exactly in int64, including both outer sums.
    angular_nodes = angular(shape)
    for w_b, w_c, w in angular_nodes:
        first = monomials((w_b * msum - means[1]) % P, left)
        second = monomials((w_c * msum - means[2]) % P, right)
        v = (((coeff @ second) % P) * first).sum(axis=0) % P
        total = (total + w * int(v @ weights)) % P
    return total, len(resist), len(angular_nodes)


if __name__ == "__main__":
    expected = (
        Fraction(4170228116417, 25758708393600),
        Fraction(-1360113875759, 21465590328000),
        Fraction(26094891211, 1022170968000),
        Fraction(5763274313, 161395416000),
        Fraction(3756403411, 151432736000),
    )
    orders = ("0146723589", "0146723589", "0146723589", "6901472358", "0691472358")
    data = fixture()
    for (shape, _, _), s, answer in zip(SHAPES, orders, expected):
        order = tuple(map(int, s))
        primitive = numerator(data, shape, order, False)
        assert primitive == numerator(data, shape, order, True)
        replay, nr, na = cubature(primitive, data[0], shape, order)
        assert replay == field(answer), (shape, replay, field(answer))
        print(
            "shape",
            shape,
            "resistance nodes",
            nr,
            "conductance nodes",
            na,
            "Wick mean",
            replay,
            "mod",
            P,
            flush=True,
        )
