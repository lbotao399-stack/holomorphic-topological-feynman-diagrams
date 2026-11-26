"""Exact finite-field product cubature for the normalized theta_166 measure."""

from collections import defaultdict
from functools import cache
from math import factorial, prod

P = 1009


def compositions(d, n):
    if n == 1:
        yield (d,)
        return
    for i in range(d + 1):
        for r in compositions(d - i, n - 1):
            yield (i, *r)


def inv(x):
    return pow(x, -1, P)


@cache
def simplex(s):
    # Dimension five. The proof by exact monomial integration below includes every orbit.
    ans = defaultdict(int)
    for i in range(s + 1):
        l = s - i
        h = 6 + 2 * l
        w = (
            (-1) ** i
            * factorial(5)
            * pow(h, 2 * s + 1, P)
            * inv(4**s * factorial(i) * factorial(6 + 2 * s - i))
        )
        for x in compositions(l, 6):
            node = tuple((1 + 2 * k) * inv(h) % P for k in x)
            ans[node] = (ans[node] + w) % P
    return dict(ans)


def smolyak(level=5):
    nodes = defaultdict(int)
    for layer, sign in [(level, 1), (level - 1, -1)]:
        for r in range(layer + 1):
            for x, wx in simplex(r).items():
                for y, wy in simplex(layer - r).items():
                    nodes[x, y] = (nodes[x, y] + sign * wx * wy) % P
    return {k: v for k, v in nodes.items() if v}


def check_simplex():
    for r in range(6):
        nodes = simplex(r)
        for d in range(2 * r + 2):
            for a in compositions(d, 6):
                if sorted(a) != list(a):
                    continue
                exact = (
                    factorial(5)
                    * prod(factorial(i) for i in a)
                    * inv(factorial(5 + d))
                    % P
                )
                assert (
                    sum(
                        w * prod(pow(x, j, P) for x, j in zip(t, a))
                        for t, w in nodes.items()
                    )
                    % P
                    == exact
                )


def solve(a, b):
    rows = [[*r, x] for r, x in zip(a, b)]
    for i in range(len(rows)):
        h = next(k for k in range(i, len(rows)) if rows[k][i] % P)
        rows[i], rows[h] = rows[h], rows[i]
        c = inv(rows[i][i])
        rows[i] = [x * c % P for x in rows[i]]
        for j in range(len(rows)):
            if j == i:
                continue
            c = rows[j][i]
            rows[j] = [(v - c * u) % P for v, u in zip(rows[j], rows[i])]
    return [row[-1] for row in rows]


def conductance():
    # w_B=r z, w_C=r(1-z) with r~Beta(1,11/2), z~Beta(1/2,1/2).
    # Search for split orthogonal polynomials in this finite field.
    def moment(n):
        return 2**n * factorial(n) * inv(prod(13 + 2 * i for i in range(n))) % P

    def monic(n):
        return [
            *solve(
                [[moment(i + j) for j in range(n)] for i in range(n)],
                [-moment(n + i) for i in range(n)],
            ),
            1,
        ]

    polys = [monic(n) for n in [5, 6, 7]]
    ev = [
        [sum(pow(x, i, P) * a for i, a in enumerate(poly)) % P for x in range(P)]
        for poly in polys
    ]

    def ab():
        for x in range(40):
            for y in range(x):
                if (ev[0][x] * ev[1][y] - ev[1][x] * ev[0][y]) % P:
                    yield solve(
                        [[ev[0][x], ev[1][x]], [ev[0][y], ev[1][y]]],
                        [-ev[2][x], -ev[2][y]],
                    )

    for a, b in ab():
        radial = [
            x for x in range(P) if (ev[2][x] + a * ev[0][x] + b * ev[1][x]) % P == 0
        ]
        if len(radial) == 7:
            break
    else:
        raise ValueError("no split radial rule")
    wr = solve(
        [[pow(t, j, P) for t in radial] for j in range(7)],
        [moment(j) for j in range(7)],
    )
    zeta = next(
        x
        for x in range(2, P)
        if pow(x, 12, P) == P - 1 and pow(x, 8, P) != 1 and pow(x, 24, P) == 1
    )
    zs = [
        (2 + pow(zeta, 2 * i + 1, P) + pow(zeta, -2 * i - 1, P)) * inv(4) % P
        for i in range(6)
    ]
    pts = [
        (r * z % P, r * (1 - z) % P, v * inv(6) % P)
        for r, v in zip(radial, wr)
        for z in zs
    ]
    for d in range(12):
        for i in range(d + 1):
            j = d - i
            exact = (
                prod(range(1, 2 * i, 2))
                * prod(range(1, 2 * j, 2))
                * inv(prod(range(13, 13 + 2 * d, 2)))
                % P
            )
            assert sum(v * pow(x, i, P) * pow(y, j, P) for x, y, v in pts) % P == exact
    return pts


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("part", choices=["rule", "conductance", "check"])
    part = parser.parse_args().part
    if part == "check":
        check_simplex()
        print(
            "resistance nodes", len(smolyak()), "conductance nodes", len(conductance())
        )
    elif part == "rule":
        for (x, y), w in smolyak().items():
            print(*x, *y, w)
    else:
        for row in conductance():
            print(*row)
