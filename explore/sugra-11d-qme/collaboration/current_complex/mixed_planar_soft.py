"""The sharp planar soft order in the marked current three-class."""

from itertools import permutations
from fractions import Fraction
from math import factorial
from sympy import Matrix, Poly, Rational, symbols

from mixed_current_lift import lift
from mixed_planar import basis, diff, linear
from mixed_planar_stream import stream


def current_coordinates(q):
    return [(k, w) for k, z in (("dx", -1), ("dy", 1)) for w in basis(5, z)[q]]


def boundaries():
    words = current_coordinates(2)
    relations = []
    for kind, z in (("dx", -1), ("dy", 1)):
        for w in basis(5, z - 2)[2]:
            col = {(kind, k): v for k, v in linear(w, "raise").items()}
            if kind == "dy":
                col["dx", w] = 1
            relations.append([col.get(i, 0) for i in words])
    result = []
    for v in Matrix(relations).nullspace():
        source = dict(zip(words, v))
        result.append(
            {
                (k, w): sum(source.get((k, t), 0) * n for t, n in diff(w).items())
                for k, w in current_coordinates(3)
            }
        )
    assert len(result) == 2
    return result


def symbol(cochain, a, b):
    # Plane waves on the two pure Hamiltonians have momenta (1,0),(0,1).
    # (a,b) is the sum of the three momenta at the insertion point.
    momenta = ((1, 0), (0, 1), (a - 1, b - 1))
    value = 0
    for w, c in cochain.items():
        for i, j, sign in ((0, 1, 1), (1, 0, -1)):
            term = c * sign
            for n, (_, h, x) in zip((i, j, 2), w):
                y = h + 2 - x
                p, q = momenta[n]
                term *= Rational(1, factorial(x) * factorial(y)) * p**x * q**y
            value += term
    return Poly(value, a, b)


def verify():
    scalar, current, _ = stream()
    a, b = symbols("a b")
    potential = symbol(scalar, a, b)
    assert potential == Poly(
        -(a * a * b * b - a * a * b + a * a - a * b * b - a * b + b * b) / 24, a, b
    )
    components = []
    for kind, multiplier in (("dx", -b), ("dy", a)):
        column = []
        for c in [*boundaries(), current]:
            column.append(symbol({w: v for (k, w), v in c.items() if k == kind}, a, b))
        assert column[-1] == potential * Poly(multiplier, a, b)
        components.append(column)
    for degree, expected in ((2, (2, 2)), (3, (2, 3))):
        rows = [
            [f.coeff_monomial(a**j * b ** (n - j)) for f in block]
            for block in components
            for n in range(degree + 1)
            for j in range(n + 1)
        ]
        m = Matrix(rows)
        assert (m[:, :2].rank(), m.rank()) == expected
    print("planar current order", 3, "cubic ranks", 2, 3)


def verify_unmarked_three():
    from block_homotopy import cochain

    a, b = symbols("a b")
    _, current, _ = stream()
    blocks = []
    for kind, weight, shear in (("dx", -1, 1 - b), ("dy", 1, a - 1)):
        columns = []
        for c in boundaries():
            columns.append(symbol({w: v for (k, w), v in c.items() if k == kind}, a, b))
        columns.append(Poly(shear, a, b))
        for name in ("A", "B"):
            values = {}
            for word in basis(5, weight)[3]:
                if word[-1][1] < 1:
                    continue
                inputs = [{(x, h + 2 - x): 1} for _, h, x in word]
                form = cochain(name, inputs)
                value = -form.get((0, 1), 0) if kind == "dx" else form.get((1, 0), 0)
                if value:
                    values[word] = value
            columns.append(symbol(values, a, b))
        columns.append(
            symbol({w: v for (k, w), v in current.items() if k == kind}, a, b)
        )
        blocks.append(columns)
    for degree, expected in ((2, (4, 4)), (3, (5, 6))):
        rows = Matrix(
            [
                [p.coeff_monomial(a**j * b ** (n - j)) for p in col]
                for col in blocks
                for n in range(degree + 1)
                for j in range(n + 1)
            ]
        )
        assert (rows[:, :-1].rank(), rows.rank()) == expected
    print("unmarked base and shear columns, cubic ranks", 5, 6)


def sign(p):
    return (-1) ** sum(a > b for i, a in enumerate(p) for b in p[i + 1 :])


def scalar_boundaries():
    words = basis(6, 0)[3]
    A = Matrix([[linear(w, "raise").get(k, 0) for k in words] for w in basis(6, -2)[3]])
    out = []
    for sol in A.nullspace():
        src = dict(zip(words, map(Fraction, sol)))
        out.append(
            {
                w: sum(src.get(k, 0) * v for k, v in diff(w).items())
                for w in basis(6, 0)[4]
            }
        )
    return out


def ordinary_wave_coefficient(w, plane):
    z = 0
    for p in permutations(plane):
        v = sign(tuple(plane.index(t) for t in p))
        for (_, h, i), (x, y) in zip(w[:3], p):
            v *= Fraction(x**i * y ** (h + 2 - i), factorial(i) * factorial(h + 2 - i))
        z += v
    return z


def verify_scalar_four():
    _, P, *_ = lift()
    cols = scalar_boundaries() + [P]
    rows = []
    for s, t in ((1, 1), (1, -1), (2, -1), (-2, -1), (1, 2)):
        plane = ((1, 0), (0, 1), (s, t))
        m = (-1 - s, -1 - t)
        pre = {w: ordinary_wave_coefficient(w, plane) for w in basis(6, 0)[4]}
        for n in range(3):
            for i in range(n + 1):
                j = n - i
                vals = []
                for col in cols:
                    vals.append(
                        sum(
                            c
                            * pre[w]
                            * Fraction(
                                m[0] ** (x - i) * m[1] ** (y - j),
                                factorial(i)
                                * factorial(j)
                                * factorial(x - i)
                                * factorial(y - j),
                            )
                            for w, c in col.items()
                            for h, x in [w[-1][1:]]
                            for y in [h + 2 - x]
                            if x >= i and y >= j
                        )
                    )
                rows.append((n, vals))
    for_degree = []
    for d in (1, 2):
        A = Matrix([r for n, r in rows if n <= d])
        for_degree.append((A[:, :-1].rank(), A.rank()))
    assert for_degree == [(4, 4), (5, 6)]
    print("four-scalar boundary and augmented ranks", for_degree)


if __name__ == "__main__":
    verify()
    verify_scalar_four()
    verify_unmarked_three()
