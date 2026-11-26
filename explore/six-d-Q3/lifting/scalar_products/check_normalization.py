"""Two-loop Ward test at general momenta in two independent sl(3) factors."""

from fractions import Fraction as F
from itertools import permutations
from random import Random
import sympy as sp
from kernel import kernel, plus, crossed


def bracket(a, b):
    return a * b - b * a


def chain(x, ghosts, y):
    for c in reversed(ghosts):
        y = bracket(c, y)
    return int(sp.trace(x * y))


def sign(mask, added):
    return (-1) ** sum(
        (added & ((1 << i) - 1)).bit_count()
        for i in range(mask.bit_length())
        if mask >> i & 1
    )


def merge(parts):
    total = 0
    factor = 1
    for mask in parts:
        if mask & total:
            return 0
        factor *= sign(total, mask)
        total |= mask
    return factor


def check_probe(momentum):
    # Grassmann bit zero is the infinitesimal Q0 parameter. Other bits label ghost waves.
    x = [sp.diag(1, -1, 0), sp.diag(0, 1, -1)]
    y = [sp.diag(1, 1, -2), sp.diag(1, -1, 0)]
    ghost = []
    for i, j in ((0, 1), (1, 2), (2, 0), (0, 1), (1, 0)):
        c = sp.zeros(3)
        c[i, j] = 1
        ghost.append(c)
    mask = [1 << (j + 1) for j in range(5)]
    p, px, py = momentum[:5], momentum[5:7], momentum[7:9]
    all_terms = []
    for sector, idx in enumerate(((0, 1, 2), (3, 4))):
        Cs = [(mask[i], p[i], ghost[i]) for i in idx]
        Cs.extend(
            (1 | mask[i] | mask[j], plus(p[i], p[j]), bracket(ghost[i], ghost[j]))
            for i in idx
            for j in idx
            if i < j
        )
        X = [(0, px[sector], x[sector])] + [
            (1 | mask[i], plus(px[sector], p[i]), bracket(ghost[i], x[sector]))
            for i in idx
        ]
        Y = [(0, py[sector], y[sector])] + [
            (1 | mask[i], plus(py[sector], p[i]), bracket(ghost[i], y[sector]))
            for i in idx
        ]
        terms = []
        for mx, ax, colx in X:
            for m1, p1, c1 in Cs:
                for m2, p2, c2 in Cs:
                    for my, _, coly in Y:
                        factor = merge([mx, m1, m2, my])
                        if not factor:
                            continue
                        color = chain(colx, (c1, c2), coly)
                        if color:
                            terms.append(
                                (mx | m1 | m2 | my, factor * color, ax, (p1, p2))
                            )
        all_terms.append(terms)
    q0 = F()
    for mask0, v0, ax, pair0 in all_terms[0]:
        for mask1, v1, bx, pair1 in all_terms[1]:
            if mask0 | mask1 != 63:
                continue
            q0 += merge([mask0, mask1]) * v0 * v1 * kernel(plus(ax, bx), [pair0, pair1])
    box = F()
    for tri in permutations((0, 1, 2)):
        darr = crossed(p[tri[0]], p[tri[1]])
        det = sum(i * j for i, j in zip(darr, p[tri[2]]))
        y0 = tuple(
            sum(F(3 - j, 4) * p[i][k] for j, i in enumerate(tri)) for k in range(3)
        )
        y1 = tuple(
            sum(F(j + 1, 4) * p[i][k] for j, i in enumerate(tri)) for k in range(3)
        )
        for pair in permutations((3, 4)):
            s = merge([mask[j] for j in (*tri, *pair)])
            color2 = chain(x[1], [ghost[j] for j in pair], y[1])
            for u, left, right, vertex in ((1, x[0], y[0], y0), (-1, y[0], x[0], y1)):
                color1 = chain(left, [ghost[j] for j in tri], right)
                position = plus(plus(px[0], px[1]), vertex)
                box += F(u * s * color1 * color2 * det, 6) * kernel(
                    position, [(p[pair[0]], p[pair[1]])]
                )
    assert 2 * q0 == box, (q0, box)
    return q0


def main():
    rng = Random(182)
    probes = [
        check_probe([tuple(rng.randrange(-3, 4) for j in range(3)) for i in range(9)])
        for _ in range(5)
    ]
    assert any(probes)
    print("unrestricted two-loop Ward residues", ", ".join(str(x) for x in probes))
    for m in range(1, 50):
        triangle = (F(1, 2), F(1, m + 1), F(1, m * (m + 1)))
        simplex = (
            F(1, 2 * (m + 2)),
            F(1, m + 1) - F(1, m + 2),
            (F(1, m) - F(2, m + 1) + F(1, m + 2)) / 2,
        )
        assert tuple((m + 2) * s for s in simplex) == triangle


if __name__ == "__main__":
    main()
