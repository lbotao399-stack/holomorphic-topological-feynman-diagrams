#!/usr/bin/env python3
"""Exact checks for the uniform Lie-algebra and first-order-cohomology proofs.

No repository imports.  Python >= 3.10 and SymPy suffice.

1. Generate each root system by integral simple reflections, determine its
   highest root and Dynkin support, and check the complete root count.
2. Independently evaluate the complete Racah--Speiser Weyl sum for A1--A4,
   B2--B4, C2--C4, D4, F4, G2 and E6 (51,840 Weyl elements in the last case).
3. Compute the exact symmetric polynomial cocycle ranks at bidegrees 0..4.
   Construct and verify every odd-monomial primitive in those slices.
4. Check the primitive's CE sign using all 27 ordered sl2 color triples.

The all-rank theorem is proved analytically in sections/lie_classification.tex;
the finite checks here are independent regression certificates.
"""

from __future__ import annotations

from collections import deque
from itertools import permutations, product
import json

import sympy as sp


def cartan(kind: str, rank: int) -> tuple[tuple[int, ...], ...]:
    """A[i,j] = <alpha_j, alpha_i^vee>."""
    a = [[2 * int(i == j) for j in range(rank)] for i in range(rank)]

    def edge(i: int, j: int, ij: int = -1, ji: int = -1) -> None:
        a[i][j], a[j][i] = ij, ji

    if kind in {"A", "B", "C"}:
        for i in range(rank - 1):
            edge(i, i + 1)
        if kind == "B":
            edge(rank - 2, rank - 1, -1, -2)
        if kind == "C":
            edge(rank - 2, rank - 1, -2, -1)
    elif kind == "D":
        for i in range(rank - 2):
            edge(i, i + 1)
        edge(rank - 3, rank - 1)
    elif kind == "E":
        edge(0, 2)
        edge(1, 3)
        for i in range(2, rank - 1):
            edge(i, i + 1)
    elif kind == "F":
        edge(0, 1)
        edge(1, 2, -1, -2)
        edge(2, 3)
    elif kind == "G":
        edge(0, 1, -3, -1)
    else:
        raise ValueError((kind, rank))
    return tuple(map(tuple, a))


def roots(a: tuple[tuple[int, ...], ...]) -> set[tuple[int, ...]]:
    rank = len(a)
    out = {tuple(int(i == j) for i in range(rank)) for j in range(rank)}
    pending = deque(out)
    while pending:
        root = pending.popleft()
        for i in range(rank):
            reflected = list(root)
            reflected[i] -= sum(a[i][j] * root[j] for j in range(rank))
            reflected = tuple(reflected)
            if reflected not in out:
                out.add(reflected)
                pending.append(reflected)
    return out


def labels(a, root):
    return tuple(sum(x * y for x, y in zip(row, root)) for row in a)


def weyl_sum(a, root_set, theta_labels):
    """Visit the whole regular Weyl orbit; do not use highest-root support."""
    rank = len(a)
    highest_regular = tuple(x + 1 for x in theta_labels)
    root_labels = {labels(a, root) for root in root_set}
    signs = {highest_regular: 1}
    pending = deque([highest_regular])
    multiplicity = 0
    contributors = []
    while pending:
        weight = pending.popleft()
        sign = signs[weight]
        difference = tuple(x - y for x, y in zip(highest_regular, weight))
        weight_multiplicity = rank if not any(difference) else int(
            difference in root_labels
        )
        multiplicity += sign * weight_multiplicity
        if weight_multiplicity:
            contributors.append((difference, sign, weight_multiplicity))
        for i in range(rank):
            reflected = tuple(weight[j] - weight[i] * a[j][i] for j in range(rank))
            if reflected in signs:
                assert signs[reflected] == -sign
            else:
                signs[reflected] = -sign
                pending.append(reflected)
    return multiplicity, len(signs), contributors


def check_roots():
    cases = (
        [("A", n) for n in range(1, 9)]
        + [(k, n) for k in ("B", "C") for n in range(2, 9)]
        + [("D", n) for n in range(4, 9)]
        + [("E", 6), ("E", 7), ("E", 8), ("F", 4), ("G", 2)]
    )
    expected_exceptional = {("E", 6): 72, ("E", 7): 126, ("E", 8): 240,
                            ("F", 4): 48, ("G", 2): 12}
    rows = []
    for kind, rank in cases:
        a = cartan(kind, rank)
        root_set = roots(a)
        if kind == "A":
            expected_roots = rank * (rank + 1)
        elif kind in {"B", "C"}:
            expected_roots = 2 * rank * rank
        elif kind == "D":
            expected_roots = 2 * rank * (rank - 1)
        else:
            expected_roots = expected_exceptional[kind, rank]
        assert len(root_set) == expected_roots, (kind, rank, len(root_set))
        positive = [v for v in root_set if min(v) >= 0]
        height = max(map(sum, positive))
        highest = [v for v in positive if sum(v) == height]
        assert len(highest) == 1
        theta_labels = labels(a, highest[0])
        assert min(theta_labels) >= 0
        support = sum(x != 0 for x in theta_labels)
        expected = 2 if kind == "A" and rank >= 2 else 1
        assert support == expected
        row = {"type": f"{kind}{rank}", "roots": len(root_set),
               "highest_root": highest[0], "dynkin_labels": theta_labels,
               "multiplicity": support}
        if rank <= 4 or (kind, rank) == ("E", 6):
            direct, order, contributors = weyl_sum(a, root_set, theta_labels)
            assert direct == support
            assert len(contributors) == 1 + rank - support
            assert all(sign == -1 and mult == 1
                       for delta, sign, mult in contributors if any(delta))
            row["full_weyl_order"] = order
            row["full_weyl_sum"] = direct
        rows.append(row)
    return rows


x, y, u, v, r, s, t = sp.symbols("x y u v r s t")


def substitute(poly, first, second):
    return sp.expand(poly.subs(dict(zip((x, y, u, v), (*first, *second))),
                              simultaneous=True))


def span(polynomials, variables):
    if not polynomials:
        return []
    monomials = sorted({m for p in polynomials
                        for m, coefficient in sp.Poly(p, *variables).terms()
                        if coefficient})
    columns = [sp.Matrix([sp.Poly(p, *variables).coeff_monomial(m)
                          for m in monomials]) for p in polynomials]
    matrix = sp.Matrix.hstack(*columns)
    _, pivots = matrix.rref()
    return [polynomials[i] for i in pivots]


def cocycle(poly):
    return sp.expand(
        substitute(poly, (x, y), (u, v))
        + substitute(poly, (x + u, y + v), (r, s))
        - substitute(poly, (u, v), (r, s))
        - substitute(poly, (x, y), (u + r, v + s))
    )


def symmetric_basis(a, b):
    momenta = [(x, y), (u, v), (-x - u, -y - v)]
    polynomials = []
    for i in range(a + 1):
        for j in range(b + 1):
            poly = x**i * u**(a - i) * y**j * v**(b - j)
            symmetrized = sp.expand(sum(
                substitute(poly, order[0], order[1])
                for order in permutations(momenta)
            ))
            if symmetrized:
                polynomials.append(symmetrized)
    return span(polynomials, (x, y, u, v))


def check_polynomial_cocycles():
    rows = []
    for a, b in product(range(5), repeat=2):
        basis = symmetric_basis(a, b)
        images = [cocycle(poly) for poly in basis]
        rank = len(span(images, (x, y, u, v, r, s)))
        cycles = len(basis) - rank
        expected_cycles = int(a + b == 0 or (a + b >= 3 and (a + b) % 2))
        assert cycles == expected_cycles, (a, b, cycles, expected_cycles)
        boundaries = int(a + b >= 3 and (a + b) % 2)
        if boundaries:
            primitive = x**a * y**b
            poly = sp.expand((x + u)**a * (y + v)**b - primitive - u**a * v**b)
            assert poly != 0
            assert cocycle(poly) == 0
            assert substitute(poly, (u, v), (x, y)) == poly
            assert substitute(poly, (x, y), (-x - u, -y - v)) == poly
            one_form_x = sp.diff(poly, u).subs({u: 0, v: 0})
            one_form_y = sp.diff(poly, v).subs({u: 0, v: 0})
            integrated = sp.integrate(
                x * one_form_x.subs({x: t*x, y: t*y}, simultaneous=True)
                + y * one_form_y.subs({x: t*x, y: t*y}, simultaneous=True),
                (t, 0, 1),
            )
            assert sp.expand(integrated - primitive) == 0
        assert cycles - boundaries == int(a + b == 0)
        rows.append({"bidegree": (a, b), "C3": len(basis), "rank_d3": rank,
                     "Z3": cycles, "B3": boundaries,
                     "H3": cycles - boundaries})
    # A symmetric quadratic invariant is a useful nonclosed control.
    assert cocycle(x*x + u*u + (-x-u)**2) != 0
    return rows


def check_ce_sign():
    e = sp.Matrix([[0, 1], [0, 0]])
    f = e.T
    h = sp.diag(1, -1)
    generators = (e, f, h)
    momenta = ((x, y), (u, v), (-x-u, -y-v))
    primitive = lambda p: p[0]**3 * p[1]**2
    polynomial = primitive((x+u, y+v)) - primitive((x,y)) - primitive((u,v))
    nonzero = 0
    for colors in product(generators, repeat=3):
        actual = 0
        for i in range(3):
            for j in range(i+1, 3):
                k = next(k for k in range(3) if k not in (i,j))
                bracket = colors[i] * colors[j] - colors[j] * colors[i]
                merged = tuple(momenta[i][a]+momenta[j][a] for a in range(2))
                actual += (-1)**(i+j) * (-primitive(merged)) * sp.trace(bracket*colors[k])
        expected = polynomial * sp.trace(colors[0]*(colors[1]*colors[2]-colors[2]*colors[1]))
        assert sp.expand(actual-expected) == 0
        nonzero += int(sp.expand(expected) != 0)
    assert nonzero == 6
    return {"ordered_color_triples": 27, "nonzero_color_triples": nonzero,
            "primitive_degree": 5, "all_CE_signs": "PASS"}


def main():
    output = {"root_systems": check_roots(),
              "polynomial_cocycles": check_polynomial_cocycles(),
              "CE_primitive": check_ce_sign()}
    output["status"] = "PASS"
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
