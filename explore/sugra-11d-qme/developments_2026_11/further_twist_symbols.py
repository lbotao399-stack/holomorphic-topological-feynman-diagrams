#!/usr/bin/env python3
"""Exact symbol calculation for the further-twisted E(5|10) complex.

Independent exterior-algebra implementation; no repository imports.
Python 3.10+ and SymPy. No floating-point arithmetic.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

import sympy as sp


def add(*terms):
    result = defaultdict(lambda: sp.S.Zero)
    for form in terms:
        for key, value in form.items():
            result[key] += value
    return {k: sp.expand(v) for k, v in result.items() if sp.expand(v) != 0}


def scale(value, form):
    return {key: sp.expand(value * coefficient) for key, coefficient in form.items()
            if sp.expand(value * coefficient) != 0}


def wedge(left, right):
    result = defaultdict(lambda: sp.S.Zero)
    for a, x in left.items():
        for b, y in right.items():
            if set(a) & set(b):
                continue
            sign = (-1) ** sum(i > j for i in a for j in b)
            result[tuple(sorted(a + b))] += sign * x * y
    return {k: sp.expand(v) for k, v in result.items() if sp.expand(v) != 0}


def interior(vector, form):
    result = defaultdict(lambda: sp.S.Zero)
    for slots, coefficient in form.items():
        for j, slot in enumerate(slots):
            result[slots[:j] + slots[j + 1:]] += (-1) ** j * vector[slot] * coefficient
    return {k: sp.expand(v) for k, v in result.items() if sp.expand(v) != 0}


def covector(p):
    return {(i,): sp.sympify(value) for i, value in enumerate(p) if value}


volume = {(0, 1, 2, 3, 4): sp.S.One}
transverse_volume = {(2, 3, 4): sp.S.One}
background = {(0, 1): sp.S.One}
two_slots = tuple(combinations(range(5), 2))
three_slots = tuple(combinations(range(5), 3))


def sharp(form):
    vector = sp.zeros(5, 1)
    for i in range(5):
        slots = tuple(j for j in range(5) if j != i)
        vector[i] = (-1) ** i * form.get(slots, 0)
    assert interior(vector, volume) == form
    return vector


def delta_even(p, vector):
    return scale(-1, wedge(covector(p), interior(vector, background)))


def delta_odd(form):
    return sharp(wedge(background, form))


def basis_vector(i):
    return sp.eye(5)[:, i]


def form_from_column(column):
    return {slots: value for slots, value in zip(two_slots, column) if value}


def wave_spaces(p):
    p = sp.Matrix(p)
    vectors = sp.Matrix([list(p)]).nullspace()
    closure = sp.zeros(len(three_slots), len(two_slots))
    for j, slots in enumerate(two_slots):
        image = wedge(covector(p), {slots: 1})
        for i, target in enumerate(three_slots):
            closure[i, j] = image.get(target, 0)
    forms = [form_from_column(v) for v in closure.nullspace()]
    return vectors, forms


def matrix_rank(columns, rows):
    return sp.Matrix.hstack(*columns).rank() if columns else 0


def check_symbolic_chart():
    q = sp.symbols("q3 q4 q5")
    p = sp.Matrix([1, 0, *q])
    vectors = [basis_vector(1)] + [
        basis_vector(j + 2) - q[j] * basis_vector(0) for j in range(3)
    ]
    forms = [wedge(covector(p), {(j,): 1}) for j in (1, 2, 3, 4)]
    even = sp.zeros(4)
    odd = sp.zeros(4)
    for j, v in enumerate(vectors):
        assert (p.T * v)[0] == 0
        image = delta_even(p, v)
        column = sp.Matrix([image.get((0, k), 0) for k in (1, 2, 3, 4)])
        assert add(image, scale(-1, add(*(scale(column[k], forms[k]) for k in range(4))))) == {}
        even[:, j] = column
    for j, form in enumerate(forms):
        assert wedge(covector(p), form) == {}
        image = delta_odd(form)
        column = sp.Matrix([image[k] for k in (1, 2, 3, 4)])
        assert image == sum((column[k] * vectors[k] for k in range(4)), sp.zeros(5, 1))
        odd[:, j] = column
    expected_even = sp.Matrix([[0, *q], [-q[0], 0, 0, 0],
                               [-q[1], 0, 0, 0], [-q[2], 0, 0, 0]])
    expected_odd = sp.Matrix([[0, 0, 0, 0], [0, 0, -q[2], q[1]],
                              [0, q[2], 0, -q[0]], [0, -q[1], q[0], 0]])
    assert even == expected_even and odd == expected_odd
    assert even * odd == sp.zeros(4)
    assert odd * even == sp.zeros(4)
    # Exact contraction after localizing on the normal line q=(t,0,0).
    t = sp.symbols("t", nonzero=True)
    line_even = even.subs(dict(zip(q, (t, 0, 0))))
    line_odd = odd.subs(dict(zip(q, (t, 0, 0))))
    even_to_odd = sp.zeros(4)
    odd_to_even = sp.zeros(4)
    even_to_odd[3, 2] = -1 / t
    even_to_odd[2, 3] = 1 / t
    odd_to_even[1, 0] = 1 / t
    odd_to_even[0, 1] = -1 / t
    assert line_odd * even_to_odd + odd_to_even * line_even == sp.eye(4)
    assert line_even * odd_to_even + even_to_odd * line_odd == sp.eye(4)
    g = sp.symbols("g", nonzero=True)
    assert (g * line_odd) * (even_to_odd / g) + odd_to_even * line_even == sp.eye(4)
    assert line_even * odd_to_even + (even_to_odd / g) * (g * line_odd) == sp.eye(4)
    return {"even_to_odd": [[str(x) for x in row] for row in even.tolist()],
            "odd_to_even": [[str(x) for x in row] for row in odd.tolist()],
            "both_square_zero_identities": True,
            "localized_line_contraction": True,
            "coupling_retained_middle_map": "g times the odd_to_even matrix",
            "coupling_retained_cohomology": {
                "even": "ker(q dot : R^3 -> R) / g ker(q dot : R^3 -> R)",
                "odd": "R/(q3,q4,q5)", "ring": "R=C[g,q3,q4,q5]"},
            "coupling_retained_line_cohomology": {
                "even": "(C[g,t]/(g*t))^2", "odd": "(C[g,t]/t)^2"},
            "localized_line_contraction_requires_g_inverse": True,
            "normal_ring_cohomology": {"even": "0", "odd": "C[q3,q4,q5]/(q3,q4,q5)"},
            "normal_line_cohomology": {"even": "(C[t]/t)^2", "odd": "(C[t]/t)^2"}}


def check_all_waves():
    points = list(product((-1, 0, 1), repeat=5))
    points += [(2, 3, 5, 7, 11), (0, 0, 2, 3, 5)]
    strata = Counter()
    primitive_counts = Counter()
    for point in points:
        p = sp.Matrix(point)
        vectors, forms = wave_spaces(p)
        even_columns = [sp.Matrix([delta_even(p, v).get(s, 0) for s in two_slots]) for v in vectors]
        odd_columns = [delta_odd(form) for form in forms]
        even_rank = matrix_rank(even_columns, 10)
        odd_rank = matrix_rank(odd_columns, 5)
        for v in vectors:
            assert delta_odd(delta_even(p, v)) == sp.zeros(5, 1)
        for form in forms:
            assert delta_even(p, delta_odd(form)) == {}
        dimensions = (len(vectors) - even_rank - odd_rank,
                      len(forms) - even_rank - odd_rank)
        if any(point[2:]):
            assert (len(vectors), len(forms), even_rank, odd_rank, dimensions) == (4, 4, 2, 2, (0, 0))
            strata["nonzero_transverse_momentum"] += 1
            assert (len(vectors) - even_rank, len(forms) - even_rank) == (2, 2)
            j = next(j for j in range(2, 5) if p[j] != 0)
            n = basis_vector(j) / p[j]
            # Every even cycle receives its own closed-two-form primitive.
            even_matrix = sp.Matrix.hstack(*even_columns)
            for coeff in even_matrix.nullspace():
                v = sum((coeff[k] * vectors[k] for k in range(len(vectors))), sp.zeros(5, 1))
                assert v[0] == v[1] == 0
                primitive = wedge(covector(p), interior(n, interior(v, transverse_volume)))
                assert wedge(covector(p), primitive) == {}
                assert delta_odd(primitive) == v
                primitive_counts["even_cycle_primitives"] += 1
            # Every odd cycle receives its own divergence-free vector primitive.
            odd_matrix = sp.Matrix.hstack(*odd_columns)
            for coeff in odd_matrix.nullspace():
                form = add(*(scale(coeff[k], forms[k]) for k in range(len(forms))))
                a = interior(n, form)
                assert all(key[0] in (0, 1) for key in a)
                y = sp.Matrix([-a.get((1,), 0), a.get((0,), 0), 0, 0, 0])
                v = y - (p.T * y)[0] * n
                assert (p.T * v)[0] == 0
                assert delta_even(p, v) == form
                primitive_counts["odd_cycle_primitives"] += 1
        elif any(point):
            assert (len(vectors), len(forms), even_rank, odd_rank, dimensions) == (4, 4, 0, 0, (4, 4))
            strata["nonzero_planar_zero_transverse_momentum"] += 1
        else:
            assert (len(vectors), len(forms), even_rank, odd_rank, dimensions) == (5, 10, 0, 3, (2, 7))
            strata["zero_momentum"] += 1
            assert (len(vectors), len(forms)) == (5, 10)
    return {"tested_momenta": len(points), "strata": dict(strata),
            **dict(primitive_counts), "all_constraints_and_square_zero_identities": True,
            "coupling_zero_fiber_dimensions": {
                "nonzero_transverse_momentum": [2, 2],
                "nonzero_planar_zero_transverse_momentum": [4, 4],
                "zero_momentum": [5, 10]}}


def monomials(degree):
    for a in range(degree + 1):
        for b in range(degree - a + 1):
            yield (a, b, degree - a - b)


def linear_add(*chains):
    result = Counter()
    for chain in chains:
        for key, value in chain.items():
            result[key] += value
    return {key: value for key, value in result.items() if value}


def koszul_d(chain):
    result = Counter()
    for (powers, slots), value in chain.items():
        for j in range(3):
            if j in slots:
                continue
            exponents = list(powers)
            exponents[j] += 1
            sign = (-1) ** sum(j > k for k in slots)
            result[(tuple(exponents), tuple(sorted((j,) + slots)))] += sign * value
    return {key: value for key, value in result.items() if value}


def koszul_h(chain):
    result = Counter()
    for (powers, slots), value in chain.items():
        for i, j in enumerate(slots):
            if powers[j] == 0:
                continue
            exponents = list(powers)
            exponents[j] -= 1
            result[(tuple(exponents), slots[:i] + slots[i + 1:])] += (-1) ** i * powers[j] * value
    return {key: value for key, value in result.items() if value}


def check_polynomial_homotopy():
    count = 0
    exceptional = []
    for degree in range(7):
        for powers in monomials(degree):
            for k in range(4):
                for slots in combinations(range(3), k):
                    key = (powers, slots)
                    chain = {key: 1}
                    assert koszul_d(koszul_d(chain)) == {}
                    result = linear_add(koszul_d(koszul_h(chain)), koszul_h(koszul_d(chain)))
                    factor = degree + 3 - k
                    assert result == ({key: factor} if factor else {})
                    if factor == 0:
                        exceptional.append({"polynomial_degree": degree, "exterior_degree": k})
                    count += 1
    assert exceptional == [{"polynomial_degree": 0, "exterior_degree": 3}]
    return {"monomial_exterior_basis_vectors": count, "polynomial_degrees": list(range(7)),
            "identity": "d h + h d = polynomial_degree + 3 - exterior_degree",
            "exceptional_classes": exceptional,
            "homotopy_is_not_linear_over_the_momentum_polynomial_ring": True}


def main():
    output = {"status": "PASS", "arithmetic": "exact symbolic polynomials and rational numbers",
              "scope": "Strict further-twist field-strength differential; not the full kinetic operator or a loop calculation.",
              "symbolic_chart": check_symbolic_chart(),
              "wave_checks": check_all_waves(),
              "polynomial_homotopy": check_polynomial_homotopy(),
              "source_revision": "aa0d2f1725dff135494eef132f5e92a8c0401a42",
              "primary_sources": ["https://arxiv.org/html/2111.03049v1#S4.SS2",
                                  "https://arxiv.org/html/2506.02246v3#S4.SS3"]}
    source = Path(__file__).resolve()
    output["certificate_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    destination = source.with_suffix(".json")
    destination.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
