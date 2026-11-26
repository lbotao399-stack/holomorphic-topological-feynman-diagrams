#!/usr/bin/env python3
"""Exact all-order maximal joining test on products of protected constituents.

No earlier repository evaluator is imported. The all-loop fan coefficient is
an explicit inherited input. The new matrix, Berezin, multiplicity, and
order-isolation calculations are independent exact checks.
"""
from collections import defaultdict
from itertools import permutations
from math import factorial
from pathlib import Path
import json

import sympy as sp


HERE = Path(__file__).resolve().parent
PINNED_COMMIT = "29ef4be3c92b8055ac51067357524c1048e6cda0"
g, tau = sp.symbols("g tau", nonzero=True)


def wedge(left, right):
    if left & right:
        return None, 0
    inversions = sum((right & ((1 << i) - 1)).bit_count()
                     for i in range(6) if left >> i & 1)
    return left | right, (-1)**inversions


def left_derivative(mask, index):
    if not mask >> index & 1:
        return None, 0
    return mask ^ (1 << index), (-1)**((mask & ((1 << index)-1)).bit_count())


def apply_component_derivative(term, component):
    """delta/dC acting on theta_1..3, beta_1..3 and one commuting b."""
    coefficient, mask, b_power = term
    if component == "b":
        if not b_power:
            return sp.S.Zero, 0, 0
        return sp.expand(-sp.I*g*b_power*coefficient), mask, b_power-1
    i = int(component)
    remaining, sign = left_derivative(mask, 3+i)
    if not sign:
        return sp.S.Zero, 0, 0
    result, insertion_sign = wedge(1 << i, remaining)
    if not insertion_sign:
        return sp.S.Zero, 0, 0
    return sp.expand(-sign*insertion_sign*coefficient), result, b_power


def berezin_assignments():
    records = []
    total_orders = 0
    for n in range(2, 15):
        m = n+1
        expected_one = -(sp.I*g*(-1))**(n-2)*factorial(n-2)
        summed = 0
        count = 0
        # The three beta slots are distinct. Repeated b derivatives are not
        # artificially relabeled; differentiating b^(n-2) supplies (n-2)!.
        for positions in permutations(range(m), 3):
            sequence = ["b"]*m
            for i, position in enumerate(positions):
                sequence[position] = str(i)
            term = (sp.S.One, 0b111000, n-2)
            for component in reversed(sequence):
                term = apply_component_derivative(term, component)
            coefficient, mask, power = term
            assert mask == 0b000111 and power == 0
            assert sp.expand(coefficient-expected_one) == 0
            summed += coefficient
            count += 1
        expected_sum = -(-sp.I*g)**(n-2)*factorial(n+1)
        assert sp.expand(summed-expected_sum) == 0
        records.append({"n": n, "distinct_beta_b_operator_orders": count,
                        "labeled_assignment_factor": factorial(n+1),
                        "berezin_sum": str(sp.factor(summed))})
        total_orders += count
    return records, total_orders


PARITY = {"c": 1, "beta1": 1, "beta2": 1, "beta3": 1,
          "gamma1": 0, "gamma2": 0, "gamma3": 0, "b": 0}


def canonical_trace(word):
    candidates = []
    for split in range(len(word)):
        left, right = word[:split], word[split:]
        sign = (-1)**(sum(PARITY[x] for x in left)*sum(PARITY[x] for x in right))
        candidates.append((right+left, sign))
    canonical = min(word for word, sign in candidates)
    signs = {sign for word, sign in candidates if word == canonical}
    return canonical, (next(iter(signs)) if len(signs) == 1 else 0)


def theta_derivative(poly, index):
    result = defaultdict(lambda: sp.S.Zero)
    for (mask, word), coefficient in poly.items():
        remaining, sign = left_derivative(mask, index)
        if sign:
            result[remaining, word] += sign*coefficient
    return dict(result)


def trace_product(left, right):
    result = defaultdict(lambda: sp.S.Zero)
    for (a, u), x in left.items():
        for (b, v), y in right.items():
            mask, sign = wedge(a, b)
            if not sign:
                continue
            sign *= (-1)**(sum(PARITY[w] for w in u)*b.bit_count())
            word, trace_sign = canonical_trace(u+v)
            if trace_sign:
                result[mask, word] += x*y*sign*trace_sign
    return {key: sp.simplify(value) for key, value in result.items() if value != 0}


def graviton_components():
    superfield = {
        (0, ("c",)): sp.I*g,
        (1, ("gamma1",)): 1, (2, ("gamma2",)): 1, (4, ("gamma3",)): 1,
        (3, ("beta3",)): -1, (5, ("beta2",)): 1, (6, ("beta1",)): -1,
        (7, ("b",)): 1/(sp.I*g),
    }
    squared = []
    for i in range(3):
        derivative = theta_derivative(superfield, i)
        squared.append(trace_product(derivative, derivative))
    checks = []
    for derivative_axis, mask, word in (
        (1, 4, ("beta1", "gamma2")),
        (2, 1, ("beta2", "gamma3")),
        (0, 2, ("beta3", "gamma1")),
    ):
        component = {w: sp.simplify(-value/2)
                     for (a, w), value in squared[derivative_axis].items() if a == mask}
        assert component == {word: 1}
        checks.append({"/dtheta": derivative_axis+1, "theta_mask": mask,
                       "normalized_trace_component": str(component)})
    even = {word: sp.simplify(sp.I*g*value/2)
            for (mask, word), value in squared[0].items() if mask == 6}
    assert even == {("b", "gamma1"): 1, ("beta2", "beta3"): -sp.I*g}
    return {"odd": checks, "even": {",".join(k): str(v) for k, v in even.items()}}


def epsilon(i, j, k):
    i, j, k = int(i), int(j), int(k)
    if len({i,j,k}) < 3:
        return 0
    return (-1)**sum(a > b for a,b in ((i,j),(i,k),(j,k)))


def matrix_color_checks():
    epsilon_matrix = sp.Matrix(3, 3, lambda d,e: epsilon(2,d,e))
    root_first = sp.Matrix([1, sp.I, 0])
    root_second = sp.Matrix([1, -sp.I, 0])
    t3 = sp.diag(sp.Rational(1,2), -sp.Rational(1,2))
    D = sp.Matrix([[0,1],[0,0]])
    E = D.T
    commutated = D
    rows = []
    for m in range(1, 25):
        commutated = t3*commutated-commutated*t3
        direct_trace = sp.trace(commutated*E)
        epsilon_value = (root_first.T*epsilon_matrix**m*root_second)[0]
        assert direct_trace == 1
        assert sp.simplify(epsilon_value - 2*sp.I**(-m)*direct_trace) == 0
        reverse = E
        for _ in range(m):
            reverse = t3*reverse-reverse*t3
        assert sp.trace(reverse*D) == (-1)**m
        rows.append({"spokes": m, "root_pair_trace": str(direct_trace),
                     "direct_epsilon_chain": str(epsilon_value),
                     "reversed_endpoint_trace": str(sp.trace(reverse*D))})
    block_tests = 0
    for rank in range(2, 8):
        for gap in (-3,-1,1,2,4):
            X = sp.diag(sp.Rational(gap,2),-sp.Rational(gap,2),*([0]*(rank-2)))
            D = sp.zeros(rank); D[0,1] = 1
            E = D.T
            nested = D
            for m in range(1, 19):
                nested = X*nested-nested*X
                binomial = sum(
                    ((-1)**j*sp.binomial(m,j)*X**(m-j)*D*X**j
                     for j in range(m+1)), sp.zeros(rank))
                assert nested == binomial
                assert sp.trace(nested*E) == gap**m
                block_tests += 1
    return rows, block_tests


def loop_isolation():
    configurations = 0
    for n in range(2, 41):
        surviving = []
        for selected_beta_pairs in range(n-1):
            scalar_letters = n+1-selected_beta_pairs
            for scalar_derivatives in range(scalar_letters+1):
                # Every beta and b must be removed at the measurement locus.
                spokes = n+1+selected_beta_pairs+scalar_derivatives
                loop = spokes-1
                theta_degree = 3+2*selected_beta_pairs+2*scalar_derivatives
                assert theta_degree == 3+2*(loop-n)
                if theta_degree == 3:
                    surviving.append((loop,selected_beta_pairs,scalar_derivatives))
                configurations += 1
        assert surviving == [(n,0,0)]
    return {"families_checked": 39, "monomial_and_derivative_counts": configurations,
            "only_surviving_configuration": "loop=n; no beta-pair term; no scalar derivative"}


def exact_coefficients():
    rows = []
    for n in range(2, 25):
        fan = (-1)**n*sp.I**(n+1)/(2*factorial(n)**3*(n+1))
        berezin = -(-sp.I*g)**(n-2)*factorial(n+1)
        endpoint_ghosts = (sp.I*g)**2
        color_conversion = tau**n*sp.I**(-n-1)
        ordered = sp.simplify(fan*berezin*endpoint_ghosts*color_conversion)
        expected = -(sp.I*g*tau)**n/(2*factorial(n)**2)
        assert sp.simplify(ordered-expected) == 0
        # Moving the two odd amplitudes past each other, reversing the dotted
        # contraction, and reversing the adjoint chain multiply to +1.
        exchange = -(-1)**n*(-1)**(n+1)
        assert exchange == 1
        polarized = sp.simplify(2*ordered)
        assert sp.simplify(polarized+(sp.I*g*tau)**n/factorial(n)**2) == 0
        at_standard_trace = sp.simplify(polarized.subs(tau,sp.Rational(1,2)))
        rows.append({"n": n,
                     "coefficient_after_removing_(32pi2)^-n": str(polarized),
                     "tau_one_half": str(at_standard_trace)})
    # The known independent Pauli Q2 formula gives the same root-wave value.
    pauli = sp.I*g**2*(-2*sp.I)/sp.Integer(32768)
    new = -(sp.I*g*sp.Rational(1,2))**2/(factorial(2)**2*32**2)
    assert sp.simplify(pauli-new) == 0
    # Separate protected two-trace Q1 test. Derive its two Berezin orders
    # from gamma_1^2 beta_1, retaining the derivative factor of two.
    one_loop_berezin = 0
    for sequence in (("gamma", "beta"), ("beta", "gamma")):
        coefficient, mask, scalar_power = sp.S.One, 1 << 3, 2
        for component in reversed(sequence):
            if component == "gamma":
                remaining, sign = wedge(0b110, mask)
                coefficient *= -scalar_power*sign
                mask, scalar_power = remaining, scalar_power-1
            else:
                remaining, sign = left_derivative(mask, 3)
                mask, insertion_sign = wedge(1, remaining)
                coefficient *= -sign*insertion_sign
        assert mask == 0b111 and scalar_power == 1 and coefficient == 2
        one_loop_berezin += coefficient
    one_loop_fan = (-1)*sp.I**2/(2*factorial(1)**3*2)
    one_loop_ordered = sp.simplify(
        one_loop_fan*one_loop_berezin*(sp.I*g)**2*tau*sp.I**(-2))
    one_loop = sp.simplify(2*one_loop_ordered)
    assert one_loop == 2*tau*g**2
    assert one_loop.subs(tau,sp.Rational(1,2)) == g**2
    return {"n_ge_two": rows,
            "one_loop_coefficient_after_removing_(32pi2)^-1": str(one_loop),
            "two_loop_pauli_control_after_removing_pi^-4": str(pauli)}


def constant_cartan_control():
    v = sp.symbols("v")
    values = []
    for n in range(2, 25):
        poly = {0b111: sp.Rational(1, 8)}
        even = {0: v/2, 0b110: -sp.I*g/2}
        for _ in range(n-2):
            result = defaultdict(lambda: sp.S.Zero)
            for a, x in poly.items():
                for b, y in even.items():
                    mask, sign = wedge(a, b)
                    if sign:
                        result[mask] += sign*x*y
            poly = {mask: sp.expand(x) for mask, x in result.items() if x}
        assert poly == {0b111: v**(n-2)/2**(n+1)}
        values.append({"n": n, "coefficient_of_chi1_chi2_chi3": str(poly[0b111])})
    return {"one_loop_pair_coefficient_of_chi1": "1/4", "higher_family": values}


def main():
    berezin, assignment_count = berezin_assignments()
    colors, matrix_count = matrix_color_checks()
    result = {
        "status": "PASS",
        "arithmetic": "exact integers, rationals and symbolic powers of i,g,tau",
        "source_revision": PINNED_COMMIT,
        "source_git_blob_sha1": {
            "sources/LLC_main.tex": "db6b3994c61a469268e57b71bf426cc1281b4a10",
            "explore/twisted-holography-N4/NOTE.md": "439f5c157fb587742c954131f4a6beaa0e784de6",
        },
        "inherited_formula": "sources/LLC_main.tex equations eq-Qnzero and eq-A-16",
        "graviton_components": graviton_components(),
        "berezin_assignments": berezin,
        "distinct_berezin_operator_orders_checked": assignment_count,
        "direct_structure_constant_chain": colors,
        "additional_exact_matrix_checks": matrix_count,
        "loop_isolation": loop_isolation(),
        "coefficients": exact_coefficients(),
        "constant_cartan_nonexactness_controls": constant_cartan_control(),
        "all_order_proof": "REPORT.md",
        "scope": "Boundary cochain differential and bulk tree targets, not a computed bulk transfer.",
    }
    path = HERE/"all_order_join.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status": "PASS",
                      "distinct_berezin_operator_orders_checked": assignment_count,
                      "matrix_checks": matrix_count+len(colors),
                      "isolated_loop_families": 39},sort_keys=True))


if __name__ == "__main__":
    main()
