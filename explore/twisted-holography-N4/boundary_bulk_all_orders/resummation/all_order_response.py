#!/usr/bin/env python3
"""Independent exact arithmetic for the all-order protected-product detector.

No imports from the boundary or bulk certificates. Gaussian rational numbers
are pairs of fractions, so no floating-point complex arithmetic enters.
The analytic proof, not the finite table, establishes the all-order statement.
"""
from fractions import Fraction
from math import factorial
from pathlib import Path
import json


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def multiply(a, b):
    return a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0]


def scale(a, x):
    return x*a[0], x*a[1]


def i_power(n):
    return ((1, 0), (0, 1), (-1, 0), (0, -1))[n % 4]


def matrix_multiply(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def commutator(a, b):
    ab, ba = matrix_multiply(a, b), matrix_multiply(b, a)
    return tuple(tuple(ab[i][j]-ba[i][j] for j in range(2)) for i in range(2))


def trace_product(a, b):
    product = matrix_multiply(a, b)
    return product[0][0] + product[1][1]


def contracted_chain(n, gap):
    """Complete two-root color factor, including exchanged odd output waves."""
    diagonal = ((Fraction(gap, 2), 0), (0, -Fraction(gap, 2)))
    plus, minus = ((0, 1), (0, 0)), ((0, 0), (1, 0))
    first, second = plus, minus
    for _ in range(n + 1):
        first, second = commutator(diagonal, first), commutator(diagonal, second)
    # Reordering the odd probes costs -1, and exchanging the n spinor
    # contractions contributes (-1)^n.
    return (trace_product(first, minus)
            + (-1)**(n+1)*trace_product(second, plus))


def source_prefactor(n, tau):
    """Raw fan phase x Berezin phase x trace conversion x spoke assignments.

    Powers g^n and pi^(-2n) are removed. The factor 2 from the ordered ghost
    probes is supplied by contracted_chain, not hidden in this prefactor.
    """
    fan_phase = scale(i_power(n+1), (-1)**n)
    berezin_phase = scale(multiply(i_power(-(n-2)), i_power(2)), -1)
    trace_phase = i_power(-(n+1))
    total = multiply(multiply(fan_phase, berezin_phase), trace_phase)
    rational = (tau**n * factorial(n+1)
                / (2 * Fraction(32)**n * factorial(n)**3 * (n+1)))
    return scale(total, rational)


def serialized(value):
    return {"real": str(value[0]), "imaginary": str(value[1])}


def main():
    rows = []
    checks = 0
    tau = Fraction(1, 2)
    for n in range(2, 25):
        actual = scale(source_prefactor(n, tau), contracted_chain(n, 1))
        expected = scale(i_power(n), -Fraction(1, 64**n*factorial(n)**2))
        assert actual == expected
        rows.append({"boundary_loop_order": n,
                     "coefficient_after_removing_g_power_and_pi_power": serialized(actual),
                     "bulk_touched_input_traces": n+1,
                     "bulk_output_traces": 1,
                     "ribbon_genus": 0,
                     "internal_color_faces": 0,
                     "fixed_tHooft_N_power": -n})
        checks += 1
        for gap in (1, 2, 3):
            assert contracted_chain(n, gap) == 2*gap**(n+1)
            for area in (-2, -1, 1, 3):
                coefficient = scale(source_prefactor(n, tau),
                                    contracted_chain(n, gap)*area**n)
                formula = scale(expected, gap**(n+1)*area**n)
                assert coefficient == formula
                checks += 1
    # Build the entire power series by the differential equation z f''+f'=f.
    # Its initial value is 1. The recurrence determines every coefficient.
    series = [Fraction(1)]
    for n in range(1, 65):
        series.append(series[-1]/(n*n))
    assert series == [Fraction(1, factorial(n)**2) for n in range(65)]
    for n in range(64):
        assert (n+1)**2*series[n+1] == series[n]
    # The subtractions at orders zero and one make division by u^2 formal.
    tail = [Fraction(0), Fraction(0), *series[2:]]
    assert tail[:2] == [0, 0]
    for n in range(2, 25):
        actual = scale(i_power(n), -tail[n]/64**n)
        recorded = rows[n-2]["coefficient_after_removing_g_power_and_pi_power"]
        assert serialized(actual) == recorded
    result = {
        "status": "PASS",
        "arithmetic": "integers and Gaussian rational numbers",
        "source_revision": "29ef4be3c92b8055ac51067357524c1048e6cda0",
        "inherited_input": "sources/LLC_main.tex, eq-Qnzero and eq-A-16",
        "independent_phase_color_and_area_checks": checks,
        "all_order_series_recurrence_checked_through": 64,
        "one_loop_protected_pair_after_removing_g_squared_over_pi_squared": "1/32",
        "rows": rows,
        "scope": "Boundary detector and required composite-map Ward sum; not a computed bulk amplitude.",
    }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status": "PASS", "phase_color_and_area_checks": checks,
                      "series_recurrence_orders": 64}))


if __name__ == "__main__":
    main()
