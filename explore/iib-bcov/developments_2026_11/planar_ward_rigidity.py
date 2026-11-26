#!/usr/bin/env python3
"""Exact, independent algebraic certificate for the planar bulk Ward tests.

This program imports no other repository module. It does not integrate a new
worldsheet or heat-kernel graph. Its only inherited analytic input is the two
previously computed reduced binary coefficients, recorded with pinned sources
in the output. All new claims checked here concern complete cyclic compositions
of those binaries and the cancellation of an unknown higher-genus coefficient.

Run with Python 3 and SymPy:
    python explore/iib-bcov/developments_2026_11/planar_ward_rigidity.py
    python explore/iib-bcov/developments_2026_11/planar_ward_rigidity.py --check

No floating-point arithmetic is used. The report contains the all-order proof;
finite checks below independently check its polynomial and differential data.
"""

from argparse import ArgumentParser
from fractions import Fraction
from functools import cache
from math import comb, factorial, gcd
from pathlib import Path
import json

import sympy as sp


PINNED_COMMIT = "aa0d2f1725dff135494eef132f5e92a8c0401a42"
REPOSITORY = "lbotao399-stack/holomorphic-topological-feynman-diagrams"
SOURCE_BLOBS = {
    "explore/iib-bcov/README.md": "20565aba20bcb73bc80ad71442bc4c649ca5a182",
    "explore/iib-bcov/planar_cubic_primary.py": "88b317a8ce0e2520f356773fc3aa5e957c0c0a4e",
    "explore/iib-bcov/operator/README.md": "a9f0affbe0656474aabf70144cfbad9adaaed462",
}

# q_g is used only in this program for the reduced coefficient
# (-1)^g (2*pi*i)^(5*g) Q_g|_B. Thus q_0=1 and q_1=0.
FERMAT_RATIONAL = sp.Rational(1163706689, 25609494822912000)
FERMAT_ZETA_COEFFICIENT = -sp.Rational(74897, 1976041267200)
OLD_FIXTURE = ((1, 0), (0, 1), (1, 2))
NEW_FIXTURE = ((1, 0), (0, 1), (2, 3))
FIXTURES = (
    OLD_FIXTURE,
    NEW_FIXTURE,
    ((1, 2), (-2, 3), (3, -1)),
    ((2, -1), (3, 4), (-1, 5)),
    ((1, -3), (2, 1), (-2, 4)),
)


def determinant(p, q):
    return p[0] * q[1] - p[1] * q[0]


def cyclic_symbol(a, b, c, inner_genus, outer_genus):
    """Complete cyclic composition for one ordered pair of vertex genera."""
    m, n = 2 * inner_genus + 1, 2 * outer_genus + 1
    return a**m * (b - c)**n + b**m * (c - a)**n + c**m * (a - b)**n


@cache
def nested_derivatives(p, q, r, inner_genus, outer_genus):
    """Evaluate nested differential operators by the Leibniz rule at zero.

    Each vertex of order m applies (d1 tensor d2 - d2 tensor d1)^m.
    The outer vertex differentiates the entire inner output. Its derivatives
    are distributed separately over the two inner waves with binomial weights.
    This evaluator never forms a momentum determinant or a fused momentum.
    """
    m, n = 2 * inner_genus + 1, 2 * outer_genus + 1
    value = 0
    for k in range(n + 1):
        for ell in range(m + 1):
            sign_and_weight = (-1)**(k + ell) * comb(n, k) * comb(m, ell)
            for u in range(n - k + 1):
                for v in range(k + 1):
                    value += (
                        sign_and_weight * comb(n - k, u) * comb(k, v)
                        * p[0]**(m - ell + u) * p[1]**(ell + v)
                        * q[0]**(ell + n - k - u) * q[1]**(m - ell + k - v)
                        * r[0]**k * r[1]**(n - k)
                    )
    return value


def cyclic_derivatives(fixture, inner_genus, outer_genus):
    p, q, r = fixture
    return sum(
        nested_derivatives(*arguments, inner_genus, outer_genus)
        for arguments in ((p, q, r), (q, r, p), (r, p, q))
    )


def complete_jacobiator(fixture, genus, coefficients, derivative_evaluator=False):
    p, q, r = fixture
    a, b, c = determinant(p, q), determinant(q, r), determinant(r, p)
    value = 0
    for inner in range(genus + 1):
        outer = genus - inner
        if derivative_evaluator:
            composition = cyclic_derivatives(fixture, inner, outer)
        else:
            composition = cyclic_symbol(a, b, c, inner, outer)
        value += coefficients[inner] * coefficients[outer] * composition
    return sp.expand(value)


def coefficient_of_unknown(genus, second_fixture=False):
    if second_fixture:
        return 2 * (1 + 2**(2 * genus + 2) + 4**(2 * genus + 1)
                    - 3**(2 * genus + 2))
    return 3**(2 * genus + 1) - 2**(2 * genus + 3) + 5


def certificate():
    a, b, c = sp.symbols("a b c")
    s, t = sp.symbols("s t")
    q = (sp.S.One, sp.S.Zero, *sp.symbols("q2:7"))
    zeta_three = sp.zeta(3)
    fermat = FERMAT_RATIONAL + FERMAT_ZETA_COEFFICIENT * zeta_three

    # Universal first-order polynomials, checked as exact polynomials in three
    # independent determinant variables, not only at the displayed fixtures.
    first_two = sp.expand(cyclic_symbol(a, b, c, 0, 2)
                          + cyclic_symbol(a, b, c, 2, 0))
    first_three = sp.expand(cyclic_symbol(a, b, c, 0, 3)
                            + cyclic_symbol(a, b, c, 3, 0))
    common = a * b * c * (a - b) * (a - c) * (b - c)
    assert sp.expand(first_two - 10 * common) == 0
    assert sp.expand(first_three - 7 * common
                     * (3*a*a + 3*b*b + 3*c*c - 2*a*b - 2*a*c - 2*b*c)) == 0

    # Independent nested differentiation checks every ordered genus split.
    differential_checks = 0
    for fixture in FIXTURES:
        p, qv, r = fixture
        determinants = determinant(p, qv), determinant(qv, r), determinant(r, p)
        for genus in range(7):
            for inner in range(genus + 1):
                outer = genus - inner
                assert cyclic_derivatives(fixture, inner, outer) == cyclic_symbol(
                    *determinants, inner, outer)
                differential_checks += 1

    # Closed integer formulas at the two distinguished fixtures, including all
    # genera through twelve as a check of the formula proved in the report.
    first_nonzero_values = {}
    for genus in range(13):
        row = {}
        for name, fixture, second in (
            ("r=(1,2)", OLD_FIXTURE, False),
            ("r=(2,3)", NEW_FIXTURE, True),
        ):
            p, qv, r = fixture
            determinants = determinant(p, qv), determinant(qv, r), determinant(r, p)
            direct = cyclic_symbol(*determinants, 0, genus)
            direct += cyclic_symbol(*determinants, genus, 0)
            assert direct == coefficient_of_unknown(genus, second)
            row[name] = str(direct)
        first_nonzero_values[str(genus)] = row

    # Evaluate the FULL genus coefficients before taking any difference.
    # In particular the genus-six row initially retains q2*q4.
    rows = {}
    for genus in (4, 5, 6):
        original = complete_jacobiator(OLD_FIXTURE, genus, q)
        additional = complete_jacobiator(NEW_FIXTURE, genus, q)
        assert original == complete_jacobiator(OLD_FIXTURE, genus, q, True)
        assert additional == complete_jacobiator(NEW_FIXTURE, genus, q, True)
        old_leading = int(original.coeff(q[genus]))
        new_leading = int(additional.coeff(q[genus]))
        divisor = gcd(old_leading, new_leading)
        weight_new, weight_old = old_leading // divisor, new_leading // divisor
        difference = sp.expand(weight_new * additional - weight_old * original)
        assert difference.coeff(q[genus]) == 0
        rows[str(genus)] = {
            "complete_r_1_2": str(original),
            "complete_r_2_3": str(additional),
            "weight_r_2_3": weight_new,
            "weight_r_1_2_subtracted": weight_old,
            "difference": str(difference),
        }

    assert rows["4"]["difference"] == "-57600*q2**2"
    assert rows["5"]["difference"] == "-5040000*q2*q3"
    difference_six = sp.sympify(rows["6"]["difference"])
    assert sp.expand(difference_six.subs(q[2], 0) + 20180160*q[3]**2) == 0

    values = {
        "genus_2_new_fixture": sp.factor(720 * s / 51840),
        "genus_3_new_fixture": sp.expand(20160 * t * fermat),
        "genus_4_difference": sp.factor(-57600 * (s / 51840)**2),
        "genus_5_difference": sp.expand(-5040000 * s / 51840 * t * fermat),
        "genus_6_difference_on_s_zero": sp.expand(-20180160 * t*t * fermat**2),
    }
    assert values["genus_4_difference"] == -s*s / 46656
    assert sp.expand(values["genus_3_new_fixture"] - t * (
        sp.Rational(1163706689, 1270312243200)
        - sp.Rational(74897, 98017920)*zeta_three)) == 0

    # Positive control: the normalized Moyal sine/hyperbolic sine coefficients.
    # q1=1/6 is chosen so qg=1/(2g+1)!, with the formal genus grading retained.
    for genus in range(9):
        control = sum(
            cyclic_symbol(a, b, c, inner, genus-inner)
            / sp.Integer(factorial(2*inner+1)*factorial(2*(genus-inner)+1))
            for inner in range(genus + 1)
        )
        assert sp.expand(control) == 0
    perturbed = [sp.Rational(1, factorial(2*genus+1)) for genus in range(3)]
    perturbed[2] += 1
    assert complete_jacobiator(NEW_FIXTURE, 2, perturbed) == 720

    # Arithmetic controls on inherited graph results. These are not fresh
    # integrations and are labelled accordingly in the certificate.
    theta = [sp.Rational(41, 90720), sp.Rational(11, 11340),
             sp.Rational(1, 2016), sp.Rational(1, 2520)]
    assert sum(theta) == sp.Rational(1, 432) == 120 * sp.Rational(1, 51840)
    double_square = -sp.Rational(3361, 27433728000)
    tetrahedron = (sp.Rational(13904672281, 177843714048000)
                   - sp.Rational(74897, 1176215040)*zeta_three)
    hexagon = -sp.Rational(1, 580608)
    assert sp.expand(double_square + tetrahedron + hexagon - 1680*fermat) == 0

    # Exact integral-test bound. For strictly decreasing x^(-3), the omitted
    # sum from k=9 is strictly greater than the integral from 9 to infinity.
    lower = sum((Fraction(1, k**3) for k in range(1, 9)), Fraction(1, 162))
    lower = sp.Rational(lower.numerator, lower.denominator)
    assert lower == sp.Rational(2136104771, 1778112000)
    ratio = -FERMAT_RATIONAL / FERMAT_ZETA_COEFFICIENT
    gap = lower - ratio
    upper = FERMAT_RATIONAL + FERMAT_ZETA_COEFFICIENT*lower
    assert gap == sp.Rational(327281302787, 133175254464000) > 0
    assert upper == -sp.Rational(327281302787, 3513622689703526400000) < 0

    # The inherited Hesse discriminant identity ensures one of the measured
    # lower coefficients is nonzero whenever the projective cubic is smooth.
    hesse_parameter = sp.symbols("b")
    hesse_s = hesse_parameter**4 - hesse_parameter
    hesse_t = 1 - 20*hesse_parameter**3 - 8*hesse_parameter**6
    assert sp.expand(hesse_t**2 - 64*hesse_s**3
                     - (1 + 8*hesse_parameter**3)**3) == 0

    return {
        "status": "PASS",
        "schema": "planar-bulk-ward-certificate-v1",
        "provenance": {
            "repository": REPOSITORY,
            "pinned_commit": PINNED_COMMIT,
            "source_git_blob_sha1": SOURCE_BLOBS,
            "reduced_coefficient_convention": "q_g = (-1)^g (2*pi*i)^(5*g) Q_g|_B; q_0=1; q_1=0",
            "inherited_analytic_input": {
                "q2": "s/51840",
                "q3_rational_times_t": str(FERMAT_RATIONAL),
                "q3_zeta3_coefficient_times_t": str(FERMAT_ZETA_COEFFICIENT),
            },
            "scope": "Exact cochain-level cyclic compositions of isolated local binary operations; no new period integration; no claim of descent to D3 charges.",
        },
        "verification": {
            "arithmetic": "integer and rational symbolic arithmetic only",
            "independent_derivative_fixture_count": len(FIXTURES),
            "independent_derivative_genus_range_inclusive": [0, 6],
            "independent_ordered_genus_split_checks": differential_checks,
            "first_nonzero_fixture_formula_genus_range_inclusive": [0, 12],
            "moyal_complete_polynomial_genus_range_inclusive": [0, 8],
            "deliberate_q2_perturbation_jacobiator_at_new_fixture": "720",
            "universal_first_order_factorizations": [str(sp.factor(first_two)), str(sp.factor(first_three))],
        },
        "momentum_fixtures": [list(map(list, fixture)) for fixture in FIXTURES],
        "first_nonzero_integer_coefficients": first_nonzero_values,
        "complete_jacobiators_and_unknown_coefficient_cancellation": rows,
        "new_exact_BCOV_values": {key: str(value) for key, value in values.items()},
        "strict_exact_fermat_sign_bound": {
            "zeta3_strict_lower_bound": str(lower),
            "rational_zero_threshold": str(ratio),
            "lower_bound_minus_threshold": str(gap),
            "fermat_coefficient_strict_upper_bound": str(upper),
        },
        "inherited_graph_arithmetic_controls_not_reintegrated": {
            "genus2_theta_entries": list(map(str, theta)),
            "genus2_theta_sum": str(sum(theta)),
            "genus3_double_square": str(double_square),
            "genus3_tetrahedron": str(tetrahedron),
            "genus3_hexagon_with_genus2_insertion": str(hexagon),
            "genus3_sum": str(sp.expand(double_square+tetrahedron+hexagon)),
        },
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_suffix(".json"))
    parser.add_argument("--check", action="store_true",
                        help="Check that the committed JSON agrees, without rewriting it.")
    args = parser.parse_args()
    result = certificate()
    if args.check:
        assert json.loads(args.output.read_text()) == result, "Certificate JSON differs."
        action = "Checked"
    else:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        action = "Wrote"
    print(f"{action} {args.output.name}: "
          f"{result['verification']['independent_ordered_genus_split_checks']} "
          "independent differential checks; "
          "Moyal polynomial control through genus 8; "
          "exact genus 4, 5, 6 Ward differences verified.")


if __name__ == "__main__":
    main()
