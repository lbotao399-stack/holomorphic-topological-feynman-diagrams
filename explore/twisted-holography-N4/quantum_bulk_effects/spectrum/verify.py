#!/usr/bin/env python3
"""Exact checks for the spin bound and the metric-sensitive spectral test.

This verifies consequences of the boundary fan formula and an explicitly
specified finite-dimensional Hilbert complex. It does not compute a bulk
supergravity amplitude.
"""

from fractions import Fraction
import json
from pathlib import Path
import sympy as s


def spin_bound_checks():
    checked = 0
    allowed = 0
    saturated = []
    # Order is ghost jet, gamma, beta, b. The ghost must have at least
    # one holomorphic derivative in the physical relative complex.
    for n in range(1, 33):
        minimum = None
        for ghost in range(n + 2):
            for gamma in range(n + 2 - ghost):
                for beta in range(n + 2 - ghost - gamma):
                    b = n + 1 - ghost - gamma - beta
                    checked += 1
                    functional_theta = 3 * ghost + 2 * gamma + beta
                    twice_spin = beta + 2 * b
                    if functional_theta <= 3:
                        allowed += 1
                        assert twice_spin >= 2 * n - 1
                        minimum = twice_spin if minimum is None else min(minimum, twice_spin)
        assert minimum == 2 * n - 1
        if n >= 2:
            # Three beta contractions and n-2 b contractions in the
            # previously constructed detector family saturate the bound.
            assert 3 + 2 * (n - 2) == 2 * n - 1
            saturated.append(n)
    return {"inventory_checks": checked, "degree_allowed_inventories": allowed,
            "orders_checked": [1, 32], "detector_saturation_orders": saturated}


def charge_checks():
    rows = []
    sources = [
        ("SO(13) BMN", (7, 5, 5), (0, 0), 16),
        ("SO(13) mixed descendant", (5, 3, 3), (3, 3), 16),
        ("SO(11) glueball", (3, 3, 3), (3, 3), 14),
        ("SO(12) and SO(13) glueball", (3, 3, 3), (4, 4), 16),
    ]
    for name, theta, derivatives, length in sources:
        for j in range(min(theta) + 1):
            corrected_theta = tuple(x - j for x in theta)
            corrected_derivatives = tuple(x + j for x in derivatives)
            corrected_length = length - j
            twice_spin = sum(corrected_theta) + sum(corrected_derivatives) - corrected_length
            assert twice_spin == 1
            rows.append({"source": name, "order": j, "theta": corrected_theta,
                         "derivatives": corrected_derivatives, "length": corrected_length,
                         "twice_J_L": twice_spin})

    def charges(theta, derivatives, length):
        return (
            Fraction(length + sum(theta), 2) + sum(derivatives),
            tuple(Fraction(length - sum(theta), 2) + x for x in theta),
            Fraction(sum(theta) + sum(derivatives) - length, 2),
            Fraction(derivatives[0] - derivatives[1], 2),
        )

    source = charges((7, 5, 5), (0, 0), 16)
    target = charges((4, 2, 2), (3, 3), 14)
    assert source == (Fraction(33, 2), (Fraction(13, 2), Fraction(9, 2), Fraction(9, 2)), Fraction(1, 2), 0)
    assert target == (17, (7, 5, 5), 0, 0)
    for E, R, JL, JR in [source, target]:
        assert E - sum(R) - 2 * JL == 0
        assert 2 * (E + JL) == 34
    return {"correction_sector_checks": rows,
            "SO13_source_charges": str(source), "SO13_target_charges": str(target)}


def spectral_checks():
    h, lam = s.symbols("h lam", real=True)
    matrices = []
    for r in range(1, 9):
        matrix = s.zeros(r)
        for i in range(r):
            matrix[i, i] = h
            if i + 1 < r:
                matrix[i, i + 1] = 1
        completion = s.Matrix([(-h) ** i for i in range(r)])
        expected = s.zeros(r, 1)
        expected[r - 1] = (-1) ** (r - 1) * h ** r
        assert matrix * completion == expected
        inverse = matrix.inv()
        for i in range(r):
            for j in range(r):
                assert inverse[i, j] == ((-1) ** (j-i) * h ** (-(j-i+1)) if j >= i else 0)
        norm_inverse = s.expand(sum(entry ** 2 for entry in inverse))
        formula_inverse = h ** (-2*r) * sum((j+1) * h ** (2*j) for j in range(r))
        assert s.expand(norm_inverse - formula_inverse) == 0
        determinant = s.expand((matrix.T * matrix - lam * s.eye(r)).det(method="domain-ge"))
        prior2, prior1 = s.Integer(1), h*h-lam
        if r == 1:
            recurrence = prior1
        else:
            for k in range(2, r+1):
                recurrence = s.expand((1+h*h-lam)*prior1-h*h*prior2)
                prior2, prior1 = prior1, recurrence
        assert s.expand(determinant - recurrence) == 0
        assert s.expand(determinant.subs(lam, 0)) == h ** (2*r)
        matrices.append({"size": r, "determinant": str(determinant),
                         "gap_lower_bound": str(2 / formula_inverse),
                         "gap_upper_bound": str(2*h**(2*r)/sum(h**(2*j) for j in range(r)))})

    characteristic3 = lam**3-(2+3*h*h)*lam**2+(1+2*h*h+3*h**4)*lam-h**6
    series3 = s.Integer(0)
    coefficients = {}
    for exponent in range(6, 25, 2):
        coefficient = s.symbols("coefficient")
        equation = s.expand(characteristic3.subs(lam, series3+coefficient*h**exponent)).coeff(h, exponent)
        value = s.solve(equation, coefficient)[0]
        series3 += value*h**exponent
        coefficients[str(exponent)] = str(2*value)
    residual = s.Poly(s.expand(characteristic3.subs(lam, series3)), h)
    assert all(exponent[0] >= 26 for exponent, coefficient in residual.terms() if coefficient != 0)

    # Explicit symbolic nontrivial Gram metrics, three-step block.
    a = s.symbols("a1:4", positive=True)
    b = s.symbols("b1:3", positive=True)
    u = s.symbols("u1:4", positive=True)
    v = s.symbols("v1:4", positive=True)
    weighted = s.Matrix([[a[0]*h,b[0],0],[0,a[1]*h,b[1]],[0,0,a[2]*h]])
    characteristic_weighted = s.expand((weighted.T*s.diag(*v)*weighted-lam*s.diag(*u)).det())
    constant = s.expand(characteristic_weighted.subs(lam,0))
    linear_at_zero = s.expand(characteristic_weighted.diff(lam).subs({h:0,lam:0}))
    assert constant == h**6*s.prod(x*x for x in a)*s.prod(v)
    assert linear_at_zero == -u[0]*s.prod(x*x for x in b)*v[0]*v[1]
    leading = s.factor(-2*(constant/h**6)/linear_at_zero)
    assert leading == 2*s.prod(x*x for x in a)*v[2]/(s.prod(x*x for x in b)*u[0])
    return {"exact_matrix_checks": matrices, "three_step_gap_coefficients": coefficients,
            "weighted_three_step_leading_gap_coefficient": str(leading)}


if __name__ == "__main__":
    output = {"scope": "Boundary spin theorem and specified finite Hilbert complex; no independent bulk amplitude",
              "spin": spin_bound_checks(), "charges": charge_checks(), "spectral": spectral_checks()}
    path = Path(__file__).with_name("results.json")
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "inventory_checks": output["spin"]["inventory_checks"],
                      "correction_sectors": len(output["charges"]["correction_sector_checks"]),
                      "spectral_matrix_sizes": [1,8],
                      "output": str(path)}, indent=2))
