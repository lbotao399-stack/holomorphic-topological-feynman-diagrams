#!/usr/bin/env python3
"""Exact checks for the c-b fan sign, auxiliary gauge defect, and N=2 check.

The paper's functional derivative is evaluated at ig=1 only to strip an
invertible normalization.  The restored g and pi powers are recorded in the
report; this program checks all Grassmann, SU(2), and rational factors.
No full absolute quantum differential or corrected brane complex is assumed.
"""
from fractions import Fraction
import json
from pathlib import Path

# Odd order theta1, theta2, theta3, c, beta1, psi, eps_plus, eps_minus.
# Polynomial keys consist of the odd bit mask and powers of gamma1 and b.
def add(*polys):
    result = {}
    for poly in polys:
        for key, value in poly.items():
            result[key] = result.get(key, Fraction(0)) + value
    return {key: value for key, value in result.items() if value}


def mul(left, right):
    result = {}
    for (mask1, gamma1, b1), value1 in left.items():
        for (mask2, gamma2, b2), value2 in right.items():
            if mask1 & mask2:
                continue
            inversions = sum((mask2 & ((1 << i) - 1)).bit_count()
                             for i in range(8) if (mask1 >> i) & 1)
            key = (mask1 | mask2, gamma1 + gamma2, b1 + b2)
            value = (-1)**inversions * value1 * value2
            result[key] = result.get(key, Fraction(0)) + value
    return {key: value for key, value in result.items() if value}


def odd_derivative(poly, index):
    result = {}
    for (mask, gamma, b), value in poly.items():
        if not ((mask >> index) & 1):
            continue
        sign = (-1)**((mask & ((1 << index) - 1)).bit_count())
        result[(mask ^ (1 << index), gamma, b)] = sign * value
    return result


def even_derivative(poly, index):
    result = {}
    for key, value in poly.items():
        power = key[index]
        if power:
            new_key = list(key)
            new_key[index] -= 1
            result[tuple(new_key)] = value * power
    return result


def monomial(mask=0, gamma=0, b=0, coefficient=1):
    return {(mask, gamma, b): Fraction(coefficient)}


def superfield_derivative(poly):
    # A.4, with ig=1: theta123 d_c - theta23 d_gamma - theta1 d_beta - d_b.
    return add(
        mul(monomial(7), odd_derivative(poly, 3)),
        mul(monomial(6, coefficient=-1), even_derivative(poly, 1)),
        mul(monomial(1, coefficient=-1), odd_derivative(poly, 4)),
        mul(monomial(coefficient=-1), even_derivative(poly, 2)),
    )


def epsilon(a, b, c):
    if len({a, b, c}) < 3:
        return 0
    return (-1)**sum(x > y for i, x in enumerate((a, b, c))
                       for y in (a, b, c)[i + 1:])


def gaussian_product(*values):
    """Multiply exact (real, imaginary) Fraction pairs, with g=1."""
    real, imaginary = Fraction(1), Fraction(0)
    for other_real, other_imaginary in values:
        other_real, other_imaginary = Fraction(other_real), Fraction(other_imaginary)
        real, imaginary = (real * other_real - imaginary * other_imaginary,
                           real * other_imaginary + imaginary * other_real)
    return real, imaginary


def main():
    gamma_beta = superfield_derivative(superfield_derivative(monomial(1 << 4, gamma=1)))
    c_b = superfield_derivative(superfield_derivative(monomial(1 << 3, b=1)))
    c_b_psi = superfield_derivative(superfield_derivative(monomial((1 << 3) | (1 << 5), b=1)))
    assert gamma_beta == monomial(7, coefficient=2)
    assert c_b == monomial(7, coefficient=-2)
    assert c_b_psi == monomial(7 | (1 << 5), coefficient=-2)

    # Restore the phase directly: fan; two assignments; 1/(ig); -ig;
    # the two ig ghost endpoints.  g=1 here only removes its known power.
    phase = gaussian_product((Fraction(-1, 128), 0), (2, 0),
                             (0, -1), (0, -1), (0, 1), (0, 1))
    assert phase == (Fraction(-1, 64), 0)
    auxiliary_phase = gaussian_product(phase, (0, 1), (2, 0))
    assert auxiliary_phase == (0, Fraction(-1, 32))

    # The root spinor invariant has two ordered odd-wave contributions.
    eps_plus, eps_minus = monomial(1 << 6), monomial(1 << 7)
    invariant = add(mul(eps_plus, eps_minus),
                    mul(monomial(coefficient=-1), mul(eps_minus, eps_plus)))
    assert invariant == monomial((1 << 6) | (1 << 7), coefficient=2)

    # SU(2), B=3, differentiated ghosts both in color 1.
    color = [sum(epsilon(a, r, 1) * epsilon(3, r, 1) for r in (1, 2, 3))
             for a in (1, 2, 3)]
    assert color == [0, 0, 1]

    # Remove g^2/pi^2 from determinant and trace answers at N=2, tau=1/2.
    determinant = Fraction(1, 32) * Fraction(1, 2) * (-4)
    known_trace = -Fraction(2, 64) * 2
    assert determinant == known_trace == -Fraction(1, 16)

    data = {
        "status": "PASS",
        "superfield_derivative_gamma_beta_theta123": 2,
        "superfield_derivative_c_b_theta123": -2,
        "superfield_derivative_c_b_psi_theta123_psi": -2,
        "exact_phase_restoration_Q1_c_b_real_imaginary": [str(x) for x in phase],
        "exact_phase_restoration_su2_auxiliary_real_imaginary": [str(x) for x in auxiliary_phase],
        "ordered_two_ghost_spinor_factor": 2,
        "su2_auxiliary_color_vector_B3_ghosts1": color,
        "Q1_c_b_coefficient_in_g2_over_pi2": "-1/64",
        "su2_auxiliary_defect_coefficient_in_i_g3_over_pi2_t3_psi": "-1/32",
        "N2_determinant_and_known_trace_coefficient_in_g2_over_pi2": str(determinant),
        "scope": "Exact evaluation of the unreduced one-loop fan; no absolute-extension or brane-QME theorem asserted."
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(data, indent=2) + "\n")
    print("PASS: three ordered Grassmann derivative checks, two exact complex phases, ghost spinor sign, SU(2) color, N=2 determinant normalization")


if __name__ == "__main__":
    main()
