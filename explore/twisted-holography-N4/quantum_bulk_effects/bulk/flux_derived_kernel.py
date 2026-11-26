#!/usr/bin/env python3
"""Exact bulk flux derived-bracket calculation; no boundary Q_n is imported."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import json


def add(*values):
    result = {}
    for value in values:
        for key, coefficient in value.items():
            result[key] = result.get(key, Fraction()) + coefficient
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def term(exponents, mask=0, coefficient=1):
    return {(tuple(exponents), mask): Fraction(coefficient)} if coefficient else {}


def scale(value, coefficient):
    return {key: coefficient * entry for key, entry in value.items() if coefficient * entry}


def wedge(left, right):
    result = {}
    for (a, ma), ca in left.items():
        for (b, mb), cb in right.items():
            if ma & mb:
                continue
            inversions = sum(bool(ma & (1 << i)) and bool(mb & (1 << j)) for i in range(5) for j in range(i))
            key = (tuple(x + y for x, y in zip(a, b)), ma | mb)
            result[key] = result.get(key, Fraction()) + (-1) ** inversions * ca * cb
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def divergence(value):
    result = {}
    for (exponents, mask), coefficient in value.items():
        for i in range(5):
            if not mask & (1 << i) or not exponents[i]:
                continue
            sign = (-1) ** ((mask & ((1 << i) - 1)).bit_count())
            powers = list(exponents)
            powers[i] -= 1
            key = (tuple(powers), mask ^ (1 << i))
            result[key] = result.get(key, Fraction()) + sign * exponents[i] * coefficient
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def bracket(left, right):
    result = {}
    for key, coefficient in left.items():
        one = {key: coefficient}
        degree = key[1].bit_count()
        result = add(result, scale(add(
            divergence(wedge(one, right)),
            scale(wedge(divergence(one), right), -1),
            scale(wedge(one, divergence(right)), -(-1) ** degree),
        ), -1 if degree % 2 == 0 else 1))
    return result


def diagonal(weights):
    values = []
    for i, weight in enumerate(weights):
        powers = [0] * 5
        powers[i] = 1
        values.append(term(powers, 1 << i, weight))
    return add(*values)


def main():
    flux = term((0, 0, -1, -1, -1), 3)
    scalar_flux = term((0, 0, -1, -1, -1))
    assert not divergence(flux)
    assert not bracket(flux, flux)

    polynomial_tests = 0
    nonzero_polynomial_tests = 0
    for a, b, c, d in product(range(6), repeat=4):
        first = term((a, b, 0, 0, 0))
        second = term((c, d, 0, 0, 0))
        actual = bracket(bracket(flux, first), second)
        determinant = a * d - b * c
        expected = term((a + c - 1, b + d - 1, -1, -1, -1), coefficient=-determinant)
        assert actual == expected, (a, b, c, d, actual, expected)
        assert bracket(bracket(flux, second), first) == scale(actual, -1)
        assert not divergence(actual)
        polynomial_tests += 1
        nonzero_polynomial_tests += bool(actual)

    assert bracket(bracket(flux, term((1, 0, 0, 0, 0))), term((0, 1, 0, 0, 0))) == scale(scalar_flux, -1)

    diagonal_tests = 0
    for weights in product(range(-2, 3), repeat=5):
        value = diagonal(weights)
        assert divergence(value) == term((0, 0, 0, 0, 0), coefficient=sum(weights))
        assert bracket(flux, value) == scale(flux, sum(weights))
        diagonal_tests += 1

    tangent_ideal_tests = 0
    for powers in product(range(3), repeat=5):
        for mask in range(32):
            actual = bracket(flux, term(powers, mask))
            assert all(output_mask & 3 for (_, output_mask) in actual)
            tangent_ideal_tests += 1

    result = {
        'status': 'PASS',
        'scope': 'Independent classical flux-derived bulk kernel; not an evaluated quantum Ward amplitude or a Q1 coefficient match.',
        'coordinate_order': ['z1', 'z2', 'w1', 'w2', 'w3'],
        'schouten_convention': '[vector, function] = vector(function)',
        'polynomial_pair_checks': polynomial_tests,
        'nonzero_polynomial_pair_checks': nonzero_polynomial_tests,
        'diagonal_scaling_checks': diagonal_tests,
        'linear_flux_tangent_ideal_checks': tangent_ideal_tests,
        'unit_example': '[[F,z1],z2] = -1/(w1*w2*w3)',
        'general_polynomial_kernel': '-N*(a*d-b*c)*z1^(a+c-1)*z2^(b+d-1)/(w1*w2*w3)',
        'plane_wave_kernel': '-N*(p1*q2-p2*q1)*exp((p+q).z)/(w1*w2*w3)',
        'normal_euler_relation': 'D(sum_i w_i d/dw_i) = 3*t + 3*N*F',
        'linear_type_obstruction': 'Every [F,regular polyvector] has at least one tangent vector, so it cannot equal the scalar glueball mode.',
        'unresolved': ['full composite/source dictionary', 'normalization relating bulk source fields to N=4 fields', 'Ward graph and external kernels', 'higher Q_n matching'],
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    print('PASS', polynomial_tests, 'polynomial pairs,', diagonal_tests, 'Euler weights,', tangent_ideal_tests, 'linear type tests')


if __name__ == '__main__':
    main()
