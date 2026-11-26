#!/usr/bin/env python3
"""Exact Laurent-polyvector and Cech source-residue check; not a Q_n match."""
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import json

# Coordinates: z1,z2,w1,w2,w3. Keys are (Laurent exponents, exterior mask).
ZERO_EXP = (0, 0, 0, 0, 0)
TANGENT = 3
TOP = 31


def clean(value):
    return {key: Fraction(coefficient) for key, coefficient in value.items() if coefficient}


def term(exponents, mask, coefficient=1):
    return clean({(tuple(exponents), mask): coefficient})


def add(*values):
    result = {}
    for value in values:
        for key, coefficient in value.items():
            result[key] = result.get(key, Fraction()) + coefficient
    return clean(result)


def scale(value, coefficient):
    return clean({key: coefficient * entry for key, entry in value.items()})


def wedge(left, right):
    result = {}
    for (exponents_a, mask_a), coefficient_a in left.items():
        for (exponents_b, mask_b), coefficient_b in right.items():
            if mask_a & mask_b:
                continue
            inversions = sum(
                bool(mask_a & (1 << a)) and bool(mask_b & (1 << b))
                for a in range(5) for b in range(a)
            )
            key = (tuple(a + b for a, b in zip(exponents_a, exponents_b)), mask_a | mask_b)
            result[key] = result.get(key, Fraction()) + (-1) ** inversions * coefficient_a * coefficient_b
    return clean(result)


def divergence(value):
    result = {}
    for (exponents, mask), coefficient in value.items():
        for coordinate in range(5):
            if not (mask & (1 << coordinate)) or not exponents[coordinate]:
                continue
            sign = (-1) ** ((mask & ((1 << coordinate) - 1)).bit_count())
            new_exponents = list(exponents)
            new_exponents[coordinate] -= 1
            key = (tuple(new_exponents), mask ^ (1 << coordinate))
            result[key] = result.get(key, Fraction()) + sign * exponents[coordinate] * coefficient
    return clean(result)


def bracket(left, right):
    result = {}
    # The convention is [vector, polyvector] = Lie derivative.
    for key, coefficient in left.items():
        one = {key: coefficient}
        degree = key[1].bit_count()
        inside = add(
            divergence(wedge(one, right)),
            scale(wedge(divergence(one), right), -1),
            scale(wedge(one, divergence(right)), -(-1) ** degree),
        )
        result = add(result, scale(inside, 1 if degree % 2 else -1))
    return result


def parity(order):
    return (-1) ** sum(order[a] > order[b] for a in range(3) for b in range(a + 1, 3))


def root_vector(i, j):
    exponents = [0] * 5
    exponents[2 + j] = 1
    return term(exponents, 1 << (2 + i))


def cartan(h):
    result = {}
    for i in range(3):
        exponents = [0] * 5
        exponents[2 + i] = 1
        result = add(result, term(exponents, 1 << (2 + i), h[i]))
    return result


def normalizable_mode(i, j, k):
    exponents_1 = [0, 0, -1, -1, -1]
    exponents_1[2 + j] = -3
    exponents_2 = [0, 0, -1, -1, -1]
    exponents_2[2 + j] = -2
    exponents_2[2 + k] = -2
    sign = parity((i, j, k))
    return add(
        term(exponents_1, TANGENT | (1 << (2 + k)), 2 * sign),
        term(exponents_2, TANGENT | (1 << (2 + j)), -sign),
    )


def serialize(value):
    return [
        {"exponents": list(exponents), "exterior_mask": mask, "coefficient": str(coefficient)}
        for (exponents, mask), coefficient in sorted(value.items())
    ]


def main():
    flux = term((0, 0, -1, -1, -1), TANGENT)
    assert not divergence(flux)
    assert not bracket(flux, flux)
    residue_key = ((0, 0, -1, -1, -1), TOP)
    completed_sources = []
    checked_residues = 0
    examples = []
    for i, j, k in permutations(range(3)):
        mode = normalizable_mode(i, j, k)
        source = root_vector(i, j)
        flux_bracket = bracket(flux, source)
        expected_exponents = [0, 0, -1, -1, -1]
        expected_exponents[2 + i] = -2
        expected_exponents[2 + j] = 0
        assert flux_bracket == term(expected_exponents, TANGENT)
        # This cochain has only the chart-pair component missing j.
        # delta has coefficient (-1)^j on that component.
        primitive = scale(flux_bracket, (-1) ** (j + 1))
        cech_delta_primitive = scale(primitive, (-1) ** j)
        assert not add(cech_delta_primitive, flux_bracket)
        assert not divergence(source)
        assert not divergence(primitive)
        assert not bracket(flux, primitive)
        assert not divergence(mode)
        assert not bracket(flux, mode)
        assert not wedge(mode, primitive)
        completed_sources.append({
            "root": [i + 1, j + 1],
            "cech_primitive_chart_pair": [a + 1 for a in range(3) if a != j],
            "flux_bracket": serialize(flux_bracket),
            "primitive": serialize(primitive),
            "full_linear_closure": True,
            "primitive_does_not_change_residue": True,
        })
        for h1 in range(-4, 5):
            for h2 in range(-4, 5):
                h = (h1, h2, -h1 - h2)
                regular = cartan(h)
                assert not divergence(regular)
                assert not bracket(flux, regular)
                # The bulk residue equals the actual Cartan root weight.
                expected = Fraction(h[j] - h[i])
                assert bracket(regular, source) == scale(source, expected)
                product = wedge(wedge(mode, source), regular)
                assert product == term(residue_key[0], residue_key[1], expected)
                residue = product.get(residue_key, Fraction())
                assert residue == expected
                checked_residues += 1
        if (i, j, k) == (0, 1, 2):
            for h in ((0, 1, -1), (1, -1, 0)):
                product = wedge(wedge(mode, source), cartan(h))
                examples.append({"cartan": list(h), "top_product": serialize(product), "normalized_residue": str(product.get(residue_key, 0))})
    result = {
        "status": "PASS",
        "scope": "Exact transverse cyclic residue of flux-closed source lifts and a normalizable mode. Not an integrated bulk amplitude and not a Q_n match.",
        "coordinate_order": ["z1", "z2", "w1", "w2", "w3"],
        "cech_delta_convention": "K23-K13+K12",
        "residue_normalization": "(2*pi*i)^(-3) times product-circle integral",
        "completed_root_sources": completed_sources,
        "cartan_residue_checks": checked_residues,
        "exact_general_result": "residue(mode_ij wedge source_ij wedge H) = h_j-h_i for h_1+h_2+h_3=0",
        "examples": examples,
        "unresolved": ["radial boundary action", "tangential external kernels", "physical source/operator normalization", "quantum composite-observable transfer", "comparison with Q_n"],
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS 6 full Cech source lifts and", checked_residues, "exact Cartan residue checks")


if __name__ == "__main__":
    main()
