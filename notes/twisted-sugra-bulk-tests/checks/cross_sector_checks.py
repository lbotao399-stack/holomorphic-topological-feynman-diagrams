#!/usr/bin/env python3
"""Independent exact arithmetic for the bulk measurement and Ward targets."""
from fractions import Fraction as F
from itertools import permutations
from math import factorial
from pathlib import Path
import json


def determinant(p, q):
    return p[0]*q[1] - p[1]*q[0]


def ordered_nested_coefficient(p, q, r, left_genus, right_genus):
    result = 0
    for first, second, third in ((p, q, r), (q, r, p), (r, p, q)):
        merged = (first[0]+second[0], first[1]+second[1])
        result += (determinant(first, second)**(2*left_genus+1)
                   * determinant(merged, third)**(2*right_genus+1))
    return result


def jacobi_tests():
    p, q = (1, 0), (0, 1)
    rows = {}
    for genus, factors in ((4, (7, 162)), (5, (32, 1389)), (6, (13, 1038))):
        probes = {}
        for name, r in (("first", (1, 2)), ("second", (2, 3))):
            linear = (ordered_nested_coefficient(p, q, r, 0, genus)
                      + ordered_nested_coefficient(p, q, r, genus, 0))
            if genus == 4:
                quadratic = ordered_nested_coefficient(p, q, r, 2, 2)
            elif genus == 5:
                quadratic = (ordered_nested_coefficient(p, q, r, 2, 3)
                             + ordered_nested_coefficient(p, q, r, 3, 2))
            else:
                quadratic = ordered_nested_coefficient(p, q, r, 3, 3)
            probes[name] = {"unknown_genus_coefficient": linear,
                            "known_product_coefficient": quadratic}
        second_weight, first_weight = factors
        assert (second_weight*probes["second"]["unknown_genus_coefficient"]
                - first_weight*probes["first"]["unknown_genus_coefficient"]) == 0
        difference = (second_weight*probes["second"]["known_product_coefficient"]
                      - first_weight*probes["first"]["known_product_coefficient"])
        assert difference == {4: -57600, 5: -5040000, 6: -20180160}[genus]
        rows[str(genus)] = {"probes": probes, "second_weight": second_weight,
                            "first_weight": first_weight,
                            "remaining_product_coefficient": difference}
    assert F(-57600, 51840**2) == -F(1, 46656)
    return rows


def genus_four_measurement():
    # Complete inherited three-mark periods; no period integration is repeated.
    periods = {(4, 0, 0): F(731, 1024192512000),
               (3, 1, 0): F(1493, 438939648000),
               (2, 2, 0): F(157, 29262643200),
               (2, 1, 1): F(433, 37623398400)}
    orbit_weights = {(4, 0, 0): 2, (3, 1, 0): -1,
                     (2, 2, 0): 1, (2, 1, 1): 0}
    weights = {}
    total = F(0)
    for levels, coefficient in orbit_weights.items():
        orbit = set(permutations(levels))
        total += len(orbit)*coefficient*periods[levels]
        for ordering in orbit:
            weights[ordering] = coefficient
    lowerings = []
    for a in range(4):
        for b in range(4-a):
            c = 3-a-b
            lower = (a, b, c)
            total_lowering = sum(weights[tuple(lower[j]+int(i==j) for j in range(3))]
                                 for i in range(3))
            assert total_lowering == 0
            lowerings.append(lower)
    assert len(weights) == 15 and len(lowerings) == 10
    assert total == -F(31, 1024192512000)
    fiber_degree = factorial(4)**3
    profile_divisor = factorial(3)**2 * factorial(4)
    elementary_output = -total*fiber_degree/profile_divisor
    assert elementary_output == F(31, 64012032000)
    return {
        "ordered_descendant_coefficients": len(weights),
        "independent_single_donor_equations": len(lowerings),
        "planar_period_combination": str(total),
        "vertical_abel_fiber_degree": fiber_degree,
        "fixed_profile_factorial_divisor": profile_divisor,
        "normalized_genus_four_elementary_response": str(elementary_output),
        "normalization": "(2*pi*i)^20 times the selected elementary scalar response",
        "lower_genus_rank_inequalities": [
            {"genus": genus, "required_vertical_marks": 2*genus+4,
             "available_vertical_complex_rank": 3*genus}
            for genus in range(4)],
        "scope": "Selected polarized flat bulk response using inherited complete periods.",
    }


def main():
    result = {"status": "PASS", "arithmetic": "exact integers and fractions",
              "jacobi": jacobi_tests(), "genus_four": genus_four_measurement()}
    destination = Path(__file__).resolve().parents[1] / "verification/cross_sector_checks.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status": "PASS",
                      "genus_four_response": result["genus_four"]["normalized_genus_four_elementary_response"],
                      "unknown_period_free_jacobi_tests": 3}))


if __name__ == "__main__":
    main()
