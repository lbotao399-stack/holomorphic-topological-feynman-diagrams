"""Independently replay the manuscript's complete finite wheel base range.

This checker imports no production graph or residue evaluator.  It checks
the closed principal coefficients against finite Taylor multiplication,
forms every convolution coefficient, retains and cancels the polynomial
part, and computes every rational and zeta(3) coefficient for q=3,...,20.
The all-order argument is in sections/periods.tex; the q>=21 bound is
checked here through its exact starting fraction and polynomial identity.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb, factorial, prod
from pathlib import Path
import argparse
import json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REFERENCE = ROOT / "explore/circuit-general/core-period-closure/certificate.json"
SECTION = HERE.parent / "sections/periods.tex"
OUTPUT = HERE.parent / "verification/wheel_residue.json"
EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def digest(value):
    return sha256(canonical(value)).hexdigest()


def harmonic(last, power=1):
    return sum((Fraction(1, k ** power) for k in range(1, last + 1)), Fraction())


def cycle_rank(selected):
    parent = list(range(4))

    def root(vertex):
        while parent[vertex] != vertex:
            vertex = parent[vertex]
        return vertex

    rank = 0
    for edge in selected:
        left, right = map(root, EDGES[edge])
        if left != right:
            parent[left] = right
            rank += 1
    return len(selected) - rank


def convergence(q):
    weights = (q - 1, 1, 1, 1, q - 1, q - 1)
    require(q >= 3 and all(1 <= a < q for a in weights), "invalid integer K4 weights")
    require(sum(weights) == 3 * q, "wrong projective homogeneity")
    margins = []
    for mask in range(1, 63):
        selected = tuple(i for i in range(6) if mask & (1 << i))
        margin = sum(weights[i] for i in selected) - q * cycle_rank(selected)
        require(margin > 0, "a K4 scalar weight is not absolutely convergent")
        margins.append(margin)
    return weights, margins


def histogram_from_original_weights(weights, q):
    spokes, rims = weights[:3], weights[3:]
    drift = 0
    nodes = []
    for i in range(3):
        drift += 2 * spokes[i] + rims[i - 1] + rims[i] - 2 * q
        defect = rims[i] - 1
        nodes.extend(range(drift - defect, drift + defect + 1, 2))
    require(drift == 0, "the radial drift does not close")
    minimum = min(nodes)
    require(all((node - minimum) % 2 == 0 for node in nodes), "node parity failed")
    counts = Counter((node - minimum) // 2 for node in nodes)
    histogram = tuple(counts[j] for j in range(max(counts) + 1))
    prefactor = prod(factorial(q - a - 1) * factorial(rho - 1)
                     for a, rho in zip(spokes, rims))
    require(histogram == (2,) * (q - 2) + (3,), "wrong occupied multiplicities")
    require(prefactor == factorial(q - 2) ** 4, "wrong complete gamma normalization")
    require(len(histogram) - 1 < q - 1, "the convolution contour strip is empty")
    return histogram, prefactor


def closed_principal_coefficients(r):
    top = Fraction(1, factorial(r) ** 2)
    H, H2 = harmonic(r), harmonic(r, 2)
    coefficients = {(r, 3): top, (r, 2): -2 * H * top,
                    (r, 1): (2 * H * H + H2) * top}
    for j in range(r):
        highest = -Fraction(1, factorial(j) ** 2 * factorial(r - j) ** 2 * (r - j))
        slope = -Fraction(3, j - r) - 2 * sum(
            (Fraction(1, j - l) for l in range(r) if l != j), Fraction())
        coefficients[j, 2] = highest
        coefficients[j, 1] = highest * slope
    return coefficients


def taylor_principal_coefficients(histogram):
    """Independent truncated products of (j-l+z)^(-N_l), degree at most two."""
    answer = {}
    for j, Nj in enumerate(histogram):
        series = [Fraction(1)] + [Fraction()] * (Nj - 1)
        for l, Nl in enumerate(histogram):
            if l == j:
                continue
            factor = [Fraction((-1) ** k * comb(Nl + k - 1, k), (j - l) ** (Nl + k))
                      for k in range(Nj)]
            changed = [Fraction()] * Nj
            for first, value in enumerate(series):
                for second, other in enumerate(factor[:Nj - first]):
                    changed[first + second] += value * other
            series = changed
        for k in range(1, Nj + 1):
            answer[j, k] = series[Nj - k]
    return answer


def multiplicity_numerator(r):
    coefficients = [0, 0, 2]
    for j in range(1, r + 1):
        changed = [0] * (len(coefficients) + 2)
        for degree, value in enumerate(coefficients):
            changed[degree] -= j * j * value
            changed[degree + 2] += value
        coefficients = changed
    return coefficients


def polynomial_quotient(coefficients, h, power):
    """Exact integer long division by the monic polynomial (s-h)^power."""
    denominator = [comb(power, k) * (-h) ** (power - k) for k in range(power + 1)]
    remainder = list(coefficients)
    if len(remainder) <= power:
        return []
    quotient = [0] * (len(remainder) - power)
    for degree in range(len(remainder) - 1, power - 1, -1):
        value = remainder[degree]
        quotient[degree - power] = value
        for k, coefficient in enumerate(denominator):
            remainder[degree - power + k] -= value * coefficient
    require(all(value == 0 for value in remainder[power:]), "polynomial division failed")
    return quotient


def coefficient_record(coefficients):
    return [[list(key), str(value)] for key, value in sorted(coefficients.items())]


def level(reference):
    q = reference["q"]
    r = q - 2
    weights, margins = convergence(q)
    histogram, prefactor = histogram_from_original_weights(weights, q)
    principal = closed_principal_coefficients(r)
    independent = taylor_principal_coefficients(histogram)
    require(principal == independent, "closed principal coefficients disagree with Taylor products")

    convolution = defaultdict(Fraction)
    for (j, k), first in principal.items():
        for (l, t), second in principal.items():
            convolution[j - l, k + t - 1] += (
                (-1) ** (sum(histogram) + t) * comb(k + t - 2, k - 1) * first * second)
    for (h, u), value in convolution.items():
        require(convolution.get((-h, u), Fraction()) == (-1) ** (u + 1) * value,
                "the complete convolution is not odd")

    M = multiplicity_numerator(r)
    denominator = factorial(2 * r + 2)
    spectral = defaultdict(Fraction)
    polynomial_part = [Fraction()] * len(M)
    for (h, pole), value in convolution.items():
        for u in range(1, pole + 1):
            order = pole - u
            derivative = sum(M[degree] * comb(degree, order) * h ** (degree - order)
                             for degree in range(order, len(M)))
            spectral[h, u] += value * Fraction(derivative, denominator)
        for degree, coefficient in enumerate(polynomial_quotient(M, h, pole)):
            polynomial_part[degree] += value * Fraction(coefficient, denominator)
    require(not any(polynomial_part), "the complete spectral polynomial part did not cancel")
    require(all(not value for (h, u), value in spectral.items() if u > 3),
            "a spectral pole of order greater than three remains")
    for (h, u), value in spectral.items():
        require(spectral.get((-h, u), Fraction()) == (-1) ** (u + 1) * value,
                "the spectral principal parts fail the oddness relation")
    simple = sum((value for (h, u), value in spectral.items() if u == 1), Fraction())
    even = sum((value for (h, u), value in spectral.items() if u == 2), Fraction())
    require(simple == 0, "a harmonic-series coefficient remains")
    require(even == 0, "an even-zeta coefficient remains")

    rational = -prefactor * sum((value * harmonic(r - h, u)
                                 for (h, u), value in spectral.items()), Fraction())
    zeta = prefactor * sum((value for (h, u), value in spectral.items() if u == 3), Fraction())
    require(zeta == Fraction(12 * factorial(r) ** 2, factorial(2 * r + 2)),
            "the uniform zeta(3) coefficient failed")
    binomial_sum = sum((-1) ** (r - k) * comb(r, k) * comb(r + k, k)
                       for k in range(r + 1))
    require(binomial_sum == 1, "the displayed coefficient-extraction identity failed")
    require(set(reference["value"]) == {"0", "3"}, "unexpected reference coefficient space")
    require(rational == Fraction(reference["value"]["0"])
            and zeta == Fraction(reference["value"]["3"]),
            "the independently computed complete value differs from the wheel certificate")
    require(rational == Fraction(-1, 6) if q == 3 else rational > 0,
            "the required finite rational-coefficient sign failed")
    return {
        "q": q, "original_weights": list(weights), "proper_subsets_checked": len(margins),
        "minimum_convergence_margin": min(margins), "occupied_multiplicities": list(histogram),
        "normalization_prefactor": str(prefactor), "contour_strip_largest_node": r,
        "principal_coefficients_checked": len(principal),
        "principal_coefficients_sha256": digest(coefficient_record(principal)),
        "spectral_coefficients_checked": len(spectral),
        "spectral_coefficients_sha256": digest(coefficient_record(spectral)),
        "polynomial_part": "0", "simple_pole_total": "0", "even_zeta_total": "0",
        "rational_coefficient": str(rational), "zeta_3_coefficient": str(zeta),
        "uniform_zeta_coefficient_binomial_sum": binomial_sum,
    }


def uniform_bound():
    trees = [tree for tree in combinations(range(6), 3) if cycle_rank(tree) == 0]
    require(len(trees) == 16, "the K4 Kirchhoff polynomial does not have sixteen terms")
    require(comb(6, 3) == 20 and Fraction(20, 6 ** 3) == Fraction(5, 54),
            "wrong simplex upper bound")
    q = 21
    ratio = Fraction(54, 5) ** q * Fraction(
        factorial(q - 2) * factorial(2 * q - 2), 12 * factorial(3 * q - 1))
    difference = 4 * ratio.numerator - 5 * ratio.denominator
    require(ratio > Fraction(5, 4) and difference > 0, "q=21 lower bound is insufficient")
    raw_difference = [36 - 10, -108 - 45, 72 - 45]
    shifted = [sum(value * comb(degree, k) * 6 ** (degree - k)
                   for degree, value in enumerate(raw_difference) if degree >= k)
               for k in range(3)]
    require(raw_difference == [26, -153, 27] and shifted == [80, 171, 27],
            "the exact all-order growth polynomial failed")
    return {
        "kirchhoff_terms": len(trees), "complete_square_free_triple_terms": comb(6, 3),
        "simplex_upper_bound": "5/54", "tail_start_q": q,
        "lower_ratio_at_21": str(ratio), "four_numerator_minus_five_denominator": str(difference),
        "ratio_difference_polynomial_low_first": raw_difference,
        "ratio_difference_after_q_equals_6_plus_x_low_first": shifted,
        "zeta3_strict_upper_bound": "5/4",
    }


def build():
    reference = json.loads(REFERENCE.read_text())
    rows = reference["finite_k4_base_range"]
    require([row["q"] for row in rows] == list(range(3, 21)), "incomplete finite base range")
    levels = [level(row) for row in rows]
    bound = uniform_bound()
    recorded_bound = reference["uniform_positivity_proof"]
    require(bound["lower_ratio_at_21"] == recorded_bound["lower_ratio_at_21"]
            and bound["four_numerator_minus_five_denominator"]
            == recorded_bound["four_numerator_minus_five_denominator"],
            "the exact uniform bound differs from its recorded value")
    section = SECTION.read_text()
    for label in ("periods:principal-coefficients", "periods:k4-lower-principal",
                  "periods:k4-finite-residues", "periods:k4-rational-algorithm",
                  "periods:ratio-growth", "periods:wheel-span"):
        require("\\label{" + label + "}" in section, "missing manuscript equation " + label)
    paths = (Path(__file__).resolve(), SECTION, REFERENCE)
    return {
        "schema": "methods-paper-independent-wheel-residue-v1",
        "production_evaluator_imported": False,
        "arithmetic": "Integer and fractions.Fraction only; all q=3,...,20 are checked.",
        "method": "Closed principal coefficients versus finite Taylor products, complete convolution, polynomial long division, and finite harmonic corrections.",
        "levels": levels, "uniform_bound": bound,
        "summary": {
            "complete_base_values": len(levels),
            "proper_subset_convergence_checks": sum(row["proper_subsets_checked"] for row in levels),
            "principal_coefficient_checks": sum(row["principal_coefficients_checked"] for row in levels),
            "spectral_coefficient_checks": sum(row["spectral_coefficients_checked"] for row in levels),
            "complete_polynomial_parts_cancelled": len(levels),
        },
        "dependencies_sha256": {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                                for path in paths},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="regenerate the complete exact certificate")
    arguments = parser.parse_args()
    result = build()
    if arguments.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        temporary = OUTPUT.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        temporary.replace(OUTPUT)
        print("Wrote the independently derived complete wheel base certificate.")
    else:
        require(OUTPUT.exists(), "missing wheel_residue.json; regenerate explicitly with --write")
        require(json.loads(OUTPUT.read_text()) == result, "exact wheel replay or source digests changed")
        print("Replayed every independent wheel coefficient and source digest.")
    for name, count in result["summary"].items():
        print(f"{name}: {count}")
    print("Uniform q>=21 bound and its positive growth polynomial: exact.")
    print("The unbounded wheel theorem is proved in periods.tex; finite ranks are not extrapolated.")
    print("canonical_certificate_sha256:", digest(result))


if __name__ == "__main__":
    main()
