#!/usr/bin/env python3
"""Exact checks for the all-loop wheel stabilization theorem.

The proof is in README.md.  The finite checks here certify its finite base
range, normalization, polynomial inequality, and independent graph reductions.
No empirical rank plateau is used as a proof of the all-loop statement.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from math import comb, factorial
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def module(relative, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


WHEEL = module("explore/circuit-general/cores/wheel.py", "closure_wheel")
GRAPH = module("explore/circuit-general/cores/gutil.py", "closure_graph")


def encoded(value):
    return {str(weight): str(coefficient)
            for weight, coefficient in sorted(value.items()) if coefficient}


def wheel_edges(n):
    # Ordered orientations used in README.md.  Scalar periods are independent
    # of these orientations; retain them explicitly for reproducibility.
    return ([(0, i) for i in range(1, n + 1)]
            + [(i, i + 1) for i in range(1, n)] + [(n, 1)])


def path_wheel(n, q, special=False):
    spoke = [q - 1] + [q - 2] * (n - 2) + [q - 1]
    if special:
        spoke[1] = q - 1
    rim = [1] * (n - 1) + [q - 1]
    assert sum(spoke) + sum(rim) == n * q
    edges = wheel_edges(n)
    graph = {tuple(sorted(edge)): power
             for edge, power in zip(edges, spoke + rim)}
    assert len(graph) == 2 * n
    return graph, spoke, rim


def parallel_factor(a, b, q):
    """Unnormalized I_q factor B(q-a,q-b); new edge power is a+b-q."""
    assert 0 < a < q and 0 < b < q and a + b > q
    return Fraction(factorial(q - a - 1) * factorial(q - b - 1),
                    factorial(2 * q - a - b - 1))


def remove_unique_vertex(graph, vertex, q):
    """Star--triangle followed by every newly formed parallel pair.

    The star has sum of Schwinger powers q.  Replace its edge opposite
    a neighbor pair by q minus that edge's power.  Star--triangle has no
    scalar prefactor in the I_q convention used here.
    """
    incident = {v if u == vertex else u: power
                for (u, v), power in graph.items() if vertex in (u, v)}
    assert len(incident) == 3 and sum(incident.values()) == q
    neighbors = sorted(incident)
    result = {edge: power for edge, power in graph.items()
              if vertex not in edge}
    factor = Fraction(1)
    merges = []
    for i, u in enumerate(neighbors):
        for v in neighbors[i + 1:]:
            opposite = next(w for w in neighbors if w not in (u, v))
            new_power = q - incident[opposite]
            assert 0 < new_power < q
            edge = (u, v)
            if edge in result:
                old_power = result[edge]
                beta = parallel_factor(old_power, new_power, q)
                factor *= beta
                result[edge] = old_power + new_power - q
                merges.append({"edge": list(edge), "powers": [old_power, new_power],
                               "new_power": result[edge], "factor": str(beta)})
            else:
                result[edge] = new_power
    return result, factor, {"vertex": vertex,
                            "star": [[v, incident[v]] for v in neighbors],
                            "parallel_merges": merges}


def rational_witness(n):
    q = n
    graph, spoke, rim = path_wheel(n, q)
    factor, steps = Fraction(1), []
    for vertex in range(2, n):
        graph, beta, step = remove_unique_vertex(graph, vertex, q)
        factor *= beta
        steps.append(step)
    assert graph == {(0, 1): 1, (0, n): 1, (1, n): q - 2}
    assert factor == Fraction(1, factorial(q - 2) ** 3)
    triangle = Fraction(factorial(q - 3), factorial(q - 1))
    value = factor * triangle
    expected = Fraction(1, (q - 1) * (q - 2) * factorial(q - 2) ** 3)
    assert value == expected
    direct = WHEEL.integrate(spoke, rim, q)
    assert direct == {0: expected}
    if n <= 7:
        assert GRAPH.convergent(tuple(wheel_edges(n)), spoke + rim, q)
    return {"loops": n, "q": q, "spokes": spoke, "rims": rim,
            "value": str(value), "unique_vertex_steps": steps,
            "parallel_factor": str(factor)}


def k4_value(q):
    spoke, rim = [q - 1, 1, 1], [1, q - 1, q - 1]
    assert GRAPH.convergent(tuple(wheel_edges(3)), spoke + rim, q)
    value = WHEEL.integrate(spoke, rim, q)
    coefficient = Fraction(12 * factorial(q - 2) ** 2, factorial(2 * q - 2))
    assert set(value) == {0, 3}
    assert value[3] == coefficient
    assert (value[0] == Fraction(-1, 6) if q == 3 else value[0] > 0)
    histogram, prefactor = WHEEL.indices(spoke, rim, q)
    assert histogram == (2,) * (q - 2) + (3,)
    assert prefactor == factorial(q - 2) ** 4
    r = q - 2
    binomial_sum = sum((-1) ** (r - k) * comb(r, k) * comb(r + k, k)
                       for k in range(r + 1))
    assert binomial_sum == 1
    return value


def k4_embedding(n):
    assert n >= 4
    q = n - 1
    graph, spoke, rim = path_wheel(n, q, special=True)
    factor, steps = Fraction(1), []
    for vertex in range(3, n):
        graph, beta, step = remove_unique_vertex(graph, vertex, q)
        factor *= beta
        steps.append(step)
    assert graph == {(0, 1): q - 1, (0, 2): 1, (0, n): 1,
                     (1, 2): 1, (2, n): q - 1, (1, n): q - 1}
    assert factor == Fraction(1, factorial(q - 2) ** 3)
    base = k4_value(q)
    direct = WHEEL.integrate(spoke, rim, q)
    assert direct == {w: factor * c for w, c in base.items()}
    if n <= 7:
        assert GRAPH.convergent(tuple(wheel_edges(n)), spoke + rim, q)
    return {"loops": n, "q": q, "spokes": spoke, "rims": rim,
            "value": encoded(direct), "unique_vertex_steps": steps,
            "parallel_factor": str(factor)}


def determinant(rows):
    rows = [[Fraction(entry) for entry in row] for row in rows]
    result = Fraction(1)
    for column in range(len(rows)):
        pivot = next((i for i in range(column, len(rows)) if rows[i][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[pivot], rows[column] = rows[column], rows[pivot]
            result *= -1
        entry = rows[column][column]
        result *= entry
        for i in range(column + 1, len(rows)):
            scale = rows[i][column] / entry
            rows[i] = [a - scale * b for a, b in zip(rows[i], rows[column])]
    return result


def spanning_control(n):
    threshold = max(3, n - 1)
    values = []
    if n == 3:
        values.append({0: Fraction(1, 2)})
    else:
        values.append({int(w): Fraction(c) for w, c in k4_embedding(n)["value"].items()})
    for q in range(2, n):
        value = WHEEL.integrate([q - 1] * n, [1] * n, q)
        assert value == WHEEL.radial_master(n, q)
        lowest_weight = 2 * n - 2 * q + 1
        assert min(value) == lowest_weight
        assert value[lowest_weight] == Fraction(2 * comb(2 * n - 2, n - 1),
                                                 factorial(2 * q - 2))
        values.append(value)
    weights = [0] + list(range(3, 2 * n - 2, 2))
    full_det = determinant([[value.get(w, 0) for w in weights] for value in values])
    assert full_det != 0
    return {"loops": n, "threshold": threshold, "weights": weights,
            "spanning_values": [encoded(value) for value in values],
            "coefficient_determinant": str(full_det),
            "lower_level_radial_values_lift_by": "multiplication by U^(threshold-q)"}


def positivity_bound():
    q = 21
    lower_ratio = (Fraction(54, 5) ** q
                   * Fraction(factorial(q - 2) * factorial(2 * q - 2),
                              12 * factorial(3 * q - 1)))
    difference = 4 * lower_ratio.numerator - 5 * lower_ratio.denominator
    assert difference > 0
    # Exact polynomial identity after setting q=6+x:
    # 36(q-1)(2q-1)-5(3q+1)(3q+2)=27x^2+171x+80.
    left = [36 * 5 * 11, 36 * (5 * 2 + 11), 36 * 2]
    right = [5 * 19 * 20, 5 * (19 * 3 + 20 * 3), 5 * 9]
    shifted = [a - b for a, b in zip(left, right)]
    assert shifted == [80, 171, 27] and all(c > 0 for c in shifted)
    return {"kirchhoff_simplex_upper_bound": "5/54",
            "zeta3_strict_upper_bound": "5/4",
            "all_q_tail_starts": q, "lower_ratio_at_21": str(lower_ratio),
            "four_numerator_minus_five_denominator": str(difference),
            "ratio_monotonicity_after_q_6_coefficients_low_first": shifted}


def build():
    finite = [{"q": q, "value": encoded(k4_value(q))}
              for q in range(3, 21)]
    rational = [rational_witness(n) for n in range(3, 13)]
    embedded = [k4_embedding(n) for n in range(4, 13)]
    spans = [spanning_control(n) for n in range(3, 13)]
    dependencies = [Path(__file__), HERE / "README.md",
                    ROOT / "explore/circuit-general/cores/wheel.py",
                    ROOT / "explore/circuit-general/cores/gutil.py"]
    return {"schema": "wheel_scalar_closure_v1",
            "theorem": "For every n>=3 and integer q>=max(3,n-1), the full convergent wheel scalar span is Q+sum_{j=2}^{n-1} Q*zeta(2j-1).",
            "normalization": "I_q=int_Omega prod(t_e^(alpha_e-1))/U^q; no Gamma normalization",
            "orientations": "spokes 0->i (i=1..n), rims i->i+1 (i=1..n-1), final rim n->1",
            "arithmetic_independence_assumed": False,
            "finite_k4_base_range": finite,
            "uniform_positivity_proof": positivity_bound(),
            "rational_witness_controls": rational,
            "k4_embedding_controls": embedded,
            "spanning_controls": spans,
            "dependencies_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                    for path in dependencies}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="replace the exact certificate")
    args = parser.parse_args()
    certificate = build()
    path = HERE / "certificate.json"
    if args.write:
        path.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
        print("Wrote", path)
    else:
        stored = json.loads(path.read_text())
        if stored != certificate:
            raise AssertionError("stored certificate differs from complete exact replay")
    print("PASS: 18 K4 finite bases, the uniform q>=21 inequality, 10 rational wheel controls,")
    print("      9 independent star--triangle embeddings, and 10 complete spanning determinants.")
    print("The all-n and all-q proof is in README.md; finite ranks are controls, not an extrapolation.")


if __name__ == "__main__":
    main()
