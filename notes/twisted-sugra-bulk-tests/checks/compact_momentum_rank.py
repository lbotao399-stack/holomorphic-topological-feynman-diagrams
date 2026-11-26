#!/usr/bin/env python3
"""Exact routing and sharp compact-momentum support tests for the 11d theta kernel.

The certificate concerns the geometric kernel. It does not infer a nonzero
complete supergravity anomaly from one geometric graph.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from random import Random
import json

import sympy as sp


def basis(columns, dimension):
    if not columns:
        return sp.zeros(dimension, 0)
    vectors = sp.Matrix.hstack(*columns).columnspace()
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(dimension, 0)


def intersection(left, right):
    """Compute intersection by solving equal-image equations over Q."""
    dimension = left.rows
    if left.cols == 0 or right.cols == 0:
        return sp.zeros(dimension, 0)
    relations = left.row_join(-right).nullspace()
    vectors = [left * relation[:left.cols, 0] for relation in relations]
    return basis(vectors, dimension)


def compact_matrix(paths, dimension):
    zero = sp.zeros(dimension, 1)
    columns = []
    for path, momenta in enumerate(paths):
        for p in momenta:
            if path == 0:
                columns.append(p.col_join(zero))
            elif path == 1:
                columns.append(p.col_join(p))
            else:
                columns.append(zero.col_join(p))
    return sp.Matrix.hstack(*columns) if columns else sp.zeros(2*dimension, 0)


def rank_certificate():
    rng = Random(20261126)
    rows = []
    for dimension in range(1, 6):
        cases = []
        unit = [sp.eye(dimension)[:, i] for i in range(dimension)]
        cases.extend([
            ([], [], []),
            (unit, [], []),
            (unit, [-v for v in unit], []),
            (unit, unit, [-2*v for v in unit]),
        ])
        for sample in range(20):
            paths = [[], [], []]
            for path in range(3):
                for _ in range(rng.randrange(dimension + 2)):
                    paths[path].append(
                        sp.Matrix([rng.randrange(-3, 4) for _ in range(dimension)])
                    )
            total = sum((v for path in paths for v in path), sp.zeros(dimension, 1))
            paths[sample % 3].append(-total)
            cases.append(tuple(paths))
        for index, paths in enumerate(cases):
            spans = [basis(path, dimension) for path in paths]
            common = intersection(intersection(spans[0], spans[1]), spans[2])
            expected = sum(span.cols for span in spans) - common.cols
            actual = compact_matrix(paths, dimension).rank()
            assert actual == expected
            nonzero_slots = sum(any(x != 0 for x in v) for path in paths for v in path)
            assert actual <= min(2*dimension, nonzero_slots)
            rows.append({
                "compact_dimension": dimension,
                "case": index,
                "path_span_dimensions": [span.cols for span in spans],
                "triple_intersection_dimension": common.cols,
                "compact_rank": actual,
                "nonzero_compact_slots": nonzero_slots,
            })
    return rows


def sharp_support_certificate():
    # Shape (1,6,6): five columns on B and five on C.
    units = [sp.eye(5)[:, i] for i in range(5)]
    paths = ([], units, [-v for v in units])
    determinant = compact_matrix(paths, 5).det()
    assert determinant == -1
    assert sum((v for path in paths for v in path), sp.zeros(5, 1)) == sp.zeros(5, 1)
    rows = []
    for compact_dimension in range(1, 6):
        compact_paths = [
            [v[:compact_dimension, 0] for v in path]
            for path in paths
        ]
        number = sum(any(x != 0 for x in v) for path in compact_paths for v in path)
        rank = compact_matrix(compact_paths, compact_dimension).rank()
        assert number == rank == 2*compact_dimension
        rows.append({"compact_dimension": compact_dimension,
                     "compact_external_slots": number,
                     "compact_rank": rank,
                     "full_holomorphic_determinant": int(determinant)})
    scalar = sp.simplify(
        sp.Rational(1, 2) * sp.gamma(sp.Rational(11, 2))
        * sp.gamma(sp.Rational(1, 2))**2
        / (sp.gamma(sp.Rational(13, 2)) * sp.factorial(5)**2)
    )
    assert scalar == sp.pi / 158400
    return {"cases": rows, "positive_scalar_integral_over_pi": "1/158400",
            "oriented_kernel_over_pi": str(Fraction(int(determinant), 158400))}


def main():
    ranks = rank_certificate()
    result = {
        "status": "PASS",
        "arithmetic": "exact rational linear algebra",
        "rank_formula_checks": len(ranks),
        "rank_formula_cases": ranks,
        "sharp_support": sharp_support_certificate(),
        "scope": "Geometric theta kernel, before supergravity color and graph sums.",
    }
    destination = Path(__file__).resolve().parents[1] / "verification/compact_momentum_rank.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "rank_formula_checks": len(ranks),
                      "sharp_support_cases": 5,
                      "oriented_theta_166_kernel_over_pi": "-1/158400"}))


if __name__ == "__main__":
    main()
