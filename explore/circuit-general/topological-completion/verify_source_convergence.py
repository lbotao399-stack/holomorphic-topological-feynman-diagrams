"""Absolute convergence before any orbit cancellation, at every source.

Full-degree source evaluations are unisolvent after graph relabelings.
Each original numerator monomial is convergent separately. Every current
numerator uses only Kirchhoff cotree monomials, so each original current
is bounded throughout the positive projective domain at fixed sources.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import sys

from flint import fmpq_mat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT/"code"), str(ROOT/"explore/circuit-general/cores")]
from htk.graph import Graph
from htk.topological_polynomial import BulkPolynomial
import cores
import gutil
from evaluate_sources import polynomial_record
from source_covariance import automorphisms, homogeneous


def complete_ranks(edges, points):
    relabeled = set()
    for point in points:
        full = list(point) + [-sum(point)]
        for permutation, *_ in automorphisms(tuple(map(tuple, edges))):
            relabeled.add(tuple(full[permutation.index(i)] for i in range(5)))
    result = {}
    for degree in (1, 4, 5):
        monomials = list(homogeneous(5, degree))
        matrix = fmpq_mat([[prod(x**a for x, a in zip(point, exponent))
                           for exponent in monomials] for point in sorted(relabeled)])
        rank = matrix.rank()
        if degree in (1, 4):
            assert rank == len(monomials), (degree, rank, len(monomials))
        result[str(degree)] = {"exact_rank": rank, "all_homogeneous_monomials": len(monomials)}
    return {"distinct_relabelled_source_points": len(relabeled), "ranks": result}


def main():
    records = []
    for name in ("prism", "K33"):
        data = json.loads((HERE/(name+"_source_covariance.json")).read_text())
        points = [tuple(row["sources"]) for row in data["samples"]]
        if name == "K33":
            points.extend([(1, 2, 3, 4, 5), (2, -3, 5, 7, -13)])
        graph = Graph(6, data["edges"])
        _, ranks = gutil.subset_table(tuple(graph.edges))
        convergence = {}
        samples = []
        for index, point in enumerate(points):
            method = "gaussian_minor" if name == "prism" else "pfaffian"
            density = BulkPolynomial(graph, [[x] for x in point], H=1, T=1, method=method)
            density.check_at((2, 3, 5, 7, 11, 13, 17, 19, 23))
            record = polynomial_record(density.numerator)
            if index < len(data["samples"]):
                saved = json.loads((HERE/"samples"/(name+"_%02d.json" % index)).read_text())
                assert record == saved["density"]
            cotrees = set(density.U.to_dict())
            assert all(coefficient == 1 for coefficient in density.U.to_dict().values())
            assert all(set(row[0].to_dict()) <= cotrees for row in density.current_numerators)
            smallest = None
            for exponent in density.numerator.to_dict():
                alpha = tuple(int(x)+1 for x in exponent)
                key = (alpha, density.power)
                if key not in convergence:
                    assert sum(alpha) == graph.loops*density.power
                    doubled_q = 2*density.power
                    assert doubled_q.denominator == 1
                    sums = [0]*(1 << graph.ne)
                    margin = None
                    for mask in range(1, (1 << graph.ne)-1):
                        low = mask & -mask
                        sums[mask] = sums[mask ^ low] + alpha[low.bit_length()-1]
                        difference = 2*sums[mask] - int(doubled_q)*ranks[mask]
                        assert difference > 0, (name, point, alpha, mask, difference)
                        margin = difference if margin is None else min(margin, difference)
                    convergence[key] = Fraction(margin, 2)
                margin = convergence[key]
                smallest = margin if smallest is None else min(smallest, margin)
            samples.append({"sources": point, "bare_numerator": record, "pole": str(density.power),
                            "smallest_strict_cut_margin": str(smallest),
                            "all_nine_original_currents_bounded_by_cotree_weights": True})
            print(name, index, "all original numerator monomials converge; all currents bounded", flush=True)
        record = {"core": name, "samples": samples,
                  "distinct_original_scalar_monomials_checked": len(convergence),
                  "every_proper_edge_subset_checked": (1 << graph.ne)-2,
                  "all_source_coefficient_convergence": complete_ranks(graph.edges, points)}
        records.append(record)
        print(name, record["all_source_coefficient_convergence"], flush=True)
    report = {"passed": True, "cores": records,
              "inputs_sha256": {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                  for path in [Path(__file__), ROOT/"code/htk/topological_polynomial.py",
                               HERE/"source_covariance.py"]}}
    (HERE/"source_convergence.json").write_text(json.dumps(report, indent=2)+"\n")
    print("PASS: arbitrary-source coefficient convergence and every finite holomorphic-current word", flush=True)


if __name__ == "__main__":
    main()
