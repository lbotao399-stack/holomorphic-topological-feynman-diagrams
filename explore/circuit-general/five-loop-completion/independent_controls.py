"""Independent original-network, ratio-integral and complete-coefficient controls.

The SP controls use the original polynomial ratio Jacobian and simplex
integrals, without scalar period tables.  The non-SP controls compare
complete polynomial/scalar dictionaries before integration and compare
every original first-current coefficient after a different cycle choice.
Source matrices of determinants4,-2,0 test both signs and degenerate frames.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

from bulk3_periods import Bulk3Periods, HERE, ROOT, current_minor, encoded, polynomial_record
from certify import ROWS, atomic_json, hashes
from chain_currents import ChainCurrentPeriods
from current_relations import CurrentRelations
from planned_cases import candidate_trees
from symmetry_controls import frame
from htk import brute, sp_reduce
from htk.graph import Graph
from htk.linalg import det
from htk.treeflow import tree_flow


def decode(value):
    return {key: Fraction(coefficient) for key, coefficient in value.items()}


def case(index):
    return json.loads((HERE / "cases" / ("bulk3_%03d.json" % index)).read_text())


def scalar_record(coefficients):
    rows = [[list(alpha), str(value)] for alpha, value in sorted(coefficients.items())]
    return {"complete_coefficients": len(rows),
            "sha256": hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()}


def sp_controls():
    controls = []
    for index in range(3):
        row, saved = ROWS[index], case(index)
        graph = Graph(row["nv"], row["edges"])
        reduction = sp_reduce(graph, D=3)
        density = reduction.density(frame(row), method="generic", check=2)
        bare = density.integrate()
        currents = reduction.currents(frame(row))
        values = [[density.integrate(value) for value in spinor] for spinor in currents]
        expected = [[decode(value) for value in spinor]
                    for spinor in saved["all_original_first_current_coefficients"]]
        encoded_values = [[({"1": value} if value else {}) for value in spinor] for spinor in values]
        if ({"1": bare} if bare else {}) != decode(saved["bare_coefficient"]) or encoded_values != expected:
            raise ArithmeticError("An independent complete SP ratio integral differs.")
        controls.append({"index": index, "method": "original_ratio_jacobian_and_simplex_integrals",
                         "density_terms": len(density.num), "independent_density_points": 2,
                         "bare_coefficient": str(bare),
                         "all_original_current_coefficients": [[str(v) for v in spinor] for spinor in values]})
    return controls


def source_controls(index, engine):
    row = ROWS[index]
    graph = engine.G
    alternate = next(tree for score, tree in candidate_trees(graph, engine.sources) if tuple(tree) != engine.tree)
    proxy = SimpleNamespace(n=5, N=16, C=graph.cycle_matrix(alternate))
    canonical = {tuple(map(int, p)): int(c) for p, c in engine._bare.terms()}
    matrices = [((2, 1, 0), (0, 1, 2), (1, 0, 1)),
                ((1, 1, 0), (1, 0, 1), (0, 1, 1)),
                ((1, 2, 0), (2, 4, 0), (0, 0, 1))]
    controls = []
    point = [Fraction(e + 2, e + 3) for e in range(16)]
    original = brute.currents(graph, point, engine.sources)
    for matrix in matrices:
        sources = [[0, 0, 0] for _ in range(graph.nv)]
        for vertex, value in zip(row["S0"][:3], matrix):
            sources[vertex] = list(value)
        sources[row["S0"][-1]] = [-sum(values) for values in zip(*matrix)]
        determinant = det(matrix)
        flow = [[int(x) for x in values] for values in tree_flow(graph, sources, alternate)]
        numerator = current_minor.numerator(proxy, sources, flow,
                                            determinant=current_minor.polynomial_det)
        observed = {tuple(map(int, p)): int(c) for p, c in numerator.terms()}
        expected = {p: int(determinant ** 5) * c for p, c in canonical.items() if determinant ** 5 * c}
        if observed != expected:
            raise ArithmeticError("The full source-transformed polynomial differs coefficientwise.")
        changed = brute.currents(graph, point, sources)
        predicted = [[sum(original[e][b] * matrix[b][a] for b in range(3)) for a in range(3)]
                     for e in range(16)]
        if changed != predicted:
            raise ArithmeticError("A full original current source direction differs.")
        if brute.density(graph, point, sources).rational() != determinant ** 5 * brute.density(
                graph, point, engine.sources).rational():
            raise ArithmeticError("The independent original density source law failed.")
        controls.append({"source_matrix": [list(values) for values in matrix],
                         "source_determinant": str(determinant),
                         "complete_transformed_numerator": polynomial_record(numerator),
                         "all_original_transformed_currents": [[str(v) for v in spinor] for spinor in changed]})
    return {"alternate_source_cycle_tree": list(alternate), "complete_source_transformations": controls}


def nonsp_controls():
    controls = []
    for index in (3, 46):
        saved, engine = case(index), Bulk3Periods(ROWS[index])
        try:
            engine.tree_polynomials()
            pushed = ChainCurrentPeriods(engine)
            coefficient_checks = []
            edges = range(16) if index == 3 else engine.G.cotree(engine.tree)
            for edge in edges:
                for component in range(3):
                    original = engine.scalar_coefficients(engine._bare * engine._currents[edge][component])
                    after_fibres = pushed.coefficients(pushed.numerator(edge, component))
                    if original != after_fibres:
                        raise ArithmeticError("The entire first-current scalar coefficient dictionary differs.")
                    coefficient_checks.append({"edge": edge, "source_direction": component,
                                               **scalar_record(original)})
            transformed = source_controls(index, engine)
            alternate = Bulk3Periods(ROWS[index], tree=transformed["alternate_source_cycle_tree"])
            try:
                if alternate._bare != engine._bare:
                    raise ArithmeticError("The alternative cycle tree changed the complete original polynomial.")
                bare = alternate.integrate()
                if encoded(bare) != saved["bare_coefficient"]:
                    raise ArithmeticError("The complete alternative-tree bare value changed.")
                if index == 3:
                    values = alternate.complete_first_currents(bare, direct=True)
                    method = "all_original_edges_integrated_directly_with_alternate_cycle_tree"
                else:
                    direct = ChainCurrentPeriods(alternate)
                    values = [[direct.integrate(edge, component) for component in range(3)] for edge in range(16)]
                    method = "all_original_edges_integrated_after_exact_fibres_with_alternate_cycle_tree"
                encoded_values = [[encoded(value) for value in spinor] for spinor in values]
                if encoded_values != saved["all_original_first_current_coefficients"]:
                    raise ArithmeticError("A complete original first-current value differs from the census.")
                relations = CurrentRelations(ROWS[index])
                relations.check(bare, values)
                flat = [value for spinor in values for value in spinor]
                if relations.reconstruct(bare, {j: flat[j] for j in relations.free}) != values:
                    raise ArithmeticError("The exact symmetry reconstruction changed the full tensor.")
                controls.append({"index": index, "complete_original_numerator": polynomial_record(engine._bare),
                                 "complete_scalar_dictionary_comparisons": coefficient_checks,
                                 "all_original_current_integral_method": method,
                                 "bare_coefficient": encoded(bare),
                                 "all_original_first_current_coefficients": encoded_values,
                                 "complete_current_relation_rows": len(relations.rows),
                                 "independent_current_entries": len(relations.free), **transformed})
            finally:
                alternate.close()
        finally:
            engine.close()
        print("Recomputed all original source directions and current vectors for", index, flush=True)
    return controls


def build():
    sp, nonsp = sp_controls(), nonsp_controls()
    bindings = hashes()
    extras = [Path(__file__), HERE / "chain_currents.py", HERE / "current_relations.py",
              HERE / "planned_cases.py", HERE / "symmetry_controls.py",
              ROOT / "code/htk/spred.py", ROOT / "code/htk/simplex.py", ROOT / "code/htk/__init__.py"]
    extras += [HERE / "cases" / ("bulk3_%03d.json" % index) for index in (0, 1, 2, 3, 46)]
    bindings.update({str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in extras})
    return {"schema": "five-loop-bulk3-independent-controls-v1",
            "independent_sp_ratio_controls": sp, "complete_nonsp_controls": nonsp,
            "summary": {"full_current_tensors_reintegrated": 5,
                        "all_original_current_coefficients_reintegrated": 240,
                        "complete_scalar_dictionaries_compared": sum(len(x["complete_scalar_dictionary_comparisons"]) for x in nonsp),
                        "all_scalar_coefficients_compared": sum(y["complete_coefficients"] for x in nonsp
                                                               for y in x["complete_scalar_dictionary_comparisons"]),
                        "complete_nonunit_or_degenerate_source_polynomials": 6,
                        "full_original_current_source_transformation_checks": 288},
            "dependencies_sha256": bindings}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build()
    path = HERE / "independent_certificate.json"
    if args.write:
        atomic_json(path, result)
    elif result != json.loads(path.read_text()):
        raise ArithmeticError("An independently recomputed complete control or source digest changed.")
    print(json.dumps(result["summary"], indent=2))
