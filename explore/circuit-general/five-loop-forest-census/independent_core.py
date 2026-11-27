"""Independent original-graph control for the sourced five-loop forest 42.

Construct a complete original source numerator with a source-flow tree and
cycle basis different from the census's selected plan.  Integrate the bare
kernel and both original edge-zero current entries directly with
CoreCurrentPeriods.  No ForestReduction call, passive contraction, sourced
ear operation, or cotree reconstruction enters these three integrals.

Finite scalar-period tables are shared with the census; the independent
parts are the original source/cycle plan, the full numerator, original
current integration, and exact grounded-network pointwise controls.
"""
import argparse
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(HERE.parent / "all-product-forests"))
sys.path.insert(0, str(HERE.parent / "symmetry"))

from catalog import forest_graph, forests
from certify_integrals import finite_record, source_assignment
from htk.forest_periods import CoreCurrentPeriods
from tree_minor import _is_tree, prepare_source_numerator


CENSUS_TREE = (1, 2, 4, 5, 8, 10)
IMPLEMENTATIONS = (
    "code/htk/forest_periods.py", "code/htk/treeflow.py",
    "explore/nonSP/verify/k4_periods/current_minor.py",
    "explore/circuit-general/symmetry/tree_minor.py",
    "explore/circuit-general/all-product-forests/catalog.py",
    "explore/circuit-general/all-product-forests/certify_integrals.py",
)


def canonical(value):
    return json.loads(json.dumps(value, sort_keys=True))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def polynomial_record(polynomial):
    terms = [[list(map(int, exponent)), str(coefficient)]
             for exponent, coefficient in sorted(polynomial.terms())]
    return {"terms": len(terms), "total_degrees": sorted({sum(exponent) for exponent, coefficient in terms}),
            "all_coefficients_sha256": digest(terms)}


def calculate():
    graph, marked = forest_graph(forests(5)[42])
    sources, superpotential = source_assignment(graph, marked)
    candidates = [tree for tree in combinations(range(graph.ne), graph.nv - 1)
                  if tree != CENSUS_TREE and _is_tree(graph, tree)]
    core = CoreCurrentPeriods(graph, sources, marked, check=3, max_terms=1000000,
                              max_period_terms=1000000, period_timeout=600)
    try:
        if core._sp is not None:
            raise ArithmeticError("the independent non-series-parallel control changed topology")
        start = perf_counter()
        plan = prepare_source_numerator(core, tree_candidates=candidates,
                                        max_tree_subsets=1000, max_plans=1)
        if tuple(plan["chosen_tree"]) == CENSUS_TREE:
            raise ArithmeticError("the independent source/cycle plan was not changed")
        print("Independent complete source numerator: %d terms; tree %s; %.3f seconds." %
              (len(core._bare), plan["chosen_tree"], perf_counter() - start), flush=True)
        source_numerator = polynomial_record(core._bare)
        # Compare every coefficient with a fresh construction in the census
        # basis.  This is a polynomial equality, beyond pointwise sampling.
        reference = CoreCurrentPeriods(graph, sources, marked, check=0, max_terms=1000000)
        try:
            reference_plan = prepare_source_numerator(reference, tree_candidates=[CENSUS_TREE],
                                                       max_tree_subsets=1, max_plans=1)
            if core._bare != reference._bare:
                raise ArithmeticError("the full source numerators in the two cycle bases differ")
        finally:
            reference.close()
        print("All %d source-numerator coefficients agree in the two original cycle bases." % len(core._bare),
              flush=True)
        core._tree_polynomials()
        # Coefficientwise original Kirchhoff laws, with the original source
        # at every vertex and both entries of every source spinor retained.
        for vertex in range(graph.nv):
            for entry in range(2):
                value = sum(((int(tail == vertex) - int(head == vertex)) * core._currents[2 * edge + entry]
                             for edge, (tail, head) in enumerate(graph.edges)), core._context.constant(0))
                if value != sources[vertex][entry] * core._U:
                    raise ArithmeticError("a full original current polynomial violates source conservation")
        current_numerators = [polynomial_record(value) for value in core._currents]
        values, statistics = [], []
        for entry in (None, 0, 1):
            exponent = [0] * (2 * graph.ne)
            if entry is not None:
                exponent[entry] = 1
            value = core.integrate({tuple(exponent): Fraction(1)}, allow_series=False)
            values.append(finite_record(value))
            statistics.append(dict(core.last_stats))
            print("Complete original-graph %s integral: %s" %
                  ("bare" if entry is None else "edge 0 entry %d" % entry, values[-1]), flush=True)
        return canonical({"schema_version": 1, "loops": 5, "index": 42,
            "scope": "Independent complete original-graph bare kernel and first-current spinor on edge 0, using an alternate complete source/cycle plan. Finite scalar-period tables are shared with the census.",
            "edges": graph.edges, "marked": marked, "full_balanced_sources": sources,
            "superpotential_vertices": superpotential,
            "census_source_flow_tree": CENSUS_TREE, "independent_source_minor_plan": plan,
            "census_basis_source_minor_plan": reference_plan,
            "full_source_numerator": source_numerator,
            "full_source_numerator_coefficients_compared": len(core._bare),
            "all_original_current_numerators": current_numerators,
            "coefficientwise_original_source_conservation_checks": 2 * graph.nv,
            "independent_original_network_density_points": core.check,
            "independent_original_network_current_entries_checked": core.check * 2 * graph.ne,
            "bare": values[0], "edge": 0, "full_first_current_spinor": values[1:],
            "integration_statistics": statistics,
            "ForestReduction_used": False, "forest_precontractions_used": False,
            "current_reconstruction_used_for_integrated_spinor": False,
            "all_values_finite": True, "residual_periods": [],
            "implementation_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                       for name in IMPLEMENTATIONS}})
    finally:
        core.close()


def compare_census(result, census_path):
    hashes = result["implementation_sha256"]
    if set(hashes) != set(IMPLEMENTATIONS):
        raise ArithmeticError("the independent control dependency list is incomplete")
    for name in IMPLEMENTATIONS:
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != hashes[name]:
            raise ArithmeticError("the independent control dependency changed: " + name)
    document = json.loads(census_path.read_text())
    record = next(row for row in document["records"] if row["index"] == 42) if "records" in document else document
    for field in ("loops", "index", "edges", "marked", "full_balanced_sources", "superpotential_vertices", "bare"):
        if result[field] != record[field]:
            raise ArithmeticError("independent original-graph control differs from census field " + field)
    if result["full_first_current_spinor"] != record["all_original_first_current_spinors"][0]:
        raise ArithmeticError("the independently integrated complete original current spinor differs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--compare-only", action="store_true",
                        help="compare the completed independent certificate with the separately completed census")
    parser.add_argument("--census", type=Path, default=HERE / "certificate.json")
    options = parser.parse_args()
    path = HERE / "independent_core_certificate.json"
    if options.compare_only:
        if options.write:
            parser.error("--compare-only does not regenerate the independent calculation")
        result = json.loads(path.read_text())
        compare_census(result, options.census)
        print("Every recorded independent source and all three complete kernels agree with the census.", flush=True)
        return
    result = calculate()
    if options.write:
        # All three complete independent integrals have succeeded. Preserve
        # that result even if the separately running census is not ready.
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif result != json.loads(path.read_text()):
        raise ArithmeticError("the independent complete original-graph certificate changed")
    compare_census(result, options.census)
    print("Independent bare kernel and complete first-current spinor agree with the census exactly.", flush=True)


if __name__ == "__main__":
    main()
