"""Compare every scalar coefficient before and after chain-first aggregation.

The default hard five-loop graph is topology 47.  The complete original
source numerator is constructed once, with independent grounded-network
checks.  Each bare/first-current input then follows both scalar reductions:
the historical monomial-by-monomial automorphism map and the new finite
sum by chain powers followed by the same map.  Complete dictionaries of
denominator powers, scalar powers and exact coefficients must be equal.
No finite period library is used by this audit.
"""
import argparse
from fractions import Fraction
import hashlib
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
from certify_integrals import source_assignment
from htk import Graph
from htk.forest_chain_periods import ChainFirstPeriods
from htk.forest_periods import CoreCurrentPeriods
from tree_minor import prepare_source_numerator

IMPLEMENTATIONS = (
    "code/htk/forest_chain_periods.py", "code/htk/forest_periods.py",
    "code/htk/series.py", "code/htk/treeflow.py",
    "explore/circuit-general/symmetry/tree_minor.py",
    "explore/nonSP/verify/k4_periods/current_minor.py",
    "explore/circuit-general/all-product-forests/catalog.py",
    "explore/circuit-general/all-product-forests/certify_integrals.py",
)


def canonical(value):
    return json.loads(json.dumps(value, sort_keys=True))


def scalar_digest(polynomials):
    values = [[power, [[list(alpha), str(coefficient)]
                       for alpha, coefficient in sorted(polynomial.items())]]
              for power, polynomial in sorted(polynomials.items())]
    return hashlib.sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()


def large_rational_regression():
    """Exact source homogeneity beyond Python's integer-string digit limit.

    Scaling every original source by 10**800 multiplies a complete three-loop
    bare kernel by 10**4800.  Direct p/q conversion must preserve every
    scalar coefficient without decimal-string parsing or a digit-limit change.
    """
    graph = Graph(5, [(0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)])
    sources = ((1, 2), (2, -3), (0, 0), (9, 6), (-12, -5))
    scale = 10 ** 800
    scaled_sources = tuple(tuple(scale * entry for entry in row) for row in sources)
    small = ChainFirstPeriods(graph, sources, (3, 4), check=0, max_terms=1000000)
    large = ChainFirstPeriods(graph, scaled_sources, (3, 4), check=0, max_terms=1000000)
    try:
        polynomial = {(0,) * (2 * graph.ne): Fraction(1)}
        original = CoreCurrentPeriods._core_polynomials(small, polynomial)
        scaled = ChainFirstPeriods._core_polynomials(large, polynomial)
        expected = {power: {alpha: coefficient * scale ** (2 * graph.loops)
                            for alpha, coefficient in terms.items()}
                    for power, terms in original.items()}
        if scaled != expected:
            raise ArithmeticError("a large exact rational source coefficient lost source homogeneity")
        bits = max(coefficient.numerator.bit_length() for terms in scaled.values() for coefficient in terms.values())
        if bits <= 15000:
            raise ArithmeticError("the exact rational regression did not exceed the default string-conversion limit")
        return {"source_scale": "10**800", "complete_source_homogeneity_degree": 2 * graph.loops,
                "complete_scalar_coefficients_checked": sum(len(terms) for terms in scaled.values()),
                "maximum_exact_numerator_bits": bits,
                "all_large_rational_coefficients_preserved": True,
                "integer_string_digit_limit_changed": False}
    finally:
        small.close()
        large.close()


def calculate(index=47):
    graph, marked = forest_graph(forests(5)[index])
    sources, superpotential = source_assignment(graph, marked)
    core = ChainFirstPeriods(graph, sources, marked, check=2, max_terms=4000000,
                             max_period_terms=1000000)
    try:
        if core._sp is not None:
            raise ArithmeticError("the independent hard-core audit requires a non-series-parallel graph")
        plan = prepare_source_numerator(core, max_tree_subsets=1000, max_plans=8)
        print("Topology %d complete source numerator: %d monomials." % (index, len(core._bare)), flush=True)
        records = []
        for entry in (None, 0, 1):
            exponent = [0] * (2 * graph.ne)
            if entry is not None:
                exponent[entry] = 1
            polynomial = {tuple(exponent): Fraction(1)}
            started = perf_counter()
            original = CoreCurrentPeriods._core_polynomials(core, polynomial)
            original_statistics = dict(core.last_stats)
            changed = ChainFirstPeriods._core_polynomials(core, polynomial)
            changed_statistics = dict(core.last_stats)
            if original != changed:
                raise ArithmeticError("chain-first aggregation changed a complete scalar coefficient dictionary")
            if core._worker is not None:
                raise ArithmeticError("the scalar dictionary audit unexpectedly queried finite period tables")
            name = "bare" if entry is None else "edge_0_entry_%d" % entry
            records.append({"name": name, "original_current_exponents": exponent,
                            "complete_scalar_dictionary_sha256": scalar_digest(original),
                            "complete_denominator_powers": sorted(original),
                            "complete_scalar_terms": sum(len(value) for value in original.values()),
                            "old_and_new_dictionaries_identical": True,
                            "original_statistics": original_statistics,
                            "chain_first_statistics": changed_statistics})
            print("%s: every one of %d scalar coefficients agrees; %.3f seconds." %
                  (name, records[-1]["complete_scalar_terms"], perf_counter() - started), flush=True)
        return canonical({"schema_version": 1, "index": index, "loops": graph.loops,
            "scope": "Equality of the entire scalar-coefficient dictionaries for the original bare kernel and both original edge-zero first-current entries on a hard five-loop forest. Both routes retain the same complete original source/current numerator; no period evaluation is used.",
            "edges": graph.edges, "marked": marked, "full_balanced_sources": sources,
            "superpotential_vertices": superpotential, "complete_source_minor_plan": plan,
            "bare_and_complete_first_current_spinor_checked": True,
            "all_current_terms_and_sources_retained": True, "finite_period_tables_used": False,
            "source_numerator_monomials": len(core._bare), "records": records,
            "total_scalar_coefficients_compared": sum(record["complete_scalar_terms"] for record in records),
            "large_rational_source_regression": large_rational_regression(),
            "implementation_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                       for name in IMPLEMENTATIONS}})
    finally:
        core.close()


def check_stored(result):
    hashes = result["implementation_sha256"]
    if set(hashes) != set(IMPLEMENTATIONS):
        raise ArithmeticError("the chain-first audit dependency list is incomplete")
    for name in IMPLEMENTATIONS:
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != hashes[name]:
            raise ArithmeticError("the chain-first audit dependency changed: " + name)
    graph, marked = forest_graph(forests(5)[47])
    sources, superpotential = source_assignment(graph, marked)
    expected = canonical({"index": 47, "loops": 5, "edges": graph.edges, "marked": marked,
                          "full_balanced_sources": sources, "superpotential_vertices": superpotential})
    if {name: result[name] for name in expected} != expected:
        raise ArithmeticError("the complete original audit graph or sources changed")
    if [row["name"] for row in result["records"]] != ["bare", "edge_0_entry_0", "edge_0_entry_1"]:
        raise ArithmeticError("the bare kernel or a complete spinor entry is missing")
    total = 0
    for index, row in enumerate(result["records"]):
        exponent = [0] * (2 * graph.ne)
        if index:
            exponent[index - 1] = 1
        count = row["complete_scalar_terms"]
        if (row["original_current_exponents"] != exponent
                or not row["old_and_new_dictionaries_identical"]
                or count != row["original_statistics"]["core_period_orbits"]
                or count != row["chain_first_statistics"]["core_period_orbits"]):
            raise ArithmeticError("a complete scalar-coefficient comparison is inconsistent")
        total += count
    if (total != result["total_scalar_coefficients_compared"]
            or not result["bare_and_complete_first_current_spinor_checked"]
            or not result["all_current_terms_and_sources_retained"]
            or result["finite_period_tables_used"]):
        raise ArithmeticError("the declared complete comparison scope changed")
    if large_rational_regression() != result["large_rational_source_regression"]:
        raise ArithmeticError("the complete large-rational source regression changed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-stored", action="store_true",
                        help="validate the recorded complete scope and dependency hashes, without repeating the hard-core dictionaries")
    options = parser.parse_args()
    path = HERE / "chain_first_audit.json"
    if options.check_stored:
        if options.write:
            parser.error("--check-stored does not regenerate the complete scalar dictionaries")
        check_stored(json.loads(path.read_text()))
        print("Recorded complete chain-first audit, exact dependencies and large-rational regression verified.", flush=True)
        return
    result = calculate()
    if options.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif result != json.loads(path.read_text()):
        raise ArithmeticError("the complete independent chain-first audit changed")
    print("Every scalar coefficient agrees for the bare kernel and complete first-current spinor.", flush=True)


if __name__ == "__main__":
    main()
