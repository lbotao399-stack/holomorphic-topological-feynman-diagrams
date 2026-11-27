"""Fresh original-graph integrations, with all edges integrated separately.

These controls share the exact numerator formula and finite period tables,
but do not use the interpolation or cotree reconstruction being checked.
The source/cycle tree differs from every interpolation construction.
"""
import argparse
import json
from time import perf_counter

from source_currents import (HERE, OPTIONS, graph, calculate_sources, digest,
                            evaluate, normalized_sources, polynomial_record)
from certify import atomic_json, dependency_hashes

ALTERNATE_TREE = (0, 2, 4, 6, 8, 10)
FRESH_SOURCES = ((0, 0), (2, 1), (1, 2), (0, 0), (0, 0), (-1, 4), (-2, -7))
DEGENERATE_SOURCES = ((0, 0), (1, 0), (2, 0), (0, 0), (0, 0), (0, 1), (-3, -1))
CASES = (("fresh_nonunit_frame", FRESH_SOURCES),
         ("degenerate_frame", DEGENERATE_SOURCES))


def compare_original_bases(sources, expected_record):
    from htk.forest_periods import CoreCurrentPeriods
    from tree_minor import prepare_source_numerator
    core = CoreCurrentPeriods(graph(), sources, **OPTIONS)
    try:
        plan = prepare_source_numerator(core, tree_candidates=[graph().spanning_tree()],
                                        max_tree_subsets=1, max_plans=1)
        observed = polynomial_record(core._bare)
        if observed != expected_record["original_numerator"]:
            raise ArithmeticError("a complete original numerator differs across source/cycle bases")
        return {"all_original_coefficients_compared": len(core._bare),
                "full_numerator": observed, "comparison_source_minor_plan": plan,
                "whole_coefficient_dictionary_equality": True}
    finally:
        core.close()


def compute(*, force=False):
    for name, sources in CASES:
        path = HERE / (name + ".json")
        hashes = dependency_hashes()
        if path.exists() and not force:
            existing = json.loads(path.read_text())
            if existing["dependency_sha256"] != hashes:
                raise ArithmeticError("an independent control has stale dependencies")
            print("Reused complete independent control", name, flush=True)
            continue
        started = perf_counter()
        def progress(edge, value):
            print("Independent %s original edge %d: complete spinor, %.2f seconds." %
                  (name, edge, perf_counter() - started), flush=True)
        result = calculate_sources(sources, direct=True, tree=ALTERNATE_TREE, progress=progress)
        basis_check = compare_original_bases(sources, result)
        atomic_json(path, {"schema_version": 1, "name": name, "result": result,
                          "independent_source_basis_comparison": basis_check,
                          "dependency_sha256": hashes,
                          "interpolation_used": False, "cotree_reconstruction_used": False})
        print("Completed independent full current control", name, flush=True)


def check(*, write=False):
    path = HERE / "certificate.json"
    document = json.loads(path.read_text())
    controls = {}
    for name, sources in CASES:
        control_path = HERE / (name + ".json")
        control = json.loads(control_path.read_text())
        if control["dependency_sha256"] != dependency_hashes():
            raise ArithmeticError("an independent control dependency changed")
        actual = control["result"]
        expected_sources = [[str(value) for value in row] for row in normalized_sources(sources)]
        if actual["full_balanced_sources"] != expected_sources:
            raise ArithmeticError("the independent control source assignment changed")
        if actual["source_minor_plan"]["chosen_tree"] != list(ALTERNATE_TREE):
            raise ArithmeticError("the independent source/cycle plan changed")
        comparison = control["independent_source_basis_comparison"]
        if comparison["comparison_source_minor_plan"]["chosen_tree"] != graph().spanning_tree():
            raise ArithmeticError("the independent comparison source/cycle plan changed")
        if comparison["full_numerator"] != actual["original_numerator"]:
            raise ArithmeticError("the complete original numerator comparison is not bound")
        expected = evaluate(document["reconstruction"], sources)
        for field in ("bare", "all_original_first_current_spinors"):
            if actual[field] != expected[field]:
                raise ArithmeticError("the complete fresh-source result differs: " + name)
        if actual["used_cotree_reconstruction"] or actual["independently_integrated_original_edges"] != list(range(11)):
            raise ArithmeticError("the purported independent full-edge control is incomplete")
        controls[name] = {"file": control_path.name, "sha256": digest(control),
                          "all_23_integrated_values_agree": True,
                          "all_original_spinor_indices_retained": True,
                          "original_source_minor_tree": ALTERNATE_TREE,
                          "whole_source_numerator_basis_check": control["independent_source_basis_comparison"]}
    if write:
        document["independent_controls"] = controls
        atomic_json(path, document)
    elif document["independent_controls"] != json.loads(json.dumps(controls)):
        raise ArithmeticError("the integrated control summaries changed")
    print("Both independent original-graph controls match every source-dependent spinor.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--compute-only", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    if args.compute or args.compute_only or args.recompute:
        compute(force=args.recompute)
    if not args.compute_only:
        check(write=args.write)


if __name__ == "__main__":
    main()
