"""Quick exact census replay: all outputs, source relations and declared inputs.

This checker does not repeat expensive scalar integration.  It regenerates
the complete graph/source relations, verifies every reconstructed original
current, rechecks independent certificates and source bindings, and writes
or compares an aggregate containing all complete answers and explicit gaps.
Run independent_controls.py for fresh complete original integrations.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from bulk3_periods import HERE, ROOT, add, encoded
from certify import ROWS, atomic_json, hashes, zero_witness
from current_relations import CurrentRelations
from htk import brute
from htk.graph import Graph
from htk.treeflow import tree_flow


EXPECTED_BARE = set(range(69))
EXPECTED_FIRST = set(range(14)) | set(range(46, 69))


def canonical(value):
    return json.loads(json.dumps(value, sort_keys=True))


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def decode(value):
    require(isinstance(value, dict), "A complete period value is not a dictionary.")
    result = {key: Fraction(coefficient) for key, coefficient in value.items()}
    require(encoded(result) == value, "A complete coefficient is noncanonical or contains a discarded zero.")
    require(set(result) <= {"1", "z3", "z5", "z3*z3"}, "An undeclared period channel appeared.")
    return result


def check_bindings(bindings):
    for name, expected in bindings.items():
        path = ROOT / name
        require(path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                "A bound exact source or artifact changed: " + name)


def check_first(data):
    row = ROWS[data["index"]]
    bare = decode(data["bare_coefficient"])
    values = [[decode(value) for value in spinor] for spinor in data["all_original_first_current_coefficients"]]
    require(len(values) == 16 and all(len(spinor) == 3 for spinor in values),
            "A complete original edge or source direction is missing.")
    relations = CurrentRelations(row)
    relations.check(bare, values)
    targets = {(entry["edge"], entry["source_direction"]): decode(entry["value"])
               for entry in data["complete_finite_targets"] if entry["edge"] is not None}
    if "complete_current_relations" in data:
        require(canonical(relations.record()) == data["complete_current_relations"],
                "A complete current reconstruction certificate changed.")
        independent = {position: targets[divmod(position, 3)] for position in relations.free}
        reconstructed = relations.reconstruct(bare, independent)
    else:
        graph = Graph(row["nv"], row["edges"])
        sources = [[Fraction(0)] * 3 for _ in range(graph.nv)]
        for a, vertex in enumerate(row["S0"][:3]):
            sources[vertex][a] = Fraction(1)
        sources[row["S0"][-1]] = [Fraction(-1)] * 3
        tree = data["source_cycle_tree"]
        cycles, flow = graph.fundamental_cycles(tree), tree_flow(graph, sources, tree)
        cotree = graph.cotree(tree)
        reconstructed = []
        for edge in range(graph.ne):
            spinor = []
            for a in range(3):
                value = add({}, bare, flow[edge][a])
                for r, chord in enumerate(cotree):
                    value = add(value, targets[chord, a], cycles[r][edge])
                spinor.append(value)
            reconstructed.append(spinor)
    require(reconstructed == values, "A complete stored original current differs from its exact reconstruction.")
    return len(relations.rows)


def check_independent():
    path = HERE / "independent_certificate.json"
    data = json.loads(path.read_text())
    check_bindings(data["dependencies_sha256"])
    require(data["summary"] == {
        "full_current_tensors_reintegrated": 5,
        "all_original_current_coefficients_reintegrated": 240,
        "complete_scalar_dictionaries_compared": 63,
        "all_scalar_coefficients_compared": 101517,
        "complete_nonunit_or_degenerate_source_polynomials": 6,
        "full_original_current_source_transformation_checks": 288}, "Independent full-control scope changed.")
    for record in data["complete_nonsp_controls"]:
        saved = json.loads((HERE / "cases" / ("bulk3_%03d.json" % record["index"])).read_text())
        require(record["bare_coefficient"] == saved["bare_coefficient"]
                and record["all_original_first_current_coefficients"] == saved["all_original_first_current_coefficients"],
                "An independently reintegrated full tensor differs from the census.")
    return data


def build():
    baseline = hashes()
    symmetry_path = HERE / "symmetry_certificate.json"
    symmetry = json.loads(symmetry_path.read_text())
    check_bindings(symmetry["dependencies_sha256"])
    witnesses = {record["index"]: record for record in symmetry["complete_bare_zeros"]}
    require(len(ROWS) == 410 and len(witnesses) == 27, "The complete source-supported graph universe changed.")
    for row in ROWS:
        edges = list(map(tuple, row["edges"]))
        degree = Counter(vertex for edge in edges for vertex in edge)
        require(row["nv"] == 12 and len(edges) == len(set(edges)) == 16
                and all(left < right for left, right in edges)
                and sorted(v for v, count in degree.items() if count == 2) == row["S0"]
                and len(row["S0"]) == 4 and all(count in (2, 3) for count in degree.values()),
                "The original simple cubic graph or complete bivalent source support changed.")
        witness = zero_witness(row)
        require((witness is not None) == (row["index"] in witnesses),
                "The exhaustive signed-automorphism classification changed.")
        if witness is not None:
            stored = witnesses[row["index"]]
            for name in ("vertex_permutation", "source_permutation", "vertex_sign", "reversed_original_edges",
                         "source_determinant_sign"):
                require(witness[name] == stored[name], "A complete zero automorphism changed.")
    independent = check_independent()
    bare_outputs, first_outputs = [], []
    family_counts, first_counts, channels = Counter(), Counter(), Counter()
    target_count = coefficient_count = relation_count = density_points = 0
    case_paths = sorted((HERE / "cases").glob("bulk3_*.json"))
    require({int(path.stem.split("_")[1]) for path in case_paths} == EXPECTED_BARE,
            "The declared finite census must contain exactly original indices0 through68.")
    for path in case_paths:
        data = json.loads(path.read_text())
        index = data["index"]
        require(index == int(path.stem.split("_")[1]), "A complete case has the wrong original index.")
        row = ROWS[index]
        require(data["ordered_oriented_edges"] == row["edges"] and data["source_vertices"] == row["S0"]
                and data["vertices"] == row["nv"] and data["core_classification"] == row["core"],
                "An original graph, edge orientation, or source support changed.")
        require(data["dependencies_sha256"] == baseline, "A finite case has stale baseline sources.")
        check_bindings(data.get("planning_dependencies_sha256", {}))
        require(data["bare_status"] == "complete", "An announced bare polynomial is incomplete.")
        bare = decode(data["bare_coefficient"])
        if index in witnesses:
            require(not bare, "An exact bare coefficient violates its signed automorphism.")
        targets = data["complete_finite_targets"]
        if data.get("method") == "complete_signed_automorphism":
            require(index in witnesses and not targets, "A zero certificate discarded unresolved finite targets.")
            require(data["independent_signed_zero_certificate_sha256"] == hashlib.sha256(symmetry_path.read_bytes()).hexdigest(),
                    "The complete signed-zero proof changed.")
        else:
            entries = [entry for entry in targets if entry["edge"] is None]
            require(len(entries) == 1 and entries[0]["value"] == data["bare_coefficient"],
                    "The complete scalar bare target is missing or inconsistent.")
            graph = Graph(row["nv"], row["edges"])
            sources = [[Fraction(0)] * 3 for _ in range(graph.nv)]
            for a, vertex in enumerate(row["S0"][:3]):
                sources[vertex][a] = Fraction(1)
            sources[row["S0"][-1]] = [Fraction(-1)] * 3
            checks = data["original_grounded_laplacian_density_checks"]
            require(len(checks) == 2, "An original-density control is missing.")
            for control in checks:
                point = list(map(Fraction, control["resistances"]))
                U = brute.symanzik_U(graph, point)
                require(str(U) == control["U"] and brute.density(graph, point, sources).rational() * U ** 8
                        == Fraction(control["complete_numerator"]), "An exact original-density control failed.")
                density_points += 1
        for entry in targets:
            decode(entry["value"])
            require(2 <= entry["denominator_power"] <= 9, "A finite scalar target has an invalid denominator.")
            require(len(entry["complete_scalar_coefficients_sha256"]) == 64,
                    "The complete scalar target digest is missing.")
            if entry["weighted_core_periods"]:
                require(entry["period_stats"]["read_only_tables"] is True,
                        "An undeclared table-building path was used.")
                require(entry["period_stats"]["q"] == entry["denominator_power"]
                        and entry["period_stats"]["weighted_periods"] == entry["weighted_core_periods"],
                        "The complete period lookup count or denominator changed.")
            target_count += 1
            coefficient_count += entry["weighted_core_periods"]
        output = {"index": index, "core": row["core"], "source_vertices": row["S0"],
                  "ordered_oriented_edges": row["edges"], "bare_coefficient": data["bare_coefficient"]}
        bare_outputs.append(output)
        family_counts[row["core"]] += 1
        channels.update(bare.keys())
        if index in EXPECTED_FIRST:
            require(data["first_current_status"] == "complete", "An announced complete current tensor is incomplete.")
            relation_count += check_first(data)
            values = data["all_original_first_current_coefficients"]
            first_outputs.append({"index": index, "source_vertices": row["S0"],
                                  "ordered_oriented_edges": row["edges"],
                                  "all_original_first_current_coefficients": values})
            first_counts[row["core"]] += 1
            for spinor in values:
                for value in spinor:
                    channels.update(decode(value).keys())
        else:
            require(data["first_current_status"] == "not_requested"
                    and "all_original_first_current_coefficients" not in data,
                    "An unavailable prism q9 tensor was misclassified.")
    for index in sorted(set(witnesses) - EXPECTED_BARE):
        row = ROWS[index]
        bare_outputs.append({"index": index, "core": row["core"], "source_vertices": row["S0"],
                             "ordered_oriented_edges": row["edges"], "bare_coefficient": {},
                             "method": "complete_signed_automorphism"})
        family_counts[row["core"]] += 1
    complete_bare = {item["index"] for item in bare_outputs}
    require(len(complete_bare) == 92 and len(first_outputs) == 37, "The declared complete result counts changed.")
    local_sources = list(HERE.glob("*.py")) + [HERE / "README.md"]
    inherited = list((ROOT / "explore/circuit-general/cores").glob("*.py"))
    inherited += list((ROOT / "explore/circuit-general/fiveloop").glob("*.py"))
    inherited += [ROOT / "code/htk/spred.py", ROOT / "code/htk/simplex.py", ROOT / "code/htk/__init__.py"]
    bindings = dict(baseline)
    bindings.update({str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in local_sources + inherited})
    artifacts = case_paths + [symmetry_path, HERE / "independent_certificate.json"]
    return {"schema": "five-loop-bulk3-complete-supported-source-results-v1",
            "source_scope": "arbitrary lambda_s0,lambda_s1,lambda_s2; lambda_s3=-sum; every other original source zero",
            "bare_reconstruction": "det(lambda_s0,lambda_s1,lambda_s2)^5 times bare_coefficient",
            "current_reconstruction": "the same determinant power times sum_j lambda_sj^a times all_original_first_current_coefficients[edge][j]",
            "summary": {"total_bulk3_graphs": 410, "complete_bare_source_polynomials": 92,
                        "new_nonsp_bare_source_polynomials": 89,
                        "complete_first_current_tensors": 37, "new_nonsp_first_current_tensors": 34,
                        "complete_original_first_current_vectors": 592,
                        "complete_first_current_source_coefficients": 1776,
                        "complete_bare_signed_zeros": 27,
                        "original_density_points_replayed": density_points,
                        "complete_integrated_relation_rows_replayed": relation_count,
                        "complete_scalar_targets": target_count,
                        "complete_scalar_coefficients": coefficient_count,
                        "period_channels": sorted(channels)},
            "complete_bare_by_core": dict(sorted(family_counts.items())),
            "complete_first_current_by_core": dict(sorted(first_counts.items())),
            "all_complete_bare_kernels": sorted(bare_outputs, key=lambda item: item["index"]),
            "all_complete_first_current_kernels": first_outputs,
            "unresolved_bare_bulk3_indices": sorted(set(range(410)) - complete_bare),
            "unresolved_first_current_bulk3_indices": sorted(set(range(410)) - EXPECTED_FIRST),
            "unchanged_d2_gaps": {"Z5": [49, 50, 51, 53, 54, 58], "K33+e": [57, 61, 65],
                                 "V7_0": [66], "V7_2": [68], "V7_3": [67]},
            "independent_control_summary": independent["summary"],
            "source_sha256": bindings,
            "artifacts_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in artifacts}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build()
    path = HERE / "certificate.json"
    if args.write:
        atomic_json(path, result)
    else:
        require(result == json.loads(path.read_text()), "The complete census aggregate or a source binding changed.")
    print(json.dumps(result["summary"], indent=2))
    print("Quick replay checked complete stored values and exact relations; it did not repeat all period integrations.")
    print("certificate_sha256:", hashlib.sha256(path.read_bytes()).hexdigest())
