"""Compare complete reconstruction with independently integrated old records."""
from fractions import Fraction
import argparse
import json

from certify import (HERE, PREVIOUS, add_records, canonical, check_conservation, digest,
                     forest_graph, forests, reconstruct_currents, source_assignment)
from htk import Graph


def check(records):
    count = 0
    for record in records:
        graph = Graph(len(record["full_balanced_sources"]), record["edges"])
        sources = [tuple(map(Fraction, row)) for row in record["full_balanced_sources"]]
        cotree = graph.cotree(graph.spanning_tree())
        integrated = {edge: record["all_original_first_current_spinors"][edge] for edge in cotree}
        result = canonical(reconstruct_currents(graph, sources, record["bare"], integrated))
        if result != record["all_original_first_current_spinors"]:
            raise ArithmeticError("source-flow reconstruction disagrees with an independently integrated edge")
        check_conservation(record)
        count += 2 * graph.ne
    return count


def check_orientation(records):
    """An odd edge permutation changes the density sign; reversal changes y."""
    count = 0
    for record in records:
        edges = list(map(tuple, record["edges"]))
        permutation = [1, 0] + list(range(2, len(edges)))
        transformed_edges, transformed_currents = [], []
        for position, old in enumerate(permutation):
            reverse = position in (0, len(edges) - 1)
            transformed_edges.append(edges[old][::-1] if reverse else edges[old])
            coefficient = 1 if reverse else -1
            transformed_currents.append([add_records(((coefficient, value),))
                                         for value in record["all_original_first_current_spinors"][old]])
        transformed = dict(record, edges=transformed_edges,
                           bare=add_records(((-1, record["bare"]),)),
                           all_original_first_current_spinors=transformed_currents)
        count += check([canonical(transformed)])
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls-only", action="store_true",
                        help="check the existing references before a new full census has been assembled")
    options = parser.parse_args()
    previous = json.loads((PREVIOUS / "integral_certificate.json").read_text())
    entries = check(previous["records"])
    print("Complete Kirchhoff reconstruction agrees with all %d independently integrated entries through four loops." % entries)
    oriented = check_orientation(previous["records"])
    print("All %d transported entries also preserve reversed currents and the odd edge-permutation sign." % oriented)
    direct_path = HERE / "direct_reference.json"
    if direct_path.exists():
        direct = json.loads(direct_path.read_text())
        if digest(direct["records"]) != direct["records_sha256"]:
            raise ArithmeticError("the complete all-edge reference digest changed")
        if (direct["graph_count"] != len(direct["records"]) or
                len({record["index"] for record in direct["records"]}) != len(direct["records"])):
            raise ArithmeticError("the direct reference has a missing or repeated record")
        entries = check(direct["records"])
        if entries != direct["first_current_entries"]:
            raise ArithmeticError("a direct-reference current entry is missing")
        print("Complete Kirchhoff reconstruction agrees with all %d independently integrated five-loop entries." % entries)
    else:
        raise FileNotFoundError("the complete per-edge five-loop reference is required")
    path = HERE / "certificate.json"
    if path.exists():
        certificate = json.loads(path.read_text())
        if len(certificate["records"]) != len(forests(5)) or digest(certificate["records"]) != certificate["records_sha256"]:
            raise ArithmeticError("the complete finite census digest changed")
        for record in certificate["records"]:
            graph, marked = forest_graph(forests(5)[record["index"]])
            sources, superpotential = source_assignment(graph, marked)
            expected = canonical({"loops": 5, "index": record["index"], "edges": graph.edges,
                                  "marked": marked, "full_balanced_sources": sources,
                                  "superpotential_vertices": superpotential})
            if {key: record[key] for key in expected} != expected:
                raise ArithmeticError("a topology or full source assignment changed")
        if len({record["index"] for record in certificate["records"]}) != 200:
            raise ArithmeticError("the topology census has a missing or repeated entry")
        if (certificate["topologies"] != 200 or certificate["bare_kernels"] != 200 or
                certificate["exact_integrated_kirchhoff_identities"] != 2800 or
                not certificate["all_values_finite"] or certificate["residual_periods"]):
            raise ArithmeticError("the complete census status or declared counts changed")
        entries = check(certificate["records"])
        if entries != 4400 or certificate["complete_first_current_entries"] != entries:
            raise ArithmeticError("a complete first-current entry is missing")
        if direct_path.exists():
            values = {record["index"]: record for record in certificate["records"]}
            for record in direct["records"]:
                if values[record["index"]] != record:
                    raise ArithmeticError("the finite census differs from a complete all-edge reference")
        print("Every topology, source, finite current entry and complete integrated source identity is retained.")
    elif not options.controls_only:
        raise FileNotFoundError("the complete 200-forest census certificate is required")


if __name__ == "__main__":
    main()
