"""Complete one full-source fixture for every admissible six-loop K5 core.

Three vertices can be distributed among the ten K5 edges in 220 ways.
The exhaustive admissibility/isomorphism check leaves exactly three shapes.
This script computes their bare T=0 kernels at the deterministic source
vectors from run_examples.sources(graph, 2, 10602); it does not claim a
symbolic result for arbitrary source vectors.

Run:
  python3 explore/circuit-general/series-fibres/run_k5_family.py --threads 4

Use --recompute to reconstruct every polynomial again (validated modular
checkpoints may be used), replay both rational fibre integrals, and repeat
the complete period contraction. Without it, existing fixture certificates
are checked for exact graph/source and current implementation agreement.
"""
import argparse
from collections import Counter
from contextlib import redirect_stdout
from fractions import Fraction
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import sys
import tempfile
import time

import networkx as nx

from run_examples import HERE, catalog, sources, subdivide
from htk.series import SeriesReduction
from core_numerator import certified_core_numerator, integrate_core, _fingerprint
from reciprocal import reciprocal_plan
from run_core import replay_fibres, certify_nonzero


CORE_EDGES = tuple(combinations(range(5), 2))
SOURCE_SEED = 10602
FIXTURES = (
    ("two_edge_path_plus_disjoint_edge", (0, 7, 9), "six_loop_k5.json"),
    ("triangle", (0, 1, 4), "six_loop_k5_triangle.json"),
    ("three_edge_path", (0, 4, 7), "six_loop_k5_three_edge_path.json"),
)
_SHAPE_BY_DEGREES = {
    (0, 0, 2, 2, 2): "triangle",
    (0, 1, 1, 2, 2): "three_edge_path",
    (1, 1, 1, 1, 2): "two_edge_path_plus_disjoint_edge",
    (0, 1, 1, 1, 3): "three_edge_star",
}


def graph_for(allocation):
    counts = Counter(allocation)
    return subdivide(5, CORE_EDGES, [1 + counts[edge] for edge in range(10)])


def _shape(allocation):
    if len(set(allocation)) == 1:
        return "three_insertions_on_one_edge"
    if len(set(allocation)) == 2:
        first, second = (CORE_EDGES[edge] for edge in sorted(set(allocation)))
        return ("two_insertions_plus_adjacent_edge" if set(first) & set(second)
                else "two_insertions_plus_disjoint_edge")
    degrees = [0] * 5
    for edge in allocation:
        for vertex in CORE_EDGES[edge]:
            degrees[vertex] += 1
    return _SHAPE_BY_DEGREES[tuple(sorted(degrees))]


def _sparsity_failure(graph):
    """Independent exhaustive H=2 induced-subgraph test, with a witness."""
    for size in range(2, graph.nv + 1):
        for vertices in combinations(range(graph.nv), size):
            selected = set(vertices)
            edges = sum(u in selected and v in selected for u, v in graph.edges)
            if edges and edges > 2 * size - 3:
                return dict(vertices=list(vertices), edges=edges, allowed_edges=2 * size - 3)
    return None


def enumerate_family():
    """Check all weak compositions of three inserted vertices and exact isomorphisms."""
    classes = []
    total, admissible, repeated = 0, 0, 0
    for allocation in combinations_with_replacement(range(10), 3):
        graph = graph_for(allocation)
        network = nx.Graph(graph.to_networkx())
        failure = _sparsity_failure(graph)
        okay = failure is None
        if okay != catalog.tight_ok(network, 2):
            raise ArithmeticError("independent admissibility tests disagree")
        total += 1
        admissible += okay
        repeated += len(set(allocation)) < 3
        if len(set(allocation)) < 3 and okay:
            raise ArithmeticError("repeated subdivisions violate the five-core-vertex bound")
        for record in classes:
            if nx.is_isomorphic(network, record["network"]):
                if okay != record["admissible"] or _shape(allocation) != record["shape"]:
                    raise ArithmeticError("isomorphic subdivisions received different classifications")
                record["labeled_allocations"] += 1
                break
        else:
            classes.append(dict(network=network, shape=_shape(allocation), admissible=okay,
                                labeled_allocations=1, representative_allocation=list(allocation),
                                representative_oriented_edges=graph.edges,
                                first_failed_induced_subgraph=failure))
    for first, second in combinations(classes, 2):
        if nx.is_isomorphic(first["network"], second["network"]):
            raise ArithmeticError("enumeration contains duplicate isomorphism classes")
    good = {record["shape"]: record["labeled_allocations"]
            for record in classes if record["admissible"]}
    expected = {"triangle": 10, "three_edge_path": 60,
                "two_edge_path_plus_disjoint_edge": 30}
    if total != 220 or repeated != 100 or admissible != 100 or good != expected or len(classes) != 7:
        raise ArithmeticError("unexpected exhaustive K5 subdivision census")
    for record in classes:
        del record["network"]
    return dict(core_vertices=5, core_edges=CORE_EDGES, inserted_vertices=3,
                candidate_allocations=total, repeated_edge_allocations=repeated,
                admissible_labeled_allocations=admissible,
                inadmissible_labeled_allocations=total - admissible,
                all_isomorphism_classes=len(classes), admissible_isomorphism_classes=3,
                admissibility="Every induced vertex subset with an edge obeys E <= 2V - 3; total E=13, V=8.",
                admissibility_derivation=[
                    "The five core vertices permit at most seven induced edges, forcing three distinct subdivided edges.",
                    "Every four-core-vertex subset must contain a subdivided edge. Among three-edge patterns this excludes exactly the star.",
                    "Smaller core subsets satisfy the bound. Adding a bivalent vertex contributes at most two edges and one vertex, preserving or weakening the bound; subsets with fewer than two core vertices are forests.",
                ],
                isomorphism_test="networkx.is_isomorphic on the complete subdivided graph",
                completeness="Suppressing the three degree-two vertices recovers K5. Every K5-core graph with V=8 and E=13 is one of these 220 allocations.",
                classes=sorted(classes, key=lambda item: item["shape"]))


def _validate_fixture(record, graph, expected_sources, *, verify_implementation=True):
    if (record.get("vertices") != graph.nv
            or [tuple(edge) for edge in record.get("original_oriented_edges", ())] != graph.edges
            or [[Fraction(value) for value in row] for row in record.get("full_balanced_sources", ())]
            != [[Fraction(value) for value in row] for row in expected_sources]):
        raise ValueError("stored fixture graph, orientations, or source vectors do not match exactly")
    if (record.get("loops") != 6 or record.get("holomorphic_dimension") != 2
            or record.get("q") != 8 or record.get("remaining_projective_integrals") != 0
            or not record.get("all_original_sources_retained")
            or not record.get("graph_admissibility_exactly_checked")):
        raise ValueError("stored fixture is not a complete admissible full-source K5 kernel")
    proof = record["polynomial_certificate"]
    if verify_implementation:
        reduction = SeriesReduction(graph, expected_sources)
        current = _fingerprint(reduction, reciprocal_plan(reduction))
        if proof.get("fingerprint") != current:
            raise ValueError("stored fixture implementation fingerprint changed; use --recompute")
    if (not proof.get("strict_uniqueness_inequality")
            or int(proof["modulus"]) <= 2 * int(proof["coefficient_bound"])
            or not all(prime.get("homogeneous_check") for prime in proof["primes"])):
        raise ValueError("stored fixture lacks the integer-CRT uniqueness certificate")
    replays = record.get("independent_rational_fibre_replays", ())
    if len(replays) != 2 or not all(replay.get("exact_equality") for replay in replays):
        raise ValueError("stored fixture lacks both independent rational fibre replays")
    table = HERE.parent / "k5" / "tables" / "K5_q8.json"
    if record["period_certificate"]["sha256"] != hashlib.sha256(table.read_bytes()).hexdigest():
        raise ValueError("stored fixture refers to a different K5 period table")
    values = {label: Fraction(value) for label, value in record["exact_value"].items()}
    if record.get("nonzero_certificate") != certify_nonzero(values):
        raise ValueError("stored nonzero certificate differs from the recalculated rational bounds")


def _compute_fixture(shape, graph, source_rows, arguments):
    reduction = SeriesReduction(graph, source_rows)
    if not catalog.tight_ok(graph.to_networkx(), 2):
        raise ArithmeticError("family fixture is not admissible")
    start = time.monotonic()
    polynomial = certified_core_numerator(reduction, threads=arguments.threads,
                                          checkpoint_dir=arguments.checkpoint_dir)
    replays = replay_fibres(reduction, polynomial)
    print(shape + ": two independent rational fibre replays passed", flush=True)
    values, period_certificate = integrate_core(polynomial)
    period_certificate["file"] = "explore/circuit-general/k5/tables/K5_q8.json"
    name = "six_loop_K5_H2" if shape == "two_edge_path_plus_disjoint_edge" else "six_loop_K5_H2_" + shape
    return dict(name=name,
                status="complete bare T=0 projective kernel",
                convention="code/htk/CONVENTIONS.md: K_Gamma(lambda; z=0)",
                loops=6, holomorphic_dimension=2, vertices=graph.nv,
                original_oriented_edges=graph.edges, full_balanced_sources=source_rows,
                source_seed=SOURCE_SEED, all_original_sources_retained=True,
                graph_admissibility_exactly_checked=True, core_edges=polynomial.edges,
                q=polynomial.q, chain_edge_ids=[chain.edges for chain in reduction.chains],
                original_determinant_size=graph.ne - 1,
                contracted_determinant_size=reduction.determinant_size,
                finite_field_cubature_nodes=27,
                exact_value={label: str(value) for label, value in sorted(values.items())},
                nonzero_certificate=certify_nonzero(values),
                polynomial_certificate=polynomial.certificate,
                independent_rational_fibre_replays=replays,
                period_certificate=period_certificate,
                remaining_projective_integrals=0,
                elapsed_seconds=round(time.monotonic() - start, 3))


def _write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, value):
        for stream in self.streams:
            stream.write(value)
        return len(value)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def run(arguments):
    census = enumerate_family()
    print("K5 allocations: 220 checked; 100 admissible; 3 admissible isomorphism classes", flush=True)
    for entry in census["classes"]:
        print("%s: %d labeled allocations; admissible=%s" %
              (entry["shape"], entry["labeled_allocations"], entry["admissible"]), flush=True)
    records = []
    for shape, selected, filename in FIXTURES:
        graph = graph_for(selected)
        source_rows = sources(graph, 2, SOURCE_SEED)
        if any(not any(row) for row in source_rows) or any(sum(row[i] for row in source_rows) for i in range(2)):
            raise ArithmeticError("source fixture must be balanced with every vertex vector nonzero")
        path = arguments.first_certificate if filename == "six_loop_k5.json" else arguments.output.parent / filename
        if path.is_file() and not arguments.recompute:
            record = json.loads(path.read_text())
            enriched = "nonzero_certificate" not in record
            if enriched:
                record["nonzero_certificate"] = certify_nonzero(record["exact_value"])
            _validate_fixture(record, graph, source_rows)
            print(shape + ": reused exact graph/source-matched certificate " + str(path), flush=True)
            mode = "validated_existing_certificate"
            if enriched or record["period_certificate"]["file"] != "explore/circuit-general/k5/tables/K5_q8.json":
                record["period_certificate"]["file"] = "explore/circuit-general/k5/tables/K5_q8.json"
                _write_json(path, record)
        elif filename == "six_loop_k5.json" and not arguments.recompute:
            raise FileNotFoundError("run run_core.py first to create " + str(path))
        else:
            print(shape + ": reconstructing the complete core polynomial", flush=True)
            record = _compute_fixture(shape, graph, source_rows, arguments)
            _validate_fixture(record, graph, source_rows, verify_implementation=False)
            _write_json(path, record)
            mode = "reconstructed_with_two_rational_replays_and_complete_period_contraction"
        record = dict(record, family_shape=shape, subdivided_core_edge_ids=list(selected),
                      subdivided_core_edges=[CORE_EDGES[edge] for edge in selected],
                      source_seed=SOURCE_SEED, certificate_use=mode,
                      fixture_file=path.name, fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        print(shape + ": " + json.dumps(record["exact_value"], sort_keys=True), flush=True)
        bound = record["nonzero_certificate"]
        print(shape + ": exact nonzero interval [" + bound["lower"] + ", " + bound["upper"] + "]", flush=True)
        records.append(record)
    result = dict(scope="All three admissible K5-core topologies at six loops and H=2; one complete bare T=0 kernel at the displayed fixed full source assignment for each topology.",
                  symbolic_arbitrary_sources=False, all_admissible_topologies_covered=True,
                  every_original_vertex_source_nonzero=True, all_three_values_proved_nonzero=True,
                  source_seed=SOURCE_SEED,
                  enumeration=census, examples=records, remaining_projective_integrals=0)
    _write_json(arguments.output, result)
    print("All three full-source fixtures are complete; remaining projective integrals: 0", flush=True)
    print("Family certificate: " + str(arguments.output), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--checkpoint-dir", type=Path,
                        default=Path(tempfile.gettempdir()) / "htk-certified-core")
    parser.add_argument("--first-certificate", type=Path, default=HERE / "six_loop_k5.json")
    parser.add_argument("--output", type=Path, default=HERE / "k5_family.json")
    parser.add_argument("--log", type=Path, default=HERE / "k5_family.log")
    parser.add_argument("--enumerate-only", action="store_true")
    parser.add_argument("--recompute", action="store_true",
                        help="reconstruct, independently replay, and integrate all three fixtures; checkpoints may be used")
    arguments = parser.parse_args()
    if arguments.threads < 1:
        parser.error("--threads must be positive")
    if arguments.enumerate_only:
        print(json.dumps(enumerate_family(), indent=2))
        return
    arguments.log.parent.mkdir(parents=True, exist_ok=True)
    with arguments.log.open("w") as log:
        with redirect_stdout(_Tee(sys.stdout, log)):
            run(arguments)


if __name__ == "__main__":
    main()
