"""Reproduce complete higher-loop cubic-forest current moments.

Write the exact fixtures, values, source polynomials and resource counters:
    python3 explore/circuit-general/forest-moments/certify_forests.py

Recompute and compare every certificate field except measured runtimes:
    python3 explore/circuit-general/forest-moments/certify_forests.py --check

The two nine-loop examples are independently integrated by the old global
ratio-chart evaluator.  The 36-loop repeated-fork value is checked against
the existing factorial formula, explicitly as a known-formula regression.
The 20-loop alternating tree uses the proved local recurrence; the separate
tests_exact.py supplies exhaustive original-chart controls through five
loops.  Every reported value is a completed integral, with rational sources
and a specified original-edge current word.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))

from htk.graph import Graph
from htk.forest_model import CubicForest, forest_graph
from htk.forest_moments import ForestMoments


def chain_tree(size, *, alternating=False, terminal="A"):
    tree = ("A", "B")
    for step in range(size - 1):
        leaf = ("A" if step % 2 else "B") if alternating else terminal
        tree = (leaf, tree)
    return tree


def generic_sources(graph, roots):
    sources = [(0, 0)] * graph.nv
    sources[graph.nv - 2] = (2, -3)
    for index, root in enumerate(roots, 1):
        sources[root] = (index, index * index + 1)
    sources[-1] = tuple(-sum(row[component] for row in sources) for component in range(2))
    return sources


def fixtures():
    graph, roots = forest_graph((chain_tree(20, alternating=True),))
    sources = [(0, 0)] * graph.nv
    sources[roots[0]], sources[-2], sources[-1] = (2, 3), (-1, 5), (-1, -8)
    yield dict(name="twenty_loop_alternating_tree", graph=graph, sources=sources,
               powers={(0, 0): 1, (graph.ne // 2, 1): 1, (graph.ne - 1, 0): 1},
               old_chart=False)

    graph, roots = forest_graph(((chain_tree(4), chain_tree(4, alternating=True)),))
    yield dict(name="nine_loop_asymmetric_branched_tree", graph=graph,
               sources=generic_sources(graph, roots),
               powers={(0, 0): 2, (graph.ne // 2, 1): 1, (graph.ne - 1, 1): 1},
               old_chart=True)

    fork = (("A", "B"), ("A", "B"))
    graph, roots = forest_graph((chain_tree(2, terminal="B"),
                                 chain_tree(4, alternating=True), fork))
    A, B = graph.nv - 2, graph.nv - 1
    model = CubicForest(graph, (A, B))
    fork_root = next(root for root in roots if model.canonical_nodes[root].encoding == fork)
    child = model.canonical_nodes[fork_root].ports[0]
    first_spoke = next(edge for edge, (left, right) in enumerate(graph.edges)
                       if left == roots[0] and right in (A, B))
    yield dict(name="nine_loop_nonparallel_three_root_forest", graph=graph,
               sources=generic_sources(graph, roots),
               powers={(graph.edges.index((child, A)), 0): 1,
                       (first_spoke, 1): 2, (graph.ne - 1, 0): 1}, old_chart=True)

    count = 12
    A, B = 3 * count, 3 * count + 1
    block_order, spokes = [], []
    sources = [(0, 0)] * (3 * count) + [(1, 0), (-1, -count)]
    for block in range(count):
        root, left, right = 3 * block, 3 * block + 1, 3 * block + 2
        sources[root] = (0, 1)
        block_order.extend(((left, root), (left, A), (left, B),
                            (right, root), (right, A), (right, B)))
        spokes.append((left, A))
    block_order.append((A, B))
    edges = sorted(tuple(sorted(edge)) for edge in block_order)
    graph = Graph(3 * count + 2, edges)
    permutation = [edges.index(tuple(sorted(edge))) for edge in block_order]
    inversions = sum(first > second for index, first in enumerate(permutation)
                     for second in permutation[index + 1:])
    factor = (-1) ** (graph.loops * (graph.loops + 1) // 2 + inversions)
    powers = {(edges.index(edge), 0): 1 for edge in spokes}
    yield dict(name="thirty_six_loop_repeated_fork_known_formula", graph=graph,
               sources=sources, powers=powers, old_chart=False,
               known_formula=dict(fork_count=count, convention_factor=factor,
                                  value=str(Fraction(1, 72 ** count * factorial(4 * count + 1)))))


def original_chart(graph, sources, powers):
    sys.path.insert(0, str(ROOT / "explore/nonSP/verify/full_ward"))
    from forest_jets import evaluate
    word = [location for location, degree in powers.items() for _ in range(degree)]
    return evaluate(tuple(graph.edges), graph.nv - 2, graph.nv - 1, sources, word)


def record(fixture):
    graph, sources, powers = fixture["graph"], fixture["sources"], fixture["powers"]
    marked = (graph.nv - 2, graph.nv - 1)
    start = time.monotonic()
    engine = ForestMoments(graph, sources, marked)
    value = engine.moment(powers)
    elapsed = time.monotonic() - start
    polynomial = engine.source_polynomial(powers)
    stats = dict(engine.last_stats)
    if not value:
        raise ArithmeticError("the selected higher-loop witness unexpectedly vanished")
    A_source = tuple(map(Fraction, sources[marked[0]]))
    reconstructed = sum((coefficient * A_source[0] ** exponent[0] * A_source[1] ** exponent[1]
                         for exponent, coefficient in polynomial.items()), Fraction(0))
    if reconstructed != value:
        raise ArithmeticError("source-polynomial evaluation does not reproduce the moment")
    checks = dict(source_polynomial_exact_equality=True)
    if fixture["old_chart"]:
        reference_start = time.monotonic()
        reference = original_chart(graph, sources, powers)
        if reference != value:
            raise ArithmeticError("independent original global chart disagrees")
        checks["independent_original_chart"] = dict(exact_equality=True, value=str(reference),
                                                    elapsed_seconds=round(time.monotonic() - reference_start, 6))
    if "known_formula" in fixture:
        known = fixture["known_formula"]
        if known["convention_factor"] * value != Fraction(known["value"]):
            raise ArithmeticError("the existing repeated-fork factorial formula disagrees")
        smaller = dict(powers)
        del smaller[next(reversed(smaller))]
        vanishing = engine.moment(smaller)
        if vanishing:
            raise ArithmeticError("omitting a required fork mark did not vanish")
        checks["known_formula_regression"] = dict(
            formula="1/(72^k * (4*k+1)!)", k=known["fork_count"],
            conversion_from_original_htk_value=known["convention_factor"],
            value=known["value"], exact_equality=True,
            one_omitted_mark_value=str(vanishing),
            source="explore/nonSP/verify/full_ward/fork_jets.py")
    source_polynomial = [dict(exponents=list(exponent), coefficient=str(coefficient))
                         for exponent, coefficient in sorted(polynomial.items())]
    result = dict(name=fixture["name"], status="complete exact current moment",
                  convention="code/htk/CONVENTIONS.md: K_Gamma with original oriented edge currents",
                  loops=graph.loops, vertices=graph.nv, original_oriented_edges=graph.edges,
                  marked_vertices=marked, component_roots=engine.model.original_roots,
                  canonical_tree_encodings=engine.model.forest,
                  full_balanced_sources=[[str(value) for value in row] for row in sources],
                  current_powers=[dict(edge=edge, oriented_endpoints=graph.edges[edge],
                                       component=component, power=power)
                                  for (edge, component), power in sorted(powers.items())],
                  total_current_degree=sum(powers.values()), exact_value=str(value),
                  A_source_polynomial=source_polynomial,
                  recognized_cubic_product_forest=True, all_original_sources_retained=True,
                  remaining_integrals=0, contraction_stats=stats, independent_checks=checks,
                  elapsed_seconds=round(elapsed, 6))
    print(result["name"], ":", result["exact_value"], flush=True)
    print("  loops=%d current_degree=%d maximum_terms=%d local_ratios=%d"
          % (result["loops"], result["total_current_degree"],
             stats["largest_polynomial_terms"], stats["largest_local_ratio_count"]), flush=True)
    return result


def certificate():
    files = (ROOT / "code/htk/forest_model.py", ROOT / "code/htk/forest_moments.py",
             HERE / "certify_forests.py")
    return dict(method="local cubic-forest boundary polynomials and exact terminal source translation",
                scope="finite original-edge current moments, H=2, T=0, cubic product forests",
                arithmetic="exact rational coefficients; timings are descriptive only",
                implementation_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                       for path in files},
                independent_low_loop_suite="explore/circuit-general/forest-moments/tests_exact.py",
                examples=[record(fixture) for fixture in fixtures()])


def without_times(value):
    if isinstance(value, dict):
        return {key: without_times(part) for key, part in value.items() if not key.endswith("_seconds")}
    if isinstance(value, list):
        return [without_times(part) for part in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "forest_certificate.json")
    arguments = parser.parse_args()
    observed = certificate()
    # Normalize tuple/list representation exactly as it appears in the JSON.
    serialized = json.dumps(observed, indent=2) + "\n"
    if arguments.check:
        expected = json.loads(arguments.output.read_text())
        if without_times(json.loads(serialized)) != without_times(expected):
            raise ArithmeticError("stored forest certificate differs from exact replay")
        print("All stored forest certificates replay exactly.", flush=True)
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized)
        print("Certificate:", arguments.output, flush=True)


if __name__ == "__main__":
    main()
