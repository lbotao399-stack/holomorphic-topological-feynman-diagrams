#!/usr/bin/env python3
"""Exhaust the tight cubic four-loop graphs independently.

No graph atlas, NetworkX, or original enumeration routine is imported.
The paper proves that suppression of the four bivalent vertices gives
a loopless biconnected cubic multigraph on six vertices, with edge
multiplicity at most two. We enumerate that finite set directly from
the vertex-degree constraints, then all four-subdivision distributions.

A graph is identified by its core endpoint pairs and path lengths.
Every automorphism preserves the six degree-three vertices; after its
action on the core is fixed, its remaining freedom permutes equal-length
parallel paths. This proves the automorphism formula used below.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement, permutations
import json
from math import factorial, prod
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PAPER = Path(__file__).resolve().parents[1]
PERMUTATIONS = tuple(permutations(range(6)))


def connected(edges, removed=None):
    vertices = set(range(6)) - {removed}
    reached = {min(vertices)}
    while True:
        enlarged = reached | {
            target
            for first, second in edges
            for source, target in ((first, second), (second, first))
            if source in reached and target in vertices
        }
        if enlarged == reached:
            return reached == vertices
        reached = enlarged


def core_signature(edges, order):
    return tuple(sorted((min(order[a], order[b]), max(order[a], order[b]))
                        for a, b in edges))


def cores():
    """Every labelled symmetric integer adjacency array of row sum three."""
    multiplicities = Counter()

    def visit(vertex, remaining, edges):
        if vertex == 6:
            assert all(value == 0 for value in remaining)
            if all(connected(edges, removed)
                   for removed in (None, 0, 1, 2, 3, 4, 5)):
                canonical = min(core_signature(edges, order)
                                for order in PERMUTATIONS)
                multiplicities[canonical] += 1
            return
        for targets in combinations_with_replacement(
                range(vertex + 1, 6), remaining[vertex]):
            counts = Counter(targets)
            if any(count > 2 or count > remaining[target]
                   for target, count in counts.items()):
                continue
            next_remaining = remaining.copy()
            next_remaining[vertex] = 0
            for target, count in counts.items():
                next_remaining[target] -= count
            visit(vertex + 1, next_remaining,
                  edges + [(vertex, target) for target in targets])

    visit(0, [3] * 6, [])
    return multiplicities


def weighted_signature(edges, lengths, order):
    return tuple(sorted((min(order[a], order[b]), max(order[a], order[b]), length)
                        for (a, b), length in zip(edges, lengths)))


def stored_signature(chains):
    edges = [(chain[0], chain[-1]) for chain in chains]
    lengths = [len(chain) - 1 for chain in chains]
    return min(weighted_signature(edges, lengths, order)
               for order in PERMUTATIONS)


def subdivide(edges, lengths):
    next_vertex = 6
    expanded = []
    for (first, last), length in zip(edges, lengths):
        path = [first] + list(range(next_vertex, next_vertex + length - 1)) + [last]
        next_vertex += length - 1
        expanded.extend(zip(path, path[1:]))
    assert next_vertex == 10 and len(expanded) == 13
    return expanded


def main():
    core_counts = cores()
    assert sum(core_counts.values()) == 550 and len(core_counts) == 5
    generated = {}
    core_rows = []
    checked_subsets = 0
    for core, labelled_count in sorted(core_counts.items()):
        automorphisms = [order for order in PERMUTATIONS
                         if core_signature(core, order) == core]
        # Orbit-stabilizer is a second check on the labelled enumeration.
        assert labelled_count * len(automorphisms) == factorial(6)
        seen = set()
        retained = 0
        for split in combinations_with_replacement(range(9), 4):
            counts = Counter(split)
            lengths = [1 + counts[edge] for edge in range(9)]
            key = min(weighted_signature(core, lengths, order)
                      for order in automorphisms)
            if key in seen:
                continue
            seen.add(key)
            expanded = subdivide(core, lengths)
            edge_masks = [(1 << first) | (1 << second)
                          for first, second in expanded]
            tight = True
            for subset in range(1 << 10):
                if subset.bit_count() < 2:
                    continue
                checked_subsets += 1
                internal_edges = sum(subset & edge == edge for edge in edge_masks)
                if 2 * internal_edges > 3 * subset.bit_count() - 4:
                    tight = False
                    break
            if not tight:
                continue
            canonical = min(weighted_signature(core, lengths, order)
                            for order in PERMUTATIONS)
            assert canonical not in generated
            identity = weighted_signature(core, lengths, tuple(range(6)))
            stabilizer = sum(weighted_signature(core, lengths, order) == identity
                             for order in automorphisms)
            path_permutations = prod(factorial(count)
                                     for count in Counter(identity).values())
            generated[canonical] = stabilizer * path_permutations
            retained += 1
        core_rows.append({
            'edges': core,
            'labelled_representatives': labelled_count,
            'vertex_automorphisms': len(automorphisms),
            'subdivision_orbits_examined': len(seen),
            'tight_subdivision_orbits': retained,
        })

    source = ROOT / 'deepen/o3/higher_orders/four-loop-graphs.json'
    supplied = json.loads(source.read_text())
    expected = {stored_signature(record['chains']): record['aut']
                for record in supplied}
    assert len(expected) == len(supplied) == 50
    assert generated == expected
    retained_counts = [row['tight_subdivision_orbits'] for row in core_rows]
    assert sorted(retained_counts) == [0, 3, 8, 11, 28]
    result = {
        'arithmetic': 'exact integer arithmetic',
        'enumeration': 'degree constraints; all vertex permutations; all subdivision distributions',
        'original_enumeration_or_graph_atlas_imported': False,
        'eligible_labelled_cores': sum(core_counts.values()),
        'eligible_unlabelled_cores': len(core_counts),
        'core_classes': core_rows,
        'retained_graphs': len(generated),
        'surviving_core_classes': sum(count > 0 for count in retained_counts),
        'induced_subset_inequalities_examined': checked_subsets,
        'all_isomorphism_classes_equal_supplied_table': True,
        'all_automorphism_orders_equal_supplied_table': True,
        'automorphism_order_histogram': dict(sorted(Counter(generated.values()).items())),
        'source_sha256': {str(source.relative_to(ROOT)): sha256(source.read_bytes()).hexdigest()},
    }
    destination = PAPER / 'verification/four_loop_census.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
