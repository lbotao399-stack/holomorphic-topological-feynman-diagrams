"""Independent exhaustive simple-graph control of the product-forest census."""
from itertools import combinations
import json

from catalog import canonical_tree, count_table, forest_graph, forests


def independent_key(vertices, edges):
    A, B = vertices - 2, vertices - 1
    adjacent = [set() for _ in range(vertices)]
    for first, second in edges:
        adjacent[first].add(second)
        adjacent[second].add(first)
    interactions = set(range(A))
    if any(len(adjacent[vertex]) < 2 for vertex in interactions):
        return None
    # Enumerate every vertex subset, independently of the production
    # maximum-connected-subtree criterion.
    for order in range(2, vertices):
        for subset in combinations(range(vertices), order):
            subset = set(subset)
            if sum(first in subset and second in subset for first, second in edges) > 2 * order - 3:
                return None
    remaining, components = set(interactions), []
    while remaining:
        seen, pending = set(), [min(remaining)]
        while pending:
            vertex = pending.pop()
            if vertex in seen:
                continue
            seen.add(vertex)
            pending.extend((adjacent[vertex] & interactions) - seen)
        remaining -= seen
        if sum(len(adjacent[vertex] & seen) for vertex in seen) != 2 * (len(seen) - 1):
            return None
        if sum(len(adjacent[vertex]) > 3 for vertex in seen) > 1:
            return None
        if sum(len(adjacent[vertex] & {A, B}) for vertex in seen) != len(seen) + 1:
            return None
        if not any(A in adjacent[vertex] for vertex in seen) or not any(B in adjacent[vertex] for vertex in seen):
            return None
        order = sorted(seen)
        labels = {vertex: index for index, vertex in enumerate(order)}
        tree = tuple(tuple(labels[neighbor] for neighbor in sorted(adjacent[vertex] & seen))
                     for vertex in order)
        spokes = tuple(int(A in adjacent[vertex]) + 2 * int(B in adjacent[vertex]) for vertex in order)
        components.append(canonical_tree(tree, spokes))
    return tuple(sorted(components))


def check():
    records = []
    for loops in range(1, 5):
        vertices = loops + 2
        marked_edge = (loops, loops + 1)
        possible = [edge for edge in combinations(range(vertices), 2) if edge != marked_edge]
        enumerated, expected, labeled, examined = set(), set(forests(loops)), 0, 0
        for physical in combinations(possible, 2 * loops):
            examined += 1
            key = independent_key(vertices, physical + (marked_edge,))
            if key is not None:
                enumerated.add(key)
                labeled += 1
        if enumerated != expected:
            raise ArithmeticError("the tree census disagrees with exhaustive original simple graphs")
        for forest in expected:
            graph, _ = forest_graph(forest)
            if independent_key(graph.nv, graph.edges) != forest:
                raise ArithmeticError("canonical tree reconstruction changed a marked topology")
        records.append({"loops": loops, "all_simple_graphs_examined": examined,
                        "labeled_forests": labeled, "unlabeled_forests": len(enumerated)})
    return {"exhaustive_simple_graph_controls": records, "counts_through_six": count_table(6)}


if __name__ == "__main__":
    print(json.dumps(check(), indent=2))
