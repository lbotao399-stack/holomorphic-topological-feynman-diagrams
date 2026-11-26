"""Cubic BF graphs marked at the operator. Generate isomorphism classes."""

from itertools import combinations
from pathlib import Path
import sys

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "code"))
from htk.catalog import sparse_ok


def canonical(graph, spokes=()):
    """Canonical decorated-tree code, or a dihedral necklace of rooted trees."""
    spokes = frozenset(spokes)

    def branch(vertex, excluded):
        children = sorted(
            branch(child, {vertex}) for child in graph[vertex] if child not in excluded
        )
        return ("1" if vertex in spokes else "0") + "(" + "".join(children) + ")"

    cycles = nx.cycle_basis(graph)
    if not cycles:
        return "T" + min(branch(center, set()) for center in nx.center(graph))
    (cycle,) = cycles
    cores = [branch(v, set(cycle)) for v in cycle]
    necklaces = [
        sequence[i:] + sequence[:i]
        for sequence in (cores, cores[::-1])
        for i in range(len(cores))
    ]
    return "C" + "".join(min(necklaces))


def interiors(size):
    cycles = set()
    for tree in nx.nonisomorphic_trees(size):
        if max(dict(tree.degree()).values()) > 3:
            continue
        yield tree, False
        for u, v in nx.non_edges(tree):
            if max(tree.degree(u), tree.degree(v)) == 3:
                continue
            core = tree.copy()
            core.add_edge(u, v)
            key = canonical(core)
            if key not in cycles:
                cycles.add(key)
                yield core, True


def graphs(loops):
    if loops < 1:
        raise ValueError("at least one loop required")
    root = 2 * loops + 1
    for interior, cycle in interiors(root):
        num_spokes = loops if cycle else loops + 1
        if num_spokes < 2:
            continue
        required = {v for v, degree in interior.degree() if degree == 1}
        choices = {v for v, degree in interior.degree() if degree == 2}
        if not len(required) <= num_spokes <= len(required) + len(choices):
            continue
        seen = set()
        for extra in combinations(sorted(choices), num_spokes - len(required)):
            spokes = required.union(extra)
            key = canonical(interior, spokes)
            if key in seen:
                continue
            seen.add(key)
            full = interior.copy()
            full.add_edges_from((root, v) for v in spokes)
            if sparse_ok(full, 3):
                yield full, root, key


if __name__ == "__main__":
    from collections import Counter

    for n in range(1, 5):
        print(n, dict(Counter(code[0] for _, _, code in graphs(n))), flush=True)
