"""Independent planar-face dual and prism-table replay of the physical five-loop core.

The face construction and direct table lookup do not import per5 or its
K5-minus-an-edge duality implementation. The established current-minor
formula provides the exact original-graph polynomial.
"""

from collections import defaultdict
from fractions import Fraction as F
import importlib.util
from itertools import combinations
import json
from math import factorial
from pathlib import Path
from types import SimpleNamespace
import sys

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))
from htk import Graph
from htk.treeflow import tree_flow

EDGES = ((0, 1), (1, 3), (0, 2), (2, 4), (0, 5), (0, 6),
         (3, 5), (3, 6), (4, 5), (4, 6), (5, 6))
SOURCES = ((0, 0), (1, 2), (3, -1), (0, 0), (0, 0), (5, 3), (-9, -4))
CORE = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 3), (1, 4),
        (2, 3), (2, 4), (3, 4))
EXPECTED = {"1": F(3327759589, 46656000), "z3": F(50942017, 324000),
            "z5": -F(45619, 180)}


def spanning_tree_masks(edges):
    vertices = len({vertex for edge in edges for vertex in edge})
    return {sum(1 << edge for edge in chosen)
            for chosen in combinations(range(len(edges)), vertices - 1)
            if nx.is_tree(nx.Graph([edges[edge] for edge in chosen]))}


def exact_original_numerator():
    path = ROOT / "explore/nonSP/verify/k4_periods/current_minor.py"
    spec = importlib.util.spec_from_file_location("_independent_sourced_current_minor", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    G = Graph(7, EDGES)
    proxy = SimpleNamespace(n=G.loops, N=G.ne, C=G.cycle_matrix())
    flow = [[int(value) for value in row] for row in tree_flow(G, SOURCES)]
    polynomial = module.numerator(proxy, SOURCES, flow)
    return {tuple(map(int, exponent)): int(coefficient)
            for exponent, coefficient in polynomial.terms()}


def replay(numerator=None, q=7, expected=None):
    """Integrate a complete original-graph polynomial on the fixed nine-edge core."""
    if numerator is None:
        assert q == 7
        numerator = exact_original_numerator()
    if expected is None and q == 7:
        expected = EXPECTED
    planar, embedding = nx.check_planarity(nx.Graph(CORE))
    assert planar
    halfedge_face, faces = {}, []
    for left in embedding:
        for right in embedding[left]:
            if (left, right) not in halfedge_face:
                face = embedding.traverse_face(left, right)
                index = len(faces)
                for i in range(len(face)):
                    halfedge_face[face[i], face[(i + 1) % len(face)]] = index
                faces.append(face)
    dual = [(halfedge_face[left, right], halfedge_face[right, left]) for left, right in CORE]
    directory = "assembly" if q == 8 else "cores"
    data = json.loads((ROOT / ("explore/circuit-general/%s/tables/prism_q%d.json"
                              % (directory, q))).read_text())
    prism = [tuple(edge) for edge in data["edges"]]
    mapping = next(nx.isomorphism.GraphMatcher(nx.Graph(dual), nx.Graph(prism)).isomorphisms_iter())
    positions = {frozenset(edge): index for index, edge in enumerate(prism)}
    edge_map = [positions[frozenset((mapping[left], mapping[right]))] for left, right in dual]
    assert sorted(edge_map) == list(range(9))

    # This coefficient identity proves U_G(1/x)*prod(x)=U_dual(x)
    # for every positive resistance assignment, without point sampling.
    primal_trees = spanning_tree_masks(CORE)
    dual_trees = spanning_tree_masks(dual)
    assert primal_trees == {511 ^ tree for tree in dual_trees}
    assert len(primal_trees) == 75

    automorphisms = []
    for image in nx.isomorphism.GraphMatcher(nx.Graph(prism), nx.Graph(prism)).isomorphisms_iter():
        automorphisms.append(tuple(positions[frozenset((image[left], image[right]))]
                                   for left, right in prism))
    table = {tuple(powers): {label: F(value) for label, value in values.items()}
             for powers, values in data["values"]}
    grouped, core_vectors = defaultdict(F), set()
    for exponent, coefficient in numerator.items():
        powers = ((exponent[0] + exponent[1] + 2, exponent[2] + exponent[3] + 2)
                  + tuple(value + 1 for value in exponent[4:]))
        assert sum(powers) == 5 * q and all(0 < power < q for power in powers)
        core_vectors.add(powers)
        weight = F(coefficient * factorial(exponent[0]) * factorial(exponent[1])
                   * factorial(exponent[2]) * factorial(exponent[3]),
                   factorial(exponent[0] + exponent[1] + 1)
                   * factorial(exponent[2] + exponent[3] + 1))
        transformed = [None] * 9
        for old, new in enumerate(edge_map):
            transformed[new] = q - powers[old]
        representative = min(tuple(transformed[i] for i in image) for image in automorphisms)
        grouped[representative] += weight

    # Check every strict edge-subset condition on every distinct core
    # weight vector before any orbit coefficient cancels.
    vectors = np.asarray(sorted(core_vectors), dtype=np.int64)
    minimum_margin = None
    for mask in range(1, 511):
        chosen = [index for index in range(9) if mask & (1 << index)]
        subgraph = nx.Graph([CORE[index] for index in chosen])
        cycles = len(chosen) - len(subgraph) + nx.number_connected_components(subgraph)
        margin = int((vectors[:, chosen].sum(axis=1) - q * cycles).min())
        assert margin > 0
        minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)

    value, orbits = defaultdict(F), 0
    for representative, coefficient in grouped.items():
        if coefficient:
            orbits += 1
            for label, period in table[representative].items():
                value[label] += coefficient * period
    value = {label: coefficient for label, coefficient in value.items() if coefficient}
    if expected is not None:
        assert value == expected
    return value, {"numerator_terms": len(numerator), "core_period_orbits": orbits,
                   "spanning_tree_monomials": len(primal_trees),
                   "minimum_convergence_margin": minimum_margin,
                   "independent_dual_edge_map": edge_map}


if __name__ == "__main__":
    value, details = replay()
    print(details)
    print({label: str(coefficient) for label, coefficient in value.items()})
