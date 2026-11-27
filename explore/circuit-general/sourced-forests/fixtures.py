"""Original graphs, source assignments and words for the sourced-forest checks."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))
from htk import Graph


def sourced_fork():
    graph = Graph(5, [(0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)],
                  "sourced-three-loop-fork")
    sources = [(1, 2), (3, -1), (0, 0), (5, 3), (-9, -4)]
    return graph, sources, (3, 4)


def paired_higher_five():
    graph = Graph(7, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (2, 6),
                      (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)],
                  "paired-higher-valence-five-loop-forest")
    sources = [(1, 2), (3, -1), (-2, 4), (0, 0), (0, 0), (5, 3), (-7, -8)]
    return graph, sources, (5, 6)


def irreducible_quartic_five():
    # Vertex0 has four contracted fields and no surviving field. The two
    # degree-two cubic sources are1,2; vertices3,4 are source-free.
    graph = Graph(7, [(0, 1), (1, 3), (0, 2), (2, 4), (0, 5), (0, 6),
                      (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)],
                  "irreducible-quartic-five-loop-product-forest")
    sources = [(0, 0), (1, 2), (3, -1), (0, 0), (0, 0), (5, 3), (-9, -4)]
    return graph, sources, (5, 6)


def nonplanar_quartic_six():
    graph = Graph(8, [(0, 1), (1, 3), (0, 2), (2, 4), (0, 5), (0, 6),
                      (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7)],
                  "nonplanar-quartic-six-loop-product-forest")
    sources = [(0, 0), (1, 2), (3, -1), (0, 0), (0, 0), (0, 0), (5, 3), (-9, -4)]
    return graph, sources, (6, 7)


def higher_paired_comb(pairs):
    if isinstance(pairs, bool) or not isinstance(pairs, int) or pairs < 1:
        raise ValueError("pairs must be a positive integer")
    edges = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 3), (2, 4), (3, 4)]
    sources = [(1, 2), (3, -1), (-2, 4), (5, 3), (-7, -8)]
    parent = 0
    for _ in range(pairs):
        edges.remove((parent, 3))
        edges.remove((parent, 4))
        first, second = len(sources), len(sources) + 1
        sources += [(0, 0), (0, 0)]
        edges += [(first, parent), (first, 3), (first, 4),
                  (second, parent), (second, 3), (second, 4)]
        parent = first
    word = {(edges.index((parent, 3)), 0): 1}
    return Graph(len(sources), edges, "sourced-paired-comb-%d" % pairs), sources, (3, 4), word
