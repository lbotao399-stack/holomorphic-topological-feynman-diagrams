"""Ordered source-complete forests for the sourced-ear reduction controls."""
from fractions import Fraction as F
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))
from htk.graph import Graph


CORE_EDGES = ((0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
GENERIC_CORE_SOURCES = ((1, 2), (3, -1), (-2, 4), (5, 3), (-7, -8))
PHYSICAL_CORE_SOURCES = ((2, -4), (1, 2), (0, 0), (5, 3), (-8, -1))


def sourced_ear_forest(ears, *, physical=False):
    """A non-SP fork with arbitrary sourced cubic ears at its W--A edge.

    In the physical assignment, vertex1 has contracted degree ears+3 and
    one surviving field, so a superpotential of degree ears+4 supports it.
    Vertex0 and every added degree-two vertex are cubic vertices with one
    surviving field. Vertex2 is cubic with all three fields contracted and
    exactly zero source. No source is imposed to vanish merely to simplify
    the computational method.
    """
    if isinstance(ears, bool) or not isinstance(ears, int) or ears < 1:
        raise ValueError("ears must be a positive integer")
    rows = [list(map(F, row)) for row in (PHYSICAL_CORE_SOURCES if physical else GENERIC_CORE_SOURCES)]
    edges = list(CORE_EDGES)
    for i in range(ears):
        vertex = len(rows)
        source = (F(i + 2), F(3 * i - 3))
        edges.extend(((1, vertex), (vertex, 3)))
        rows.append(list(source))
        for a in range(2):
            rows[4][a] -= source[a]
    graph = Graph(len(rows), edges, "%s-sourced-ear-%d" % ("physical" if physical else "generic", ears))
    metadata = {
        "interaction_edges_form_a_tree": True,
        "higher_vertex": 1,
        "higher_vertex_contracted_degree": ears + 3,
        "higher_superpotential_degree": ears + 4 if physical else None,
        "higher_vertex_surviving_fields": 1 if physical else None,
        "cubic_vertex_contracted_and_surviving_fields":
            {0: (2, 1), 2: (3, 0), **{v: (2, 1) for v in range(5, 5 + ears)}} if physical else None,
        "physical_source_pattern": physical,
        "scope": "An individual graph kernel; no Wick/color/operator sum is asserted.",
    }
    return graph, tuple(map(tuple, rows)), (3, 4), metadata
