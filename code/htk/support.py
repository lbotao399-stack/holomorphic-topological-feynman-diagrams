"""Exact graph-rank bounds for holomorphic Schwinger numerator support.

For P = U**(n+H) D and every nonempty edge set S,
    degree_S(P) <= n * rank(C[:, S]) - 1,
    rank(C[:, S]) = len(S) + 1 - components(G without S).
All source components are retained.  See the proof in
explore/circuit-general/series-fibres/DERIVATION.md.

Only flats of the cycle-column matroid are needed.  A closure is computed by
adding bridges of G without S, including isolated vertices in component
counts.  Enumeration has a hard budget; any returned partial constraint set
is still valid and may safely be used for interpolation.
"""
from dataclasses import dataclass

__all__ = ["NumeratorSupport", "cycle_column_rank", "numerator_support"]


def _components_and_bridges(G, removed):
    adjacency = [[] for _ in range(G.nv)]
    for e, (u, v) in enumerate(G.edges):
        if removed & (1 << e):
            continue
        adjacency[u].append((v, e))
        adjacency[v].append((u, e))
    entry = [-1] * G.nv
    low = [-1] * G.nv
    bridges = 0
    clock = 0

    def visit(vertex, parent_edge):
        nonlocal clock, bridges
        entry[vertex] = low[vertex] = clock
        clock += 1
        for following, edge in adjacency[vertex]:
            if edge == parent_edge:
                continue
            if entry[following] < 0:
                visit(following, edge)
                low[vertex] = min(low[vertex], low[following])
                if low[following] > entry[vertex]:
                    bridges |= 1 << edge
            else:
                low[vertex] = min(low[vertex], entry[following])

    components = 0
    for vertex in range(G.nv):
        if entry[vertex] < 0:
            components += 1
            visit(vertex, -1)
    return components, bridges


def cycle_column_rank(G, edges):
    """Rank of the indicated columns of any full cycle matrix of connected G."""
    selected = set(edges)
    if any(not isinstance(e, int) or not 0 <= e < G.ne for e in selected):
        raise ValueError("edge IDs must be integers in range")
    if not G.is_connected():
        raise ValueError("a connected graph is required")
    mask = sum(1 << e for e in selected)
    components, _ = _components_and_bridges(G, mask)
    return len(selected) + 1 - components


@dataclass(frozen=True)
class NumeratorSupport:
    """Constraints use original edge IDs, excluding the fixed gauge edge.

    Each (edges, bound) means sum of exponents on edges <= bound.
    ``complete`` records whether all matroid flats were enumerated; an
    incomplete list loses only possible optimizations, never correctness.
    """
    gauge_edge: int
    constraints: tuple
    flat_count: int
    complete: bool
    zero: bool


def numerator_support(G, gauge_edge=0, *, max_flats=4096):
    """Find certified subset bounds without enumerating a full exponent cube.

    Individual bounds n-1 and total degree n*n-1 are implicit in the caller's
    lower set.  Redundant flat inequalities are omitted here.  The algorithm
    terminates at max_flats rather than silently spending exponential time.
    """
    if not isinstance(max_flats, int) or max_flats < 1:
        raise ValueError("max_flats must be a positive integer")
    if not isinstance(gauge_edge, int) or not 0 <= gauge_edge < G.ne:
        raise ValueError("gauge_edge must index an edge")
    if G.loops < 1 or not G.is_connected():
        raise ValueError("a connected graph with a loop is required")
    ground = ((1 << G.ne) - 1) ^ (1 << gauge_edge)
    _, initial_bridges = _components_and_bridges(G, 0)
    if initial_bridges:
        return NumeratorSupport(gauge_edge, (), 1, True, True)
    n = G.loops
    seen = {0}
    queue = [0]
    constraints = {}
    complete = True
    for flat in queue:
        components, _ = _components_and_bridges(G, flat)
        rank = flat.bit_count() + 1 - components
        bound = n * rank - 1
        if flat and bound < min((n - 1) * flat.bit_count(), n * n - 1):
            constraints[flat] = bound
        missing = ground ^ flat
        while missing:
            bit = missing & -missing
            missing ^= bit
            candidate = flat | bit
            _, bridges = _components_and_bridges(G, candidate)
            closure = candidate | (bridges & ground)
            if closure in seen:
                continue
            if len(seen) >= max_flats:
                complete = False
                continue
            seen.add(closure)
            queue.append(closure)
    ordered = tuple((tuple(e for e in range(G.ne) if mask & (1 << e)), bound)
                    for mask, bound in sorted(constraints.items(),
                                              key=lambda item: (item[0].bit_count(), item[0])))
    return NumeratorSupport(gauge_edge, ordered, len(seen), complete, False)
