"""Exact series-chain contraction of the full holomorphic current determinant.

The scalar density is the one in htk/CONVENTIONS.md, Section 2 (T = 0).
Every source and every original edge current is retained.  Bivalent chains
produce constant source columns; these are eliminated once, before any
Schwinger-point evaluation.  A graph with m maximal chains needs an (m-1)
square determinant, instead of the original (ne-1) square determinant.

This is an algebraic contraction, not a period evaluator.  ``fibre_density``
provides the exact polynomial in chain fractions at fixed core resistances,
including the positive series-coordinate Jacobian.  It is suitable for exact
simplex cubature.  Arbitrary finite current monomials are supported.
"""
from dataclasses import dataclass
from fractions import Fraction as F

from . import brute
from .linalg import det, inverse, perm_sign
from .treeflow import tree_flow

__all__ = ["SeriesChain", "series_chains", "SeriesReduction"]


@dataclass(frozen=True)
class SeriesChain:
    """A maximal path through bivalent vertices, in traversal order.

    ``signs[j]`` compares traversal with the stored orientation of ``edges[j]``.
    The endpoints may coincide (a closed chain).  A pure polygon is cut at
    vertex zero, so it is represented by one closed chain.
    """

    edges: tuple
    vertices: tuple
    signs: tuple

    @property
    def length(self):
        return len(self.edges)


def series_chains(G):
    """Partition all original edge IDs into deterministic maximal chains.

    Parallel edges, reversed orientations and pure cycles are supported.
    The graph must be connected and have at least one edge.
    """
    if not G.ne or not G.is_connected():
        raise ValueError("series contraction needs a connected graph with edges")
    adjacency = [[] for _ in range(G.nv)]
    for e, (u, v) in enumerate(G.edges):
        adjacency[u].append((e, v))
        adjacency[v].append((e, u))
    anchors = {v for v in range(G.nv) if len(adjacency[v]) != 2}
    if not anchors:
        anchors = {0}
    used, chains = set(), []
    for start in sorted(anchors):
        for first, _ in adjacency[start]:
            if first in used:
                continue
            vertex, edge = start, first
            edges, vertices, signs = [], [start], []
            while True:
                if edge in used:
                    raise ValueError("chain traversal encountered an already used edge")
                used.add(edge)
                a, b = G.edges[edge]
                following = b if a == vertex else a
                edges.append(edge)
                signs.append(1 if a == vertex else -1)
                vertices.append(following)
                if following in anchors:
                    break
                next_edges = [(e, w) for e, w in adjacency[following] if e != edge]
                if len(next_edges) != 1:
                    raise ValueError("chain interior is not bivalent")
                vertex, edge = following, next_edges[0][0]
            chains.append(SeriesChain(tuple(edges), tuple(vertices), tuple(signs)))
    if len(used) != G.ne:
        raise ValueError("chain partition did not cover all edges")
    return tuple(chains)


def _contract_constant_columns(columns, dimension):
    """Factor det[columns, variable] by exact, recorded row elimination.

    Returns (rank, multiplier, remaining rows of the elimination matrix).
    If the fixed columns are dependent, every completed determinant is zero.
    No source direction or current component is dropped from the identity.
    """
    count = len(columns)
    rows = [[F(columns[j][i]) for j in range(count)] for i in range(dimension)]
    transform = [[F(int(i == j)) for j in range(dimension)] for i in range(dimension)]
    multiplier = F(1)
    for column in range(count):
        pivot = next((i for i in range(column, dimension) if rows[i][column]), None)
        if pivot is None:
            # The determinant is already identically zero, but later fixed
            # columns may still increase the matrix rank.  Report the true
            # rank rather than only the length of the first independent run.
            rank = column
            for later in range(column + 1, count):
                pivot = next((i for i in range(rank, dimension) if rows[i][later]), None)
                if pivot is None:
                    continue
                rows[rank], rows[pivot] = rows[pivot], rows[rank]
                value = rows[rank][later]
                for i in range(rank + 1, dimension):
                    if rows[i][later]:
                        factor = rows[i][later] / value
                        for j in range(later, count):
                            rows[i][j] -= factor * rows[rank][j]
                rank += 1
            return rank, F(0), ()
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            transform[column], transform[pivot] = transform[pivot], transform[column]
            multiplier = -multiplier
        value = rows[column][column]
        multiplier *= value
        for i in range(column + 1, dimension):
            if not rows[i][column]:
                continue
            factor = rows[i][column] / value
            for j in range(column, count):
                rows[i][j] -= factor * rows[column][j]
            transform[i] = [a - factor * b for a, b in zip(transform[i], transform[column])]
    return count, multiplier, tuple(tuple(row) for row in transform[count:])


class SeriesReduction:
    """Compress a T=0 current determinant, with exact rational sources.

    ``lam`` has the conventions of ``htk.brute.sources``. ``gauge`` is a chain
    index: the determinant uses the projective covector which equals one on
    that chain and zero on the other edges.  It divides by that chain's total
    resistance, so all gauge choices give exactly the same ``density(t)``.

    Public data: ``chains``, ``fixed_columns``, ``fixed_rank``,
    ``determinant_size`` and ``zero``.  ``zero`` certifies dependence of the
    constant source columns; other source-dependent zeros can still occur.
    ``fibre_degree(powers)`` is the safe total-degree bound m-1+|powers| in
    all chain-fraction variables together.  It does not concern core times.
    """

    def __init__(self, G, lam, gauge=0):
        self.G = G
        self.chains = series_chains(G)
        self.n = G.loops
        self.lam = brute.sources(G, lam)
        self.H = len(self.lam[0]) if self.lam else 0
        if any(len(row) != self.H for row in self.lam):
            raise ValueError("source vectors must have a common dimension")
        if self.H < 1 or self.n * self.H != G.ne - 1:
            raise ValueError("requires H >= 1 and n H = ne-1 (T = 0)")
        self.full_sources = self.lam + [[-sum(row[i] for row in self.lam)
                                        for i in range(self.H)]]
        if not isinstance(gauge, int) or not 0 <= gauge < len(self.chains):
            raise ValueError("gauge must index a maximal chain")
        self.gauge = gauge
        self.C = G.cycle_matrix()
        self.y0 = tree_flow(G, self.full_sources)
        self.core_cycles = []
        for chain in self.chains:
            column = tuple(chain.signs[0] * self.C[r][chain.edges[0]]
                           for r in range(self.n))
            if any(tuple(sign * self.C[r][edge] for r in range(self.n)) != column
                   for edge, sign in zip(chain.edges, chain.signs)):
                raise ArithmeticError("cycle columns disagree along a series chain")
            self.core_cycles.append(column)

        # Group columns in path order, then replace each group by its weighted
        # mean and successive differences.  That replacement has determinant 1.
        edge_order = [edge for chain in self.chains for edge in chain.edges]
        base_positions, fixed_positions = [], []
        position = 0
        fixed = []
        for chain, cycle in zip(self.chains, self.core_cycles):
            base_positions.append(position)
            for j, vertex in enumerate(chain.vertices[1:-1], start=1):
                fixed_positions.append(position + j)
                fixed.append(tuple(cycle[r] * self.full_sources[vertex][i]
                                   for i in range(self.H) for r in range(self.n)))
            position += chain.length
        omitted = base_positions[gauge]
        minor_positions = [i for i in range(G.ne) if i != omitted]
        ordered_positions = fixed_positions + [i for i in base_positions if i != omitted]
        minor_order = [minor_positions.index(i) for i in ordered_positions]
        self.column_sign = perm_sign(edge_order) * (-1) ** omitted * perm_sign(minor_order)
        self.fixed_columns = tuple(fixed)
        self.fixed_rank, multiplier, self.elimination_rows = _contract_constant_columns(
            fixed, G.ne - 1)
        self.determinant_factor = self.column_sign * multiplier
        self.zero = multiplier == 0
        self.determinant_size = len(self.chains) - 1
        self._core_key = None
        self._inverse = None
        self._U = None

    def _check_powers(self, powers):
        if powers is None:
            return {}
        out = dict(powers)
        for key, power in out.items():
            if (not isinstance(key, tuple) or len(key) != 2
                    or any(not isinstance(x, int) for x in key)
                    or not 0 <= key[0] < self.G.ne or not 0 <= key[1] < self.H
                    or not isinstance(power, int) or power < 0):
                raise ValueError("powers must map (edge, component) to nonnegative integers")
        return out

    def fibre_degree(self, powers=None):
        """Safe total polynomial degree at fixed positive core resistances."""
        return self.determinant_size + sum(self._check_powers(powers).values())

    def core_resistances(self, t):
        t = tuple(F(x) for x in t)
        if len(t) != self.G.ne or any(x < 0 for x in t):
            raise ValueError("need one nonnegative time per original edge")
        resistances = tuple(sum(t[e] for e in chain.edges) for chain in self.chains)
        if any(x <= 0 for x in resistances):
            raise ValueError("each chain must have positive total resistance")
        return resistances

    def fractions(self, t):
        t = tuple(F(x) for x in t)
        resistances = self.core_resistances(t)
        return tuple(tuple(t[e] / resistance for e in chain.edges)
                     for chain, resistance in zip(self.chains, resistances))

    def times(self, resistances, fractions):
        resistances = tuple(F(x) for x in resistances)
        fractions = tuple(tuple(F(x) for x in row) for row in fractions)
        if len(resistances) != len(self.chains) or len(fractions) != len(self.chains):
            raise ValueError("need one resistance and fraction tuple per chain")
        if any(x <= 0 for x in resistances):
            raise ValueError("core resistances must be positive")
        out = [F(0)] * self.G.ne
        for chain, resistance, row in zip(self.chains, resistances, fractions):
            if len(row) != chain.length or any(x < 0 for x in row) or sum(row) != 1:
                raise ValueError("chain fractions must be nonnegative and sum to one")
            for edge, fraction in zip(chain.edges, row):
                out[edge] = resistance * fraction
        return out

    def _core_inverse(self, resistances):
        key = tuple(resistances)
        if key != self._core_key:
            matrix = [[sum(resistance * cycle[r] * cycle[s]
                           for resistance, cycle in zip(key, self.core_cycles))
                       for s in range(self.n)] for r in range(self.n)]
            self._inverse = inverse(matrix)
            self._U = det(matrix)
            self._core_key = key
        return self._inverse, self._U

    def currents(self, t):
        """All original oriented edge currents; only an n by n loop solve.

        At fixed core resistances the inverse and U are cached, and the
        currents are affine in all chain fractions.  Zero fractions use the
        polynomial boundary value; no inverse zero edge time is formed.
        """
        t = tuple(F(x) for x in t)
        resistances = self.core_resistances(t)
        inverse_matrix, _ = self._core_inverse(resistances)
        rhs = [[-sum(self.C[r][e] * t[e] * self.y0[e][i] for e in range(self.G.ne))
                for i in range(self.H)] for r in range(self.n)]
        loop = [[sum(inverse_matrix[r][s] * rhs[s][i] for s in range(self.n))
                 for i in range(self.H)] for r in range(self.n)]
        return [[self.y0[e][i] + sum(self.C[r][e] * loop[r][i] for r in range(self.n))
                 for i in range(self.H)] for e in range(self.G.ne)]

    def density(self, t, powers=None):
        """Exact D(t) times a specified monomial in original edge currents."""
        powers = self._check_powers(powers)
        t = tuple(F(x) for x in t)
        resistances = self.core_resistances(t)
        if self.zero:
            return F(0)
        currents = self.currents(t)
        columns = []
        for j, (chain, cycle, resistance) in enumerate(
                zip(self.chains, self.core_cycles, resistances)):
            if j == self.gauge:
                continue
            mean = [sum(t[e] * sign * currents[e][i]
                        for e, sign in zip(chain.edges, chain.signs)) / resistance
                    for i in range(self.H)]
            column = [cycle[r] * mean[i] for i in range(self.H) for r in range(self.n)]
            columns.append([sum(a * b for a, b in zip(row, column))
                            for row in self.elimination_rows])
        rows = [[columns[j][i] for j in range(self.determinant_size)]
                for i in range(self.determinant_size)]
        result = (self.determinant_factor * det(rows)
                  / (resistances[self.gauge] * self._U ** self.H))
        for (edge, component), power in powers.items():
            result *= currents[edge][component] ** power
        return result

    def fibre_density(self, resistances, fractions, powers=None):
        """D(t) prod_chain R_chain^(length-1), times current insertions.

        Fractions use ordinary simplex Lebesgue measure (volume 1/(L-1)!).
        The returned coefficient uses the original graph's density sign.
        Fixing one core resistance to one gives a projective core chart.
        """
        t = self.times(resistances, fractions)
        result = self.density(t, powers)
        for chain, resistance in zip(self.chains, self.core_resistances(t)):
            result *= resistance ** (chain.length - 1)
        return result
