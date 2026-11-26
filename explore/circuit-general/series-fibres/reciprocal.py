"""Reciprocal core interpolation and rigorous integer-CRT certificates.

Let P=U**(n+H) D, and integrate P over all chain fractions, WITHOUT the
series Jacobian.  Its reciprocal core polynomial is
    Q(g) = prod_c g_c**(n-1) * integral P(t_cj=x_cj/g_c) dx.
It has total degree (n-1)*(V_core-2), individual degrees <= n-1, and
    degree_T Q <= n * graphic_rank(T) - len(T)
for every proper core edge set T.  Only proper graphic flats are needed.

The coefficient certificate uses integral tree-current coefficients,
permanent bounds for a determinant numerator, and the elementary
multivariate Mahler inequality.  It gives an unconditional CRT stopping
threshold, not a stabilization heuristic. See RECIPROCAL.md for proofs.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from math import comb, factorial, gcd
from operator import index
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ASSEMBLY = HERE.parent / "assembly"
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ASSEMBLY))

from htk.graph import Graph
from htk.linalg import det
from htk.series import SeriesReduction
from htk.treeflow import tree_flow
from htk import brute
from numerator import LowerSet

__all__ = ["ReciprocalPlan", "reciprocal_plan", "reciprocal_gauge_counts",
           "integer_certificate"]


def _integer(value, name, minimum=0):
    if isinstance(value, bool):
        raise TypeError("%s must be an integer" % name)
    try:
        value = index(value)
    except TypeError as error:
        raise TypeError("%s must be an integer" % name) from error
    if value < minimum:
        raise ValueError("%s must be at least %d" % (name, minimum))
    return value


def _lcm(a, b):
    return a // gcd(a, b) * b


def _core(reduction):
    if not isinstance(reduction, SeriesReduction):
        raise TypeError("reduction must be a SeriesReduction")
    endpoints = [(chain.vertices[0], chain.vertices[-1]) for chain in reduction.chains]
    if any(u == v for u, v in endpoints):
        raise ValueError("reciprocal core planning requires a loopless series core")
    vertices = sorted({v for edge in endpoints for v in edge})
    if len(vertices) < 2:
        raise ValueError("reciprocal core planning needs at least two core vertices")
    labels = {v: i for i, v in enumerate(vertices)}
    edges = tuple((labels[u], labels[v]) for u, v in endpoints)
    core = Graph(len(vertices), edges)
    if not core.is_connected() or core.loops != reduction.n:
        raise ArithmeticError("series core does not preserve connectivity and loop number")
    return core


def _graphic_closure(core, mask):
    parent = list(range(core.nv))

    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    rank = 0
    for edge, (u, v) in enumerate(core.edges):
        if mask & (1 << edge):
            u, v = find(u), find(v)
            if u != v:
                parent[u] = v
                rank += 1
    closure = sum(1 << edge for edge, (u, v) in enumerate(core.edges)
                  if find(u) == find(v))
    return closure, rank


def _core_constraints(core, max_subsets):
    if (1 << core.ne) > max_subsets:
        raise ValueError("complete reciprocal support enumeration exceeds max_subsets")
    n = core.loops
    degree = (n - 1) * (core.nv - 2)
    full = (1 << core.ne) - 1
    flats = {}
    for mask in range(1, full):
        closure, rank = _graphic_closure(core, mask)
        if closure == full:
            continue
        flats[closure] = rank
    constraints = []
    for mask, rank in flats.items():
        edges = tuple(edge for edge in range(core.ne) if mask & (1 << edge))
        bound = n * rank - len(edges)
        if bound < min((n - 1) * len(edges), degree):
            constraints.append((edges, bound))
    return tuple(sorted(constraints, key=lambda item: (len(item[0]), item[0]))), len(flats) + 2


def _dehomogenized_constraints(constraints, variables, per_edge, degree):
    location = {edge: axis for axis, edge in enumerate(variables)}
    result = {}
    for edges, bound in constraints:
        # Dropping the nonnegative gauge exponent preserves the inequality.
        selected = tuple(location[edge] for edge in edges if edge in location)
        if bound < min(per_edge * len(selected), degree):
            result[selected] = min(result.get(selected, bound), bound)
    return tuple(sorted(result.items(), key=lambda item: (len(item[0]), item[0])))


@dataclass(frozen=True)
class ReciprocalPlan:
    """Full exponents and core edges follow the chain order of the reduction.

    LowerSet axes follow ``variable_edges``. Some grid coefficients must be
    zero by gauge-containing homogeneous bounds; ``support_mask`` identifies
    these exactly. ``sample_points`` gives positive integer conductances,
    not resistances. The omitted conductance equals one.
    """

    lower_set: object
    gauge: int
    total_degree: int
    per_edge_degree: int
    core_edges: tuple
    core_vertices: int
    variable_edges: tuple
    chain_lengths: tuple
    constraints: tuple
    flat_count: int
    loop_number: int
    core_spanning_trees: int
    grid_laplacian_bound: int

    def exponents(self):
        exponents = np.zeros((self.lower_set.M, len(self.core_edges)), dtype=np.int64)
        exponents[:, self.variable_edges] = self.lower_set.dig
        exponents[:, self.gauge] = self.total_degree - exponents.sum(axis=1)
        return exponents

    def _exponents(self, exponents):
        if exponents is None:
            return self.exponents()
        exponents = np.asarray(exponents)
        if (exponents.ndim != 2 or exponents.shape[1] != len(self.core_edges)
                or not np.issubdtype(exponents.dtype, np.integer)):
            raise ValueError("exponents must be an integer array with one column per chain")
        return exponents

    def support_mask(self, exponents=None):
        exponents = self._exponents(exponents)
        good = ((exponents >= 0).all(axis=1)
                & (exponents <= self.per_edge_degree).all(axis=1)
                & (exponents.sum(axis=1) == self.total_degree))
        for selected, bound in self.constraints:
            good &= exponents[:, selected].sum(axis=1) <= bound
        return good

    def period_powers(self, exponents=None):
        """Original resistance-period powers, including the series Jacobian.

        Validate nonzero coefficients with support_mask before period lookup;
        this function preserves the alignment of the full interpolation grid.
        """
        exponents = self._exponents(exponents)
        return self.per_edge_degree - exponents + np.array(self.chain_lengths, dtype=np.int64)

    def sample_points(self, start=0, stop=None):
        start = _integer(start, "start")
        stop = self.lower_set.M if stop is None else _integer(stop, "stop")
        if not 0 <= start <= stop <= self.lower_set.M:
            raise ValueError("invalid interpolation sample range")
        points = np.ones((stop - start, len(self.core_edges)), dtype=np.int64)
        points[:, self.variable_edges] += self.lower_set.dig[start:stop].astype(np.int64)
        return points


def reciprocal_plan(reduction, *, gauge=None, max_points=2000000, max_subsets=65536):
    """Plan the complete reciprocal support for a connected loopless core.

    The graph need not be K5; it may be any loopless series core that fits
    the explicit support and grid budgets. This bare-numerator planner does
    not include current insertions. Exhausting a budget raises an exception.
    """
    max_points = _integer(max_points, "max_points", 1)
    max_subsets = _integer(max_subsets, "max_subsets", 1)
    core = _core(reduction)
    gauge = reduction.gauge if gauge is None else _integer(gauge, "gauge")
    if not 0 <= gauge < core.ne:
        raise ValueError("gauge must index a core edge")
    constraints, flat_count = _core_constraints(core, max_subsets)
    degree = (core.loops - 1) * (core.nv - 2)
    per_edge = core.loops - 1
    variables = tuple(edge for edge in range(core.ne) if edge != gauge)
    lowered = _dehomogenized_constraints(constraints, variables, per_edge, degree)
    lower_set = LowerSet(core.ne - 1, per_edge, degree,
                         constraints=lowered, max_points=max_points)
    tree_count = int(det(brute.laplacian(core, [1] * core.ne)))
    laplacian_bound = tree_count * (per_edge + 1) ** (core.nv - 1)
    return ReciprocalPlan(lower_set, gauge, degree, per_edge, tuple(core.edges), core.nv,
                          variables, tuple(chain.length for chain in reduction.chains),
                          constraints, flat_count, core.loops, tree_count, laplacian_bound)


def reciprocal_gauge_counts(reduction, *, max_points=2000000, max_subsets=65536):
    """Exact grid counts for all gauges, without retaining multiple grids."""
    core = _core(reduction)
    counts = []
    for gauge in range(core.ne):
        plan = reciprocal_plan(reduction, gauge=gauge,
                               max_points=max_points, max_subsets=max_subsets)
        counts.append(plan.lower_set.M)
    return tuple(counts)


def _simplex_denominator(dimension, degree, max_moment_terms):
    if dimension == 0:
        return 1
    terms = comb(degree + dimension, dimension)
    if terms > max_moment_terms:
        raise ValueError("exact simplex denominator enumeration exceeds max_moment_terms")
    denominator = 1
    for positions in combinations_with_replacement(range(dimension + 1), degree):
        powers = [0] * (dimension + 1)
        for position in positions:
            powers[position] += 1
        numerator = 1
        for power in powers:
            numerator *= factorial(power)
        moment_denominator = factorial(degree + dimension) // numerator
        denominator = _lcm(denominator, moment_denominator)
    return denominator


def _tree_current_norms(G, sources, max_tree_candidates):
    candidates = comb(G.ne, G.nv - 1)
    if candidates > max_tree_candidates:
        raise ValueError("tree-current certificate exceeds max_tree_candidates")
    H = len(sources[0])
    heights = [[0] * H for _ in G.edges]
    count = 0
    for tree in combinations(range(G.ne), G.nv - 1):
        try:
            currents = tree_flow(G, sources, tree)
        except ValueError:
            continue
        count += 1
        for edge in tree:
            for component in range(H):
                value = currents[edge][component]
                if value.denominator != 1:
                    raise ArithmeticError("scaled tree current is not integral")
                heights[edge][component] += abs(value.numerator)
    if not count:
        raise ArithmeticError("connected graph has no spanning tree")
    return heights, count, candidates


def _permanent(matrix):
    size = len(matrix)
    values = [0] * (1 << size)
    values[0] = 1
    for mask in range(1, len(values)):
        row = mask.bit_count() - 1
        remaining, value = mask, 0
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            column = bit.bit_length() - 1
            if matrix[row][column]:
                value += matrix[row][column] * values[mask ^ bit]
        values[mask] = value
    return values[-1]


def integer_certificate(reduction, *, max_tree_candidates=1000000,
                        max_permanent_size=20, max_moment_terms=100000):
    """Proven coefficient denominator and height for the bare reciprocal Q.

    Returns a JSON-serializable dictionary. ``denominator * Q`` has integer
    coefficients bounded in absolute value by ``coefficient_bound``. CRT
    reconstruction is unique once its modulus EXCEEDS ``crt_modulus_must_exceed``.
    Rational sources are supported by one common scaling denominator, whose
    homogeneous power n*H is included in the returned denominator.
    """
    max_tree_candidates = _integer(max_tree_candidates, "max_tree_candidates", 1)
    max_permanent_size = _integer(max_permanent_size, "max_permanent_size", 1)
    max_moment_terms = _integer(max_moment_terms, "max_moment_terms", 1)
    core = _core(reduction)
    G, n, H = reduction.G, reduction.n, reduction.H
    if G.ne - 1 > max_permanent_size:
        raise ValueError("determinant height certificate exceeds max_permanent_size")
    source_denominator = 1
    for row in reduction.full_sources:
        for value in row:
            source_denominator = _lcm(source_denominator, F(value).denominator)
    scaled = [[int(F(value) * source_denominator) for value in row]
              for row in reduction.full_sources]
    simplex_denominators = [_simplex_denominator(chain.length - 1, n - 1, max_moment_terms)
                            for chain in reduction.chains]
    moment_denominator = 1
    for denominator in simplex_denominators:
        moment_denominator *= denominator
    heights, trees, candidates = _tree_current_norms(G, scaled, max_tree_candidates)
    permanents = []
    for gauge in range(G.ne):
        columns = [edge for edge in range(G.ne) if edge != gauge]
        matrix = [[abs(reduction.C[row][edge]) * heights[edge][component]
                   for edge in columns] for component in range(H) for row in range(n)]
        permanents.append(_permanent(matrix))
    gauge = min(range(G.ne), key=permanents.__getitem__)
    mahler_factor = 1 << (G.ne * (n - 1))
    bound = moment_denominator * mahler_factor * permanents[gauge]
    core_trees = int(det(brute.laplacian(core, [1] * core.ne)))
    return {
        "scope": "bare T=0 reciprocal core polynomial; complete sources retained",
        "denominator": moment_denominator * source_denominator ** (n * H),
        "moment_denominator": moment_denominator,
        "simplex_denominators": simplex_denominators,
        "source_denominator": source_denominator,
        "source_homogeneous_degree": n * H,
        "scaled_integer_sources": scaled,
        "tree_candidates": candidates,
        "spanning_trees": trees,
        "tree_current_coefficient_l1": heights,
        "determinant_size": G.ne - 1,
        "determinant_permanent_by_gauge": permanents,
        "height_gauge_edge": gauge,
        "mahler_factor": mahler_factor,
        "coefficient_bound": bound,
        "coefficient_bound_bits": bound.bit_length(),
        "crt_modulus_must_exceed": 2 * bound,
        "core_spanning_trees": core_trees,
        "grid_conductance_max": n,
        "grid_laplacian_bound": core_trees * n ** (core.nv - 1),
        "proof": "integer polynomial divisibility; exact Dirichlet denominators; positive moment contraction; multivariate Mahler coefficient bound; permanent of exact tree-current coefficient norms",
    }
