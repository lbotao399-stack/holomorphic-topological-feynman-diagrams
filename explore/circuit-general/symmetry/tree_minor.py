"""Bounded tree/basis planning for the exact source-complete linear minor.

The current_minor determinant and CoreCurrentPeriods arithmetic are reused
unchanged.  A spanning tree selects BOTH a particular source flow and its
positively oriented fundamental cycle basis.  This is an exact coordinate
change, not a source restriction.  The original edge order is never changed.

A small linear-matrix support score ranks candidates; it is a heuristic, not
an optimality assertion.  Each attempted determinant has the caller's exact
polynomial term bound.  Failure tries another complete plan, and exhausting
all plans raises rather than returning a truncated or zero polynomial.
"""
from fractions import Fraction
from itertools import combinations, islice
from math import comb, lcm, prod
from numbers import Integral
from types import SimpleNamespace

from flint import fmpq, fmpq_mpoly_ctx

from htk import brute
from htk.forest_periods import CoreCurrentPeriods, ForestPeriodBudgetExceeded, _research_module
from htk.treeflow import tree_flow


def _integer(value, name, minimum=1):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError('%s must be an integer >= %d' % (name, minimum))
    return int(value)


def _is_tree(graph, tree):
    if len(tree) != graph.nv - 1 or len(set(tree)) != len(tree):
        return False
    parent = list(range(graph.nv))
    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex
    for edge in tree:
        u, v = map(find, graph.edges[edge])
        if u == v:
            return False
        parent[u] = v
    return True


def prepare_source_numerator(core, *, tree_candidates=None, max_tree_subsets=1000,
                             max_plans=8):
    """Populate a CoreCurrentPeriods numerator cache using bounded exact plans.

    The return value is a deterministic statistics dictionary.  The complete
    rational numerator is ``core._bare`` and remains compatible with all its
    original-current and series-chain methods.  No period worker is opened.
    Only connected tight H=2 graphs with at least one loop are used.

    With no explicit candidates, inspect at most ``max_tree_subsets`` edge
    subsets, always also including the default tree.  Partial tree scans are
    valid: they affect only the chance of finding a cheap coordinate system.
    With explicit candidates, reject invalid trees and reject a list longer
    than that bound.  Try at most ``max_plans`` ranked trees, each under the
    existing ``core.max_terms`` and ``core.max_matrix_entries`` bounds.
    ``core.check`` controls independent exact Kirchhoff density checks.
    """
    if not isinstance(core, CoreCurrentPeriods):
        raise TypeError('core must be a CoreCurrentPeriods instance')
    max_tree_subsets = _integer(max_tree_subsets, 'max_tree_subsets')
    max_plans = _integer(max_plans, 'max_plans')
    if core._bare is not None:
        return {'cached': True, 'numerator_terms': len(core._bare)}
    graph = core.G
    n = graph.loops
    if n < 1:
        raise ValueError('the linear-minor planner requires at least one loop')
    if (n * n - 1) ** 2 > core.max_matrix_entries:
        raise ForestPeriodBudgetExceeded('linear determinant exceeds max_matrix_entries')
    module = _research_module('explore/nonSP/verify/k4_periods/current_minor.py',
                              '_htk_planned_current_minor')
    scales = [lcm(*(row[i].denominator for row in core.sources)) for i in range(2)]
    scaled = [[int(row[i] * scales[i]) for i in range(2)] for row in core.sources]
    context = fmpq_mpoly_ctx.get(tuple('t%d' % e for e in range(graph.ne)))
    denominator = prod(scale ** n for scale in scales)
    candidates = {tuple(graph.spanning_tree())}
    scanned = 0
    if tree_candidates is None:
        for tree in islice(combinations(range(graph.ne), graph.nv - 1), max_tree_subsets):
            scanned += 1
            if _is_tree(graph, tree):
                candidates.add(tree)
        exhaustive = scanned == comb(graph.ne, graph.nv - 1)
    else:
        candidates.clear()
        for tree in islice(iter(tree_candidates), max_tree_subsets + 1):
            scanned += 1
            if scanned > max_tree_subsets:
                raise ForestPeriodBudgetExceeded('explicit trees exceed max_tree_subsets')
            tree = tuple(_integer(edge, 'tree edge', 0) for edge in tree)
            if any(edge >= graph.ne for edge in tree) or not _is_tree(graph, tree):
                raise ValueError('an explicit candidate is not a spanning tree')
            candidates.add(tuple(sorted(tree)))
        if not candidates:
            raise ValueError('at least one explicit tree candidate is required')
        exhaustive = False
    ranked = []
    for tree in sorted(candidates):
        score = []
        def capture(matrix, ctx):
            score.extend((sum(len(value) for row in matrix for value in row),
                          sum(bool(value) for row in matrix for value in row)))
            return ctx.constant(0)
        flow = [[int(value) for value in row] for row in tree_flow(graph, scaled, tree)]
        proxy = SimpleNamespace(n=n, N=graph.ne, C=graph.cycle_matrix(tree))
        module.numerator(proxy, scaled, flow, determinant=capture)
        # No callback means that the constant source/cycle columns are
        # rank-deficient, proving the entire numerator is exactly zero.
        ranked.append((tuple(score) if score else (0, 0), tree, not score))
    ranked.sort()
    failed = []
    for score, tree, certified_zero in ranked[:max_plans]:
        flow = [[int(value) for value in row] for row in tree_flow(graph, scaled, tree)]
        proxy = SimpleNamespace(n=n, N=graph.ne, C=graph.cycle_matrix(tree))
        try:
            polynomial = module.numerator(proxy, scaled, flow, determinant=core._determinant)
            core._bounded(polynomial)
        except ForestPeriodBudgetExceeded:
            failed.append(tree)
            continue
        value = context.from_dict({tuple(map(int, exponent)): fmpq(int(coefficient), denominator)
                                   for exponent, coefficient in polynomial.terms()})
        for index in range(core.check):
            point = [Fraction(edge + index + 2, edge + 1) for edge in range(graph.ne)]
            U = brute.symanzik_U(graph, point)
            expected = brute.density(graph, point, core.sources).rational() * U ** (n + 2)
            observed = Fraction(str(value(*(fmpq(t.numerator, t.denominator) for t in point))))
            if observed != expected:
                raise ArithmeticError('planned linear numerator failed exact original-density check')
        # Do not install an incomplete or independently inconsistent cache.
        core._context, core._bare = context, value
        return {'cached': False, 'scanned_subsets': scanned,
                'distinct_valid_trees': len(candidates), 'exhaustive_tree_scan': exhaustive,
                'max_plans': max_plans, 'attempted_plans': len(failed) + 1,
                'budget_failed_trees': failed, 'chosen_tree': tree,
                'linear_matrix_terms': score[0], 'linear_matrix_nonzero_entries': score[1],
                'constant_rank_zero': certified_zero,
                'source_component_scales': scales, 'numerator_terms': len(value),
                'independent_density_checks': core.check,
                'max_terms_per_attempt': core.max_terms,
                'max_matrix_entries': core.max_matrix_entries}
    raise ForestPeriodBudgetExceeded('all %d planned trees exceeded max_terms=%d; no numerator cached' %
                                     (min(len(ranked), max_plans), core.max_terms))
