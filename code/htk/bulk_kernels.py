"""One exact evaluator for admissible bulk kernels in any (H,T).

Routes: exact bulk vanishing, complete current/source polynomials when
available, series-parallel integration, general polynomial numerator plus
read-only scalar periods, and a complete convergent rational series.
Sources and finite original holomorphic-current words are retained in all
routes.  Zero topological source is required; defects are outside this API.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import factorial, prod
from pathlib import Path
from types import MappingProxyType, SimpleNamespace
import json
import subprocess
import sys
import weakref

from . import brute
from .bulk_series import ExactBulkSeries
from .core_sector_periods import _cycle_count
from .forest_cubes import _integer, _rational
from .forest_periods import (_PeriodWorker, _stop_process, MissingForestPeriod,
                             ForestPeriodBudgetExceeded)
from .graph import Graph
from .series import series_chains
from .simplex import PiNumber
from .spred import sp_reduce, NotSeriesParallel
from .topological_polynomial import BulkPolynomial, rational, source_vectors

__all__ = ["evaluate_bulk", "BulkKernelValue"]


@dataclass(frozen=True)
class BulkKernelValue:
    coefficients: object
    route: str
    evidence: object

    def as_sympy(self):
        import sympy as sp
        result = sp.Integer(0)
        for label, coefficient in self.coefficients.items():
            term = sp.Rational(coefficient.numerator, coefficient.denominator)
            if label != "1":
                for factor in label.split("*"):
                    if factor.startswith("sqrtpi^"):
                        term *= sp.pi**sp.Rational(int(factor[7:]), 2)
                    elif factor.startswith("pi^"):
                        term *= sp.pi**int(factor[3:])
                    elif factor.startswith("z") and factor[1:].isdigit():
                        term *= sp.zeta(int(factor[1:]))
                    else:
                        raise ValueError("unknown exact scalar label: " + factor)
            result += term
        return result


def _value(coefficients, route, **evidence):
    return BulkKernelValue(MappingProxyType({k: rational(v) for k, v in coefficients.items() if v}),
                           route, MappingProxyType(evidence))


class _BulkPeriodWorker(_PeriodWorker):
    def __init__(self, timeout):
        self.timeout = timeout
        path = Path(__file__).with_name("_bulk_period_worker.py")
        self.process = subprocess.Popen([sys.executable, "-u", str(path)], stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                        text=True, bufsize=1)
        self._finalizer = weakref.finalize(self, _stop_process, self.process)


def _saved_source_polynomial(graph, sources, word, H, T):
    if (H, T, graph.nv, graph.ne) != (1, 1, 6, 9) or len(word) > 1:
        return None
    import networkx as nx
    if len(set(tuple(sorted(e)) for e in graph.edges)) != graph.ne:
        return None
    root = Path(__file__).resolve().parents[2]
    for name in ("prism", "K33"):
        path = root / "explore/circuit-general/topological-completion" / (name+"_complete_polynomials.json")
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        canonical = [tuple(e) for e in data["edges"]]
        mapping = next(nx.algorithms.isomorphism.GraphMatcher(
            nx.Graph(graph.edges), nx.Graph(canonical)).isomorphisms_iter(), None)
        if mapping is None:
            continue
        indices = {tuple(sorted(e)): i for i, e in enumerate(canonical)}
        edge_map = [indices[tuple(sorted((mapping[a], mapping[b])))] for a, b in graph.edges]
        orientation = (-1)**sum(a > b for i, a in enumerate(edge_map) for b in edge_map[i+1:])
        full = [row[0] for row in sources] + [-sum(row[0] for row in sources)]
        canonical_source = [Fraction(0)]*6
        for old, new in mapping.items():
            canonical_source[new] = full[old]
        output = 0
        if word:
            edge, component = word[0]
            output = 1+edge_map[edge]
            a, b = graph.edges[edge]
            if (mapping[a], mapping[b]) != canonical[edge_map[edge]]:
                orientation *= -1
        coefficients = {}
        for exponent, label, scalar in data["bare_and_all_nine_current_polynomials"][output]:
            coefficient = orientation*Fraction(scalar)*prod(x**p for x, p in zip(canonical_source[:5], exponent))
            coefficients[label] = coefficients.get(label, Fraction(0)) + coefficient
        return _value(coefficients, "complete source polynomial", core=name,
                      all_original_sources_retained=True, all_original_current_slots_retained=True)
    return None


def _scalar_components(density, word):
    numerator, power = density.inserted(word)
    graph = density.graph
    chains = series_chains(graph)
    anchors = sorted({chain.vertices[0] for chain in chains} | {chain.vertices[-1] for chain in chains})
    labels = {v: i for i, v in enumerate(anchors)}
    edges = [(labels[chain.vertices[0]], labels[chain.vertices[-1]]) for chain in chains]
    coefficients = {}
    for exponent, coefficient in numerator.to_dict().items():
        alpha = tuple(sum(int(exponent[e])+1 for e in chain.edges) for chain in chains)
        beta_factor = Fraction(1)
        for chain, total in zip(chains, alpha):
            beta_factor *= Fraction(prod(factorial(int(exponent[e])) for e in chain.edges), factorial(total-1))
        coefficients[alpha] = coefficients.get(alpha, Fraction(0)) + rational(coefficient)*beta_factor
    return edges, {a: c for a, c in coefficients.items() if c}, power, chains


def _check_admissible(graph, dimension):
    """All Laman inequalities, after exact contraction of bivalent chains.

    An offending subgraph can be stripped of its degree-one vertices;
    the remaining subgraph consists of whole maximal chains.  For a
    connected subgraph its inequality is E >= dimension*loops+1.  Testing
    every chain subset includes each connected component separately.
    A bridgeless graph of at most five loops has at most twelve chains
    unless it is a polygon (one chain), irrespective of dimension.
    """
    cycles = graph.fundamental_cycles()
    if any(all(not row[e] for row in cycles) for e in range(graph.ne)):
        raise ValueError("a positive-loop admissible tight graph cannot have a bridge")
    chains = series_chains(graph)
    if len(chains) > 20:
        raise ValueError("admissibility verification exceeds the explicit twenty-chain limit")
    endpoints = [(chain.vertices[0], chain.vertices[-1]) for chain in chains]
    core = SimpleNamespace(nv=graph.nv, edges=endpoints)
    for mask in range(1, 1 << len(chains)):
        selected = tuple(i for i in range(len(chains)) if mask & (1 << i))
        cycles = _cycle_count(core, selected)
        if cycles and sum(chains[i].length for i in selected) < dimension*cycles+1:
            raise ValueError("the graph violates a local Laman inequality on original chains %s" % (selected,))
    return True


def evaluate_bulk(graph, sources=None, *, H=None, T=0, word=(), route="auto",
                  period_timeout=180, polynomial_options=None, series_options=None):
    """Evaluate the full bulk kernel or an original holomorphic-current word.

    A finite answer is a BulkKernelValue.  Missing finite period data returns
    an ExactBulkSeries with an explicit reason; ``terms`` and ``interval``
    compute its rational coefficients and rigorous remainders.  ``route``
    may force ``series`` for an independent integral control.  All sources
    must be exact rationals; full and non-root source conventions are those
    of htk.brute.sources.  No dimension, source, or current is silently cut.
    """
    T = _integer(T, "topological dimension")
    if route not in ("auto", "periods", "series"):
        raise ValueError("route must be auto, periods, or series")
    lam = source_vectors(graph, sources, H=H)
    inferred = len(lam[0]) if lam else 0
    H = inferred if H is None else _integer(H, "holomorphic dimension")
    if H != inferred or any(len(row) != H for row in lam):
        raise ValueError("source dimensions disagree with H")
    if not graph.is_connected() or graph.ne != (H+T)*graph.loops+1:
        raise ValueError("a connected tight bulk graph is required")
    word = tuple((_integer(e, "original edge"), _integer(a, "spinor index")) for e, a in word)
    if any(e >= graph.ne or a >= H for e, a in word):
        raise ValueError("current slot outside the original graph/source vectors")
    if graph.loops == 0:
        current = brute.currents(graph, [1], lam)
        return _value({"1": prod(current[e][a] for e, a in word)}, "single edge")
    _check_admissible(graph, H+T)
    if T >= 2 or (T == 1 and graph.loops % 2):
        return _value({}, "exact bulk vanishing", reason="rho squared" if T >= 2 else "odd Gaussian moment",
                      original_current_slots=word, zero_topological_momentum=True)
    if route == "auto":
        saved = _saved_source_polynomial(graph, lam, word, H, T)
        if saved is not None:
            return saved
        try:
            reduction = sp_reduce(graph, D=H+T)
        except NotSeriesParallel:
            reduction = None
        if reduction is not None:
            result = reduction.integrate(lam, T=T, powers=dict(Counter(word)), check=2)
            if isinstance(result, PiNumber):
                coefficients = {"1" if not k else "pi^%d" % k: v for k, v in result.c.items()}
            else:
                coefficients = {"1": rational(result)}
            return _value(coefficients, "series-parallel", independent_density_controls=2,
                          all_original_current_slots_retained=True)
    density = BulkPolynomial(graph, lam, H=H, T=T, **(polynomial_options or {}))
    if not density.numerator:
        return _value({}, "exact zero polynomial", complete_polynomial_circuit_identities=True)
    density.check_at(tuple(Fraction(e+2, e+1) for e in range(graph.ne)))
    edges, coefficients, power, chains = _scalar_components(density, word)
    if not coefficients:
        return _value({}, "exact series-fibre cancellation", complete_original_current_word=word)
    if any(a == b for a, b in edges):
        # A pure polygon contracts to one scalar loop; its projective
        # zero-dimensional integral is one.  Other self loops are not
        # silently assigned a period.
        if len(edges) == 1 and len(coefficients) == 1:
            return _value({"1": next(iter(coefficients.values()))}, "complete series fibres")
        raise ValueError("the contracted scalar graph contains a nontrivial self loop")
    reason = "complete rational series explicitly requested"
    if route != "series":
        worker = _BulkPeriodWorker(_integer(period_timeout, "period_timeout", 1))
        try:
            value, stats = worker.integrate(edges, str(power), coefficients)
            return _value(value, "exact scalar periods", **stats,
                          full_graph_density_control=True, all_source_fibres_integrated=True)
        except (MissingForestPeriod, ForestPeriodBudgetExceeded) as error:
            reason = str(error)
        finally:
            worker.close()
    core = Graph(1+max(v for edge in edges for v in edge), edges)
    from flint import fmpq_mpoly_ctx
    context = fmpq_mpoly_ctx.get(tuple("t%d" % e for e in range(core.ne)))
    from .forest_cube_expansion import _q
    numerator = context.from_dict({tuple(a-1 for a in alpha): _q(c) for alpha, c in coefficients.items()})
    scalar = SimpleNamespace(graph=core, symbolic_sources=False,
                             inserted=lambda _: (numerator, power))
    result = ExactBulkSeries(scalar, reason=reason, **(series_options or {}))
    result.original_graph = graph
    result.original_sources = tuple(map(tuple, lam))
    result.original_word = word
    result.series_chains = chains
    return result
