"""Exact integration of current polynomials on fully sourced forest cores.

The fast path reuses the generic series-parallel density.  Otherwise the
linear current determinant from ``explore/nonSP/verify/k4_periods/current_minor.py``
constructs the complete Schwinger numerator, and exact tree flows insert
every requested current.  Series beta factors and graph symmetries are
combined before a read-only period lookup.  Missing finite period data
produces an explicit convergent-series object, never a zero or a truncated
number presented as the answer.  Its rational terms and remainder bounds
reuse ``explore/nonSP/verify/hepp_series/period.py``.

The two small research modules have isolated import names.  The larger
period stack runs in a private worker, so its historical top-level module
names cannot replace modules in the calling process.  No table is built.
"""
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import importlib.util
from itertools import combinations
import json
from math import comb, factorial, lcm, prod
from numbers import Integral
from pathlib import Path
import selectors
import subprocess
import sys
from types import MappingProxyType, SimpleNamespace
import weakref

from . import brute
from .forest_contractions import normalize_forest_sources
from .graph import Graph
from .series import series_chains
from .spred import NotSeriesParallel, plan_generic, sp_reduce
from .treeflow import tree_flow

__all__ = ["CoreCurrentPeriods", "ExactPeriod", "ExactForestSeries", "SeriesTerm",
           "ForestPeriodBudgetExceeded", "MissingForestPeriod"]


class ForestPeriodBudgetExceeded(RuntimeError):
    """An explicit exact-computation budget was exceeded."""


class MissingForestPeriod(NotImplementedError):
    """The finite period library has no value; an exact series is available."""


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError("%s must be an integer at least %d" % (name, minimum))
    return int(value)


def _fraction(value):
    if hasattr(value, "p") and hasattr(value, "q"):
        return Fraction(int(value.p), int(value.q))
    return Fraction(value)


def _rational(value):
    if isinstance(value, float):
        raise TypeError("exact rational coefficients are required")
    return Fraction(value)


@lru_cache(maxsize=None)
def _research_module(relative, name):
    path = Path(__file__).resolve().parents[2] / relative
    if not path.is_file():
        raise FileNotFoundError("the exact forest backend needs this repository module: %s" % path)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ExactPeriod(Mapping):
    """A finite exact rational combination of the explicitly named zeta products.

    Labels are ``'1'``, ``'z3'``, ``'z5'`` or products such as ``'z3*z3'``.
    Equality compares this expression; it makes no independence assertion
    about distinct products of zeta values.
    """
    def __init__(self, coefficients):
        clean = {}
        for label, value in coefficients.items():
            if label != "1" and (not isinstance(label, str) or any(
                    not part.startswith("z") or not part[1:].isdigit()
                    or int(part[1:]) < 2 for part in label.split("*"))):
                raise ValueError("unknown exact period label: %r" % label)
            value = _rational(value)
            if value:
                clean[label] = value
        self.coefficients = MappingProxyType(clean)

    def __getitem__(self, key):
        return self.coefficients[key]

    def __iter__(self):
        return iter(self.coefficients)

    def __len__(self):
        return len(self.coefficients)

    def __repr__(self):
        return "ExactPeriod(%r)" % dict(self.coefficients)

    def __eq__(self, other):
        if isinstance(other, Mapping):
            return dict(self.coefficients) == {k: Fraction(v) for k, v in other.items() if v}
        if isinstance(other, (Integral, Fraction)):
            return all(label == "1" for label in self) and self.coefficients.get("1", 0) == other
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, (Integral, Fraction)):
            other = {"1": Fraction(other)}
        if not isinstance(other, Mapping):
            return NotImplemented
        value = dict(self.coefficients)
        for label, coefficient in other.items():
            value[label] = value.get(label, Fraction(0)) + Fraction(coefficient)
        return _finite(value)

    __radd__ = __add__

    def __mul__(self, scalar):
        if not isinstance(scalar, (Integral, Fraction)):
            return NotImplemented
        return _finite({label: scalar * value for label, value in self.items()})

    __rmul__ = __mul__

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + (-other)

    def as_sympy(self):
        import sympy
        return sum((sympy.Rational(value.numerator, value.denominator) *
                    (1 if label == "1" else prod(sympy.zeta(int(part[1:]))
                                                for part in label.split("*")))
                    for label, value in self.items()), sympy.Integer(0))


def _finite(coefficients):
    coefficients = {key: Fraction(value) for key, value in coefficients.items() if value}
    if all(key == "1" for key in coefficients):
        return coefficients.get("1", Fraction(0))
    return ExactPeriod(coefficients)


@dataclass(frozen=True)
class SeriesTerm:
    order: int
    term: Fraction
    partial_sum: Fraction
    error_bound: Fraction

    @property
    def lower(self):
        return self.partial_sum - self.error_bound

    @property
    def upper(self):
        return self.partial_sum + self.error_bound


class ExactForestSeries:
    """An exact infinite rational series with a certified remainder after every order.

    ``components`` is a finite sequence of ``(q, numerator)`` on ``edges``.
    The coefficient dictionary contains ordinary monomial exponents.  The
    integral has already been replaced by the explicit convergent series
    implemented in ``terms``.  A returned partial sum alone is never called
    the exact value.  Construction is lazy; the edge budget is enforced
    before the subset recursion is allocated.
    """
    def __init__(self, edges, components, *, max_edges=18, max_terms=1000000, reason=None,
                 constant=0):
        self.edges = tuple(tuple(_integer(vertex, "edge vertex") for vertex in edge) for edge in edges)
        if any(len(edge) != 2 for edge in self.edges):
            raise ValueError("every graph edge has two endpoints")
        self.components = tuple((_integer(q, "denominator power", 1),
                                 MappingProxyType({tuple(_integer(value, "Schwinger exponent") for value in a):
                                                   _rational(c) for a, c in numerator.items() if c}))
                                for q, numerator in components)
        if any(len(exponent) != len(self.edges) for _, numerator in self.components for exponent in numerator):
            raise ValueError("one Schwinger exponent is required per core edge")
        self.max_edges = _integer(max_edges, "max_series_edges", 1)
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.reason = reason
        self.constant = _rational(constant)
        self._engines = None

    def __repr__(self):
        return "ExactForestSeries(edges=%d, components=%d, numerator_terms=%d)" % (
            len(self.edges), len(self.components), sum(len(p) for _, p in self.components))

    def _setup(self):
        if len(self.edges) > self.max_edges:
            raise ForestPeriodBudgetExceeded("exact series subset recursion exceeds max_series_edges=%d; "
                                             "raise this explicit budget to allocate 2**edges subsets" % self.max_edges)
        if self._engines is None:
            module = _research_module("explore/nonSP/verify/hepp_series/period.py",
                                      "_htk_forest_hepp_series")
            self._engines = tuple(module.HeppSeries(self.edges, q, dict(numerator))
                                  for q, numerator in self.components)
        return self._engines

    def terms(self, last, *, method="gegenbauer"):
        """Yield exact rational coefficients, partial sums and rigorous absolute tails."""
        last = _integer(last, "series order")
        engines = self._setup()
        # PU^j has at most len(P)*comb(trees+j-1,j) expanded terms. This is
        # a conservative pre-allocation bound, not a numerical truncation.
        for engine in engines:
            predicted = len(engine.numerator) * comb(engine.tree_count + last - 1, last)
            if predicted > self.max_terms:
                raise ForestPeriodBudgetExceeded("series expansion exceeds max_terms=%d; "
                                                 "increase the explicit exact polynomial budget" % self.max_terms)
        streams = [engine.terms(last, method=method, permutations="graph") for engine in engines]
        for order in range(last + 1):
            rows = [next(stream) for stream in streams]
            yield SeriesTerm(order, sum((_fraction(row[1]) for row in rows),
                                       self.constant if order == 0 else Fraction(0)),
                             sum((_fraction(row[2]) for row in rows), self.constant),
                             sum((_fraction(row[3]) for row in rows), Fraction(0)))

    def interval(self, last, *, method="gegenbauer"):
        row = None
        for row in self.terms(last, method=method):
            pass
        return row.lower, row.upper


def _stop_process(process):
    if process.poll() is None:
        process.terminate()
        try:
            process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
    for stream in (process.stdin, process.stdout):
        if stream is not None:
            stream.close()


class _PeriodWorker:
    def __init__(self, timeout):
        self.timeout = timeout
        path = Path(__file__).with_name("_forest_period_worker.py")
        self.process = subprocess.Popen([sys.executable, "-u", str(path)], stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                        text=True, bufsize=1)
        self._finalizer = weakref.finalize(self, _stop_process, self.process)

    def close(self):
        self._finalizer()

    def integrate(self, edges, q, coefficients):
        if self.process.poll() is not None:
            raise ArithmeticError("isolated finite-period worker is no longer running")
        request = {"edges": edges, "q": q,
                   "coefficients": [[list(a), str(c)] for a, c in coefficients.items()]}
        self.process.stdin.write(json.dumps(request, separators=(",", ":")) + "\n")
        self.process.stdin.flush()
        with selectors.DefaultSelector() as ready:
            ready.register(self.process.stdout, selectors.EVENT_READ)
            if not ready.select(self.timeout):
                self.close()
                raise ForestPeriodBudgetExceeded("finite period lookup exceeded period_timeout=%s seconds" % self.timeout)
        line = self.process.stdout.readline()
        if not line:
            raise ArithmeticError("isolated finite-period worker exited without a result")
        answer = json.loads(line)
        if "missing" in answer:
            raise MissingForestPeriod(answer["missing"])
        if "error" in answer:
            raise ArithmeticError("finite-period worker: " + answer["error"])
        return {label: Fraction(value) for label, value in answer["coefficients"].items()}, answer["stats"]


class CoreCurrentPeriods:
    """Reusable source-complete integration after local forest contractions.

    Input polynomials map edge-major exponent tuples of length ``2*G.ne``
    to rational coefficients.  Component zero precedes component one for
    every edge.  All original orientation and contraction signs must already
    be present in the polynomial.  Rational source rows use htk conventions.
    Budgets limit allocations, not the mathematical family or loop order.
    """
    def __init__(self, G, sources, marked=None, *, check=2, max_terms=1000000,
                 max_tree_candidates=1000000, max_period_terms=250000,
                 max_series_edges=18, period_timeout=120, max_matrix_entries=1000000):
        self.G = self.graph = Graph(G.nv, G.edges, G.name)
        if not self.G.is_connected() or self.G.ne != 2 * self.G.loops + 1:
            raise ValueError("a connected tight H=2 graph is required")
        self.sources = normalize_forest_sources(self.G, sources)
        self.marked = None if marked is None else tuple(marked)
        self.check = _integer(check, "check count")
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.max_tree_candidates = _integer(max_tree_candidates, "max_tree_candidates", 1)
        self.max_period_terms = _integer(max_period_terms, "max_period_terms", 1)
        self.max_series_edges = _integer(max_series_edges, "max_series_edges", 1)
        self.max_matrix_entries = _integer(max_matrix_entries, "max_matrix_entries", 1)
        self.period_timeout = _integer(period_timeout, "period_timeout", 1)
        self._sp = self._density = self._bare = self._currents = self._U = self._worker = None
        self._context = self._core = self._chains = self._permutations = None
        self._current_powers = {}
        self._last_cancellations = {}
        self._sp_density_factors = None
        self.last_stats = {}
        try:
            self._sp = sp_reduce(self.G, D=2, plan=plan_generic(self.G))
        except NotSeriesParallel:
            pass

    def close(self):
        if self._worker is not None:
            self._worker.close()
            self._worker = None

    def _polynomial(self, polynomial):
        if not isinstance(polynomial, Mapping):
            raise TypeError("current polynomial must map full exponent tuples to rational coefficients")
        out = {}
        for exponent, coefficient in polynomial.items():
            exponent = tuple(_integer(value, "current exponent") for value in exponent)
            if len(exponent) != 2 * self.G.ne:
                raise ValueError("a current exponent must have length 2*number_of_edges")
            coefficient = _rational(coefficient)
            if coefficient:
                out[exponent] = out.get(exponent, Fraction(0)) + coefficient
        if len(out) > self.max_terms:
            raise ForestPeriodBudgetExceeded("input current polynomial exceeds max_terms")
        return {exponent: coefficient for exponent, coefficient in out.items() if coefficient}

    def _bounded(self, polynomial):
        if len(polynomial) > self.max_terms:
            raise ForestPeriodBudgetExceeded("exact polynomial exceeds max_terms=%d" % self.max_terms)
        return polynomial

    def _determinant(self, matrix, ctx):
        """Budgeted version of current_minor's block/Bareiss determinant."""
        import networkx as nx
        size = len(matrix)
        if size == 0:
            return ctx.constant(1)
        if size == 1:
            return self._bounded(matrix[0][0])
        support = nx.Graph()
        support.add_nodes_from(range(2 * size))
        support.add_edges_from((i, size + j) for i, row in enumerate(matrix)
                               for j, value in enumerate(row) if value)
        matching = nx.bipartite.maximum_matching(support, range(size))
        if len(matching) != 2 * size:
            return ctx.constant(0)
        permutation = [matching[i] - size for i in range(size)]
        directed = nx.DiGraph()
        directed.add_nodes_from(range(size))
        directed.add_edges_from((i, matching[j]) for i in range(size) for j in support[i])
        blocks = list(nx.strongly_connected_components(directed))
        if len(blocks) > 1:
            sign = (-1) ** sum(x > y for i, x in enumerate(permutation) for y in permutation[i + 1:])
            result = sign * ctx.constant(1)
            for block in blocks:
                indices = sorted(block)
                result = self._bounded(result * self._determinant(
                    [[matrix[i][permutation[j]] for j in indices] for i in indices], ctx))
                if not result:
                    break
            return result
        rows = [list(row) for row in matrix]
        previous, sign = ctx.constant(1), 1
        while len(rows) > 1:
            row_counts = [sum(bool(value) for value in row) for row in rows]
            column_counts = [sum(bool(value) for value in column) for column in zip(*rows)]
            chosen = min(((len(value), row_counts[i] * column_counts[j], i, j)
                          for i, row in enumerate(rows) for j, value in enumerate(row) if value), default=None)
            if chosen is None:
                return ctx.constant(0)
            _, _, i, j = chosen
            if i:
                rows[0], rows[i] = rows[i], rows[0]
                sign = -sign
            if j:
                for row in rows:
                    row[0], row[j] = row[j], row[0]
                sign = -sign
            pivot = rows[0][0]
            following = []
            for row in rows[1:]:
                line = []
                for x, y in zip(row[1:], rows[0][1:]):
                    numerator = self._bounded(pivot * x - row[0] * y)
                    quotient, remainder = divmod(numerator, previous)
                    if remainder:
                        raise ArithmeticError("non-exact Bareiss division in current numerator")
                    line.append(self._bounded(quotient))
                following.append(line)
            previous, rows = pivot, following
        return self._bounded(sign * rows[0][0])

    def _source_numerator(self):
        if self._bare is not None:
            return
        from flint import fmpq, fmpq_mpoly_ctx
        n = self.G.loops
        if (n * n - 1) ** 2 > self.max_matrix_entries:
            raise ForestPeriodBudgetExceeded("linear determinant exceeds max_matrix_entries")
        module = _research_module("explore/nonSP/verify/k4_periods/current_minor.py",
                                  "_htk_forest_current_minor")
        scales = [lcm(*(row[i].denominator for row in self.sources)) for i in range(2)]
        scaled = [[int(row[i] * scales[i]) for i in range(2)] for row in self.sources]
        flow = [[int(value) for value in row] for row in tree_flow(self.G, scaled)]
        proxy = SimpleNamespace(n=n, N=self.G.ne, C=self.G.cycle_matrix())
        polynomial = module.numerator(proxy, scaled, flow, determinant=self._determinant)
        self._bounded(polynomial)
        denominator = prod(value ** n for value in scales)
        self._context = fmpq_mpoly_ctx.get(tuple("t%d" % edge for edge in range(self.G.ne)))
        self._bare = self._context.from_dict({tuple(map(int, exponent)): fmpq(int(coefficient), denominator)
                                             for exponent, coefficient in polynomial.terms()})
        for index in range(self.check):
            point = [Fraction(edge + index + 2, edge + 1) for edge in range(self.G.ne)]
            U = brute.symanzik_U(self.G, point)
            expected = brute.density(self.G, point, self.sources).rational() * U ** (n + 2)
            actual = _fraction(self._bare(*(fmpq(value.numerator, value.denominator) for value in point)))
            if actual != expected:
                raise ArithmeticError("linear numerator failed independent Kirchhoff density check")

    def _tree_polynomials(self):
        if self._currents is not None:
            return
        from flint import fmpq
        self._source_numerator()
        candidates = comb(self.G.ne, self.G.nv - 1)
        if candidates > self.max_tree_candidates:
            raise ForestPeriodBudgetExceeded("tree current construction exceeds max_tree_candidates")
        scalar, currents = {}, [[{} for _ in range(2)] for _ in self.G.edges]
        for tree in combinations(range(self.G.ne), self.G.nv - 1):
            parent = list(range(self.G.nv))
            def find(vertex):
                while parent[vertex] != vertex:
                    parent[vertex] = parent[parent[vertex]]
                    vertex = parent[vertex]
                return vertex
            valid = True
            for edge in tree:
                left, right = (find(v) for v in self.G.edges[edge])
                if left == right:
                    valid = False
                    break
                parent[left] = right
            if not valid:
                continue
            selected = set(tree)
            exponent = tuple(int(edge not in selected) for edge in range(self.G.ne))
            scalar[exponent] = 1
            flow = tree_flow(self.G, self.sources, tree)
            for edge in tree:
                for component in range(2):
                    if flow[edge][component]:
                        currents[edge][component][exponent] = flow[edge][component]
        def convert(values):
            return self._context.from_dict({exponent: fmpq(value.numerator, value.denominator)
                                             for exponent, value in values.items()})
        self._U = convert({exponent: Fraction(value) for exponent, value in scalar.items()})
        self._currents = tuple(convert(values) for row in currents for values in row)
        for index in range(self.check):
            point = [Fraction(edge + 3, edge + index + 2) for edge in range(self.G.ne)]
            U = brute.symanzik_U(self.G, point)
            actual = brute.currents(self.G, point, self.sources)
            for edge in range(self.G.ne):
                for component in range(2):
                    observed = _fraction(self._currents[2 * edge + component](
                        *(fmpq(value.numerator, value.denominator) for value in point)))
                    if observed != U * actual[edge][component]:
                        raise ArithmeticError("tree numerator failed independent original-current check")

    def _chain_data(self):
        if self._chains is not None:
            return
        self._chains = series_chains(self.G)
        vertices = sorted({vertex for chain in self._chains for vertex in (chain.vertices[0], chain.vertices[-1])})
        self._core = tuple((vertices.index(chain.vertices[0]), vertices.index(chain.vertices[-1]))
                           for chain in self._chains)
        import networkx as nx
        incidence = nx.Graph()
        for edge, ends in enumerate(self._core):
            incidence.add_node(("edge", edge), kind="edge")
            for vertex in ends:
                incidence.add_node(("vertex", vertex), kind="vertex")
                incidence.add_edge(("edge", edge), ("vertex", vertex))
        same = nx.isomorphism.categorical_node_match("kind", None)
        matcher = nx.isomorphism.GraphMatcher(incidence, incidence, node_match=same)
        permutations = []
        for mapping in matcher.isomorphisms_iter():
            permutations.append(tuple(mapping["edge", edge][1] for edge in range(len(self._core))))
            # Any set of actual automorphisms gives valid equality-preserving
            # aggregation; exhaustiveness affects speed only.
            if len(permutations) == 4096:
                break
        self._permutations = tuple(permutations)

    def _schwinger(self, polynomial):
        self._source_numerator()
        if not self._bare or not polynomial:
            return {}
        from flint import fmpq, fmpq_mpoly_ctx
        grouped = {}
        for exponent, coefficient in polynomial.items():
            degree = sum(exponent)
            grouped.setdefault(degree, {})[exponent] = coefficient
        # Kirchhoff identities and local current brackets create exact U
        # factors only after all monomials have been added. Cancel them here,
        # before a large bare multiplication or any weighted-period lookup.
        reduced = {}
        self._last_cancellations = {}
        current_context = fmpq_mpoly_ctx.get(tuple("y%d" % index for index in range(2 * self.G.ne)))
        for degree, terms in grouped.items():
            q = self.G.loops + 2 + degree
            cancelled = 0
            if degree:
                self._tree_polynomials()
            expression = current_context.from_dict({exponent: fmpq(coefficient.numerator, coefficient.denominator)
                                                    for exponent, coefficient in terms.items()})
            content, factors = expression.factor()
            value = self._context.constant(content)
            for factor, multiplicity in factors:
                mapped = self._context.constant(0)
                for exponent, coefficient in factor.terms():
                    term = self._context.constant(coefficient)
                    for location, power in enumerate(exponent):
                        power = int(power)
                        if not power:
                            continue
                        key = location, power
                        if key not in self._current_powers:
                            self._current_powers[key] = self._bounded(self._currents[location] ** power)
                        term = self._bounded(term * self._current_powers[key])
                    mapped = self._bounded(mapped + term)
                if not mapped:
                    value = mapped
                    break
                factor_cancelled = 0
                if self._U:
                    while True:
                        quotient, remainder = divmod(mapped, self._U)
                        if remainder:
                            break
                        mapped = quotient
                        factor_cancelled += 1
                multiplicity = int(multiplicity)
                cancelled += factor_cancelled * multiplicity
                q -= factor_cancelled * multiplicity
                value = self._bounded(value * self._bounded(mapped ** multiplicity))
            if not value:
                continue
            if degree and self._U:
                while True:
                    quotient, remainder = divmod(value, self._U)
                    if remainder:
                        break
                    value = quotient
                    q -= 1
                    cancelled += 1
            self._last_cancellations[degree] = cancelled
            reduced[q] = self._bounded(reduced.get(q, self._context.constant(0)) + value)
        result = {}
        for q, value in reduced.items():
            numerator = self._bounded(self._bare * value)
            if not numerator:
                continue
            if self._U:
                while True:
                    quotient, remainder = divmod(numerator, self._U)
                    if remainder:
                        break
                    numerator = quotient
                    q -= 1
            result[q] = self._bounded(result.get(q, self._context.constant(0)) + numerator)
        return {q: value for q, value in result.items() if value}

    def _core_polynomials(self, polynomial):
        schwinger = self._schwinger(polynomial)
        self._chain_data()
        result = {}
        for q, numerator in schwinger.items():
            values = {}
            for exponent, coefficient in numerator.terms():
                exponent = tuple(map(int, exponent))
                alpha = tuple(sum(exponent[edge] + 1 for edge in chain.edges) for chain in self._chains)
                factor = _fraction(coefficient) * Fraction(prod(factorial(k) for k in exponent),
                                                          prod(factorial(k - 1) for k in alpha))
                orbit = min(tuple(alpha[edge] for edge in permutation) for permutation in self._permutations)
                values[orbit] = values.get(orbit, Fraction(0)) + factor
            values = {alpha: coefficient for alpha, coefficient in values.items() if coefficient}
            if len(values) > self.max_period_terms:
                raise ForestPeriodBudgetExceeded("weighted core periods exceed max_period_terms")
            result[q] = values
        self.last_stats = {"backend": "schwinger_numerator", "input_current_terms": len(polynomial),
                           "bare_numerator_terms": len(self._bare),
                           "inserted_numerator_terms": sum(len(value) for value in schwinger.values()),
                           "core_edges": len(self._core), "used_automorphisms": len(self._permutations),
                           "core_period_orbits": sum(len(value) for value in result.values()),
                           "denominator_powers": sorted(result), "independent_checks": self.check}
        self.last_stats["current_U_cancellations"] = dict(self._last_cancellations)
        return result

    def _integrate_sp(self, polynomial):
        """Factor in current variables before substituting the ratio currents."""
        from flint import fmpq, fmpq_mpoly_ctx
        if self._density is None:
            self._density = self._sp.density(self.sources, T=0, method="generic", check=self.check)
        currents = tuple(value for row in self._sp.currents(self.sources) for value in row)
        ring = currents[0].ring
        context = fmpq_mpoly_ctx.get(tuple("y%d" % index for index in range(2 * self.G.ne)))
        expression = context.from_dict({exponent: fmpq(coefficient.numerator, coefficient.denominator)
                                        for exponent, coefficient in polynomial.items()})
        content, factorization = expression.factor()
        content = _fraction(content)
        factors = [ring(content.numerator) / content.denominator]
        def ratio_factors(value):
            coefficient, split = value.factor_list()
            result = [ring(coefficient)]
            for factor, multiplicity in split:
                result.extend([factor] * int(multiplicity))
            return result
        if self._sp_density_factors is None:
            self._sp_density_factors = tuple(part for factor in self._density.factors
                                            for part in ratio_factors(factor))
        for factor, multiplicity in factorization:
            mapped = ring.zero
            for exponent, coefficient in factor.terms():
                coefficient = _fraction(coefficient)
                term = ring(coefficient.numerator) / coefficient.denominator
                for value, power in zip(currents, exponent):
                    power = int(power)
                    if power:
                        term = self._bounded(term * value ** power)
                mapped = self._bounded(mapped + term)
            # Preserve repeated factors for the existing support-aware
            # simplex integrator; it multiplies only the factors needed by
            # the next eliminated ratio block.
            factors.extend(ratio_factors(mapped) * int(multiplicity))
        if self._density.T == 0 and not any(self._density.exps) and self._sp.blocks:
            value = self._density._integrate_polynomial([*self._sp_density_factors, *factors])
            if any(any(exponent) for exponent in value):
                raise ArithmeticError("series-parallel integration retained a ratio variable")
            value = Fraction(int(value.LC.numerator), int(value.LC.denominator))
        else:
            insertion = ring.one
            for factor in factors:
                insertion = self._bounded(insertion * factor)
            value = Fraction(self._density.integrate(insertion))
        self.last_stats = {"backend": "series_parallel_generic", "input_current_terms": len(polynomial),
                           "independent_checks": self.check,
                           "current_factors": sum(int(power) for _, power in factorization),
                           "distinct_current_factors": len(factorization),
                           "ratio_factors": len(factors),
                           "density_factors": len(self._sp_density_factors),
                           "substituted_factor_terms": sum(len(factor) for factor in factors)}
        return value

    def _as_series(self, polynomials, reason=None):
        components = [(q, {tuple(power - 1 for power in alpha): coefficient
                            for alpha, coefficient in values.items()})
                      for q, values in polynomials.items() if values]
        self.last_stats["backend"] = "exact_convergent_series"
        if reason:
            self.last_stats["finite_period_unavailable"] = reason
        return ExactForestSeries(self._core, components, max_edges=self.max_series_edges,
                                 max_terms=self.max_terms, reason=reason)

    def series(self, polynomial):
        """Return the complete explicit rational series, not a partial sum."""
        polynomial = self._polynomial(polynomial)
        if not polynomial or self._sp is not None:
            value = self.integrate(polynomial)
            self.last_stats["backend"] = "exact_terminating_series"
            return ExactForestSeries((), (), constant=value,
                                     max_edges=self.max_series_edges, max_terms=self.max_terms)
        return self._as_series(self._core_polynomials(polynomial))

    def integrate(self, polynomial, *, allow_series=True):
        """Complete finite period when known, otherwise its exact convergent series."""
        polynomial = self._polynomial(polynomial)
        if not polynomial:
            self.last_stats = {"backend": "exact_zero", "input_current_terms": 0}
            return Fraction(0)
        if self._sp is not None:
            return self._integrate_sp(polynomial)
        polynomials = self._core_polynomials(polynomial)
        if not any(polynomials.values()):
            self.last_stats["backend"] = "exact_zero"
            return Fraction(0)
        if self._worker is None:
            self._worker = _PeriodWorker(self.period_timeout)
        total, records = {}, []
        try:
            for q, coefficients in polynomials.items():
                value, record = self._worker.integrate(self._core, q, coefficients)
                records.append(record)
                for label, coefficient in value.items():
                    total[label] = total.get(label, Fraction(0)) + coefficient
        except MissingForestPeriod as error:
            if not allow_series:
                raise
            return self._as_series(polynomials, str(error))
        except (ForestPeriodBudgetExceeded, ArithmeticError, BrokenPipeError):
            # A timeout or failed process must not poison the reusable
            # source/numerator cache with references to closed pipes.
            self.close()
            raise
        self.last_stats["backend"] = "finite_periods"
        self.last_stats["period_records"] = records
        return _finite(total)
