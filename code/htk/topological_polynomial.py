"""Exact bulk Pfaffian numerators, with every source and original current.

This is a direct exterior-algebra calculation.  It does not fit samples.
The topological Wick contractions and the holomorphic minors are summed
before exact polynomial division.  In particular, no individual divergent
scalar monomial is integrated before the complete graph orbit is collected.

The API accepts rational source vectors, or ``sources=None, H=...`` for
independent formal sources at every non-root vertex.  The root source is
minus their sum.  It returns the full numerator of the reduced htk density,
and the numerator of every original edge current.  Only bulk, zero
topological source, holomorphic-current insertions are in scope.
"""

from fractions import Fraction
from itertools import combinations
from math import comb, lcm
from types import SimpleNamespace

from flint import fmpq, fmpq_mpoly_ctx

from . import brute
from .forest_cubes import _integer, _rational
from .linalg import det, inverse
from .forest_periods import _research_module
from .treeflow import tree_flow

__all__ = ["BulkPolynomial", "BulkPolynomialBudgetExceeded"]


class BulkPolynomialBudgetExceeded(RuntimeError):
    """A stated complete-polynomial budget was exceeded; no answer is returned."""


def rational(value):
    if hasattr(value, "p"):
        return Fraction(int(value.p), int(value.q))
    return _rational(value)


def _q(value):
    value = rational(value)
    return fmpq(value.numerator, value.denominator)


def _sign(order):
    return (-1) ** sum(a > b for i, a in enumerate(order) for b in order[i + 1:])


def source_vectors(graph, sources, H=None):
    """Validate every supplied vertex and component before normalization.

    In particular the conserved root row is checked in full, and unknown
    dictionary vertices cannot be silently discarded by the older helper.
    """
    dimension = None if H is None else _integer(H, "holomorphic dimension")
    if sources is None:
        if dimension not in (None, 0):
            raise ValueError("exact rational source vectors are required when H is positive")
        return [[] for _ in range(graph.nv-1)]
    if isinstance(sources, dict):
        supplied = {}
        for vertex, row in sources.items():
            vertex = _integer(vertex, "original source vertex")
            if vertex >= graph.nv:
                raise ValueError("source vertex outside the original graph")
            supplied[vertex] = list(row)
        rows = list(supplied.values())
    else:
        supplied = [list(row) for row in sources]
        if len(supplied) not in (graph.nv-1, graph.nv):
            raise ValueError("need one source vector per original non-root vertex, or per vertex")
        rows = supplied
    if dimension is None:
        dimension = len(rows[0]) if rows else 0
    for row in rows:
        if len(row) != dimension:
            raise ValueError("every supplied source row, including the root, must have H components")
        for value in row:
            _rational(value)
    return [[_rational(value) for value in row]
            for row in brute.sources(graph, supplied, H=dimension)]


class BulkPolynomial:
    """Construct ``numerator / U**power`` and all ``current_numerators / U``.

    ``T=0`` uses the holomorphic determinant; ``T=1`` uses its complete
    exterior product with the Gaussian Pfaffian.  At zero topological
    momentum, odd-loop T=1 densities and positive-loop T>=2 densities
    vanish identically.  The latter vanishing applies to the stated bulk
    sector, not to defects or inserted topological currents.
    """

    def __init__(self, graph, sources=None, *, H=None, T=1,
                 max_terms=1_000_000, max_total_terms=3_000_000,
                 max_minors=100_000, check=True, method="auto"):
        self.graph = graph
        self.T = _integer(T, "topological dimension")
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.max_total_terms = _integer(max_total_terms, "max_total_terms", 1)
        self.max_minors = _integer(max_minors, "max_minors", 1)
        if not graph.is_connected() or graph.loops < 1:
            raise ValueError("a connected graph with at least one loop is required")
        self.symbolic_sources = sources is None
        if sources is None:
            if H is None:
                raise ValueError("H is required for formal sources")
            self.H = _integer(H, "holomorphic dimension")
            normalized = None
        else:
            normalized = source_vectors(graph, sources, H=H)
            self.H = len(normalized[0])
            if any(len(row) != self.H for row in normalized):
                raise ValueError("every source must have the same dimension")
            if H is not None and _integer(H, "holomorphic dimension") != self.H:
                raise ValueError("H disagrees with the source vectors")
            normalized = [[_rational(x) for x in row] for row in normalized]
        if graph.ne != (self.H + self.T) * graph.loops + 1:
            raise ValueError("the graph does not satisfy N-1=(H+T)*loops")
        self.source_order = tuple((v, a) for v in range(graph.nv - 1)
                                  for a in range(self.H)) if sources is None else ()
        names = tuple("t%d" % e for e in range(graph.ne)) + tuple(
            "lambda%d_%d" % (v, a) for v, a in self.source_order)
        self.context = fmpq_mpoly_ctx.get(names, ordering="lex")
        self.names = names
        variables = self.context.gens()
        self.times = variables[:graph.ne]
        if sources is None:
            self.sources = tuple(tuple(variables[graph.ne + v*self.H + a]
                                      for a in range(self.H))
                                 for v in range(graph.nv - 1))
        else:
            self.sources = tuple(tuple(self.context.constant(_q(x)) for x in row)
                                 for row in normalized)
        self.numeric_sources = normalized
        if method not in ("auto", "pfaffian", "gaussian_minor", "linear_minor"):
            raise ValueError("unknown exact density algorithm")
        if method == "auto":
            if normalized is None:
                method = "pfaffian"
            else:
                if self.T == 0:
                    method = "linear_minor"
                elif self.T == 1 and graph.loops in (2, 4) and self.H >= 2:
                    method = "gaussian_minor"
                else:
                    method = "pfaffian"
        if method in ("gaussian_minor", "linear_minor") and normalized is None:
            raise ValueError("the linear minor requires exact rational source values")
        if method == "gaussian_minor" and (self.T != 1 or graph.loops not in (2, 4)):
            raise ValueError("Gaussian polarization is implemented at T=1, loops=2 or 4")
        if method == "linear_minor" and self.T != 0:
            raise ValueError("the holomorphic linear minor requires T=0")
        self.method = method
        self.cycles = graph.cycle_matrix()
        self.minors_used = 0
        self.divisions = 0
        self._build_currents()
        self._build_density()
        if check:
            self.check_polynomial_identities()

    def _bounded(self, polynomial):
        if len(polynomial) > self.max_terms:
            raise BulkPolynomialBudgetExceeded("complete polynomial exceeds max_terms")
        return polynomial

    def _buffer(self, polynomials):
        if sum(len(p) for p in polynomials) > self.max_total_terms:
            raise BulkPolynomialBudgetExceeded("complete polynomial buffer exceeds max_total_terms")

    def _sum(self, terms):
        result = self.context.constant(0)
        for term in terms:
            result = self._bounded(result + term)
        return result

    def _product(self, terms):
        result = self.context.constant(1)
        for term in terms:
            result = self._bounded(result * term)
        return result

    def determinant(self, matrix):
        """Fraction-free determinant, checking every exact Bareiss division."""
        rows = [list(row) for row in matrix]
        if not rows:
            return self.context.constant(1)
        self._buffer(p for row in rows for p in row)
        if len(rows) > 6:
            import networkx as nx
            size = len(rows)
            support = nx.Graph()
            support.add_nodes_from(range(2*size))
            support.add_edges_from((i, size+j) for i, row in enumerate(rows)
                                   for j, p in enumerate(row) if p)
            matching = nx.bipartite.maximum_matching(support, range(size))
            if len(matching) != 2*size:
                return self.context.constant(0)
            permutation = [matching[i]-size for i in range(size)]
            directed = nx.DiGraph()
            directed.add_nodes_from(range(size))
            directed.add_edges_from((i, matching[j]) for i in range(size) for j in support[i])
            blocks = list(nx.strongly_connected_components(directed))
            if len(blocks) > 1:
                result = self.context.constant(_sign(permutation))
                for block in blocks:
                    indices = sorted(block)
                    result = self._bounded(result * self.determinant(
                        [[rows[i][permutation[j]] for j in indices] for i in indices]))
                return result
        previous = self.context.constant(1)
        orientation = 1
        while len(rows) > 1:
            row_counts = [sum(bool(p) for p in row) for row in rows]
            column_counts = [sum(bool(p) for p in column) for column in zip(*rows)]
            candidates = [(len(p), row_counts[i]*column_counts[j], i, j)
                          for i, row in enumerate(rows) for j, p in enumerate(row) if p]
            if not candidates:
                return self.context.constant(0)
            _, _, i, j = min(candidates)
            if i:
                rows[0], rows[i] = rows[i], rows[0]
                orientation *= -1
            if j:
                for row in rows:
                    row[0], row[j] = row[j], row[0]
                orientation *= -1
            pivot = rows[0][0]
            following = []
            for row in rows[1:]:
                entries = []
                for entry, top in zip(row[1:], rows[0][1:]):
                    numerator = self._bounded(pivot * entry - row[0] * top)
                    quotient, remainder = divmod(numerator, previous)
                    if remainder:
                        raise ArithmeticError("nonexact Bareiss division")
                    entries.append(self._bounded(quotient))
                following.append(entries)
            self._buffer(p for row in following for p in row)
            rows, previous = following, pivot
        return orientation * rows[0][0]

    def _build_currents(self):
        graph, cycles, times = self.graph, self.cycles, self.times
        n = graph.loops
        self.matrix = [[self._sum(cycles[r][e] * cycles[s][e] * times[e]
                                 for e in range(graph.ne))
                        for s in range(n)] for r in range(n)]
        self.U = self.determinant(self.matrix)
        self.adjugate = [[(-1) ** (r+s) * self.determinant(
            [[self.matrix[i][j] for j in range(n) if j != r]
             for i in range(n) if i != s]) for s in range(n)] for r in range(n)]
        tree = graph.spanning_tree()
        incidence = graph.incidence()
        # Solve only the unimodular tree incidence matrix.  This step keeps
        # all non-root sources, including those on subdivided edges.
        tree_inverse = inverse([[incidence[e][v] for e in tree]
                                for v in range(graph.nv - 1)])
        particular = [[self.context.constant(0) for _ in range(self.H)]
                      for _ in graph.edges]
        for k, e in enumerate(tree):
            for a in range(self.H):
                particular[e][a] = self._sum(_q(tree_inverse[k][v]) * self.sources[v][a]
                                             for v in range(graph.nv - 1))
        self.current_numerators = [[self.context.constant(0) for _ in range(self.H)]
                                   for _ in graph.edges]
        for a in range(self.H):
            forcing = [self._sum(cycles[r][e] * times[e] * particular[e][a]
                                 for e in range(graph.ne)) for r in range(n)]
            correction = [-self._sum(self.adjugate[r][s] * forcing[s] for s in range(n))
                          for r in range(n)]
            for e in range(graph.ne):
                self.current_numerators[e][a] = self._bounded(
                    self.U * particular[e][a]
                    + self._sum(cycles[r][e] * correction[r] for r in range(n)))
        self._buffer(p for row in self.current_numerators for p in row)

    def _wick(self, slots, response):
        if not slots:
            return self.context.constant(1)
        first, rest = slots[0], slots[1:]
        return self._sum(response[first][e] * self._wick(rest[:j] + rest[j+1:], response)
                         for j, e in enumerate(rest))

    def _holomorphic_minor(self, edges):
        if self.H == 1:
            integer_minor = int(det([[row[e] for e in edges] for row in self.cycles]))
            if not integer_minor:
                return self.context.constant(0)
            return integer_minor * self._product(self.current_numerators[e][0] for e in edges)
        return self.determinant([[self.cycles[r][e] * self.current_numerators[e][a]
                                  for e in edges]
                                 for a in range(self.H) for r in range(self.graph.loops)])

    def _build_density(self):
        graph, n = self.graph, self.graph.loops
        self.numerator = self.context.constant(0)
        self.power = Fraction(0)
        if self.T >= 2 or (self.T == 1 and n % 2):
            self.zero_reason = "rho squared" if self.T >= 2 else "odd Gaussian moment"
            return
        self.zero_reason = None
        if self.method in ("gaussian_minor", "linear_minor"):
            self._build_from_linear_minor()
            return
        remaining = tuple(range(graph.ne - 1))
        if self.T == 0:
            raw = self._holomorphic_minor(remaining)
            integer_power = self.H * (n + 1)
            scalar = (-1) ** (graph.ne - 1)
        else:
            if comb(graph.ne - 1, n) > self.max_minors:
                raise BulkPolynomialBudgetExceeded("complete exterior expansion exceeds max_minors")
            response = [[self._sum(self.cycles[r][e] * self.adjugate[r][s] * self.cycles[s][f]
                                   for r in range(n) for s in range(n))
                         for f in range(graph.ne)] for e in range(graph.ne)]
            self._buffer(p for row in response for p in row)
            raw = self.context.constant(0)
            for topological in combinations(remaining, n):
                top_minor = int(det([[row[e] for e in topological] for row in self.cycles]))
                if not top_minor:
                    continue
                holomorphic = tuple(e for e in remaining if e not in topological)
                hol = self._holomorphic_minor(holomorphic)
                if not hol:
                    continue
                self.minors_used += 1
                raw = self._bounded(raw + _sign(topological + holomorphic) * top_minor
                                    * self._wick(topological, response) * hol)
            integer_power = self.H * (n + 1) + n // 2
            scalar = Fraction((-1) ** (graph.ne - 1 + n // 2), 2 ** (n // 2))
        raw = self._bounded(_q(scalar) * raw)
        raw, remainder = divmod(raw, self.times[-1])
        if remainder:
            raise ArithmeticError("the complete form is not Euler-horizontal")
        while raw and integer_power:
            quotient, remainder = divmod(raw, self.U)
            if remainder:
                break
            raw = self._bounded(quotient)
            integer_power -= 1
            self.divisions += 1
        self.numerator = raw
        self.power = Fraction(integer_power) + Fraction(self.T, 2) if raw else Fraction(0)
        if raw:
            degree = graph.loops * self.power - graph.ne
            if any(sum(map(int, ex[:graph.ne])) != degree for ex in raw.to_dict()):
                raise ArithmeticError("density has incorrect projective homogeneity")

    def _linear_minor(self, extra=None):
        module = _research_module("explore/nonSP/verify/k4_periods/current_minor.py",
                                  "_htk_bulk_current_minor")
        scale = lcm(*(x.denominator for row in self.numeric_sources for x in row))
        sources = [[int(scale * x) for x in row] for row in self.numeric_sources]
        if extra is not None:
            for row, value in zip(sources, extra):
                row.append(int(value))
        sources.append([-sum(row[a] for row in sources)
                        for a in range(self.H + (extra is not None))])
        flow = [[int(x) for x in row] for row in tree_flow(self.graph, sources)]
        proxy = SimpleNamespace(n=self.graph.loops, N=self.graph.ne, C=self.cycles)
        numerator = module.numerator(proxy, sources, flow, edge_times=self.times,
                                     determinant=lambda matrix, context: self.determinant(matrix))
        return self._bounded(numerator / (scale ** (self.H * self.graph.loops)))

    def _build_from_linear_minor(self):
        """Gaussian averaging of one extra, complete source direction.

        If xi_e are independent with variance 1/(2t_e), their vertex
        sources B^t xi have covariance Lambda/2.  Solving Kirchhoff's
        equations gives precisely the zero-mean topological-current
        covariance.  Thus the T=1 form is sqrt(U) times the Gaussian
        average of the H+1 holomorphic form, with its n-row block moved
        to the front.  Quadratic/quartic polarization evaluates this
        Gaussian expectation exactly at n=2/4, without sampling.
        """
        n = self.graph.loops
        if self.T == 0:
            numerator = self._linear_minor()
            power = Fraction(n + self.H)
        else:
            incidence = self.graph.incidence()
            diagonals = [self._linear_minor(row) for row in incidence]
            self._buffer(diagonals)
            numerator = self.context.constant(0)
            def divide(polynomial, monomial):
                quotient, remainder = divmod(polynomial, monomial)
                if remainder:
                    raise ArithmeticError("Gaussian source polarization left a time denominator")
                return self._bounded(quotient)
            if n == 2:
                for e, value in enumerate(diagonals):
                    numerator = self._bounded(numerator + divide(value, 2*self.times[e]))
            else:
                for e, value in enumerate(diagonals):
                    numerator = self._bounded(numerator + divide(3*value, 4*self.times[e]**2))
                for e, f in combinations(range(self.graph.ne), 2):
                    plus = [a+b for a, b in zip(incidence[e], incidence[f])]
                    minus = [a-b for a, b in zip(incidence[e], incidence[f])]
                    mixed = self._bounded(self._linear_minor(plus) + self._linear_minor(minus)
                                          - 2*diagonals[e] - 2*diagonals[f])
                    numerator += divide(mixed, 8*self.times[e]*self.times[f])
                    numerator = self._bounded(numerator)
            numerator *= (-1) ** (self.H * n*n)
            power = Fraction(n + self.H) + Fraction(1, 2)
        while numerator and power >= 1:
            quotient, remainder = divmod(numerator, self.U)
            if remainder:
                break
            numerator = self._bounded(quotient)
            power -= 1
            self.divisions += 1
        self.numerator = numerator
        self.power = power if numerator else Fraction(0)
        if numerator and any(sum(map(int, ex[:self.graph.ne])) != n*power-self.graph.ne
                             for ex in numerator.to_dict()):
            raise ArithmeticError("Gaussian numerator has incorrect projective homogeneity")

    def check_polynomial_identities(self):
        """Check full adjugate, Kirchhoff and voltage identities symbolically."""
        graph, n, cycles = self.graph, self.graph.loops, self.cycles
        for r in range(n):
            for s in range(n):
                if self._sum(self.matrix[r][k] * self.adjugate[k][s] for k in range(n)) != (
                        self.U if r == s else 0):
                    raise ArithmeticError("incorrect adjugate")
        incidence = graph.incidence()
        for a in range(self.H):
            for v in range(graph.nv - 1):
                if self._sum(incidence[e][v] * self.current_numerators[e][a]
                             for e in range(graph.ne)) != self.U * self.sources[v][a]:
                    raise ArithmeticError("polynomial Kirchhoff identity failed")
            for row in cycles:
                if self._sum(row[e] * self.times[e] * self.current_numerators[e][a]
                             for e in range(graph.ne)):
                    raise ArithmeticError("polynomial voltage identity failed")
        return True

    def inserted(self, word=()):
        """Complete original-current word; each entry is (edge, spinor index)."""
        result = self.numerator
        power = self.power
        for edge, component in word:
            edge = _integer(edge, "original edge")
            component = _integer(component, "spinor index")
            if edge >= self.graph.ne or component >= self.H:
                raise ValueError("current slot outside the original graph/source vector")
            result = self._bounded(result * self.current_numerators[edge][component])
            power += 1
        return result, power

    def coefficients(self, word=()):
        """Separate time and source exponents without dropping any monomial."""
        numerator, power = self.inserted(word)
        result = {}
        for exponent, coefficient in numerator.to_dict().items():
            time = tuple(int(x) for x in exponent[:self.graph.ne])
            source = tuple(int(x) for x in exponent[self.graph.ne:])
            result.setdefault(time, {})[source] = rational(coefficient)
        return result, power

    def check_at(self, times, sources=None):
        """Independent grounded-Laplacian/Pfaffian control at exact inputs."""
        values = tuple(_rational(x) for x in times)
        if len(values) != self.graph.ne or min(values) <= 0:
            raise ValueError("need one positive rational time per original edge")
        lam = self.numeric_sources if sources is None else source_vectors(self.graph, sources, H=self.H)
        if lam is None:
            raise ValueError("a source specialization is required for this control")
        substitutions = {"t%d" % e: _q(t) for e, t in enumerate(values)}
        if self.symbolic_sources:
            substitutions.update({"lambda%d_%d" % (v, a): _q(lam[v][a])
                                  for v, a in self.source_order})
        def value(p):
            constant = p.subs(substitutions).to_dict()
            if any(any(e) for e in constant):
                raise ValueError("source specialization is incomplete")
            return sum((rational(c) for c in constant.values()), Fraction(0))
        actual_currents = brute.currents(self.graph, values, lam)
        u = value(self.U)
        for e in range(self.graph.ne):
            for a in range(self.H):
                if value(self.current_numerators[e][a]) != u * actual_currents[e][a]:
                    raise ArithmeticError("original-current control failed")
        if self.T <= 1:
            actual = brute.density(self.graph, values, lam, T=self.T, check_basic=True)
            if self.numerator:
                integer_power = self.power - Fraction(self.T, 2)
                if value(self.numerator) != actual.coeff * u ** int(integer_power):
                    raise ArithmeticError("original-graph density control failed")
            elif actual.coeff:
                raise ArithmeticError("reported zero density is nonzero")
        return True
