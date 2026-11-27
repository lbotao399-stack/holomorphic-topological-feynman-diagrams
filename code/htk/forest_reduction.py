"""Composed exact reduction of complete product-forest current polynomials.

All vertices, sources, edge orientations and current indices use the input
htk graph.  Kirchhoff polynomial identities are simplified before choosing
cubic messages, source-free stars, sourced ears, or a complete residual
period/series.  Inhomogeneous current polynomials are recombined before
period lookup.  The old scalar-word dispatch and its certificates remain
unchanged; import this additive interface from htk.forest_reduction.
"""
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
import gzip
import hashlib
from itertools import product
import json
from operator import index
from pathlib import Path
from types import MappingProxyType

from .forest_model import CubicForest, UnsupportedForestSources, recognize_product_forest
from .forest_moments import ForestMoments
from .forest_contractions import ForestContractions, normalize_forest_sources
from .forest_ears import ForestEarMoments, find_ear_bundles
from .forest_periods import (CoreCurrentPeriods, ExactForestSeries, ExactPeriod,
                             ForestPeriodBudgetExceeded, MissingForestPeriod,
                             _finite, _research_module)
from .graph import Graph
from .linalg import perm_sign
from .treeflow import tree_flow

__all__ = ['ForestReduction', 'ForestCurrentTensor', 'ForestReductionBudgetExceeded',
           'reduce_product_forest']

_ROOT = Path(__file__).resolve().parents[2]
_CERTIFICATE = 'explore/circuit-general/forest-finite-jets/source_certificate.json'
_CERTIFICATE_SHA256 = 'cca48d6331144f7604fda126240100d464e5034ff01ba16f133fb4d735232b55'
_CERTIFICATE_CACHE = None


class ForestReductionBudgetExceeded(RuntimeError):
    """A complete polynomial or tensor exceeds an explicit resource bound."""


def _integer(value, name, minimum=0):
    if isinstance(value, bool):
        raise TypeError('%s must be an integer' % name)
    try:
        value = index(value)
    except TypeError as error:
        raise TypeError('%s must be an integer' % name) from error
    if value < minimum:
        raise ValueError('%s must be at least %d' % (name, minimum))
    return value


def _rational(value):
    if isinstance(value, (float, complex, bool)):
        raise TypeError('current coefficients must be exact rational numbers')
    return Fraction(value)


def _bracket(left, right):
    return left[0] * right[1] - left[1] * right[0]


def _word(graph, powers):
    if powers is not None and not isinstance(powers, Mapping):
        raise TypeError('current powers must map (original edge, component) to nonnegative integers')
    exponent = [0] * (2 * graph.ne)
    for location, power in ({} if powers is None else powers).items():
        if not isinstance(location, tuple) or len(location) != 2:
            raise ValueError('a current location is (original edge, component)')
        edge, component = (_integer(location[0], 'current edge'),
                           _integer(location[1], 'current component'))
        if edge >= graph.ne or component >= 2:
            raise ValueError('current location lies outside the original graph')
        exponent[2 * edge + component] = _integer(power, 'current power')
    return {tuple(exponent): Fraction(1)}


def _verified_certificate():
    """Pin the proved content and check every recorded dependency, without solving.

    A cache is reused only while all dependency metadata are unchanged.
    Changed content or an absent file disables this optional finite shortcut;
    it cannot silently supply coefficients from a stale certificate.
    """
    global _CERTIFICATE_CACHE
    try:
        path = _ROOT / _CERTIFICATE
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != _CERTIFICATE_SHA256:
            return None, 'certificate_content_changed'
        data = json.loads(raw)
        expected = {}
        for name in ('implementation_sha256', 'period_table_sha256', 'reference_file_sha256'):
            expected.update(data[name])
        witness = path.parent / data['witness_file']
        if not witness.resolve().is_relative_to(_ROOT):
            return None, 'certificate_path_invalid'
        paths = {(_ROOT / name): digest for name, digest in expected.items()}
        paths[witness] = data['witness_file_sha256']
        metadata = []
        for filename in sorted(paths):
            if not filename.resolve().is_relative_to(_ROOT):
                return None, 'certificate_path_invalid'
            stat = filename.stat()
            metadata.append((str(filename), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns))
        metadata = tuple(metadata)
        if _CERTIFICATE_CACHE is not None and _CERTIFICATE_CACHE[0] == metadata:
            return _CERTIFICATE_CACHE[1], 'verified'
        for filename, expected_digest in paths.items():
            if hashlib.sha256(filename.read_bytes()).hexdigest() != expected_digest:
                return None, 'certificate_dependency_changed:' + str(filename.relative_to(_ROOT))
        if hashlib.sha256(gzip.decompress(witness.read_bytes())).hexdigest() != data['witness_uncompressed_sha256']:
            return None, 'certificate_witness_changed'
        if (not data['zero_residuals_for_every_coefficient'] or
                data['shared_elimination_statistics']['residual_periods'] or
                len(data['finite_coefficients_in_basis_order']) != len(data['source_reconstruction']['basis'])):
            return None, 'certificate_incomplete'
        _CERTIFICATE_CACHE = metadata, data
        return data, 'verified'
    except (OSError, KeyError, ValueError, TypeError):
        return None, 'certificate_unavailable_or_invalid'


def _certified_bare(graph, sources, marked):
    """Transport the certified source family, including the edge-order sign."""
    if graph.nv != 7 or graph.ne != 11:
        return None, 'different_topology'
    import networkx as nx
    data, status = _verified_certificate()
    if data is None:
        return None, status
    source_record = data['source_reconstruction']
    edges = tuple(map(tuple, source_record['original_graph']['edges']))
    reference_marked = tuple(source_record['marked'])
    first, second = nx.Graph(), nx.Graph()
    for vertex in range(graph.nv):
        role = 1 if vertex == marked[0] else 2 if vertex == marked[1] else 0
        first.add_node(vertex, role=role)
        role = 1 if vertex == reference_marked[0] else 2 if vertex == reference_marked[1] else 0
        second.add_node(vertex, role=role)
    first.add_edges_from(graph.edges)
    second.add_edges_from(edges)
    match = nx.isomorphism.GraphMatcher(first, second,
                 node_match=nx.isomorphism.categorical_node_match('role', 0))
    positions = {frozenset(edge): index for index, edge in enumerate(edges)}
    isomorphic = False
    for mapping in match.isomorphisms_iter():
        isomorphic = True
        mapped = [None] * graph.nv
        for vertex, row in enumerate(sources):
            mapped[mapping[vertex]] = row
        if any(any(mapped[vertex]) for vertex in source_record['structurally_zero_source_vertices']):
            continue
        permutation = [positions[frozenset((mapping[u], mapping[v]))] for u, v in graph.edges]
        sign = perm_sign(permutation)
        difference = tuple(a - b for a, b in zip(mapped[5], mapped[6]))
        factors = (_bracket(mapped[1], mapped[3]), _bracket(mapped[3], difference),
                   _bracket(mapped[1], difference))
        values = {}
        for (i, j), coefficients in zip(source_record['basis'], data['finite_coefficients_in_basis_order']):
            weight = sign * factors[0] ** (5 - i - j) * factors[1] ** i * factors[2] ** j
            for label, coefficient in coefficients.items():
                values[label] = values.get(label, Fraction()) + weight * Fraction(coefficient)
        return _finite(values), {'status': 'verified', 'certificate': _CERTIFICATE,
                                'certificate_sha256': _CERTIFICATE_SHA256,
                                'original_to_certificate_vertices': tuple(mapping[v] for v in range(graph.nv)),
                                'edge_permutation_sign': sign}
    return None, 'source_support_outside_certificate' if isomorphic else 'different_topology'


@dataclass(frozen=True)
class ForestCurrentTensor:
    """All ordered free-index entries; no contraction or endpoint factor.

    Each entry key is a complete tuple of zero-based spinor coordinates.
    The two coordinates are storage for each full two-dimensional spinor
    slot. Every one of 2**len(indices) entries is present, including zeros.
    Original edge orientations are retained. These are htk current tensors;
    physical endpoint derivative factors must be applied by the caller.
    """
    edge_slots: tuple
    indices: tuple
    entries: object

    @property
    def rank(self):
        return len(self.indices)

    def __getitem__(self, coordinates):
        return self.entries[tuple(coordinates)]


class ForestReduction:
    """Reuse one exact reduction for arbitrary finite original-current polynomials.

    Finite values are Fraction or ExactPeriod. An unavailable finite period
    gives the complete ExactForestSeries when allow_series=True. The finite
    source certificate applies only to a whole constant insertion after
    Kirchhoff simplification; no finite zeta term is incorrectly placed in
    the rational constant field of a series.
    """
    def __init__(self, G, sources, marked, *, check=2, max_terms=250000,
                 max_states=100000, max_source_samples=10000,
                 max_tree_candidates=1000000, max_period_terms=250000,
                 max_series_edges=18, period_timeout=120, max_matrix_entries=1000000,
                 max_tensor_entries=65536, allow_series=True, use_certificates=True,
                 max_minor_tree_subsets=1000, max_minor_plans=8):
        self.structure = recognize_product_forest(G, marked)
        self.graph = self.G = Graph(G.nv, G.edges, G.name)
        self.marked = self.structure.marked
        self.sources = normalize_forest_sources(self.graph, sources)
        self.check = _integer(check, 'check', 0)
        self.max_terms = _integer(max_terms, 'max_terms', 1)
        self.max_states = _integer(max_states, 'max_states', 1)
        self.max_source_samples = _integer(max_source_samples, 'max_source_samples', 1)
        self.max_tensor_entries = _integer(max_tensor_entries, 'max_tensor_entries', 1)
        self.max_minor_tree_subsets = _integer(max_minor_tree_subsets, 'max_minor_tree_subsets', 1)
        self.max_minor_plans = _integer(max_minor_plans, 'max_minor_plans', 1)
        if not isinstance(allow_series, bool) or not isinstance(use_certificates, bool):
            raise TypeError('allow_series and use_certificates must be booleans')
        self.allow_series, self.use_certificates = allow_series, use_certificates
        self.options = dict(check=self.check, max_terms=self.max_terms,
             max_tree_candidates=_integer(max_tree_candidates, 'max_tree_candidates', 1),
             max_period_terms=_integer(max_period_terms, 'max_period_terms', 1),
             max_series_edges=_integer(max_series_edges, 'max_series_edges', 1),
             period_timeout=_integer(period_timeout, 'period_timeout', 1),
             max_matrix_entries=_integer(max_matrix_entries, 'max_matrix_entries', 1))
        self._forest = self._contractions = self._core = self._ears = None
        self._linear_currents = None
        self._results, self._series_results, self._cube_engines = {}, {}, {}
        self.method, self.last_stats = 'unprepared', {}
        if self.structure.cubic:
            try:
                CubicForest(self.graph, self.marked).normalize_sources(self.sources)
            except UnsupportedForestSources:
                pass
            else:
                self._forest = ForestMoments(self.graph, self.sources, self.marked,
                                              max_terms=self.max_terms, max_states=self.max_states)
                self.method = 'cubic_boundary_moments'

    def close(self):
        if self._core is not None:
            self._core.close()
        if self._ears is not None:
            self._ears.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()

    def _polynomial(self, polynomial):
        if polynomial is None:
            return _word(self.graph, None)
        if not isinstance(polynomial, Mapping):
            raise TypeError('a current polynomial maps full exponent tuples to rational coefficients')
        result = {}
        for exponent, coefficient in polynomial.items():
            if not isinstance(exponent, tuple) or len(exponent) != 2 * self.graph.ne:
                raise ValueError('a current exponent tuple has length 2*original_edges')
            exponent = tuple(_integer(power, 'current exponent') for power in exponent)
            coefficient = _rational(coefficient)
            if coefficient:
                result[exponent] = result.get(exponent, Fraction()) + coefficient
        result = {exponent: coefficient for exponent, coefficient in result.items() if coefficient}
        if len(result) > self.max_terms:
            raise ForestReductionBudgetExceeded('input current polynomial exceeds max_terms')
        return result

    def normal_form(self, polynomial):
        """Exact current polynomial modulo all original Kirchhoff source equations.

        Tree currents are replaced by their source flow plus fundamental-cycle
        currents. Independent cotree currents retain their original edge slots.
        The returned polynomial uses the same length-2*ne exponent convention.
        """
        from flint import fmpq, fmpq_mpoly_ctx
        polynomial = self._polynomial(polynomial)
        if not polynomial or all(not any(exponent) for exponent in polynomial):
            return polynomial
        context = fmpq_mpoly_ctx.get(tuple('y%d' % i for i in range(2 * self.graph.ne)))
        if self._linear_currents is None:
            tree = tuple(self.graph.spanning_tree())
            cycles = self.graph.fundamental_cycles(tree)
            cotree = self.graph.cotree(tree)
            flow = tree_flow(self.graph, self.sources, tree)
            variables = context.gens()
            self._linear_currents = tuple(
                context.constant(fmpq(flow[edge][component].numerator, flow[edge][component].denominator))
                + sum((cycles[cycle][edge] * variables[2 * original + component]
                       for cycle, original in enumerate(cotree)), context.constant(0))
                for edge in range(self.graph.ne) for component in range(2))
        expression = context.from_dict({exponent: fmpq(value.numerator, value.denominator)
                                        for exponent, value in polynomial.items()})
        def multiply(left, right):
            # A conservative pre-allocation bound. Cancellations could make
            # a larger raw product small; deferring this optional normal
            # form retains the complete original polynomial in that case.
            if len(left) * len(right) > self.max_terms:
                raise ForestReductionBudgetExceeded('Kirchhoff product expansion exceeds max_terms')
            value = left * right
            if len(value) > self.max_terms:
                raise ForestReductionBudgetExceeded('Kirchhoff normal form exceeds max_terms')
            return value

        def power(value, exponent):
            result = context.constant(1)
            while exponent:
                if exponent & 1:
                    result = multiply(result, value)
                exponent >>= 1
                if exponent:
                    value = multiply(value, value)
            return result

        # Factoring a sparse polynomial of enormous degree need not be
        # cheaper than reduction. High degrees use bounded monomial
        # substitution directly; smaller bracket factors are kept intact.
        if len(polynomial) <= 64 and max(map(sum, polynomial), default=0) <= 32:
            scalar, factors = expression.factor()
        else:
            scalar, factors = fmpq(1), ((expression, 1),)
        powers_cache = {}
        mapped = []
        for factor, multiplicity in factors:
            value = context.constant(0)
            for exponent, coefficient in factor.terms():
                term = context.constant(coefficient)
                for location, degree in enumerate(exponent):
                    degree = int(degree)
                    if not degree:
                        continue
                    key = location, degree
                    if key not in powers_cache:
                        powers_cache[key] = power(self._linear_currents[location], degree)
                    term = multiply(term, powers_cache[key])
                value += term
                if len(value) > self.max_terms:
                    raise ForestReductionBudgetExceeded('Kirchhoff factor exceeds max_terms')
            if not value:
                return {}
            mapped.append((value, int(multiplicity)))
        result = context.constant(scalar)
        for value, multiplicity in mapped:
            result = multiply(result, power(value, multiplicity))
        return {tuple(map(int, exponent)): Fraction(str(value)) for exponent, value in result.terms() if value}

    def _prepare(self):
        if self._core is not None or self._ears is not None:
            return
        original = CoreCurrentPeriods(self.graph, self.sources, self.marked, **self.options)
        if original._sp is not None:
            self._core = original
            self.method = 'series_parallel_generic'
            return
        original.close()
        contractions = ForestContractions(self.graph, self.sources, self.marked,
                                           check=self.check, max_terms=self.max_terms)
        graph, sources, marked = (contractions.reduced_graph, contractions.reduced_sources,
                                  contractions.reduced_marked)
        core = CoreCurrentPeriods(graph, sources, marked, **self.options)
        ears, installed = None, False
        try:
            if core._sp is None and find_ear_bundles(graph, marked):
                ears = ForestEarMoments(graph, sources, marked,
                                       max_source_samples=self.max_source_samples, **self.options)
            self._contractions, self._core, self._ears = contractions, core, ears
            installed = True
        finally:
            if not installed:
                core.close()
        self.method = ('sourcefree_stars_then_sourced_ears' if contractions.steps and ears is not None
                       else 'sourced_ears' if ears is not None
                       else 'sourcefree_stars' if contractions.steps else 'core_periods')

    def _pushforward(self, polynomial):
        if self._contractions is None:
            return polynomial, {}
        groups = {}
        for exponent, value in polynomial.items():
            groups.setdefault(sum(exponent), {})[exponent] = value
        result, records = {}, []
        for degree in sorted(groups):
            pushed = self._contractions.contract(polynomial=groups[degree])
            records.append(pushed.stats)
            for exponent, value in pushed.polynomial.items():
                result[exponent] = result.get(exponent, Fraction()) + value
        result = {exponent: value for exponent, value in result.items() if value}
        if len(result) > self.max_terms:
            raise ForestReductionBudgetExceeded('recombined passive polynomial exceeds max_terms')
        return result, {'steps': self._contractions.steps, 'homogeneous_groups': tuple(records),
                        'complete_output_terms': len(result)}

    def _normalize_request(self, polynomial):
        original = self._polynomial(polynomial)
        # Cubic messages already avoid global loop-coordinate expansion.
        # A single nonconstant word has no cross-word cancellation to expose.
        if self._forest is not None or len(original) <= 1:
            return original, 'not_needed'
        try:
            return self.normal_form(original), 'exact'
        except ForestReductionBudgetExceeded:
            # This optional coordinate simplification may be more expensive
            # than local reduction. Retain the entire original polynomial.
            return original, 'deferred_at_budget'

    def _minor_plan(self):
        if self._core._sp is not None or self._core._bare is not None or self._core.G.loops < 4:
            return None
        module = _research_module('explore/circuit-general/symmetry/tree_minor.py',
                                   '_htk_uniform_forest_tree_minor')
        return module.prepare_source_numerator(self._core, max_tree_subsets=self.max_minor_tree_subsets,
                                               max_plans=self.max_minor_plans)

    def integrate(self, polynomial=None, *, allow_series=None):
        """Integrate the complete finite polynomial, including mixed degrees."""
        if allow_series is None:
            allow_series = self.allow_series
        if not isinstance(allow_series, bool):
            raise TypeError('allow_series must be boolean')
        polynomial, normalized = self._normalize_request(polynomial)
        key = tuple(sorted(polynomial.items()))
        if key in self._results:
            value, statistics = self._results[key]
            if isinstance(value, ExactForestSeries) and not allow_series:
                raise MissingForestPeriod(value.reason or 'the complete finite period is unavailable')
            self.last_stats = dict(statistics)
            return value
        statistics = {'all_original_sources_retained': True, 'kirchhoff_normalization': normalized,
                      'current_terms_after_normalization': len(polynomial)}
        if not polynomial:
            value = Fraction()
            statistics['backend'] = 'exact_kirchhoff_zero'
        elif self._forest is not None:
            value = Fraction()
            for exponent, coefficient in polynomial.items():
                word = {(i // 2, i % 2): power for i, power in enumerate(exponent) if power}
                value += coefficient * self._forest.moment(word)
            statistics.update(backend='cubic_boundary_moments', shared_subtree_messages=True)
        else:
            constant = polynomial.get((0,) * (2 * self.graph.ne)) if len(polynomial) == 1 else None
            certified, status = (None, 'not_requested')
            if constant is not None and self.use_certificates:
                certified, status = _certified_bare(self.graph, self.sources, self.marked)
            statistics['finite_certificate'] = status
            if certified is not None:
                value = constant * certified
                self.method = 'certified_five_loop_source_polynomial'
                statistics['backend'] = self.method
            else:
                self._prepare()
                pushed, contraction = self._pushforward(polynomial)
                statistics['contraction'] = contraction
                if not pushed:
                    value = Fraction()
                    statistics['backend'] = 'exact_pushforward_zero'
                elif self._ears is not None:
                    value = self._ears.integrate(pushed, allow_series=allow_series)
                    statistics.update(self._ears.last_stats)
                else:
                    if len(pushed) == 1 and (0,) * (2 * self._core.G.ne) in pushed and self.use_certificates:
                        certified, status = _certified_bare(self._core.G, self._core.sources,
                                                           self._core.marked)
                        statistics['reduced_finite_certificate'] = status
                    else:
                        certified = None
                    if certified is not None:
                        value = pushed[(0,) * (2 * self._core.G.ne)] * certified
                        statistics['backend'] = 'certified_five_loop_source_polynomial'
                    else:
                        statistics['tree_minor_plan'] = self._minor_plan()
                        value = self._core.integrate(pushed, allow_series=allow_series)
                        statistics.update(self._core.last_stats)
        self._results[key] = value, dict(statistics)
        self.last_stats = dict(statistics)
        return value

    def moment(self, powers=None, *, allow_series=None):
        return self.integrate(_word(self.graph, powers), allow_series=allow_series)

    def moments(self, words, *, allow_series=None):
        return tuple(self.moment(word, allow_series=allow_series) for word in words)

    def series(self, polynomial=None, *, powers=None):
        """Complete exact rational series; finite zeta shortcuts are bypassed."""
        if powers is not None:
            if polynomial is not None:
                raise ValueError('supply a complete polynomial or a current word, not both')
            polynomial = _word(self.graph, powers)
        polynomial, normalized = self._normalize_request(polynomial)
        key = tuple(sorted(polynomial.items()))
        if key not in self._series_results:
            if not polynomial or self._forest is not None:
                value = self.integrate(polynomial)
                series = ExactForestSeries((), (), constant=value,
                           max_edges=self.options['max_series_edges'], max_terms=self.max_terms)
                statistics = dict(self.last_stats, backend='exact_terminating_series')
            else:
                self._prepare()
                pushed, contraction = self._pushforward(polynomial)
                if self._ears is not None:
                    series = self._ears.series(polynomial=pushed)
                    statistics = dict(self._ears.last_stats, contraction=contraction)
                else:
                    plan = self._minor_plan() if pushed else None
                    series = self._core.series(pushed)
                    statistics = dict(self._core.last_stats, contraction=contraction, tree_minor_plan=plan)
                statistics.update(all_original_sources_retained=True, kirchhoff_normalization=normalized)
            self._series_results[key] = series, statistics
        value, statistics = self._series_results[key]
        self.last_stats = dict(statistics)
        return value

    def cube_series(self, powers=None, *, polynomial=None, **budgets):
        """Direct-current cube series on the ORIGINAL graph and complete sources.

        This explicit universal route constructs no global Schwinger
        numerator. It is separate from .series(), which uses the reduced
        core's existing tropical/Hepp series. No local current polynomial
        or reduced source assignment is substituted for the original input.
        """
        if powers is not None and polynomial is not None:
            raise ValueError('supply a complete polynomial or a current word, not both')
        from .forest_cubes import ForestCubes
        options = dict(check=self.check)
        options.update(budgets)
        key = tuple(sorted(options.items()))
        if key not in self._cube_engines:
            self._cube_engines[key] = ForestCubes(self.graph, self.sources, self.marked, **options)
        engine = self._cube_engines[key]
        if polynomial is not None:
            result = engine.series(polynomial=self._polynomial(polynomial))
        else:
            # Validate the public word convention before handing off.
            _word(self.graph, powers)
            result = engine.series(powers)
        self.last_stats = dict(getattr(engine, 'last_stats', {}),
                               backend='direct_current_cube_series',
                               original_graph=True, all_original_sources_retained=True)
        return result

    def tensor(self, edges, *, indices=None, allow_series=None):
        """Return every entry of the ordered current spinor tensor on edges."""
        edges = tuple(_integer(edge, 'tensor edge') for edge in edges)
        if any(edge >= self.graph.ne for edge in edges):
            raise ValueError('tensor edge lies outside the original graph')
        if len(edges) >= self.max_tensor_entries.bit_length():
            raise ForestReductionBudgetExceeded('the complete spinor tensor exceeds max_tensor_entries')
        if indices is None:
            indices = tuple('a%d' % (i + 1) for i in range(len(edges)))
        else:
            indices = tuple(indices)
        if (len(indices) != len(edges) or any(not isinstance(name, str) or not name for name in indices)
                or len(set(indices)) != len(indices)):
            raise ValueError('provide one distinct free-index name per tensor slot')
        entries = {}
        for components in product((0, 1), repeat=len(edges)):
            word = {}
            for edge, component in zip(edges, components):
                word[edge, component] = word.get((edge, component), 0) + 1
            entries[components] = self.moment(word, allow_series=allow_series)
        return ForestCurrentTensor(edges, indices, MappingProxyType(entries))


def reduce_product_forest(G, sources, marked, polynomial=None, **options):
    """Evaluate one complete current polynomial and release any period workers."""
    with ForestReduction(G, sources, marked, **options) as reduction:
        return reduction.integrate(polynomial)
