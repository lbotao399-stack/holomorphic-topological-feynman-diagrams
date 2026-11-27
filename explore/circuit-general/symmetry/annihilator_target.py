"""Source-target reduction by verified polynomial annihilators of U.

For homogeneous polynomial vector fields a_i and multiplier b satisfying
sum_i a_i partial_i U = b U, the existing syz5.relation supplies the complete
projective total-derivative identity, including coordinate faces.  We ask
only for identities touching the target, then eliminate its full signed
linear combination and retain an exact original-row witness.

This is the parametric-annihilator method of Bitoun, Bogner, Klausen and
Panzer, arXiv:1712.09215, applied to a complete source numerator.  The local
identity and convergence/face convention are inherited from fiveloop/syz5.py.
No period table or polynomial-annihilator basis is built automatically.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
from types import MappingProxyType

import engine
from family_target import (FamilyTarget, FamilyTargetBudgetExceeded,
                           _integer, _labels, _rational)
from relation_witness import WitnessResult, eliminate_with_witness, verify_row_witness
import syz5


class UnavailableAnnihilators(LookupError):
    """No verified polynomial archive for the requested ordered core."""


@dataclass(frozen=True)
class Annihilator:
    degree: int
    vector_fields: tuple
    scalar_multiplier: object


def _graph(edges, max_vertices=8, max_edges=14):
    maximum_vertices = _integer(max_vertices, 'max_vertices', 2)
    maximum_edges = _integer(max_edges, 'max_edges', 1)
    edges = tuple(tuple(edge) for edge in edges)
    if not edges or len(edges) > maximum_edges:
        raise FamilyTargetBudgetExceeded('input core exceeds max_edges')
    if any(len(edge) != 2 for edge in edges):
        raise ValueError('each edge must have two endpoints')
    edges = tuple(tuple(_integer(vertex, 'vertex', 0) for vertex in edge) for edge in edges)
    if any(first == second for first, second in edges):
        raise ValueError('the input core must have no self edge')
    if len({frozenset(edge) for edge in edges}) != len(edges):
        raise ValueError('the input core must be simple')
    if len({vertex for edge in edges for vertex in edge}) > maximum_vertices:
        raise FamilyTargetBudgetExceeded('input core exceeds max_vertices')
    if not engine.nx.is_connected(engine.nx.Graph(edges)):
        raise ValueError('the input core must be connected')
    return edges


def _polynomial(rows, variables, degree):
    result = {}
    for powers, coefficient in rows:
        powers = tuple(_integer(power, 'polynomial exponent', 0) for power in powers)
        if len(powers) != variables or sum(powers) != degree:
            raise ValueError('annihilator polynomial has the wrong homogeneous degree')
        result[powers] = result.get(powers, Fraction()) + _rational(coefficient)
    return MappingProxyType({powers: value for powers, value in result.items() if value})


def _normalized_generator(generator, variables):
    degree = _integer(generator.degree, 'annihilator degree', 1)
    if len(generator.vector_fields) != variables:
        raise ValueError('one vector field is required per edge')
    return Annihilator(degree,
                       tuple(_polynomial(vector.items(), variables, degree)
                             for vector in generator.vector_fields),
                       _polynomial(generator.scalar_multiplier.items(), variables, degree - 1))


def load_annihilators(path, edges, *, max_archive_bytes=2000000,
                     max_generators=128, max_generator_terms=10000):
    """Load a rational JSON archive and transport it by exact graph isomorphism.

    Loading does not establish the differential identity; callers must use
    verify_annihilators, as AnnihilatorTarget and the witness verifier do.
    JSON is used here rather than executing or unpickling external data.
    """
    edges = _graph(edges)
    path = Path(path)
    maximum_bytes = _integer(max_archive_bytes, 'max_archive_bytes', 1)
    maximum_generators = _integer(max_generators, 'max_generators', 1)
    maximum_terms = _integer(max_generator_terms, 'max_generator_terms', 1)
    if path.stat().st_size > maximum_bytes:
        raise FamilyTargetBudgetExceeded('annihilator archive exceeds max_archive_bytes')
    archive = json.loads(path.read_text(encoding='utf-8'))
    if archive.get('schema_version') != 1:
        raise ValueError('unknown annihilator archive schema')
    original_edges = _graph(archive['edges'])
    matcher = engine.nx.isomorphism.GraphMatcher(engine.nx.Graph(original_edges), engine.nx.Graph(edges))
    mapping = next(matcher.isomorphisms_iter(), None)
    if mapping is None:
        raise UnavailableAnnihilators('the archive graph is not isomorphic to the input core')
    if len(archive['generators']) > maximum_generators:
        raise FamilyTargetBudgetExceeded('annihilator archive exceeds max_generators')
    indices = {frozenset(edge): index for index, edge in enumerate(edges)}
    permutation = [indices[frozenset((mapping[u], mapping[v]))] for u, v in original_edges]

    def transport(polynomial):
        result = {}
        for powers, coefficient in polynomial.items():
            transformed = [0] * len(edges)
            for old, new in enumerate(permutation):
                transformed[new] = powers[old]
            result[tuple(transformed)] = coefficient
        return MappingProxyType(result)

    result = []
    terms = 0
    for row in archive['generators']:
        degree = _integer(row['degree'], 'annihilator degree', 1)
        if len(row['vector_fields']) != len(edges):
            raise ValueError('one vector field is required per edge')
        terms += sum(len(polynomial) for polynomial in row['vector_fields']) + len(row['scalar_multiplier'])
        if terms > maximum_terms:
            raise FamilyTargetBudgetExceeded('annihilator archive exceeds max_generator_terms')
        vector = [None] * len(edges)
        for old, new in enumerate(permutation):
            vector[new] = transport(_polynomial(row['vector_fields'][old], len(edges), degree))
        scalar = transport(_polynomial(row['scalar_multiplier'], len(edges), degree - 1))
        result.append(Annihilator(degree, tuple(vector), scalar))
    return tuple(result)


def verify_annihilators(edges, generators, *, max_polynomial_products=5000000):
    """Verify every coefficient of sum a_i partial_i U - b U is zero."""
    edges = _graph(edges)
    maximum = _integer(max_polynomial_products, 'max_polynomial_products', 1)
    polynomial = syz5.upoly(edges)
    derivatives = []
    for index in range(len(edges)):
        derivative = {}
        for powers, coefficient in polynomial.items():
            if powers[index]:
                shifted = list(powers)
                shifted[index] -= 1
                derivative[tuple(shifted)] = coefficient * powers[index]
        derivatives.append(derivative)
    products = 0
    for generator in generators:
        generator = _normalized_generator(generator, len(edges))
        if len(generator.vector_fields) != len(edges):
            raise ValueError('annihilator vector field has the wrong edge count')
        cost = sum(len(vector) * len(derivative)
                   for vector, derivative in zip(generator.vector_fields, derivatives))
        cost += len(generator.scalar_multiplier) * len(polynomial)
        if products + cost > maximum:
            raise FamilyTargetBudgetExceeded('annihilator verification exceeds max_polynomial_products')
        difference = Counter()
        for vector, derivative in zip(generator.vector_fields, derivatives):
            for first, coefficient in vector.items():
                for second, value in derivative.items():
                    difference[tuple(a + b for a, b in zip(first, second))] += coefficient * value
        for first, coefficient in generator.scalar_multiplier.items():
            for second, value in polynomial.items():
                difference[tuple(a + b for a, b in zip(first, second))] -= coefficient * value
        if any(difference.values()):
            raise ArithmeticError('a supplied polynomial does not annihilate the graph denominator')
        products += cost
    return {'generators': len(generators), 'kirchhoff_terms': len(polynomial),
            'verified_polynomial_products': products}


def _items(edges, generator, beta, q):
    if len(beta) != len(edges) or any(not isinstance(power, int) or isinstance(power, bool)
                                     or power < 1 for power in beta):
        raise ValueError('annihilator base powers must be positive integers')
    loops = engine.gutil.cycle_rank(edges)
    if sum(beta) != loops * q - generator.degree + 1:
        raise ValueError('annihilator base has the wrong projective degree')
    return syz5.relation(edges, (generator.vector_fields, generator.scalar_multiplier), list(beta), q)


class AnnihilatorTarget:
    """Generate complete identities touching a whole scalar target.

    Candidate/storage limits stop generation explicitly.  Every already
    registered input coefficient remains in the target; rejected identities
    are not partially used.  The subsequent exact solver returns a finite
    value only if all residual coefficients cancel.
    """

    def __init__(self, edges, terms, q, generators, *, max_input_terms=10000,
                 max_variables=30000, max_rows=100000,
                 max_row_coefficients=3000000, max_images=256,
                 max_polynomial_products=5000000, max_generator_terms=10000):
        self.edges = _graph(edges)
        generators = tuple(generators)
        maximum_terms = _integer(max_generator_terms, 'max_generator_terms', 1)
        if sum(sum(len(vector) for vector in generator.vector_fields)
               + len(generator.scalar_multiplier) for generator in generators) > maximum_terms:
            raise FamilyTargetBudgetExceeded('generators exceed max_generator_terms')
        self.generators = tuple(_normalized_generator(generator, len(self.edges)) for generator in generators)
        self.verification = verify_annihilators(self.edges, self.generators,
                                                max_polynomial_products=max_polynomial_products)
        self.plan = FamilyTarget(self.edges, terms, q, max_input_terms=max_input_terms,
                                 max_variables=max_variables, max_rows=max_rows,
                                 max_row_coefficients=max_row_coefficients)
        self.plan.system.counts['annihilator'] = 0
        self.maximum_images = _integer(max_images, 'max_images', 1)
        if len(self.plan.system.fam[0]['autos']) > self.maximum_images:
            raise FamilyTargetBudgetExceeded('graph automorphisms exceed max_images')
        self.row_origins = []
        self.generation = {'candidate_bases': 0, 'convergence_rejections': 0,
                           'classification_rejections': 0, 'accepted_rows': 0,
                           'generation_stopped_by': None}
        self.generated = False

    def generate(self, *, max_candidates=100000, include_symmetries=True):
        maximum = _integer(max_candidates, 'max_candidates', 0)
        if not isinstance(include_symmetries, bool):
            raise ValueError('include_symmetries must be boolean')
        if self.generated:
            raise ValueError('candidate generation has already run on this target')
        self.generated = True
        shifts = []
        for generator in self.generators:
            current = set(generator.scalar_multiplier)
            for index, vector in enumerate(generator.vector_fields):
                for powers in vector:
                    shifted = list(powers)
                    shifted[index] -= 1
                    current.add(tuple(shifted))
            shifts.append(tuple(sorted(current)))
        system = self.plan.system
        seen = set()
        for column in self.plan.input_columns:
            target = system.reps[column][1]
            images = {tuple(target[index] for index in permutation)
                      for permutation in system.fam[0]['autos']} if include_symmetries else {target}
            for powers in sorted(images):
                for identifier, generator in enumerate(self.generators):
                    for shift in shifts[identifier]:
                        beta = tuple(power - change for power, change in zip(powers, shift))
                        key = (identifier, beta)
                        if min(beta) < 1 or key in seen:
                            continue
                        if len(seen) >= maximum:
                            self.generation['generation_stopped_by'] = 'max_candidates'
                            return MappingProxyType(dict(self.generation))
                        seen.add(key)
                        self.generation['candidate_bases'] += 1
                        items = _items(self.edges, generator, beta, self.plan.q)
                        if items is None:
                            self.generation['convergence_rejections'] += 1
                            continue
                        before = len(system.rows)
                        try:
                            accepted = system.add_relation(items, 'annihilator')
                        except FamilyTargetBudgetExceeded as error:
                            self.generation['generation_stopped_by'] = str(error)
                            return MappingProxyType(dict(self.generation))
                        if not accepted:
                            self.generation['classification_rejections'] += 1
                        elif len(system.rows) > before:
                            self.row_origins.append((identifier, beta))
                            self.generation['accepted_rows'] += 1
        return MappingProxyType(dict(self.generation))

    def reduce(self, **budgets):
        result = eliminate_with_witness(self.plan.system.rows, self.plan.target,
                                        [powers for _, powers in self.plan.system.reps],
                                        self.edges, self.plan.q, **budgets)
        statistics = dict(result.statistics)
        statistics.update(self.generation)
        statistics.update(self.verification)
        statistics.update(input_terms=self.plan.input_terms,
                          input_orbits=len(self.plan.input_columns), read_only_tables=True)
        return WitnessResult(result.finite_coefficients, result.residuals,
                             MappingProxyType(statistics), result.row_weights)

    def witness_records(self, result):
        """Return (generator index, positive base powers, exact row weight)."""
        return tuple((self.row_origins[row][0], self.row_origins[row][1], weight)
                     for row, weight in result.row_weights.items())


def verify_annihilator_witness(edges, terms, q, generators, witness, *,
                               max_witness_rows=20000, max_variables=30000,
                               max_row_coefficients=3000000):
    """Rebuild only the witnessed complete identities and sum them exactly.

    Returns the complete finite scalar dictionary.  If even one target
    coefficient is left, or a boundary cannot be evaluated, verification
    raises; no search, elimination, or unrequested table build occurs.
    """
    maximum = _integer(max_witness_rows, 'max_witness_rows', 1)
    plan = AnnihilatorTarget(edges, terms, q, generators, max_variables=max_variables,
                             max_rows=maximum, max_row_coefficients=max_row_coefficients)
    system = plan.plan.system
    weights = {}
    finite = {}
    for count, (identifier, beta, coefficient) in enumerate(witness, 1):
        if count > maximum:
            raise FamilyTargetBudgetExceeded('witness exceeds max_witness_rows')
        identifier = _integer(identifier, 'generator index', 0)
        if identifier >= len(plan.generators):
            raise ValueError('witness generator index is out of range')
        coefficient = _rational(coefficient)
        if not coefficient:
            raise ValueError('a witness row weight must be nonzero')
        beta = tuple(beta)
        items = _items(plan.edges, plan.generators[identifier], beta, plan.plan.q)
        if items is None:
            raise ArithmeticError('a witness derivative is not absolutely integrable')
        before = len(system.rows)
        if not system.add_relation(items, 'annihilator'):
            raise ArithmeticError('a complete witness boundary period is unavailable')
        if len(system.rows) != before + 1:
            raise ArithmeticError('a witness row is empty')
        weights[before] = coefficient
        for label, value in system.rows[before][1].items():
            finite[label] = finite.get(label, Fraction()) + coefficient * value
    finite = {label: value for label, value in finite.items() if value}
    verify_row_witness(system.rows, plan.plan.target, {}, finite, weights)
    return MappingProxyType(finite)
