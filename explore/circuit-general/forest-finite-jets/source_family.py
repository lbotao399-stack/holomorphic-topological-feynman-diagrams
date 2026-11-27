"""Exact covariant source reconstruction on the five-loop forest topology.

Only vertices 1, 3, 5 and 6 carry sources; all three independent spinors
are arbitrary and vertex 6 enforces conservation.  The full bare kernel,
not a selected source coefficient or current projection, is reconstructed.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

from flint import fmpq, fmpq_mat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'code'))
sys.path.insert(0, str(HERE.parent / 'symmetry'))
from htk import Graph
from htk.forest_periods import CoreCurrentPeriods
from tree_minor import prepare_source_numerator

EDGES = ((0, 1), (1, 2), (0, 3), (0, 4), (2, 5), (2, 6),
         (3, 5), (3, 6), (4, 5), (4, 6), (5, 6))
MARKED = (5, 6)
# The additional zero follows structurally from the 3<->4 automorphism.  Keeping
# its candidate coefficient here provides an additional exact zero control.
BASIS = tuple((i, degree - i) for degree in (1, 3, 5)
              for i in range(degree + 1) if i != 5)
POINTS = ((1, 1), (1, -1), (3, 1), (1, 0), (2, 0), (0, 1),
          (0, 2), (0, 3), (1, 2), (2, 1), (2, 2))
TREE_CANDIDATES = ((0, 3, 4, 6, 8, 10), (0, 3, 5, 7, 9, 10), (1, 2, 5, 7, 9, 10))
SOURCE_OPTIONS = {'check': 1, 'max_terms': 750000, 'max_period_terms': 200000}
PLANNER_OPTIONS = {'max_tree_subsets': 3, 'max_plans': 3}


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def polynomial_record(terms):
    rows = sorted([[list(map(int, powers)), str(Fraction(str(coefficient)))]
                   for powers, coefficient in terms if coefficient])
    return {'nonzero_terms': len(rows), 'sha256': digest(rows)}


def bracket(first, second):
    return first[0] * second[1] - first[1] * second[0]


def sources_at(u, v):
    marked = (Fraction(u - 1, 2), Fraction(v - 1, 2))
    return ((0, 0), (1, 0), (0, 0), (0, 1), (0, 0), marked,
            (-1 - marked[0], -1 - marked[1]))


def exact_scalar_numerator(sources, *, exhaustive=False):
    """Return every scalar coefficient and its original-numerator digest."""
    core = CoreCurrentPeriods(Graph(7, EDGES), sources, **SOURCE_OPTIONS)
    try:
        plan = prepare_source_numerator(core, **({'max_tree_subsets': 1000, 'max_plans': 8}
                                                if exhaustive else
                                                dict(PLANNER_OPTIONS, tree_candidates=TREE_CANDIDATES)))
        original = polynomial_record(core._bare.terms())
        components = core._core_polynomials({(0,) * (2 * core.G.ne): Fraction(1)})
        if tuple(components) != (7,):
            raise ArithmeticError('the complete bare source family changed denominator level')
        values = components[7]
        return values, {'full_balanced_sources': [[str(value) for value in row] for row in core.sources],
                        'original_numerator': original, 'scalar_numerator': polynomial_record(values.items()),
                        'tree_plan': plan, 'scalar_edges': tuple(core._core)}
    finally:
        core.close()


def reconstruct_source_coefficients(*, progress=None):
    matrix = fmpq_mat([[u ** i * v ** j for i, j in BASIS] for u, v in POINTS])
    inverse = matrix.inv()
    samples, records = [], []
    for point in POINTS:
        sample, record = exact_scalar_numerator(sources_at(*point))
        records.append(dict(record, interpolation_point=point))
        samples.append(sample)
        if progress:
            progress(point, record)
    if any(row['scalar_edges'] != records[0]['scalar_edges'] for row in records):
        raise ArithmeticError('source sampling changed the ordered scalar core')
    support = sorted(set().union(*(set(sample) for sample in samples)))
    targets = []
    for k, (i, j) in enumerate(BASIS):
        target = {}
        for powers in support:
            coefficient = (-1) ** i * sum((Fraction(str(inverse[k, n])) * sample.get(powers, 0)
                                           for n, sample in enumerate(samples)), Fraction())
            if coefficient:
                target[powers] = coefficient
        targets.append(target)
    if targets[BASIS.index((0, 5))]:
        raise ArithmeticError('the source-free 3<->4 antisymmetry zero failed')
    rank_matrix = fmpq_mat([[fmpq(str(target.get(powers, 0))) for powers in support] for target in targets])
    _, rank = rank_matrix.rref()
    if rank != 10:
        raise ArithmeticError('the ten nonzero full source numerators are no longer independent')
    metadata = {'original_graph': {'vertices': 7, 'edges': EDGES, 'loops': 5},
                'marked': MARKED, 'source_support': (1, 3, 5, 6),
                'structurally_zero_source_vertices': (0, 2, 4),
                'source_degree': 10, 'bracket_degree': 5, 'basis': BASIS,
                'basis_meaning': '<lambda_1,lambda_3>^(5-i-j) '
                                '<lambda_3,lambda_5-lambda_6>^i '
                                '<lambda_1,lambda_5-lambda_6>^j',
                'sample_points': POINTS, 'samples': records,
                'source_options': SOURCE_OPTIONS, 'planner_options': PLANNER_OPTIONS,
                'explicit_tree_candidates': TREE_CANDIDATES,
                'scalar_edges': records[0]['scalar_edges'], 'q': 7,
                'source_coefficient_rank': rank, 'union_scalar_support': len(support),
                'coefficient_numerators': [polynomial_record(target.items()) for target in targets]}
    return tuple(targets), metadata


def evaluate_source_numerators(targets, sources):
    first, second = sources[1], sources[3]
    difference = tuple(a - b for a, b in zip(sources[5], sources[6]))
    factors = (bracket(first, second), bracket(second, difference), bracket(first, difference))
    result = {}
    for (i, j), target in zip(BASIS, targets):
        weight = factors[0] ** (5 - i - j) * factors[1] ** i * factors[2] ** j
        for powers, coefficient in target.items():
            result[powers] = result.get(powers, Fraction()) + weight * coefficient
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def fresh_source_control(targets):
    # A fresh frame point (u,v)=(3,2), followed by the invertible source
    # transformation with columns (2,1),(1,2).  det=3 is retained exactly.
    sources = ((0, 0), (2, 1), (0, 0), (1, 2), (0, 0),
               (Fraction(5, 2), 2), (Fraction(-11, 2), -5))
    actual, record = exact_scalar_numerator(sources, exhaustive=True)
    expected = evaluate_source_numerators(targets, sources)
    if actual != expected:
        raise ArithmeticError('covariant reconstruction changed a complete fresh source coefficient')
    return dict(record, equality='every complete scalar numerator coefficient',
                reconstructed_numerator=polynomial_record(expected.items()))


def five_loop_source_kernel(sources):
    """Evaluate the certified full bare polynomial for exact supported sources.

    The graph and original edge order are EDGES above.  Source rows follow
    the existing htk convention; arbitrary balanced rational spinors at
    1,3,5,6 are allowed, including all degenerate specializations.  Sources
    at 0,2,4 must vanish.  The result is Fraction or ExactPeriod, never an
    unresolved integral.  Run certify_sources.py to independently replay
    the fixed coefficient archive rather than trusting its stored contents.
    """
    from htk.forest_contractions import normalize_forest_sources
    from htk.forest_periods import ExactPeriod
    sources = normalize_forest_sources(Graph(7, EDGES), sources)
    if any(any(sources[vertex]) for vertex in (0, 2, 4)):
        raise ValueError('this certified source family requires zero sources at vertices 0,2,4')
    certificate = json.loads((HERE / 'source_certificate.json').read_text())
    if (certificate.get('schema_version') != 1 or
            certificate['source_reconstruction']['basis'] != [list(pair) for pair in BASIS] or
            certificate['source_reconstruction']['original_graph']['edges'] != [list(edge) for edge in EDGES] or
            len(certificate['finite_coefficients_in_basis_order']) != len(BASIS) or
            not certificate['zero_residuals_for_every_coefficient']):
        raise ValueError('the stored certificate does not describe this complete source family')
    difference = tuple(a - b for a, b in zip(sources[5], sources[6]))
    factors = (bracket(sources[1], sources[3]), bracket(sources[3], difference),
               bracket(sources[1], difference))
    values = {}
    for (i, j), coefficients in zip(BASIS, certificate['finite_coefficients_in_basis_order']):
        weight = factors[0] ** (5 - i - j) * factors[1] ** i * factors[2] ** j
        for label, coefficient in coefficients.items():
            values[label] = values.get(label, Fraction()) + weight * Fraction(coefficient)
    values = {label: value for label, value in values.items() if value}
    return values.get('1', Fraction()) if all(label == '1' for label in values) else ExactPeriod(values)
