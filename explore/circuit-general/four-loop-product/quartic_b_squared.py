"""Complete four-loop quartic abelian product with two original b fields.

Interaction: -i*g*c*sum(q_i*gamma_i*beta_i) + gamma1*gamma2*gamma3*gamma4.
Inputs: beta1*(partial_{+dot a} beta2)*gamma5 and b**2*beta3.
The spectator gamma5 is free and neutral. Every symbolic cubic charge
monomial, original Wick sign, surviving source, and free spinor is retained.

At fourth order, field counting requires exactly one W and three matter
vertices. The two distinguishable b occurrences represent the derivative
multiplicity of the actual unnormalized input b**2. All four W placements
are checked separately. No extra input factorial is divided out.

The complete degree-nine source polynomial has spectator degree at most
five. The 21-node triangular grid gives an exact all-source certificate;
covariant restoration uses the original three source spinors and brackets.
The complete coefficient of hbar**4 has physical multiplier
    +i*g**3/(1536*(8*pi**2)**4)
when the displayed sources are replaced by their physical derivatives.
This already includes lambda=-partial/2 and the differentiated head edge.

Default verification reassembles every saved graph value and verifies all
source coefficients, matching records, hashes, and printed metadata.
--recompute repeats every exact geometric integral. --recompute --write
regenerates the complete certificate. No period table is built.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'code'))
sys.path.insert(0, str(HERE.parent / 'nonrenormalizable-product'))
import quartic_wick as previous
from htk.graph import Graph
from htk.forest_periods import CoreCurrentPeriods

CERTIFICATE = HERE / 'four_loop_certificate.json'
NODES = tuple((a, b) for a in range(6) for b in range(6 - a))
CONTROLS = (((1, 2), (3, -1), (-2, 4)),
            ((2, -1), (1, 3), (4, -2)),
            ((1, 2), (2, 4), (3, -1)),
            ((1, 2), (0, 0), (-2, 4)))


def admissibility(edges):
    if len(set(edges)) != len(edges):
        return {'reason': 'parallel edges'}
    for size in range(2, 7):
        for vertices in combinations(range(6), size):
            count = sum(a in vertices and b in vertices for a, b in edges)
            if count > 2 * size - 3:
                return {'reason': 'Laman subset', 'vertices': vertices, 'edges': count}
    return None


def fields(flavors, W_vertex=0):
    matter = tuple(v for v in range(4) if v != W_vertex)
    result = []
    for vertex in range(4):
        if vertex == W_vertex:
            result += [('gamma', flavor, vertex, 0) for flavor in range(1, 5)]
        else:
            flavor = flavors[matter.index(vertex)]
            result += [('c', 0, vertex, 0), ('gamma', flavor, vertex, 0), ('beta', flavor, vertex, 0)]
    return tuple(result + [('beta', 1, 4, 0), ('beta', 2, 4, 0), ('gamma5', 5, 4, 0),
                           ('b', 0, 5, 0), ('b', 0, 5, 1), ('beta', 3, 5, 0)])


@lru_cache(None)
def wick_matchings(W_vertex=0):
    if isinstance(W_vertex, bool) or W_vertex not in range(4):
        raise ValueError('W_vertex must be one of the four original interaction slots')
    valid, rejected = [], []
    matter = tuple(v for v in range(4) if v != W_vertex)
    for flavors in product(range(1, 5), repeat=3):
        original = fields(flavors, W_vertex)
        antifields = tuple(field for field in original if field[0] in ('b', 'beta'))
        ordinary = tuple(field for field in original if field[0] in ('c', 'gamma'))

        def match(index, used, pairs):
            if index == len(antifields):
                yield pairs
                return
            left = antifields[index]
            for right in ordinary:
                if right in used or left[2] == right[2]:
                    continue
                if left[0] == 'b' and right[0] != 'c':
                    continue
                if left[0] == 'beta' and (right[0] != 'gamma' or left[1] != right[1]):
                    continue
                yield from match(index + 1, used | {right}, pairs + ((left, right),))

        for pairs in match(0, set(), ()):
            edges = tuple(sorted(tuple(sorted((a[2], b[2]))) for a, b in pairs) + [(4, 5)])
            record = {'W_vertex': W_vertex, 'matter_vertices': matter, 'matter_flavors': flavors,
                      'pairs': pairs, 'ordered_edges': edges}
            failure = admissibility(edges)
            if failure is not None:
                rejected.append(record | {'vanishing_witness': failure})
                continue
            remaining, sign, trace = list(original), 1, []
            for edge in reversed(edges[:-1]):
                pair = next((a, b) for a, b in pairs if tuple(sorted((a[2], b[2]))) == edge)
                odd = next(field for field in pair if field[0] in ('beta', 'c'))
                preceding = sum(field[0] in ('beta', 'c') for field in remaining[:remaining.index(odd)])
                pairing_sign = -1 if odd[0] == 'beta' else 1
                sign *= (-1) ** preceding * pairing_sign
                trace.append((edge, odd, preceding, pairing_sign))
                for field in pair:
                    remaining.remove(field)
            if Counter(field[0] for field in remaining) != Counter({'gamma': 1, 'c': 1, 'gamma5': 1}):
                raise ArithmeticError('an original surviving field was lost')
            if next(field[1] for field in remaining if field[0] == 'gamma') != 4:
                raise ArithmeticError('the surviving superpotential slot changed')
            edge = next(tuple(sorted((a[2], b[2]))) for a, b in pairs if ('beta', 2, 4, 0) in (a, b))
            if edge[1] != 4:
                raise ArithmeticError('the differentiated input was not kept at its original head')
            valid.append(record | {'wick_sign': sign, 'left_derivative_trace': tuple(trace),
                                   'surviving_fields': tuple(remaining), 'current_edge': edges.index(edge),
                                   'derivative_factor': -2})
    return tuple(valid), tuple(rejected)


def canonical_matching(record):
    candidates = []
    for order in permutations(record['matter_vertices']):
        mapping = {record['W_vertex']: 0, **{vertex: index for index, vertex in enumerate(order, 1)}, 4: 4, 5: 5}
        oriented = tuple((mapping[a], mapping[b]) for a, b in record['ordered_edges'])
        edges = tuple(sorted(tuple(sorted(edge)) for edge in oriented))
        permutation = tuple(edges.index(tuple(sorted(edge))) for edge in oriented)
        locations = tuple(mapping[next(field[2] for field in record['surviving_fields'] if field[0] == species)]
                          for species in ('gamma', 'c', 'gamma5'))
        edge = record['current_edge']
        orientation = 1 if oriented[edge][0] < oriented[edge][1] else -1
        sign = record['wick_sign'] * previous._parity(permutation) * orientation
        candidates.append(((edges, locations), permutation[edge], sign))
    return min(candidates)


@lru_cache(None)
def grouped_matchings(W_vertex=0):
    records, rejected = wick_matchings(W_vertex)
    result = defaultdict(lambda: defaultdict(lambda: defaultdict(Fraction)))
    for record in records:
        key, edge, sign = canonical_matching(record)
        result[key][tuple(sorted(record['matter_flavors']))][edge] += sign
    clean = {key: {charges: {edge: coefficient for edge, coefficient in weights.items() if coefficient}
                   for charges, weights in sorted(charge_terms.items())}
             for key, charge_terms in sorted(result.items())}
    return clean, records, rejected


def source_rows(locations, vectors):
    if len(vectors) != 3 or any(len(vector) != 2 for vector in vectors):
        raise ValueError('all three surviving source spinors are required')
    if any(isinstance(value, (float, complex, bool)) for vector in vectors for value in vector):
        raise TypeError('exact rational source values are required')
    rows = [[Fraction(0), Fraction(0)] for _ in range(6)]
    for vertex, vector in zip(locations, vectors):
        for coordinate in range(2):
            rows[vertex][coordinate] += Fraction(vector[coordinate])
    for coordinate in range(2):
        rows[5][coordinate] = -sum(rows[vertex][coordinate] for vertex in range(5))
    return tuple(map(tuple, rows))


def assemble_graphs(graphs):
    result = defaultdict(lambda: [defaultdict(Fraction), defaultdict(Fraction)])
    for row in graphs:
        for charge, weights in row['symbolic_charge_current_terms']:
            charge = tuple(charge)
            result[charge]
            for edge, weight in weights.items():
                for coordinate, value in enumerate(row['full_first_current_moments'][str(edge)]):
                    for label, coefficient in value.items():
                        result[charge][coordinate][label] += Fraction(weight) * Fraction(coefficient)
    return {charge: tuple({label: coefficient for label, coefficient in vector.items() if coefficient}
                          for vector in spinor) for charge, spinor in sorted(result.items())}


def evaluate(vectors, *, progress=False):
    groups, _, _ = grouped_matchings()
    graphs = []
    for number, ((edges, locations), charges) in enumerate(groups.items(), 1):
        rows = source_rows(locations, vectors)
        engine = CoreCurrentPeriods(Graph(6, edges), rows, marked=(4, 5), check=1)
        moments = {}
        try:
            for edge in sorted({edge for weights in charges.values() for edge in weights}):
                spinor = []
                for coordinate in range(2):
                    exponent = [0] * 18
                    exponent[2 * edge + coordinate] = 1
                    value = engine.integrate({tuple(exponent): Fraction(1)}, allow_series=False)
                    spinor.append(previous._periods(value))
                moments[str(edge)] = tuple(spinor)
        finally:
            engine.close()
        graphs.append({'ordered_edges': edges, 'survivor_vertices_gamma4_c_gamma5': locations,
                       'full_balanced_sources': rows,
                       'symbolic_charge_current_terms': [(charge, weights) for charge, weights in charges.items()],
                       'full_first_current_moments': moments})
        if progress:
            print('  complete graph', number, '/', len(groups), flush=True)
    return {'source_spinors_gamma4_c_gamma5': tuple(map(tuple, vectors)), 'graphs': graphs,
            'symbolic_charge_coefficients': [(charge, spinor) for charge, spinor in assemble_graphs(graphs).items()]}


def covariant_restore(polynomials):
    first, second = map(dict, polynomials)
    if first.get((0, 5), 0) or second.get((5, 0), 0):
        raise ArithmeticError('the complete interpolant violates spinor covariance')
    result = []
    for a in range(5):
        b = 4 - a
        coefficient = first.pop((a + 1, b), 0)
        if coefficient != second.pop((a, b + 1), 0):
            raise ArithmeticError('the full degree-five part is not along the spectator source')
        if coefficient:
            result.append((2, (0, a, b), coefficient))
    for vector, polynomial in enumerate((first, second)):
        for (a, b), coefficient in sorted(polynomial.items()):
            if a + b > 4:
                raise ArithmeticError('a nonpolynomial covariant term remained')
            if coefficient:
                result.append((vector, (4 - a - b, a, b), coefficient))
    return tuple(result)


def restore_samples(samples):
    values = {tuple(sample['source_spinors_gamma4_c_gamma5'][2]): assemble_graphs(sample['graphs']) for sample in samples}
    if set(values) != set(NODES):
        raise ArithmeticError('the full degree-five interpolation lower set is required')
    charges = sorted(next(iter(values.values())))
    if any(sorted(value) != charges for value in values.values()):
        raise ArithmeticError('a degenerate source assignment lost a symbolic charge monomial')
    labels = sorted({label for value in values.values() for spinor in value.values() for vector in spinor for label in vector} | {'z3', 'z5'})
    result = []
    for charge in charges:
        channels = {}
        for label in labels:
            polynomials = tuple(previous.interpolate({node: value[charge][coordinate].get(label, 0)
                                                       for node, value in values.items()}, degree=5)
                                for coordinate in range(2))
            restored = covariant_restore(polynomials)
            if restored:
                channels[label] = restored
        result.append((charge, channels))
    return result


def specialize(covariants, charges=(1, -1, 0, 0)):
    if len(charges) != 4:
        raise ValueError('four interacting matter charges are required')
    result = defaultdict(lambda: defaultdict(Fraction))
    for monomial, channels in covariants:
        multiplier = Fraction(1)
        for flavor in monomial:
            multiplier *= Fraction(charges[flavor - 1])
        for label, terms in channels.items():
            for vector, powers, coefficient in terms:
                result[label][vector, tuple(powers)] += multiplier * Fraction(coefficient)
    return {label: tuple((vector, powers, coefficient) for (vector, powers), coefficient in sorted(terms.items()) if coefficient)
            for label, terms in result.items() if any(terms.values())}


def implementation_hashes():
    paths = [Path(__file__), Path(previous.__file__), ROOT / 'code/htk/forest_periods.py',
             ROOT / 'code/htk/_forest_period_worker.py', ROOT / 'code/htk/spred.py',
             ROOT / 'explore/nonSP/verify/k4_periods/current_minor.py']
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    return value


def certify(*, progress=False):
    groups, records, rejected = grouped_matchings()
    if (len(records), len(rejected), len(groups)) != (324, 732, 14):
        raise ArithmeticError('the complete original matching census changed')
    for slot in (1, 2, 3):
        if grouped_matchings(slot)[0] != groups:
            raise ArithmeticError('the four superpotential placements differ')
    samples = []
    for index, node in enumerate(NODES, 1):
        samples.append(evaluate(((1, 0), (0, 1), node)))
        if progress:
            print('Complete four-loop source grid', index, '/', len(NODES), flush=True)
    covariants = restore_samples(samples)
    if any(set(channels) - {'1'} for _, channels in covariants):
        raise ArithmeticError('a nonrational term survives in a symbolic charge coefficient')
    controls = []
    for vectors in CONTROLS:
        direct = evaluate(vectors)
        restored = {tuple(charge): tuple({'1': entry} if entry else {} for entry in previous.evaluate_covariant(channels.get('1', ()), vectors))
                    for charge, channels in covariants}
        if restored != assemble_graphs(direct['graphs']):
            raise ArithmeticError('an independent full-source replay failed')
        controls.append(direct)
        if progress:
            print('Complete independent source control', len(controls), '/', len(CONTROLS), flush=True)
    return jsonable({'schema': 1, 'scope': 'complete g^3*hbar^4 coefficient, all four symbolic charges and all three surviving source spinors',
                     'inputs': ['beta1*(partial_{+dot a}beta2)*gamma5', 'b**2*beta3'],
                     'superpotential': 'gamma1*gamma2*gamma3*gamma4',
                     'free_neutral_spectator': 'gamma5, absent from the interaction',
                     'original_field_order': 'fields(matter_flavors,W_vertex), with the two b occurrences separately indexed 0 and1',
                     'matching_counts': {'valid': len(records), 'excluded': len(rejected), 'canonical_geometries': len(groups)},
                     'complete_matching_ledger': records, 'excluded_matchings_with_witnesses': rejected,
                     'W_placements_checked': [0, 1, 2, 3],
                     'source_polynomial_proof': {'total_degree': 9, 'spectator_degree_bound': 5, 'complete_lower_set': NODES,
                                                 'bound': 'rank at most four for the spectator-source part of the eight-row current determinant, plus one inserted current',
                                                 'restoration': 'det(M)^4*M and polynomial continuation to dependent source spinors'},
                     'samples': samples, 'complete_covariant_charge_polynomial': covariants,
                     'covariant_bracket_order': ['<lambda_gamma4,lambda_c>', '<lambda_gamma5,lambda_c>', '<lambda_gamma4,lambda_gamma5>'],
                     'covariant_vector_order': ['lambda_gamma4', 'lambda_c', 'lambda_gamma5'],
                     'physical_charge_specialization': [1, -1, 0, 0],
                     'complete_physical_covariant_spinor': specialize(covariants),
                     'independent_full_source_replays': controls,
                     'normalization': {'W_placements_divided_by_4_factorial': '1/6', 'graph_sign': 1, 'manuscript_kernel_divided_by_htk': 1,
                                       'three_matter_phases': 'i', 'head_derivative': -2, 'source_substitution': '(-1/2)^9',
                                       'physical_multiplier': '+i*g^3/(1536*(8*pi^2)^4)', 'input_b_squared_factorial_divided_out': False},
                     'all_symbolic_charge_period_terms_cancel': True,
                     'implementation_sha256': implementation_hashes()})


def verify_saved(data):
    for relative, expected in data['implementation_sha256'].items():
        if hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() != expected:
            raise ArithmeticError('a certified implementation changed: ' + relative)
    groups, records, rejected = grouped_matchings()
    if jsonable(records) != data['complete_matching_ledger'] or jsonable(rejected) != data['excluded_matchings_with_witnesses']:
        raise ArithmeticError('an original complete Wick ledger changed')
    for slot in (1, 2, 3):
        if grouped_matchings(slot)[0] != groups:
            raise ArithmeticError('a superpotential placement changed')
    for sample in data['samples'] + data['independent_full_source_replays']:
        vectors = sample['source_spinors_gamma4_c_gamma5']
        if len(sample['graphs']) != len(groups):
            raise ArithmeticError('an original graph was omitted')
        for graph in sample['graphs']:
            key = (tuple(map(tuple, graph['ordered_edges'])), tuple(graph['survivor_vertices_gamma4_c_gamma5']))
            if jsonable([(charge, weights) for charge, weights in groups[key].items()]) != graph['symbolic_charge_current_terms']:
                raise ArithmeticError('a complete color/current coefficient changed')
            if jsonable(source_rows(key[1], vectors)) != graph['full_balanced_sources']:
                raise ArithmeticError('an original source was changed')
        assembled = assemble_graphs(sample['graphs'])
        if jsonable(list(assembled.items())) != sample['symbolic_charge_coefficients']:
            raise ArithmeticError('a saved complete geometric sum differs')
    covariants = restore_samples(data['samples'])
    if jsonable(covariants) != data['complete_covariant_charge_polynomial']:
        raise ArithmeticError('the exact complete source interpolation changed')
    if any(set(channels) - {'1'} for _, channels in covariants):
        raise ArithmeticError('a nonrational charge coefficient was discarded')
    if jsonable(specialize(covariants)) != data['complete_physical_covariant_spinor']:
        raise ArithmeticError('the physical charge specialization changed')
    for sample in data['independent_full_source_replays']:
        vectors = sample['source_spinors_gamma4_c_gamma5']
        restored = {tuple(charge): tuple({'1': entry} if entry else {} for entry in previous.evaluate_covariant(channels.get('1', ()), vectors))
                    for charge, channels in covariants}
        if restored != assemble_graphs(sample['graphs']):
            raise ArithmeticError('a full-source polynomial replay disagrees')
    print('Verified all 324 matchings, 732 exclusions, 14 geometries, and every covariant charge coefficient.', flush=True)


def write_certificate(data):
    # Each complete ledger/sample remains explicit while avoiding whitespace-only bulk.
    compact_keys = ('complete_matching_ledger', 'excluded_matchings_with_witnesses', 'samples', 'independent_full_source_replays')
    ordinary = {key: value for key, value in data.items() if key not in compact_keys}
    text = json.dumps(ordinary, indent=2)[:-2]
    for key in compact_keys:
        text += ',\n  ' + json.dumps(key) + ': [\n'
        text += ',\n'.join('    ' + json.dumps(row, separators=(',', ':')) for row in data[key])
        text += '\n  ]'
    CERTIFICATE.write_text(text + '\n}\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recompute', action='store_true')
    parser.add_argument('--write', action='store_true')
    options = parser.parse_args()
    if options.write and not options.recompute:
        parser.error('--write requires --recompute')
    start = perf_counter()
    if options.recompute:
        data = certify(progress=True)
        if options.write:
            write_certificate(data)
            print('Wrote the complete four-loop product certificate.', flush=True)
        elif json.loads(CERTIFICATE.read_text()) != data:
            raise ArithmeticError('a complete regenerated certificate field differs')
        else:
            print('Every complete recomputation field agrees.', flush=True)
    else:
        verify_saved(json.loads(CERTIFICATE.read_text()))
    print('Elapsed %.3f seconds.' % (perf_counter() - start), flush=True)


if __name__ == '__main__':
    main()
