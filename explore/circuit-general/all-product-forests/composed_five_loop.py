"""Complete physical five-loop scalar and free-spinor composition certificate.

Run without arguments to reconstruct and replay every exact field.  --write
replaces the certificate after all five independently composed integrals pass.
The full sources are fixed explicitly; this is not an all-source formula or a
Wick/color sum.  Spinors retain every entry and are also expressed in the
independent original source basis lambda_0, lambda_3.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))

from htk import Graph
from htk.forest_contractions import ForestContractions
from htk.forest_ears import ForestEarMoments
from htk.forest_model import recognize_product_forest
from htk.forest_periods import ExactPeriod
from htk.forest_reduction import ForestReduction


EDGES = ((0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4),
         (5, 3), (6, 1), (6, 5), (6, 3))
SOURCES = ((1, 2), (3, -1), (0, 0), (5, 3), (-11, -1), (2, -3), (0, 0))
MARKED = (3, 4)
CURRENT_EDGE = EDGES.index((6, 1))
CERTIFICATE = HERE / 'composed_five_loop_certificate.json'
BUDGETS = {'check': 2, 'max_terms': 250000, 'max_source_samples': 10000,
           'max_period_terms': 250000, 'max_matrix_entries': 1000000,
           'period_timeout': 120}


def native(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [native(item) for item in value]
    return value


def finite(value):
    if isinstance(value, (int, Fraction)):
        return {'1': str(value)} if value else {}
    if not isinstance(value, ExactPeriod):
        raise ArithmeticError('the complete five-loop expression retained an unresolved period')
    return {label: str(coefficient) for label, coefficient in sorted(value.items()) if coefficient}


def add(records):
    result = {}
    for weight, coefficients in records:
        for label, coefficient in coefficients.items():
            result[label] = result.get(label, Fraction()) + weight * Fraction(coefficient)
    return {label: str(value) for label, value in sorted(result.items()) if value}


def source_basis(entries):
    """An invertible exact two-spinor change of basis, retaining both entries."""
    first, second = SOURCES[0], SOURCES[3]
    determinant = first[0] * second[1] - first[1] * second[0]
    if not determinant:
        raise ArithmeticError('the requested original source basis is singular')
    coefficients = [add(((Fraction(second[1], determinant), entries[0]),
                          (Fraction(-second[0], determinant), entries[1]))),
                    add(((Fraction(-first[1], determinant), entries[0]),
                          (Fraction(first[0], determinant), entries[1])))]
    for entry in range(2):
        if add(((first[entry], coefficients[0]), (second[entry], coefficients[1]))) != entries[entry]:
            raise ArithmeticError('the full source-spinor reconstruction failed')
    return {'lambda_0,+dot_a': coefficients[0], 'lambda_3,+dot_a': coefficients[1]}


def polynomial_record(polynomial):
    return [[list(exponent), str(value)] for exponent, value in sorted(polynomial.items()) if value]


def insertions(graph):
    zero = (0,) * (2 * graph.ne)
    result = [('bare', {zero: Fraction(1)})]
    for entry in range(2):
        exponent = [0] * (2 * graph.ne)
        exponent[2 * CURRENT_EDGE + entry] = 1
        result.append(('current_' + str(entry), {tuple(exponent): Fraction(1)}))
    for entry in range(2):
        exponent = [0] * (2 * graph.ne)
        exponent[2 * CURRENT_EDGE + entry] = 1
        result.append(('mixed_' + str(entry), {zero: Fraction(SOURCES[0][entry]),
                                               tuple(exponent): Fraction(-2, 7)}))
    return result


def implementation_hashes():
    htk = ('forest_reduction', 'forest_contractions', 'forest_ears', 'forest_periods',
           '_forest_period_worker', 'forest_model', 'forest_moments', 'graph', 'linalg',
           'brute', 'treeflow', 'star', 'spred')
    paths = [ROOT / ('code/htk/' + name + '.py') for name in htk]
    paths += list((ROOT / 'explore/circuit-general/cores').glob('*.py'))
    paths += [ROOT / ('explore/circuit-general/fiveloop/' + name)
              for name in ('per5.py', 'msys5.py', 'skel5.json')]
    paths += [ROOT / 'explore/nonSP/verify/k4_periods/current_minor.py',
              ROOT / 'explore/circuit-general/symmetry/tree_minor.py']
    paths.append(Path(__file__).resolve())
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(set(paths))}


def calculate():
    graph = Graph(7, EDGES)
    structure = recognize_product_forest(graph, MARKED)
    if structure.higher_valence_vertices != ((1,),) or structure.degree_two_vertices != ((0, 5),):
        raise ArithmeticError('the physical superpotential/source topology changed')
    if any(SOURCES[vertex] != (0, 0) for vertex in (2, 6)):
        raise ArithmeticError('a source-free matter/gauge vertex was sourced')
    requests = insertions(graph)
    direct_values, direct_stats = {}, {}
    with ForestReduction(graph, SOURCES, MARKED, allow_series=False, use_certificates=False,
                         **BUDGETS) as reduction:
        direct_values['bare'] = finite(reduction.moment())
        direct_stats['bare'] = native(reduction.last_stats)
        tensor = reduction.tensor((CURRENT_EDGE,), indices=('dot_a',), allow_series=False)
        if set(tensor.entries) != {(0,), (1,)}:
            raise ArithmeticError('the complete free-spinor tensor was not retained')
        for entry in range(2):
            direct_values['current_' + str(entry)] = finite(tensor[(entry,)])
            reduction.moment({(CURRENT_EDGE, entry): 1}, allow_series=False)
            direct_stats['current_' + str(entry)] = native(reduction.last_stats)
        for name, polynomial in requests[3:]:
            direct_values[name] = finite(reduction.integrate(polynomial, allow_series=False))
            direct_stats[name] = native(reduction.last_stats)
        method = reduction.method
        if method != 'sourcefree_stars_then_sourced_ears':
            raise ArithmeticError('the intended exact composition was bypassed')

    # This reference does not use ForestReduction, its Kirchhoff normal form,
    # certificate shortcut, or dispatcher.  Existing local identities are
    # composed directly, retaining every homogeneous group and coefficient.
    contractions = ForestContractions(graph, SOURCES, MARKED,
                                       check=BUDGETS['check'], max_terms=BUDGETS['max_terms'])
    independent = {}
    with ForestEarMoments(contractions.reduced_graph, contractions.reduced_sources,
                          contractions.reduced_marked, **BUDGETS) as ears:
        for name, polynomial in requests:
            groups = {}
            for exponent, coefficient in polynomial.items():
                groups.setdefault(sum(exponent), {})[exponent] = coefficient
            pushed, records = {}, []
            for degree, group in sorted(groups.items()):
                contraction = contractions.contract(polynomial=group)
                records.append({'degree': degree, 'statistics': native(contraction.stats),
                                'complete_polynomial': polynomial_record(contraction.polynomial)})
                for exponent, coefficient in contraction.polynomial.items():
                    pushed[exponent] = pushed.get(exponent, Fraction()) + coefficient
            pushed = {exponent: coefficient for exponent, coefficient in pushed.items() if coefficient}
            value = finite(ears.integrate(pushed, allow_series=False))
            if value != direct_values[name]:
                raise ArithmeticError('independent complete composition failed: ' + name)
            independent[name] = {'value': value, 'homogeneous_pushforwards': records,
                                 'complete_recombined_polynomial': polynomial_record(pushed),
                                 'ear_statistics': native(ears.last_stats)}
        ear_steps = native(ears.steps)
        if len(ears._moves) != 1 or ears._moves[0].k != 1:
            raise ArithmeticError('the recorded single-ear affine source transfer changed')
        move = ears._moves[0]
        constant = move.shifted_sources(contractions.reduced_sources, (Fraction(0),))
        endpoint = move.shifted_sources(contractions.reduced_sources, (Fraction(1),))
        slopes = tuple(tuple(b - a for a, b in zip(left, right)) for left, right in zip(constant, endpoint))
        terminal = {'vertices': ears.reduced_graph.nv, 'edges': native(ears.reduced_graph.edges),
                    'marked': native(ears.reduced_marked),
                    'source_transfer': 'constant + ear_ratio * slope, integrated exactly; no midpoint source substitution',
                    'full_source_constants': native(constant), 'full_source_slopes': native(slopes)}

    current = [direct_values['current_' + str(entry)] for entry in range(2)]
    mixed = [direct_values['mixed_' + str(entry)] for entry in range(2)]
    for entry in range(2):
        expected = add(((SOURCES[0][entry], direct_values['bare']), (Fraction(-2, 7), current[entry])))
        if mixed[entry] != expected:
            raise ArithmeticError('complete mixed-degree spinor linearity failed')
    return {'schema_version': 1,
            'scope': 'Complete fixed-source five-loop physical product-forest scalar, original-edge current spinor, and mixed-degree spinor; one degree-five superpotential at vertex1, no source restriction beyond the explicitly recorded physical assignment; not a full Wick/color sum.',
            'convention': 'htk H=2,T=0 geometric kernel in the original ordered/oriented edges, before physical endpoint factors or displacement factorials',
            'original_graph': {'vertices': graph.nv, 'edges': native(EDGES), 'marked': native(MARKED)},
            'full_balanced_original_sources': native(SOURCES),
            'superpotential': {'vertex': 1, 'contracted_valence': 4, 'degree': 5},
            'degree_two_source_vertices': [0, 5], 'source_free_internal_vertices': [2, 6],
            'budgets': BUDGETS, 'method': method,
            'bare_kernel': direct_values['bare'],
            'original_current_spinor': {'oriented_edge': [6, 1], 'original_edge_index': CURRENT_EDGE,
                'free_index': '+dot_a', 'full_entries': current,
                'complete_original_source_basis': source_basis(current)},
            'mixed_degree_spinor': {'insertion': 'lambda_0,+dot_a - (2/7) y_(6,1),+dot_a',
                'free_index': '+dot_a', 'full_entries': mixed,
                'complete_original_source_basis': source_basis(mixed)},
            'source_basis_determinant': '-7',
            'all_original_sources_retained': True, 'every_free_spinor_entry_retained': True,
            'exact_comparison_integrals': 5, 'finite_unresolved_remainder': [],
            'original_request_polynomials': {name: polynomial_record(value) for name, value in requests},
            'dispatcher_statistics': direct_stats,
            'independent_existing_rule_composition': independent,
            'passive_stages': native(contractions.steps),
            'passive_remainder': {'vertices': contractions.reduced_graph.nv,
                'edges': native(contractions.reduced_graph.edges), 'marked': native(contractions.reduced_marked),
                'full_sources': native(contractions.reduced_sources)},
            'sourced_ear_stages': ear_steps, 'terminal_affine_source_graph': terminal,
            'implementation_sha256': implementation_hashes()}


def replay(path=CERTIFICATE):
    result = calculate()
    if result != json.loads(Path(path).read_text()):
        raise ArithmeticError('the complete physical five-loop composition certificate differs')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    arguments = parser.parse_args()
    start = perf_counter()
    if arguments.write:
        result = calculate()
        CERTIFICATE.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    else:
        result = replay()
    print('All five complete physical five-loop integrals agree; unresolved remainder is empty.')
    print('Bare:', result['bare_kernel'])
    print('Full current spinor in original source basis:',
          result['original_current_spinor']['complete_original_source_basis'])
    print('Full mixed spinor in original source basis:',
          result['mixed_degree_spinor']['complete_original_source_basis'])
    print('Complete %s finished in %.3f s.' % ('generation' if arguments.write else 'replay', perf_counter() - start))


if __name__ == '__main__':
    main()
