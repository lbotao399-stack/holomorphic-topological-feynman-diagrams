"""Verify the complete finite five-loop value through an exact row witness.

Default execution reconstructs the original source numerator, verifies every
archived polynomial annihilator, rebuilds the witnessed complete identities,
and sums them exactly.  --recompute repeats bounded candidate generation and
elimination too.  --recompute --write creates the compact certificate.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'code'))
sys.path.insert(0, str(HERE.parent / 'symmetry'))

from htk import ExactForestSeries, Graph, ProductForestMoments
from annihilator_target import (AnnihilatorTarget, load_annihilators,
                                 verify_annihilator_witness, verify_annihilators)


IDENTITY_PATH = HERE.parent / 'forest-target-reduction/dimension_certificate.json'
# Immutable mathematical input from commit 2aa881f.  The whole 17-term
# expression is used, including the previously evaluated finite contribution.
IDENTITY_SHA256 = 'fa765ab202321680e708120baafd524d91fb58abce7f113a43ef2af68e3d9d12'
ARCHIVE = HERE / 'annihilators.json'
CERTIFICATE = HERE / 'finite_certificate.json'
CONSTRUCTION = {'max_input_terms': 10000, 'max_variables': 30000,
                'max_rows': 100000, 'max_row_coefficients': 3000000,
                'max_images': 256, 'max_polynomial_products': 5000000,
                'max_generator_terms': 10000}
GENERATION = {'max_candidates': 100000, 'include_symmetries': True}
ELIMINATION = {'max_operations': 10000000, 'max_active_coefficients': 3000000,
               'max_pivots': 30000, 'max_history_entries': 2000000}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest_terms(terms):
    rows = sorted([[list(map(int, powers)), str(Fraction(str(value)))]
                   for powers, value in terms if value])
    return hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()


def input_identity():
    if sha(IDENTITY_PATH) != IDENTITY_SHA256:
        raise ArithmeticError('the immutable complete input identity changed')
    identity = json.loads(IDENTITY_PATH.read_text())
    edges = tuple(map(tuple, identity['scalar_core']['edges']))
    q = identity['scalar_core']['q']
    terms = [(tuple(powers), Fraction(coefficient))
             for powers, coefficient in identity['all_residual_scalar_periods']]
    if q != 7 or len(terms) != 17:
        raise ArithmeticError('the complete five-loop input identity changed scope')
    return identity, edges, q, terms


def check_original_sources(identity):
    graph_data = identity['original_graph']
    graph = Graph(graph_data['vertices'], graph_data['edges'], 'finite-five-loop-source-check')
    with ProductForestMoments(graph, identity['full_balanced_sources'], identity['marked'],
                              **identity['source_options']) as forest:
        series = forest.series()
        if not isinstance(series, ExactForestSeries) or len(series.components) != 1 or series.constant:
            raise ArithmeticError('the original complete kernel changed scalar components')
        q, numerator = series.components[0]
        if q != identity['scalar_core']['q'] or list(map(list, series.edges)) != identity['scalar_core']['edges']:
            raise ArithmeticError('the ordered original scalar core changed')
        original_digest = digest_terms(forest._core._bare.terms())
        scalar_digest = digest_terms(numerator.items())
    if original_digest != identity['original_numerator']['sha256']:
        raise ArithmeticError('an original source-numerator coefficient changed')
    if scalar_digest != identity['scalar_core']['sha256']:
        raise ArithmeticError('a complete scalar-numerator coefficient changed')
    return {'original_numerator_sha256': original_digest,
            'scalar_numerator_sha256': scalar_digest,
            'original_terms': identity['original_numerator']['nonzero_terms'],
            'scalar_terms': identity['scalar_core']['nonzero_terms']}


def total_value(identity, remainder):
    value = {label: Fraction(coefficient) for label, coefficient in identity['finite_coefficients'].items()}
    for label, coefficient in remainder.items():
        value[label] = value.get(label, Fraction()) + coefficient
    return {label: coefficient for label, coefficient in value.items() if coefficient}


def positive_interval(value):
    # Integral test after the first two positive terms:
    # 1+1/8+1/18 <= zeta(3) <= 1+1/8+1/8.
    if set(value) != {'1', 'z3'} or value['z3'] >= 0:
        raise ArithmeticError('the expected two-coefficient finite value changed')
    zeta_lower, zeta_upper = Fraction(85, 72), Fraction(5, 4)
    lower = value['1'] + value['z3'] * zeta_upper
    upper = value['1'] + value['z3'] * zeta_lower
    if not 0 < lower <= upper:
        raise ArithmeticError('the exact interval does not certify positivity')
    return {'zeta3_lower': str(zeta_lower), 'zeta3_upper': str(zeta_upper),
            'kernel_lower': str(lower), 'kernel_upper': str(upper),
            'strictly_positive': True,
            'proof': 'Sum 1/k^3 for k=1,2, then bound the remaining sum between '
                     'the integrals from 3 and from 2 to infinity; multiply by the '
                     'negative zeta3 coefficient and add the exact rational coefficient.'}


def implementation_hashes():
    paths = {Path(__file__).resolve(), ROOT / 'code/htk/_forest_period_worker.py',
             ROOT / 'explore/nonSP/verify/k4_periods/current_minor.py'}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT):
                paths.add(path)
    return {str(path.relative_to(ROOT)): sha(path) for path in sorted(paths)}


def recompute():
    started = perf_counter()
    identity, edges, q, terms = input_identity()
    source_check = check_original_sources(identity)
    print('Reconstructed every original source coefficient: %d original, %d scalar (%.3f s).' %
          (source_check['original_terms'], source_check['scalar_terms'], perf_counter() - started), flush=True)
    generators = load_annihilators(ARCHIVE, edges)
    plan = AnnihilatorTarget(edges, terms, q, generators, **CONSTRUCTION)
    plan.generate(**GENERATION)
    print('Generated %d complete annihilator rows on %d scalar variables (%.3f s).' %
          (len(plan.plan.system.rows), len(plan.plan.system.reps), perf_counter() - started), flush=True)
    result = plan.reduce(**ELIMINATION)
    if not result.complete:
        raise ArithmeticError('finite closure failed; %d exact residuals remain: %s' %
                              (len(result.residuals), dict(result.statistics)))
    witness = plan.witness_records(result)
    checked_remainder = verify_annihilator_witness(edges, terms, q, generators, witness)
    if dict(checked_remainder) != dict(result.value):
        raise ArithmeticError('independent witness reconstruction disagrees with elimination')
    value = total_value(identity, checked_remainder)
    positivity = positive_interval(value)
    print('Complete finite closure: %s; %d independently summable witness rows (%.3f s).' %
          ({label: str(coefficient) for label, coefficient in value.items()}, len(witness),
           perf_counter() - started), flush=True)
    archive = json.loads(ARCHIVE.read_text())
    cache = ROOT / archive['source']['path']
    if sha(cache) != archive['source']['sha256']:
        raise ArithmeticError('the historical annihilator cache changed')
    data = {
        'schema_version': 1,
        'scope': identity['scope'],
        'convention': 'original htk ordered oriented graph; full balanced sources; bare current word',
        'original_graph': identity['original_graph'],
        'full_balanced_sources': identity['full_balanced_sources'],
        'marked': identity['marked'],
        'input_identity': {'path': str(IDENTITY_PATH.relative_to(ROOT)), 'sha256': IDENTITY_SHA256,
                           'complete_residual_periods_used': len(terms),
                           'previous_finite_coefficients': identity['finite_coefficients']},
        'scalar_core': identity['scalar_core'],
        'fresh_original_source_check': source_check,
        'annihilator_archive': {'path': str(ARCHIVE.relative_to(ROOT)), 'sha256': sha(ARCHIVE),
                                'historical_source': archive['source'],
                                'exact_polynomial_verification': plan.verification},
        'construction_options': CONSTRUCTION,
        'generation_options': GENERATION,
        'elimination_options': ELIMINATION,
        'statistics': dict(result.statistics),
        'remainder_finite_coefficients': {label: str(coefficient)
                                          for label, coefficient in sorted(checked_remainder.items())},
        'complete_kernel_coefficients': {label: str(coefficient) for label, coefficient in sorted(value.items())},
        'complete_finite_reduction': True,
        'all_residual_scalar_periods': [],
        'positive_rational_interval': positivity,
        'witness_convention': 'Each row [generator index, positive beta vector, rational weight] '
                              'multiplies the COMPLETE syz5.relation, including every coordinate face. '
                              'Their exact sum is the entire archived 17-term scalar target.',
        'witness': [[identifier, list(beta), str(coefficient)]
                    for identifier, beta, coefficient in witness],
    }
    data['implementation_sha256'] = implementation_hashes()
    paths = [HERE.parent / 'cores/tables/prism_q7.json', HERE.parent / 'cores/tables/K33_q7.json',
             HERE.parent / 'k5/tables/K5_q7.json']
    data['period_table_sha256'] = {str(path.relative_to(ROOT)): sha(path) for path in paths}
    return data


def verify_saved(data):
    identity, edges, q, terms = input_identity()
    for collection in ('implementation_sha256', 'period_table_sha256'):
        for relative, expected in data[collection].items():
            if sha(ROOT / relative) != expected:
                raise ArithmeticError('a certified implementation/table file changed: %s' % relative)
    if sha(ARCHIVE) != data['annihilator_archive']['sha256']:
        raise ArithmeticError('the certified polynomial archive changed')
    source_check = check_original_sources(identity)
    if source_check != data['fresh_original_source_check']:
        raise ArithmeticError('the fresh complete source reconstruction differs')
    generators = load_annihilators(ARCHIVE, edges)
    if verify_annihilators(edges, generators) != data['annihilator_archive']['exact_polynomial_verification']:
        raise ArithmeticError('the exact polynomial-verification record changed')
    witness = [(identifier, tuple(beta), Fraction(coefficient))
               for identifier, beta, coefficient in data['witness']]
    remainder = verify_annihilator_witness(edges, terms, q, generators, witness)
    value = total_value(identity, remainder)
    if {label: str(coefficient) for label, coefficient in remainder.items()} != data['remainder_finite_coefficients']:
        raise ArithmeticError('the witnessed finite remainder changed')
    if {label: str(coefficient) for label, coefficient in value.items()} != data['complete_kernel_coefficients']:
        raise ArithmeticError('the complete finite kernel changed')
    if positive_interval(value) != data['positive_rational_interval']:
        raise ArithmeticError('the rigorous positivity certificate changed')
    if not data['complete_finite_reduction'] or data['all_residual_scalar_periods']:
        raise ArithmeticError('the finite certificate incorrectly retains a period')
    print('Verified all %d complete row identities and every original source coefficient.' % len(witness), flush=True)
    print('Complete kernel: %s.' % data['complete_kernel_coefficients'], flush=True)
    print('Strictly positive rational interval: [%s, %s].' %
          (data['positive_rational_interval']['kernel_lower'],
           data['positive_rational_interval']['kernel_upper']), flush=True)


def write_certificate(data):
    # Keep one full witness row per line; every coefficient remains visible.
    metadata = dict(data)
    witness = metadata.pop('witness')
    text = json.dumps(metadata, indent=2, sort_keys=True)
    text = text[:-2] + ',\n  "witness": [\n'
    text += ',\n'.join('    ' + json.dumps(row, separators=(',', ':')) for row in witness)
    text += '\n  ]\n}\n'
    CERTIFICATE.write_text(text, encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recompute', action='store_true')
    parser.add_argument('--write', action='store_true')
    options = parser.parse_args()
    if options.write and not options.recompute:
        parser.error('--write requires --recompute')
    started = perf_counter()
    if options.recompute:
        data = json.loads(json.dumps(recompute()))
        if options.write:
            write_certificate(data)
            print('Wrote the complete finite certificate.', flush=True)
        elif data != json.loads(CERTIFICATE.read_text()):
            raise ArithmeticError('a bounded-recomputation certificate field differs')
        else:
            print('Every bounded-recomputation certificate field replays.', flush=True)
    else:
        verify_saved(json.loads(CERTIFICATE.read_text()))
    print('Elapsed %.3f seconds.' % (perf_counter() - started), flush=True)


if __name__ == '__main__':
    main()
