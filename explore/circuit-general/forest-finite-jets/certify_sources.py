"""Reconstruct and certify the full covariant five-loop bare source kernel.

Default regenerates the original source numerators and every inherited row,
then directly sums the stored original-row witnesses.  --write additionally
runs the shared exact elimination and writes a deterministic certificate.
No full scalar period level or table is constructed.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter

from source_family import (HERE, ROOT, BASIS, EDGES, MARKED, digest, polynomial_record,
                           reconstruct_source_coefficients, evaluate_source_numerators,
                           fresh_source_control, sources_at)
from htk.forest_periods import _research_module
from family_target import FamilyTarget, _labels
from annihilator_target import load_annihilators, AnnihilatorTarget
from batch_witness import eliminate_batch_with_witness
from relation_witness import verify_row_witness

CERTIFICATE = HERE / 'source_certificate.json'
WITNESS = HERE / 'source_witnesses.json.gz'
CONSTRUCTION = {'max_input_terms': 20000, 'max_variables': 70000,
                'max_rows': 260000, 'max_row_coefficients': 5000000}
ELIMINATION = {'max_operations': 100000000, 'max_active_coefficients': 5000000,
               'max_pivots': 70000, 'max_history_entries': 6000000}
REFERENCE_FILES = ('forest-target-reduction/dimension_certificate.json',
                   'forest-finite-target/annihilators.json',
                   'forest-finite-target/finite_certificate.json')


def hashes(paths):
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(set(paths))}


def implementation_hashes():
    paths = {Path(__file__).resolve(),
             ROOT / 'explore/nonSP/verify/k4_periods/current_minor.py',
             HERE.parent / 'forest-target-reduction/certify_target.py'}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT):
                paths.add(path)
    return hashes(paths)


def row_digest(representatives, rows):
    value = hashlib.sha256()
    value.update(json.dumps(representatives, separators=(',', ':')).encode())
    for row, finite in rows:
        value.update(b'\n')
        record = [sorted([[int(index), str(Fraction(coefficient))] for index, coefficient in row.items()]),
                  sorted([[label, str(Fraction(coefficient))] for label, coefficient in _labels(finite).items()])]
        value.update(json.dumps(record, separators=(',', ':')).encode())
    return value.hexdigest()


def inherited_system(targets, edges, progress):
    source = evaluate_source_numerators(targets, sources_at(1, 1))
    previous = json.loads((HERE.parent / REFERENCE_FILES[0]).read_text())
    exponent_record = polynomial_record((tuple(power - 1 for power in alpha), coefficient)
                                        for alpha, coefficient in source.items())
    if exponent_record != {key: previous['scalar_core'][key] for key in ('nonzero_terms', 'sha256')}:
        raise ArithmeticError('the old full source target no longer matches its immutable input')
    plan = FamilyTarget(edges, sorted(source.items()), 7, **CONSTRUCTION)
    plan.seed_finite(search_states=40)
    progress('Regenerated %d complete finite seeds.' % plan.system.counts['finite_seed'])
    plan.close_local(x_edge_limits={9: 4})
    progress('Closed %d local scalar variables.' % len(plan.system.reps))
    old = _research_module('explore/circuit-general/forest-target-reduction/certify_target.py',
                           '_htk_finite_source_dimension_reference')
    plan.add_dimension_rows(old.dimension_inputs(edges), 5)
    plan.close_local(x_edge_limits={9: 4})
    progress('Regenerated %d local/dimension rows.' % len(plan.system.rows))
    generators = load_annihilators(HERE.parent / REFERENCE_FILES[1], edges)
    remainder = [(tuple(powers), Fraction(coefficient))
                 for powers, coefficient in previous['all_residual_scalar_periods']]
    annihilator = AnnihilatorTarget(edges, remainder, 7, generators)
    generated = dict(annihilator.generate(max_candidates=100000, include_symmetries=True))
    if generated['generation_stopped_by'] is not None:
        raise ArithmeticError('the bounded inherited annihilator generation stopped')
    progress('Regenerated %d complete annihilator rows.' % len(annihilator.plan.system.rows))
    representatives, indices, rows = [], {}, []
    for system in (plan.system, annihilator.plan.system):
        remap = []
        for _, powers in system.reps:
            if powers not in indices:
                indices[powers] = len(representatives)
                representatives.append(powers)
            remap.append(indices[powers])
        rows.extend(({remap[index]: Fraction(value) for index, value in row.items()}, _labels(finite))
                    for row, finite in system.rows)
    converted = []
    for target in targets:
        converted.append({indices[powers]: value for powers, value in target.items()})
    metadata = {'variables': len(representatives), 'rows': len(rows),
                'rows_sha256': row_digest(representatives, rows),
                'local_relation_counts': dict(plan.system.counts),
                'annihilator_generation': generated,
                'annihilator_polynomial_verification': dict(annihilator.verification),
                'construction_options': CONSTRUCTION,
                'finite_seed_search_states': 40, 'vertex_row_edge_power_limits': {'9': 4},
                'lower_dimension': 5, 'dimension_rows': 687,
                'all_coordinate_faces_retained': True, 'read_only_period_tables': True}
    return rows, converted, representatives, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    options = parser.parse_args()
    started = perf_counter()
    def progress(message):
        print('%s (%.3f s)' % (message, perf_counter() - started), flush=True)
    targets, source_metadata = reconstruct_source_coefficients(progress=lambda point, record:
        progress('Source %s: %d original / %d scalar terms.' %
                 (point, record['original_numerator']['nonzero_terms'],
                  record['scalar_numerator']['nonzero_terms'])))
    control = fresh_source_control(targets)
    progress('Fresh general-frame source: every one of %d scalar coefficients agrees.' %
             control['scalar_numerator']['nonzero_terms'])
    edges = tuple(source_metadata['scalar_edges'])
    rows, converted, representatives, relations = inherited_system(targets, edges, progress)
    previous = None if options.write else json.loads(CERTIFICATE.read_text())
    if options.write:
        results = eliminate_batch_with_witness(rows, converted, representatives, edges, 7, **ELIMINATION)
        progress('Shared elimination statistics: %s' % dict(results[0].statistics))
        if any(not result.complete for result in results):
            raise ArithmeticError('the shared source reduction retained at least one scalar period: %s' %
                                  [dict(result.statistics) for result in results])
        weights = [[[int(index), str(weight)] for index, weight in result.row_weights.items()]
                   for result in results]
        payload = json.dumps(weights, separators=(',', ':')).encode()
        WITNESS.write_bytes(gzip.compress(payload, mtime=0))
        finite = [{label: str(value) for label, value in result.value.items()} for result in results]
        statistics = dict(results[0].statistics)
        statistics.pop('witness_rows')
        statistics.pop('target_residual_periods')
    else:
        if hashlib.sha256(WITNESS.read_bytes()).hexdigest() != previous['witness_file_sha256']:
            raise ArithmeticError('the complete witness archive changed')
        weights = json.loads(gzip.decompress(WITNESS.read_bytes()))
        if len(weights) != len(targets):
            raise ArithmeticError('a source coefficient witness is missing')
        finite = previous['finite_coefficients_in_basis_order']
        statistics = previous['shared_elimination_statistics']
        for target, witness, values in zip(converted, weights, finite):
            verify_row_witness(rows, target, {}, {label: Fraction(value) for label, value in values.items()},
                               {int(index): Fraction(weight) for index, weight in witness})
    progress('All eleven original-row identities verified with zero residual.')
    # Independent old specialization: lambda_5=0, all earlier nonzero
    # sources retained.  In the centered basis every surviving covariant
    # contributes; this comparison is not one selected coefficient.
    old_value = {}
    for (i, _), values in zip(BASIS, finite):
        for label, value in values.items():
            old_value[label] = old_value.get(label, Fraction()) + (-1) ** i * Fraction(value)
    old_value = {label: value for label, value in old_value.items() if value}
    expected_old = {'1': Fraction(1187, 435456000), 'z3': Fraction(-816, 435456000)}
    if old_value != expected_old:
        raise ArithmeticError('the complete previously certified source specialization changed')
    data = {'schema_version': 1,
            'scope': 'full bare H=2,T=0 kernel; arbitrary balanced spinors at vertices 1,3,5,6; '
                     'sources at 0,2,4 vanish structurally; no current insertion',
            'source_reconstruction': source_metadata,
            'fresh_source_control': control, 'relations': relations,
            'finite_coefficients_in_basis_order': finite,
            'zero_residuals_for_every_coefficient': True,
            'shared_elimination_options': ELIMINATION,
            'shared_elimination_statistics': statistics,
            'witness_row_counts': [len(witness) for witness in weights],
            'witness_file': WITNESS.name,
            'witness_file_sha256': hashlib.sha256(WITNESS.read_bytes()).hexdigest(),
            'witness_uncompressed_sha256': hashlib.sha256(gzip.decompress(WITNESS.read_bytes())).hexdigest(),
            'old_complete_specialization': {label: str(value) for label, value in sorted(old_value.items())},
            'normalization': 'I_q=integral_Omega product(t_e^(alpha_e-1))/U^q; no Gamma(alpha) '
                             'or completion factor; all source/orientation signs retained',
            'reference_file_sha256': hashes(HERE.parent / name for name in REFERENCE_FILES),
            'period_table_sha256': hashes(HERE.parent / name for name in
                 ('cores/tables/prism_q7.json', 'cores/tables/K33_q7.json', 'k5/tables/K5_q7.json',
                  'fiveloop/tables/K33+e_q5.json')),
            'implementation_sha256': implementation_hashes()}
    data = json.loads(json.dumps(data))
    if options.write:
        CERTIFICATE.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
        progress('Wrote source_certificate.json and exact source_witnesses.json.gz.')
    elif data != previous:
        raise ArithmeticError('a deterministic source certificate field differs')
    else:
        progress('Every deterministic full-source certificate field replays.')


if __name__ == '__main__':
    main()
