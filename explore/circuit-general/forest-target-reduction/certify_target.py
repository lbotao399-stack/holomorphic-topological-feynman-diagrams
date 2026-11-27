"""Reconstruct and reduce every coefficient of the missing five-loop forest.

Default: replay the compact exact certificate.  --write: create it.  The
optional --dimension-shift additionally replays the bounded unsuccessful
q=5 -> q=7 lift; it has its own certificate and never changes the fast one.
No giant numerator, equation matrix, or intermediate row dump is stored.
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
from family_target import FamilyTarget, IncompleteFamilyTarget
import engine


EDGES = ((0, 1), (1, 2), (0, 3), (0, 4), (2, 5), (2, 6),
         (3, 5), (3, 6), (4, 5), (4, 6), (5, 6))
SOURCES = ((0, 0), (1, 0), (0, 0), (0, 1), (0, 0), (0, 0), (-1, -1))
MARKED = (5, 6)
SOURCE_OPTIONS = {'check': 2, 'max_terms': 1000000}
CONSTRUCTION_OPTIONS = {'max_input_terms': 10000, 'max_variables': 70000,
                        'max_rows': 250000, 'max_row_coefficients': 5000000}
ELIMINATION_OPTIONS = {'max_operations': 100000000,
                       'max_active_coefficients': 5000000, 'max_pivots': 70000}
INPUT_DIGESTS = {
    'original': '0598c5a835d1a0f2f6612f170bae9deb0d915a14a78c6644963d9a384ad1c50e',
    'scalar': '393711b3394c9642862a0bee38f24b39e8a02f6859d5ca443335a00065b9acc3',
}


def digest(rows):
    return hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode('utf-8')).hexdigest()


def numerator_record(terms):
    rows = [[list(map(int, powers)), str(Fraction(str(coefficient)))]
            for powers, coefficient in terms if coefficient]
    rows.sort()
    return {'nonzero_terms': len(rows), 'sha256': digest(rows)}


def dimension_inputs(edges):
    """All 687 q5 orbits whose U² shifts stay in the sixteen source blocks."""
    path = HERE.parent / 'fiveloop/tables/K33+e_q5.json'
    table = json.loads(path.read_text())
    old_edges = tuple(map(tuple, table['edges']))
    mapping = next(engine.nx.isomorphism.GraphMatcher(
        engine.nx.Graph(old_edges), engine.nx.Graph(edges)).isomorphisms_iter())
    edge_index = {frozenset(edge): index for index, edge in enumerate(edges)}
    permutation = [edge_index[frozenset((mapping[u], mapping[v]))] for u, v in old_edges]
    # This ordered graph is recorded explicitly in the certificate.  These
    # column sums are invariants of its diamond rule, not discarded sources.
    columns = ((0, 3, 4), (1, 5, 6), (2, 7, 8))
    expected = ((0, 1), (0, 2), (0, 3), (1, 4), (1, 5),
                (2, 4), (2, 5), (3, 4), (3, 5), (4, 5))
    if tuple(edges) != expected:
        raise ArithmeticError('the scalar edge order changed; map source-support invariants explicitly')
    automorphisms = engine.gutil.edge_automorphisms(edges)
    result = []
    for powers, value in table['values']:
        transformed = [0] * len(edges)
        for old, new in enumerate(permutation):
            transformed[new] = powers[old]
        if transformed[9] == 3 and all(sum(transformed[index] for index in column) <= 8
                                       for column in columns):
            result.append((engine.gutil.orbit_rep(tuple(transformed), automorphisms), value))
    result.sort()
    if len(result) != 687:
        raise ArithmeticError('the complete lower-level source-block selection changed')
    return result


def implementation_hashes():
    paths = {Path(__file__).resolve(), ROOT / 'code/htk/_forest_period_worker.py',
             ROOT / 'explore/nonSP/verify/k4_periods/current_minor.py'}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT):
                paths.add(path)
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def calculate(with_dimension=False):
    started = perf_counter()
    graph = Graph(7, EDGES, 'source-complete-K33-plus-edge-target')
    with ProductForestMoments(graph, SOURCES, MARKED, **SOURCE_OPTIONS) as forest:
        series = forest.moment()
        if not isinstance(series, ExactForestSeries) or series.reason != 'K33+e at q=7 not tabulated':
            raise ArithmeticError('the original finite-table status changed')
        if len(series.components) != 1 or series.constant:
            raise ArithmeticError('the original complete kernel changed scalar components')
        q, numerator = series.components[0]
        original = numerator_record(forest._core._bare.terms())
        scalar = numerator_record(numerator.items())
        if original['sha256'] != INPUT_DIGESTS['original'] or scalar['sha256'] != INPUT_DIGESTS['scalar']:
            raise ArithmeticError('the complete original numerator differs from the immutable input reference')
        edges = tuple(series.edges)
        terms = [(tuple(int(power) + 1 for power in exponent), coefficient)
                 for exponent, coefficient in sorted(numerator.items())]
    print('Reconstructed all %d original and %d scalar coefficients (%.3f s).' %
          (original['nonzero_terms'], scalar['nonzero_terms'], perf_counter() - started), flush=True)
    system = FamilyTarget(edges, terms, q, **CONSTRUCTION_OPTIONS)
    system.seed_finite(search_states=40)
    print('Finite input seeds: %d (%.3f s).' %
          (system.system.counts['finite_seed'], perf_counter() - started), flush=True)
    system.close_local(x_edge_limits={9: 4})
    before_dimension = {'variables': len(system.system.reps), 'rows': len(system.system.rows)}
    if with_dimension:
        system.add_dimension_rows(dimension_inputs(edges), 5)
        system.close_local(x_edge_limits={9: 4})
    print('Exact local system: %d variables and %d rows (%.3f s).' %
          (len(system.system.reps), len(system.system.rows), perf_counter() - started), flush=True)
    result = system.reduce(**ELIMINATION_OPTIONS)
    print('Whole target reduced to %d residual periods (%.3f s).' %
          (len(result.residuals), perf_counter() - started), flush=True)
    rejected = False
    try:
        result.value
    except IncompleteFamilyTarget:
        rejected = True
    residuals = sorted([[list(term.powers), str(term.coefficient)] for term in result.residuals])
    data = {
        'schema_version': 1,
        'scope': 'H=2,T=0; full listed balanced sources; bare current word',
        'original_graph': {'vertices': 7, 'loops': graph.loops, 'edges': EDGES},
        'full_balanced_sources': SOURCES, 'marked': MARKED,
        'source_options': SOURCE_OPTIONS,
        'original_numerator': original,
        'scalar_core': {'edges': edges, 'q': q, **scalar},
        'construction_options': CONSTRUCTION_OPTIONS,
        'finite_seed_search_states_per_completed_orbit': 40,
        'vertex_identity_edge_power_limits': {'9': 4},
        'elimination_options': ELIMINATION_OPTIONS,
        'dimension_shift': with_dimension,
        'local_system_before_dimension_shift': before_dimension,
        'finite_coefficients': {label: str(value) for label, value in sorted(result.finite_coefficients.items())},
        'all_residual_scalar_periods': residuals,
        'residual_sha256': digest(residuals),
        'statistics': dict(result.statistics),
        'complete_finite_reduction': result.complete,
        'incomplete_finite_value_rejected': rejected,
        'kernel_nonvanishing_proved': False,
        'identity': ('original full sourced kernel = finite_coefficients + sum of EVERY '
                     '[powers, coefficient] row in all_residual_scalar_periods, with '
                     'I_q = integral_Omega product(t_e^(powers_e-1))/U_G^q on scalar_core.edges. '
                     'No Gamma(powers) factor and no completed-Pi normalization.'),
        'meaning': ('This is the entire kernel, including every unresolved coefficient. '
                    'The finite contribution alone is not the value. No zero/nonzero '
                    'conclusion for the original kernel follows from this decomposition.'),
        'input_digest_encoding': ('SHA256 of UTF-8 compact JSON of sorted '
                                  '[[full exponent list, reduced Fraction string], ...].'),
    }
    data['implementation_sha256'] = implementation_hashes()
    tables = [HERE.parent / 'cores/tables/prism_q7.json',
              HERE.parent / 'cores/tables/K33_q7.json',
              HERE.parent / 'k5/tables/K5_q7.json']
    if with_dimension:
        tables.append(HERE.parent / 'fiveloop/tables/K33+e_q5.json')
    data['period_table_sha256'] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                   for path in tables}
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--dimension-shift', action='store_true')
    options = parser.parse_args()
    started = perf_counter()
    result = json.loads(json.dumps(calculate(options.dimension_shift)))
    name = 'dimension_certificate.json' if options.dimension_shift else 'target_certificate.json'
    path = HERE / name
    if options.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
        print('Wrote %s.' % name, flush=True)
    else:
        if result != json.loads(path.read_text()):
            raise ArithmeticError('a deterministic whole-target certificate field differs')
        print('Every exact certificate field replays.', flush=True)
    print('Finite contribution: %s.' % result['finite_coefficients'], flush=True)
    print('All residuals retained: %d; complete finite result: %s.' %
          (len(result['all_residual_scalar_periods']), result['complete_finite_reduction']), flush=True)
    print('Elapsed %.3f seconds.' % (perf_counter() - started), flush=True)


if __name__ == '__main__':
    main()
