"""Replay the full missing five-loop forest through bounded conformal reduction.

The default reconstructs every original/source-complete coefficient, reduces
the whole scalar combination, and compares every deterministic certificate
field.  --write refreshes forest_combination_certificate.json.  The compact
record hashes the complete sorted numerator and residual expression; neither
is omitted from the computation or replaced by a numerical sample.

This is an explicitly incomplete finite reduction.  The known finite part
plus every recorded-by-digest residual is the whole expression.  The existing
public forest backend still supplies its complete exact convergent series.
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

from htk import ExactForestSeries, Graph, ProductForestMoments
from combination import IncompleteCombination, reduce_period_combination


EDGES = ((0, 1), (1, 2), (0, 3), (0, 4), (2, 5), (2, 6),
         (3, 5), (3, 6), (4, 5), (4, 6), (5, 6))
SOURCES = ((0, 0), (1, 0), (0, 0), (0, 1), (0, 0), (0, 0), (-1, -1))
MARKED = (5, 6)
SOURCE_OPTIONS = {'check': 2, 'max_terms': 1000000}
REDUCTION_OPTIONS = {'max_states': 30000, 'max_nodes': 50000,
                     'max_input_terms': 10000, 'max_edges': 10,
                     'max_vertices': 8, 'use_inverse_uniqueness': True}


def _digest(rows):
    encoded = json.dumps(rows, separators=(',', ':'), ensure_ascii=True).encode('utf-8')
    return hashlib.sha256(encoded).hexdigest()


def numerator_digest(terms):
    rows = sorted((tuple(map(int, exponent)), str(Fraction(str(coefficient))))
                  for exponent, coefficient in terms if coefficient)
    return {'nonzero_terms': len(rows),
            'sha256': _digest([[list(exponent), coefficient] for exponent, coefficient in rows])}


def residual_digest(terms):
    rows = []
    for term in terms:
        weighted = []
        for first, second, weight in term.weights:
            if Fraction(weight).denominator != 1:
                raise ArithmeticError('an integer uniqueness move produced a fractional edge weight')
            weighted.append([int(first), int(second), int(weight)])
        rows.append([term.q, term.vertices, weighted, str(term.coefficient)])
    rows.sort()
    return {'nonzero_terms': len(rows), 'sha256': _digest(rows)}


def implementation_hashes():
    # Include the actual repository modules loaded in this fresh CLI process.
    # The determinant is loaded under an isolated name without sys.modules,
    # and the finite-period worker runs in its own process, so add both files.
    paths = {HERE / 'certify_forest_combination.py',
             ROOT / 'code/htk/_forest_period_worker.py',
             ROOT / 'explore/nonSP/verify/k4_periods/current_minor.py'}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if not filename:
            continue
        path = Path(filename).resolve()
        if path.suffix == '.py' and path.is_relative_to(ROOT):
            paths.add(path)
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def calculate():
    start = perf_counter()
    graph = Graph(7, EDGES, 'missing-K33-plus-edge-five-loop-product-forest')
    with ProductForestMoments(graph, SOURCES, MARKED, **SOURCE_OPTIONS) as forest:
        series = forest.moment()
        if not isinstance(series, ExactForestSeries):
            raise ArithmeticError('the missing-period fixture no longer selects the exact series')
        if series.reason != 'K33+e at q=7 not tabulated':
            raise ArithmeticError('the finite-period availability of this fixture changed')
        if forest.last_stats['contraction']['contraction_steps'] != 0:
            raise ArithmeticError('the irreducible sourced fixture acquired a passive contraction')
        if len(series.components) != 1 or series.constant:
            raise ArithmeticError('the bare fixture no longer has one scalar denominator power')
        q, numerator = series.components[0]
        if q != 7:
            raise ArithmeticError('the expected scalar denominator power changed')
        bare = numerator_digest(forest._core._bare.terms())
        core = numerator_digest(numerator.items())
        terms = [(series.edges, tuple(exponent + 1 for exponent in powers), coefficient)
                 for powers, coefficient in sorted(numerator.items())]
        print('Reconstructed every source-complete coefficient: %d original, %d core terms (%.3f s).'
              % (bare['nonzero_terms'], core['nonzero_terms'], perf_counter() - start), flush=True)
        reduction_start = perf_counter()
        result = reduce_period_combination(terms, q, **REDUCTION_OPTIONS)
        print('Reduced the complete combination: %d residuals, %d expanded states (%.3f s).'
              % (len(result.residuals), result.statistics['expanded_states'],
                 perf_counter() - reduction_start), flush=True)
        if result.complete:
            raise ArithmeticError('this recorded incomplete finite reduction has become complete')
        try:
            _ = result.value
        except IncompleteCombination as error:
            incomplete = {'rejected': True, 'exception': type(error).__name__, 'message': str(error)}
        else:
            raise ArithmeticError('an incomplete reduction exposed its finite contribution as the value')
        data = {
            'schema_version': 1,
            'scope': 'H=2,T=0, listed full balanced sources, bare current word',
            'convention': 'htk kernel with the listed original edge order and orientations',
            'original_graph': {'vertices': graph.nv, 'loops': graph.loops, 'edges': graph.edges},
            'full_balanced_sources': SOURCES,
            'marked': MARKED,
            'current_powers': [],
            'source_reconstruction_options': SOURCE_OPTIONS,
            'original_bare_numerator': {'denominator_power': graph.loops + 2, **bare},
            'scalar_core': {'edges': series.edges, 'denominator_power': q, **core},
            'original_finite_period_unavailable_reason': series.reason,
            'reduction_options': REDUCTION_OPTIONS,
            'finite_coefficients': {label: str(value)
                                    for label, value in sorted(result.finite_coefficients.items())},
            'complete_residual_expression': residual_digest(result.residuals),
            'statistics': dict(result.statistics),
            'complete_finite_reduction': result.complete,
            'finite_value_access_control': incomplete,
            'kernel_nonvanishing_proved': False,
            'meaning': ('The finite contribution plus EVERY residual completed graph period is '
                        'the entire original kernel. The finite contribution alone is not the '
                        'answer. Unresolved coefficients are retained exactly; no zero or '
                        'nonzero conclusion is drawn for the kernel.'),
            'numerator_digest_encoding': (
                'UTF-8 compact JSON of [[full exponent tuple as list, reduced Fraction string], ...], '
                'sorted lexicographically; all nonzero coefficients included. Bare variables '
                'follow original edges; scalar variables follow scalar_core.edges.'),
            'residual_digest_encoding': (
                'UTF-8 compact JSON of sorted [q, vertex count, [[first,second,integer weight], ...], '
                'reduced Fraction coefficient string] rows. Every residual is coefficient times '
                'the complete weighted graph period Pi in symmetry/engine.py normalization. '
                'The graph and every signed edge weight are included, not a graph hash alone.'),
        }
    table_paths = (
        ROOT / 'explore/circuit-general/cores/tables/prism_q7.json',
        ROOT / 'explore/circuit-general/cores/tables/K33_q7.json',
        ROOT / 'explore/circuit-general/k5/tables/K5_q7.json',
    )
    data['implementation_sha256'] = implementation_hashes()
    data['period_table_sha256'] = {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in table_paths}
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='refresh the deterministic certificate')
    arguments = parser.parse_args()
    start = perf_counter()
    result = json.loads(json.dumps(calculate()))
    path = HERE / 'forest_combination_certificate.json'
    if arguments.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        print('Wrote the complete-expression reduction certificate.', flush=True)
    else:
        expected = json.loads(path.read_text(encoding='utf-8'))
        if result != expected:
            raise ArithmeticError('a deterministic whole-forest reduction certificate field differs')
        print('Every exact whole-forest reduction certificate field replays.', flush=True)
    print('Finite contribution: %s.' % result['finite_coefficients'], flush=True)
    print('Every residual remains: %d; finite value access is rejected.'
          % result['complete_residual_expression']['nonzero_terms'], flush=True)
    print('Elapsed %.3f seconds.' % (perf_counter() - start), flush=True)


if __name__ == '__main__':
    main()
