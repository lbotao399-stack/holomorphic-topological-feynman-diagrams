"""Optional controlled default-versus-planned exact numerator benchmark.

The optimized run is always performed.  --baseline additionally reruns the
unchanged default source numerator and compares every original polynomial
and every scalar-core coefficient; that baseline can take minutes.
"""
import argparse
from fractions import Fraction
import json
from time import perf_counter

from source_family import HERE, EDGES, sources_at, polynomial_record
from htk import Graph
from htk.forest_periods import CoreCurrentPeriods
from tree_minor import prepare_source_numerator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', action='store_true')
    options = parser.parse_args()
    sources = sources_at(3, 1)
    planned = CoreCurrentPeriods(Graph(7, EDGES), sources, check=0, max_terms=3000000)
    baseline = None
    try:
        started = perf_counter()
        plan = prepare_source_numerator(planned, tree_candidates=((1, 2, 5, 7, 9, 10),),
                                        max_tree_subsets=1, max_plans=1)
        numerator_seconds = perf_counter() - started
        scalar = planned._core_polynomials({(0,) * 22: Fraction(1)})[7]
        print('Planned numerator %.3f s; full scalar aggregation %.3f s; %d / %d terms.' %
              (numerator_seconds, perf_counter() - started, len(planned._bare), len(scalar)))
        certificate = json.loads((HERE / 'source_certificate.json').read_text())
        reference = next(row for row in certificate['source_reconstruction']['samples']
                         if row['interpolation_point'] == [3, 1])
        assert polynomial_record(planned._bare.terms()) == reference['original_numerator']
        assert polynomial_record(scalar.items()) == reference['scalar_numerator']
        if options.baseline:
            baseline = CoreCurrentPeriods(Graph(7, EDGES), sources, check=0, max_terms=3000000)
            started = perf_counter()
            baseline._source_numerator()
            print('Default numerator %.3f s.' % (perf_counter() - started))
            assert baseline._bare == planned._bare
            assert baseline._core_polynomials({(0,) * 22: Fraction(1)})[7] == scalar
            print('Every original and scalar-core coefficient agrees exactly.')
        else:
            print('Both full coefficient digests equal the source certificate.')
    finally:
        planned.close()
        if baseline is not None:
            baseline.close()


if __name__ == '__main__':
    main()
