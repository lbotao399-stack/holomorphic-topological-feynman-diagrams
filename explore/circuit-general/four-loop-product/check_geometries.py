"""Compare every generic-source graph moment with passive forest contraction."""
from collections import Counter
from fractions import Fraction
import json
from time import perf_counter

import quartic_b_squared as calculation
from htk import Graph, ProductForestMoments


def main():
    start = perf_counter()
    data = json.loads(calculation.CERTIFICATE.read_text())
    sample = data['independent_full_source_replays'][0]
    methods = Counter()
    checks = 0
    for row in sample['graphs']:
        graph = Graph(6, tuple(map(tuple, row['ordered_edges'])))
        sources = tuple(tuple(map(Fraction, vector)) for vector in row['full_balanced_sources'])
        with ProductForestMoments(graph, sources, marked=(4, 5), check=2, allow_series=False) as plan:
            methods[plan.method] += 1
            for edge, spinor in row['full_first_current_moments'].items():
                for coordinate, expected in enumerate(spinor):
                    actual = calculation.previous._periods(plan.moment({(int(edge), coordinate): 1}))
                    if actual != {label: Fraction(value) for label, value in expected.items()}:
                        raise ArithmeticError('an independently contracted original current disagrees')
                    checks += 1
    if methods['sourcefree_contractions'] == 0:
        raise ArithmeticError('no independent passive contraction was exercised')
    print('Every original graph current agrees:', checks, 'complete moments;', dict(methods))
    print('Elapsed %.3f seconds.' % (perf_counter() - start))


if __name__ == '__main__':
    main()
