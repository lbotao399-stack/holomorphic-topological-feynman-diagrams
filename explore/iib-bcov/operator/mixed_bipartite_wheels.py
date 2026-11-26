"""Four-wheel cuts of the q=13 one-star mixed Ward chamber."""
from argparse import ArgumentParser
from fractions import Fraction as F
import gzip
import json
from math import factorial, prod
from pathlib import Path
from subprocess import run
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from quartic_prism_completion import completed, prism_charts
from quartic_prism_physical_ibp import compile_program, zeta_interval
from quartic_wheel_exact import fields, field_table, reconstruct

TOTAL = (F(612680153354189825651120954254291, 3629389754900829974326272000000000),
         F(247723220447053749700079347, 548660582751448219852800000),
         -F(2783935140019231, 4056996571992000))


def evaluate(complete, wheel, scratch, q, inputs):
    output = run([complete, str(q), '--wheel-charts'], input=inputs, text=True,
                 check=True, capture_output=True)
    cuts = {}
    for line in output.stdout.splitlines():
        row = line.split()
        assert len(row) == 13
        seed = tuple(map(int, row[:9]))
        assert seed not in cuts
        cuts[seed] = (F(row[9]), int(row[11]), F(row[12]))
    keys = scratch/'keys'
    keys.write_text(''.join(f'{k}\n' for k in sorted({k for _, k, _ in cuts.values()})))
    primes = fields(q)
    tables = [field_table(wheel, keys, scratch, None, q, p) for p in primes]
    values_file = scratch/'wheel.tsv'
    reconstruct(q, primes, tables, keys, values_file)
    values = {}
    for line in values_file.read_text().splitlines():
        k, *columns = line.split()
        values[int(k)] = tuple(map(F, columns))
    return cuts, values


def lower_dimension(complete, wheel, scratch):
    reference = json.loads((HERE.parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    known = {}
    inputs = []
    for ex, coeffs in zip(reference['exponents'], reference['coefficients'], strict=True):
        seed = tuple(a+1 for a in ex)
        key = min(prism_charts(completed(seed, 6), 6))
        norm = F(factorial(5), prod(map(factorial, ex)))
        expected = tuple(norm*F(str(c)) for c in coeffs)
        assert known.setdefault(key, expected) == expected
        inputs.append(' '.join(map(str, (*seed, 1, 0))))
    cuts, values = evaluate(complete, wheel, scratch, 6, '\n'.join(inputs))
    for key, (_, wheel_key, scale) in cuts.items():
        assert tuple(scale*x for x in values[wheel_key]) == known[key]
    print('All three independent q=6 coefficients verified on wheel components:', len(cuts))


def partial_fractions(values):
    from quartic_bipartite_wheel import period
    # Choose the least expensive keys by the number of Gegenbauer shifts.
    def powers(key):
        a = [(key >> (4*i)) & 15 for i in reversed(range(9))]
        return (a[4], a[7], a[8], a[5]), (a[3], a[1], a[6], a[2])
    for key in sorted(values, key=lambda k: (prod(13-x for x in powers(k)[1]), k))[:12]:
        spokes, ring = powers(key)
        integrated = period(spokes, ring, 13)
        assert tuple(F(str(integrated.get(k, 0))) for k in (0, 3, 5)) == values[key]
    print('All three q=13 coefficients verified by direct rational partial fractions on 12 keys.')


def physical(complete, wheel, scratch, check):
    with gzip.open(HERE/'mixed_bipartite_one_star.components.gz', 'rt') as table:
        rows = [r[:11] for line in table if (r := line.split())[11] == '0']
    assert len(rows) == 35563
    cuts, values = evaluate(complete, wheel, scratch, 13,
                            ''.join(' '.join(r)+'\n' for r in rows))
    assert len(cuts) == len(values) == 10940
    if check:
        partial_fractions(values)
    sums = [F(), F(), F()]
    for weight, key, scale in cuts.values():
        for i, value in enumerate(values[key]):
            sums[i] += weight*scale*value
    assert tuple(sums) == TOTAL
    lo = hi = sums[0]
    for weight, power in zip(sums[1:], (3, 5)):
        low, high = sorted(weight*x for x in zeta_interval(power, terms=96))
        lo += low
        hi += high
    print('10940 wheel components:', *sums)
    print('Open components:', len(rows)-len(cuts))
    print('10^20 times the subtotal lies between', (10**20*lo)//1, 'and', -((-10**20*hi)//1))


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Compare FORCER and direct partial fractions')
    args = parser.parse_args()
    with TemporaryDirectory() as work:
        scratch = Path(work)
        complete, wheel = scratch/'complete', scratch/'wheel'
        compile_program('quartic_prism_completion.cpp', complete)
        compile_program('quartic_wheel_bank.cpp', wheel)
        if args.check:
            lower_dimension(complete, wheel, scratch)
        physical(complete, wheel, scratch, args.check)
