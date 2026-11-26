"""Pair descendant graph weights with the three certified positive K3,3 period tables."""
from argparse import ArgumentParser
from fractions import Fraction
import json
from pathlib import Path
from subprocess import run
import sys

IIB = Path(__file__).resolve().parents[1] / 'iib-bcov'
sys.path.insert(0, str(IIB))
from quartic_scalar_response import BIPARTITE


def contract(native, tables, banks):
    values = {p: [] for p in banks}
    for i, expected in enumerate(BIPARTITE):
        result = run([native, tables / f'channel_{i}.tsv', tables / 'anchor.tsv',
                      IIB / 'quartic_k33.weights.gz', str(i), *banks],
                     capture_output=True, text=True, check=True)
        control, *rows = result.stdout.splitlines()
        assert tuple(map(Fraction, control.split()[-2:])) == expected
        assert len(rows) == len(banks)
        for p, line in zip(banks, rows):
            tag, source, count, first, second = line.split('\t')
            assert tag == 'bank' and p == Path(source) and Fraction(second) == 0
            if values[p] and int(count) != values[p][0][0]:
                raise ValueError(f'Inconsistent number of moments for {p}')
            values[p].append((int(count), str(Fraction(first))))
    return {p.name: {'orbits': rows[0][0], 'coefficients': [c for _,c in rows]}
            for p, rows in values.items()}


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('native', type=Path, help='compiled quartic_bipartite_totals.cpp')
    parser.add_argument('tables', type=Path, help='output of quartic_bipartite_replay.py')
    parser.add_argument('banks', type=Path, nargs='+')
    args = parser.parse_args()
    print(json.dumps(contract(args.native.resolve(), args.tables, args.banks), indent=2))
