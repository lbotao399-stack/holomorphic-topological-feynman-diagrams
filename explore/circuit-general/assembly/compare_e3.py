"""Compare assembled bare kernels (results/*.jsonl) with atlas/e3_*.log (fastker, same seeds).
Standard-source records (bulk3 --std): K_e3 = c det(lambda_1, lambda_2, lambda_3)^4 with the e3 sources.
usage: python3 compare_e3.py results/FILE.jsonl [...]"""
import ast
import glob
import json
import os
import re
import sys
from fractions import Fraction as F

import sets

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = os.path.join(HERE, '..', 'atlas')
PAT = re.compile(r'^(\S+) #(\d+) SP=(\w+) fastker K = (\{.*?\}) \(')


def e3_values():
    out = {}
    for f in (glob.glob(os.path.join(ATLAS, 'e3_bulk3_s*.log')) + glob.glob(os.path.join(ATLAS, 'e3_twistor4_s*.log'))
              + [os.path.join(ATLAS, 'e3_twistor3.log')]):
        for line in open(f):
            m = PAT.match(line)
            if not m:
                continue
            name = 'bulk3' if m.group(1).startswith('bulk3') else m.group(1)
            out[(name, int(m.group(2)))] = {k: F(v) for k, v in ast.literal_eval(m.group(4)).items()}
    return out


def det3(L):
    return (L[0][0] * (L[1][1] * L[2][2] - L[1][2] * L[2][1]) - L[0][1] * (L[1][0] * L[2][2] - L[1][2] * L[2][0])
            + L[0][2] * (L[1][0] * L[2][1] - L[1][1] * L[2][0]))


def main(files):
    ref = e3_values()
    its = {name: dict(sets.items(name)) for name in ('bulk3',)}
    n = match = 0
    degenerate = []
    for fn in files:
        for line in open(fn):
            r = json.loads(line)
            key = (r['set'], r['idx'])
            if key not in ref or 'bare' not in r:
                continue
            val = {k: F(v) for k, v in r['bare'].items()}
            if 'legs' in r:
                G, O, S0 = its[r['set']][r['idx']]
                Gh, lab, lam, full = sets.instance(r['set'], r['idx'], G, O, S0, 3)
                d4 = F(det3([full[v] for v in r['legs']])) ** 4
                val = {k: d4 * v for k, v in val.items() if d4 * v}
                if d4 == 0:
                    degenerate.append(r['idx'])
            n += 1
            ok = val == ref[key]
            match += ok
            if not ok:
                print('MISMATCH', key, val, ref[key])
    print('compared %d bare kernels with e3 logs: %d equal; e3 sources with det(lambda_1,lambda_2,lambda_3) = 0: %s' % (n, match, degenerate))


if __name__ == '__main__':
    main(sys.argv[1:])
