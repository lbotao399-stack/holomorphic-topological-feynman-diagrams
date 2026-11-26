"""Compare the triangle-rule prism periods with the FORCER bank of the nonSP run
(1481 monomials t^k / U^6, prism edges 01,02,03,14,15,23,25,34,45)."""
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

import prop
from laurent import as_zeta_dict

BANK = Path(__file__).resolve().parents[2] / 'nonSP/verify/prism_forcer/periods.json'
PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]


def main(limit=None, cuts=(0, 1), stride=1):
    data = json.loads(BANK.read_text())
    rows = list(zip(data['exponents'], data['coefficients']))[::stride]
    if limit:
        rows = rows[:limit]
    t0 = time.time()
    bad = 0
    for i, (ex, co) in enumerate(rows):
        powers = [k + 1 for k in ex]
        ref = {'1': F(co[0]), 'z3': F(co[1]), 'z5': F(co[2])}
        ref = {k: v for k, v in ref.items() if v}
        for cut in cuts:
            got = as_zeta_dict(prop.cut_period(PRISM, powers, 6, cut=cut))
            if got != ref:
                bad += 1
                print('MISMATCH', ex, 'cut', cut, got, ref, flush=True)
    print('checked %d monomials x %d cuts, mismatches %d, %.1fs' % (len(rows), len(cuts), bad, time.time() - t0))


if __name__ == '__main__':
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else None
    stride = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    main(lim, stride=stride)
