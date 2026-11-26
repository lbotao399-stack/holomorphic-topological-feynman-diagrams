"""Closure of the regulator-free relation systems (all convergent orbits determined) for the
prism and K33 at q = 3,4,5 (and 7/2), and the prism system at q = 6 against the 1481-entry
FORCER bank of explore/nonSP/verify/prism_forcer/periods.json."""
import json
import os
from fractions import Fraction as F

import cores
import gutil
import msys

for name in ('prism', 'K33'):
    for q in (3, 4, 5):
        S = msys.MultiSystem([cores.FAMILIES[name]], q)
        S.build(verbose=False)
        solved, piv = S.solve(verbose=False)
        print('%-5s q=%d: convergent orbits %d, rows %d %s, determined %d' % (name, q, len(S.reps), len(S.rows), S.counts, len(solved)), flush=True)

S = msys.MultiSystem([cores.PRISM], 6)
S.build(verbose=False)
solved, piv = S.solve(verbose=False)
print('prism q=6: convergent orbits %d, determined %d' % (len(S.reps), len(solved)), flush=True)
bank = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../nonSP/verify/prism_forcer/periods.json')
data = json.load(open(bank))
BANK_PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]
ok = bad = 0
for ex, co in zip(data['exponents'], data['coefficients']):
    pw = [k + 1 for k in ex]
    val = cores.period(BANK_PRISM, pw, 6)
    ref = {k: F(v) for k, v in zip(('1', 'z3', 'z5'), co) if F(v)}
    if val == ref:
        ok += 1
    else:
        bad += 1
print('FORCER bank (q=6, 1481 monomials): agree %d, disagree %d' % (ok, bad), flush=True)
