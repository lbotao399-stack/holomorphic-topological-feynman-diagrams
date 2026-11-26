"""Compare the tower K5 tables with K5 values computed elsewhere in the project (independent methods).
Sources: graphical/REPORT.md 4/6/7 (graphical functions), symmetry/REPORT.md 3.3 (engine), fiveloop/REPORT.md 5
(Y -> K33+e)."""
import json
import os
from fractions import Fraction as F

from k5tables import canon5

HERE = os.path.dirname(os.path.abspath(__file__))
REF = [
    # (q, alpha, value, source)
    (3, (1, 1, 2, 2, 2, 2, 2, 2, 2, 2), {'z3': F(5, 3), 'z5': F(-5, 3)}, 'graphical 6 (6 plans); fiveloop 5'),
    (3, (1, 2, 2, 2, 2, 2, 2, 1, 2, 2), {'z3': F(-2), 'z5': F(5, 2)}, 'fiveloop 5 (Y->K33+e)'),
    (7, (4, 2, 5, 5, 5, 6, 5, 2, 5, 3), {'1': F(3251321, 64665216000), 'z3': F(52643, 449064000), 'z5': F(-23, 124740)},
     'graphical 7 / symmetry 3.3'),
    (8, (5, 5, 6, 1, 6, 4, 5, 6, 5, 5), {'1': F(105943, 1569209241600), 'z3': F(-17, 302702400)}, 'graphical 7 / symmetry 3.3'),
    (8, (4, 4, 6, 7, 6, 1, 5, 5, 5, 5), {'1': F(-62311, 224172748800), 'z3': F(1, 4324320)}, 'graphical 7 / symmetry 3.3'),
    (8, (5, 6, 4, 1, 6, 5, 6, 7, 2, 6), {'1': F(263371, 7472424960000), 'z3': F(-1, 34594560)}, 'graphical 7'),
    (8, (3, 3, 7, 5, 5, 4, 7, 3, 5, 6), {'1': F(-8711086259, 3923023104000), 'z3': F(-26946407, 4540536000),
                                         'z5': F(5689, 630630)}, 'graphical 7'),
    (8, (2, 5, 7, 6, 6, 5, 4, 7, 4, 2), {'1': F(378433667, 4707627724800), 'z3': F(21719, 100900800),
                                         'z5': F(-5, 15288)}, 'graphical 7'),
    (8, (5, 5, 7, 1, 7, 2, 7, 6, 4, 4), {'1': F(-735555239, 201180672000), 'z3': F(-17797783, 1816214400),
                                         'z5': F(3755, 252252)}, 'graphical 7'),
]
tabs = {}
for q, a, v, src in REF:
    if q not in tabs:
        p = os.path.join(HERE, 'tables', 'K5_q%d.json' % q)
        tabs[q] = {tuple(r[0]): (r[1], {k: F(x) for k, x in r[2].items()}) for r in json.load(open(p))['values']} \
            if os.path.exists(p) else None
    t = tabs[q]
    if t is None:
        print('q=%d %s: no tower table' % (q, a))
        continue
    kd, val = t[canon5(a)]
    print('q=%d %s [%s]: tower %s ; reference (%s) %s' % (q, a, kd, 'EQUAL' if val == v else 'DIFFERENT', src, v))
