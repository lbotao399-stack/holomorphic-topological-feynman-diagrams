"""High-level samples: prism at q = 6, 7 against cores/tables; K33 at q = 6 (completion) against cores/tables;
Z5 at q = 5 against the fiveloop relation table if present; one Z5 orbit at q = 7 (timing)."""
import json
import os
import sys
import time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'cores'))
import seq
import mzv
from validate import to_labels
from coverage import cores_list

E = dict(cores_list())


def run(name, q, a, ref=None):
    t = time.time()
    try:
        v, pl = seq.period(list(E[name]), list(a), q)
    except Exception as e:
        print('%s q=%d %s: %s %s' % (name, q, a, type(e).__name__, e), flush=True)
        return
    tag = ''
    if ref is not None:
        tag = ' ref:%s' % ('OK' if to_labels(v) == ref else 'MISMATCH %s' % ref)
    print('%s q=%d %s: %s  (%.1fs, %s, %d appends)%s' % (name, q, list(a), mzv.tostr(v), time.time() - t,
                                                          pl['route'], pl['appends'], tag), flush=True)


for name, q in (('prism', 6), ('prism', 7), ('K33', 6)):
    tab = json.load(open(os.path.join(HERE, '..', 'cores', 'tables', '%s_q%d.json' % (name, q))))['values']
    picks = [tab[0], tab[len(tab) // 2], tab[-1]]
    for r, val in picks:
        c = seq.plans_for(list(E[name]), list(r), q)
        if not c:
            print('%s q=%d %s: no plan' % (name, q, r))
            continue
        run(name, q, r, {k: F(x) for k, x in val.items()})
run('Z5', 7, (4, 3, 3, 4, 3, 4, 3, 3, 4, 4))
