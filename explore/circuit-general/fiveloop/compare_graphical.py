"""Exact comparison of the five-loop tables (tables/) with the independent graphical-function tables
(../graphical/tables/, graphical functions in d = 2q, edge order of skel5.json)."""
import glob
import json
import os
from fractions import Fraction as F

import per5
import gutil

HERE = os.path.dirname(os.path.abspath(__file__))
for path in sorted(glob.glob(os.path.join(HERE, '..', 'graphical', 'tables', '*.json'))):
    d = json.load(open(path))
    name, q = d['core'], int(d['q'])
    if name not in per5.SKEL:
        continue
    E = per5.SKEL[name]
    if [tuple(e) for e in d['edges']] != list(E):
        print(name, q, 'edge order differs; skipped')
        continue
    autos = gutil.edge_automorphisms(tuple(E))
    ok = bad = miss = 0
    for a, val in d['values']:
        v = {k: F(x) for k, x in val.items() if F(x)}
        try:
            mine = per5.period(list(E), list(a), q)
        except Exception:
            miss += 1
            continue
        if mine == v:
            ok += 1
        else:
            bad += 1
    print('%-6s q=%d: graphical %d entries; agree %d, differ %d, not available here %d' % (name, q, len(d['values']), ok, bad, miss), flush=True)
