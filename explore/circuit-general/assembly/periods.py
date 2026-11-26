"""Core periods I_q(G; alpha) = int prod t^(alpha-1) / U^q (projective), exact.

Adapted copy of cores/cores.py (table / period): tables of prism and K33 orbits are read from
cores/tables (q <= 7) or assembly/tables (built here, q >= 8) and built on demand with the
relation system cores/msys.py (imported unchanged).  Wheels W_n by cores/wheel.py, any q;
series-parallel pieces by msys.MultiSystem.classify.

API
  table(name, q)                      {orbit rep: {label: Fraction}}
  period(edges, powers, q)            {label: Fraction}, labels '1', 'z3', 'z5', ...
  CoreEvaluator(edges, q)             cached, vectorizable lookups for one core graph
"""
import json
import os
import sys
import time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CORES = os.path.abspath(os.path.join(HERE, '..', 'cores'))
if CORES not in sys.path:
    sys.path.insert(0, CORES)

import gutil   # noqa: E402
import msys    # noqa: E402

TABLES_OUT = os.path.join(HERE, 'tables')
TABLES_IN = [TABLES_OUT, os.path.join(CORES, 'tables')]
PRISM = ((0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5))
K33 = ((0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5))
FAMILIES = {'prism': PRISM, 'K33': K33}

_systems = {}
_tables = {}


def system(q):
    if q not in _systems:
        _systems[q] = msys.MultiSystem([PRISM, K33], q)
    return _systems[q]


def table(name, q, verbose=True):
    q = int(q)
    key = (name, q)
    if key in _tables:
        return _tables[key]
    for d in TABLES_IN:
        path = os.path.join(d, '%s_q%d.json' % (name, q))
        if os.path.exists(path):
            raw = json.load(open(path))
            tab = {tuple(r): {k: F(v) for k, v in val.items()} for r, val in raw['values']}
            _tables[key] = tab
            return tab
    t = time.time()
    S = msys.MultiSystem([FAMILIES[name]], q)
    S.build(verbose=verbose)
    solved, piv = S.solve(verbose=verbose)
    if len(solved) != len(S.reps):
        und = [S.reps[j] for j in range(len(S.reps)) if j not in solved]
        raise ArithmeticError('%s q=%d: %d undetermined orbits, e.g. %s' % (name, q, len(und), und[:3]))
    tab = {S.reps[j][1]: v for j, v in solved.items()}
    os.makedirs(TABLES_OUT, exist_ok=True)
    path = os.path.join(TABLES_OUT, '%s_q%d.json' % (name, q))
    json.dump({'core': name, 'edges': FAMILIES[name], 'q': str(q), 'seconds': round(time.time() - t, 1),
               'values': [[list(r), {k: str(x) for k, x in v.items()}] for r, v in sorted(tab.items())]},
              open(path, 'w'))
    _tables[key] = tab
    return tab


def period(edges, powers, q):
    """exact period; raises ValueError if divergent / unrecognized"""
    S = system(q)
    res = S.classify(list(edges), list(powers), q)
    if res is None:
        raise ValueError('divergent or unrecognized core for powers %s' % (powers,))
    const, unk = res
    out = dict(const)
    for j, s in unk:
        f, rep = S.reps[j]
        val = table('prism' if f == 0 else 'K33', q)[tuple(rep)]
        for k, v in val.items():
            out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


class CoreEvaluator:
    """Period lookups for one core graph at level q, memoized per power vector."""

    def __init__(self, edges, q):
        self.edges = tuple(tuple(e) for e in edges)
        self.q = q
        self.cache = {}
        self.ncalls = 0

    def __call__(self, powers):
        key = tuple(int(x) for x in powers)
        v = self.cache.get(key)
        if v is None:
            self.ncalls += 1
            v = period(self.edges, key, self.q)
            self.cache[key] = v
        return v
