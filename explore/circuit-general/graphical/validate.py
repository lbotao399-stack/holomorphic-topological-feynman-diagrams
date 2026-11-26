"""Graphical-function periods vs independent project values (cores.period: wheel.py radial algorithm,
msys relation tables) for all convergent weight orbits.   usage: python3 validate.py GRAPH Q [max_orbits]"""
import os
import sys
import time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'cores'))
import gutil
import cores
import seq
import mzv

PRISM = ((0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5))


def wheel(n):
    return tuple([(0, i) for i in range(1, n + 1)] + [(i, i % n + 1) for i in range(1, n + 1)])


K33 = ((0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5))
GRAPHS = {'K4': wheel(3), 'W4': wheel(4), 'W5': wheel(5), 'prism': PRISM, 'Z5': cores.zigzag(5), 'K33': K33}


def to_labels(val):
    out = {}
    for m, c in val.items():
        if not m:
            out['1'] = out.get('1', 0) + c
        elif len(m) == 1:
            out[m[0]] = out.get(m[0], 0) + c
        else:
            out['*'.join(m)] = c
    return out


def orbits(edges, q):
    n = gutil.cycle_rank(edges)
    autos = gutil.edge_automorphisms(tuple(tuple(e) for e in edges))
    reps = set()
    for a in gutil.compositions(n * q, len(edges), 1, q - 1):
        if gutil.convergent(edges, a, q):
            reps.add(gutil.orbit_rep(a, autos))
    return sorted(reps)


if __name__ == '__main__':
    name, q = sys.argv[1], int(sys.argv[2])
    mx = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
    route = sys.argv[4] if len(sys.argv) > 4 else 'any'
    E = GRAPHS[name]
    reps = orbits(E, q)
    print('%s q=%d: %d convergent orbits' % (name, q, len(reps)), flush=True)
    ok = bad = 0
    for a in reps[:mx]:
        t = time.time()
        if route == 'completion':
            cands = [p for p in seq.plans_for(list(E), list(a), q) if p['route'] != 'direct']
            if not cands:
                print('  %s  no completion route' % (a,))
                continue
            pl = cands[0]
            A = seq.run(pl)
            pref = F(1, seq.gamma_int(q))
            for x in a:
                pref *= seq.gamma_int(q - x)
            val = mzv.scal(A, pref)
        else:
            try:
                val, pl = seq.period(list(E), list(a), q)
            except ValueError as e:
                print('  %s  %s' % (a, e))
                continue
        mine = to_labels(val)
        t1 = time.time() - t
        ref = cores.period(E, a, q)
        ref = {k: F(v) for k, v in ref.items() if v}
        same = mine == ref
        ok += same
        bad += not same
        print('  %s  %s  %s  (%.1fs, %s)' % (a, 'OK ' if same else 'MISMATCH', mzv.tostr(val), t1, pl['route']), flush=True)
        if not same:
            print('      ref', ref)
    print('%s q=%d: %d agree, %d disagree' % (name, q, ok, bad))
