"""Exact checks of engine values on five- and six-loop cores, by identities NOT used by the engine.
 (A) dimension shift  I_{q-1}(G; a) = sum_C I_q(G; a + 1_C)  (C cotrees), every term convergent and
     engine-evaluated, both sides through independent chains;
 (B) phi^4 limit q = 2, a = 1:  P(Z5) = 441/8 z7, P(W5) = 70 z7, P(K4+K4) = P(K33+e) = 36 z3^2,
     from level-3 engine values through (A);
 (C) Z5 at q = 4: engine values vs the orbits determined by the msys relation system."""
import json
import os
import sys
import time
from fractions import Fraction as F

import engine as en
import gutil
import msys
import cores

HERE = os.path.dirname(os.path.abspath(__file__))
SK = {r['name']: [tuple(e) for e in r['edges']] for r in json.load(open(os.path.join(HERE, '..', 'fiveloop', 'skel5.json')))}
SK['K5'] = [(i, j) for i in range(5) for j in range(i + 1, 5)]


def vadd(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}


def dim_shift(name, G, q, maxn=None):
    """all convergent orbits at level q-1 (entries >= 1) whose cotree shifts are all engine-evaluated"""
    n = gutil.cycle_rank(G)
    autos = gutil.edge_automorphisms(tuple(G))
    cot = gutil.cotrees(tuple(G))
    reps = set()
    for c in gutil.compositions(n * (q - 1), len(G), 1, q - 2):
        if gutil.convergent(tuple(G), list(c), q - 1):
            reps.add(gutil.orbit_rep(c, autos))
    ok = bad = skipped = 0
    example = None
    for al in sorted(reps):
        lhs = en.period(G, al, q - 1)
        if lhs is None:
            skipped += 1
            continue
        tot = {}
        good = True
        for C in cot:
            b = list(al)
            for e in C:
                b[e] += 1
            if not gutil.convergent(tuple(G), b, q):
                good = False
                break
            r = en.period(G, b, q)
            if r is None:
                good = False
                break
            tot = vadd(tot, r[0])
        if not good:
            skipped += 1
            continue
        if tot == lhs[0]:
            ok += 1
            if example is None:
                example = (al, lhs)
        else:
            bad += 1
            print('MISMATCH', name, q, al, en.vstr(lhs[0]), '|', en.vstr(tot), flush=True)
        if maxn and ok + bad >= maxn:
            break
    return len(reps), ok, bad, skipped, example


if __name__ == '__main__':
    t0 = time.time()
    print('(A) dimension shift q-1 -> q', flush=True)
    # cores with full coverage (every cotree term evaluated); partial cores: check_routes.py
    for name in ['K5-e', 'W5', 'K4+K4', 'V7_1', 'V8_0']:
        for q in (3, 4, 5):
            t = time.time()
            if name == 'K5' and q > 4:
                continue
            nr, ok, bad, sk, ex = dim_shift(name, SK[name], q)
            print('  %-6s q=%d->%d: level-%d orbits %d, checked %d, mismatches %d, not all terms evaluated %d  (%.1fs)'
                  % (name, q - 1, q, q - 1, nr, ok, bad, sk, time.time() - t), flush=True)
            if ex:
                print('        e.g. I_%d%s = %s   chain %s' % (q - 1, ex[0], en.vstr(ex[1][0]), ' '.join(ex[1][1])), flush=True)
    print('(B) q = 2, unit weights, from level 3', flush=True)
    ref = {'Z5': {('z7',): F(441, 8)}, 'W5': {('z7',): F(70)}, 'K4+K4': {('z3', 'z3'): F(36)},
           'K33+e': {('z3', 'z3'): F(36)}}
    for name, R in ref.items():
        G = SK[name]
        tot = {}
        good = True
        for C in gutil.cotrees(tuple(G)):
            b = [1] * len(G)
            for e in C:
                b[e] += 1
            r = en.period(G, b, 3) if gutil.convergent(tuple(G), b, 3) else None
            if r is None:
                good = False
                break
            tot = vadd(tot, r[0])
        direct = en.period(G, [1] * len(G), 2)
        print('  %-6s sum_C I_3(1+1_C) = %s  (all terms evaluated: %s)  expected %s  match %s;  engine at q=2 directly: %s'
              % (name, en.vstr(tot) if good else '-', good, en.vstr(R), good and tot == R,
                 en.vstr(direct[0]) if direct else None), flush=True)
    print('(C) Z5 q=4 vs msys-determined orbits', flush=True)
    Z5 = cores.zigzag(5)
    S = msys.MultiSystem([Z5], 4)
    S.lower = cores.lower_core_value
    S.build(verbose=False)
    solved, piv = S.solve(verbose=False)
    ok = bad = ev = 0
    for j, (f, rep) in enumerate(S.reps):
        r = en.period(Z5, rep, 4)
        if r is None:
            continue
        ev += 1
        if j in solved:
            if r[0] == en.from_cores(solved[j]):
                ok += 1
            else:
                bad += 1
                print('MISMATCH Z5', rep, flush=True)
    print('  Z5 q=4: orbits %d, msys-determined %d, engine-evaluated %d, compared %d, mismatches %d'
          % (len(S.reps), len(solved), ev, ok + bad, bad), flush=True)
    print('total %.1fs' % (time.time() - t0))
