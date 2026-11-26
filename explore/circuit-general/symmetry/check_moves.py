"""Exact check of the single moves on completed graphs of tabulated cores.
For each convergent orbit (G, a, q) of K4, W4, W5, prism, K33 and every applicable move
G -> G' (M1 star->triangle, M2 triangle->star, FT Fourier/planar duality of a decompletion, TW twist at a
4-vertex cut, PR product at a 3-vertex cut): compare Pi(G) and factor * Pi(G') whenever both sides
have a tabulated decompletion (terminal evaluation only, no search)."""
import sys
import time
from collections import Counter
from fractions import Fraction as F

import engine as en
import gutil
from check_completion import K4, W4, W5, PRISM, K33

W6 = [(0, i) for i in range(1, 7)] + [(i, i % 6 + 1) for i in range(1, 7)]
CASES = [('K4', K4, (3, 4, 5, 6)), ('W4', W4, (3, 4, 5)), ('W5', W5, (3, 4)), ('W6', W6, (3,)),
         ('prism', PRISM, (3, 4, 5)), ('K33', K33, (3, 4, 5))]
if len(sys.argv) > 1:
    CASES = [c for c in CASES if c[0] in sys.argv[1:]]

tot = Counter()
for name, G, qs in CASES:
    n = gutil.cycle_rank(G)
    autos = gutil.edge_automorphisms(tuple(G))
    for q in qs:
        t = time.time()
        cnt = Counter()
        bad = Counter()
        reps = set()
        for c in gutil.compositions(n * q, len(G), 1, q - 1):
            if gutil.convergent(tuple(G), list(c), q):
                reps.add(gutil.orbit_rep(c, autos))
        for al in sorted(reps):
            V, E, inf = en.complete(G, al, q)
            r0 = en.terminal(V, E, q)
            if r0 is None:
                continue
            v0 = r0[0]
            for V2, E2, f2, desc in en.successors(V, E, q, use_m2=True, split_m2=True):
                kind = desc[:2]
                if not en.balanced(V2, E2, q):
                    bad[kind + '-unbalanced'] += 1
                    continue
                r2 = en.terminal(V2, E2, q)
                if r2 is None:
                    continue
                cnt[kind] += 1
                if en.canon(V2, E2) != en.canon(V, E):
                    cnt[kind + '-nontrivial'] += 1
                if en.vscale(r2[0], f2) != v0:
                    bad[kind] += 1
                    if bad[kind] <= 3:
                        print('MISMATCH', name, q, al, desc, en.vstr(v0), '|', en.vstr(en.vscale(r2[0], f2)), flush=True)
            ps = en.product_split(V, E, q)
            if ps is not None:
                (P1, P2), S = ps
                r1 = en.terminal(P1[0], P1[1], q)
                r2 = en.terminal(P2[0], P2[1], q)
                if r1 is not None and r2 is not None:
                    cnt['PR'] += 1
                    if en.vmul(r1[0], r2[0]) != v0:
                        bad['PR'] += 1
                        if bad['PR'] <= 3:
                            print('MISMATCH PR', name, q, al, S, flush=True)
        print('%-5s q=%d orbits %d  checks %s  mismatches %s  (%.1fs)' % (name, q, len(reps), dict(cnt), dict(bad),
                                                                         time.time() - t), flush=True)
        tot.update(cnt)
print('total checks', dict(tot))
