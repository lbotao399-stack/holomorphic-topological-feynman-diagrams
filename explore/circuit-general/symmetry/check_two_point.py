"""Two-point (numerator) formula: with c at infinity and a, b at 0, 1 the edge ab drops out;
b_ab = -k <= 0:  Pi(G) = Gamma(q+k)/prod_{H'} Gamma(b_e) sum_m c_m I_{q+k}(H'; q+k-b-m),
(Psi^{ab}_{H'})^k = sum_m c_m t^m (2-forests separating a, b),  H' = G - c - ab.
Check: completed graphs of K4, W4, W5, prism, K33 having an edge with b <= 0; compare with a clean
decompletion (engine.terminal)."""
import time
from collections import Counter

import engine as en
import gutil
from check_completion import K4, W4, W5, PRISM, K33

tot = Counter()
for name, G, qs in [('K4', K4, (3, 4, 5, 6)), ('W4', W4, (3, 4, 5)), ('W5', W5, (3, 4)),
                    ('prism', PRISM, (3, 4, 5)), ('K33', K33, (3, 4, 5))]:
    n = gutil.cycle_rank(G)
    autos = gutil.edge_automorphisms(tuple(G))
    for q in qs:
        t = time.time()
        reps = set()
        for c in gutil.compositions(n * q, len(G), 1, q - 1):
            if gutil.convergent(tuple(G), list(c), q):
                reps.add(gutil.orbit_rep(c, autos))
        cnt = Counter()
        for al in sorted(reps):
            V, E, inf = en.complete(G, al, q)
            r0 = en.terminal(V, E, q)
            if r0 is None:
                continue
            for c in range(V):
                for (a, b), w in E.items():
                    if c in (a, b) or w > 0:
                        continue
                    if any(x < 0 for (i, j), x in E.items() if c not in (i, j) and (i, j) != (a, b)):
                        continue
                    r = en.pi_two_point(V, E, c, (a, b), q, kmax=6)
                    if r is None:
                        continue
                    k = -w
                    cnt['k=%d' % k] += 1
                    if r[0] != r0[0]:
                        cnt['mismatch'] += 1
                        print('MISMATCH', name, q, al, c, (a, b), en.vstr(r[0]), '|', en.vstr(r0[0]), flush=True)
        print('%-5s q=%d orbits %d  two-point evaluations compared %s  (%.1fs)' % (name, q, len(reps), dict(cnt), time.time() - t), flush=True)
        tot.update(cnt)
print('total', dict(tot))
