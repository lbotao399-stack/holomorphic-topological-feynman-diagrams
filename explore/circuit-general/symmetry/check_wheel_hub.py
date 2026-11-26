"""Wheel hub <-> infinity:  I_q(W_n; s, r)/prod Gamma(q - s_i) = I_q(W_n; s', r)/prod Gamma(q - s'_i),
s'_i = 2q - s_i - r_{i-1} - r_i  (rim r_i joins rim vertices i, i+1).  Direct check with wheel.py."""
from fractions import Fraction as F
from math import factorial

import engine as en
import gutil
import wheel

for n, qs in ((3, (3, 4, 5, 6, 7)), (4, (3, 4, 5, 6)), (5, (3, 4, 5)), (6, (3, 4))):
    E = [(0, i) for i in range(1, n + 1)] + [(i, i % n + 1) for i in range(1, n + 1)]
    for q in qs:
        ok = bad = 0
        for c in gutil.compositions(n * q, 2 * n, 1, q - 1):
            if not gutil.convergent(tuple(E), list(c), q):
                continue
            s, r = c[:n], c[n:]
            sp = [2 * q - s[i] - r[i - 1] - r[i] for i in range(n)]
            if min(sp) < 1 or max(sp) > q - 1 or not gutil.convergent(tuple(E), list(sp) + list(r), q):
                continue
            a = wheel.integrate(s, r, q)
            b = wheel.integrate(sp, r, q)
            fa = F(1)
            for x in s:
                fa *= factorial(q - x - 1)
            fb = F(1)
            for x in sp:
                fb *= factorial(q - x - 1)
            if {k: v / fa for k, v in a.items()} == {k: v / fb for k, v in b.items()}:
                ok += 1
            else:
                bad += 1
        print('W%d q=%d: pairs %d, mismatches %d' % (n, q, ok + bad, bad), flush=True)
