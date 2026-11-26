"""Check of the Weyl-formula multiplicity (used in t7_sizes.py / analyze.py) against a direct computation:
dim of the joint kernel of the raising operators E_{i,i+1} on C_quot slot (w, t=(k..k), N)
== sum_s sgn(s) dim C_{t + rho - s rho}.  Writes check_weyl.log."""
import itertools
from cohom import Space
from jets import Echelon, vec
from t4_homotopy import gl_action
from t7_sizes import dimC_N, weyl_mult

out = open('check_weyl.log', 'w')


def log(s):
    print(s, flush=True)
    out.write(s + '\n')


for n, w, t, Gmaxh in ((2, 1, (-1, -1), 6), (2, 1, (-2, -2), 8), (2, 0, (0, 0), 6), (3, 1, (0, 0, 0), 6), (3, 0, (0, 0, 0), 6)):
    S = Space(n, w, Gmaxh)
    A = S.A
    sl = S.slot(t)
    for N in sorted(sl):
        src = sl[N]
        # joint kernel: stack images in the direct sum of the target slots
        tg = []
        for i in range(n - 1):
            tt = tuple(t[j] + (j == i) - (j == i + 1) for j in range(n))
            tg.append(S.slot(tt).get(N, []))
        offs = [0]
        for L in tg:
            offs.append(offs[-1] + len(L))
        idx = {}
        for i, L in enumerate(tg):
            for j, m in enumerate(L):
                idx[(i, m)] = offs[i] + j
        E = Echelon()
        ops = [gl_action(A, i, i + 1) for i in range(n - 1)]
        for m in src:
            v = {}
            for i, op in enumerate(ops):
                img = A.apply_der({m: 1}, op, False)
                for mm, c in img.items():
                    v[idx[(i, mm)]] = v.get(idx[(i, mm)], 0) + c
            E.insert(v)
        ker = len(src) - E.rank
        wm = weyl_mult(n, w, t, N, 'C')
        log('n=%d w=%d t=%s N=%2d: dim C=%6d, dim ker(E_{i,i+1}) = %6d, Weyl formula = %6d, equal: %s' % (n, w, t, N, len(src), ker, wm, ker == wm))
