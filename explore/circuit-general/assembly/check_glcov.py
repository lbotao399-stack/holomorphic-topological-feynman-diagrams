"""GL(H) covariance of assembled kernels and moments (exact):  lambda -> g lambda, g in GL(H, Z),
  K(g lambda) = det(g)^n K(lambda),  m_e(g lambda) = det(g)^n g m_e(lambda),
  m_ef(g lambda) = det(g)^n (g (x) g) m_ef(lambda)     [drv: y -> g y pointwise, D -> det(g)^n D].
usage: python3 check_glcov.py SET ORDER IDX [IDX ...]"""
import random
import sys
import time
from fractions import Fraction as F

import sympy

import sets
import assembly
from run_batch import HSET


def rand_g(H, rng):
    while True:
        g = [[rng.randint(-2, 2) for _ in range(H)] for _ in range(H)]
        d = sympy.Matrix(g).det()
        if d != 0 and abs(d) <= 3:
            return g, int(d)


def main(argv):
    name, order = argv[0], int(argv[1])
    H = HSET[name]
    its = dict(sets.items(name))
    for idx in [int(x) for x in argv[2:]]:
        G, O, S0 = its[idx]
        Gh, lab, lam, full = sets.instance(name, idx, G, O, S0, H)
        rng = random.Random(100 + idx)
        g, d = rand_g(H, rng)
        full2 = [[sum(g[i][j] * v[j] for j in range(H)) for i in range(H)] for v in full]
        t = time.time()
        r1 = assembly.kernel_moments(Gh.nv, Gh.edges, full, order=order, verbose=False)
        r2 = assembly.kernel_moments(Gh.nv, Gh.edges, full2, order=order, verbose=False)
        n = Gh.ne - Gh.nv + 1
        dn = F(d) ** n
        labs = set()
        for r in (r1, r2):
            for v in r.values():
                labs |= set(v)
        bad = cnt = 0
        for lb in labs:
            cnt += 1
            bad += r2['bare'].get(lb, 0) != dn * r1['bare'].get(lb, 0)
        for e in range(Gh.ne):
            for i in range(H):
                for lb in labs:
                    cnt += 1
                    rhs = dn * sum(g[i][a] * r1.get(('y', e, a), {}).get(lb, 0) for a in range(H))
                    bad += r2.get(('y', e, i), {}).get(lb, 0) != rhs
        if order >= 2:
            def m2(r, e, a, f, b):
                k = ('yy', e, a, f, b) if (e, a) <= (f, b) else ('yy', f, b, e, a)
                return r.get(k, {})
            for e in range(Gh.ne):
                for f in range(e, Gh.ne):
                    for i in range(H):
                        for j in range(H):
                            for lb in labs:
                                cnt += 1
                                rhs = dn * sum(g[i][a] * g[j][b] * m2(r1, e, a, f, b).get(lb, 0) for a in range(H) for b in range(H))
                                bad += m2(r2, e, i, f, j).get(lb, 0) != rhs
        print('%s #%d order %d: det g = %d, %d values compared, %d mismatches; bare %s (%.1fs)' % (
            name, idx, order, d, cnt, bad, {k: str(v) for k, v in r1['bare'].items()}, time.time() - t), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
