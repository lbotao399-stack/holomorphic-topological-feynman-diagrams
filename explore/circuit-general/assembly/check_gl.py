"""GL(3) covariance check of assembled first moments: at the e3 sources lambda and at standard sources
(legs e_1, e_2, e_3, root leg -e_1-e_2-e_3):  K(lambda) = det^4 K(std),  m_e(lambda) = det^4 sum_a m_e^a(std) lambda_a."""
import sys
import time
from fractions import Fraction as F

import sets
import assembly
from run_batch import std_sources

idxs = [int(x) for x in sys.argv[1:]] or [12, 40]
its = dict(sets.items('bulk3'))
for idx in idxs:
    G, O, S0 = its[idx]
    Gh, lab, lam, full = sets.instance('bulk3', idx, G, O, S0, 3)
    legs, fstd = std_sources(Gh, lab, S0, O, 3)
    t = time.time()
    r1 = assembly.kernel_moments(Gh.nv, Gh.edges, full, order=1, verbose=False)
    r2 = assembly.kernel_moments(Gh.nv, Gh.edges, fstd, order=1, verbose=False)
    L = [full[v] for v in legs]
    det = (L[0][0] * (L[1][1] * L[2][2] - L[1][2] * L[2][1]) - L[0][1] * (L[1][0] * L[2][2] - L[1][2] * L[2][0])
           + L[0][2] * (L[1][0] * L[2][1] - L[1][1] * L[2][0]))
    d4 = F(det) ** 4
    bad = 0
    labs = {'1', 'z3', 'z5'}
    for lb in labs:
        if r1['bare'].get(lb, 0) != d4 * r2['bare'].get(lb, 0):
            bad += 1
    for e in range(Gh.ne):
        for i in range(3):
            for lb in labs:
                lhs = r1.get(('y', e, i), {}).get(lb, 0)
                rhs = d4 * sum(r2.get(('y', e, a), {}).get(lb, 0) * L[a][i] for a in range(3))
                if lhs != rhs:
                    bad += 1
    print('bulk3 #%d det %d: bare %s, std c = %s; GL mismatches %d (%.1fs)' % (
        idx, det, {k: str(v) for k, v in r1['bare'].items()}, {k: str(v) for k, v in r2['bare'].items()}, bad, time.time() - t), flush=True)
