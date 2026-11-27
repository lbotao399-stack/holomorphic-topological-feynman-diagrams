"""Families so_{2n}, so_{2n+1}, sp_{2n} (lie.so_tr / sp_tr, kappa = tr_V, diagonal Cartan basis eps_i):
Cartan detector o'|_C with h_v supported on eps_1..eps_3 and fixed momenta, as a function of n.  It is a polynomial
in n (root sums over the remaining indices); degree found by interpolation and verified on the other n.
No root in F_p (two primes)  =>  o'|_C != 0, hence [O_2] != 0, for every n in the family.
Same for sl_N (lie.sl_tr, h on the first 3 simple coroots) as a control of the formal result."""
import numpy as np
from common import Log, PRIMES
import lie
import cartan_lib as cl
import polymod as fp

log = Log(__file__)


def value(L, p, hsmall, mom):
    det_ = cl.Detector(L, p)
    r = len(L.cartan)
    hc = np.zeros((8, r), dtype=np.int64)
    hc[:, :hsmall.shape[1]] = hsmall
    return det_.o_prime(hc, mom)


rng = np.random.default_rng(123)
hsmall = rng.integers(-4, 5, size=(8, 3))
mom = rng.integers(-4, 5, size=(7, 4)); mom = np.vstack([mom, -mom.sum(axis=0)])
fams = {
    'so(2n)':   (lambda n: lie.so_tr(2 * n), range(4, 13)),
    'so(2n+1)': (lambda n: lie.so_tr(2 * n + 1), range(4, 13)),
    'sp(2n)':   (lambda n: lie.sp_tr(n), range(4, 13)),
    'sl(n+1)':  (lambda n: lie.sl_tr(n + 1), range(4, 13)),
}
from common import crt_pair, ratrec
for name, (mk, ns) in fams.items():
    ns = list(ns)
    roots = []
    for p in PRIMES[:4]:
        vals = [value(mk(n), p, hsmall, mom) for n in ns]
        deg = None
        for d in range(0, 6):
            poly = fp.interp(ns[:d + 1], vals[:d + 1], p)
            if all(fp.peval(poly, n, p) == v for n, v in zip(ns, vals)):
                deg = d
                break
        nr = fp.n_roots_mod(poly, p) if deg is not None else None
        log('%-9s p = %d: values at n = %s nonzero: %s; polynomial degree %s (verified on all %d n); roots in F_p: %s'
            % (name, p, ns, all(vals), deg, len(ns), nr))
        if deg == 1:
            roots.append((-poly[0] * pow(poly[1], p - 2, p) % p, p))
    if len(roots) == 4:
        v, m = roots[0]
        for a, q in roots[1:]:
            v, m = crt_pair(v, m, a, q)
        log('%-9s linear polynomial in n: root (CRT over 4 primes, rational reconstruction) = %s' % (name, ratrec(v, m)))
log.close()
