"""All-entries checks (dense, modulo p < 2^19) of the cochain identities, for all 6 QB coefficients mu:
   F133^mu + F223^mu - (1/2) beta^mu = 0,   F223^mu - x beta^mu - a G_kk^mu = 0   (Killing form).
"""
import sys, time, itertools
import numpy as np
from fractions import Fraction as Fr
import vlie, vpoint, vb_dense as vd
from va_points import CLAIM

def alt_mu(Dn, W, mu):
    p = Dn.p
    out = np.zeros_like(W)
    for s in vpoint.PERMS:
        sg, c = vpoint.CQ[s]
        if c[mu]:
            out = np.fmod(out + ((sg * c[mu]) % p) * vd.permT(W, s), p)
    return out

def beta_mu(Dn, mu):
    p = Dn.p
    out = np.zeros_like(Dn.FF)
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        c = vpoint.quad(vpoint.br(k, l), vpoint.br(k, l))
        if c[mu]:
            out = np.fmod(out + ((-1) ** (i + j) * c[mu] % p) * vd.permT(Dn.FF, (i, j, k, l)), p)
    return out

for p in vd.P19[:int(sys.argv[1])]:
    for key in sys.argv[2:]:
        t0 = time.time()
        L = vlie.build(key); d = L.d
        V, K, Kinv = L.mod(p)
        Cp = np.zeros((d, d, d), dtype=np.int64); np.add.at(Cp, (L.A, L.B, L.Cc), V); Cp %= p
        Dn = vd.Dense(d, Cp, K, Kinv, p); Dn.build()
        KKl = np.fmod(np.einsum('ab,ce->abce', Dn.K, Dn.K), p)
        h = vd.r(Fr(1, 2), p)
        x, a = CLAIM.get(key, (None, None))
        ok1 = ok2 = True; nz = []
        for mu in range(6):
            F1 = alt_mu(Dn, Dn.W133, mu); F2 = alt_mu(Dn, Dn.W223, mu); B = beta_mu(Dn, mu)
            ok1 &= not np.fmod(F1 + F2 + (p - h) * B, p).any()
            if x is not None:
                G = alt_mu(Dn, KKl, mu)
                ok2 &= not np.fmod(F2 + (p - vd.r(x, p)) * B + (p - vd.r(a, p)) * G, p).any()
            nz.append(int(np.count_nonzero(F2)))
        print('%-16s p=%d d=%d  all entries, 6 coefficients: F133+F223 = beta/2: %s ; F223 = x beta + a G_kk (x,a)=(%s,%s): %s ; nonzero entries of F223^mu: %s  (%.1fs)'
              % (L.name, p, d, ok1, x, a, ok2 if x is not None else '-', nz, time.time() - t0), flush=True)
