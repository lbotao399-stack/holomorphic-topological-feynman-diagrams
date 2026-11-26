"""C6 lemma: B(u,v) = det[(u;0),(0;v),(Lam_A;0),(Lam_B;Lam_B),(0;Lam_C)] on span(Lam_P) x span(Lam_P):
zero for P = A, C; for P = B antisymmetric and in general nonzero (only the symmetric part vanishes)."""
import sys, random
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/verify")
sys.dont_write_bytecode = True
from vcore import *
rng = random.Random(8)
for H, pqr in [(3, (2, 2, 3)), (4, (2, 3, 4)), (4, (3, 3, 3)), (5, (3, 4, 4)), (6, (4, 4, 5))]:
    th = theta_std(*pqr)
    lam = rand_sources(th.nv, th.Np, H, rng)
    def Bf(u, v):
        z = [0] * H
        cols = [col(u, z, H), col(z, v, H)] + [col(lam[w], z, H) for w in th.interior["A"]] + \
               [col(lam[w], lam[w], H) for w in th.interior["B"]] + [col(z, lam[w], H) for w in th.interior["C"]]
        return det([[c[i] for c in cols] for i in range(2 * H)])
    out = {}
    for P in "ABC":
        I = th.interior[P]
        vals = [(Bf(lam[a], lam[b]), Bf(lam[b], lam[a])) for a in I for b in I]
        out[P] = ("all zero" if all(x == 0 and y == 0 for x, y in vals) else
                  "antisymmetric, max |B| = %s" % max(abs(x) for x, _ in vals) if all(x == -y for x, y in vals) else "neither")
    print("H=%d %s:" % (H, pqr), out)
