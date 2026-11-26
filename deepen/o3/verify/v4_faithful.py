"""Function-space dimension of fixed-point-free trace monomials on sl_N^k (and all trace monomials on gl_N^k)."""
import time
from fractions import Fraction as Fr
import vla
from vfix import Faithful
t0 = time.time()
for (N, sl, k, mode) in [(3, True, 3, 'all'), (3, True, 4, 'all'), (3, True, 5, 'all'),
                         (4, True, 4, 'random'), (4, True, 5, 'random'), (5, True, 4, 'random'), (5, True, 5, 'random'),
                         (3, False, 4, 'random'), (5, False, 5, 'random')]:
    F = Faithful(N, sl, k, mode=mode)
    extra = ""
    if F.rank < len(F.perms) and mode != 'all':
        rel = F.antisym_relation()
        # check the relation exactly on all evaluation tuples and that rank + 1 = #perms
        ok = all(sum(c * F.E[i][j] for i, c in rel.items()) == 0 for j in range(F.ntuples))
        extra = "; antisymmetrizer relation vanishes on all %d tuples: %s; #perms - 1 = %d" % (F.ntuples, ok, len(F.perms) - 1)
    print("%s%d k=%d (%s, %d tuples): #taus %d, rank %d%s" % ('sl' if sl else 'gl', N, k, mode, F.ntuples, len(F.perms), F.rank, extra), flush=True)
print("time %.1fs" % (time.time() - t0))
