"""C2/C7 for unsorted path lengths: theta_std(perm of (p,q,r)) (B = the shared path of C3 can be any of the three)."""
import sys, random, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/verify")
sys.dont_write_bytecode = True
from vcore import *
rng = random.Random(99)
ok, n = True, 0
for T, H in [(0, 2), (0, 3), (0, 4), (1, 1), (1, 2), (1, 3), (1, 4)]:
    for pqr in theta_list(H + T):
        for perm in sorted(set(itertools.permutations(pqr))):
            th = theta_std(*perm)
            for _ in range(2):
                t = rand_t(th.ne, rng)
                lam = rand_sources(th.nv, th.Np, H, rng)
                bf, U, _ = density_bf(th, t, lam, H, T, rng=rng)
                cl, _ = density_closed(th, t, lam, H, T)
                ok &= (bf == cl)
            Bf = incidence_full(th.nv, th.edges)
            d = det([[Bf[e][v] for e in range(th.ne)] for v in range(th.nv) if v != th.Np] + th.C3())
            ok &= ((1 if d > 0 else -1) == sg(perm[1] + 1))
            n += 1
print("unsorted (p,q,r) orders: %d graphs, brute = closed density and sgn det = (-1)^(q+1): %s" % (n, ok))
