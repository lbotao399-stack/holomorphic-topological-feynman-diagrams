"""(c) n = 6, (H,T) = (1,1): engine kernels at integer sources for SP graphs with p != 3 and a few with p = 3.
Each density self-checked by the engine (check=1); pointwise check with my brute force is in c_n6_check."""
import sys, random, time, json
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
import enumerate_graphs as eg
import c_pi

which = [int(x) for x in sys.argv[1].split(",")]
gs = [Gx for Gx in eg.enumerate_tight(2, 6) if eg.is_sp(Gx)]
rng = random.Random(61)
for k in which:
    G, _ = htk.from_networkx(gs[k], None)
    red = spred.sp_reduce(G, prefer="large")
    plan = [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]
    L, p, expU = c_pi.loops_and_p(G.ne, plan)
    for trial in range(2):
        lam = [[Fr(rng.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]))] for _ in range(G.nv - 1)]
        t0 = time.time()
        dens = red.density(lam, T=1, check=1)
        K = dens.integrate()
        print("graph %d p=%d degrees %s lambda %s: K = %s  (#terms of P: %d, %.1fs)" % (
            k, p, sorted(G.degree(v) for v in range(G.nv)), [int(l[0]) for l in lam], K, len(dens.num.terms()), time.time() - t0), flush=True)
