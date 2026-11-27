"""Closed-form rooted theta moments (thetak.py) vs htk:
(a) all H = 1,2,3 rooted shapes: K(lam;0) and moments prod y_e^{p_e} (root edges, p in {0,1,2}^H) at
    random rational sources vs htk.kernel_exact(powers=...) (series-parallel reduction, brute-force self-check);
(b) H = 1, 2: full symbolic Taylor polynomial to order 2 in z on the root edges vs htk sp_reduce(...).taylor.
Log: check_thetak.log"""
import sys, random, time, itertools
from fractions import Fraction as Fr
sys.path.insert(0, ".")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import thetak
import htk
from htk import graph, sp_reduce
import sympy as sp
from math import factorial

log = open("check_thetak.log", "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()

for H in (1, 2, 3):
    shapes = thetak.rooted_shapes(H)
    out("H=%d: theta arms %s; rooted shapes (arms, root, |Aut_O|): %s" % (H, thetak.theta_arms(H), shapes))
    for arms, rs, aut in shapes:
        t0 = time.time()
        T = thetak.RootedTheta(H, arms, rs)
        s, r = T.fix_sign(seed=7, check_moments=0)
        G = T.G
        rng = random.Random(3)
        nchk = 0
        for trial in range(3 if H < 3 else 2):
            lamv = {u: [Fr(rng.randint(-6, 6)) for _ in range(H)] for u in T.surv}
            lam = [lamv.get(v, [Fr(0)] * H) for v in range(T.nv - 1)]
            key = tuple(tuple(rng.randint(0, 2 if H < 3 else 1) for _ in range(H)) for _ in T.root_edges)
            mom = T.moment(key)
            val = Fr(0)
            for beta, c in mom.items():
                t = c
                for k, u in enumerate(T.surv):
                    for i in range(H):
                        t *= lamv[u][i] ** beta[k][i]
                val += t
            powers = {(e, i): p[i] for e, p in zip(T.root_edges, key) for i in range(H) if p[i]}
            K = htk.kernel_exact(G, lam, T=1, powers=powers or None, check=1)
            ok = (K.c.get(1, Fr(0)) == val) and all(k == 1 or v == 0 for k, v in K.c.items())
            nchk += ok
            if not ok:
                out("   MISMATCH", arms, rs, key, K.c, val)
        out("  %s @ %s: |Aut_O| = %d, sign s = %+d, edges %s; moments equal to htk at %d/%d random (lam, p)  [%.1fs]"
            % (arms, rs, aut, s, T.edges, nchk, 3 if H < 3 else 2, time.time() - t0))
    if H == 3:
        continue
    # (b) symbolic Taylor to order 2
    for arms, rs, aut in shapes:
        t0 = time.time()
        T = thetak.RootedTheta(H, arms, rs)
        T.fix_sign(seed=7, check_moments=0)
        G = graph.Graph(T.nv, T.edges, "t")
        lam = [[sp.Symbol("l%d_%d" % (v, i)) if v in T.surv else sp.Integer(0) for i in range(H)] for v in range(T.nv - 1)]
        z = {e: [sp.Symbol("z%d_%d" % (e, i)) for i in range(H)] for e in T.root_edges}
        Ks = sp_reduce(G, prefer="large").taylor(lam, z, 2, T=1, check=1)
        Kh = sp.expand(sum(Ks) / sp.pi)
        # rebuild from moments: K = sum_p M_p prod (-z)^p / p!
        mine = 0
        for tot in range(3):
            for key in itertools.product(*[[p for p in itertools.product(range(3), repeat=H)] for _ in T.root_edges]):
                if sum(sum(p) for p in key) != tot:
                    continue
                mom = T.moment(key)
                zf = 1
                for e, p in zip(T.root_edges, key):
                    for i in range(H):
                        zf *= (-z[e][i]) ** p[i] / sp.Integer(factorial(p[i]))
                for beta, c in mom.items():
                    t = sp.Rational(c.numerator, c.denominator) * zf
                    for k, u in enumerate(T.surv):
                        for i in range(H):
                            t *= lam[u][i] ** beta[k][i]
                    mine += t
        d = sp.expand(Kh - mine)
        out("  (b) H=%d %s @ %s: htk Taylor (order <= 2 in z at O) - closed form = %s  [%.1fs]" % (H, arms, rs, d, time.time() - t0))
