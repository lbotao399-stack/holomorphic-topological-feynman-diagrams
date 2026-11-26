"""Miscellaneous checks.
 (1) dictionary forms: alt((1/3) p1^1 p2^2 d) = D;  alt((1/3)((p1^1)^2 (p2^2)^2 - p1^1 p1^2 p2^1 p2^2) f) = u2f;
     alt((1/3)((p1^1)^3 (p2^2)^3 - 3 (p1^1)^2 p1^2 p2^1 (p2^2)^2) d) = D3   (gl_N and sl_N).
 (2) tree graph (one heat edge, vertices a = v0, b = v1 = root, externals x1,x2 at v0, x3,x4 at v1): own Wick sign
     sigma, w_tree = sigma a(L0) kappa^{-1} b(L1) with q_(e,0) = -y_e = -(p1+p2), q_(e,1) = +y_e; check
     (1/4) alt(w_tree) = {a,b} for all 36 pairs of local orders, (a,b) = (D,u2f), (D,D), (s0,u2f), (S1,u2f) (gl_N).
 (3) vanishing (engine): polygon(4) at (H,T) = (2,1) with current insertions; tight three-loop graphs (catalog).
 (4) d D3 explicit formula of the explorer: dD3(x0..x3) = sum_{i<j} (-1)^{i+j} [p_k p_l]^3 tr([X_i,X_j]{X_k,X_l}).
"""
import itertools, random, sys, time
from fractions import Fraction as Fr
import numpy as np
from vlib import *
import vexp

t0 = time.time()
for sl in (False, True):
    tag = 'sl' if sl else 'gl'
    def mono3(a0, b0, a1, b1):
        e = [0] * 6
        e[0], e[3], e[1], e[4] = a0, b0, a1, b1
        return tuple(e)
    PD = {mono3(1, 0, 0, 1): Fr(1, 3)}
    Pu = {mono3(2, 0, 0, 2): Fr(1, 3), mono3(1, 1, 1, 1): Fr(-1, 3)}
    P3 = {mono3(3, 0, 0, 3): Fr(1, 3), mono3(2, 1, 1, 2): Fr(-1)}
    for name, P, T, ref in [("D", PD, 'd', cubic('d', 1, sl)), ("u2f", Pu, 'f', cubic('f', 2, sl)), ("D3", P3, 'd', cubic('d', 3, sl))]:
        w = from_raw(3, [(LP.const(c), words, P) for c, words in TENS3[T]], sl)
        print("%s: alt(P W) = %s: %s" % (tag, name, is_zero(lin((1, alt(w, 3, sl)), (-1, ref)))))

# (2) tree
SL = False
D = cubic('d', 1, SL); u2f = cubic('f', 2, SL); s0 = cubic('f', 0, SL); S1 = cubic('s1', 1, SL)
alg = vexp.Algebra(3, SL)
random.seed(11)
pts = []
for _ in range(2):
    Xs = [alg.random_element(2) for _ in range(4)]
    ps = [(random.randint(-3, 3), random.randint(-3, 3)) for _ in range(3)]
    ps.append((-sum(p[0] for p in ps), -sum(p[1] for p in ps)))
    pts.append((Xs, ps))
src0 = [('e', 0), ('x', 1), ('x', 2)]
src1 = [('e', 1), ('x', 3), ('x', 4)]
target = [('e', 0), ('e', 1), ('x', 1), ('x', 2), ('x', 3), ('x', 4)]
for (an, a), (bn, b) in [(("D", D), ("u2f", u2f)), (("D", D), ("D", D)), (("s0", s0), ("u2f", u2f)), (("S1", S1), ("u2f", u2f))]:
    br = bracket(a, 3, b, 3, SL)
    allok = True
    vals = []
    for L0 in itertools.permutations(src0):
        for L1 in itertools.permutations(src1):
            pos = {h: i for i, h in enumerate(target)}
            sigma = seq_sign([pos[h] for h in list(L0) + list(L1)])
            for Xs, ps in pts:
                # (1/4) alt(w_tree): sum over sigma in S4 of sgn * w_tree(x_sigma)
                tot = Fr(0)
                for perm in itertools.permutations(range(4)):
                    sg = perm_sign(perm)
                    X = [Xs[perm[i]] for i in range(4)]
                    P = [ps[perm[i]] for i in range(4)]
                    ye = (P[0][0] + P[1][0], P[0][1] + P[1][1])
                    def arg(h, e_mat):
                        if h[0] == 'x':
                            return X[h[1] - 1], P[h[1] - 1]
                        if h[1] == 0:
                            return e_mat, (-ye[0], -ye[1])
                        return e_mat, ye
                    s = Fr(0)
                    for e in range(alg.dim):
                        A_args = [arg(h, alg.basis[e]) for h in L0]
                        av = vexp.evaluate(a, 3, [m for m, _ in A_args], [p for _, p in A_args], 3)
                        if av == 0:
                            continue
                        for e2 in range(alg.dim):
                            if alg.kinv[e][e2] == 0:
                                continue
                            B_args = [arg(h, alg.basis[e2]) for h in L1]
                            bv = vexp.evaluate(b, 3, [m for m, _ in B_args], [p for _, p in B_args], 3)
                            s += av * alg.kinv[e][e2] * bv
                    tot += sg * sigma * s
                ref = vexp.evaluate(br, 4, Xs, ps, 3)
                vals.append((tot / 4, ref))
                allok = allok and (tot / 4 == ref)
    print("tree (%s,%s): (1/4) alt(sigma a(L0) b(L1)) == {a,b} for 36 local-order pairs x 2 points: %s  (e.g. %s vs %s)"
          % (an, bn, allok, vals[0][0], vals[0][1]), flush=True)

# (3) vanishing
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
import sympy as sp
G = htk.polygon(4)
lam = [[Fr(1), Fr(2)], [Fr(-1), Fr(3)], [Fr(2), Fr(-5, 3)]]
red = htk.sp_reduce(G, prefer="large")
for pw in [None, {(0, 0): 1}, {(0, 0): 1, (1, 1): 1}, {(0, 0): 2, (2, 1): 2}, {(1, 0): 1, (3, 1): 1, (2, 0): 1, (0, 1): 1}]:
    try:
        K = red.integrate(lam, 1, powers=pw, check=2)
        print("polygon(4), (H,T)=(2,1), insertion %s: K = %s" % (pw, K))
    except Exception as ex:
        print("polygon(4) insertion %s: engine error %s" % (pw, ex))
try:
    gs = htk.catalog.tight_graphs(3, 3)
    print("tight three-loop graphs at D = 3: %d" % len(gs))
    nz = 0
    for i, g in enumerate(gs):
        lamg = [[Fr(random.randint(-4, 4)), Fr(random.randint(-4, 4))] for _ in range(g.nv - 1)]
        t = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(g.ne)]
        dens = htk.brute.density(g, t, lamg, T=1)
        z = (dens.coeff == 0)
        nz += (not z)
        print("  graph %d (nv %d, ne %d): brute T=1 density at a random rational t zero: %s" % (i, g.nv, g.ne, z))
    print("  nonzero densities: %d" % nz)
except Exception as ex:
    print("three-loop check error:", ex)

# (4) explicit formula for d D3
for sl in (False, True):
    alg = vexp.Algebra(3, sl)
    D3 = cubic('d', 3, sl)
    dD3 = d(D3, 3, sl)
    ok = True
    for _ in range(4):
        Xs = [alg.random_element(2) for _ in range(4)]
        ps = [(random.randint(-3, 3), random.randint(-3, 3)) for _ in range(3)]
        ps.append((-sum(p[0] for p in ps), -sum(p[1] for p in ps)))
        tot = Fr(0)
        for i in range(4):
            for j in range(i + 1, 4):
                k, l = [m for m in range(4) if m not in (i, j)]
                b = ps[k][0] * ps[l][1] - ps[k][1] * ps[l][0]
                C = Xs[i].dot(Xs[j]) - Xs[j].dot(Xs[i])
                A = Xs[k].dot(Xs[l]) + Xs[l].dot(Xs[k])
                tot += (-1) ** (i + j) * Fr(b) ** 3 * np.trace(C.dot(A))
        ok = ok and (tot == vexp.evaluate(dD3, 4, Xs, ps, 3))
    print("%s_3: d D3 = sum_{i<j} (-1)^{i+j} [p_k p_l]^3 tr([X_i,X_j]{X_k,X_l}) at 4 points: %s" % ('sl' if sl else 'gl', ok))
print("time %.1fs" % (time.time() - t0))
