"""Convention checks of the Q_2 graph rule (log: check_conventions.log):
 (a) edge order / orientation independence: kernels recomputed by htk for a random permutation and
     random reversal of the edges of each bitriangle; pairs extracted in that order/orientation;
 (b) CS vertex: Q_2 identical with kappa != 0 and kappa = 0;
 (c) tadpole colour weight for U(1): Q2(B B)|tad and Q2(B dB)|tad versus sum q^3."""
import sys, random, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import sympy as sp
from fractions import Fraction as Fr
from htk import graph, sp_reduce
import ht3d
from ht3d import Engine, mono, pstr, psum, pscale, load_kernels, padd
from theories import *

M = 2
rng = random.Random(7)
l = sp.symbols("l0 l1 l2"); zs = sp.symbols("z0:5")
pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
alt = {}
for missing in pairs:
    edges = [p for p in pairs if p != missing]
    perm = list(range(5)); rng.shuffle(perm)
    rev = [rng.random() < 0.5 for _ in range(5)]
    aedges = [(edges[k][1], edges[k][0]) if rev[i] else edges[k] for i, k in enumerate(perm)]
    G = graph.Graph(4, aedges, "alt")
    Ks = sp_reduce(G, prefer="large").taylor([[l[0]], [l[1]], [l[2]]], {e: [zs[e]] for e in range(5)}, M, T=1, check=2)
    alt[tuple(edges)] = (tuple(aedges), [sp.expand(K) for K in Ks])
    print("bitriangle", edges, "-> alternative order/orientation", aedges)

def moments(Ks):
    mom = {}
    for K in Ks:
        K = sp.expand(K / sp.pi)
        for zmon, coeff in sp.Poly(K, *zs).terms():
            sgn = (-1) ** sum(zmon); fac = 1
            for p in zmon:
                fac *= sp.factorial(p)
            d = {}
            for lm, c in sp.Poly(sp.expand(coeff), *l).terms():
                d[lm] = Fr(int(sp.numer(c)), int(sp.denom(c))) * sgn * int(fac)
            mom[zmon] = d
    return mom

altmom = {k: (ae, moments(Ks)) for k, (ae, Ks) in alt.items()}

class AltEngine(Engine):
    def Q2_mono(self, O, only=None):
        th = self.th; out = {}
        Il = th.I_list
        for edges, (aedges, mom_tab) in altmom.items():
            for (m0, c0), (m1, c1), (m2, c2) in itertools.product(Il, Il, Il):
                mons = [m0, m1, m2, O]
                cands = []
                ok = True
                for (v, w) in aedges:
                    cl = [(i, j) for i, x in enumerate(mons[v]) for j, y in enumerate(mons[w]) if th.eta(x, y)]
                    if not cl:
                        ok = False; break
                    cands.append(cl)
                if not ok:
                    continue
                for choice in itertools.product(*cands):
                    used = set(); good = True
                    for (v, w), (i, j) in zip(aedges, choice):
                        if (v, i) in used or (w, j) in used:
                            good = False; break
                        used.add((v, i)); used.add((w, j))
                    if not good:
                        continue
                    pexp, sgn = [], 1
                    for (v, w), (i, j) in zip(aedges, choice):
                        dx, dy = mons[v][i][2], mons[w][j][2]
                        pexp.append(dx + dy); sgn *= (-1) ** dx
                    mom = mom_tab.get(tuple(pexp))
                    if not mom:
                        continue
                    legpairs = [((v, i), (w, j)) for (v, w), (i, j) in zip(aedges, choice)]
                    t = self._term(mons, [c0, c1, c2], aedges, legpairs, None, {k: v * sgn for k, v in mom.items()})
                    for m, c in t.items():
                        padd(out, m, c)
        return {m: c / 6 for m, c in out.items() if c != 0}

K4 = load_kernels(4)
print("(a) edge order/orientation independence")
for th in [su2_adj(1), u1([1, 2, -3], W={((0, 0), (1, 0), (2, 0)): 1}, kappa=Fr(1, 2))]:
    E, Ea = Engine(th, K4), AltEngine(th, K4)
    A = E.A
    pool = [('c', a, 0) for a in range(th.G)] + [('B', a, d) for a in range(th.G) for d in (0, 1)]
    for k, m in enumerate(th.multiplets):
        pool += [('p%d' % k, i, 0) for i in range(m[0])] + [('s%d' % k, i, d) for i in range(m[0]) for d in (0, 1)]
    ntest = nnz = nbad = 0
    for deg in (2, 3):
        for combo in itertools.combinations_with_replacement(pool, deg):
            if sum(L[2] for L in combo) > M:
                continue
            O = mono(A, *combo)
            if not O:
                continue
            a, b = E.Q2(O), Ea.Q2(O)
            ntest += 1; nnz += bool(a)
            if psum(a, pscale(b, -1)):
                nbad += 1
                if nbad < 3:
                    print("   MISMATCH", pstr(O), pstr(a), "|", pstr(b))
    print("   %s: %d monomials, Q2 != 0 on %d, mismatches %d" % (th.name, ntest, nnz, nbad))

print("(b) CS vertex contributes nothing to Q_2")
for th0, th1 in [(u1([1, 2, -3], W={((0, 0), (1, 0), (2, 0)): 1}), u1([1, 2, -3], W={((0, 0), (1, 0), (2, 0)): 1}, kappa=Fr(5, 3))),
                 (su2_adj(1), su2_adj(1, kappa=Fr(2, 7)))]:
    E0, E1 = Engine(th0, K4), Engine(th1, K4)
    A = E0.A
    pool = [('c', a, d) for a in range(th0.G) for d in (0, 1)] + [('B', a, d) for a in range(th0.G) for d in (0, 1)]
    for k, m in enumerate(th0.multiplets):
        pool += [('p%d' % k, i, 0) for i in range(m[0])] + [('s%d' % k, i, 0) for i in range(m[0])]
    ntest = nbad = 0
    for deg in (2, 3):
        for combo in itertools.combinations_with_replacement(pool, deg):
            O = mono(A, *combo)
            if not O:
                continue
            ntest += 1
            if psum(E0.Q2(O), pscale(E1.Q2(O), -1)):
                nbad += 1
    print("   %s kappa=%s vs 0: %d monomials, differences %d" % (th1.name, th1.kappa, ntest, nbad))

print("(c) U(1) tadpole weight")
for qs in [[1], [2], [1, 1], [1, 2, -3], [1, -1], [3, -1, -1, -1]]:
    th = u1(qs)
    E = Engine(th, K4); A = E.A
    s3 = sum(Fr(q) ** 3 for q in qs)
    r1 = E.Q2(mono(A, ('B', 0, 0), ('B', 0, 0)), 'tad')
    r2 = E.Q2(mono(A, ('B', 0, 0), ('B', 0, 1)), 'tad')
    r3 = E.Q2(mono(A, ('B', 0, 1), ('B', 0, 1)), 'tad')
    r4 = E.Q2(mono(A, ('B', 0, 0), ('B', 0, 2)), 'tad')
    print("   charges %s, sum q^3 = %s: Q2(B B)|tad = %s ; Q2(B dB)|tad = %s ; Q2(dB dB)|tad = %s ; Q2(B ddB)|tad = %s"
          % (qs, s3, pstr(r1), pstr(r2), pstr(r3), pstr(r4)))

print("(d) GL(2): tadpole weight vs symmetrized trace, and {Q0,Q2} = 0 on cubic/quadratic monomials")
import numpy as np
for nf, na in [(1, 0), (1, 1), (2, 0)]:
    th = gl_theory(2, nf, na)
    E = Engine(th, K4); A = E.A
    rhos = [mlt[1] for mlt in th.multiplets]
    def tr3(a, b, c):
        s = 0
        for rho in rhos:
            Ma = np.zeros((2, 2), dtype=object); Mb = np.zeros((2, 2), dtype=object); Mc = np.zeros((2, 2), dtype=object)
            for (i, j), v in rho[a].items(): Ma[i, j] = Fr(v)
            for (i, j), v in rho[b].items(): Mb[i, j] = Fr(v)
            for (i, j), v in rho[c].items(): Mc[i, j] = Fr(v)
            s += np.trace(Ma.dot(Mb).dot(Mc))
        return Fr(s)
    nbad = ntest = 0
    for a in range(4):
        for b in range(a, 4):
            r = E.Q2(mono(A, ('B', a, 0), ('B', b, 0)), 'tad')
            pred = {}
            for c in range(4):
                val = Fr(1, 3) * (tr3(a, b, c) + tr3(b, a, c))
                if val:
                    s, m = A.canon([('c', c, 2)])
                    pred[m] = s * val
            ntest += 1
            if psum(r, pscale(pred, -1)):
                nbad += 1
                if nbad < 4:
                    print("   B%d B%d: engine %s ; (1/3)(tr rho_a rho_b rho_c + tr rho_b rho_a rho_c) d^2 c^c = %s" % (a, b, pstr(r), pstr(pred)))
    pool = [('c', a, 0) for a in range(4)] + [('B', a, 0) for a in range(4)]
    for k, m in enumerate(th.multiplets):
        pool += [('p%d' % k, i, 0) for i in range(2)] + [('s%d' % k, i, 0) for i in range(2)]
    nz = bad = tot = 0
    for deg in (2, 3):
        for combo in itertools.combinations_with_replacement(pool, deg):
            O = mono(A, *combo)
            if not O:
                continue
            tot += 1
            q2 = E.Q2(O); nz += bool(q2)
            if psum(E.Q0(q2), E.Q2(E.Q0(O))):
                bad += 1
    print("   %s: tadpole formula mismatches %d/%d ; {Q0,Q2}: %d monomials, Q2 != 0 on %d, failures %d" % (th.name, nbad, ntest, tot, nz, bad))
