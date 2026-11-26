"""Labelling test of the sign rule
    sign(W^pattern_Gamma relative to I_Gamma) = koszul(Gamma) * (-1)^{#spokes oriented O -> v},
koszul = sign of the permutation taking the vertex fields (vertex order, legs in the pattern order)
to (tree-edge pairs (tail field, head field) in edge order, spoke fields in delta order, outputs in
vertex order).  For random relabellings (vertex permutation, edge reversals, edge order) of graphs
the kernel sign is recomputed by brute force and the relation is solved:
    n = 2 (graphs A, B relabelled):  kappa_A, kappa_B   predicted  rule / 2
    n = 3 (argument 'n3'):           kappa_G1..G5       predicted  rule / |Aut|."""
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

import opcheck as oc
import q3lib
from q3lib import Reduction, Brute, LAM, R, random_point, lam_sub

rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)


def relabel(name, new):
    g = q3lib.ALL[name]
    root = g['root']
    E = g['edges']
    N = len(E)
    perm = list(range(root))
    rng.shuffle(perm)
    vmap = {v: perm[v] for v in range(root)}
    vmap[root] = root
    rev = [rng.random() < 0.5 for _ in range(N)]
    order = list(range(N))
    rng.shuffle(order)                      # new edge k = old edge order[k]
    E2 = []
    for k in range(N):
        a, b = E[order[k]]
        a, b = vmap[a], vmap[b]
        if rev[order[k]]:
            a, b = b, a
        E2.append((a, b))
    inv = {old: k for k, old in enumerate(order)}
    sp2 = []
    for e, s in g['spokes']:
        k = inv[e]
        sp2.append((k, +1 if E2[k][0] == root else -1))
    mv2 = [(m[0],) + tuple(vmap[x] for x in m[1:]) for m in g['moves']]
    q3lib.ALL[new] = dict(edges=E2, spokes=sp2, moves=mv2, root=root)
    # kernel sign by brute force
    red = Reduction(new)
    B = Brute(red.edges, red.cotree(), root)
    nr = 3 * len(red.steps)
    dom = red.domain()
    facs = red.factors()
    dt = [[sp.diff(te, r) for r in R[:nr]] for te in red.t]
    for _ in range(3):
        lam = [[Fr(rng.randint(-5, 5)) for _ in range(3)] for _ in range(root)]
        pt = random_point(dom, rng)
        subs = {r: sp.Rational(pt[r].numerator, pt[r].denominator) for r in R[:nr]}
        t = [Fr(int(sp.Rational(te.subs(subs)).p), int(sp.Rational(te.subs(subs)).q)) for te in red.t]
        Jt = [[Fr(int(sp.Rational(d.subs(subs)).p), int(sp.Rational(d.subs(subs)).q)) for d in row] for row in dt]
        om = q3lib.fdet([t] + [[Jt[e][k] for e in range(N)] for k in range(nr)])
        y, D = B.dy_dt(t, lam)
        rows = [[sum((D[c][a][j] * Jt[j][kk] for j in range(N)), Fr(0)) for kk in range(nr)] for c in B.cot for a in range(3)]
        val = B.epsT * B.detBT * q3lib.fdet(rows) * (1 if om > 0 else -1)
        dp = sp.Rational(sp.prod([f_.subs(lam_sub(lam)).subs(subs) for f_ in facs]))
        if dp != 0:
            s = val / Fr(int(dp.p), int(dp.q))
            assert s in (1, -1)
            oc.SIGNS[new] = int(s)
            break
    oc.AUT[new] = oc.AUT[name]
    return new


def rule(name):
    G = oc.op(name)
    red = G.red
    E = red.edges
    root = red.root
    fields = [(v, typ, x) for v, L in enumerate(G.legs) for (typ, x) in L]
    spoke_e = {e for e, s in red.spokes}
    target = []
    for e in range(len(E)):
        if e in spoke_e:
            continue
        a, b = E[e]
        target.append(next(i for i, fl in enumerate(fields) if fl[0] == a and fl[1] == 'R' and fl[2] == e))
        target.append(next(i for i, fl in enumerate(fields) if fl[0] == b and fl[1] == 'R' and fl[2] == e))
    target += [next(i for i, fl in enumerate(fields) if fl[1] == 'A' and fl[2] == j) for j in range(G.m)]
    target += [next(i for i, fl in enumerate(fields) if fl[1] == 'D' and fl[0] == v) for v in G.outs]
    nout = sum(1 for e, s in red.spokes if E[e][0] == root)
    return oc.perm_sign(target) * (-1) ** nout


def solve(cols, rhs_vals):
    ks = sp.symbols('k1:%d' % (len(cols) + 1))
    M = sp.Matrix([[sp.Rational(c[i].numerator, c[i].denominator) for c in cols] for i in range(len(rhs_vals))])
    rhs = sp.Matrix([-sp.Rational(x.numerator, x.denominator) / 4 for x in rhs_vals])
    return M.rank(), M.row_join(rhs).rank(), sp.linsolve((M, rhs), *ks)


n, f, fd = oc.fdense('so5')
mode = sys.argv[1] if len(sys.argv) > 1 else 'n2'
print('rule on the reference labellings:', {g: rule(g) for g in ['box', 'A', 'B', 'G1', 'G2', 'G3', 'G4', 'G5']})
if mode == 'n2':
    for trial in range(6):
        t0 = time.time()
        a = relabel('A', 'A_%d' % trial)
        b = relabel('B', 'B_%d' % trial)
        rs = np.random.RandomState(trial)
        Ov = [rs.randint(-2, 3, size=n).astype(np.int64) for _ in range(3)]
        tds = [oc.TestDataV(n, 700 + s, n=5) for s in range(4)]
        types = [3, 4, 2]
        E11 = oc.evaluate3(['box', 'box'], types, f, Ov, tds)[0]
        EA = [x + y for x, y in zip(oc.evaluate3([a, 'Q0'], types, f, Ov, tds)[0], oc.evaluate3(['Q0', a], types, f, Ov, tds)[0])]
        EB = [x + y for x, y in zip(oc.evaluate3([b, 'Q0'], types, f, Ov, tds)[0], oc.evaluate3(['Q0', b], types, f, Ov, tds)[0])]
        rk, rka, sol = solve([EA, EB], E11)
        print('trial %d: A edges %s spokes %s (I-sign %+d), B edges %s spokes %s (I-sign %+d); rank %d/%d, kappa = %s, '
              'rule/2 = (%s, %s)  [%.0fs]' % (trial, q3lib.ALL[a]['edges'], q3lib.ALL[a]['spokes'], oc.SIGNS[a],
                                             q3lib.ALL[b]['edges'], q3lib.ALL[b]['spokes'], oc.SIGNS[b], rk, rka, sol,
                                             Fr(rule(a), 2), Fr(rule(b), 2), time.time() - t0), flush=True)
else:
    t0 = time.time()
    names = [relabel(g, g + '_r') for g in ['G1', 'G2', 'G3', 'G4', 'G5']]
    for g in names:
        print(g, 'edges', q3lib.ALL[g]['edges'], 'spokes', q3lib.ALL[g]['spokes'], 'I-sign', oc.SIGNS[g], 'rule', rule(g))
    rngn = np.random.RandomState(17)
    Ov = [rngn.randint(-2, 3, size=n).astype(np.int64) for _ in range(4)]
    tds = [oc.TestDataV(n, 800 + s, n=6) for s in range(6)]
    types = [3, 4, 2, 2]
    cols = []
    for g in names:
        cols.append([x + y for x, y in zip(oc.evaluate3([g, 'Q0'], types, f, Ov, tds)[0], oc.evaluate3(['Q0', g], types, f, Ov, tds)[0])])
        print(g, 'done [%.0fs]' % (time.time() - t0), flush=True)
    Rhs = [Fr(0)] * len(tds)
    for a, b in [('box', 'A'), ('A', 'box'), ('box', 'B'), ('B', 'box')]:
        Rhs = [x + y for x, y in zip(Rhs, oc.evaluate3([a, b], types, f, Ov, tds)[0])]
    rk, rka, sol = solve(cols, Rhs)
    print('n = 3 relabelled: rank %d/%d, kappa = %s, rule/|Aut| = %s  [%.0fs]'
          % (rk, rka, sol, tuple(Fr(rule(g), oc.AUT[g]) for g in names), time.time() - t0))
