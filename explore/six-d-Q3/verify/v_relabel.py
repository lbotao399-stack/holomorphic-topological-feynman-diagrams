"""C9: the sign rule sigma_Gamma = sgn(pi) (-1)^{#spokes O->v}.
(a) sigma on the reference labellings (box, A, B, G1..G5).
(b) labelling independence tested at operator level, without nilpotency: for random relabellings
    Gamma' (vertex permutation, edge reversals, edge-order permutation) the kernel ratio
    eps = I_Gamma'/I_Gamma is computed by brute force (own Kirchhoff layer, same Schwinger point),
    the operator hat Q^I_Gamma' is built from the relabelled leg pattern (T, output order, delta order),
    and sigma_Gamma' <Phi, hat Q^I_Gamma' O> = sigma_Gamma <Phi, hat Q^I_Gamma O> is checked exactly on
    test operators (with and without an input derivative)."""
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import vk
import vops
from vops import Op, pair, so_n, build_spec, GRAPH_SPECS

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
NREL = int(sys.argv[2]) if len(sys.argv) > 2 else 4
dim, f = so_n(5)


def sigma_rule(edges):
    verts = sorted({x for e in edges for x in e if x != 'O'})
    deg = {v: sum(v in e for e in edges) for v in verts}
    L1 = []
    for v in verts:
        ins = [i for i, e in enumerate(edges) if e[1] == v]
        outs = [i for i, e in enumerate(edges) if e[0] == v]
        L1 += [('e', i, 'h') for i in ins] + ([('out', v)] if deg[v] == 2 else []) + [('e', i, 't') for i in outs]
    spokes = [i for i, e in enumerate(edges) if 'O' in e]
    L2 = []
    for i, e in enumerate(edges):
        if i not in spokes:
            L2 += [('e', i, 't'), ('e', i, 'h')]
    for i in spokes:
        L2.append(('e', i, 'h') if edges[i][0] == 'O' else ('e', i, 't'))
    L2 += [('out', v) for v in verts if deg[v] == 2]
    assert sorted(L1) == sorted(L2)
    p = [L1.index(x) for x in L2]
    s = vk.perm_sign(p)
    nout = sum(1 for i in spokes if edges[i][0] == 'O')
    return s * (-1) ** nout


def density_t(edges, verts, t, dt, lam, tree):
    net = vk.Network(edges, verts)
    y, D = net.dy_dt(t, lam)
    sT, cot = net.tree_sign(tree)
    k = len(dt[0])
    m = len(t)
    rows = [[sum(D[c][a][j] * dt[j][kk] for j in range(m)) for kk in range(k)] for c in cot for a in range(3)]
    orient = vk.det([t] + [[dt[j][kk] for j in range(m)] for kk in range(k)])
    return sT * vk.det(rows) * (1 if orient > 0 else -1)


def relabel(name):
    base = GRAPH_SPECS[name]
    E = base['edges']
    verts = base['verts']
    newlab = list(range(1, len(verts) + 1))
    rng.shuffle(newlab)
    vmap = dict(zip(verts, newlab))
    vmap['O'] = 'O'
    order = list(range(len(E)))
    rng.shuffle(order)                       # new edge k = old edge order[k]
    rev = [rng.random() < 0.5 for _ in E]
    newE = []
    for k in order:
        a, b = E[k]
        a, b = vmap[a], vmap[b]
        newE.append((b, a) if rev[k] else (a, b))
    old_spokes = base['spokes']
    new_spokes = [i for i, e in enumerate(newE) if 'O' in e]
    smap = [old_spokes.index(order[i]) for i in new_spokes]
    return newE, vmap, order, smap


def eps_bruteforce(name, newE, vmap, order):
    g = vk.GRAPHS[name]
    verts = g['verts']
    nverts = sorted(vmap[v] for v in verts)
    vals = []
    while len(vals) < 3:
        r = vk.random_point(name, rng)
        t, dt = vk.chart(name, r)
        lam = {v: tuple(Fr(rng.randint(-5, 5)) for _ in range(3)) for v in verts}
        tree = vk.Network(g['edges'], verts).random_tree(rng)
        d0 = density_t(g['edges'], verts, t, dt, lam, tree)
        t2 = [t[order[k]] for k in range(len(order))]
        dt2 = [dt[order[k]] for k in range(len(order))]
        lam2 = {vmap[v]: lam[v] for v in verts}
        net2 = vk.Network(newE, nverts)
        tree2 = net2.random_tree(rng)
        d1 = density_t(newE, nverts, t2, dt2, lam2, tree2)
        if d0 == 0:
            assert d1 == 0
            continue
        vals.append(d1 / d0)
    assert len(set(vals)) == 1 and abs(vals[0]) == 1, vals
    return int(vals[0])


TESTS = {1: (['g1', 'g2'], ['c'] * 3), 2: (['g1', 'g2', 'b'], ['c'] * 4), 3: (['g1', 'g2', 'b', 'b'], ['c'] * 5)}
LOOPS = {'box': 1, 'A': 2, 'B': 2, 'G1': 3, 'G2': 3, 'G3': 3, 'G4': 3, 'G5': 3}

t00 = time.time()
print('(a) sigma on the reference labellings:', {nm: sigma_rule(GRAPH_SPECS[nm]['edges']) for nm in LOOPS})
cnt = 0
okall = True
neq = nnz = 0
for name in LOOPS:
    Ot, Ft = TESTS[LOOPS[name]]
    for k in range(NREL):
        newE, vmap, order, smap = relabel(name)
        eps = eps_bruteforce(name, newE, vmap, order)
        nm2 = '%s_r%d' % (name, k)
        GRAPH_SPECS[nm2] = build_spec(dict(edges=newE, doms=GRAPH_SPECS[name]['doms'], base=name,
                                           vmap={v: vmap[v] for v in GRAPH_SPECS[name]['verts']}, smap=smap, eps=eps))
        s0, s1 = sigma_rule(GRAPH_SPECS[name]['edges']), sigma_rule(newE)
        O = [(('O', i), t) for i, t in enumerate(Ot)]
        cols = {('O', i): [rng.randint(-2, 2) for _ in range(dim)] for i in range(len(Ot))}
        res = []
        for jets in ({}, {('O', 0): [rng.randrange(3)]}):
            F = (Ft, np.array([[rng.randint(-2, 2) for _ in range(dim)] for _ in Ft], dtype=float),
                 [[rng.randint(-4, 4) for _ in range(3)] for _ in Ft])
            v0 = pair([Op(name, 'd')], O, cols, jets, F, dim, f, {})
            v1 = pair([Op(nm2, 'd')], O, cols, jets, F, dim, f, {})
            ok = (s1 * v1 == s0 * v0)
            okall &= ok
            neq += 1
            nnz += (v0 != 0)
            res.append('jets %s: <Q_G O> = %s, <Q_G\' O> = %s, sigma\' v\' == sigma v: %s' % (jets, v0, v1, ok))
        cnt += 1
        print('%s relabelled %d: edges %s, eps = I\'/I = %+d, sigma = %+d, sigma\' = %+d; %s'
              % (name, k, newE, eps, s0, s1, '; '.join(res)), flush=True)
print('%d relabellings, %d comparisons, sigma\' v\' == sigma v in all: %s, comparisons with v != 0: %d  [%.0fs]'
      % (cnt, neq, okall, nnz, time.time() - t00))
