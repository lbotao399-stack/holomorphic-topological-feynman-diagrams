"""Task 4: elimination of a degree-4 vertex in the D = 2, n = 6 graphs without circuit chart (I-core K5).
S = the four edges at v (k = 4, m = 1, n_S = 0): response Lr = rank-structured mesh c_ab = g_a g_b / G,
transfer h_a = g_a / G, fiber dimension 0 (rho_S = 1).  Checks (exact):
  (i)  Omega_Gamma = eps Omega_out (mode R) at random points, sourced and unsourced v;
  (ii) mode P (dh terms dropped) for unsourced v (dh lambda_v = 0);
  (iii) Lr = g g^T / G - diag on the terminals, c_ab c_cd = c_ac c_bd."""
import json
import os
import random
from fractions import Fraction as Fr

import locality as L
from htk.graph import Graph
from check_locality import rand_lam

cen = json.load(open(os.path.join(L.HERE, '../atlas/census.json')))
gs = cen['data']['2,6']['graphs']
irr = [g for g in gs if not g['sp'] and not g['B'] and g['k_gen'] is None]
print('D=2 n=6 graphs without reduction:', len(irr))
rng = random.Random(4)
res = {}
for gi, g in enumerate(irr[:8]):
    E = [tuple(e) for e in g['edges']]
    nv = 1 + max(max(e) for e in E)
    G = Graph(nv, E)
    deg4 = [v for v in range(nv - 1) if G.degree(v) == 4]
    v = deg4[0]
    S = [e for e, (a, b) in enumerate(E) if v in (a, b)]
    sp_ = L.Split(G, S)
    for sourced in (True, False):
        for mode in ('R', 'P'):
            t = L.rand_times(G.ne, rng)
            lam = rand_lam(G, 2, rng)
            if not sourced:
                lam[v] = [0, 0]
            Db, Df, _ = L.factor_T0(sp_, t, lam, mode)
            ok = 'zero' if Db == 0 and Df == 0 else ('+1' if Db == Df else ('-1' if Db == -Df else 'FAIL'))
            key = ('v sourced' if sourced else 'v unsourced', mode)
            res.setdefault(key, []).append(ok)
    # rank structure of the response
    t = L.rand_times(G.ne, rng)
    K = sp_.kron(t)
    gv = {sp_.dS[a]: None for a in range(4)}
    gg = []
    for a in range(4):
        e = [e for e in S if sp_.dS[a] in G.edges[e]][0]
        gg.append(1 / Fr(t[e]))
    Gs = sum(gg)
    okr = all(-K['Lr'][a][b] == gg[a] * gg[b] / Gs for a in range(4) for b in range(4) if a != b)
    okh = all(K['E'][0][a] == gg[a] / Gs for a in range(4))
    res.setdefault('Lr = g g^T/G, h = g/G', []).append(okr and okh)
for k, v in res.items():
    print('  %-28s : %s' % (k, v))
