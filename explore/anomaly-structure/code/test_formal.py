"""formal_gl.py against explicit components: colour weights (gl2, gl3, sl2 via restriction) and
the CE differential (evaluation of formal d omega on gl3 vs explicit d of the explicit tensor)."""
import sys, itertools
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
from htk import theta, polygon
import lie, cochain as cc, formal_gl as fg

graphs = {'box': (polygon(4), [0, 1, 2, 3])}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, interior = theta(*pqr)
    graphs['theta%d%d%d' % pqr] = (G, [v for P in interior for v in P])
for name, (G, ext) in graphs.items():
    Wf = fg.colour_weight_formal(G.nv, G.edges, ext)
    print(name, 'formal colour weight:', {t: d for t, d in sorted(Wf.items())})
    for g in [lie.gl(2), lie.gl(3)]:
        We, sc = cc.colour_weight(g, G.nv, G.edges, ext)
        Wx = fg.tensor_from_formal(Wf, g, 4)
        print('   %s explicit == formal: %s' % (g.name, np.array_equal(We * 1, Wx * sc)))

# formal d vs explicit d on gl3, H = 2, weight (2,2): use theta133 with a test polynomial
mom = cc.Momenta(6, 2)
P = mom.p
G, ext = graphs['theta133']
Wf = fg.colour_weight_formal(G.nv, G.edges, ext)
Ptest = (P[0][0] * P[1][1] - P[0][1] * P[1][0]) * (P[2][0] * P[0][1] + 2 * P[1][0] * P[2][1])
wf = fg.alternate_formal(mom, Ptest, Wf, 4)
dwf = fg.d_formal(mom, wf, 4)
for g in [lie.gl(2), lie.gl(3)]:
    N = g.mats[0].shape[0]
    We, sc = cc.colour_weight(g, G.nv, G.edges, ext)
    we = cc.alternate(mom, Ptest, We, 4)
    dwe = cc.d_CE(mom, g, we, 4)
    # evaluate formal on basis (einsum)
    num = fg.at_N(dwf, N)
    by_m = {}
    for (tau, M), v in num.items():
        by_m.setdefault(M, {})[tau] = {0: v}
    ok = True
    for M in set(by_m) | set(dwe):
        T = fg.tensor_from_formal(by_m.get(M, {}), g, 5)
        E = dwe.get(M, np.zeros_like(T))
        ok &= np.array_equal(T * sc, E)
    print('formal d == explicit d on', g.name, ':', ok, flush=True)
