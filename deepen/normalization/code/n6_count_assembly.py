"""(d) assembly arithmetic of O_2, cross-check of the msym-odd kernels with htk, (e) counting for D = 2, 3.

(e) Anomaly graphs (one heat edge, cubic vertices, external legs = fields A of total degree 1), tight:
    |I| = Dn + 1, |V| = (D-1)n + 2.  External legs E = 3|V| - 2|I|.  The edge forms fill the relative positions,
    so the external block must supply a D-form (tau, zb directions) at the common point: ghost number gh = E - D.
    Independent parity count of the graph functional: |V| D (vertex measures Omega_v^H with int over 2H+T dims)
    + (D-1)(|I|-1) (P) + D (K) + E (fields), mod 2.
"""
import sys
import itertools
import sympy as sp
import networkx as nx
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import catalog, Graph

print("(e) counting")
for D in (2, 3):
    for n in range(1, 7):
        I_, V_ = D * n + 1, (D - 1) * n + 2
        E = 3 * V_ - 2 * I_
        gh = E - D
        par = (V_ * D + (D - 1) * (I_ - 1) + D + E) % 2
        print("  D=%d (H,T)=(%d,1) n=%d: |I|=%2d |V|=%2d  E=3|V|-2|I|=%2d  holomorphic derivatives %d (n per component)  "
              "gh=E-D=%2d  parity(count)=%d  T=1 kernel %s"
              % (D, D - 1, n, I_, V_, E, (D - 1) * n, gh, par, "zero (rho_n = 0, n odd)" if n % 2 else "nonzero possible"))
print("  formula: E = (D-3)n + 4, gh = (D-3)(n-1) + 1")
print("  tight graphs with all degrees <= 3 (htk.catalog, subgraph bound), external legs = sum(3 - deg):")
for D, n in [(2, 1), (2, 2), (3, 1), (3, 2)]:
    gs = catalog.enumerate_tight(D, n)
    cub = [G for G in gs if max(dict(G.degree()).values()) <= 3]
    desc = []
    for G in cub:
        degs = sorted(dict(G.degree()).values())
        desc.append((degs, sum(3 - d for d in degs)))
    print("    D=%d n=%d: %d tight graphs, %d with max degree <= 3: %s" % (D, n, len(gs), len(cub), desc))

print("(d) msym-odd kernels cross-checked with htk ((H,T) = (2,1))")
p = {k: [sp.Symbol('%s_1' % k), sp.Symbol('%s_2' % k)] for k in ['A1', 'B1', 'B2', 'C1', 'C2', 'l1', 'l2', 'l4', 'l5', 'u', 'w', 'a', 'b']}


def br(x, y):
    return x[0] * y[1] - x[1] * y[0]


def rootN(pp, q, r):
    nv = pp + q + r - 1
    N_, Np = nv - 1, nv - 2
    edges, nxt, paths = [], 0, []
    for n_ in (pp, q, r):
        prev, ints = N_, []
        for j in range(n_):
            cur = Np if j == n_ - 1 else nxt
            if j < n_ - 1:
                ints.append(nxt)
                nxt += 1
            edges.append((prev, cur))
            prev = cur
        paths.append(ints)
    return Graph(nv, edges, "rootN"), paths


G, paths = rootN(1, 3, 3)
lam = {v: [0, 0] for v in range(G.nv - 1)}
for v, nm in zip(paths[1] + paths[2], ['B1', 'B2', 'C1', 'C2']):
    lam[v] = p[nm]
K = htk.kernel_exact(G, lam, T=1)
print("  rootN(1,3,3): K_htk - (pi/20)[B1B2][C1C2] =", sp.simplify(K - sp.pi / 20 * br(p['B1'], p['B2']) * br(p['C1'], p['C2'])))
G, paths = rootN(2, 2, 3)
lam = {v: [0, 0] for v in range(G.nv - 1)}
for v, nm in zip(paths[0] + paths[1] + paths[2], ['A1', 'B1', 'C1', 'C2']):
    lam[v] = p[nm]
K = htk.kernel_exact(G, lam, T=1)
print("  rootN(2,2,3): K_htk + (pi/30)[A1B1][C1C2] =", sp.simplify(K + sp.pi / 30 * br(p['A1'], p['B1']) * br(p['C1'], p['C2'])))
GA = htk.graph_A()
lam = {0: p['l1'], 1: p['l2'], 2: [0, 0], 3: p['l4'], 4: p['l5']}
KA = htk.kernel_exact(GA, lam, T=1)
GB = htk.graph_B()
lam = {0: p['u'], 1: p['w'], 2: [0, 0], 3: p['a'], 4: p['b']}
KB = htk.kernel_exact(GB, lam, T=1)
print("  graph_A: K_htk = %s ; graph_B: K_htk = %s" % (sp.factor(KA), sp.factor(KB)))
print("    msym-odd theta-kernels claim: K_A = (pi/20)(-<12><45>), K_B = -(pi/30)(-<uw><ab>):",
      sp.simplify(KA + sp.pi / 20 * br(p['l1'], p['l2']) * br(p['l4'], p['l5'])) == 0,
      sp.simplify(KB - sp.pi / 30 * br(p['u'], p['w']) * br(p['a'], p['b'])) == 0)

print("(d) assembly of O_2 (exact arithmetic)")
pi = sp.pi
hb = sp.Symbol('hbar', positive=True)
U = sp.Symbol('U')            # universal sign (convention-dependent), U^2 = 1
vs = -1                       # vs_ext * s1 * s2 at n = 2 (n3_extraction.log)
Nnorm = U * vs * hb ** 2 / (64 * pi ** 5)
K133, K223 = -pi / 20, pi / 30
Aut = 4
c = sp.simplify(Nnorm / Aut * pi / 60)
print("  raw G_Gamma = vs (8 pi^{5/2})^{-2} K_Gamma = %s K_Gamma" % sp.simplify(vs / (8 * pi ** sp.Rational(5, 2)) ** 2))
print("  Nnorm = U * vs * hbar^2/(64 pi^5) =", Nnorm)
print("  O_2 = Nnorm/4 (K_133 W_133 + K_223 W_223), K_133 = -3 pi/60 [..], K_223 = 2 pi/60 [..]  ->  O_2 = c int [d1d2][d3d4] (-3 W_133 + 2 W_223),"
      " c = pi Nnorm/240 =", c, "; check pi*Nnorm/240 - c =", sp.simplify(pi * Nnorm / 240 - c))
# class: F_133 = -F_223 mod exact (anomaly-structure c7, confirmed)
F133, F223 = sp.symbols('F133 F223')
cls = sp.expand((c * (-3 * F133 + 2 * F223)).subs(F133, -F223))
print("  [O_2] = c(-3[F_133] + 2[F_223]) with [F_133] = -[F_223]:", sp.factor(cls), " = (pi Nnorm/48)[F_223]:",
      sp.simplify(cls - pi * Nnorm / 48 * F223) == 0)
cls_opp = sp.expand((c * (-3 * F133 - 2 * F223)).subs(F133, -F223))
print("  with the opposite relative sign: ", sp.factor(cls_opp), "= (pi Nnorm/240)[F_223]:",
      sp.simplify(cls_opp - pi * Nnorm / 240 * F223) == 0)
print("  |c| = %s, |c'| = %s" % (sp.simplify(c.subs(U, 1) * vs), sp.simplify(pi * Nnorm.subs(U, 1) * vs / 48)))
