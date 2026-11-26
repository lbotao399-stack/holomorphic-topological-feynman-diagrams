"""(d) Theta(J_1) piece by piece: for each graph G and J_1 vertex v,
    |Aut G|^{-1} alt(K_G^(v)/pi W_G^(v,d)) = d( x [p1p2]^3 d + y [p1p2]^3 s1 )   (gl_N; sl_N: y-term restricts to 0),
x, y polynomials in N (exact solve in the SL2-invariant, C-odd 3-cochains of weight (3,3))."""
from common import *

KD = load('kernels.pkl')
AUT = {(1, 3, 3): 4, (2, 2, 3): 4}
deg = (3, 3)
sl = False
k4 = spaces.keys_of(4, deg, sl)
u3 = ce.poly_mul(ce.poly_mul(ce.bracket_poly(0, 1, 3), ce.bracket_poly(0, 1, 3)), ce.bracket_poly(0, 1, 3))
Ts = {'d': {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}},
      's1': {((0,), (1, 2)): {0: Fr(1)}, ((0, 2), (1,)): {0: Fr(1)}, ((0, 1), (2,)): {0: Fr(1)}}}
Ps = [ce.canon_raw(ce.raw_from(u3, Ts[T]), 3, sl=sl) for T in ('d', 's1')]
dPs = [spaces.const_part(ce.d_ce(P_, 3, sl)) for P_ in Ps]
_, M = la.matrix_from(dPs, k4)
cc = la.certified_rank(M)
assert cc['rank'] == 2
Nn = sp.Symbol('N')


def prim(c):
    pw = sorted({p for d in c.values() for p in d})
    xs = [0, 0]
    for p in pw:
        x = {kk: Fr(d[p]) for kk, d in c.items() if d.get(p)}
        ok, sol = la.in_span(M, la.column(x, k4), cc)
        assert ok
        for j in range(2):
            xs[j] += sp.Rational(sol[j].numerator, sol[j].denominator) * Nn ** p
    return [sp.factor(x) for x in xs]


tot = [0, 0]
for pqr in AUT:
    G, ext, lam = graph_data(pqr)
    sub = [0, 0]
    for v in range(G.nv):
        Wv = colour_weight(G.nv, G.edges, ext, {v: 'd'})
        piece = ce.scale(alt_cochain(KD[pqr]['Kv'][v], Wv), Fr(1, AUT[pqr]))
        x = prim(piece)
        sub = [sub[j] + x[j] for j in range(2)]
        print("theta%s, J_1 at v=%d: d( %s [p1p2]^3 d + %s [p1p2]^3 s1 )" % (pqr, v, x[0], x[1]))
    print("theta%s total: d( %s [p1p2]^3 d + %s [p1p2]^3 s1 )" % (pqr, sp.factor(sub[0]), sp.factor(sub[1])))
    tot = [tot[j] + sub[j] for j in range(2)]
print("theta' total: d( %s [p1p2]^3 d + %s [p1p2]^3 s1 )" % (sp.factor(tot[0]), sp.factor(tot[1])))
