"""(d)-(f) Formal sl_N, weight (3,3), details:
 (1) theta' split by graph and by vertex class (trivalent N, N' / external vertices): closedness, exactness;
 (2) controls on the closedness test: kernels with the end sign of internal derivatives dropped (+y on both ends),
     single-graph sign flips; closedness of the results;
 (3) explicit primitives: theta' = d eta_th,  {D, iota} = d eta_br,  and i_3 = a eps pi nu xi with
     d xi = {D, iota} + theta'  (H^3_(3,3) = 0, C^2_(3,3) = 0: xi unique);
 (4) N-dependence of the primitives."""
import time
from common import *

t0 = time.time()
sl = True
deg = (3, 3)
H = load('h4w3_sl.pkl')
C3, B4, k4 = H['C3'], H['B4'], H['k4']
_, MB = la.matrix_from(B4, k4)
cB = la.certified_rank(MB)
KD = load('kernels.pkl')
AUT = {(1, 3, 3): 4, (2, 2, 3): 4}
D = J_cochain('d', sl)


def npow_parts(c):
    pw = sorted({p for d in c.values() for p in d})
    return {p: {kk: d[p] for kk, d in c.items() if d.get(p)} for p in pw}


def primitive(c):
    """eta (Npoly coefficients) with d eta = c, exact solve per N-power; None if not exact."""
    out = {}
    for p, x in npow_parts(c).items():
        ok, sol = la.in_span(MB, la.column(x, k4), cB)
        if not ok:
            return None
        for j, cj in enumerate(sol):
            if cj:
                for kk, val in C3[j].items():
                    ce.add_term(out, kk, {p: Fr(val)}, cj)
    return out


def closed(c):
    return is_zero(ce.d_ce(c, 4, sl))


def piece(pqr, vs, flip_end=False):
    G, ext, lam = graph_data(pqr)
    out = {}
    for v in vs:
        if flip_end:
            ins = insertion_nosign(G, ext, v)
            Kp = kernel_over_pi(G, lam, ins)
        else:
            Kp = KD[pqr]['Kv'][v]
        Wv = colour_weight(G.nv, G.edges, ext, {v: 'd'})
        out = ce.lin((1, out), (Fr(1, AUT[pqr]), alt_cochain(Kp, Wv, sl=True)))
    return out


def insertion_nosign(G, ext, v):
    local = local_orders(G.nv, G.edges, ext)[v]
    h1, h2 = local[0], local[1]

    def q(h, y, R):
        if h[0] == 'x':
            i = ext.index(h[1])
            return [R.from_expr(PSYM[i][0]), R.from_expr(PSYM[i][1])]
        return [y[h[1]][0], y[h[1]][1]]

    def ins(y):
        R = y[0][0].ring
        a, b = q(h1, y, R), q(h2, y, R)
        return a[0] * b[1] - a[1] * b[0]
    return ins


print("(1) pieces of theta' (sl_N): closed / exact")
pieces = {}
for pqr in AUT:
    for nm, vs in (("trivalent {0,5}", [0, 5]), ("external {1,2,3,4}", [1, 2, 3, 4]), ("all", range(6))):
        c = piece(pqr, vs)
        pieces[(pqr, nm)] = c
        cl = closed(c)
        ex = primitive(c) is not None if cl else None
        print("  theta%s, vertices %s: closed %s, exact %s" % (pqr, nm, cl, ex), flush=True)
th_tot = ce.lin((1, pieces[((1, 3, 3), "all")]), (1, pieces[((2, 2, 3), "all")]))
print("  sum == theta' of o3_theta:", is_zero(ce.lin((1, th_tot), (-1, H['thp']))))
print("(2) controls")
c = ce.lin((1, pieces[((1, 3, 3), "all")]), (-1, pieces[((2, 2, 3), "all")]))
print("  theta'_133 - theta'_223: closed %s" % closed(c))
c = ce.lin((1, piece((1, 3, 3), range(6), True)), (1, piece((2, 2, 3), range(6), True)))
print("  end sign dropped (+y on both ends): closed %s" % closed(c), flush=True)
c = ce.lin((1, pieces[((1, 3, 3), "trivalent {0,5}")]), (-1, pieces[((1, 3, 3), "external {1,2,3,4}")]),
           (1, pieces[((2, 2, 3), "trivalent {0,5}")]), (-1, pieces[((2, 2, 3), "external {1,2,3,4}")]))
print("  trivalent - external: closed %s" % closed(c), flush=True)

print("(3) primitives")
eta_th = primitive(H['thp'])
print("  theta' = d eta_th:", eta_th is not None and is_zero(ce.lin((1, ce.d_ce(eta_th, 3, sl)), (-1, H['thp']))))
print("  eta_th =", show(eta_th, 3))
eta_br = primitive(H['br'])
print("  {D, iota} = d eta_br:", eta_br is not None and is_zero(ce.lin((1, ce.d_ce(eta_br, 3, sl)), (-1, H['br']))))
print("  eta_br =", show(eta_br, 3))
xi = ce.lin((1, eta_th), (1, eta_br))
print("  xi = eta_th + eta_br, d xi == {D, iota} + theta':", is_zero(ce.lin((1, ce.d_ce(xi, 3, sl)), (-1, H['o3p']))))
print("  xi =", show(xi, 3))
# express in the 6-dim basis C^3_(3,3)
print("  C^3_(3,3) basis:")
for i, c in enumerate(C3):
    print("   c%d = %s" % (i, show(c, 3)))
save(dict(eta_th=eta_th, eta_br=eta_br, xi=xi), 'primitives_sl.pkl')
print("time %.1fs" % (time.time() - t0))
