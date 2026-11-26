"""Common definitions for task o3 (order hbar^3 of commutative CS on R x C^2, (H,T) = (2,1)).

Graphs: htk.theta(p,q,r) engine labelling (N = 0, N' = root = 5, ext = interior vertices in path order).
Local order L'_v (normalization b5/b12, flagship Prop. 3.5): [x_v if v external] + incident half-edges
sorted by (edge, end).  Half-edge (e,s): s = 0 tail e(0), s = 1 head e(1).
Momentum flowing into a vertex through a half-edge (holomorphic derivative acting on that leg):
    external x_v: lambda_v (= p of the external leg);  (e,0): -y_e;  (e,1): +y_e        (o1 check)
Vertex tensors in local order (words = cycles of traces over the local positions 0,1,2):
    f = tr(012) - tr(021),  d = tr(012) + tr(021),  s1 = tr(0)tr(12) + tr(1)tr(02) + tr(2)tr(01),  s0 = tr0 tr1 tr2.
Colour weight W = sgn(L -> L') prod_v T_v(local order) prod_e kappa^{-1}, kappa = tr on gl_N (double line):
    tr(A X_e) tr(B X^e) = tr(AB),  tr(X_e A X^e B) = tr A tr B,  tr 1 = N.
"""
import sys, itertools, pickle, os
from fractions import Fraction as Fr
import sympy as sp

P_ROOT = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + ''
sys.path.insert(0, P_ROOT + '/code')
sys.path.insert(0, P_ROOT + '/deepen/deform/code')
sys.path.insert(0, P_ROOT + '/explore/anomaly-structure/code')
import htk
from htk import theta
import ce, la, spaces

OUT = P_ROOT + '/deepen/o3'
DATA = OUT + '/code/data'
os.makedirs(DATA, exist_ok=True)

PSYM = [[sp.Symbol('p%d_%d' % (i, c)) for c in range(2)] for i in range(4)]

from ribbon import TENS, canon_cycle, canon_tau, colour_weight, local_orders, orientation_sign, perm_sign


def graph_data(pqr):
    G, interior = theta(*pqr)
    ext = [v for Pp in interior for v in Pp]
    lam = {v: (list(PSYM[ext.index(v)]) if v in ext else [0, 0]) for v in range(G.nv - 1)}
    return G, ext, lam


def insertion_bracket(G, ext, v, which=(0, 1), lamnum=None):
    """callable y -> [q_h q_h'] (h, h' = local positions `which` at v) as a ring element.
    lamnum: numeric sources {vertex: [l1, l2]} (default: symbolic PSYM)."""
    local = local_orders(G.nv, G.edges, ext)[v]
    h1, h2 = local[which[0]], local[which[1]]

    def q(h, y, R):
        if h[0] == 'x':
            if lamnum is not None:
                return [R.ground_new(lamnum[h[1]][0]), R.ground_new(lamnum[h[1]][1])]
            i = ext.index(h[1])
            return [R.from_expr(PSYM[i][0]), R.from_expr(PSYM[i][1])]
        e, s = h[1], h[2]
        sg = -1 if s == 0 else 1
        return [sg * y[e][0], sg * y[e][1]]

    def ins(y):
        R = y[0][0].ring
        a, b = q(h1, y, R), q(h2, y, R)
        return a[0] * b[1] - a[1] * b[0]
    return ins


def kernel_over_pi(G, lam, insertion=None):
    """K / pi as a sympy polynomial (rational coefficients) in the p's; asserts the pi^1 form."""
    red = htk.sp_reduce(G, prefer="large")          # as htk.kernel_exact for T = 1
    K = sp.expand(red.integrate(lam, 1, insertion=insertion, check=2))
    Kp = sp.expand(K / sp.pi)
    assert not Kp.has(sp.pi), K
    return Kp


def poly_to_raw(Kp):
    """sympy polynomial in PSYM -> dict fullmono (slot-major over 4 slots) -> Fraction."""
    gens = [PSYM[i][c] for i in range(4) for c in range(2)]
    if Kp == 0:
        return {}
    Pl = sp.Poly(Kp, *gens)
    out = {}
    for mon, c in Pl.terms():
        out[tuple(mon)] = Fr(int(sp.Rational(c).p), int(sp.Rational(c).q))
    return out


def alt_cochain(Kp, W, sl=False):
    return ce.alt_raw(ce.raw_from(poly_to_raw(Kp), W), 4, sl=sl)


def J_cochain(name, sl=False):
    """[p1p2] T for T in {d, s1, s0} (canonical, alternating; = D, S1, S0 of deform)."""
    u = ce.bracket_poly(0, 1, 3)
    Ts = {'d': {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}},
          's1': {((0,), (1, 2)): {0: Fr(1)}, ((0, 2), (1,)): {0: Fr(1)}, ((0, 1), (2,)): {0: Fr(1)}},
          's0': {((0,), (1,), (2,)): {0: Fr(1)}}}
    return ce.canon_raw(ce.raw_from(u, Ts[name]), 3, sl=sl)


def f_bracket_sq(sl=False):
    """[p1p2]^2 f, f = tr(123) - tr(132) (= [q0q1]^2 f), canonical 3-cochain of weight (2,2)."""
    u = ce.bracket_poly(0, 1, 3)
    u2 = ce.poly_mul(u, u)
    return ce.canon_raw(ce.raw_from(u2, {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(-1)}}), 3, sl=sl)


def save(obj, name):
    pickle.dump(obj, open(DATA + '/' + name, 'wb'))


def load(name):
    return pickle.load(open(DATA + '/' + name, 'rb'))


def is_zero(c):
    return ce.clean(c) == {}


def show(c, k, maxterms=None):
    terms = []
    items = sorted(c.items(), key=lambda x: (x[0][0], x[0][1]))
    for (t, m), v in items:
        mon = "".join(("p%d%s" % (s + 1, "ab"[cc])) + ("^%d" % m[2 * s + cc] if m[2 * s + cc] > 1 else "")
                      for s in range(k - 1) for cc in range(2) if m[2 * s + cc])
        vv = {p: str(x) for p, x in v.items()} if isinstance(v, dict) else str(v)
        terms.append("%s %s %s" % (vv, mon or "1", ce.tau_str(t)))
    if maxterms and len(terms) > maxterms:
        terms = terms[:maxterms] + ["... (%d terms)" % len(items)]
    return " + ".join(terms)


def inv_span(k, sl, deg):
    """SL2-invariant spanning set at bidegree (w,w): alt(B tau), B = products of w brackets [p_i p_j] of the k
    momenta (raw, slot-major), tau over trace monomials; returns N-independent canonical cochains."""
    assert deg[0] == deg[1]
    w = deg[0]
    pairs = list(itertools.combinations(range(k), 2))
    out = []
    for tup in itertools.combinations_with_replacement(pairs, w):
        poly = {tuple([0] * (2 * k)): 1}
        for (i, j) in tup:
            poly = ce.poly_mul(poly, ce.bracket_poly(i, j, k))
        for tau in spaces.taus(k, sl):
            c = ce.clean(ce.alt_raw({(tau, m): {0: Fr(v)} for m, v in poly.items()}, k, sl=sl))
            if c:
                out.append(spaces.const_part(c))
    return out


def rev_tau(t):
    return ce.canon_tau([tuple(reversed(c)) for c in t])


def Cop(w, k):
    """charge conjugation X -> -X^T on a k-cochain (trace basis)."""
    s = (-1) ** k
    out = {}
    for (t, m), d in w.items():
        ce.add_term(out, (rev_tau(t), m), d if isinstance(d, dict) else {0: Fr(d)}, s)
    return out
