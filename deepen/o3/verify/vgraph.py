"""Graph data, Wick sign, colour factors (trace calculus) and J_1-vertex kernels for the o3 verification.
Own conventions (independent of the explorer's code):
  half-edge (e,0) at e(0) (tail), (e,1) at e(1) (head), ('x',v) external leg at v.
  local order (choice 'A'): external leg last, internal half-edges sorted by (edge, end) DESCENDING;
              (choice 'B'): external first, internal ascending (used for the covariance check).
  Wick sign sigma = sgn(source -> target), source = concatenation of local orders in vertex order,
              target = (e0,0),(e0,1),(e1,0),(e1,1),..., x_{v1},..,x_{v4} (externals in path order).
  Momentum into the vertex through a half-edge: (e,0): -y_e, (e,1): +y_e, external: lambda_v,
              lambda_v written through Kirchhoff, lambda_v = sum_{e(0)=v} y_e - sum_{e(1)=v} y_e.
"""
import sys
from fractions import Fraction as Fr
import sympy as sp
P_ROOT = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + ''
sys.path.insert(0, P_ROOT + '/code')
import htk
from vlib import LP, ONE, perm_from_words, has_fixed, seq_sign, add_to

PS = [[sp.Symbol('q%d_%d' % (i, c)) for c in range(2)] for i in range(4)]


def theta_data(pqr):
    G, interior = htk.theta(*pqr)
    ext = [v for P in interior for v in P]
    lam = {v: (list(PS[ext.index(v)]) if v in ext else [0, 0]) for v in range(G.nv - 1)}
    return G, ext, lam


def local_orders(G, ext, choice='A'):
    out = []
    for v in range(G.nv):
        hs = [('e', e, s) for e, ab in enumerate(G.edges) for s in (0, 1) if ab[s] == v]
        if choice == 'A':
            hs = sorted(hs, key=lambda h: (h[1], h[2]), reverse=True)
            if v in ext:
                hs = hs + [('x', v)]
        elif choice == 'B':
            hs = sorted(hs, key=lambda h: (h[1], h[2]))
            if v in ext:
                hs = [('x', v)] + hs
        elif choice == 'C':     # rotate A by one
            hs = sorted(hs, key=lambda h: (h[1], h[2]), reverse=True)
            if v in ext:
                hs = hs + [('x', v)]
            hs = hs[1:] + hs[:1]
        out.append(hs)
    return out


def wick_sign(G, ext, loc):
    src = [h for hs in loc for h in hs]
    tgt = []
    for e in range(len(G.edges)):
        tgt += [('e', e, 0), ('e', e, 1)]
    tgt += [('x', v) for v in ext]
    pos = {h: i for i, h in enumerate(tgt)}
    return seq_sign([pos[h] for h in src])


TENS = {'f': [(1, (0, 1, 2)), (-1, (0, 2, 1))],
        'd': [(1, (0, 1, 2)), (1, (0, 2, 1))],
        's1': [(1, ((0,), (1, 2))), (1, ((1,), (0, 2))), (1, ((2,), (0, 1)))],
        's0': [(1, ((0,), (1,), (2,)))]}


def _tens_words(T, hs):
    out = []
    for c, w in TENS[T]:
        if isinstance(w[0], tuple):
            out.append((c, [[hs[i] for i in ww] for ww in w]))
        else:
            out.append((c, [[hs[i] for i in w]]))
    return out


def colour(G, ext, loc, vtypes, sl):
    """sum over contractions of prod_v T_v(local order) prod_e kappa^{-1}; returns {perm on ext slots: LP}."""
    terms = [(ONE, [])]
    for v in range(G.nv):
        T = vtypes.get(v, 'f')
        new = []
        for c, ws in terms:
            for cc, words in _tens_words(T, loc[v]):
                new.append((c * cc, ws + [list(w) for w in words]))
        terms = new
    Ninv = LP.N(-1)
    NN = LP.N(1)
    for e in range(len(G.edges)):
        h, hp = ('e', e, 0), ('e', e, 1)
        new = []
        for c, ws in terms:
            i = next(k for k, w in enumerate(ws) if h in w)
            j = next(k for k, w in enumerate(ws) if hp in w)
            if i != j:
                w1 = ws[i]; a = w1.index(h); A = w1[a + 1:] + w1[:a]      # tr(X_e A)
                w2 = ws[j]; b = w2.index(hp); B = w2[b + 1:] + w2[:b]     # tr(X^e B)
                rest = [w for k, w in enumerate(ws) if k not in (i, j)]
                # gl: tr(A B)
                if A + B:
                    new.append((c, rest + [A + B]))
                else:
                    new.append((c * NN, rest))
                if sl:
                    cc = -c * Ninv
                    r2 = list(rest)
                    for part in (A, B):
                        if part:
                            r2.append(part)
                        else:
                            cc = cc * NN
                    new.append((cc, r2))
            else:
                w = ws[i]; a = w.index(h); w = w[a:] + w[:a]       # (X_e, B..., X^e, A...)
                b = w.index(hp)
                B = w[1:b]; A = w[b + 1:]
                rest = [x for k, x in enumerate(ws) if k != i]
                cc, r2 = c, list(rest)
                for part in (A, B):
                    if part:
                        r2.append(part)
                    else:
                        cc = cc * NN
                new.append((cc, r2))
                if sl:
                    AB = A + B
                    if AB:
                        new.append((-c * Ninv, rest + [AB]))
                    else:
                        new.append((-c * Ninv * NN, rest))
        terms = new
    out = {}
    for c, ws in terms:
        words = [[ext.index(h[1]) for h in w] for w in ws]
        perm = perm_from_words(words, len(ext))
        if sl and has_fixed(perm):
            continue
        add_to(out, perm, c)
    return out


def q_of(h, y, G):
    if h[0] == 'x':
        v = h[1]
        q = [0, 0]
        for e, (a, b) in enumerate(G.edges):
            for i in range(2):
                if a == v:
                    q[i] = q[i] + y[e][i]
                if b == v:
                    q[i] = q[i] - y[e][i]
        return q
    e, s = h[1], h[2]
    sg = -1 if s == 0 else 1
    return [sg * y[e][0], sg * y[e][1]]


def insertion(G, loc, v, pair=(0, 1)):
    h1, h2 = loc[v][pair[0]], loc[v][pair[1]]

    def ins(y):
        a, b = q_of(h1, y, G), q_of(h2, y, G)
        return a[0] * b[1] - a[1] * b[0]
    return ins


_red = {}


def kernel_over_pi(G, lam, ins=None, check=2):
    key = (tuple(G.edges), G.nv)
    if key not in _red:
        _red[key] = htk.sp_reduce(G, prefer="large")
    K = sp.expand(_red[key].integrate(lam, 1, insertion=ins, check=check))
    Kp = sp.expand(K / sp.pi)
    assert not Kp.has(sp.pi), K
    return Kp


def poly_to_full(Kp, k=4):
    """sympy polynomial in PS -> full-momentum dict {(a_0..a_3, b_0..b_3): Fr}."""
    gens = [PS[i][0] for i in range(k)] + [PS[i][1] for i in range(k)]
    if Kp == 0:
        return {}
    out = {}
    for mon, c in sp.Poly(Kp, *gens).terms():
        c = sp.Rational(c)
        out[tuple(mon)] = Fr(int(c.p), int(c.q))
    return out
