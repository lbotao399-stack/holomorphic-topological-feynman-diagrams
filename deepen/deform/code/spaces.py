"""Spaces of invariant cochains (trace basis, formal N), obstruction cochains, exact bases."""
import sys
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import ce, la


def taus(k, sl):
    t = ce.all_taus(k)
    return [x for x in t if not ce.has_singleton(x)] if sl else t


def spanning(k, deg, sl):
    """alt(tau (x) mu) over taus and canonical monomials of bidegree deg (N-independent)."""
    out = []
    for tau in taus(k, sl):
        for mu in ce.monos_canonical(k, deg):
            c = ce.alt_raw({(tau, ce.full_of_canon(mu)): {0: Fr(1)}}, k, sl=sl)
            if c:
                out.append(c)
    return out


def const_part(c):
    """N-independent cochain -> dict key -> Fraction (asserts N-independence)."""
    out = {}
    for k, d in c.items():
        assert set(d) <= {0}, d
        out[k] = d[0]
    return out


def keys_of(k, deg, sl):
    return [(t, m) for t in taus(k, sl) for m in ce.monos_canonical(k, deg)]


def basis(k, deg, sl):
    """exact basis of C^k_inv(deg): list of cochains (N-independent), certified."""
    sp = [const_part(c) for c in spanning(k, deg, sl)]
    keys = keys_of(k, deg, sl)
    if not sp:
        return [], keys
    _, M = la.matrix_from(sp, keys)
    cert = la.certified_rank(M)
    B = [sp[j] for j in cert['pivcols']]
    return B, keys


def to_np(c):
    return {k: {0: Fr(v)} for k, v in c.items()}


def obstruction_cochains():
    """F_133, F_223 (formal gl_N, Npoly coefficients) and the engine kernels (60/pi) K."""
    import htk
    import formal_gl as fg
    from c_obstruction import kernel_poly
    import sympy as sp_
    out = {}
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        G, ext, K = kernel_poly(pqr)
        W = fg.colour_weight_formal(G.nv, G.edges, ext)
        Wm = {tau: {p: Fr(v) for p, v in d.items()} for tau, d in W.items()}
        poly = ce.poly_mul(ce.bracket_poly(0, 1, 4), ce.bracket_poly(2, 3, 4))
        F = ce.alt_raw(ce.raw_from(poly, Wm), 4)
        # kernel check: (60/pi) K == c [p0 p1][p2 p3]
        P = [[sp_.Symbol('p%d_%d' % (i, c)) for c in range(2)] for i in range(6)]
        br = lambda a, b: a[0] * b[1] - a[1] * b[0]
        Kc = sp_.simplify(K / (br(P[0], P[1]) * br(P[2], P[3])))
        out[pqr] = dict(W=Wm, F=F, Kcoef=Kc)
    return out
