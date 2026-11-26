"""Colour weights of oriented graphs and Chevalley-Eilenberg cochains of g[[z_1..z_H]] modulo
translations (momentum conservation).

Colour weight (orientation rule).  Graph: vertices 0..nv-1 in label order, oriented edges, external
legs at the vertices ext = [v_1, ..., v_k] (one leg each).  Half-edges:
    L  = (x_{v_1}, ..., x_{v_k}, (e_0,out), (e_0,in), (e_1,out), (e_1,in), ...)
    L' = concatenation over v = 0..nv-1 of [x_v (if v external)] + [half-edges of the edges at v,
         sorted by (edge, end)],
    W[a_1..a_k] = sgn(L -> L') * sum prod_v fl[half-edges at v in L'-order] prod_e kinv[(e,out),(e,in)],
external index a_i on x_{v_i}.  W changes sign under a transposition of two vertex labels and under
the reversal of an edge, as K_Gamma does for D = 3; K_Gamma W is independent of these choices.

Cochains.  An alternating k-cochain omega(X_1 e^{p_1 z}, ..., X_k e^{p_k z}) =
    sum_m  m(p_1, ..., p_{k-1}) T_m[a_1..a_k]  X_1^{a_1} ... X_k^{a_k}
in canonical form: p_k = -(p_1 + ... + p_{k-1}) substituted, m monomials, T_m integer tensors.
CE differential:
    (d omega)(x_0..x_k) = sum_{i<j} (-1)^{i+j} omega([x_i, x_j], x_0, .., ^i, .., ^j, .., x_k),
    [X e^{pz}, Y e^{qz}] = [X, Y] e^{(p+q)z}.
"""
import itertools
import numpy as np
import sympy as sp

LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def colour_weight(lie, nv, edges, ext):
    L = [('x', v) for v in ext]
    for e in range(len(edges)):
        L += [('e', e, 0), ('e', e, 1)]
    local = []
    for v in range(nv):
        hs = [('x', v)] if v in ext else []
        hs += sorted([('e', e, s) for e, ab in enumerate(edges) for s in (0, 1) if ab[s] == v],
                     key=lambda h: (h[1], h[2]))
        local.append(hs)
    Lp = [h for hs in local for h in hs]
    assert sorted(Lp, key=str) == sorted(L, key=str)
    pos = {h: i for i, h in enumerate(L)}
    sgn = perm_sign([pos[h] for h in Lp])
    # letters
    let = {h: LETTERS[i] for i, h in enumerate(L)}
    ops, terms = [], []
    for v in range(nv):
        hs = local[v]
        if len(hs) != 3:
            raise ValueError("vertex %d is not trivalent (%d half-edges)" % (v, len(hs)))
        terms.append(''.join(let[h] for h in hs))
        ops.append(lie.fl)
    for e in range(len(edges)):
        terms.append(let[('e', e, 0)] + let[('e', e, 1)])
        ops.append(lie.kinv)
    out = ''.join(let[('x', v)] for v in ext)
    W = np.einsum(','.join(terms) + '->' + out, *ops, optimize='greedy')
    return sgn * W, lie.kinv_scale ** len(edges)


# ------------------------------------------------------------------ polynomials
class Momenta:
    def __init__(self, kmax, H):
        self.H = H
        self.p = [[sp.Symbol('p%d_%d' % (i, c)) for c in range(H)] for i in range(kmax)]

    def gens(self, k):
        return [self.p[i][c] for i in range(k - 1) for c in range(self.H)]

    def canon(self, expr, k):
        """dict monomial -> rational, after p_{k-1} = -(p_0 + ... + p_{k-2})."""
        sub = {self.p[k - 1][c]: -sum(self.p[i][c] for i in range(k - 1)) for c in range(self.H)}
        e = sp.expand(sp.sympify(expr).subs(sub, simultaneous=True))
        if e == 0:
            return {}
        P = sp.Poly(e, *self.gens(k))
        return dict(P.terms())

    def mono_expr(self, m, k, args):
        """monomial m (exponents over gens(k)) evaluated at the k-1 momentum vectors args."""
        out = sp.Integer(1)
        idx = 0
        for i in range(k - 1):
            for c in range(self.H):
                if m[idx]:
                    out *= args[i][c] ** m[idx]
                idx += 1
        return out


def transpose_for(sigma):
    """axes such that np.transpose(W, axes)[a_1..a_k] = W[a_{sigma(1)}, .., a_{sigma(k)}]."""
    return list(np.argsort(sigma))


def alternate(mom, P, W, k):
    """canonical form of sum_sigma sgn(sigma) P(p_sigma(1..k)) W[a_sigma(1..k)].
    P: sympy expression in mom.p[0..k-1]; W: integer numpy tensor of rank k."""
    out = {}
    for sigma in itertools.permutations(range(k)):
        s = perm_sign(sigma)
        sub = {mom.p[i][c]: mom.p[sigma[i]][c] for i in range(k) for c in range(mom.H)}
        Ps = sp.sympify(P).subs(sub, simultaneous=True)
        Ws = np.transpose(W, transpose_for(sigma))
        for m, coef in mom.canon(Ps, k).items():
            q = sp.Rational(coef)
            assert q.q == 1, "non-integer coefficient; rescale P"
            out[m] = out.get(m, 0) + s * int(q.p) * Ws
    return {m: T for m, T in out.items() if T.any()}


def d_CE(mom, lie, omega, k):
    """omega: canonical k-cochain dict; returns canonical (k+1)-cochain dict."""
    out = {}
    P = mom.p
    for i, j in itertools.combinations(range(k + 1), 2):
        rest = [r for r in range(k + 1) if r not in (i, j)]
        args = [[P[i][c] + P[j][c] for c in range(mom.H)]] + [[P[r][c] for c in range(mom.H)] for r in rest]
        sgn = (-1) ** (i + j)
        for m, T in omega.items():
            # colour: C[a_i, a_j, a_rest...] = sum_e fu[a_i, a_j, e] T[e, rest]
            C = np.tensordot(lie.fu, T, axes=([2], [0]))           # axes (a_i, a_j, r_1..r_{k-1})
            # reorder to (a_0..a_k): current axis order is [i, j] + rest
            cur = [i, j] + rest
            C = np.transpose(C, list(np.argsort(cur)))
            expr = mom.mono_expr(m, k, args[:k - 1])
            for M, coef in mom.canon(expr, k + 1).items():
                q = sp.Rational(coef)
                assert q.q == 1
                out[M] = out.get(M, 0) + sgn * int(q.p) * C
    return {m: T for m, T in out.items() if np.any(T)}


def is_zero(cochain):
    return all(not T.any() for T in cochain.values())


def flatten(cochain, monos, dim, k):
    """vector over (monomial, tensor entries) in the order monos."""
    parts = []
    for m in monos:
        T = cochain.get(m)
        parts.append(np.zeros(dim ** k, dtype=np.int64) if T is None else T.reshape(-1))
    return np.concatenate(parts)


def all_monos(mom, k, degs):
    """canonical monomials over gens(k) with degree degs[c] in component c."""
    H = mom.H
    nvar = k - 1
    per = []
    for c in range(H):
        per.append([e for e in itertools.product(range(degs[c] + 1), repeat=nvar) if sum(e) == degs[c]])
    out = []
    for combo in itertools.product(*per):
        m = [0] * (nvar * H)
        for c in range(H):
            for i in range(nvar):
                m[i * H + c] = combo[c][i]
        out.append(tuple(m))
    return sorted(out)


def exact_rank(A):
    """exact rank of an integer matrix A (numpy int64, rows x cols) via the Gram matrix A^T A."""
    A = A[np.any(A != 0, axis=1)]
    if A.shape[0] == 0:
        return 0
    A = A.copy()
    for j in range(A.shape[1]):
        gcol = int(np.gcd.reduce(np.abs(A[:, j]))) if A[:, j].any() else 0
        if gcol > 1:
            A[:, j] //= gcol
    mx = int(np.abs(A).max())
    if A.shape[0] * mx * mx < 2 ** 62:
        G = (A.T @ A).astype(object)
    else:
        Ao = A.astype(object)
        G = Ao.T.dot(Ao)
    from fractions import Fraction as Fr
    M = [[Fr(int(x)) for x in row] for row in G.tolist()]
    n = len(M)
    r = 0
    for c in range(n):
        p = next((i for i in range(r, n) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    return r
