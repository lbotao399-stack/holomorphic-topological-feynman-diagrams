"""Exact weights of the theta kernels on C^4 by Cauchy-Binet (independent of the engine's density layer).

In ratio coordinates xi (sections of kernels_c4.py, currents y_e = sum_v resp[e][v](xi) lam_v):
    g(xi) = orientation * det[ d y_{c,a} / d xi_j ]_{(c,a); j},   rows c in the cotree (outer), a = 1..4 (inner),
    det = sum_{S1, S2} <S1><S2> det[ beta_{c1, S1} ; beta_{c2, S2} ],   beta_{c,v,j} = d resp[c][v] / d xi_j,
    <S> = det(lam_{s_1}, .., lam_{s_4}) (S sorted).
Weights: W(ks)(lam) = int g prod_{j,a} (-s_j y_{e_j, a})^{k_{j,a}}  (e_j spoke of input j, z_{e_j} = s_j z_j),
which equals d_z^{ks} K(lam; z) at z = 0 (S2 eq-moment).  Supported: total order of ks <= 1.
"""
import sys, itertools, time, pickle
from fractions import Fraction as Fr
from math import factorial
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from htk import sp_reduce
from graphs_c4 import theta_O, THETAS
import brute_c4 as bf

H = 4


def bracket_monomials(S):
    """<S> as {(v, a)-exponent dict as sorted tuple of ((v,a),1)): sign}."""
    out = {}
    for perm in itertools.permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        key = tuple(sorted((S[b], perm[b]) for b in range(4)))
        out[key] = out.get(key, 0) + (-1 if inv % 2 else 1)
    return out


class ThetaWeights:
    def __init__(self, pqr, log=print):
        self.pqr = pqr
        G, spokes, outputs, paths = theta_O(*pqr)
        self.G, self.spokes, self.outputs, self.paths = G, spokes, outputs, paths
        O, N = G.root, 0
        P = [p_[0] for p_ in paths]
        plan = [(O, N, [P[0], P[1], P[2]])] if pqr == (3, 3, 3) else [(O, N, [P[0], P[1]]), (O, N, [[G.ne], P[2]])]
        self.red = sp_reduce(G, plan=plan)
        self.coords = list(self.red.coords)
        self.cs = [sp.Symbol(c) for c in self.coords]
        self.blocks = [list(b.indep) for b in self.red.blocks]
        self.cot = bf.cotree_of(G.nv, G.edges)
        R = sp.QQ.old_poly_ring(*self.cs) if False else None
        self.dom = sp.QQ[tuple(self.cs)]
        self.resp = {e: {v: sp.Poly(self.red.resp[e][v].as_expr(), *self.cs) for v in outputs} for e in range(G.ne)}
        t0 = time.time()
        self.dets = {}
        subsets = list(itertools.combinations(outputs, 4))
        c1, c2 = self.cot
        beta = {(c, v): [self.resp[c][v].diff(x) for x in self.cs] for c in self.cot for v in outputs}
        for S1 in subsets:
            for S2 in subsets:
                rows = [beta[(c1, v)] for v in S1] + [beta[(c2, v)] for v in S2]
                M = DomainMatrix([[self.dom.from_sympy(p.as_expr()) for p in r] for r in rows], (8, 8), self.dom)
                d = M.det()
                if d != self.dom.zero:
                    self.dets[(S1, S2)] = sp.Poly(self.dom.to_sympy(d), *self.cs)
        self.orientation = self.red.orientation
        log('theta%s: cotree %s, %d nonzero Cauchy-Binet minors of %d, orientation %+d (%.1fs)'
            % (pqr, self.cot, len(self.dets), len(subsets) ** 2, self.orientation, time.time() - t0))
        self._mono = {}

    def integrate(self, P):
        pos = {c: i for i, c in enumerate(self.coords)}
        tot = Fr(0)
        for mon, c in P.terms():
            w = Fr(1)
            for blk in self.blocks:
                a = [mon[pos[x]] for x in blk]
                num = 1
                for x in a:
                    num *= factorial(x)
                w *= Fr(num, factorial(len(blk) + sum(a)))
            tot += Fr(int(c.p), int(c.q)) * w
        return tot

    def density_value(self, xi, lam):
        """g at numeric xi (dict coord -> Fraction) and lam (dict v -> 4-vector), from the Cauchy-Binet form."""
        sub = tuple(sp.Rational(xi[c].numerator, xi[c].denominator) for c in self.coords)
        val = Fr(0)
        for (S1, S2), P in self.dets.items():
            b1 = sp.Matrix([lam[v] for v in S1]).det()
            b2 = sp.Matrix([lam[v] for v in S2]).det()
            pv = sp.Rational(P.eval(sub)) if len(self.cs) else P
            val += Fr(int(pv.p), int(pv.q)) * int(b1) * int(b2)
        return self.orientation * val

    def weights(self, ks):
        """{m: Fraction}, m = tuple over outputs of 4-tuples of exponents of lam_{v,a}."""
        order = sum(sum(k) for k in ks)
        coef = {}
        if order == 0:
            for key, P in self.dets.items():
                coef[key + (None,)] = self.orientation * self.integrate(P)
        elif order == 1:
            j = next(i for i, k in enumerate(ks) if any(k)) + 1
            a = next(i for i in range(H) if ks[j - 1][i])
            e, s = next((e, s) for e, (jj, s) in self.spokes.items() if jj == j)
            for key, P in self.dets.items():
                for v in self.outputs:
                    Q = P * self.resp[e][v]
                    val = self.orientation * (-s) * self.integrate(Q)
                    if val:
                        coef[key + ((v, a),)] = val
        else:
            raise NotImplementedError
        idx = {v: i for i, v in enumerate(self.outputs)}
        out = {}
        for (S1, S2, extra), c in coef.items():
            if c == 0:
                continue
            for k1, s1 in self._bm(S1).items():
                for k2, s2 in self._bm(S2).items():
                    m = [[0] * H for _ in self.outputs]
                    for (v, aa) in k1 + k2 + ((extra,) if extra else ()):
                        m[idx[v]][aa] += 1
                    key = tuple(tuple(r) for r in m)
                    out[key] = out.get(key, 0) + c * s1 * s2
        return {m: w for m, w in out.items() if w != 0}

    def _bm(self, S):
        if S not in self._mono:
            self._mono[S] = bracket_monomials(S)
        return self._mono[S]
