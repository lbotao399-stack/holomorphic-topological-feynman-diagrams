"""Locality of Ohm's law at the level of the full kernel form (T = 0, 1).

Network Gamma (htk Graph), subnetwork S = set of edges, terminals dS = V(S) cap V(O), O = the
other edges, internal vertices S0 = V(S) - dS.  At a point t (exact Fractions):

  Kron data   L_S Laplacian of S (conductances 1/t_e), vertices ordered (dS, S0);
              Lr = L_dd - L_d0 L_00^{-1} L_0d  (k x k response),  E = -L_00^{-1} L_0d  (harmonic
              extension, m x k),  h = E^T  (source transfer:  lambda'_dS = lambda_dS + h lambda_S0).
  S forms     psi^i_r(S) = sum_{e in S} C^S_re y^i_e dt_e  (cycles of S),  U_S = det C^S t C^S^T.
  outer forms cycles of Gamma not in S, S-part replaced by the electrical flow f^(r) with the same
              terminal divergence d^(r):
              psi^i_r(out) = sum_{e in O} C_re y^i_e dt_e + sum_{e in S} f^(r)_e y^i_e dt_e
              base form:  sum_{e in O} C_re y^i_e dt_e + d^(r).dZ.j^i + d^(r).Z.dh.lambda^i_S0,
              Z = Lr^+ (grounded at the last terminal), j = div of the harmonic part of y on S.
  T = 1       Gaussian free field phi ~ N(0, 1/2 Lambda^{-1}), w_e = (B phi)_e / t_e; Markov
              decomposition phi_S0 = E phi_dS + phi^D, phi^D ~ N(0, 1/2 L_00^{-1}) independent of
              the reduced field on V(O) cup dS, whose law is N(0, 1/2 Lambda'^{-1}).

Differential forms: dicts {bitmask over Gamma edges: Fraction} (htk.linalg conventions).
"""
import os
import random
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '../../../code'))
if CODE not in sys.path:
    sys.path.insert(0, CODE)

from htk.graph import Graph                                   # noqa: E402
from htk import brute                                         # noqa: E402
from htk.linalg import (det, solve, inverse, one_form, wedge, wedge_all, form_add,   # noqa: E402
                        form_scale, basic_coefficient)


# ---------------------------------------------------------------------- small helpers

def zeros(r, c):
    return [[Fr(0)] * c for _ in range(r)]


def mm(A, B):
    if not A or not B:
        return [[Fr(0)] * (len(B[0]) if B else 0) for _ in A]
    Bt = list(zip(*B))
    return [[sum((a * b for a, b in zip(row, col)), Fr(0)) for col in Bt] for row in A]


def tr(A):
    return [list(r) for r in zip(*A)] if A else []


def madd(A, B, s=1):
    return [[a + s * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def form1(coef):
    """{edge: coefficient} -> one-form"""
    return {1 << e: Fr(c) for e, c in coef.items() if c != 0}


def fadd(*fs):
    out = {}
    for f in fs:
        out = form_add(out, f)
    return out


def spanning_tree(nv_list, edges, order):
    """Kruskal over edge indices in `order`; returns set of tree edges"""
    parent = {v: v for v in nv_list}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    tree = []
    for e in order:
        a, b = edges[e]
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
            tree.append(e)
    return tree


def fundamental_cycles(edges, tree, cot, ne):
    """rows over all ne edges; cycle runs along cotree edge c and back through the tree"""
    adj = {}
    for e in tree:
        a, b = edges[e]
        adj.setdefault(a, []).append((b, e, +1))
        adj.setdefault(b, []).append((a, e, -1))
    C = []
    for c in cot:
        a, b = edges[c]
        prev = {b: None}
        stack = [b]
        while stack:
            v = stack.pop()
            if v == a:
                break
            for (u, e, s) in adj.get(v, []):
                if u not in prev:
                    prev[u] = (v, e, s)
                    stack.append(u)
        row = [0] * ne
        row[c] = 1
        v = a
        while prev[v] is not None:
            u, e, s = prev[v]
            row[e] += s
            v = u
        C.append(row)
    return C


# ---------------------------------------------------------------------- the split

class Split:
    def __init__(self, G, S_edges, terminal_order=None):
        self.G = G
        self.S = sorted(S_edges)
        sS = set(self.S)
        self.O = [e for e in range(G.ne) if e not in sS]
        VS = set(v for e in self.S for v in G.edges[e])
        VO = set(v for e in self.O for v in G.edges[e])
        term = sorted(VS & VO)
        if terminal_order is not None:
            assert sorted(terminal_order) == term
            term = list(terminal_order)
        self.dS = term
        self.S0 = sorted(VS - VO)
        self.VO = sorted(VO)
        self.k, self.m = len(self.dS), len(self.S0)
        self.Sv = self.dS + self.S0                      # S vertex order: terminals, internal
        self.sidx = {v: i for i, v in enumerate(self.Sv)}
        if G.root in self.S0:
            raise ValueError('root must not be internal to S')
        # cycle bases: S cycles (Kruskal on S), outer cycles (cotree edges in O)
        allv = list(range(G.nv))
        tree = spanning_tree(allv, G.edges, self.S + self.O)
        ts = set(tree)
        self.treeS = [e for e in self.S if e in ts]
        if len(self.treeS) != len(self.Sv) - 1:
            raise ValueError('S not connected')
        cotS = [e for e in self.S if e not in ts]
        cotO = [e for e in self.O if e not in ts]
        self.CS = fundamental_cycles(G.edges, tree, cotS, G.ne)
        self.Cout_raw = fundamental_cycles(G.edges, tree, cotO, G.ne)
        self.nS, self.nout = len(self.CS), len(self.Cout_raw)
        # terminal divergence of the S-part of each outer cycle
        self.dout = []
        for row in self.Cout_raw:
            d = [Fr(0)] * self.k
            for e in self.S:
                a, b = G.edges[e]
                if row[e]:
                    if a in self.sidx and self.sidx[a] < self.k:
                        d[self.sidx[a]] += row[e]
                    if b in self.sidx and self.sidx[b] < self.k:
                        d[self.sidx[b]] -= row[e]
            self.dout.append(d)

    # ------------------------------------------------------------------ Kron data at t
    def kron(self, t):
        G = self.G
        nS = len(self.Sv)
        L = zeros(nS, nS)
        for e in self.S:
            a, b = G.edges[e]
            i, j = self.sidx[a], self.sidx[b]
            c = 1 / Fr(t[e])
            L[i][i] += c
            L[j][j] += c
            L[i][j] -= c
            L[j][i] -= c
        k, m = self.k, self.m
        Ldd = [r[:k] for r in L[:k]]
        Ld0 = [r[k:] for r in L[:k]]
        L0d = [r[:k] for r in L[k:]]
        L00 = [r[k:] for r in L[k:]]
        if m:
            L00i = inverse(L00)
            X = mm(L00i, L0d)                            # m x k
            Lr = madd(Ldd, mm(Ld0, X), -1)
            E = [[-x for x in r] for r in X]
        else:
            L00i, X, E, Lr = [], [], [], Ldd
        Zg = inverse([r[:k - 1] for r in Lr[:k - 1]]) if k > 1 else []
        return dict(L=L, Ldd=Ldd, Ld0=Ld0, L0d=L0d, L00=L00, L00i=L00i, X=X, E=E, Lr=Lr, Zg=Zg)

    def dkron(self, t, K, e):
        """derivatives d/dt_e of Lr, E, Zg for an S-edge e (exact)"""
        G = self.G
        k, m = self.k, self.m
        a, b = G.edges[e]
        i, j = self.sidx[a], self.sidx[b]
        c = -1 / Fr(t[e]) ** 2
        nS = len(self.Sv)
        dL = zeros(nS, nS)
        dL[i][i] += c
        dL[j][j] += c
        dL[i][j] -= c
        dL[j][i] -= c
        dLdd = [r[:k] for r in dL[:k]]
        dLd0 = [r[k:] for r in dL[:k]]
        dL0d = [r[:k] for r in dL[k:]]
        dL00 = [r[k:] for r in dL[k:]]
        X = K['X']
        if m:
            XT = tr(X)
            dLr = madd(madd(madd(dLdd, mm(dLd0, X), -1), mm(XT, dL0d), -1), mm(mm(XT, dL00), X), 1)
            dX = mm(K['L00i'], madd(dL0d, mm(dL00, X), -1))          # d(L00^-1 L0d)
            dE = [[-x for x in r] for r in dX]
        else:
            dLr, dE = dLdd, []
        Zg = K['Zg']
        dLrg = [r[:k - 1] for r in dLr[:k - 1]]
        dZg = [[-x for x in r] for r in mm(mm(Zg, dLrg), Zg)] if k > 1 else []
        return dLr, dE, dZg

    # ------------------------------------------------------------------ currents
    def flow_from_div(self, t, K, d):
        """electrical flow in S with terminal divergence d (sum d = 0), as {edge: value}"""
        k = self.k
        Zg = K['Zg']
        phig = [sum((Zg[a][b] * d[b] for b in range(k - 1)), Fr(0)) for a in range(k - 1)]
        phid = phig + [Fr(0)]
        phi0 = [sum((K['E'][v][a] * phid[a] for a in range(k)), Fr(0)) for v in range(self.m)]
        phi = phid + phi0
        out = {}
        for e in self.S:
            a, b = self.G.edges[e]
            out[e] = (phi[self.sidx[a]] - phi[self.sidx[b]]) / Fr(t[e])
        return out, phid


def lam_full(G, lam):
    """sources over all vertices (H-vectors of Fractions), root = -sum"""
    ln = brute.sources(G, lam)
    H = len(ln[0]) if ln else 0
    return ln + [[-sum((v[i] for v in ln), Fr(0)) for i in range(H)]], H


# ---------------------------------------------------------------------- T = 0 factorization

def factor_T0(sp, t, lam, mode='R'):
    """returns (D_brute, D_factorized) at the point t.
    mode 'R': outer forms from response differentials dZ, dh (any S, any k).
    mode 'P': outer forms without the dh term (k = 2, tight S: Gamma' kernel at fixed sources)."""
    G = sp.G
    t = [Fr(x) for x in t]
    D_b = brute.density(G, t, lam).coeff
    lf, H = lam_full(G, lam)
    y = brute.currents(G, t, lam)
    K = sp.kron(t, )
    k, m = sp.k, sp.m
    # S forms
    MS = [[sum((Fr(ra[e] * rb[e]) * t[e] for e in sp.S), Fr(0)) for rb in sp.CS] for ra in sp.CS]
    US = det(MS)
    Sforms = [form1({e: sp.CS[r][e] * y[e][i] for e in sp.S}) for i in range(H) for r in range(sp.nS)]
    # harmonic decomposition on S: j = div(y_S) + h lambda_S0  (div = net current from terminal into S)
    divS = [[Fr(0)] * H for _ in range(k)]
    for e in sp.S:
        a, b = G.edges[e]
        for i in range(H):
            if sp.sidx[a] < k:
                divS[sp.sidx[a]][i] += y[e][i]
            if sp.sidx[b] < k:
                divS[sp.sidx[b]][i] -= y[e][i]
    lam0 = [lf[v] for v in sp.S0]
    hl = [[sum((K['E'][v][a] * lam0[v][i] for v in range(m)), Fr(0)) for i in range(H)] for a in range(k)]
    j = [[divS[a][i] + hl[a][i] for i in range(H)] for a in range(k)]
    # derivatives of Zg and h on S edges
    dZ, dh = {}, {}
    for e in sp.S:
        _dLr, dE, dZg = sp.dkron(t, K, e)
        dZ[e] = dZg
        dh[e] = dE
    Zg = K['Zg']
    # outer forms (base) and M_out
    outforms_i = []
    Mout = zeros(sp.nout, sp.nout)
    for r, row in enumerate(sp.Cout_raw):
        d = sp.dout[r]
        for s, row2 in enumerate(sp.Cout_raw):
            d2 = sp.dout[s]
            val = sum((Fr(row[e] * row2[e]) * t[e] for e in sp.O), Fr(0))
            val += sum((d[a] * Zg[a][b] * d2[b] for a in range(k - 1) for b in range(k - 1)), Fr(0))
            Mout[r][s] = val
    for i in range(H):
        for r, row in enumerate(sp.Cout_raw):
            d = sp.dout[r]
            coef = {e: row[e] * y[e][i] for e in sp.O if row[e]}
            for e in sp.S:
                v = sum((d[a] * dZ[e][a][b] * j[b][i] for a in range(k - 1) for b in range(k - 1)), Fr(0))
                if mode == 'R' and m:
                    dhl = [sum((dh[e][v2][a] * lam0[v2][i] for v2 in range(m)), Fr(0)) for a in range(k)]
                    v += sum((d[a] * Zg[a][b] * dhl[b] for a in range(k - 1) for b in range(k - 1)), Fr(0))
                if v:
                    coef[e] = coef.get(e, 0) + v
            outforms_i.append((i, r, form1(coef)))
    Uout = det(Mout)
    # order: for each i: S cycles then outer cycles (block per component), as in the brute ordering
    forms = []
    for i in range(H):
        forms += [Sforms[i * sp.nS + r] for r in range(sp.nS)]
        forms += [f for (ii, r, f) in outforms_i if ii == i]
    W = wedge_all(forms)
    D_f = basic_coefficient(W, t) / (US * Uout) ** H
    return D_b, D_f, dict(US=US, Uout=Uout)


# ---------------------------------------------------------------------- T = 1: Gaussian Wick

class GField:
    """odd field  sum_v F[v] g_v  (+ mean, a 1-form polynomial in a formal scalar w:
    mean = {power: form}); g = Gaussian vector with covariance Cov."""
    __slots__ = ('F', 'mean')

    def __init__(self, F, mean=None):
        self.F = F                    # {var: form}
        self.mean = mean or {}        # {power: form}


def cov_forms(X, Y, Cov):
    out = {}
    for u, fu in X.F.items():
        for v, fv in Y.F.items():
            c = Cov.get((u, v), 0)
            if c:
                out = form_add(out, form_scale(wedge(fu, fv), c))
    return out


def poly_wedge(A, B):
    """A, B: {power: form}"""
    out = {}
    for pa, fa in A.items():
        for pb, fb in B.items():
            w = wedge(fa, fb)
            if w:
                out[pa + pb] = form_add(out.get(pa + pb, {}), w)
                if not out[pa + pb]:
                    del out[pa + pb]
    return out


def wick(fields, Cov):
    """< prod_r X_r >  (ordered product of odd fields) as {power of w: form}"""
    if not fields:
        return {0: {0: Fr(1)}}
    X0, rest = fields[0], fields[1:]
    out = poly_wedge(X0.mean, wick(rest, Cov)) if X0.mean else {}
    for jx in range(len(rest)):
        A = cov_forms(X0, rest[jx], Cov)
        if not A:
            continue
        sub = wick(rest[:jx] + rest[jx + 1:], Cov)
        term = poly_wedge({0: A}, sub)
        if jx % 2:
            term = {p: form_scale(f, -1) for p, f in term.items()}
        for p, f in term.items():
            out[p] = form_add(out.get(p, {}), f)
            if not out[p]:
                del out[p]
    return out


def gff_cov(G, t, verts=None):
    """covariance 1/2 Lambda^{-1} of the GFF on the non-root vertices of G (dict over pairs)"""
    Li = inverse(brute.laplacian(G, t))
    g = G.nv - 1
    return {(u, v): Li[u][v] / 2 for u in range(g) for v in range(g)}


def factor_T1(sp, t, lam, mode='R'):
    """T = 1 (zero topological sources).  Brute density vs the factorized Gaussian:
    reduced GFF on V(O) cup dS (covariance 1/2 Lambda'^{-1}) and Dirichlet GFF on S0.
    mode 'R': outer topological forms include the Dirichlet term d.Z.dh.xi, xi = L_00 phi^D;
    mode 'P': k = 2, tight S, the Dirichlet parts of the outer forms dropped.
    Returns (Density_brute, coeff, radicand) with D = coeff * radicand^(-1/2)."""
    G = sp.G
    t = [Fr(x) for x in t]
    Db = brute.density(G, t, lam, T=1)
    lf, H = lam_full(G, lam)
    y = brute.currents(G, t, lam) if H else [[] for _ in range(G.ne)]
    K = sp.kron(t)
    k, m = sp.k, sp.m
    MS = [[sum((Fr(ra[e] * rb[e]) * t[e] for e in sp.S), Fr(0)) for rb in sp.CS] for ra in sp.CS]
    US = det(MS)
    Zg = K['Zg']
    # reduced network Gamma' (outer vertices + terminals, Laplacian = O part + Lr)
    Vp = sorted(set(sp.VO) | set(sp.dS))
    if G.root not in Vp:
        raise ValueError('root must be in V(O) or dS')
    Vp_nr = [v for v in Vp if v != G.root]
    ip = {v: i for i, v in enumerate(Vp_nr)}
    g = len(Vp_nr)
    Lam = zeros(g, g)
    for e in sp.O:
        a, b = G.edges[e]
        c = 1 / t[e]
        for u, w in ((a, b), (b, a)):
            if u != G.root:
                Lam[ip[u]][ip[u]] += c
                if w != G.root:
                    Lam[ip[u]][ip[w]] -= c
    for a in range(k):
        for b in range(k):
            u, w = sp.dS[a], sp.dS[b]
            if u != G.root and w != G.root:
                Lam[ip[u]][ip[w]] += K['Lr'][a][b]
    Li = inverse(Lam)
    Cov = {}
    for u in Vp_nr:
        for w in Vp_nr:
            Cov[('R', u), ('R', w)] = Li[ip[u]][ip[w]] / 2
    for a in range(m):
        for b in range(m):
            Cov[('D', a), ('D', b)] = K['L00i'][a][b] / 2
    # vertex potential of the topological field as linear form in Gaussian variables:
    # phi_v (v in Vp) = R_v ; phi_{S0[a]} = sum_b E[a][b] phi_dS[b] + D_a
    def pot(v):
        if v == G.root:
            return {}
        if v in ip:
            return {('R', v): Fr(1)}
        a = sp.S0.index(v)
        out = {('D', a): Fr(1)}
        for b in range(k):
            u = sp.dS[b]
            if u != G.root and K['E'][a][b]:
                out[('R', u)] = out.get(('R', u), 0) + K['E'][a][b]
        return out

    def wcur(e):
        a, b = G.edges[e]
        out = {}
        for var, c in pot(a).items():
            out[var] = out.get(var, 0) + c / t[e]
        for var, c in pot(b).items():
            out[var] = out.get(var, 0) - c / t[e]
        return out
    W = {e: wcur(e) for e in range(G.ne)}

    def field_from(coef_edges):
        """sum_e coef_e w_e dt_e  -> GField"""
        F = {}
        for e, c in coef_edges.items():
            for var, x in W[e].items():
                F[var] = form_add(F.get(var, {}), {1 << e: c * x})
        return GField({v: f for v, f in F.items() if f})
    tauS = [field_from({e: sp.CS[r][e] for e in sp.S if sp.CS[r][e]}) for r in range(sp.nS)]
    # outer topological forms: O part + d.dZ.j^tau (j^tau = L phi_dS... harmonic flow) + Dirichlet
    dZ, dh = {}, {}
    for e in sp.S:
        _dLr, dE, dZg = sp.dkron(t, K, e)
        dZ[e], dh[e] = dZg, dE
    # j^tau = terminal injection of the harmonic part = Lr phi_dS (linear in R variables)
    jtau = []
    for a in range(k):
        lin = {}
        for b in range(k):
            u = sp.dS[b]
            if u != G.root and K['Lr'][a][b]:
                lin[('R', u)] = lin.get(('R', u), 0) + K['Lr'][a][b]
        jtau.append(lin)
    # xi = L00 phi^D
    xi = [{('D', b): K['L00'][a][b] for b in range(m) if K['L00'][a][b]} for a in range(m)]
    tauout = []
    for r, row in enumerate(sp.Cout_raw):
        d = sp.dout[r]
        F = {}
        for e in sp.O:
            if row[e]:
                for var, x in W[e].items():
                    F[var] = form_add(F.get(var, {}), {1 << e: row[e] * x})
        for e in sp.S:
            for b in range(k - 1):
                cb = sum((d[a] * dZ[e][a][b] for a in range(k - 1)), Fr(0))
                if cb:
                    for var, x in jtau[b].items():
                        F[var] = form_add(F.get(var, {}), {1 << e: cb * x})
            if mode == 'R' and m:
                # d.Z.dh.xi
                for v2 in range(m):
                    cv = sum((d[a] * Zg[a][b] * dh[e][v2][b] for a in range(k - 1) for b in range(k - 1)), Fr(0))
                    if cv:
                        for var, x in xi[v2].items():
                            F[var] = form_add(F.get(var, {}), {1 << e: cv * x})
        tauout.append(GField({v: f for v, f in F.items() if f}))
    # holomorphic forms (as in factor_T0, mode R or P)
    Mout = zeros(sp.nout, sp.nout)
    for r, row in enumerate(sp.Cout_raw):
        for s, row2 in enumerate(sp.Cout_raw):
            val = sum((Fr(row[e] * row2[e]) * t[e] for e in sp.O), Fr(0))
            val += sum((sp.dout[r][a] * Zg[a][b] * sp.dout[s][b] for a in range(k - 1) for b in range(k - 1)), Fr(0))
            Mout[r][s] = val
    Uout = det(Mout)
    hol = []
    if H:
        divS = [[Fr(0)] * H for _ in range(k)]
        for e in sp.S:
            a, b = G.edges[e]
            for i in range(H):
                if sp.sidx[a] < k:
                    divS[sp.sidx[a]][i] += y[e][i]
                if sp.sidx[b] < k:
                    divS[sp.sidx[b]][i] -= y[e][i]
        lam0 = [lf[v] for v in sp.S0]
        hl = [[sum((K['E'][v][a] * lam0[v][i] for v in range(m)), Fr(0)) for i in range(H)] for a in range(k)]
        j = [[divS[a][i] + hl[a][i] for i in range(H)] for a in range(k)]
        for i in range(H):
            hol += [form1({e: sp.CS[r][e] * y[e][i] for e in sp.S}) for r in range(sp.nS)]
            for r, row in enumerate(sp.Cout_raw):
                d = sp.dout[r]
                coef = {e: row[e] * y[e][i] for e in sp.O if row[e]}
                for e in sp.S:
                    v = sum((d[a] * dZ[e][a][b] * j[b][i] for a in range(k - 1) for b in range(k - 1)), Fr(0))
                    if mode == 'R' and m:
                        dhl = [sum((dh[e][v2][a] * lam0[v2][i] for v2 in range(m)), Fr(0)) for a in range(k)]
                        v += sum((d[a] * Zg[a][b] * dhl[b] for a in range(k - 1) for b in range(k - 1)), Fr(0))
                    if v:
                        coef[e] = coef.get(e, 0) + v
                hol.append(form1(coef))
    top = wick(tauS + tauout, Cov).get(0, {})
    Wf = wedge(top, wedge_all(hol)) if hol else top
    U = US * Uout
    coeff = basic_coefficient(Wf, t) / U ** H
    return Db, coeff, U


def same_density(Db, coeff, rad):
    """Db = c U^-1/2 with U = Db.radicand; compare sign and square, return ratio sign or None"""
    if Db.radicand != rad:
        return None
    if Db.coeff == coeff:
        return 1
    if Db.coeff == -coeff:
        return -1
    return None


# ---------------------------------------------------------------------- k = 2: chart-level form

class TwoTerminal:
    """k = 2 split with S standalone (root q = second terminal) and Gamma' = O + e* (p -> q, last).
    K_Gamma = kappa * int_{t_S} int_{t'} D_S(t_S; sigma, lambda_S0) D_Gamma'(t'; lambda'(t_S)),
    sigma = y'_{e*} - h_p lambda_S0 = current from p into S,
    kappa = D_Gamma |J_Phi| / (D_S D')  (pointwise identity D_Gamma J_Phi = eps D_S D'; kappa = eps sgn J_Phi)."""

    def __init__(self, G, S_edges, p=None, q=None):
        self.sp = sp = Split(G, S_edges)
        if sp.k != 2:
            raise ValueError('not two-terminal')
        if p is None:
            p, q = sp.dS
        self.p, self.q = p, q
        self.G = G
        # S standalone: internal vertices, p, q(root)
        sv = sp.S0 + [p, q]
        self.Sidx = {v: i for i, v in enumerate(sv)}
        self.GS = Graph(len(sv), [(self.Sidx[G.edges[e][0]], self.Sidx[G.edges[e][1]]) for e in sp.S])
        # Gamma': vertices of O plus terminals in original order (root kept last)
        vp = sorted(set(sp.VO) | {p, q})
        if vp[-1] != G.root:
            raise ValueError('root of Gamma must survive in Gamma\'')
        self.Pidx = {v: i for i, v in enumerate(vp)}
        self.vp = vp
        self.Gp = Graph(len(vp), [(self.Pidx[G.edges[e][0]], self.Pidx[G.edges[e][1]]) for e in sp.O]
                        + [(self.Pidx[p], self.Pidx[q])])

    def reff(self, tS):
        """R_eff of S (tS over S edges in order) and harmonic measure h_p(v), v in S0"""
        GS = self.GS
        L = brute.laplacian(GS, tS)                         # q grounded
        Li = inverse(L)
        ip = self.Sidx[self.p]
        R = Li[ip][ip]
        # harmonic measure of p seen from internal v: potential at v for unit potential at p,
        # q grounded, p held: h_v = Li[v][p] / Li[p][p]
        h = [Li[self.Sidx[v]][ip] / R for v in self.sp.S0]
        return R, h, Li

    def point(self, tS, tp, lam):
        """Gamma times from (t_S, t'), with R_eff(t_S) rescaled to t'_{e*}"""
        R, h, _ = self.reff(tS)
        c = Fr(tp[-1]) / R
        t = [None] * self.G.ne
        for j, e in enumerate(self.sp.S):
            t[e] = Fr(tS[j]) * c
        for j, e in enumerate(self.sp.O):
            t[e] = Fr(tp[j])
        return t

    def sources_prime(self, tS, lam):
        G, sp = self.G, self.sp
        lf, H = lam_full(G, lam)
        R, h, _ = self.reff(tS)
        out = {}
        for v in self.vp:
            out[v] = list(lf[v])
        s = [sum((h[a] * lf[v][i] for a, v in enumerate(sp.S0)), Fr(0)) for i in range(H)]
        tot0 = [sum((lf[v][i] for v in sp.S0), Fr(0)) for i in range(H)]
        out[self.p] = [out[self.p][i] + s[i] for i in range(H)]
        out[self.q] = [out[self.q][i] + tot0[i] - s[i] for i in range(H)]
        lamp = [out[v] for v in self.vp]
        return lamp, s, H

    def check_point(self, tS, tp, lam, T=0):
        """returns (D_Gamma * J_Phi, D_S * D_Gamma') at the point (exact)"""
        G, sp = self.G, self.sp
        tS = [Fr(x) for x in tS]
        tp = [Fr(x) for x in tp]
        t = self.point(tS, tp, lam)
        lamp, s, H = self.sources_prime(tS, lam)
        Dp = brute.density(self.Gp, tp, lamp, T=T)
        yp = brute.currents(self.Gp, tp, lamp)
        J = yp[-1]
        sigma = [J[i] - s[i] for i in range(H)]
        lf, _ = lam_full(G, lam)
        lamS = [lf[v] for v in sp.S0] + [sigma]
        DS = brute.density(self.GS, tS, lamS, T=T) if T == 0 else None
        DG = brute.density(G, t, lam, T=T)
        # Jacobian of Phi on the charts t_S[-1] fixed, t'[-1] fixed:
        # Phi^* iota_E vol_Gamma (u, v) = det[t; dPhi u; dPhi v],  (iota vol_S ^ iota vol')(u, v) = det[tS; u] det[t'; v]
        R, h, Li = self.reff(tS)
        ip = self.Sidx[self.p]
        # dR/dt_j = a_j^2 (a = unit flow p -> q)
        phi = [Li[x][ip] for x in range(len(Li))] + [Fr(0)]
        a = [(phi[self.GS.edges[j][0]] - phi[self.GS.edges[j][1]]) / tS[j] for j in range(len(tS))]
        NS, Np = len(tS), len(tp)
        rows = [t]
        c = tp[-1] / R
        for jj in range(NS - 1):
            v = [Fr(0)] * G.ne
            for j2, e in enumerate(sp.S):
                v[e] = c * ((1 if j2 == jj else 0) - tS[j2] * a[jj] ** 2 / R)
            rows.append(v)
        for jj in range(Np - 1):
            v = [Fr(0)] * G.ne
            v[sp.O[jj]] = Fr(1)
            rows.append(v)
        num = det(rows)
        den = det([tS] + [[Fr(int(i == jj)) for i in range(NS)] for jj in range(NS - 1)]) * \
            det([tp] + [[Fr(int(i == jj)) for i in range(Np)] for jj in range(Np - 1)])
        JPhi = num / den
        return DG, JPhi, DS, Dp, sigma


def rand_times(n, rng, lo=1, hi=9):
    return [Fr(rng.randint(lo, hi), rng.randint(1, 5)) for _ in range(n)]


# ---------------------------------------------------------------------- T = 1, k = 2: S alone x Gamma' alone

def cycle_forms_tau(G, t, C, Wfun):
    """psi^tau_r = sum_e C_re w_e dt_e as GFields; Wfun(e) = {var: coeff} (w_e linear in Gaussians)"""
    out = []
    for row in C:
        F = {}
        for e, c in enumerate(row):
            if c:
                for var, x in Wfun(e).items():
                    F[var] = form_add(F.get(var, {}), {1 << e: c * x})
        out.append(GField({v: f for v, f in F.items() if f}))
    return out


def density_T1_insert(G, t, lam, f, jmax):
    """[coefficients c_j], D^{[w_f^j]} = c_j U^{-1/2}, D^{[F]} = U^{-H-1/2} <F prod psi^tau> ^ prod psi^i / iota_E vol,
    GFF phi ~ N(0, 1/2 Lambda^{-1}), w_e = (B phi)_e / t_e, insertion w_f^j (j = 0..jmax)."""
    from math import factorial
    t = [Fr(x) for x in t]
    C = G.cycle_matrix()
    g = G.nv - 1
    Li = inverse(brute.laplacian(G, t))
    Cov = {(u, v): Li[u][v] / 2 for u in range(g) for v in range(g)}

    def Wfun(e):
        a, b = G.edges[e]
        out = {}
        if a < g:
            out[a] = out.get(a, 0) + 1 / t[e]
        if b < g:
            out[b] = out.get(b, 0) - 1 / t[e]
        return out
    taus = cycle_forms_tau(G, t, C, Wfun)
    wf = Wfun(f)
    # <exp(mu w) prod X> = exp(mu^2 V/2) < prod (X + mu Cov(w, X)) >
    V = sum((wf[u] * wf[v] * Cov[(u, v)] for u in wf for v in wf), Fr(0))
    fields = []
    for X in taus:
        m1 = {}
        for u, cu in wf.items():
            for v, fv in X.F.items():
                c = cu * Cov[(u, v)]
                if c:
                    m1 = form_add(m1, form_scale(fv, c))
        fields.append(GField(X.F, {1: m1} if m1 else {}))
    poly = wick(fields, Cov)                               # {power of mu: form}
    lamn = brute.sources(G, lam)
    H = len(lamn[0]) if lamn else 0
    y = brute.currents(G, t, lam) if H else None
    hol = [form1({e: C[r][e] * y[e][i] for e in range(G.ne)}) for i in range(H) for r in range(len(C))]
    holw = wedge_all(hol) if hol else {0: Fr(1)}
    U = det(brute.period_matrix(G, t, C))
    coeffs = []
    for j in range(jmax + 1):
        # j! [mu^j] exp(mu^2 V/2) poly
        tot = {}
        for l in range(j // 2 + 1):
            p = j - 2 * l
            if p in poly:
                tot = form_add(tot, form_scale(poly[p], (V / 2) ** l / factorial(l)))
        tot = form_scale(tot, factorial(j))
        w = wedge(tot, holw)
        coeffs.append(basic_coefficient(w, t) / U ** H if w else Fr(0))
    return coeffs, U


def rho_S_T1(tt, tS, lamS):
    """rho_S^{(j)}: U_S^{-H-1/2} < prod_r sum_e C_re (a_e omega + w^D_e) dt_e >_D ^ prod psi^i(S),
    coefficient of omega^j, on the S standalone graph (q root), at the times tS.
    Dirichlet GFF on S0 with covariance 1/2 L_00^{-1} (p, q grounded).  Returns ([c_j], U_S)."""
    GS = tt.GS
    t = [Fr(x) for x in tS]
    C = GS.cycle_matrix()
    m = len(tt.sp.S0)
    g = GS.nv - 1                     # non-root: S0 (0..m-1), p (m)
    Lfull = brute.laplacian(GS, t)
    L00 = [r[:m] for r in Lfull[:m]]
    L00i = inverse(L00) if m else []
    Cov = {(u, v): L00i[u][v] / 2 for u in range(m) for v in range(m)}
    Li = inverse(Lfull)
    phi = [Li[x][m] for x in range(g)] + [Fr(0)]             # unit injection at p
    R = Li[m][m]
    a = [(phi[GS.edges[e][0]] - phi[GS.edges[e][1]]) / t[e] / R * R for e in range(GS.ne)]
    # a: currents for unit current p -> q
    a = [(phi[GS.edges[e][0]] - phi[GS.edges[e][1]]) / t[e] for e in range(GS.ne)]
    fields = []
    for row in C:
        F, mean = {}, {}
        for e, c in enumerate(row):
            if not c:
                continue
            mean = form_add(mean, {1 << e: c * a[e]})
            u, w = GS.edges[e]
            for vv, s in ((u, 1), (w, -1)):
                if vv < m:
                    F[vv] = form_add(F.get(vv, {}), {1 << e: c * s / t[e]})
        fields.append(GField({v: f for v, f in F.items() if f}, {1: mean} if mean else {}))
    poly = wick(fields, Cov)
    lamn = brute.sources(GS, lamS)
    H = len(lamn[0]) if lamn else 0
    y = brute.currents(GS, t, lamS) if H else None
    hol = [form1({e: C[r][e] * y[e][i] for e in range(GS.ne)}) for i in range(H) for r in range(len(C))]
    holw = wedge_all(hol) if hol else {0: Fr(1)}
    U = det(brute.period_matrix(GS, t, C))
    out = []
    for j in range(len(C) + 1):
        w = wedge(poly.get(j, {}), holw)
        out.append(basic_coefficient(w, t) / U ** H if w else Fr(0))
    return out, U


def check_point_T1(tt, tS, tp_O, lam):
    """chart-level T = 1 identity at a point with t'_{e*} = R_eff(t_S) (scale factor 1):
    D_Gamma J_Phi = kappa sum_j rho_S^{(j)}(t_S) D'^{[w_{e*}^j]}(t').  Returns (lhs, [rho_j], [D'_j], U's)."""
    tS = [Fr(x) for x in tS]
    R, h, Li = tt.reff(tS)
    tp = [Fr(x) for x in tp_O] + [R]
    G = tt.G
    t = tt.point(tS, tp, lam)
    lamp, s, H = tt.sources_prime(tS, lam)
    yp = brute.currents(tt.Gp, tp, lamp) if H else None
    sigma = [yp[-1][i] - s[i] for i in range(H)]
    lf, _ = lam_full(G, lam)
    lamS = [lf[v] for v in tt.sp.S0] + [sigma] if H else None
    rho, US = rho_S_T1(tt, tS, lamS)
    Dps, Up = density_T1_insert(tt.Gp, tp, lamp, tt.Gp.ne - 1, len(rho) - 1)
    DG = brute.density(G, t, lam, T=1)
    # Jacobian (same as T = 0 chart)
    ip = tt.Sidx[tt.p]
    phi = [Li[x][ip] for x in range(len(Li))] + [Fr(0)]
    a = [(phi[tt.GS.edges[j][0]] - phi[tt.GS.edges[j][1]]) / tS[j] for j in range(len(tS))]
    NS, Np = len(tS), len(tp)
    rows = [t]
    for jj in range(NS - 1):
        v = [Fr(0)] * G.ne
        for j2, e in enumerate(tt.sp.S):
            v[e] = (1 if j2 == jj else 0) - tS[j2] * a[jj] ** 2 / R
        rows.append(v)
    for jj in range(Np - 1):
        v = [Fr(0)] * G.ne
        v[tt.sp.O[jj]] = Fr(1)
        rows.append(v)
    num = det(rows)
    den = det([tS] + [[Fr(int(i == jj)) for i in range(NS)] for jj in range(NS - 1)]) * \
        det([tp] + [[Fr(int(i == jj)) for i in range(Np)] for jj in range(Np - 1)])
    JPhi = num / den
    return DG, JPhi, rho, US, Dps, Up
