"""Graphical functions in d = 2 lam + 2 (lam = q - 1) for sequential graphs.

Normalization: A_G(0,1,z) = int prod_int d^d x / pi^{d/2} prod_e |x_e|^{-2 beta_e};  z = complex coordinate,
x = z, y = zbar, w = x - y.  M_lam = Delta_d / 4 = dd - lam T,  dd = d_z d_zbar,  T = w^{-1}(d_z - d_zbar).
Identities (checked, opcheck2.py): [dd, T] = T^2,  M_{mu+1} T = T M_mu,  M_1 = w^{-1} dd w.
Representation V_d:  f = sum_{k=0}^{lam-1} T^k (w^{-1} Phi_k),  Phi_k pole-free SVF.
Appending a weight-lam edge (beta = lam):  M_lam f_new = -f / Gamma(lam)
   X_k = M_1^{-1}[(lam-k) X_{k-1} - A_k / Gamma(lam)],  M_1^{-1}(w^{-1} Psi) = w^{-1} dd^{-1} Psi.
Edge factor r (rational, symmetric):  r T^k A = sum_j binom(k,j) (-1)^{k-j} T^j((T^{k-j} r) A).
Kernel of M_lam (symmetric, rational):  T^{lam-1}(w^{-1}(phi(x) - phi(y))).
Period:  int d^d z / pi^{d/2} f = (-1)^lam lam!/(2 lam)! int_C w^{2 lam} f d^2z / pi.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial

import mzv
import sv
import expand
from rat import R, ONE, factor_pow, pmul
from sv import SVF

W1 = R({(1, 0): 1, (0, 1): -1})          # w = x - y
WINV = R({(0, 0): 1}, (0, 0, 0, 0, 1))  # 1/w


def T_svf(f):
    return (f.dz() - f.dzb()).mulrat(WINV)


def T_rat(r):
    return (r.dx() - r.dy()) * WINV


def edge_factor(b0, b1):
    """|z|^{-2 b0} |1-z|^{-2 b1}  (negative weights: numerators)"""
    num = {(0, 0): 1}
    if b0 < 0:
        num = pmul(num, {(-b0, -b0): 1})
    if b1 < 0:
        num = pmul(num, pmul(factor_pow(2, -b1), factor_pow(3, -b1)))
    return R(num, (max(b0, 0), max(b0, 0), max(b1, 0), max(b1, 0), 0))


class VD:
    def __init__(self, lam, Phi):
        self.lam = lam
        self.Phi = Phi

    @staticmethod
    def rational(lam, r):
        Phi = [SVF() for _ in range(lam)]
        Phi[0] = SVF.rat(r * W1)
        return VD(lam, Phi)

    def explicit(self):
        f = SVF()
        for k in range(self.lam - 1, -1, -1):
            f = f + self.Phi[k].mulrat(WINV)
            if k:
                f = T_svf(f)
        return f

    def mul_edge(self, r):
        lam = self.lam
        Tr = [r]
        for m in range(1, lam):
            t = T_rat(Tr[-1])
            assert t.den[4] == 0, 'T^m r has a (z - zbar) pole'
            Tr.append(t)
        new = []
        for j in range(lam):
            acc = SVF()
            for k in range(j, lam):
                if self.Phi[k].iszero():
                    continue
                acc = acc + self.Phi[k].mulrat(Tr[k - j]).scale(comb(k, j) * (-1) ** (k - j))
            new.append(acc)
        return VD(lam, new)

    def symmetrize(self):
        return VD(self.lam, [(p - p.swap()).scale(F(1, 2)) for p in self.Phi])

    def swap01(self):
        return VD(self.lam, [expand.svf_1mz(p).scale(-1) for p in self.Phi])

    def maxweight(self):
        return max(p.maxweight() for p in self.Phi)

    def size(self):
        return sum(sum(len(d) for d in p.t.values()) for p in self.Phi)


def ddinv(psi):
    return psi.int_z().int_zb()


def append_particular(f):
    lam = f.lam
    g = F(1, factorial(lam - 1))
    X = []
    prev = SVF()
    for k in range(lam):
        psi = prev.scale(lam - k) - f.Phi[k].scale(g)
        Xk = ddinv(psi)
        X.append(Xk)
        prev = Xk
    return VD(lam, X)


# ---------------- kernel fixing ----------------
def _phi_basis(J0, J1, Jinf):
    out = []
    for j in range(1, J0 + 1):
        out.append(('x^-%d' % j, R({(0, 0): 1}, (j, 0, 0, 0, 0)) - R({(0, 0): 1}, (0, j, 0, 0, 0))))
    for j in range(1, J1 + 1):
        out.append(('(x-1)^-%d' % j, R({(0, 0): 1}, (0, 0, j, 0, 0)) - R({(0, 0): 1}, (0, 0, 0, j, 0))))
    for j in range(1, Jinf + 1):
        out.append(('x^%d' % j, R({(j, 0): 1}) - R({(0, j): 1})))
    return out


def _local_conditions(fexp_pf, E, point, M, D):
    """fexp_pf: expansion (dict) of the pole-free G = w^E f in the local variable at `point`.
    Returns (dict of forbidden coefficients {key: mzv}, list of divisibility remainders)."""
    if point == 'inf':
        shifted = {(l, m + E, mb + E): c for (l, m, mb), c in fexp_pf.items()}
        q, rems = expand.homog_divide(shifted, -1, E)
        out = {}
        for (l, m, mb), c in q.items():
            if m + mb > D:
                continue
            if m < -M or mb < -M:
                out[('f', l, m, mb)] = c
        return out, rems
    lin = 1 if point == '0' else -1
    q, rems = expand.homog_divide(fexp_pf, lin, E)
    out = {}
    for (l, m, mb), c in q.items():
        if m + mb > D:
            continue
        if m < M or mb < M:
            out[('f', l, m, mb)] = c
    return out, rems


def _expansions(G, D, E):
    """expansions of pole-free SVF G at 0, 1, inf, complete to the degrees needed"""
    return {'0': expand.expand0(G, D + E),
            '1': expand.expand0(expand.svf_1mz(G), D + E),
            'inf': expand.expand0(expand.svf_inv(G), D + 3 * E)}


def wpole_order(f):
    return max((r.den[4] for w, m, r in f.terms()), default=0)


def fix_kernel(vd, bounds, D=2, J=None, verbose=False):
    """add symmetric rational kernel elements T^{lam-1}(w^{-1}(phi(x)-phi(y))) so that f satisfies
    the local bounds  bounds = {'0': M0, '1': M1, 'inf': Minf}  and is regular on the diagonal."""
    lam = vd.lam
    f = vd.explicit()
    E = max(wpole_order(f), 2 * lam - 1)
    if J is None:
        J0 = max((r.den[0] for w, m, r in f.terms()), default=0) + 2
        J1 = max((r.den[2] for w, m, r in f.terms()), default=0) + 2
        Ji = max((max((i + j for (i, j) in r.num), default=0) - sum(r.den) for w, m, r in f.terms()), default=0) + 2 * lam + 2
        Ji = max(Ji, 2)
    else:
        J0 = J1 = Ji = J
    basis = _phi_basis(J0, J1, Ji)
    wE = R(factor_pow(4, E))
    kers = []
    for name, eta in basis:
        k = SVF.rat(eta * WINV)
        for _ in range(lam - 1):
            k = T_svf(k)
        kers.append(k)
    rows = {}
    rhs = {}
    pts = ('0', '1', 'inf')

    def conds(fun):
        G = fun.mulrat(wE)
        assert not G.has_w_poles()
        ex = _expansions(G, D, E)
        out = {}
        for p in pts:
            fb, rems = _local_conditions(ex[p], E, p, bounds[p], D)
            for key, c in fb.items():
                out[(p,) + key] = c
            for key, c in rems:
                out[(p, 'rem') + key] = c
        return out
    base = conds(f)
    for key, c in base.items():
        rhs[key] = mzv.scal(c, -1)
    for j, k in enumerate(kers):
        for key, c in conds(k).items():
            assert set(c) <= {()}, 'kernel element with MZV coefficient'
            rows.setdefault(key, {})[j] = c.get((), F(0))
    keys = set(rows) | set(rhs)
    # Gaussian elimination: unknowns c_j (MZV-valued)
    eqs = []
    for key in keys:
        eqs.append((dict(rows.get(key, {})), dict(rhs.get(key, {}))))
    piv = {}
    for row, b in eqs:
        row = {j: v for j, v in row.items() if v}
        for j in list(row):
            if j in piv and j in row:
                c = row.pop(j)
                prow, pb = piv[j]
                for jj, vv in prow.items():
                    if jj != j:
                        row[jj] = row.get(jj, 0) - c * vv
                        if not row[jj]:
                            del row[jj]
                b = mzv.add(b, pb, -c)
        if not row:
            if b:
                raise ArithmeticError('kernel fixing inconsistent (%d basis elements); residual %s' % (
                    len(basis), mzv.tostr(b)))
            continue
        p = min(row)
        cp = row[p]
        row = {j: v / cp for j, v in row.items()}
        b = mzv.scal(b, 1 / cp)
        for j2 in list(piv):
            r2, b2 = piv[j2]
            if p in r2:
                c = r2.pop(p)
                for jj, vv in row.items():
                    if jj != p:
                        r2[jj] = r2.get(jj, 0) - c * vv
                        if not r2[jj]:
                            del r2[jj]
                piv[j2] = (r2, mzv.add(b2, b, -c))
        piv[p] = (row, b)
    free = [j for j in range(len(kers)) if j not in piv]
    for j in free:
        if not kers[j].iszero():
            raise ArithmeticError('kernel direction %s undetermined by the conditions' % basis[j][0])
    sol = {}
    for j, (row, b) in piv.items():
        if any(jj != j for jj in row):
            # set free variables to zero
            pass
        sol[j] = b
    if verbose:
        print('    kernel: %d basis, %d conditions, %d fixed, free %s' % (
            len(kers), len(keys), len(piv), [basis[j][0] for j in free]))
    newPhi = list(vd.Phi)
    add = SVF()
    for j, c in sol.items():
        if c:
            add = add + SVF.rat(basis[j][1]).mulmzv(c)
    newPhi[lam - 1] = newPhi[lam - 1] + add
    return VD(lam, newPhi), [basis[j][0] for j in free]


def period(vd):
    """int d^d z / pi^{d/2} f(z)"""
    lam = vd.lam
    f = vd.explicit()
    phi = f.mulrat(R(factor_pow(4, 2 * lam)))
    if phi.has_w_poles():
        raise sv.Obstruction('period integrand has (z-zbar) poles')
    val = expand.integrate_C(phi)
    return mzv.scal(val, F((-1) ** lam * factorial(lam), factorial(2 * lam)))


# ---------------- power counting (beta units) ----------------
def bounds_for(edges, lam):
    """edges: list of (u, v, beta) with external vertices '0', '1', 'z'; others internal."""
    verts = set()
    for u, v, b in edges:
        verts |= {u, v}
    internal = sorted(verts - {'0', '1', 'z'}, key=str)
    out = {}
    for s in ('0', '1'):
        best = None
        for k in range(len(internal) + 1):
            for V in combinations(internal, k):
                S = set(V) | {s, 'z'}
                val = sum(b for u, v, b in edges if u in S and v in S) - (lam + 1) * len(V)
                best = val if best is None else max(best, val)
        out[s] = -best
    best = None
    for k in range(len(internal) + 1):
        for V in combinations(internal, k):
            S = set(V) | {'z'}
            val = sum(b for u, v, b in edges if u in S or v in S) - (lam + 1) * len(V)
            best = val if best is None else min(best, val)
    out['inf'] = -best
    return out


def bounds_path(edges, lam):
    """same as bounds_for, by dynamic programming when the internal vertices form a path"""
    verts = set()
    for u, v, b in edges:
        verts |= {u, v}
    internal = [v for v in verts if v not in ('0', '1', 'z')]
    if len(internal) <= 10:
        return bounds_for(edges, lam)
    adj = {v: {} for v in verts}
    for u, v, b in edges:
        adj[u][v] = adj[u].get(v, 0) + b
        adj[v][u] = adj[v].get(u, 0) + b
    inn = {v: [u for u in adj[v] if u in adj and u not in ('0', '1', 'z')] for v in internal}
    if any(len(x) > 2 for x in inn.values()):
        return bounds_for(edges, lam)
    ends = [v for v in internal if len(inn[v]) <= 1]
    if len(ends) != 2 and len(internal) > 1:
        return bounds_for(edges, lam)
    order = [ends[0]]
    while len(order) < len(internal):
        order.append(next(u for u in inn[order[-1]] if u not in order))
    out = {}
    for s in ('0', '1'):
        const = adj['z'].get(s, 0)
        # dp over order: state = previous vertex in V?
        best = {False: 0, True: None}
        prev = None
        for v in order:
            gain = adj[v].get(s, 0) + adj[v].get('z', 0) - (lam + 1)
            link = adj[v].get(prev, 0) if prev is not None else 0
            nb = {False: max(x for x in best.values() if x is not None),
                  True: max(x for x in (best[False] + gain if best[False] is not None else None,
                                        best[True] + gain + link if best[True] is not None else None) if x is not None)}
            best, prev = nb, v
        out[s] = -(const + max(x for x in best.values() if x is not None))
    const = sum(adj['z'].values())
    best = {False: 0, True: None}
    prev = None
    for v in order:
        hubs = adj[v].get('0', 0) + adj[v].get('1', 0)
        link = adj[v].get(prev, 0) if prev is not None else 0
        # edges counted iff touching V u {z}: edge (prev, v) counted if prev in V or v in V; (v, z) already in const
        gainV = hubs + link - (lam + 1)
        nb = {False: min(x for x in (best[False], (best[True] + link) if best[True] is not None else None) if x is not None),
              True: min(x for x in (best[False] + gainV if best[False] is not None else None,
                                    best[True] + gainV if best[True] is not None else None) if x is not None)}
        best, prev = nb, v
    # last vertex in order adjacent to the next path vertex? all path edges handled; edges to z are in const
    out['inf'] = -(const + min(x for x in best.values() if x is not None))
    return out
