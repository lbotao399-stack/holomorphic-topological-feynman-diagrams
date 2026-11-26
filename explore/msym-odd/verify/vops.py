"""Independent two-loop operators Qhat_Gamma and the closure test {Q0hat, sum_G k_G Qhat_G} = 0.

Qhat_Gamma O = int d^k theta  sum_{l1,l2,l3}  G_Gamma(A_1,A_2,A_3; j_1,j_2,j_3)  delta_{l1} delta_{l2} delta_{l3} O,
G_Gamma = sum_m W_Gamma(j)[m] T_Gamma(A; m),
W_Gamma(j)[m] = coefficient of prod_v lam_v^{m_v} in  d_{z_1}^{j_1} d_{z_2}^{j_2} d_{z_3}^{j_3} K_Gamma(lam; z)|_{z=0} / pi,
K_Gamma(lam; z) = int s (D_Gamma/2) prod_P w_P^{H+1/2-n_P} exp(sum_v lam_v . X_v(z)) (density verified
pointwise in kernels_verify.py), X_v = sigma z_{P_j} + (1 - sigma) sum_P w_P z_P,
so d_{z_j}^a X_v^a = c_{vj} = sigma^P_k [P = P_j] + (1 - sigma^P_k) w_{P_j}.
T_Gamma: colour tensor with output superfields; S2 eq-res-TA / eq-res-TB for graph data 'S2A', 'S2B',
arm convention for theta_rootN data (explorer's opsk.py docstring):
  arm_P[R] = f_{A_P D_1 R_1} C^{D_1}(v_1) f_{R_1 D_2 R_2} C^{D_2}(v_2) ...,  T = f_{R_A R_B R_C} arm_A arm_B arm_C.
"""
from fractions import Fraction as Fr
from itertools import permutations
import random

from vkern import Theta, theta_rootN, S2_A, S2_B, dirichlet
from vkern import perm_sign as psign


# ---------------------------------------------------------------- small multivariate polynomials
def padd(a, b, s=1):
    r = dict(a)
    for k, v in b.items():
        x = r.get(k, 0) + s * v
        if x:
            r[k] = x
        else:
            r.pop(k, None)
    return r


def pmul(a, b):
    r = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            x = r.get(k, 0) + va * vb
            if x:
                r[k] = x
            else:
                r.pop(k, None)
    return r


class KernelOp:
    def __init__(self, kind, H, pqr=None, mutate=None):
        self.H = H
        if kind == 'S2A':
            nv, ed, pa = S2_A; self.spoke = [1, 0, 2]; self.outv = [0, 1, 3, 4]
        elif kind == 'S2B':
            nv, ed, pa = S2_B; self.spoke = [0, 1, 2]; self.outv = [0, 1, 3, 4]
        else:
            nv, ed, pa = theta_rootN(*pqr); self.spoke = [0, 1, 2]
            self.outv = None
        self.TH = Theta(nv, ed, pa, H)
        TH = self.TH
        sg = TH.check(npts=4, seed=99, zero_Np=True)
        assert len(sg) == 1
        self.s = sg.pop()
        # spoke j sits on the first edge of path spoke[j]; its displacement in the direction from N is +z_j
        # (S2 data: z_1, z_2 on edges oriented away from O, -z_3 on (5->O) resp. (b->O), oriented toward O)
        if kind in ('S2A', 'S2B'):
            assert [TH.sg[P][0] for P in range(3)] == [1, 1, -1]
        if self.outv is None:
            self.outv = [v for P in range(3) for v in TH.inter[P]]
        assert sorted(self.outv) == sorted(v for P in range(3) for v in TH.inter[P])
        self.nout = len(self.outv)
        self.pos = {}                          # output vertex -> (P, k) (k 1-based along P from N)
        for P in range(3):
            for k, v in enumerate(TH.inter[P], start=1):
                self.pos[v] = (P, k)
        self.arms = [[self.outv.index(v) for v in TH.inter[P]] for P in range(3)]
        nl = self.nout * H
        self.iw = nl
        self.ix = [nl + 3 + sum(TH.nP[:P]) for P in range(3)]
        self.nvar = nl + 3 + sum(TH.nP)
        self.mutate = mutate
        self._W = {}
        # D_Gamma(lam) as a polynomial in the output lam (lam_{N'} = 0)
        A, B, C = TH.inter
        cols = A + B + C
        rows = []
        for i in range(H):
            rows.append([('l', v, i) for v in A] + [('l', v, i) for v in B] + [None] * len(C))
        for i in range(H):
            rows.append([None] * len(A) + [('l', v, i) for v in B] + [('l', v, i) for v in C])
        n = len(rows)
        D = {}
        for perm in permutations(range(n)):
            mon = [0] * self.nvar
            ok = True
            for i in range(n):
                ent = rows[i][perm[i]]
                if ent is None:
                    ok = False; break
                mon[self.lv(ent[1], ent[2])] += 1
            if ok:
                D = padd(D, {tuple(mon): Fr(psign(perm))})
        self.D = D
        # c_{vj}
        self.c = {}
        for v in self.outv:
            P, k = self.pos[v]
            for j in range(3):
                Pj = self.spoke[j]
                poly = {}
                # sigma = sum_{i > k} x^P_i  (1-based i = k+1..n_P -> 0-based k..n_P-1)
                sig = {}
                for i0 in range(k, TH.nP[P]):
                    mon = [0] * self.nvar; mon[self.ix[P] + i0] = 1
                    sig = padd(sig, {tuple(mon): Fr(1)})
                wj = [0] * self.nvar; wj[self.iw + Pj] = 1
                one = {tuple([0] * self.nvar): Fr(1)}
                wpoly = {tuple(wj): Fr(1)}
                if P == Pj:
                    poly = padd(poly, sig)
                poly = padd(poly, pmul(padd(one, sig, -1), wpoly))
                self.c[(v, j)] = poly

    def lv(self, v, i):
        return self.outv.index(v) * self.H + i

    def weights(self, js):
        """js: three H-tuples (jets of the contracted letters 1, 2, 3) -> {m: Fraction}, m = tuple of
        output jets (H-tuples) in the order of self.outv."""
        if js in self._W:
            return self._W[js]
        H, TH = self.H, self.TH
        expr = {k: v * Fr(self.s, 2) for k, v in self.D.items()}
        for j in range(3):
            for a in range(H):
                if js[j][a]:
                    lin = {}
                    for v in self.outv:
                        mon = [0] * self.nvar; mon[self.lv(v, a)] = 1
                        lin = padd(lin, pmul({tuple(mon): Fr(1)}, self.c[(v, j)]))
                    for _ in range(js[j][a]):
                        expr = pmul(expr, lin)
        out = {}
        nl = self.nout * H
        for mon, cf in expr.items():
            q, kk = dirichlet([TH.alpha[P] + mon[self.iw + P] for P in range(3)])
            assert kk == 2
            for P in range(3):
                q *= dirichlet(list(mon[self.ix[P]:self.ix[P] + TH.nP[P]]))[0]
            m = tuple(tuple(mon[o * H:(o + 1) * H]) for o in range(self.nout))
            out[m] = out.get(m, 0) + q * cf
        out = {m: v for m, v in out.items() if v}
        if self.mutate == 'double_z' and any(any(x) for x in js):
            out = {m: 2 * v for m, v in out.items()}
        self._W[js] = out
        return out


class TwoLoop:
    def __init__(self, alg, q0, ops, tensor):
        """ops: {name: KernelOp}; tensor: {name: 'S2A' | 'S2B' | 'arm'}."""
        self.S, self.q0 = alg, q0
        self.f, self.fA = q0.f, q0.fA
        self.ops, self.tensor = ops, tensor
        self._col = {}
        self._G = {}
        self._prod = {}

    # colour tensors: {tuple of output colours D (in output order): coefficient}
    def colour(self, name, As):
        key = (name, As)
        if key in self._col:
            return self._col[key]
        f, fA = self.f, self.fA
        A1, A2, A3 = As
        out = {}
        kind = self.tensor[name]
        if kind == 'S2A':
            for (D1, R1, v1) in fA.get(A1, []):
                for (D2, R2, v2) in fA.get(R1, []):
                    for (B, R3, v3) in fA.get(R2, []):
                        if B != A2:
                            continue
                        for (D4, R4, v4) in fA.get(R3, []):
                            for (D5, B5, v5) in fA.get(R4, []):
                                if B5 != A3:
                                    continue
                                k = (D1, D2, D4, D5)
                                out[k] = out.get(k, 0) + v1 * v2 * v3 * v4 * v5
        elif kind == 'S2B':
            for (Du, Ru, v1) in fA.get(A1, []):
                for (Dw, Rw, v2) in fA.get(A2, []):
                    for (Rw_, Rq, v3) in fA.get(Ru, []):
                        if Rw_ != Rw:
                            continue
                        for (Da, Rp, v4) in fA.get(Rq, []):
                            for (Db, B5, v5) in fA.get(Rp, []):
                                if B5 != A3:
                                    continue
                                k = (Du, Dw, Da, Db)
                                out[k] = out.get(k, 0) + v1 * v2 * v3 * v4 * v5
        else:
            op = self.ops[name]
            arms = []
            for P in range(3):
                vec = {((), As[P]): Fr(1)}
                for _ in op.arms[P]:
                    nv = {}
                    for (Ds, R), c in vec.items():
                        for (D, Rn, v) in fA.get(R, []):
                            k = (Ds + (D,), Rn)
                            nv[k] = nv.get(k, 0) + c * v
                    vec = nv
                arms.append(vec)
            for (DA, RA), ca in arms[0].items():
                for (DB, RB), cb in arms[1].items():
                    for (DC, RC), cc in arms[2].items():
                        v = f.get((RA, RB, RC))
                        if v:
                            k = DA + DB + DC
                            out[k] = out.get(k, 0) + ca * cb * cc * v
        out = {k: v for k, v in out.items() if v}
        self._col[key] = out
        return out

    def prod(self, Ds, ms, d):
        """theta-degree-d part of C^{D_1}(m_1) ... C^{D_r}(m_r)."""
        key = (Ds, ms, d)
        if key in self._prod:
            return self._prod[key]
        S = self.S
        p = {((), ()): Fr(1)}
        for D, m in zip(Ds, ms):
            q = S.csuper(D, m)
            r = {}
            for (o1, e1), c1 in p.items():
                n1 = S.theta_part(o1)
                for (o2, e2), c2 in q.items():
                    if n1 + S.theta_part(o2) > d:
                        continue
                    mg = S.merge(o1, o2)
                    if mg is None:
                        continue
                    o, sg = mg
                    k = (o, tuple(sorted(e1 + e2)))
                    r[k] = r.get(k, 0) + sg * c1 * c2
            p = {k: v for k, v in r.items() if v}
        p = {k: v for k, v in p.items() if S.theta_part(k[0]) == d}
        self._prod[key] = p
        return p

    def Gdec(self, name, As, js, d):
        key = (name, As, js, d)
        if key in self._G:
            return self._G[key]
        S = self.S
        W = self.ops[name].weights(js)
        col = self.colour(name, As)
        G = {}
        if W and col:
            for m, w in W.items():
                for Ds, c in col.items():
                    G = S.add(G, S.scale(self.prod(Ds, m, d), w * c))
        res = S.decompose(G)
        self._G[key] = res
        return res

    def top(self, Gd, X):
        S = self.S
        k = S.k
        full = S.full
        r = {}
        for (o, e), c in X.items():
            n = S.theta_part(o)
            Sp = o[:n]
            Sc = tuple(i for i in full if i not in Sp)
            g = Gd.get(Sc)
            if not g:
                continue
            s1 = psign(Sc + Sp)
            rest = o[n:]
            for (go, ge), gc in g.items():
                mg = S.merge(go, rest)
                if mg is None:
                    continue
                oo, s2 = mg
                sg = s1 * s2 * (-1 if (len(go) * n) & 1 else 1)
                key = (oo, tuple(sorted(ge + e)))
                v = r.get(key, 0) + sg * gc * c
                if v:
                    r[key] = v
                else:
                    r.pop(key, None)
        return r

    def Q(self, name, p):
        S = self.S
        ls = sorted(l for l in S.letters(p) if l >= 2 * S.k)
        r = {}
        for l3 in ls:
            S3_, A3, j3 = S.info[l3]
            d3 = S.sder(p, l3)
            if not d3:
                continue
            for l2 in ls:
                S2_, A2, j2 = S.info[l2]
                d23 = S.sder(d3, l2)
                if not d23:
                    continue
                for l1 in ls:
                    S1_, A1, j1 = S.info[l1]
                    d = len(S1_) + len(S2_) + len(S3_) - 2 * S.k
                    if d < 0:
                        continue
                    d123 = S.sder(d23, l1)
                    if not d123:
                        continue
                    Gd = self.Gdec(name, (A1, A2, A3), (j1, j2, j3), d)
                    if Gd:
                        r = S.add(r, self.top(Gd, d123))
        return r

    def anti(self, name, p):
        S = self.S
        return S.add(self.q0(self.Q(name, p)), self.Q(name, self.q0(p)))


# ---------------------------------------------------------------- exact rank / null space (few columns)
def rank_null(rows, ncol):
    basis = []      # list of (pivot, row) in reduced form
    for row in rows:
        r = list(row)
        for piv, b in basis:
            if r[piv]:
                f = r[piv] / b[piv]
                r = [x - f * y for x, y in zip(r, b)]
        nz = [i for i in range(ncol) if r[i]]
        if nz:
            piv = nz[0]
            # reduce existing rows
            nb = []
            for p2, b in basis:
                if b[piv]:
                    f = b[piv] / r[piv]
                    b = [x - f * y for x, y in zip(b, r)]
                nb.append((p2, b))
            basis = nb + [(piv, r)]
    rank = len(basis)
    pivs = {p for p, _ in basis}
    free = [i for i in range(ncol) if i not in pivs]
    null = []
    for fcol in free:
        v = [Fr(0)] * ncol
        v[fcol] = Fr(1)
        for p, b in basis:
            v[p] = -b[fcol] / b[p]
        null.append(v)
    return rank, null
