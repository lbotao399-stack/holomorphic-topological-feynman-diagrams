"""C10 (closed forms of hat Q_1), C11 ({hatQ0, hatQ1} = 0), C12 (hatQ1^2 + sum_G k_G {hatQ0, hatQ_G} = 0).
Own implementation in components (no theta letter), independent of ops_c41.py; theta weights from v3_thetas.py.

Letters: c^A_k (odd), gamma^A_k (even), k in N^4 (holomorphic derivatives, not divided by k!).
Superfield C = c + theta gamma.  With delta_j := delta/delta C^{A_j}(theta) = theta d/dc_j - d/dgamma_j (claim C1):
  delta_1 delta_2          = -theta (d_c1 d_g2 + d_g1 d_c2) + d_g1 d_g2,
  delta_1 delta_2 delta_3  =  theta (d_c1 d_g2 d_g3 + d_g1 d_c2 d_g3 + d_g1 d_g2 d_c3) - d_g1 d_g2 d_g3,
  F(theta) = F_0 + theta F_1 (product of odd superfields; F_1 = sum_q (-1)^{q-1} [q-th c -> gamma]),
  int dtheta F delta_1 delta_2         = -F_0 (d_c1 d_g2 + d_g1 d_c2) + F_1 d_g1 d_g2,
  int dtheta F delta_1 delta_2 delta_3 =  F_0 (d_c1 d_g2 d_g3 + d_g1 d_c2 d_g3 + d_g1 d_g2 d_c3) - F_1 d_g1 d_g2 d_g3.
hat Q_0 C^A = (1/2) f_ABC C^B C^C, odd, Q_0 theta = -theta Q_0  =>  hatQ0 c^A = (1/2) f_ABC c^B c^C,
  hatQ0 gamma^A = f_ABC c^B gamma^C ; commutes with d_a.
Pentagon (C10): F = sum f_{A1 D1 R1} f_{R1 D2 R2} f_{R2 D3 R3} f_{R3 D4 A2} eps^{abcd} sum_m W1(k1,k2;m)
  d^{m_1+e_a} C^{D1} d^{m_2+e_b} C^{D2} d^{m_3+e_c} C^{D3} d^{m_4+e_d} C^{D4},
  W1(k1,k2;m) = [u^m] int_{1>r1>..>r4>0} prod_a (sum_q r_q u_qa)^{k1_a} (sum_q (1-r_q) u_qa)^{k2_a}.
Theta (C10): F = sum_D T_G(A;D) sum_m W_G(k;m) prod_{v in vertex order} d^{m_v} C^{D_v}, T_G = arms f_{A_j D R}...f_{.. D S_j},
  then f_{S1 S2 S3}; W_G = Taylor weights of K_G(lam; z) (own, v3_thetas.py).
"""
import sys, time, pickle, itertools
from fractions import Fraction as Fr
from math import comb

H = 4
LOG = None


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if LOG:
        LOG.write(s + '\n'); LOG.flush()


# ------------------------------------------------------------------ letters and monomials
def mk(t, A, k):
    return ((((k[0] * 16 + k[1]) * 16 + k[2]) * 16 + k[3]) * 64 + A) * 2 + t


def ltyp(L):
    return L & 1


def lcol(L):
    return (L >> 1) & 63


def lk(L):
    x = L >> 7
    return (x >> 12, (x >> 8) & 15, (x >> 4) & 15, x & 15)


def sort_sign(seq):
    s = list(seq)
    sign = 1
    for i in range(1, len(s)):
        j = i
        while j > 0 and s[j - 1] > s[j]:
            s[j - 1], s[j] = s[j], s[j - 1]
            sign = -sign
            j -= 1
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            return None, 0
    return tuple(s), sign


def merge_sign(a, b):
    """sorted tuples a, b of odd letters: (sorted a+b, sign of the shuffle) or (None, 0)."""
    if not a:
        return b, 1
    if not b:
        return a, 1
    out = []
    i = j = inv = 0
    la, lb = len(a), len(b)
    while i < la and j < lb:
        x, y = a[i], b[j]
        if x < y:
            out.append(x); i += 1
        elif x > y:
            out.append(y); j += 1; inv += la - i
        else:
            return None, 0
    out.extend(a[i:]); out.extend(b[j:])
    return tuple(out), (-1 if inv & 1 else 1)


def addto(r, key, v):
    w = r.get(key)
    if w is None:
        r[key] = v
    else:
        w = w + v
        if w:
            r[key] = w
        else:
            del r[key]


def mul_left(F, O2, E2, coef, out):
    """out += coef * F . (O2, E2)."""
    for (O1, E1), c in F.items():
        O, s = merge_sign(O1, O2)
        if O is None:
            continue
        E = tuple(sorted(E1 + E2)) if E1 and E2 else (E1 or E2)
        addto(out, (O, E), coef * c * s)


def distinct_even(E):
    """(letter, multiplicity, E with one copy removed) for the distinct letters of E."""
    out = []
    seen = set()
    for i, e in enumerate(E):
        if e in seen:
            continue
        seen.add(e)
        out.append((e, E.count(e), E[:i] + E[i + 1:]))
    return out


def poly_mono(O, E, c=1):
    O2, s = sort_sign(O)
    return {(O2, tuple(sorted(E))): Fr(c * s)} if O2 is not None else {}


def pmul(p, q):
    r = {}
    for (O1, E1), c1 in p.items():
        for (O2, E2), c2 in q.items():
            O, s = merge_sign(O1, O2)
            if O is None:
                continue
            addto(r, (O, tuple(sorted(E1 + E2))), c1 * c2 * s)
    return r


def padd(p, q, s=1):
    r = dict(p)
    for k, v in q.items():
        addto(r, k, s * v)
    return r


# ------------------------------------------------------------------ structure constants
def su2():
    f = {}
    for (a, b, c) in itertools.permutations(range(3)):
        f[(a, b, c)] = Fr(perm_sign((a, b, c)))
    return 3, f


def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def so(N):
    """basis T_{ij} = E_ij - E_ji (i<j); [T_B, T_C] = sum_A f_ABC T_A; f totally antisymmetric (asserted)."""
    idx = [(i, j) for i in range(N) for j in range(i + 1, N)]
    import numpy as np
    def T(i, j):
        M = np.zeros((N, N), dtype=int); M[i, j] = 1; M[j, i] = -1
        return M
    gens = [T(*p) for p in idx]
    f = {}
    for B, TB in enumerate(gens):
        for C, TC in enumerate(gens):
            Z = TB @ TC - TC @ TB
            for A, (i, j) in enumerate(idx):
                if Z[i, j]:
                    f[(A, B, C)] = Fr(int(Z[i, j]))
    for (A, B, C), v in f.items():
        assert f.get((B, A, C), 0) == -v and f.get((A, C, B), 0) == -v and f.get((C, B, A), 0) == -v
    # Jacobi
    n = len(idx)
    for A, B, C, D in itertools.product(range(n), repeat=4):
        s = sum(f.get((A, B, E), 0) * f.get((E, C, D), 0) + f.get((A, C, E), 0) * f.get((E, D, B), 0)
                + f.get((A, D, E), 0) * f.get((E, B, C), 0) for E in range(n))
        assert s == 0
    return n, f


# ------------------------------------------------------------------ theory
EPS4 = [(p, perm_sign(p)) for p in itertools.permutations(range(4))]


def pent_moment(a, b, c, d):
    return Fr(1, (d + 1) * (c + d + 2) * (b + c + d + 3) * (a + b + c + d + 4))


class Theory:
    def __init__(self, n, f, theta_weights=None):
        self.n, self.f = n, f
        self.fby = {}
        for (A, B, C), v in f.items():
            self.fby.setdefault(A, []).append((B, C, v))
        self._q0 = {}
        self._W1 = {}
        self._Fp = {}
        self._pc = {}
        self._tc = {}
        self._Ft = {}
        self.tw = theta_weights

    # ---- Q0
    def q0_letter(self, L):
        if L in self._q0:
            return self._q0[L]
        A, k, t = lcol(L), lk(L), ltyp(L)
        r = {}
        splits = []
        for kp in itertools.product(*[range(x + 1) for x in k]):
            km = tuple(k[a] - kp[a] for a in range(H))
            b = 1
            for a in range(H):
                b *= comb(k[a], kp[a])
            splits.append((kp, km, b))
        for (B, C, v) in self.fby.get(A, []):
            for kp, km, b in splits:
                if t == 0:
                    O, s = sort_sign((mk(0, B, kp), mk(0, C, km)))
                    if O is not None:
                        addto(r, (O, ()), v * b * s / 2)
                else:
                    addto(r, ((mk(0, B, kp),), (mk(1, C, km),)), v * b)
        self._q0[L] = r
        return r

    def Q0(self, p):
        out = {}
        for (O, E), coef in p.items():
            for i, o in enumerate(O):
                s = -1 if i & 1 else 1
                pre, post = O[:i], O[i + 1:]
                for (Oq, Eq), cq in self.q0_letter(o).items():
                    On, sg = sort_sign(pre + Oq + post)
                    if On is None:
                        continue
                    En = tuple(sorted(E + Eq)) if Eq else E
                    addto(out, (On, En), coef * s * sg * cq)
            sp_ = -1 if len(O) & 1 else 1
            for e, nmult, Er in distinct_even(E):
                for (Oq, Eq), cq in self.q0_letter(e).items():
                    On, sg = merge_sign(O, Oq)
                    if On is None:
                        continue
                    En = tuple(sorted(Er + Eq))
                    addto(out, (On, En), coef * sp_ * nmult * sg * cq)
        return out

    # ---- pentagon
    def W1(self, k1, k2):
        key = (k1, k2)
        if key in self._W1:
            return self._W1[key]
        # polynomial in r (4 exps) and u (16 exps): dict ((r exps), (u exps)) -> int
        poly = {((0,) * 4, (0,) * 16): 1}

        def mulf(poly, a, kind):
            new = {}
            for (re, ue), c in poly.items():
                for q in range(4):
                    ue2 = list(ue); ue2[4 * q + a] += 1; ue2 = tuple(ue2)
                    re2 = list(re); re2[q] += 1; re2 = tuple(re2)
                    if kind == 1:            # r_q u_qa
                        new[(re2, ue2)] = new.get((re2, ue2), 0) + c
                    else:                    # (1 - r_q) u_qa
                        new[(re, ue2)] = new.get((re, ue2), 0) + c
                        new[(re2, ue2)] = new.get((re2, ue2), 0) - c
            return {k_: v for k_, v in new.items() if v}
        for a in range(H):
            for _ in range(k1[a]):
                poly = mulf(poly, a, 1)
            for _ in range(k2[a]):
                poly = mulf(poly, a, 2)
        out = {}
        for (re, ue), c in poly.items():
            m = tuple(tuple(ue[4 * q + a] for a in range(H)) for q in range(4))
            out[m] = out.get(m, 0) + c * pent_moment(*re)
        out = {m: v for m, v in out.items() if v}
        self._W1[key] = out
        return out

    def pent_colour(self, A1, A2):
        key = (A1, A2)
        if key in self._pc:
            return self._pc[key]
        r = {}
        for (D1, R1, v1) in self.fby.get(A1, []):
            for (D2, R2, v2) in self.fby.get(R1, []):
                for (D3, R3, v3) in self.fby.get(R2, []):
                    for (D4, B, v4) in self.fby.get(R3, []):
                        if B == A2:
                            k = (D1, D2, D3, D4)
                            r[k] = r.get(k, 0) + v1 * v2 * v3 * v4
        r = {k: v for k, v in r.items() if v}
        self._pc[key] = r
        return r

    def Fpent(self, A1, A2, k1, k2):
        key = (A1, A2, k1, k2)
        if key in self._Fp:
            return self._Fp[key]
        F0, F1 = {}, {}
        W = self.W1(k1, k2)
        col = self.pent_colour(A1, A2)
        for D, cD in col.items():
            for m, w in W.items():
                for perm, se in EPS4:
                    ks = [tuple(m[q][a] + (1 if a == perm[q] else 0) for a in range(H)) for q in range(4)]
                    cs = [mk(0, D[q], ks[q]) for q in range(4)]
                    O, s = sort_sign(cs)
                    wc = w * cD * se
                    if O is not None:
                        addto(F0, (O, ()), wc * s)
                    for q in range(4):
                        O1, s1 = sort_sign(cs[:q] + cs[q + 1:])
                        if O1 is None:
                            continue
                        addto(F1, (O1, (mk(1, D[q], ks[q]),)), wc * s1 * (-1 if q & 1 else 1))
        self._Fp[key] = (F0, F1)
        return F0, F1

    def Q1(self, p):
        out = {}
        for (O, E), coef in p.items():
            dE = distinct_even(E)
            # F1 d_g1 d_g2
            for g2, n2, E2 in dE:
                for g1, n1, E12 in distinct_even(E2):
                    F0, F1 = self.Fpent(lcol(g1), lcol(g2), lk(g1), lk(g2))
                    mul_left(F1, O, E12, coef * n1 * n2, out)
            # -F0 d_c1 d_g2
            for g2, n2, E2 in dE:
                for i, c1 in enumerate(O):
                    s = -1 if i & 1 else 1
                    F0, F1 = self.Fpent(lcol(c1), lcol(g2), lk(c1), lk(g2))
                    mul_left(F0, O[:i] + O[i + 1:], E2, -coef * n2 * s, out)
            # -F0 d_g1 d_c2
            for i, c2 in enumerate(O):
                s = -1 if i & 1 else 1
                for g1, n1, E1 in dE:
                    F0, F1 = self.Fpent(lcol(g1), lcol(c2), lk(g1), lk(c2))
                    mul_left(F0, O[:i] + O[i + 1:], E1, -coef * s * n1, out)
        return out

    # ---- thetas
    def theta_colour(self, name, As):
        key = (name, As)
        if key in self._tc:
            return self._tc[key]
        arms = self.tw[name]['arms']
        vecs = []
        for j in range(3):
            vec = {As[j]: {(): Fr(1)}}
            for _v in arms[j]:
                new = {}
                for R, dd in vec.items():
                    for (D, Rn, val) in self.fby.get(R, []):
                        tgt = new.setdefault(Rn, {})
                        for Ds, c in dd.items():
                            k = Ds + (D,)
                            tgt[k] = tgt.get(k, 0) + c * val
                vec = new
            vecs.append(vec)
        out = {}
        for S1, d1 in vecs[0].items():
            for (S2, S3, val) in self.fby.get(S1, []):
                if S2 in vecs[1] and S3 in vecs[2]:
                    for a, ca in d1.items():
                        for b, cb in vecs[1][S2].items():
                            for c, cc in vecs[2][S3].items():
                                k = a + b + c
                                out[k] = out.get(k, 0) + ca * cb * cc * val
        out = {k: v for k, v in out.items() if v}
        self._tc[key] = out
        return out

    def theta_W(self, name, ks):
        w = self.tw[name]
        nz = [(j, a) for j in range(3) for a in range(H) if ks[j][a]]
        if not nz:
            return w['W0']
        if len(nz) == 1 and ks[nz[0][0]][nz[0][1]] == 1:
            return w['W1'][nz[0]]
        raise NotImplementedError(ks)

    def Ftheta(self, name, As, ks, part):
        key = (name, As, ks, part)
        if key in self._Ft:
            return self._Ft[key]
        F = {}
        col = self.theta_colour(name, As)
        if col:
            W = self.theta_W(name, ks)
            for D, cD in col.items():
                for m, w in W.items():
                    cs = [mk(0, D[q], m[q]) for q in range(6)]
                    wc = w * cD
                    if part == 0:
                        O, s = sort_sign(cs)
                        if O is not None:
                            addto(F, (O, ()), wc * s)
                    else:
                        for q in range(6):
                            O1, s1 = sort_sign(cs[:q] + cs[q + 1:])
                            if O1 is None:
                                continue
                            addto(F, (O1, (mk(1, D[q], m[q]),)), wc * s1 * (-1 if q & 1 else 1))
        self._Ft[key] = F
        return F

    def Qtheta(self, name, p):
        out = {}
        for (O, E), coef in p.items():
            dE = distinct_even(E)
            # F0 d_c1 d_g2 d_g3
            for g3, n3, E3 in dE:
                for g2, n2, E23 in distinct_even(E3):
                    for i, c1 in enumerate(O):
                        s = -1 if i & 1 else 1
                        F = self.Ftheta(name, (lcol(c1), lcol(g2), lcol(g3)), (lk(c1), lk(g2), lk(g3)), 0)
                        mul_left(F, O[:i] + O[i + 1:], E23, coef * n3 * n2 * s, out)
            # F0 d_g1 d_c2 d_g3
            for g3, n3, E3 in dE:
                for i, c2 in enumerate(O):
                    s = -1 if i & 1 else 1
                    for g1, n1, E13 in distinct_even(E3):
                        F = self.Ftheta(name, (lcol(g1), lcol(c2), lcol(g3)), (lk(g1), lk(c2), lk(g3)), 0)
                        mul_left(F, O[:i] + O[i + 1:], E13, coef * n3 * s * n1, out)
            # F0 d_g1 d_g2 d_c3
            for i, c3 in enumerate(O):
                s = -1 if i & 1 else 1
                for g2, n2, E2 in dE:
                    for g1, n1, E12 in distinct_even(E2):
                        F = self.Ftheta(name, (lcol(g1), lcol(g2), lcol(c3)), (lk(g1), lk(g2), lk(c3)), 0)
                        mul_left(F, O[:i] + O[i + 1:], E12, coef * s * n2 * n1, out)
            # -F1 d_g1 d_g2 d_g3
            for g3, n3, E3 in dE:
                for g2, n2, E23 in distinct_even(E3):
                    for g1, n1, E123 in distinct_even(E23):
                        F = self.Ftheta(name, (lcol(g1), lcol(g2), lcol(g3)), (lk(g1), lk(g2), lk(g3)), 1)
                        mul_left(F, O, E123, -coef * n3 * n2 * n1, out)
        return out


def letter(t, A, k=(0, 0, 0, 0)):
    return mk({'c': 0, 'g': 1}[t], A, k)


def op(*ls):
    O = [L for L in ls if ltyp(L) == 0]
    E = [L for L in ls if ltyp(L) == 1]
    # the operator is the product in the order given
    seq_odd = [L for L in ls if ltyp(L) == 0]
    return poly_mono(seq_odd, E)


def e(a, n=1):
    return tuple(n if i == a else 0 for i in range(H))


def name_letter(L):
    t = 'c' if ltyp(L) == 0 else 'g'
    k = lk(L)
    d = ''.join('d%d' % (a + 1) * k[a] for a in range(H))
    return d + t + str(lcol(L))
