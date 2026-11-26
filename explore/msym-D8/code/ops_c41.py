"""Q_0, Q_1, Q_2 of holomorphic Chern-Simons on C^{4|1} as operators on polynomials in the letters.
Adapted from S2 code/salg.py, ops.py, ops2.py (one odd coordinate theta, four holomorphic directions).

Letters (t, A, k1, k2, k3, k4): t = 0: theta (A = 1); t = 1: c (odd); t = 3: gamma (even).
Code letters: C^A(theta)|_{dxbar=0} = c^A + theta gamma^A,  c = i sfc (S2 convention, eq-res-letters).
Polynomial: dict {(odd_tuple, even_tuple): Fraction}, tuples sorted.

delta C^B(theta') / delta C^A(theta) = delta_A^B (theta - theta')    (S2 eq-fd with one theta)
  =>  delta/delta C^A(theta) = theta d_L/dc^A - d/dgamma^A   (an even operator).
hat Q_0 C^A = (1/2) f_{ABC} C^B C^C,   Q_0 = i hat Q_0,   Q_0 theta = - theta Q_0.
hat Q_1 = int dtheta int d^8z_1 d^8z_2 int_{1>r_1>..>r_4>0} d^4r eps^{abcd} f_{A1 D1 R1} f_{R1 D2 R2} f_{R2 D3 R3} f_{R3 D4 A2}
          d_a C^{D1}(x_1) d_b C^{D2}(x_2) d_c C^{D3}(x_3) d_d C^{D4}(x_4) delta_{A1}(z_1) delta_{A2}(z_2),
          x_j = r_j z_1 + (1 - r_j) z_2.
hat Q_Gamma (theta graphs, O at a pole) = int dtheta int prod_j d^8 z_j int_{D_Gamma} J_Gamma(d) T_Gamma
          delta_{A1}(z_1) delta_{A2}(z_2) delta_{A3}(z_3),
  T_Gamma: arm j (path P_j without O, from the spoke towards N): f_{A_j D R} f_{R D' R'} ... ending in S_j;
           at N: f_{S_1 S_2 S_3};  outputs C^{D_v}(x + X_v) in the order of the vertex labels.
int dtheta theta X = X.
"""
from fractions import Fraction as Fr
from itertools import product, permutations
import sympy as sp

H = 4
TH = (0, 1, 0, 0, 0, 0)
ODD_T = (0, 1)


def is_odd(l):
    return l[0] in ODD_T


def const(x=1):
    return {((), ()): Fr(x)}


def gen(l):
    return {((l,), ()): Fr(1)} if is_odd(l) else {((), (l,)): Fr(1)}


def add(p, q, s=1):
    r = dict(p)
    for m, c in q.items():
        v = r.get(m, 0) + s * c
        if v:
            r[m] = v
        else:
            r.pop(m, None)
    return r


def add_into(r, q, s=1):
    for m, c in q.items():
        v = r.get(m, 0) + s * c
        if v:
            r[m] = v
        else:
            r.pop(m, None)
    return r


def scale(p, s):
    return {m: c * s for m, c in p.items()} if s else {}


def _merge_odd(a, b):
    if not a:
        return b, 1
    if not b:
        return a, 1
    out, i, j, inv = [], 0, 0, 0
    la, lb = len(a), len(b)
    while i < la and j < lb:
        if a[i] < b[j]:
            out.append(a[i]); i += 1
        elif a[i] > b[j]:
            out.append(b[j]); j += 1
            inv += la - i
        else:
            return None
    out.extend(a[i:]); out.extend(b[j:])
    return tuple(out), (-1 if inv & 1 else 1)


def mul(p, q):
    r = {}
    for (o1, e1), c1 in p.items():
        for (o2, e2), c2 in q.items():
            mo = _merge_odd(o1, o2)
            if mo is None:
                continue
            o, s = mo
            key = (o, tuple(sorted(e1 + e2)))
            v = r.get(key, 0) + s * c1 * c2
            if v:
                r[key] = v
            else:
                r.pop(key, None)
    return r


def lderiv(p, l):
    r = {}
    if is_odd(l):
        for (o, e), c in p.items():
            if l in o:
                i = o.index(l)
                key = (o[:i] + o[i + 1:], e)
                r[key] = r.get(key, 0) + (-c if i & 1 else c)
    else:
        for (o, e), c in p.items():
            n = e.count(l)
            if n:
                i = e.index(l)
                key = (o, e[:i] + e[i + 1:])
                r[key] = r.get(key, 0) + n * c
    return {m: c for m, c in r.items() if c}


def letters_of(p):
    s = set()
    for (o, e) in p:
        s.update(o); s.update(e)
    return s


def shift(l, a):
    k = list(l[2:]); k[a] += 1
    return (l[0], l[1]) + tuple(k)


def dpar(p, a):
    r = {}
    for l in letters_of(p):
        if l[0] == 0:
            continue
        r = add(r, mul(gen(shift(l, a)), lderiv(p, l)))
    return r


_CS = {}


def csuper(A, k=(0, 0, 0, 0)):
    key = (A, tuple(k))
    if key not in _CS:
        r = gen((1, A) + tuple(k))
        r = add(r, mul(gen(TH), gen((3, A) + tuple(k))))
        _CS[key] = r
    return _CS[key]


def theta_coeff(p, which):
    want = {'0': (), '1': (TH,)}[which]
    n = len(want)
    r = {}
    for (o, e), c in p.items():
        nth = sum(1 for l in o if l[0] == 0)
        if nth == n and o[:n] == want:
            r[(o[n:], e)] = c
    return r


def THETA(t):
    """theta factor of delta/delta C(theta) on the component of type t (S2 salg.THETA with one theta)."""
    if t == 1:
        return gen(TH)
    if t == 3:
        return const(-1)
    raise ValueError(t)


def sderiv(p, l):
    d = lderiv(p, l)
    return mul(THETA(l[0]), d) if d else {}


def eps4():
    out = {}
    for perm in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        out[perm] = -1 if inv % 2 else 1
    return out


EPS4 = eps4()


def structure_constants(name):
    if name == 'su2':
        f = {}
        for (a, b, c), v in {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}.items():
            f[(a, b, c)] = Fr(v)
        return 3, f
    if name.startswith('so') and name[2:].isdigit():
        N = int(name[2:])
        idx = {}
        for i in range(N):
            for j in range(i + 1, N):
                idx[(i, j)] = len(idx)

        def Lmat(i, j):
            M = [[0] * N for _ in range(N)]
            M[i][j] = 1; M[j][i] = -1
            return M

        def comm(X, Y):
            return [[sum(X[a][k] * Y[k][b] - Y[a][k] * X[k][b] for k in range(N)) for b in range(N)] for a in range(N)]
        f = {}
        gens = {v: Lmat(*k) for k, v in idx.items()}
        for B in gens:
            for C in gens:
                Z = comm(gens[B], gens[C])
                for (i, j), A in idx.items():
                    if Z[i][j]:
                        f[(A, B, C)] = Fr(Z[i][j])
        for (A, B, C), v in f.items():
            assert f.get((B, A, C), 0) == -v and f.get((A, C, B), 0) == -v
        return len(idx), f
    raise ValueError(name)


def L(t, A, k=(0, 0, 0, 0)):
    return gen(({'c': 1, 'g': 3}[t], A) + tuple(k))


# ------------------------------------------------------------------ pentagon moments

def pent_moment(a, b, c, d):
    """int_{1>r1>r2>r3>r4>0} r1^a r2^b r3^c r4^d."""
    return Fr(1, (d + 1) * (c + d + 2) * (b + c + d + 3) * (a + b + c + d + 4))


class Theory:
    def __init__(self, n, f):
        self.n = n
        self.f = f
        self.fby = {}
        for (A, B, C), v in f.items():
            self.fby.setdefault(A, []).append((B, C, v))
        self._q0 = {}
        self._T1 = {}
        self._W1 = {}

    def q0_letter(self, l):
        if l in self._q0:
            return self._q0[l]
        t, A, k = l[0], l[1], l[2:]
        if any(k):
            a = next(i for i in range(H) if k[i])
            kk = list(k); kk[a] -= 1
            r = dpar(self.q0_letter((t, A) + tuple(kk)), a)
        else:
            R = {}
            for (B, C, v) in self.fby.get(A, []):
                R = add(R, scale(mul(csuper(B), csuper(C)), v / 2))
            # hatQ0 (c + th g) = hatQ0 c - th hatQ0 g
            r = {1: theta_coeff(R, '0'), 3: scale(theta_coeff(R, '1'), -1)}[t]
        self._q0[l] = r
        return r

    def Q0(self, p):
        r = {}
        for l in letters_of(p):
            if l[0] == 0:
                continue
            d = lderiv(p, l)
            if d:
                add_into(r, mul(self.q0_letter(l), d))
        return r

    def W1(self, k1, k2):
        """{(m1..m4): Fraction}: d_{z1}^{k1} d_{z2}^{k2} at z = 0 of int d^4r exp(sum_j lam_j . x_j)."""
        key = (k1, k2)
        if key in self._W1:
            return self._W1[key]
        s = sp.symbols('s1:5')
        u = [[sp.Symbol('u%d%d' % (j, a)) for a in range(H)] for j in range(4)]
        expr = sp.Integer(1)
        for a in range(H):
            expr *= sum(s[j] * u[j][a] for j in range(4)) ** k1[a]
            expr *= sum((1 - s[j]) * u[j][a] for j in range(4)) ** k2[a]
        allu = [u[j][a] for j in range(4) for a in range(H)]
        poly = sp.Poly(sp.expand(expr), *allu)
        out = {}
        for mon, coeff in poly.terms():
            cp = sp.Poly(coeff, *s)
            val = Fr(0)
            for (a_, b_, c_, d_), cc in cp.terms():
                val += Fr(int(cc.p), int(cc.q)) * pent_moment(a_, b_, c_, d_)
            if val != 0:
                m = tuple(tuple(mon[H * j + a] for a in range(H)) for j in range(4))
                out[m] = val
        self._W1[key] = out
        return out

    def T1(self, A1, A2, m):
        key = (A1, A2, m)
        if key in self._T1:
            return self._T1[key]
        r = {}
        for (D1, R1, v1) in self.fby.get(A1, []):
            for (D2, R2, v2) in self.fby.get(R1, []):
                for (D3, R3, v3) in self.fby.get(R2, []):
                    for (D4, B, v4) in self.fby.get(R3, []):
                        if B != A2:
                            continue
                        for (a, b, c, d), e in EPS4.items():
                            js = [tuple(m[q][i] + (i == x) for i in range(H)) for q, x in enumerate((a, b, c, d))]
                            term = mul(mul(mul(csuper(D1, js[0]), csuper(D2, js[1])), csuper(D3, js[2])), csuper(D4, js[3]))
                            add_into(r, scale(term, v1 * v2 * v3 * v4 * e))
        self._T1[key] = r
        return r

    def G1(self, A1, A2, k1, k2):
        r = {}
        for m, w in self.W1(k1, k2).items():
            add_into(r, scale(self.T1(A1, A2, m), w))
        return r

    def Q1(self, p):
        r = {}
        ls = [l for l in letters_of(p) if l[0] != 0]
        for l2 in ls:
            d2 = sderiv(p, l2)
            if not d2:
                continue
            for l1 in ls:
                d12 = sderiv(d2, l1)
                if not d12:
                    continue
                G = self.G1(l1[1], l2[1], l1[2:], l2[2:])
                if G:
                    add_into(r, theta_coeff(mul(G, d12), '1'))
        return r


# ------------------------------------------------------------------ two loops

class TwoLoop:
    """data[name] = dict(outputs, arms, K (callable ks -> {m: Fraction}))"""

    def __init__(self, theory, kdata):
        self.th = theory
        self.fby = theory.fby
        self.kd = kdata
        self._T = {}

    def arm_vec(self, A, arm_outputs, ms):
        """vec[S] = sum f_{A D1 R1} f_{R1 D2 R2} ... f_{R D_last S} C^{D1} ... C^{D_last} (superfields)."""
        vec = {A: const(1)}
        for v, mv in zip(arm_outputs, ms):
            new = {}
            for R, p in vec.items():
                for (D, Rn, val) in self.fby.get(R, []):
                    t = scale(mul(p, csuper(D, mv)), val)
                    if Rn in new:
                        add_into(new[Rn], t)
                    else:
                        new[Rn] = t
            vec = {k: v_ for k, v_ in new.items() if v_}
        return vec

    def T(self, name, As, m):
        key = (name, As, m)
        if key in self._T:
            return self._T[key]
        arms = self.kd[name]['arms']           # list of 3 lists of output vertices (from spoke to N)
        outputs = self.kd[name]['outputs']
        pos = {v: i for i, v in enumerate(outputs)}
        vecs = [self.arm_vec(As[j], arms[j], [m[pos[v]] for v in arms[j]]) for j in range(3)]
        res = {}
        for S1, p1 in vecs[0].items():
            for (S2, S3, val) in self.fby.get(S1, []):
                if S2 in vecs[1] and S3 in vecs[2]:
                    add_into(res, scale(mul(mul(p1, vecs[1][S2]), vecs[2][S3]), val))
        self._T[key] = res
        return res

    def G(self, name, As, ks):
        r = {}
        for m, w in self.kd[name]['K'](ks).items():
            t = self.T(name, As, m)
            if t:
                add_into(r, scale(t, w))
        return r

    def Q2(self, name, p):
        r = {}
        ls = [l for l in letters_of(p) if l[0] != 0]
        for l3 in ls:
            d3 = sderiv(p, l3)
            if not d3:
                continue
            for l2 in ls:
                d23 = sderiv(d3, l2)
                if not d23:
                    continue
                for l1 in ls:
                    d123 = sderiv(d23, l1)
                    if not d123:
                        continue
                    G = self.G(name, (l1[1], l2[1], l3[1]), (l1[2:], l2[2:], l3[2:]))
                    if G:
                        add_into(r, theta_coeff(mul(G, d123), '1'))
        return r


# ------------------------------------------------------------------ fast two-loop operator (same result as TwoLoop.Q2)

def _sort_odd(seq):
    """Sorted tuple and sign of the permutation; None if a letter repeats."""
    s = list(seq)
    n = len(s)
    sign = 1
    for i in range(1, n):
        j = i
        while j > 0 and s[j - 1] > s[j]:
            s[j - 1], s[j] = s[j], s[j - 1]
            sign = -sign
            j -= 1
        if j > 0 and s[j - 1] == s[j]:
            return None, 0
    for i in range(1, n):
        if s[i - 1] == s[i]:
            return None, 0
    return tuple(s), sign


class FastTwoLoop(TwoLoop):
    """int dtheta (G_0 + theta G_1)(d_0 + theta d_1) = G_1 d_0 + G_0 d_1 (G_0 is even: six odd letters)."""

    def __init__(self, theory, kdata):
        super().__init__(theory, kdata)
        self._CA = {}
        self._G = {}

    def colour(self, name, As):
        key = (name, As)
        if key in self._CA:
            return self._CA[key]
        arms = self.kd[name]['arms']
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
        self._CA[key] = out
        return out

    def Gparts(self, name, As, ks):
        key = (name, As, ks)
        if key in self._G:
            return self._G[key]
        CA = self.colour(name, As)
        W = self.kd[name]['K'](ks)
        G0, G1 = {}, {}
        for m, w in W.items():
            for D, cf in CA.items():
                cs = [(1, D[i]) + m[i] for i in range(len(D))]
                wc = w * cf
                s0, sg = _sort_odd(cs)
                if s0 is not None:
                    k0 = (s0, ())
                    v = G0.get(k0, 0) + wc * sg
                    if v:
                        G0[k0] = v
                    else:
                        G0.pop(k0, None)
                for q in range(len(D)):
                    s1, sg1 = _sort_odd(cs[:q] + cs[q + 1:])
                    if s1 is None:
                        continue
                    k1 = (s1, ((3, D[q]) + m[q],))
                    v = G1.get(k1, 0) + wc * sg1 * (-1 if q % 2 else 1)
                    if v:
                        G1[k1] = v
                    else:
                        G1.pop(k1, None)
        self._G[key] = (G0, G1)
        return G0, G1

    def Q2(self, name, p):
        r = {}
        ls = [l for l in letters_of(p) if l[0] != 0]
        for l3 in ls:
            d3 = sderiv(p, l3)
            if not d3:
                continue
            for l2 in ls:
                d23 = sderiv(d3, l2)
                if not d23:
                    continue
                for l1 in ls:
                    d123 = sderiv(d23, l1)
                    if not d123:
                        continue
                    d0 = theta_coeff(d123, '0')
                    d1 = theta_coeff(d123, '1')
                    G0, G1 = self.Gparts(name, (l1[1], l2[1], l3[1]), (l1[2:], l2[2:], l3[2:]))
                    if d0:
                        add_into(r, mul(G1, d0))
                    if d1:
                        add_into(r, mul(G0, d1))
        return r
