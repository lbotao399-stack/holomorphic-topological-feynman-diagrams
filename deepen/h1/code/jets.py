"""Translation-invariant local functionals of commutative CS on R x C^2 in jet coordinates.

Generators (jet coordinates of the superfield A^a = sum_T A^a_T dx_T):
    gen = (a, T, al)
    a  : colour index 0..dim g - 1
    T  : bitmask of form directions, bit0 = dtau, bit1 = dzb1, bit2 = dzb2  (dx_T ordered tau < zb1 < zb2)
    al : 5-tuple of derivative orders in directions (tau, z1, z2, zb1, zb2)
    parity(gen) = (1 - |T|) mod 2 = ghost number mod 2  (A^a_T dx_T odd).
M = free supercommutative algebra on generators (reduced: no constant term).
Q  : odd derivation, Q(gen) = d^al [ coefficient of dx_T in  dA^a + (1/2) f^a_bc A^b A^c ],
     dA = sum_{u in tau,zb1,zb2} dx_u d_u A.
D_u: even derivation (total derivative), D_u (a,T,al) = (a,T,al+e_u).
eta_v: odd derivation induced by A -> iota_{d/dx_v} A.
V  = M / sum_u D_u M  (translation-invariant local functionals  int dvol m(jets)).
Exact integer arithmetic throughout (Python int).
"""
from itertools import combinations
from functools import lru_cache

FORMDIR = (0, 3, 4)          # direction index u of the form directions tau, zb1, zb2
BIT_OF_U = {0: 1, 3: 2, 4: 4}
NDIR = 5


def popcount(T):
    return bin(T).count("1")


def gpar(gen):
    return (1 - popcount(gen[1])) % 2


def ghost(gen):
    return 1 - popcount(gen[1])


# ---------------------------------------------------------------- monomials
def mono_mul(m1, m2):
    """product of sorted monomials; returns (sign, monomial) or (0, None)."""
    if not m1:
        return 1, m2
    if not m2:
        return 1, m1
    seq = list(m1) + list(m2)
    # sign of sorting, counting only odd-odd transpositions (stable insertion sort)
    sign = 1
    arr = []
    for x in seq:
        px = gpar(x)
        j = len(arr)
        while j > 0 and arr[j - 1] > x:
            if px and gpar(arr[j - 1]):
                sign = -sign
            j -= 1
        arr.insert(j, x)
    for i in range(len(arr) - 1):
        if arr[i] == arr[i + 1] and gpar(arr[i]):
            return 0, None
    return sign, tuple(arr)


def mono_par(m):
    return sum(gpar(x) for x in m) % 2


def padd(P, m, c):
    if c:
        v = P.get(m, 0) + c
        if v:
            P[m] = v
        else:
            P.pop(m, None)


def pmul(P1, P2):
    R = {}
    for m1, c1 in P1.items():
        for m2, c2 in P2.items():
            s, m = mono_mul(m1, m2)
            if s:
                padd(R, m, s * c1 * c2)
    return R


# ---------------------------------------------------------------- derivations
def derivation(P, act, odd):
    """apply derivation defined on generators by act(gen) -> poly; odd = True for odd derivations."""
    R = {}
    for m, c in P.items():
        pre = 1
        for i, x in enumerate(m):
            img = act(x)
            if img:
                left = m[:i]
                right = m[i + 1:]
                sg = pre if odd else 1
                for mi, ci in img.items():
                    s1, t = mono_mul(left, mi)
                    if not s1:
                        continue
                    s2, t2 = mono_mul(t, right)
                    if not s2:
                        continue
                    padd(R, t2, sg * s1 * s2 * c * ci)
            if gpar(x):
                pre = -pre
    return R


def D_gen(u):
    def act(g):
        al = list(g[2]); al[u] += 1
        return {((g[0], g[1], tuple(al)),): 1}
    return act


def D(P, u):
    return derivation(P, D_gen(u), False)


def Dpow(P, al):
    for u in range(NDIR):
        for _ in range(al[u]):
            P = D(P, u)
    return P


# ---------------------------------------------------------------- forms
def form_sign(T1, T2):
    """dx_T1 dx_T2 = sign dx_{T1|T2} (0 if overlap); bits ordered tau < zb1 < zb2."""
    if T1 & T2:
        return 0
    s = 0
    for i in range(3):
        if T1 >> i & 1:
            for j in range(i):
                if T2 >> j & 1:
                    s += 1
    return -1 if s % 2 else 1


def R_mul(X, Y):
    """X, Y: dict (mono, T) -> c, element = c * mono * dx_T.  (m1 dx_T1)(m2 dx_T2) = (-1)^{|T1| par(m2)} m1 m2 dx_T1 dx_T2."""
    Z = {}
    for (m1, T1), c1 in X.items():
        for (m2, T2), c2 in Y.items():
            fs = form_sign(T1, T2)
            if not fs:
                continue
            s, m = mono_mul(m1, m2)
            if not s:
                continue
            sg = fs * s * (-1 if (popcount(T1) * mono_par(m2)) % 2 else 1)
            key = (m, T1 | T2)
            v = Z.get(key, 0) + sg * c1 * c2
            if v:
                Z[key] = v
            else:
                Z.pop(key, None)
    return Z


def superfield(a, al=(0, 0, 0, 0, 0)):
    return {(((a, T, tuple(al)),), T): 1 for T in range(8)}


def R_d(X):
    """d = sum_u dx_u D_u :  d(m dx_T) = sum_u (-1)^{par m} (D_u m) dx_u dx_T."""
    Z = {}
    for (m, T), c in X.items():
        for u in FORMDIR:
            b = BIT_OF_U[u]
            fs = form_sign(b, T)
            if not fs:
                continue
            Dm = D({m: 1}, u)
            sg = fs * (-1 if mono_par(m) else 1)
            for mm, cc in Dm.items():
                key = (mm, T | b)
                v = Z.get(key, 0) + sg * c * cc
                if v:
                    Z[key] = v
                else:
                    Z.pop(key, None)
    return Z


def R_coeff(X, T):
    return {m: c for (m, TT), c in X.items() if TT == T}


class Theory:
    """f[a][b] = dict c -> structure constant f^c_{ab} (integers), [t_a, t_b] = sum_c f^c_ab t_c."""

    def __init__(self, dim, f, name):
        self.dim, self.f, self.name = dim, f, name
        self._Q0 = {}
        self._Q = {}
        self._build()

    def _build(self):
        for a in range(self.dim):
            A = superfield(a)
            QA = R_d(A)
            for b in range(self.dim):
                for c in range(b + 1, self.dim):
                    fbc = self.f[b][c].get(a, 0)
                    if fbc:
                        prod = R_mul(superfield(b), superfield(c))
                        for k, v in prod.items():
                            nv = QA.get(k, 0) + fbc * v
                            if nv:
                                QA[k] = nv
                            else:
                                QA.pop(k, None)
            for T in range(8):
                self._Q0[(a, T)] = R_coeff(QA, T)

    def Q_gen(self, g):
        if g not in self._Q:
            self._Q[g] = Dpow(self._Q0[(g[0], g[1])], g[2])
        return self._Q[g]

    def Q(self, P):
        return derivation(P, self.Q_gen, True)


def eta_gen(v_bit):
    """iota_{d/dx_v} A = sum_T (-1)^{par A_T} A_T iota_v(dx_T);  iota_v(dx_T) = (-1)^{#bits of T below v} dx_{T - v}.
    component of dx_{T'} (v not in T'): (-1)^{par A_{T'+v}} (-1)^{#bits of T' below v} A_{T'+v}."""
    def act(g):
        a, T, al = g
        if T & v_bit:
            return {}
        Tn = T | v_bit
        below = popcount(T & (v_bit - 1))
        sg = (-1) ** ((1 - popcount(Tn)) % 2 + below)
        return {((a, Tn, al),): sg}
    return act


def eta(P, v_bit):
    return derivation(P, eta_gen(v_bit), True)


# ---------------------------------------------------------------- Lie algebras
def abelian(n):
    return Theory(n, [[{} for _ in range(n)] for _ in range(n)], "abelian%d" % n)


def sl2():
    # basis e=0, f=1, h=2 ; [e,f]=h, [h,e]=2e, [h,f]=-2f
    f = [[{} for _ in range(3)] for _ in range(3)]
    def setb(a, b, d):
        f[a][b] = dict(d)
        f[b][a] = {k: -v for k, v in d.items()}
    setb(0, 1, {2: 1})
    setb(2, 0, {0: 2})
    setb(2, 1, {1: -2})
    return Theory(3, f, "sl2")


def gl2():
    # basis e=0, f=1, h=2, u=3 (identity matrix, central)
    f = [[{} for _ in range(4)] for _ in range(4)]
    def setb(a, b, d):
        f[a][b] = dict(d)
        f[b][a] = {k: -v for k, v in d.items()}
    setb(0, 1, {2: 1})
    setb(2, 0, {0: 2})
    setb(2, 1, {1: -2})
    th = Theory(4, f, "gl2")
    th.charge_fn = lambda g: {0: 2, 1: -2, 2: 0, 3: 0}[g[0]]
    return th


def sl2_charge(g):
    """ad(h)-charge of the dual coordinate (cochain charge = minus element charge; sign irrelevant)."""
    return {0: 2, 1: -2, 2: 0}[g[0]]


# ---------------------------------------------------------------- weights and enumeration
def gen_cost(g):
    """(c_tau, c_zb1, c_zb2, hz1, hz2): c_v = al_v + [v in T]."""
    a, T, al = g
    return (al[0] + (T & 1), al[3] + (T >> 1 & 1), al[4] + (T >> 2 & 1), al[1], al[2])


def mono_weight(m):
    """functional weight w = (c_tau - 1, c_zb1 - 1, c_zb2 - 1) and holomorphic bidegree (a, b)."""
    c = [0] * 5
    for g in m:
        gc = gen_cost(g)
        for i in range(5):
            c[i] += gc[i]
    return (c[0] - 1, c[1] - 1, c[2] - 1), (c[3], c[4])


def mono_ghost(m):
    return sum(ghost(g) for g in m)


def gens_within(dim, budget, holo_only=False):
    """all generators with cost <= budget componentwise. budget = (b_tau, b_zb1, b_zb2, a, b)."""
    out = []
    Ts = [0] if holo_only else range(8)
    for a in range(dim):
        for T in Ts:
            tb = (T & 1, T >> 1 & 1, T >> 2 & 1)
            if any(tb[i] > budget[i] for i in range(3)):
                continue
            rng = [range(budget[0] - tb[0] + 1), range(budget[3] + 1), range(budget[4] + 1),
                   range(budget[1] - tb[1] + 1), range(budget[2] - tb[2] + 1)]
            if holo_only:
                rng[0] = rng[3] = rng[4] = range(1)
            for a0 in rng[0]:
                for a1 in rng[1]:
                    for a2 in rng[2]:
                        for a3 in rng[3]:
                            for a4 in rng[4]:
                                out.append((a, T, (a0, a1, a2, a3, a4)))
    return sorted(out)


def enumerate_monos(dim, budget, holo_only=False, max_arity=None, charge=None, charge_fn=None):
    """all monomials (multisets of generators, odd ones distinct) with total cost == budget exactly.
    returns dict ghost -> list of monomials."""
    gens = gens_within(dim, budget, holo_only)
    costs = [gen_cost(g) for g in gens]
    pars = [gpar(g) for g in gens]
    res = {}
    target = tuple(budget)

    def rec(start, rem, cur):
        if all(r == 0 for r in rem):
            m = tuple(cur)
            if m and (charge is None or sum(charge_fn(g) for g in m) == charge):
                res.setdefault(mono_ghost(m), []).append(m)
            # zero-cost generators may still be appended; continue
        if max_arity is not None and len(cur) >= max_arity:
            return
        for i in range(start, len(gens)):
            c = costs[i]
            if any(c[j] > rem[j] for j in range(5)):
                continue
            nrem = tuple(rem[j] - c[j] for j in range(5))
            cur.append(gens[i])
            # odd generators at most once: next start i+1; even: i
            rec(i + 1 if pars[i] else i, nrem, cur)
            cur.pop()

    rec(0, target, [])
    for k in res:
        res[k] = sorted(set(res[k]))
    return res
