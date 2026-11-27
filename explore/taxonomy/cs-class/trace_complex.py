"""Translation-coinvariant CE complex of sl_N[[z_1..z_H]] in the trace basis (formal N >= number of inputs).

Cochain (k inputs x_i = X_i e^{p_i.z}, sum p = 0):  omega = sum_pi c_pi(p) T_pi(X),
T_pi = prod over cycles (a -> pi(a) -> ...) of tr(X_a X_pi(a) ...), pi fixed-point free (tr X = 0).
Alternation: c_{s pi s^-1}(p) = sgn(s) c_pi(s^* p), (s^* p)_a = p_{s(a)}.  A cochain is stored by its canonical
polynomials c_{pi_lambda} (pi_lambda = consecutive cycles, lambda decreasing), polynomial = dict
{(I, J): coefficient} meaning sum coef [p_I][p_J], [p_I] = det(p_{I_1},..,p_{I_H}), I sorted.
Differential: (d omega)(x_0..x_k) = sum_{i<j} (-1)^{i+j} omega([x_i,x_j], x_0..^i..^j..x_k).
"""
from itertools import combinations, permutations, product
from fractions import Fraction as Fr
from functools import lru_cache
import numpy as np


# ------------------------------------------------------------------ permutations
def sgn(perm):
    perm = list(perm)
    seen, s = [False] * len(perm), 1
    for i in range(len(perm)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def cycles(succ):
    seen, out = set(), []
    for s in range(len(succ)):
        if s in seen:
            continue
        c, x = [], s
        while x not in seen:
            seen.add(x); c.append(x); x = succ[x]
        out.append(c)
    return out


def cycle_type(succ):
    return tuple(sorted((len(c) for c in cycles(succ)), reverse=True))


@lru_cache(None)
def canonical(lam):
    succ, start = [], 0
    for l in lam:
        for i in range(l):
            succ.append(start + (i + 1) % l)
        start += l
    return tuple(succ)


@lru_cache(None)
def centralizer(lam):
    """all z with z pi_lam z^{-1} = pi_lam (as tuples z[a])."""
    blocks, start = [], 0
    for l in lam:
        blocks.append(list(range(start, start + l)))
        start += l
    k = start
    # group blocks by length
    bylen = {}
    for b in blocks:
        bylen.setdefault(len(b), []).append(b)
    groups = list(bylen.values())
    out = []

    def rec(gi, z):
        if gi == len(groups):
            out.append(tuple(z))
            return
        bl = groups[gi]
        l = len(bl[0])
        for perm in permutations(range(len(bl))):
            for rots in product(range(l), repeat=len(bl)):
                z2 = list(z)
                for src, dst_i, r in zip(bl, perm, rots):
                    dst = bl[dst_i]
                    for t in range(l):
                        z2[src[t]] = dst[(t + r) % l]
                rec(gi + 1, z2)
    rec(0, [None] * k)
    return out


def conjugator(succ):
    """sigma with sigma pi_lam sigma^{-1} = succ (sigma[canonical letter] = letter of succ)."""
    cs = sorted(cycles(succ), key=lambda c: -len(c))
    lam = tuple(len(c) for c in cs)
    sigma = []
    for c in cs:
        sigma.extend(c)
    return lam, tuple(sigma)


def compose(a, b):
    """(a o b)(x) = a(b(x))."""
    return tuple(a[b[x]] for x in range(len(b)))


def inverse(a):
    out = [None] * len(a)
    for i, x in enumerate(a):
        out[x] = i
    return tuple(out)


def merge(succ, a, b):
    """merge adjacent letters a -> b into letter 0; other letters in increasing order -> 1..k-1."""
    k = len(succ)
    rest = [x for x in range(k) if x not in (a, b)]
    newid = {x: i + 1 for i, x in enumerate(rest)}
    newid[a] = 0
    newid[b] = 0
    out = [None] * (k - 1)
    for x in rest:
        y = succ[x]
        out[newid[x]] = newid[y]
    out[0] = newid[succ[b]]
    return tuple(out), [a, b], rest


# ------------------------------------------------------------------ polynomials in minors
def poly_pullback(P, sigma, H):
    """Q(p) = P(sigma^* p), (sigma^* p)_a = p_{sigma(a)}."""
    out = {}
    for (I, J), c in P.items():
        I2 = [sigma[i] for i in I]
        J2 = [sigma[j] for j in J]
        sI = _sort_sign(I2)
        sJ = _sort_sign(J2)
        if sI is None or sJ is None:
            continue
        key = _key(tuple(sorted(I2)), tuple(sorted(J2)))
        out[key] = out.get(key, 0) + c * sI[0] * sJ[0]
    return {k: v for k, v in out.items() if v}


def _sort_sign(seq):
    if len(set(seq)) < len(seq):
        return None
    inv = sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] > seq[j])
    return (-1 if inv % 2 else 1,)


def _key(I, J):
    return (I, J) if I <= J else (J, I)


def poly_add(P, Q, c=1):
    out = dict(P)
    for k, v in Q.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v}


# ------------------------------------------------------------------ alternating cochains
class Cochain:
    """canonical data {lam: poly}; k = number of inputs; H = holomorphic dimension."""

    def __init__(self, k, H, data=None):
        self.k, self.H = k, H
        self.data = data or {}

    def add_term(self, succ, P, coef=1):
        """add alt(coef * P * T_succ)."""
        lam, sig = conjugator(succ)       # sig pi_lam sig^-1 = succ
        s1 = inverse(sig)                  # s1 succ s1^-1 = pi_lam
        acc = self.data.get(lam, {})
        for z in centralizer(lam):
            s = compose(z, s1)
            Q = poly_pullback(P, s, self.H)
            acc = poly_add(acc, Q, coef * sgn(s))
        if acc:
            self.data[lam] = acc
        elif lam in self.data:
            del self.data[lam]

    def scaled(self, c):
        return Cochain(self.k, self.H, {l: {k: c * v for k, v in P.items()} for l, P in self.data.items()})

    def plus(self, other, c=1):
        d = dict(self.data)
        for l, P in other.data.items():
            d[l] = poly_add(d.get(l, {}), P, c)
            if not d[l]:
                del d[l]
        return Cochain(self.k, self.H, d)


def fpf_types(k):
    """partitions of k with all parts >= 2, decreasing."""
    out = []

    def rec(n, mx, cur):
        if n == 0:
            out.append(tuple(cur)); return
        for l in range(min(n, mx), 1, -1):
            rec(n - l, l, cur + [l])
    rec(k, k, [])
    return out


# ------------------------------------------------------------------ evaluation
def sample_points(k, H, M, rng, rmax=2):
    pts = rng.integers(-rmax, rmax + 1, size=(M, k - 1, H))
    last = -pts.sum(axis=1, keepdims=True)
    return np.concatenate([pts, last], axis=1).astype(np.int64)


def minors_table(pts, H):
    """dict sorted H-tuple -> vector of minors over points (int64, exact)."""
    M, k, _ = pts.shape
    tab = {(): np.ones(M, dtype=np.int64)}
    for I in combinations(range(k), H):
        sub = pts[:, list(I), :]                    # M x H x H, rows = vectors
        tab[I] = _det_int(sub)
    return tab


def _det_int(sub):
    """exact integer determinant of M x n x n arrays (n <= 4) by Leibniz."""
    M, n, _ = sub.shape
    out = np.zeros(M, dtype=np.int64)
    for perm in permutations(range(n)):
        s = sgn(perm)
        term = np.ones(M, dtype=np.int64)
        for r in range(n):
            term = term * sub[:, r, perm[r]]
        out += s * term
    return out


def det_int_general(sub):
    """exact integer determinants for n up to 5 (Leibniz, 120 terms)."""
    return _det_int(sub)


def eval_poly(P, tab, p, coef_mod=None):
    """sum coef [I][J] at all points, mod p (coefficients Fractions/ints)."""
    first = next(iter(tab.values()))
    out = np.zeros(first.shape[0], dtype=np.int64)
    for (I, J), c in P.items():
        cm = _mod(c, p)
        if cm == 0:
            continue
        v = (tab[I] % p) * (tab[J] % p) % p
        out = (out + cm * v) % p
    return out


def _mod(c, p):
    if isinstance(c, Fr):
        return (c.numerator % p) * pow(c.denominator % p, p - 2, p) % p
    return int(c) % p


def coordinates(om, pts_by_type, tabs_by_type, p):
    """vector of c_{pi_mu}(q) for each target type mu (dict mu -> points)."""
    parts = []
    for mu in sorted(pts_by_type):
        P = om.data.get(mu, {})
        parts.append(eval_poly(P, tabs_by_type[mu], p) if P else np.zeros(pts_by_type[mu].shape[0], dtype=np.int64))
    return np.concatenate(parts)


class DiffPlan:
    """precomputed merged/pulled-back sample points for d : C^k -> C^{k+1} at coordinates (mu, point)."""

    def __init__(self, k, H, M, rng, types=None, rmax=2):
        self.k, self.H = k, H
        self.types = sorted(types or fpf_types(k + 1))
        self.pts = {mu: sample_points(k + 1, H, M, rng, rmax) for mu in self.types}
        self.tabs = {mu: minors_table(self.pts[mu], H) for mu in self.types}
        self.items = []          # per mu: list of (lam, sign, tab of pulled back merged points)
        for mu in self.types:
            pi = canonical(mu)
            lst = []
            for a in range(k + 1):
                b = pi[a]
                if b == a or pi[b] == a:
                    continue
                s_ab = (1 if a < b else -1) * (-1) ** (a + b)
                rho, ab, rest = merge(pi, a, b)
                lam, sig0 = conjugator(rho)         # sig0 pi_lam sig0^-1 = rho
                if min(lam) < 2:
                    continue
                q = self.pts[mu]
                merged = np.concatenate([(q[:, a, :] + q[:, b, :])[:, None, :], q[:, rest, :]], axis=1)
                pulled = merged[:, list(sig0), :]   # (sig0^* q')_c = q'_{sig0(c)}
                lst.append((lam, s_ab * sgn(sig0), minors_table(pulled, H)))
            self.items.append((mu, lst))

    def apply(self, om, p):
        """coordinates of d(om) (om a k-cochain)."""
        parts = []
        for mu, lst in self.items:
            n = self.pts[mu].shape[0]
            v = np.zeros(n, dtype=np.int64)
            for lam, s, tab in lst:
                P = om.data.get(lam)
                if P:
                    v = (v + s * eval_poly(P, tab, p)) % p
            parts.append(v % p)
        return np.concatenate(parts)


# ------------------------------------------------------------------ bracket of symmetric-minor cochains
def _words(succ):
    return [tuple(c) for c in cycles(succ)]


def half_bracket_sym(Tdict, H, c, N, sl=True):
    """1/2 {j,j} for j = c [p_1..p_H] T(X_0..X_H), T a symmetric invariant given in the trace basis
    {succ on H+1 letters: coef}.  Bracket (flagship): {a,b} = alt(gamma_0)/((m-1)!(n-1)!),
    gamma_0 = a(x_1..x_{m-1}, X_e) kappa^{ee'} b(X_e', y_1..y_{n-1}); here
    gamma_0 = c^2 [p_x1..p_xH][p_y1..p_yH] sum_e T(X_1..X_H, T_e) T(T^e, Y_1..Y_H).
    Letters of the result: x's 0..H-1, y's H..2H-1.  sl: Casimir of sl_N (kappa = tr)."""
    from math import factorial
    k = 2 * H
    om = Cochain(k, H)
    pref = Fr(c * c, 2 * factorial(H) ** 2)
    acc = {}
    E1, E2 = ('e', 1), ('e', 2)
    for s1, c1 in Tdict.items():
        w1 = [tuple(E1 if a == H else ('x', a) for a in w) for w in _words(s1)]
        for s2, c2 in Tdict.items():
            w2 = [tuple(E2 if a == 0 else ('y', a - 1) for a in w) for w in _words(s2)]
            i1 = next(i for i, w in enumerate(w1) if E1 in w)
            i2 = next(i for i, w in enumerate(w2) if E2 in w)
            u = w1[i1]; ku = u.index(E1); U = u[ku + 1:] + u[:ku]
            v = w2[i2]; kv = v.index(E2); V = v[kv + 1:] + v[:kv]
            rest = [w for i, w in enumerate(w1) if i != i1] + [w for i, w in enumerate(w2) if i != i2]
            terms = [(rest + [U + V], Fr(1))]
            if sl:
                terms.append((rest + [U, V], Fr(-1, N)))
            for words, f in terms:
                if any(len(w) == 0 for w in words):
                    # tr(1) = N
                    nE = sum(1 for w in words if len(w) == 0)
                    f = f * N ** nE
                    words = [w for w in words if w]
                if sl and any(len(w) == 1 for w in words):
                    continue
                succ = [None] * k
                for w in words:
                    for i, a in enumerate(w):
                        b = w[(i + 1) % len(w)]
                        ia = a[1] if a[0] == 'x' else H + a[1]
                        ib = b[1] if b[0] == 'x' else H + b[1]
                        succ[ia] = ib
                succ = tuple(succ)
                acc[succ] = acc.get(succ, 0) + pref * c1 * c2 * f
    P = {(tuple(range(H)), tuple(range(H, 2 * H))): 1}
    for succ, cf in acc.items():
        if cf:
            om.add_term(succ, P, cf)
    return om


def sym_trace_dict(lengths, sl=True):
    """polarization of prod_l tr(X^l) (lengths, sum = n) as {succ on n letters: coef} (symmetric)."""
    from math import factorial
    n = sum(lengths)
    acc = {}
    base, start = [], 0
    for l in lengths:
        base.append(list(range(start, start + l)))
        start += l
    for perm in permutations(range(n)):
        succ = [None] * n
        for b in base:
            for i in range(len(b)):
                succ[perm[b[i]]] = perm[b[(i + 1) % len(b)]]
        succ = tuple(succ)
        acc[succ] = acc.get(succ, 0) + Fr(1, factorial(n))
    return {k: v for k, v in acc.items() if v}
