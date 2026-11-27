"""Cartan detector.  For commuting inputs h_v e^{p_v.z} (h_v in a Cartan subalgebra, sum p = 0) every wedge is a
CE cycle of g[[z]] (all brackets vanish), so an exact cochain vanishes on it.  Values (mod p):

  o'|_C = sum_G (s c_1/|Aut|)_G sigma_G sum_{s in S_8} sgn(s) D_theta,G(p_s) C_std,G(h_s)
        = sum_G (...) a! b! c! sum_{ordered set partitions (S_A,S_B,S_C)} sgn(s_0) D_theta(p_{s_0}) W_G(h; S_A,S_B,S_C),
  W_G(h; S) = sum_{x,y,z} f_xyz f^xyz prod_{v in S_A} x(h_v) prod_{v in S_B} y(h_v) prod_{v in S_C} z(h_v)
  (weight basis; R_P acts diagonally; D_theta alternates and W is symmetric within each path).

  1/2{j_T,j_T}|_C = c^2/(2 (4!)^2) sum_s sgn(s) [p_s0..p_s3][p_s4..p_s7] T'(h_s),
  T'(h) = sum_{mu,nu} prod_{a<4} mu(h_a) prod_{a>=4} nu(h_a) (mu,nu)_kappa  for T = Sym tr_V X^5 (weights mu of V).
"""
import numpy as np
from itertools import combinations, permutations
from fractions import Fraction as Fr
from common import parity

SHAPES = [(1, 5, 5), (2, 4, 5), (3, 3, 5), (3, 4, 4)]
KC = {(1, 5, 5): Fr(1, 5184), (2, 4, 5): Fr(-1, 9072), (3, 3, 5): Fr(1, 10080), (3, 4, 4): Fr(-1, 22680)}
AUT = {(1, 5, 5): 4, (2, 4, 5): 2, (3, 3, 5): 4, (3, 4, 4): 4}
SIGMA = -1
DTH = {(1, 5, 5): [((0, 1, 2, 3), (4, 5, 6, 7), 1)], (2, 4, 5): [((0, 1, 2, 3), (4, 5, 6, 7), 1)],
       (3, 3, 5): [((0, 1, 2, 3), (4, 5, 6, 7), 1)],
       (3, 4, 4): [((0, 1, 2, 3), (4, 5, 6, 7), 1), ((0, 1, 2, 4), (3, 5, 6, 7), -1), ((0, 1, 3, 4), (2, 5, 6, 7), 1)]}


def fmod(c, p):
    c = Fr(c)
    return (c.numerator % p) * pow(c.denominator % p, p - 2, p) % p


def det4_mod(vs, p):
    from common import det_int
    return det_int([list(map(int, v)) for v in vs]) % p


def structure_triples(L, p):
    """arrays (x, y, z, w) with w = f_xyz f^xyz mod p, f^xyz = kappa^{xa} kappa^{yb} kappa^{zc} f_abc."""
    import scipy.sparse as sps
    d = L.dim
    a, b, c, v = L.f_low
    G = L.G.tocsc()
    s = L.s
    # f^{xyz} s^3 = sum G[x,a] G[y,b] G[z,c] f_abc ; G pairs weight w with -w
    up = {}
    Gcols = {k: (G.indices[G.indptr[k]:G.indptr[k + 1]], G.data[G.indptr[k]:G.indptr[k + 1]]) for k in range(d)}
    for aa, bb, cc, vv in zip(a, b, c, v):
        xa, ga = Gcols[aa]
        yb, gb = Gcols[bb]
        zc, gc = Gcols[cc]
        for x, g1 in zip(xa, ga):
            for y, g2 in zip(yb, gb):
                for z, g3 in zip(zc, gc):
                    key = (int(x), int(y), int(z))
                    up[key] = up.get(key, 0) + int(g1) * int(g2) * int(g3) * int(vv)
    low = {(int(x), int(y), int(z)): int(w) for x, y, z, w in zip(a, b, c, v)}
    inv_s3 = pow(s ** 3 % p, p - 2, p)
    X, Y, Z, W = [], [], [], []
    for key, lv in low.items():
        uv = up.get(key, 0)
        if uv:
            X.append(key[0]); Y.append(key[1]); Z.append(key[2])
            W.append((lv % p) * (uv % p) % p * inv_s3 % p)
    return np.array(X), np.array(Y), np.array(Z), np.array(W, dtype=np.int64)


def set_partitions(shape):
    a, b, c = (n - 1 for n in shape)
    out = []
    for SA in combinations(range(8), a):
        rest = [x for x in range(8) if x not in SA]
        for SB in combinations(rest, b):
            SC = tuple(x for x in rest if x not in SB)
            seq = list(SA) + list(SB) + list(SC)
            out.append((SA, SB, SC, parity(seq), seq))
    return out


class Detector:
    def __init__(self, L, p):
        self.L, self.p = L, p
        self.wt = L.wt.astype(np.int64)            # d x r
        self.trip = structure_triples(L, p)
        self.parts = {sh: set_partitions(sh) for sh in SHAPES}

    def sample(self, rng, cmax=5, pmax=4):
        r = self.wt.shape[1]
        hc = rng.integers(-cmax, cmax + 1, size=(8, r))
        mom = rng.integers(-pmax, pmax + 1, size=(7, 4))
        mom = np.vstack([mom, -mom.sum(axis=0)])
        return hc, mom

    def o_prime(self, hc, mom):
        p = self.p
        alpha = (self.wt @ hc.T) % p              # d x 8 : alpha[x, v] = x(h_v)
        # products over subsets
        prod = {(): np.ones(self.wt.shape[0], dtype=np.int64)}
        for k in range(1, 9):
            for S in combinations(range(8), k):
                prod[S] = prod[S[:-1]] * alpha[:, S[-1]] % p
        X, Y, Z, W = self.trip
        total = 0
        for sh in SHAPES:
            a, b, c = (n - 1 for n in sh)
            from math import factorial
            coef = fmod(KC[sh] / AUT[sh] * SIGMA * factorial(a) * factorial(b) * factorial(c), p)
            acc = 0
            for SA, SB, SC, sg, seq in self.parts[sh]:
                q = [mom[seq[s]] for s in range(8)]
                D = 0
                for I, J, e in DTH[sh]:
                    D += e * det4_mod([q[i] for i in I], p) * det4_mod([q[j] for j in J], p)
                D %= p
                if D == 0:
                    continue
                w = int(((W * prod[SA][X] % p) * prod[SB][Y] % p * prod[SC][Z] % p).sum() % p)
                acc = (acc + sg * D * w) % p
            total = (total + coef * acc) % p
        return total

    def half_bracket(self, weights_V, hc, mom, c=120):
        """1/2{j,j}|_C for j = c [p1..p4] Sym tr_V X^5; weights_V: list of integer vectors mu(h_i) (Cartan basis)."""
        p = self.p
        L = self.L
        r = self.wt.shape[1]
        # kappa restricted to the Cartan block, inverse
        K = L.kap.toarray()
        cart = L.cartan
        Kh = [[Fr(int(K[i, j])) for j in cart] for i in cart]
        import sympy as sp
        Khi = sp.Matrix(Kh).inv()
        mus = np.array(weights_V, dtype=np.int64)          # nV x r
        # inner products (mu, nu) = mu^T Kh^{-1} nu
        ip = np.zeros((len(mus), len(mus)), dtype=np.int64)
        for i in range(len(mus)):
            for j in range(len(mus)):
                val = sum(Fr(int(mus[i, a])) * Fr(Khi[a, b]) * Fr(int(mus[j, b])) for a in range(r) for b in range(r))
                ip[i, j] = fmod(val, p)
        mh = (mus @ hc.T) % p                              # nV x 8
        total = 0
        for S1 in combinations(range(8), 4):
            S2 = tuple(x for x in range(8) if x not in S1)
            sg = parity(list(S1) + list(S2))
            D = det4_mod([mom[i] for i in S1], p) * det4_mod([mom[i] for i in S2], p) % p
            u = np.ones(len(mus), dtype=np.int64)
            for v in S1:
                u = u * mh[:, v] % p
            w = np.ones(len(mus), dtype=np.int64)
            for v in S2:
                w = w * mh[:, v] % p
            Tp = int((u @ (ip % p) % p) @ w % p) if False else int(((u[:, None] * ip % p) * w[None, :] % p).sum() % p)
            total = (total + sg * D * Tp) % p
        from math import factorial
        pref = fmod(Fr(c * c, 2 * factorial(4) ** 2) * factorial(4) ** 2, p)
        # Sym tr_V X^5 polarization: T(h_0..h_4) = sum_mu prod mu(h_i)
        return total * pref % p


def sym_tensor_from_poly(Pfun, r, p):
    """symmetric 5-tensor T on C^r (mod p) with T(x,x,x,x,x) = Pfun(x) (Pfun: integer polynomial function),
    by polarization: T(e_i0..e_i4) = (1/5!) sum_{S subset} (-1)^{5-|S|} P(sum_{s in S} e_{i_s})."""
    from itertools import product as iprod
    from math import factorial
    inv120 = pow(120, p - 2, p)
    T = {}
    for idx in iprod(range(r), repeat=5):
        key = tuple(sorted(idx))
        if key in T:
            continue
        tot = 0
        for k in range(1, 6):
            for S in combinations(range(5), k):
                x = [0] * r
                for s in S:
                    x[idx[s]] += 1
                tot += (-1) ** (5 - k) * Pfun(x)
        T[key] = (tot % p) * inv120 % p
    return T


def half_bracket_tensor(L, T, hc, mom, p, c=120):
    """1/2{j,j}|_C for j = c [p1..p4] T, T symmetric invariant restricted to h (dict sorted index -> value mod p)."""
    from math import factorial
    import sympy as sp
    r = len(L.cartan)
    K = L.kap.toarray()
    cart = L.cartan
    Khi = sp.Matrix([[int(K[i, j]) for j in cart] for i in cart]).inv()
    Kinv = [[fmod(Fr(int(sp.Rational(Khi[i, j]).p), int(sp.Rational(Khi[i, j]).q)), p) for j in range(r)] for i in range(r)]
    # full tensor array
    import itertools
    Tarr = np.zeros((r,) * 5, dtype=np.int64)
    for idx in itertools.product(range(r), repeat=5):
        Tarr[idx] = T[tuple(sorted(idx))]
    Kinv = np.array(Kinv, dtype=np.int64)

    hcs = np.asarray(hc, dtype=np.int64)          # small signed entries: tensordot stays below 2^63

    def contract4(S):
        A = Tarr
        for v in S:
            A = np.tensordot(hcs[v], A, axes=([0], [0])) % p
        return [int(x) for x in A]   # vector over the last index (python ints)
    total = 0
    for S1 in combinations(range(8), 4):
        S2 = tuple(x for x in range(8) if x not in S1)
        sg = parity(list(S1) + list(S2))
        D = det4_mod([mom[i] for i in S1], p) * det4_mod([mom[i] for i in S2], p) % p
        u = contract4(S1)
        w = contract4(S2)
        Tp = sum(u[i] * int(Kinv[i, j]) * w[j] for i in range(r) for j in range(r)) % p
        total = (total + sg * D * Tp) % p
    pref = fmod(Fr(c * c, 2), p)
    return total * pref % p
