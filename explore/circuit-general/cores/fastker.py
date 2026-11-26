"""Vectorized modular kernel evaluation for large T=0 graphs (e.g. D=3 four loops, N=13).

U^(n+H) D = det[C_re Y^i_e]_{(i,r), e != 0} / (t_0 U^(n(H-1))),  Y = U y,
evaluated mod p on the grid {t = (1, 1+e): e in [0,n-1]^(N-1), |e| <= n^2-1}, Newton-
interpolated mod p on the dense cube, converted to monomials, pushed forward along series
paths to the core (p unchanged: P_orig = P_core * prod Gamma(a_e)/Gamma(a_path)) and contracted
with the exact core tables (cores.period) reduced mod p.  CRT over primes and rational
reconstruction of the final coefficients (labels 1, z3, z5).  Bare kernels only.
"""
import os
import sys
import time
from fractions import Fraction as F
from math import comb, factorial

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../../code'))
from htk import brute                      # noqa: E402
from htk.graph import Graph                # noqa: E402

import cores                               # noqa: E402
import modsolve                            # noqa: E402

PRIMES = [2147483629, 2147483587, 2147483579, 2147483563, 2147483549, 2147483543, 2147483497, 2147483489]


def modinv_vec(a, p):
    """a^(p-2) mod p elementwise (int64, p < 2^31)"""
    r = np.ones_like(a)
    b = a % p
    e = p - 2
    while e:
        if e & 1:
            r = r * b % p
        b = b * b % p
        e >>= 1
    return r


def det_mod(A, p):
    """batched determinant mod p of A (B, k, k) int64; destroys A"""
    B, k, _ = A.shape
    det = np.ones(B, dtype=np.int64)
    idx = np.arange(B)
    for c in range(k):
        col = A[:, c:, c]
        nz = col != 0
        has = nz.any(axis=1)
        piv = c + np.argmax(nz, axis=1)
        det[~has] = 0
        swap = piv != c
        if swap.any():
            rows_c = A[idx, c, :].copy()
            A[idx, c, :] = A[idx, piv, :]
            A[idx, piv, :] = rows_c
            det[swap] = (-det[swap]) % p
        d = A[:, c, c] % p
        det = det * d % p
        inv = modinv_vec(np.where(d == 0, 1, d), p)
        if c + 1 < k:
            f = A[:, c + 1:, c] * inv[:, None] % p
            A[:, c + 1:, :] = (A[:, c + 1:, :] - f[:, :, None] * A[:, c, None, :]) % p
    return det % p


def solve_mod(Mm, R, p):
    """batched solve Mm x = R mod p (Mm: (B,n,n), R: (B,n,h)); returns x and det"""
    B, n, _ = Mm.shape
    A = np.concatenate([Mm, R], axis=2) % p
    det = np.ones(B, dtype=np.int64)
    idx = np.arange(B)
    for c in range(n):
        nz = A[:, c:, c] != 0
        piv = c + np.argmax(nz, axis=1)
        swap = piv != c
        if swap.any():
            rows_c = A[idx, c, :].copy()
            A[idx, c, :] = A[idx, piv, :]
            A[idx, piv, :] = rows_c
            det[swap] = (-det[swap]) % p
        d = A[:, c, c] % p
        det = det * d % p
        inv = modinv_vec(d, p)
        A[:, c, :] = A[:, c, :] * inv[:, None] % p
        for r in range(n):
            if r == c:
                continue
            f = A[:, r, c].copy()
            A[:, r, :] = (A[:, r, :] - f[:, None] * A[:, c, :]) % p
    return A[:, :, n:], det


class FastKernel:
    def __init__(self, nv, edges, lam):
        self.G = Graph(nv, list(edges))
        self.edges = list(edges)
        G = self.G
        self.n, self.N = G.loops, G.ne
        lam_nr = brute.sources(G, lam)
        self.H = len(lam_nr[0])
        if self.n * self.H != self.N - 1:
            raise ValueError('not a T=0 holomorphic graph')
        for v in lam_nr:
            for x in v:
                if F(x).denominator != 1:
                    raise ValueError('integer sources required')
        self.lam_nr = [[int(x) for x in v] for v in lam_nr]
        self.C = np.array(G.cycle_matrix(), dtype=np.int64)          # n x N
        # particular flow y0 on the Kruskal spanning tree (B^T y0 = lambda)
        lam_full = self.lam_nr + [[-sum(v[i] for v in self.lam_nr) for i in range(self.H)]]
        import numer
        Y = numer.tree_current_polys(G, lam_full)
        # one tree: pick the cotree monomial that appears; the flow for that tree is y0
        cots = set()
        for e in range(self.N):
            for i in range(self.H):
                cots |= set(Y[e][i])
        cot0 = sorted(cots)[0]
        self.y0 = np.array([[Y[e][i].get(cot0, 0) for i in range(self.H)] for e in range(self.N)], dtype=np.int64)
        # series paths of the core
        self.paths = self.series_paths()

    def series_paths(self):
        """list of edge groups (paths through bivalent vertices) and the core edge list"""
        from collections import defaultdict
        edges = self.edges
        deg = defaultdict(int)
        for u, v in edges:
            deg[u] += 1
            deg[v] += 1
        groups = [[i] for i in range(len(edges))]
        ends = [list(e) for e in edges]
        changed = True
        while changed:
            changed = False
            for v in list(deg):
                if deg[v] != 2:
                    continue
                gi = [g for g in range(len(groups)) if v in ends[g]]
                if len(gi) != 2:
                    continue
                a, b = gi
                ea = ends[a][0] if ends[a][1] == v else ends[a][1]
                eb = ends[b][0] if ends[b][1] == v else ends[b][1]
                groups[a] = groups[a] + groups[b]
                ends[a] = [ea, eb]
                del groups[b]
                del ends[b]
                deg[v] = 0
                changed = True
                break
        return groups, [tuple(e) for e in ends]

    def samples_mod(self, grid, p, batch=40000):
        """U^(n+H) D mod p at t = (1, 1+e) for the rows of grid"""
        n, N, H = self.n, self.N, self.H
        C = self.C
        out = np.zeros(len(grid), dtype=np.int64)
        for s in range(0, len(grid), batch):
            e = grid[s:s + batch].astype(np.int64)
            B = len(e)
            t = np.concatenate([np.ones((B, 1), dtype=np.int64), 1 + e], axis=1)     # B x N
            Mm = np.einsum('re,be,se->brs', C, t, C) % p                               # B x n x n
            rhs = -np.einsum('re,be,ei->bri', C, t, self.y0) % p                       # B x n x H
            x, U = solve_mod(Mm, rhs, p)                                                # ell = x
            # Y = U y = U y0 + C^T (U ell)
            Uell = x * U[:, None, None] % p
            Y = (U[:, None, None] * self.y0[None, :, :] + np.einsum('re,bri->bei', C, Uell)) % p   # B x N x H
            # Psi rows (i, r), columns e != 0
            Psi = np.einsum('re,bei->bire', C, Y).reshape(B, H * n, N)[:, :, 1:] % p
            d = det_mod(Psi.copy(), p)
            # divide by t_0 U^(n(H-1)), t_0 = 1
            Uk = np.ones(B, dtype=np.int64)
            for _ in range(n * (H - 1)):
                Uk = Uk * U % p
            out[s:s + B] = d * modinv_vec(Uk, p) % p
        return out


def grid_points(nvar, deg, total):
    """all e in [0,deg]^nvar with sum <= total, as int8 array, plus flat cube indices"""
    pts = np.array(np.meshgrid(*[np.arange(deg + 1)] * nvar, indexing='ij'), dtype=np.int8).reshape(nvar, -1).T
    keep = pts.sum(axis=1) <= total
    return pts[keep], np.nonzero(keep)[0]


def interpolate_mod(vals_cube, nvar, deg, p, valid=None):
    """Newton divided differences on t = 1..deg+1 per axis (dense cube), restriction of the
    Newton coefficients to the downward-closed set, then monomial basis"""
    A = vals_cube.reshape([deg + 1] * nvar).astype(np.int64)
    invf = [pow(factorial(k), p - 2, p) for k in range(deg + 1)]
    for ax in range(nvar):
        A = np.moveaxis(A, ax, 0)
        S = [A[j].copy() for j in range(deg + 1)]
        newA = np.empty_like(A)
        for k in range(deg + 1):
            acc = np.zeros_like(S[0])
            for j in range(k + 1):
                acc = (acc + ((-1) ** (k - j) * comb(k, j)) % p * S[j]) % p
            newA[k] = acc * invf[k] % p
        A = np.moveaxis(newA, 0, ax)
    if valid is not None:
        A = np.where(valid.reshape(A.shape), A, 0)
    # falling factorial (t-1)(t-2)...(t-k) -> monomials in t
    falling = [[1]]
    for k in range(1, deg + 1):
        prev = falling[-1]
        cur = [0] * (k + 1)
        for i, c in enumerate(prev):
            cur[i + 1] += c
            cur[i] -= k * c
        falling.append(cur)
    for ax in range(nvar):
        A = np.moveaxis(A, ax, 0)
        newA = np.zeros_like(A)
        for k in range(deg + 1):
            for pw, c in enumerate(falling[k]):
                if c:
                    newA[pw] = (newA[pw] + (c % p) * A[k]) % p
        A = np.moveaxis(newA, 0, ax)
    return A.reshape(-1)


def kernel(nv, edges, lam, primes=PRIMES, verbose=True, min_primes=2):
    fk = FastKernel(nv, edges, lam)
    n, N, H = fk.n, fk.N, fk.H
    q = n + H
    deg, total = n - 1, n * n - 1
    t0 = time.time()
    grid, flat = grid_points(N - 1, deg, total)
    if verbose:
        print('  grid %d points (cube %d)' % (len(grid), (deg + 1) ** (N - 1)), flush=True)
    groups, core_edges = fk.paths
    # monomial exponents on the valid part of the cube (k_1..k_{N-1}); k_0 = total - sum
    cube_valid = grid.astype(np.int16)                    # same order as flat (row-major cube)
    ksum = cube_valid.sum(axis=1).astype(np.int16)
    alpha = np.concatenate([(total - ksum)[:, None], cube_valid], axis=1).astype(np.int16) + 1
    core_alpha = np.stack([alpha[:, g].sum(axis=1) for g in groups], axis=1).astype(np.int16)
    base = int(alpha.max()) * N + 2
    key = np.zeros(len(alpha), dtype=np.int64)
    for j in range(core_alpha.shape[1]):
        key = key * base + core_alpha[:, j].astype(np.int64)
    valid = np.zeros((deg + 1) ** (N - 1), dtype=bool)
    valid[flat] = True
    results = []
    rec = None
    for ip, p in enumerate(primes):
        tp = time.time()
        vals = np.zeros((deg + 1) ** (N - 1), dtype=np.int64)
        vals[flat] = fk.samples_mod(grid, p)
        coef = interpolate_mod(vals, N - 1, deg, p, valid)[flat]
        # series factor mod p
        fact = [1]
        for i in range(1, int(alpha.max()) * N + 2):
            fact.append(fact[-1] * i % p)
        factn = np.array(fact, dtype=np.int64)
        sf = np.ones(len(alpha), dtype=np.int64)
        for g in groups:
            if len(g) > 1:
                for e in g:
                    sf = sf * factn[alpha[:, e].astype(np.int64) - 1] % p
                ap = alpha[:, g].sum(axis=1).astype(np.int64)
                sf = sf * modinv_vec(factn[ap - 1], p) % p
        w = coef * sf % p
        nzm = w != 0
        uk, inv = np.unique(key[nzm], return_inverse=True)
        acc = np.zeros(len(uk), dtype=object)
        wv = w[nzm]
        # accumulate mod p
        sums = np.zeros(len(uk), dtype=np.int64)
        np.add.at(sums, inv, wv)          # may exceed int64? wv < 2^31, counts < 2^32 -> < 2^63 ok
        sums %= p
        # first index of each key for its core alpha
        first = np.zeros(len(uk), dtype=np.int64)
        idxs = np.nonzero(nzm)[0]
        first[inv] = idxs
        tot = {}
        for j in range(len(uk)):
            ca = [int(x) for x in core_alpha[first[j]]]
            val = cores.period(core_edges, ca, q)
            for lab, x in val.items():
                tot[lab] = (tot.get(lab, 0) + int(sums[j]) * (x.numerator % p) * pow(x.denominator % p, p - 2, p)) % p
        results.append((p, tot))
        if verbose:
            print('  prime %d: %d monomials, %d core configurations  %.1fs' % (ip + 1, int(nzm.sum()), len(uk), time.time() - tp), flush=True)
        # CRT + reconstruction
        labels = sorted({l for _, t in results for l in t})
        M = 1
        accv = {l: 0 for l in labels}
        for pp, t in results:
            for l in labels:
                a = accv[l]
                x = t.get(l, 0)
                accv[l] = a + M * ((x - a) * pow(M, -1, pp) % pp)
            M *= pp
        new = {}
        okr = True
        for l in labels:
            r = modsolve.ratrec(accv[l], M)
            if r is None:
                okr = False
                break
            if r:
                new[l] = r
        if okr and rec == new and ip + 1 >= min_primes:
            if verbose:
                print('  stable after %d primes, total %.1fs' % (ip + 1, time.time() - t0), flush=True)
            return new
        rec = new if okr else None
    raise ArithmeticError('reconstruction did not stabilize')
