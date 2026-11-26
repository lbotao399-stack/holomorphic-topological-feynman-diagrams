"""Bare Schwinger numerator  P(t) = U^(n+H) D(t)  (T = 0) modulo p.

P is homogeneous of degree n^2 - 1, degree <= n - 1 in every t_e [src cores/fastker.py, cores/numer.py].
Samples at t = (1, 1 + e), on a downward-closed set with |e| <= n^2 - 1,
individual degrees <= n-1 and, by default, the certified subset-rank bounds
sum_{j in S} e_j <= n rank(C_S)-1.  Set rank_bounds=False for the old grid.
(vectorized: currents by M ell = -C t y0, P = det[C_re Y^i_e]_{(i,r), e != 0} / (t_0 U^(n(H-1)))),
Newton divided differences and conversion to monomials directly on L (sparse lower-set
interpolation; L is downward closed, so every shifted index stays in L).
"""
import os
import random
import sys
from math import comb, factorial

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..', 'code')))
from htk import brute            # noqa: E402
from htk.graph import Graph      # noqa: E402
from htk.treeflow import tree_flow as particular_tree_flow  # noqa: E402
from htk.support import numerator_support  # noqa: E402

from modlin import det_mod, solve_mod, modinv_vec   # noqa: E402


def tree_flow(G, lam_full):
    """integer particular flow y0 (N x H) supported on the Kruskal spanning tree: B^T y0 = lambda"""
    return np.array(particular_tree_flow(G, lam_full), dtype=np.int64)


class LowerSet:
    """points e in [0,deg]^nvar with |e| <= total, ordered by cube index sum_k e_k (deg+1)^k;
    optional constraints (variable indices, upper bound) further restrict support.
    Generated digit by digit, without allocating a full cube.  Every restriction
    is downward closed, so the Newton transform is unchanged.
    """

    def __init__(self, nvar, deg, total, *, constraints=(), max_points=10000000):
        if any(not isinstance(x, int) or x < 0 for x in (nvar, deg, total)):
            raise ValueError('nonnegative integer dimensions and degrees are required')
        if max_points is not None and (not isinstance(max_points, int) or max_points < 1):
            raise ValueError('max_points must be positive or None')
        self.nvar, self.deg, self.total = nvar, deg, total
        b = deg + 1
        self.K = b ** nvar
        if self.K > np.iinfo(np.int64).max:
            raise ValueError('cube indices exceed int64; use series-fibre integration')
        parsed = []
        for indices, bound in constraints:
            indices = tuple(sorted(set(indices)))
            if (not isinstance(bound, int) or any(not isinstance(j, int) or not 0 <= j < nvar
                                                 for j in indices)):
                raise ValueError('invalid lower-set constraint')
            parsed.append((indices, bound))
        self.constraints = tuple(parsed)
        impossible = any(bound < 0 for _, bound in self.constraints)
        # digits from the most significant (variable nvar-1) down, so rows come out sorted by cube index
        digit_dtype = np.int8 if deg <= 127 else np.int64
        pts = np.zeros((0 if impossible else 1, 0), dtype=digit_dtype)
        sums = np.zeros(len(pts), dtype=np.int64)
        for k in range(nvar):
            axis = nvar - 1 - k
            maximum = np.minimum(deg, total - sums)
            for indices, bound in self.constraints:
                if axis not in indices:
                    continue
                columns = [nvar - 1 - j for j in indices if j > axis]
                partial = pts[:, columns].sum(axis=1, dtype=np.int64) if columns else 0
                maximum = np.minimum(maximum, bound - partial)
            counts = np.maximum(0, maximum + 1)
            length = int(counts.sum())
            if max_points is not None and length > max_points:
                raise MemoryError('lower set exceeds max_points=%d; use series-fibre integration '
                                  'or explicitly raise the allocation budget' % max_points)
            parent = np.repeat(np.arange(len(pts), dtype=np.int64), counts)
            digit = np.arange(length, dtype=np.int64) - np.repeat(np.cumsum(counts) - counts, counts)
            new = np.empty((length, k + 1), dtype=digit_dtype)
            new[:, :k] = pts[parent]
            new[:, k] = digit
            pts = new
            sums = sums[parent] + digit
        if any(not indices and bound < 0 for indices, bound in self.constraints):
            pts = pts[:0]
        # pts columns are variables nvar-1 ... 0 (most significant first); reverse to e_0 .. e_{nvar-1}
        self.dig = np.ascontiguousarray(pts[:, ::-1])
        self.M = len(self.dig)
        cube = np.zeros(self.M, dtype=np.int64)
        for k in range(nvar - 1, -1, -1):
            cube = cube * b + self.dig[:, k].astype(np.int64)
        self.cube = cube
        assert self.M < 2 or (np.diff(cube) > 0).all()
        # A dense lookup may otherwise allocate more memory than the sparse grid.
        if self.K <= 16000000 and self.K <= max(1024, 8 * self.M):
            pos = np.full(self.K, -1, dtype=np.int32)
            pos[cube] = np.arange(self.M, dtype=np.int32)
            self._pos = pos
        else:
            self._pos = None
        self.stride = [b ** k for k in range(nvar)]

    def pos_of(self, ci):
        if self._pos is not None:
            return self._pos[ci]
        return np.searchsorted(self.cube, ci)

    @property
    def pos(self):
        return self

    def __getitem__(self, ci):
        return self.pos_of(ci)

    def newton_to_monomial(self, vals, p):
        """vals: samples f(1 + e) mod p on the lower set -> monomial coefficients in t (same order)"""
        deg, b = self.deg, self.deg + 1
        f = vals.astype(np.int64) % p
        invf = [pow(factorial(k), p - 2, p) for k in range(b)]
        def sel_axis(a):
            col = self.dig[:, a]
            return [np.nonzero(col == d)[0] for d in range(b)]
        # forward differences along every axis
        for a in range(self.nvar):
            s = self.stride[a]
            sa = sel_axis(a)
            old = f.copy()
            for d in range(1, b):
                I = sa[d]
                if not len(I):
                    continue
                ci = self.cube[I]
                acc = old[I].copy()
                for j in range(d):
                    c = ((-1) ** (d - j) * comb(d, j)) % p
                    acc = (acc + c * old[self.pos[ci - (d - j) * s]]) % p
                f[I] = acc * invf[d] % p
            del old, sa
        # Newton basis prod_{j=1}^{k} (t - j) -> monomials
        falling = [[1]]
        for k in range(1, b):
            prev = falling[-1]
            cur = [0] * (k + 1)
            for i, c in enumerate(prev):
                cur[i + 1] += c
                cur[i] -= k * c
            falling.append(cur)
        for a in range(self.nvar):
            s = self.stride[a]
            sa = sel_axis(a)
            new = np.zeros_like(f)
            for k in range(b):
                I = sa[k]
                if not len(I):
                    continue
                ci = self.cube[I]
                src = f[I]
                for pw in range(k + 1):
                    c = falling[k][pw] % p
                    if not c:
                        continue
                    tgt = I if pw == k else self.pos[ci - (k - pw) * s]
                    new[tgt] = (new[tgt] + c * src) % p
            f = new
        return f


class BareNumerator:
    """P = U^(n+H) D of an htk graph with integer sources, modulo primes."""

    def __init__(self, G, lam_full, *, rank_bounds=True, max_flats=4096, max_points=10000000):
        self.G = G
        self.n, self.N = G.loops, G.ne
        self.H = len(lam_full[0])
        if self.n * self.H != self.N - 1:
            raise ValueError('not a T = 0 holomorphic graph')
        self.lam_full = [[int(x) for x in v] for v in lam_full]
        self.C = np.array(G.cycle_matrix(), dtype=np.int64)
        self.y0 = tree_flow(G, self.lam_full)
        self.deg = self.n - 1
        self.total = self.n * self.n - 1
        self.support = numerator_support(G, max_flats=max_flats) if rank_bounds else None
        self.max_points = max_points
        self._ls = None

    @property
    def lowerset(self):
        if self._ls is None:
            constraints = []
            if self.support is not None:
                constraints = [(tuple(e - 1 for e in edges), bound)
                               for edges, bound in self.support.constraints]
                if self.support.zero:
                    constraints.append((tuple(range(self.N - 1)), -1))
            self._ls = LowerSet(self.N - 1, self.deg, self.total,
                                constraints=constraints, max_points=self.max_points)
        return self._ls

    def exponents(self):
        """full homogeneous exponent vectors: (total - |e|, e), without int8 overflow."""
        ls = self.lowerset
        e0 = self.total - ls.dig.astype(np.int64).sum(axis=1)
        return np.concatenate([e0[:, None], ls.dig], axis=1)

    def samples_mod(self, tpts, p, batch=20000):
        """P(t) mod p at the rows of tpts (B x N int64 positive)"""
        n, N, H = self.n, self.N, self.H
        C = self.C
        out = np.zeros(len(tpts), dtype=np.int64)
        for s in range(0, len(tpts), batch):
            t = tpts[s:s + batch] % p
            B = len(t)
            Mm = np.einsum('re,be,se->brs', C, t, C) % p
            rhs = -np.einsum('re,be,ei->bri', C, t, self.y0) % p
            x, U = solve_mod(Mm, rhs, p)
            if np.any(U == 0) or np.any(t[:, 0] == 0):
                raise ArithmeticError('singular modular sample: U or the gauge time is zero; '
                                      'choose another prime or sampling point')
            Uell = x * U[:, None, None] % p
            Y = (U[:, None, None] * self.y0[None, :, :] + np.einsum('re,bri->bei', C, Uell)) % p
            Psi = np.einsum('re,bei->bire', C, Y).reshape(B, H * n, N)[:, :, 1:] % p
            d = det_mod(Psi.copy(), p)
            den = t[:, 0] % p
            Uk = np.ones(B, dtype=np.int64)
            for _ in range(n * (H - 1)):
                Uk = Uk * U % p
            out[s:s + B] = d * modinv_vec(den * Uk % p, p) % p
        return out

    def samples_fast(self, tpts, p):
        if np.any(np.asarray(tpts)[:, 0] % p == 0):
            raise ArithmeticError('singular modular sample: the gauge time is zero')
        try:
            import fastsample
        except ImportError:
            return self.samples_mod(tpts, p)
        out = fastsample.samples(np.ascontiguousarray(tpts, dtype=np.int64), self.C, self.y0, self.n, self.H, p)
        if (out < 0).any():
            raise ArithmeticError('U = 0 mod p at a sample point')
        return out

    def coefficients_mod(self, p, check=True, chunk=400000):
        """monomial coefficients (aligned with exponents()) mod p"""
        ls = self.lowerset
        vals = np.empty(ls.M, dtype=np.int64)
        for s in range(0, ls.M, chunk):
            d = ls.dig[s:s + chunk].astype(np.int64)
            tpts = np.concatenate([np.ones((len(d), 1), dtype=np.int64), 1 + d], axis=1)
            vals[s:s + chunk] = self.samples_fast(tpts, p)
        coef = ls.newton_to_monomial(vals, p)
        if check:
            rng = random.Random(p)
            t = np.array([[rng.randint(2, p - 1) for _ in range(self.N)]], dtype=np.int64)
            direct = int(self.samples_mod(t, p)[0])
            if int(self.samples_fast(t, p)[0]) != direct:
                raise ArithmeticError('numba sampler differs from numpy sampler')
            ex = self.exponents()
            val = np.ones(ls.M, dtype=np.int64)
            for e in range(self.N):
                pw = [1]
                for _ in range(self.total):
                    pw.append(pw[-1] * int(t[0, e]) % p)
                val = val * np.array(pw, dtype=np.int64)[ex[:, e].astype(np.int64)] % p
            interp = int((val * coef % p).sum() % p)
            if interp != direct:
                raise ArithmeticError('interpolation check failed mod %d' % p)
        return coef
