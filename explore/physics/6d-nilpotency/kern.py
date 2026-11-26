"""Kernels of the rooted C^{3|2} graphs, evaluated modulo a prime and vectorized over momentum assignments.

I_Gamma = (-1)^{n(n+1)/2} K_Gamma (htk CONVENTIONS section 3, S2 eq-kernel-def); K_Gamma = int g(xi) d xi in the
spred ratio chart.  Spokes are oriented (v, root); B^T y = lambda, the root source is -sum_v lambda_v, and the momentum of spoke j is
    p_j = y_{e_j} = sum_v c_{jv} lambda_v,   c_{jv} = - X[v][e_j],   sum_j c_{jv} = 1   (htk displacement X),
and an input jet d_a on the letter bound to spoke j multiplies the integrand by (p_j)_a (vops.py convention).

XP: polynomial in the integration coordinates with coefficient vectors (one entry per assignment), mod P.
Exponents are packed BITS bits per coordinate into one int64 key.
"""
import os
import sys
import random
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(PROOT, 'code'))
sys.path.insert(0, HERE)

import sympy as sp                                    # noqa: E402
from htk import Graph, SPReduction, brute             # noqa: E402
from rooted_trees import graphs                       # noqa: E402
from wick_weights import wick_sign, automorphisms     # noqa: E402

P = 2147483647
BITS = 5
MASK = (1 << BITS) - 1
NFACT = 600
FACT = INVFACT = None


def set_prime(p):
    global P, FACT, INVFACT
    P = p
    f = [1] * NFACT
    for i in range(1, NFACT):
        f[i] = f[i - 1] * i % P
    FACT = np.array(f, dtype=np.int64)
    INVFACT = np.array([pow(x, P - 2, P) for x in f], dtype=np.int64)


set_prime(P)


def modq(x):
    """Fraction / int / sympy QQ -> residue."""
    if hasattr(x, 'numerator') and hasattr(x, 'denominator'):
        num, den = int(x.numerator), int(x.denominator)
    else:
        num, den = int(x), 1
    if den % P == 0:
        raise ZeroDivisionError('denominator divisible by P')
    return num % P * pow(den, P - 2, P) % P


def ratrec(a, p=None):
    """Rational reconstruction of a residue (|num|, den < sqrt(p/2))."""
    p = p or P
    a %= p
    r0, r1, s0, s1 = p, a, 0, 1
    bound = int((p / 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


# =============================================================================== XP polynomials
class XP:
    __slots__ = ('k', 'C')

    def __init__(self, k, C):
        self.k, self.C = k, C          # k: int64 (M,), sorted unique; C: int64 (M, L)

    @property
    def L(self):
        return self.C.shape[1]

    def __len__(self):
        return len(self.k)


def _combine(k, C):
    if len(k) == 0:
        return XP(k, C)
    order = np.argsort(k, kind='stable')
    k = k[order]
    C = C[order]
    starts = np.flatnonzero(np.concatenate(([True], k[1:] != k[:-1])))
    if len(starts) < len(k):
        C = np.add.reduceat(C, starts, axis=0)
        k = k[starts]
    C %= P
    keep = np.any(C != 0, axis=1)
    if not keep.all():
        k, C = k[keep], C[keep]
    return XP(k, C)


def xp_const(v, L):
    v = np.broadcast_to(np.asarray(v, dtype=np.int64) % P, (L,)).copy()
    return _combine(np.zeros(1, dtype=np.int64), v[None, :])


def xp_zero(L):
    return XP(np.zeros(0, dtype=np.int64), np.zeros((0, L), dtype=np.int64))


def xp_add(*xs):
    if not xs:
        raise ValueError('xp_add of nothing')
    x0 = xs[0]
    xs = [x for x in xs if len(x)]
    if not xs:
        return x0
    return _combine(np.concatenate([x.k for x in xs]), np.concatenate([x.C for x in xs]))


def xp_scale(A, v):
    """v: residue scalar or (L,) array."""
    v = np.asarray(v, dtype=np.int64) % P
    return _combine(A.k.copy(), A.C * v % P)


def maxdeg(k):
    if len(k) == 0:
        return 0
    m = 0
    kk = k.copy()
    while np.any(kk):
        m = max(m, int((kk & MASK).max()))
        kk >>= BITS
    return m


def xp_mul(A, B, budget=2_000_000):
    if len(A) == 0 or len(B) == 0:
        return xp_zero(A.L)
    if maxdeg(A.k) + maxdeg(B.k) > MASK:
        raise OverflowError('exponent overflow in XP product')
    if len(A) < len(B):
        A, B = B, A
    L = A.L
    nb = len(B)
    chunk = max(1, budget // max(1, nb * L))
    run = None
    for i0 in range(0, len(A), chunk):
        a_k = A.k[i0:i0 + chunk]
        a_C = A.C[i0:i0 + chunk]
        kk = (a_k[:, None] + B.k[None, :]).ravel()
        CC = (a_C[:, None, :] * B.C[None, :, :]) % P
        r = _combine(kk, CC.reshape(-1, L))
        run = r if run is None else _combine(np.concatenate([run.k, r.k]), np.concatenate([run.C, r.C]))
    return run


def decode(k, nc):
    return np.array([(k >> (BITS * j)) & MASK for j in range(nc)], dtype=np.int64).T.reshape(len(k), nc)


def matmul_mod(W, B):
    """(W @ B) mod P for residues W (n x m), B (m x L), m < 2^16."""
    lo = B & 0xFFFF
    hi = B >> 16
    r = (W @ lo) % P
    r = (r + ((W @ hi) % P << 16)) % P
    return r


def full_integral(A, B, blocks, nc):
    """int A * B over all blocks (B may be None), without forming the product."""
    ea = decode(A.k, nc)
    allv = sorted(j for b in blocks for j in b)
    if B is None:
        w = np.ones(len(A), dtype=np.int64)
        for j in allv:
            w = w * FACT[ea[:, j]] % P
        for b in blocks:
            w = w * INVFACT[len(b) + ea[:, b].sum(axis=1)] % P
        return (A.C * w[:, None] % P).sum(axis=0) % P
    eb = decode(B.k, nc)
    out = np.zeros(A.L, dtype=np.int64)
    budget = 2_000_000
    ca = max(1, budget // max(1, len(B)))
    for i0 in range(0, len(A), ca):
        ea_ = ea[i0:i0 + ca]
        step = min(1 << 15, max(1, budget // len(ea_)))
        for j0 in range(0, len(B), step):
            j1 = j0 + step
            eb_ = eb[j0:j1]
            W = np.ones((len(ea_), len(eb_)), dtype=np.int64)
            for j in allv:
                W = W * FACT[ea_[:, j][:, None] + eb_[:, j][None, :]] % P
            for b in blocks:
                tot = len(b) + ea_[:, b].sum(axis=1)[:, None] + eb_[:, b].sum(axis=1)[None, :]
                W = W * INVFACT[tot] % P
            WB = matmul_mod(W, B.C[j0:j1])
            out = (out + (A.C[i0:i0 + ca] * WB % P).sum(axis=0)) % P
    return out


def xp_integrate(A, block):
    """int over the simplex {b_j >= 0, sum_j b_j <= 1} of the coordinates in block (list of indices):
    prod_j b_j^{a_j} -> prod a_j! / (m + sum a_j)!."""
    if len(A) == 0:
        return A
    m = len(block)
    w = np.ones(len(A), dtype=np.int64)
    tot = np.full(len(A), m, dtype=np.int64)
    k = A.k.copy()
    for j in block:
        e = (A.k >> (BITS * j)) & MASK
        w = w * FACT[e] % P
        tot += e
        k &= ~np.int64(MASK << (BITS * j))
    w = w * INVFACT[tot] % P
    return _combine(k, A.C * w[:, None] % P)


def xp_value(A):
    """constant term as (L,) array (all coordinates must be integrated)."""
    if len(A) == 0:
        return np.zeros(A.L, dtype=np.int64)
    if np.any(A.k != 0):
        raise ValueError('non-constant XP')
    return A.C[0] % P


def support(A):
    s = set()
    if len(A) == 0:
        return s
    kk = np.bitwise_or.reduce(A.k)
    j = 0
    while kk:
        if kk & MASK:
            s.add(j)
        kk >>= BITS
        j += 1
    return s


def integrate_factors(factors, blocks, L):
    """int prod(factors) over the product of simplices; greedy order as htk spred, with the last two
    non-constant factors integrated bilinearly over all remaining blocks."""
    if any(len(f) == 0 for f in factors):
        return np.zeros(L, dtype=np.int64)
    nc = 1 + max([j for b in blocks for j in b] + [0])
    active = [(f, support(f)) for f in factors]
    blocks = [list(b) for b in blocks]
    const = np.ones(L, dtype=np.int64)
    while True:
        nonconst = [(f, s) for f, s in active if s]
        for f, s in active:
            if not s:
                const = const * xp_value(f) % P
        active = nonconst
        if len(active) <= 2:
            if not active:
                w = 1
                for b in blocks:
                    w = w * int(INVFACT[len(b)]) % P
                return const * w % P
            A = active[0][0]
            B = active[1][0] if len(active) == 2 else None
            return const * full_integral(A, B, blocks, nc) % P
        if not blocks:
            raise ValueError('unintegrated coordinates')

        def cost(block):
            t = 1
            for f, s in active:
                if s.intersection(block):
                    t *= len(f)
            return t, -len(block)
        block = min(blocks, key=cost)
        blocks.remove(block)
        dep = [(f, s) for f, s in active if s.intersection(block)]
        active = [(f, s) for f, s in active if not s.intersection(block)]
        if not dep:
            const = const * int(INVFACT[len(block)]) % P
            continue
        dep.sort(key=lambda fs: len(fs[0]))
        comb = dep[0][0]
        sup = set(dep[0][1])
        for f, s in dep[1:]:
            comb = xp_mul(comb, f)
            sup |= s
        r = xp_integrate(comb, block)
        if len(r) == 0:
            return np.zeros(L, dtype=np.int64)
        active.append((r, sup - set(block)))


# =============================================================================== symbolic -> XP
def pack(exps):
    key = 0
    for j, e in enumerate(exps):
        if e:
            if e > MASK:
                raise OverflowError('exponent too large')
            key |= int(e) << (BITS * j)
    return key


class SymPoly:
    """Polynomial in (coords, lambda components) stored as xi-keys, lambda exponents and residues."""

    def __init__(self, poly, nc, lam_index):
        """poly: sympy PolyElement in ring (coords + parameter symbols); lam_index: parameter name -> index."""
        names = [str(g) for g in poly.ring.gens]
        pidx = [lam_index[nm] for nm in names[nc:]]
        nl = len(lam_index)
        keys, lexp, coef = [], [], []
        for mon, c in poly.terms():
            keys.append(pack(mon[:nc]))
            e = [0] * nl
            for i, x in zip(pidx, mon[nc:]):
                e[i] = x
            lexp.append(e)
            coef.append(modq(c))
        self.keys = np.array(keys, dtype=np.int64)
        self.lexp = np.array(lexp, dtype=np.int64).reshape(len(keys), nl)
        self.coef = np.array(coef, dtype=np.int64)
        self.nl = nl

    def evaluate(self, lam, L):
        """lam: (nl, L) residues -> XP in the coordinates."""
        T = len(self.keys)
        if T == 0:
            return xp_zero(L)
        vals = np.broadcast_to(self.coef[:, None], (T, L)).copy()
        maxe = int(self.lexp.max()) if self.lexp.size else 0
        for i in range(self.nl):
            col = self.lexp[:, i]
            if not col.any():
                continue
            pw = [np.ones(L, dtype=np.int64)]
            for _ in range(maxe):
                pw.append(pw[-1] * lam[i] % P)
            pw = np.array(pw)
            vals = vals * pw[col] % P
        return _combine(self.keys.copy(), vals)


# =============================================================================== graph kernels
class GraphKernel:
    def __init__(self, n, number, edges):
        self.n, self.number, self.edges = n, number, [tuple(e) for e in edges]
        self.name = 'G%d.%d' % (n, number)
        root = 2 * n + 1
        self.root = root
        self.verts = list(range(root))
        self.spokes = [i for i, e in enumerate(self.edges) if root in e]
        assert all(self.edges[i][1] == root for i in self.spokes)
        deg = {v: sum(v in e for e in self.edges) for v in self.verts}
        self.outputs = [v for v in self.verts if deg[v] == 2]
        legs = {}
        for v in self.verts:
            ins = [i for i, e in enumerate(self.edges) if e[1] == v]
            outs = [i for i, e in enumerate(self.edges) if e[0] == v]
            L = [('e', i) for i in ins] + ([('out', v)] if deg[v] == 2 else []) + [('e', i) for i in outs]
            assert len(L) == 3
            legs[v] = L
        self.legs = legs
        self.kappa = Fr(wick_sign(self.edges, root), automorphisms(self.edges, root))
        self.sign_n = (-1) ** (n * (n + 1) // 2)
        self.G = Graph(root + 1, self.edges)
        self.red = SPReduction(self.G)
        self.coords = list(self.red.coords)
        self.nc = len(self.coords)
        if self.nc * BITS > 63:
            raise OverflowError('too many coordinates for the key packing')
        cidx = {c: i for i, c in enumerate(self.coords)}
        self.blocks = [[cidx[b] for b in blk.indep] for blk in self.red.blocks]
        # lambda parameters: output v, component a
        self.lam_names = ['l%d_%d' % (v, a) for v in self.outputs for a in range(3)]
        self.lam_index = {nm: i for i, nm in enumerate(self.lam_names)}
        self.lam_sym = {v: [sp.Symbol('l%d_%d' % (v, a)) for a in range(3)] if v in self.outputs else [0, 0, 0]
                        for v in self.verts}
        self._dens = None
        self._moments = None
        self._symk = {}

    # ---------------------------------------------------------------- symbolic density
    def density_factors(self):
        if self._dens is None:
            red = self.red
            if red.tight_plan:
                dens = red.density(self.lam_sym, check=1)
                ring = dens.ring
                self._ring = ring
                self._dens = [SymPoly(f, self.nc, self.lam_index) for f in dens.factors]
                self._dens_sym = list(dens.factors)
            else:
                from blocktri import BlockTriangular
                bt = BlockTriangular(self)
                self._ring = bt.ring
                self._dens = bt
                self._dens_sym = None
        return self._dens

    def ring(self):
        self.density_factors()
        return self._ring

    def moment_polys(self):
        """c[j][v] (sympy PolyElements in the coordinate ring): p_j = sum_v c[j][v] lambda_v."""
        if self._moments is None:
            X = self.red.displacement()
            self._moments = [{v: -X[v][e] for v in self.outputs} for e in self.spokes]
            # translation invariance: sum_j c[j][v] = 1
            for v in self.outputs:
                s = sum((self._moments[j][v] for j in range(len(self.spokes))), self.red.ring.zero)
                assert s == self.red.ring.one, 'moment maps are not convex'
            # sources at non-output vertices vanish
        return self._moments

    def moment_xp(self, lam, L):
        """p_j components as XP for constant lam: dict v -> (3, L) residues.  Returns [j][a] -> XP."""
        mom = self.moment_polys()
        out = []
        for j in range(len(self.spokes)):
            comps = []
            for a in range(3):
                keys, vals = [], []
                for v in self.outputs:
                    for mon, c in mom[j][v].terms():
                        keys.append(pack(mon[:self.nc]))
                        vals.append(modq(c) * lam[v][a] % P)
                if keys:
                    comps.append(_combine(np.array(keys, dtype=np.int64), np.array(vals, dtype=np.int64)))
                else:
                    comps.append(xp_zero(L))
            out.append(comps)
        return out

    def density_xp(self, lam, L):
        """list of XP factors of the density for constant lam: dict v -> (3, L) residues."""
        dens = self.density_factors()
        arr = np.zeros((len(self.lam_names), L), dtype=np.int64)
        for v in self.outputs:
            for a in range(3):
                arr[self.lam_index['l%d_%d' % (v, a)]] = lam[v][a]
        if isinstance(dens, list):
            return [f.evaluate(arr, L) for f in dens]
        return dens.evaluate(arr, L)

    def kernel_numeric(self, lam, L, extra=()):
        """I_Gamma with insertions: int g * prod(extra) for constant output momenta lam (v -> (3, L))."""
        facs = self.density_xp(lam, L) + list(extra)
        val = integrate_factors(facs, self.blocks, L)
        return val * (self.sign_n % P) % P

    # ---------------------------------------------------------------- symbolic integrated kernel
    def kernel_symbolic(self, jets):
        """jets: tuple over spokes of tuples of components.  Returns (lexp (T, nl), coef (T,)) of
        I_Gamma(lambda) with the input moments prod_j prod_{a in jets[j]} (p_j)_a."""
        jets = tuple(tuple(x) for x in jets)
        if jets in self._symk:
            return self._symk[jets]
        dens = self.red.density(self.lam_sym, check=0) if self.red.tight_plan else None
        if dens is None:
            raise NotImplementedError('symbolic kernel for a non-tight plan')
        ring = dens.ring
        X = self.red.displacement()
        names = [str(g) for g in ring.gens]
        gens = dict(zip(names, ring.gens))
        extra = []
        for j, comps in enumerate(jets):
            e = self.spokes[j]
            for a in comps:
                pj = ring.zero
                for v in self.outputs:
                    pj = pj - X[v][e].set_ring(ring) * gens['l%d_%d' % (v, a)]
                extra.append(pj)
        res = dens._integrate_polynomial(list(dens.factors) + extra)
        sp_ = SymPoly(res, self.nc, self.lam_index)
        assert not sp_.keys.any()
        out = (sp_.lexp, sp_.coef * (self.sign_n % P) % P)
        self._symk[jets] = out
        return out


def all_graphs(n):
    return [GraphKernel(n, k, e) for k, e in enumerate(graphs(n), 1)]


# =============================================================================== evaluation of a symbolic kernel
def eval_symbolic(lexp, coef, args, L):
    """sum_t coef_t prod_i args[i]^lexp[t,i]; args[i] is an (L,) residue array or an XP.  Returns XP."""
    T, nl = lexp.shape
    const = [i for i in range(nl) if not isinstance(args[i], XP)]
    xpv = [i for i in range(nl) if isinstance(args[i], XP)]
    vals = np.broadcast_to(coef[:, None] % P, (T, L)).copy()
    maxe = int(lexp.max()) if lexp.size else 0
    for i in const:
        col = lexp[:, i]
        if not col.any():
            continue
        pw = [np.ones(L, dtype=np.int64)]
        for _ in range(maxe):
            pw.append(pw[-1] * args[i] % P)
        vals = vals * np.array(pw)[col] % P
    if not xpv:
        return xp_const(vals.sum(axis=0) % P, L) if T else xp_zero(L)
    pat = lexp[:, xpv]
    keys = {}
    for t in range(T):
        keys.setdefault(tuple(pat[t]), []).append(t)
    powc = {}

    def power(i, e):
        if (i, e) not in powc:
            powc[(i, e)] = xp_const(1, L) if e == 0 else xp_mul(power(i, e - 1), args[i])
        return powc[(i, e)]
    parts = []
    for pt, ts in keys.items():
        c = vals[ts].sum(axis=0) % P
        if not c.any():
            continue
        term = None
        for i, e in zip(xpv, pt):
            if e:
                term = power(i, e) if term is None else xp_mul(term, power(i, e))
        if term is None:
            term = xp_const(1, L)
        parts.append(xp_scale(term, c))
    return xp_add(*parts) if parts else xp_zero(L)


# =============================================================================== joint (outer x inner) evaluation
def xp_shift(A, shift):
    """re-index the coordinates of A by +shift (A must not use coordinates >= 63 // BITS - shift)."""
    return XP(A.k << np.int64(BITS * shift), A.C)


def _as_xp(x, L):
    return x if isinstance(x, XP) else xp_const(x, L)


def sympoly_mixed(sp_, args, L, shift):
    """SymPoly in (own coordinates, lambda) evaluated with lambda components args[i] that are (L,) arrays
    or XP in outer coordinates; own coordinates are shifted by `shift`.  Returns XP."""
    T = len(sp_.keys)
    if T == 0:
        return xp_zero(L)
    const = [i for i in range(sp_.nl) if not isinstance(args[i], XP)]
    xpv = [i for i in range(sp_.nl) if isinstance(args[i], XP)]
    vals = np.broadcast_to(sp_.coef[:, None], (T, L)).copy()
    maxe = int(sp_.lexp.max()) if sp_.lexp.size else 0
    for i in const:
        col = sp_.lexp[:, i]
        if not col.any():
            continue
        pw = [np.ones(L, dtype=np.int64)]
        for _ in range(maxe):
            pw.append(pw[-1] * args[i] % P)
        vals = vals * np.array(pw)[col] % P
    keys = sp_.keys << np.int64(BITS * shift)
    groups = {}
    pat = sp_.lexp[:, xpv] if xpv else np.zeros((T, 0), dtype=np.int64)
    for t in range(T):
        groups.setdefault(tuple(pat[t]), []).append(t)
    powc = {}

    def power(i, e):
        if (i, e) not in powc:
            powc[(i, e)] = xp_const(1, L) if e == 0 else xp_mul(power(i, e - 1), args[i])
        return powc[(i, e)]
    parts = []
    for pt, ts in groups.items():
        ts = np.array(ts)
        inner = _combine(keys[ts].copy(), vals[ts])
        if len(inner) == 0:
            continue
        term = inner
        for i, e in zip(xpv, pt):
            if e:
                term = xp_mul(term, power(i, e))
        parts.append(term)
    return xp_add(*parts) if parts else xp_zero(L)


def density_mixed(gk, lam, L, shift):
    """density factors of gk (tight plan) with output momenta lam[v][a] = array or XP (outer coordinates)."""
    dens = gk.density_factors()
    if not isinstance(dens, list):
        raise NotImplementedError('mixed evaluation of a block-triangular density')
    args = []
    for nm in gk.lam_names:
        v, a = map(int, nm[1:].split('_'))
        args.append(lam[v][a])
    return [sympoly_mixed(f, args, L, shift) for f in dens]


def moments_mixed(gk, lam, L, shift):
    """[j][a] -> XP of p_j components with lam[v][a] = array or XP (outer coordinates)."""
    mom = gk.moment_polys()
    out = []
    for j in range(len(gk.spokes)):
        comps = []
        for a in range(3):
            acc = None
            for v in gk.outputs:
                keys, vals = [], []
                for mon, c in mom[j][v].terms():
                    keys.append(pack(mon[:gk.nc]))
                    vals.append(modq(c))
                if not keys:
                    continue
                cv = _combine(np.array(keys, dtype=np.int64) << np.int64(BITS * shift),
                              np.broadcast_to(np.array(vals, dtype=np.int64)[:, None], (len(keys), L)).copy())
                term = xp_mul(cv, _as_xp(lam[v][a], L))
                acc = term if acc is None else xp_add(acc, term)
            comps.append(acc if acc is not None else xp_zero(L))
        out.append(comps)
    return out
