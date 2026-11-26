"""Dimension counts (generating functions) for the Budzik-Gaiotto meson complex of bg_complex_k and for
the U(N) gauge-theory cochains, per sector, without enumerating monomials.

Conventions as in bg_complex_k: generators T_w (graded-cyclic, nonvanishing), M_w (w may be empty),
P_w, MM_w, R_w (w nonempty) with flavours (a, b); parity |w|, |w|, 1+|w|, 1+|w|, |w|; ghost number
gh(w), gh(w), gh(w)-1, gh(w)-1, gh(w)-2; charges Lc, Rc.  For a word of weight omega and length l,
parity = l + Theta(omega) mod 2 and gh = l - Theta(omega).

chain_dims(alg, lam, k, charges, ghs) -> {(charge, g): dim V_{charge, g}} (float, exact below 2^53).
"""
from fractions import Fraction
from math import comb
import numpy as np


def _wle(a, b):
    return all(x <= y for x, y in zip(a, b))


def word_class_counts(alg, lam):
    """{(omega, length): number of nonempty words}, omega <= lam."""
    lam = tuple(lam)
    letters = [(alg.w[i]) for i in range(alg.n) if _wle(alg.w[i], lam)]
    cnt = {((0,) * len(lam), 0): 1}
    # DP over total weight in a graded order: process by length
    frontier = dict(cnt)
    out = {}
    while frontier:
        new = {}
        for (om, l), c in frontier.items():
            for wi in letters:
                no = tuple(x + y for x, y in zip(om, wi))
                if _wle(no, lam):
                    key = (no, l + 1)
                    new[key] = new.get(key, 0) + c
        for key, c in new.items():
            out[key] = out.get(key, 0) + c
        frontier = new
    return out


def _phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def cyclic_class_counts(alg, lam):
    """{(omega, length): number of nonvanishing graded-cyclic words}
    = coefficients of -sum_d phi(d)/d log(1 - fB(x^d) - (-1)^(d+1) fF(x^d))."""
    lam = tuple(lam)
    dim = len(lam)
    fB, fF = {}, {}
    for i in range(alg.n):
        wi = tuple(alg.w[i])
        if _wle(wi, lam):
            (fF if alg.par[i] else fB)[(wi, 1)] = (fF if alg.par[i] else fB).get((wi, 1), 0) + 1
    out = {}
    maxlen = sum(lam)
    for d in range(1, maxlen + 1):
        g = {}
        for (wi, l), c in fB.items():
            key = (tuple(d * x for x in wi), d * l)
            if _wle(key[0], lam):
                g[key] = g.get(key, 0) + c
        sF = 1 if (d + 1) % 2 == 0 else -1     # (-1)^(d+1)
        for (wi, l), c in fF.items():
            key = (tuple(d * x for x in wi), d * l)
            if _wle(key[0], lam):
                g[key] = g.get(key, 0) + sF * c
        g = {k: v for k, v in g.items() if v}
        if not g:
            continue
        # -log(1 - g) = sum_m g^m / m
        pw = {((0,) * dim, 0): 1}
        m = 0
        while True:
            m += 1
            npw = {}
            for (a, la), ca in pw.items():
                for (b, lb), cb in g.items():
                    o = tuple(x + y for x, y in zip(a, b))
                    if _wle(o, lam):
                        key = (o, la + lb)
                        npw[key] = npw.get(key, 0) + ca * cb
            pw = {k: v for k, v in npw.items() if v}
            if not pw:
                break
            for key, c in pw.items():
                out[key] = out.get(key, 0) + Fraction(_phi(d) * c, d * m)
    res = {}
    for key, v in out.items():
        if v:
            assert v.denominator == 1, (key, v)
            res[key] = int(v)
    return res


def _theta(om):
    return sum(om[:3])


class _Series:
    """Dense series: for each weight omega <= lam an array over (gh, B, l_1.., r_1..) offsets.
    k = 1: charges (B) only (Lc = Rc = B).  k = 2: (B, Lc1-Lc2, Rc1-Rc2)."""

    def __init__(self, lam, k, gh_rng, B_rng, l_rng):
        self.lam = tuple(lam)
        self.k = k
        self.g0, self.g1 = gh_rng
        self.B0, self.B1 = B_rng
        self.L = l_rng
        shape = [self.g1 - self.g0 + 1, self.B1 - self.B0 + 1]
        if k == 2:
            shape += [2 * l_rng + 1, 2 * l_rng + 1]
        self.shape = tuple(shape)
        self.weights = []
        rng = [range(x + 1) for x in self.lam]
        import itertools
        for om in itertools.product(*rng):
            self.weights.append(om)
        self.weights.sort(key=lambda o: (sum(o), o))
        self.A = {om: np.zeros(self.shape) for om in self.weights}
        z = (0,) * len(self.lam)
        idx = [0 - self.g0, 0 - self.B0] + ([self.L, self.L] if k == 2 else [])
        self.A[z][tuple(idx)] = 1.0

    def _shift(self, arr, dg, dB, dl=0, dr=0, times=1):
        out = np.zeros_like(arr)
        sh = [dg * times, dB * times] + ([dl * times, dr * times] if self.k == 2 else [])
        src = []
        dst = []
        for s, n in zip(sh, arr.shape):
            if abs(s) >= n:
                return None
            if s >= 0:
                src.append(slice(0, n - s))
                dst.append(slice(s, n))
            else:
                src.append(slice(-s, n))
                dst.append(slice(0, n + s))
        out[tuple(dst)] = arr[tuple(src)]
        return out

    def mult_class(self, om, odd, c, shifts):
        """Multiply by prod_{s in shifts} (1 - X_s)^(-c) (even) or (1 + X_s)^c (odd), X_s = x^om * shift s."""
        for (dg, dB, dl, dr) in shifts:
            self._mult_one(om, odd, c, dg, dB, dl, dr)

    def _mult_one(self, om, odd, c, dg, dB, dl, dr):
        lam = self.lam
        nmax = min((lam[i] // om[i]) for i in range(len(om)) if om[i] > 0)
        coefs = []
        for n in range(1, nmax + 1):
            coefs.append(comb(c, n) if odd else comb(n + c - 1, n))
        # process weights in decreasing order so that sources are old values
        for w in reversed(self.weights):
            acc = None
            for n, cf in enumerate(coefs, start=1):
                if cf == 0:
                    continue
                src = tuple(x - n * y for x, y in zip(w, om))
                if min(src) < 0:
                    break
                a = self.A[src]
                if not a.any():
                    continue
                s = self._shift(a, dg, dB, dl, dr, n)
                if s is None:
                    continue
                if acc is None:
                    acc = cf * s
                else:
                    acc += cf * s
            if acc is not None:
                self.A[w] = self.A[w] + acc


def chain_dims(alg, lam, k, ghs, N, extra_charges=(), B_pad=12, l_rng=16, verbose=False):
    """dim V_{charge, g} of the meson complex chains (bg_complex_k.chains conventions) in sector lam.
    Charges: T = (N..N; N..N) and extra_charges given as (Lc tuple, Rc tuple).  k in (1, 2).
    Returns {(Lc, Rc, g): dim}."""
    lam = tuple(lam)
    wc = word_class_counts(alg, lam)
    cc = cyclic_class_counts(alg, lam)
    tcl = []
    for (om, l), c in sorted(cc.items()):
        th = _theta(om)
        tcl.append((om, l - th, (l + th) & 1, c))
    wcl = []
    for (om, l), c in sorted(wc.items()):
        th = _theta(om)
        wcl.append((om, l - th, (l + th) & 1, c))
    return chain_dims_classes(lam, k, ghs, N, tcl, wcl, extra_charges, B_pad, l_rng)


def chain_dims_classes(lam, k, ghs, N, tclasses, wclasses, extra_charges=(), B_pad=12, l_rng=16):
    """Same count for arbitrary generator classes.
    tclasses: [(omega, gh, parity, count)] trace generators (no charge);
    wclasses: [(omega, gh_w, parity_w, count)] nonempty meson words; the mesons M, P, MM, R on a word
    have parity p_w, 1+p_w, 1+p_w, p_w and gh gh_w, gh_w - 1, gh_w - 1, gh_w - 2; M_0 is added."""
    lam = tuple(lam)
    gpos = sum(max(0, g) * c for (om, g, p, c) in wclasses if g > 0)
    ntot = sum(lam[3:])
    glo = min(ghs) - max(ntot, 0) - 0
    ghi = max(ntot, 0)
    S = _Series(lam, k, (glo, ghi), (-B_pad, k * N + B_pad), l_rng)
    for (om, gh, odd, c) in tclasses:
        S.mult_class(om, odd, c, [(gh, 0, 0, 0)])
    for (om, gw, pw, c) in wclasses:
        for t in range(1, 5):
            odd = pw if t in (1, 4) else 1 - pw
            gh = gw - (0, 0, 1, 1, 2)[t]
            shifts = []
            if k == 1:
                dB = {1: 1, 2: 0, 3: 0, 4: -1}[t]
                shifts.append((gh, dB, 0, 0))
            else:
                for a in range(2):
                    for b in range(2):
                        sa = 1 if a == 0 else -1
                        sb = 1 if b == 0 else -1
                        if t == 1:      # M: Lc_a +1, Rc_b +1
                            shifts.append((gh, 1, sa, sb))
                        elif t == 2:    # P: Lc_a +1, Lc_b -1
                            shifts.append((gh, 0, sa - sb, 0))
                        elif t == 3:    # MM: Rc_a -1, Rc_b +1
                            shifts.append((gh, 0, 0, sb - sa))
                        else:           # R: Rc_a -1, Lc_b -1
                            shifts.append((gh, -1, -sb, -sa))
            S.mult_class(om, odd, c, shifts)
    out = {}
    A = S.A[lam]
    targets = [((N,) * k, (N,) * k)] + [tuple(x) for x in extra_charges]
    for (Lc, Rc) in targets:
        Btot = sum(Lc)
        assert Btot == sum(Rc)
        for g in ghs:
            gi = g - S.g0
            tot = 0.0
            for Bnz in range(S.B0, S.B1 + 1):
                n0 = Btot - Bnz
                if n0 < 0:
                    continue
                bi = Bnz - S.B0
                if k == 1:
                    tot += A[gi, bi]
                else:
                    lt = Lc[0] - Lc[1]
                    rt = Rc[0] - Rc[1]
                    for l in range(-S.L, S.L + 1):
                        dl = lt - l
                        for r in range(-S.L, S.L + 1):
                            dr = rt - r
                            v = A[gi, bi, l + S.L, r + S.L]
                            if v == 0:
                                continue
                            tot += v * _m0_count(n0, dl, dr)
            out[(tuple(Lc), tuple(Rc), g)] = tot
    return out


def small_model_classes(lam):
    """Generator classes of the HA^! model in a sector with n+ = n- = 0:
    HA^! = C[x1,x2,x3] (x_I = [gamma^I], gh 0, even); T generators = reduced HC(C[x]):
    HC_n at weight a has dim binom(s(a)-1, n), gh -n, parity n mod 2 (s(a) = #nonzero entries)."""
    import itertools
    assert lam[3] == 0 and lam[4] == 0
    tcl, wcl = [], []
    for a in itertools.product(*[range(x + 1) for x in lam[:3]]):
        if sum(a) == 0:
            continue
        om = tuple(a) + (0, 0)
        s = sum(1 for x in a if x)
        for n in range(0, s):
            c = comb(s - 1, n)
            if c:
                tcl.append((om, -n, n & 1, c))
        wcl.append((om, 0, 0, 1))
    return tcl, wcl


_M0 = {}


def _m0_count(n0, dl, dr):
    key = (n0, dl, dr)
    if key in _M0:
        return _M0[key]
    c = 0
    for m11 in range(n0 + 1):
        for m12 in range(n0 + 1 - m11):
            for m21 in range(n0 + 1 - m11 - m12):
                m22 = n0 - m11 - m12 - m21
                if (m11 + m12) - (m21 + m22) == dl and (m11 + m21) - (m12 + m22) == dr:
                    c += 1
    _M0[key] = c
    return c


def chain_dims_general(alg, lam, k, ghs, N, charges=None):
    """Exact dim V_{charge, g} for any k by a sparse DP over generator classes (small sectors only).
    charges: list of (Lc tuple, Rc tuple); default [T].  Returns {(Lc, Rc, g): dim} (Python ints)."""
    import itertools
    lam = tuple(lam)
    if charges is None:
        charges = [((N,) * k, (N,) * k)]
    wc = word_class_counts(alg, lam)
    cc = cyclic_class_counts(alg, lam)
    units = sum(lam)
    gmin = min(ghs) - sum(lam[3:])
    # state key: (omega, Lc, Rc, gh)
    S = {((0,) * len(lam), (0,) * k, (0,) * k, 0): 1}

    def mult(om, odd, c, dLc, dRc, dg):
        nonlocal S
        nmax = min((lam[i] // om[i]) for i in range(len(om)) if om[i] > 0)
        coefs = [comb(c, n) if odd else comb(n + c - 1, n) for n in range(1, nmax + 1)]
        new = dict(S)
        for (w, Lc, Rc, g), v in S.items():
            for n, cf in enumerate(coefs, start=1):
                if not cf:
                    continue
                nw = tuple(x + n * y for x, y in zip(w, om))
                if not all(a <= b for a, b in zip(nw, lam)):
                    break
                ng = g + n * dg
                if ng < gmin:
                    break
                nL = tuple(x + n * y for x, y in zip(Lc, dLc))
                nR = tuple(x + n * y for x, y in zip(Rc, dRc))
                if max(nL + nR) > N + 1 + units or min(nL + nR) < -units - 1:
                    continue
                key = (nw, nL, nR, ng)
                new[key] = new.get(key, 0) + v * cf
        S = new

    for (om, l), c in sorted(cc.items()):
        th = _theta(om)
        mult(om, (l + th) & 1, c, (0,) * k, (0,) * k, l - th)
    for (om, l), c in sorted(wc.items()):
        th = _theta(om)
        pw = (l + th) & 1
        gw = l - th
        for t in range(1, 5):
            odd = pw if t in (1, 4) else 1 - pw
            gh = gw - (0, 0, 1, 1, 2)[t]
            for a in range(k):
                for b in range(k):
                    dL = [0] * k
                    dR = [0] * k
                    if t == 1:
                        dL[a] += 1; dR[b] += 1
                    elif t == 2:
                        dL[a] += 1; dL[b] -= 1
                    elif t == 3:
                        dR[a] -= 1; dR[b] += 1
                    else:
                        dR[a] -= 1; dL[b] -= 1
                    mult(om, odd, c, tuple(dL), tuple(dR), gh)
    out = {}
    for (Lc_t, Rc_t) in charges:
        for g in ghs:
            tot = 0
            for (w, Lc, Rc, gg), v in S.items():
                if w != lam or gg != g:
                    continue
                r = [x - y for x, y in zip(Lc_t, Lc)]
                cc_ = [x - y for x, y in zip(Rc_t, Rc)]
                if min(r + cc_) < 0 or sum(r) != sum(cc_):
                    continue
                tot += v * _contingency(tuple(r), tuple(cc_))
            out[(tuple(Lc_t), tuple(Rc_t), g)] = tot
    return out


_CT = {}


def _contingency(r, c):
    """Number of nonnegative integer matrices with row sums r and column sums c."""
    key = (r, c)
    if key in _CT:
        return _CT[key]
    if len(r) == 1:
        return 1 if r[0] == sum(c) else 0
    tot = 0
    def rec(j, rem, row):
        nonlocal tot
        if j == len(c) - 1:
            if rem <= c[j]:
                nc = tuple(c[i] - row[i] for i in range(len(c) - 1)) + (c[j] - rem,)
                tot_add = _contingency(r[1:], nc)
                tot_list.append(tot_add)
            return
        for x in range(min(rem, c[j]) + 1):
            row.append(x)
            rec(j + 1, rem - x, row)
            row.pop()
    tot_list = []
    rec(0, r[0], [])
    tot = sum(tot_list)
    _CT[key] = tot
    return tot
