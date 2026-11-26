"""Slot spaces of C_quot (w = 0, 1), exact ranks of Q, cohomology of C_quot and of V = C_quot / sum_k D_k C_quot.

Slot = (w, t, N), N = #X + #beta (Q raises N by 1, parity = N mod 2), G = -w - sum(t) (redundant).
Exact arithmetic: Echelon (integer fraction-free elimination with content removal) from jets.py.
"""
import itertools, time
from collections import defaultdict
from jets import Alg, Echelon, EchelonModP, vec, multi_indices, x_letters_quot


class Space:
    def __init__(self, n, w, Gmax_high):
        self.n = n
        self.w = w
        self.A = Alg(n, quot=True, sX=-1)
        A = self.A
        levels = x_letters_quot(A, Gmax_high + 1)
        self.levels = levels
        low = levels[-1] + levels[0]
        self.low_idx = defaultdict(list)
        for k in range(len(low) + 1):
            for c in itertools.combinations(low, k):
                self.low_idx[A.weight(c)].append(c)
        self.high_idx = defaultdict(list)
        highlev = [l for l in sorted(levels) if l >= 1 and l <= Gmax_high]
        letters = [(l, x) for l in highlev for x in levels[l]]
        self.Gmax_high = Gmax_high

        def rec(start, budget, chosen):
            self.high_idx[A.weight(tuple(chosen))].append(tuple(chosen))
            for j in range(start, len(letters)):
                l, x = letters[j]
                if l <= budget:
                    chosen.append(x)
                    rec(j + 1, budget - l, chosen)
                    chosen.pop()
        rec(0, Gmax_high, [])
        # w letters
        self.wl = []
        if w == 1:
            for m in range(Gmax_high + n + 3):
                for a in multi_indices(n, m):
                    self.wl.append(A.L('b', 0, a))
                    for i in range(n):
                        self.wl.append(A.L('g', i, a))
        self.cache = {}

    def G_of(self, t):
        return -self.w - sum(t)

    def slot(self, t):
        """dict N -> sorted list of monomials of weight t."""
        r = self.cache.get(t)
        if r is not None:
            return r
        A = self.A
        n = self.n
        G = self.G_of(t)
        out = defaultdict(list)
        wls = self.wl if self.w == 1 else [None]
        for wl in wls:
            if wl is None:
                Gw, tw = 0, (0,) * n
            else:
                Gw, tw = A.Gw[wl], A.tw[wl]
            if G - Gw < -n:
                continue
            for (Gl, tl), lows in self.low_idx.items():
                Gh = G - Gw - Gl
                if Gh < 0:
                    continue
                assert Gh <= self.Gmax_high, (G, Gw, Gl, Gh)
                th = tuple(t[j] - tw[j] - tl[j] for j in range(n))
                highs = self.high_idx.get((Gh, th))
                if not highs:
                    continue
                for lo in lows:
                    for hi in highs:
                        m = lo + hi + ((wl,) if wl is not None else ())
                        m = tuple(sorted(m))
                        out[A.Ndeg(m)].append(m)
        r = {N: sorted(v) for N, v in out.items()}
        self.cache[t] = r
        return r

    def drop(self, t):
        self.cache.pop(t, None)


class Cohom:
    def __init__(self, space, verbose=False):
        self.S = space
        self.A = space.A
        self.rQ = {}
        self.rI = {}
        self.rQbar = {}
        self.verbose = verbose
        self.timing = {}

    def basis(self, t, N):
        return self.S.slot(t).get(N, [])

    def rankQ(self, t, N):
        key = (t, N)
        if key in self.rQ:
            return self.rQ[key]
        src = self.basis(t, N)
        tgt = self.basis(t, N + 1)
        if not src or not tgt:
            self.rQ[key] = 0
            return 0
        idx = {m: j for j, m in enumerate(tgt)}
        t0 = time.time()
        E = Echelon()
        A = self.A
        for m in src:
            E.insert(vec(A.Qmono(m), idx))
        self.rQ[key] = E.rank
        self.timing[('Q',) + key] = time.time() - t0
        return E.rank

    def HC(self, t):
        sl = self.S.slot(t)
        if not sl:
            return {}
        Ns = range(min(sl), max(sl) + 1)
        out = {}
        for N in Ns:
            d = len(sl.get(N, []))
            h = d - self.rankQ(t, N) - self.rankQ(t, N - 1)
            out[N] = (d, h)
        return out

    def shiftw(self, t, k, s=1):
        return tuple(t[j] + s * (j == k) for j in range(len(t)))

    def V_target(self, t, N):
        """echelon in C^N_t of the D-images, then the Q-images of C^{N-1}_t."""
        key = (t, N)
        if key in self.rI:
            return
        A = self.A
        tgt = self.basis(t, N)
        idx = {m: j for j, m in enumerate(tgt)}
        E = Echelon()
        t0 = time.time()
        if tgt:
            for k in range(self.S.n):
                for m in self.basis(self.shiftw(t, k), N):
                    E.insert(vec(A.Dq({m: 1}, k), idx))
        rI = E.rank
        if tgt:
            for m in self.basis(t, N - 1):
                E.insert(vec(A.Qmono(m), idx))
        self.rI[key] = rI
        self.rQbar[(t, N - 1)] = E.rank - rI
        self.timing[('V',) + key] = time.time() - t0

    def HV(self, t):
        sl = self.S.slot(t)
        Ns = set(sl)
        for k in range(self.S.n):
            Ns |= set(self.S.slot(self.shiftw(t, k)))
        if not Ns:
            return {}
        out = {}
        for N in range(min(Ns) - 1, max(Ns) + 2):
            self.V_target(t, N)
        for N in range(min(Ns), max(Ns) + 1):
            d = len(sl.get(N, []))
            dV = d - self.rI[(t, N)]
            h = dV - self.rQbar[(t, N)] - self.rQbar.get((t, N - 1), 0)
            out[N] = (d, dV, h)
        return out


def koszul_dimV(space, t, N):
    """sum_K (-1)^|K| dim C^N_{t+e_K} (equals dim V^N_t if C is free over C[D] in this slot)."""
    n = space.n
    s = 0
    for r in range(n + 1):
        for K in itertools.combinations(range(n), r):
            tt = tuple(t[j] + (j in K) for j in range(n))
            s += (-1) ** r * len(space.slot(tt).get(N, []))
    return s


def all_weights(space, G):
    """all t with sum(t) = -w-G occurring in C_quot (w fixed by the space)."""
    A = space.A
    n = space.n
    if not hasattr(space, 'high_by_G'):
        space.high_by_G = defaultdict(set)
        for (Gh, th) in space.high_idx:
            space.high_by_G[Gh].add(th)
    ts = set()
    wls = space.wl if space.w == 1 else [None]
    for wl in wls:
        if wl is None:
            Gw, tw = 0, (0,) * n
        else:
            Gw, tw = A.Gw[wl], A.tw[wl]
        if G - Gw < -n:
            continue
        for (Gl, tl) in space.low_idx:
            Gh = G - Gw - Gl
            if Gh < 0 or Gh > space.Gmax_high:
                continue
            for th in space.high_by_G.get(Gh, ()):
                ts.add(tuple(tw[j] + tl[j] + th[j] for j in range(n)))
    return sorted(ts)
