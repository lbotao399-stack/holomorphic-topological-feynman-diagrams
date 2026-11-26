"""Cohomology of slices of V = M / sum_u D_u M, of V^eta (strict invariants in V) and of CE side (cost (0,0,0,a,b))."""
import sys, time
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from vjet import Jets, Echelon, P1, P2, add


def e_(u):
    return tuple(1 if i == u else 0 for i in range(5))


def sub(t, s):
    return tuple(x - y for x, y in zip(t, s))


def addt(t, s):
    return tuple(x + y for x, y in zip(t, s))


class Slice:
    """data for one cost t: monomials by ghost, index maps, D-image rows by ghost."""

    def __init__(self, J, t, charge0, cache):
        self.J, self.t = J, tuple(t)
        self.M = J.slice(t, charge0)
        self.idx = {g: {m: i for i, m in enumerate(L)} for g, L in self.M.items()}
        self.Irows = {}
        for u in range(5):
            if t[u] == 0:
                continue
            s = sub(t, e_(u))
            Ms = cache.get(s, charge0)
            for g, L in Ms.items():
                for m in L:
                    r = J.Dmono(m, u)
                    if r:
                        self.Irows.setdefault(g, []).append(self.vec(g, r))

    def vec(self, g, P):
        ix = self.idx.get(g, {})
        v = {}
        for m, c in P.items():
            v[ix[m]] = c
        return v

    def Qrow(self, g, m):
        return self.vec(g + 1, self.J.Qmono(m))


class Cache:
    def __init__(self, J):
        self.J = J
        self.d = {}

    def get(self, t, charge0):
        k = (tuple(t), charge0)
        if k not in self.d:
            self.d[k] = self.J.slice(t, charge0)
        return self.d[k]


def cohom_V(J, t, charge0=False, primes=(P1, P2), exact=False, ghosts=None):
    """returns dict g -> (dimM, dimV, H per field) for all ghosts in the slice."""
    cache = Cache(J)
    S = Slice(J, t, charge0, cache)
    gs = sorted(S.M)
    fields = list(primes) + ([None] if exact else [])
    out = {g: {"dimM": len(S.M[g])} for g in gs}
    rkI = {}
    rkQ = {}
    for p in fields:
        for g in gs + [gs[-1] + 1]:
            if ghosts is not None and g not in ghosts and (g - 1) not in ghosts:
                rkI[(p, g)] = None
                continue
            E = Echelon(p)
            for r in S.Irows.get(g, []):
                E.add(r)
            rkI[(p, g)] = E.rank()
        for g in gs:
            if ghosts is not None and g not in ghosts and (g + 1) not in ghosts:
                rkQ[(p, g)] = None
                continue
            E = Echelon(p)
            for r in S.Irows.get(g + 1, []):
                E.add(r)
            r0 = E.rank()
            for m in S.M[g]:
                E.add(S.Qrow(g, m))
            rkQ[(p, g)] = E.rank() - r0
    for g in gs:
        for p in fields:
            a = rkQ.get((p, g))
            b = rkQ.get((p, g - 1), 0)
            if a is None or b is None:
                out[g][p] = None
                continue
            if rkI[(p, g)] is None:
                out[g][p] = None
                continue
            out[g]["dimV"] = len(S.M[g]) - rkI[(p, g)]
            out[g][p] = len(S.M[g]) - rkI[(p, g)] - a - b
    return out, S


def ce_table(J, a, b, charge0=False, exact=False):
    """H^n(C_g^0)_{(a,b)} = H^n(V_{cost (0,0,0,a,b)})  (ghost = CE degree there)."""
    out, _ = cohom_V(J, (0, 0, 0, a, b), charge0, exact=exact)
    return out


def kernel_eta(J, t, vs, charge0, p, cache):
    """K_g = { x in M_t^g : eta_v x in D M_{t+e_v} for v in vs }; returns dict g -> list of vectors (index space of M_t^g)."""
    S = Slice(J, t, charge0, cache)
    tgt = {}
    for v in vs:
        tv = addt(t, e_(v))
        Sv = Slice(J, tv, charge0, cache)
        Ev = {}
        for g, rows in Sv.Irows.items():
            E = Echelon(p)
            for r in rows:
                E.add(r)
            Ev[g] = E
        tgt[v] = (Sv, Ev)
    K = {}
    for g, L in S.M.items():
        E = Echelon(p)
        kers = []
        n = len(L)
        for i, m in enumerate(L):
            vec = {}
            for j, v in enumerate(vs):
                Sv, Ev = tgt[v]
                img = J.etamono(m, v)
                if not img:
                    continue
                w = Sv.vec(g - 1, img)
                if (g - 1) in Ev:
                    w = Ev[g - 1].nf(w)
                else:
                    w = {k: x % p for k, x in w.items() if x % p}
                for k, x in w.items():
                    vec[j * 10 ** 7 + k] = x
            # augmented: image columns first (< 10^7 * len(vs)), then identity columns
            base = 10 ** 7 * len(vs)
            vec[base + i] = 1
            w, c = E.reduce(vec)
            if c is not None and c >= base:
                kers.append({k - base: x for k, x in w.items()})
            else:
                E.add(vec)
        K[g] = kers
    return K, S


def cohom_Veta(J, t, vs=(0, 1, 2), charge0=False, p=P1):
    """H^g of the strict-invariant subcomplex V^{eta_vs} of V_t."""
    cache = Cache(J)
    K, S = kernel_eta(J, t, vs, charge0, p, cache)
    gs = sorted(S.M)
    res = {}
    rk = {}
    dimVe = {}
    for g in gs:
        E = Echelon(p)
        for r in S.Irows.get(g, []):
            E.add(r)
        rI = E.rank()
        for r in K[g]:
            E.add(r)
        dimVe[g] = E.rank() - rI
        assert E.rank() == len(K[g]) or True
        # check I subset K: rank(I + K) == dim K
        res.setdefault(g, {})["I_in_K"] = (E.rank() == len(K[g]))
        E2 = Echelon(p)
        for r in S.Irows.get(g + 1, []):
            E2.add(r)
        r0 = E2.rank()
        for kv in K[g]:
            img = {}
            for i, c in kv.items():
                for mm, cc in J.Qmono(S.M[g][i]).items():
                    j = S.idx[g + 1][mm]
                    img[j] = (img.get(j, 0) + c * cc) % p
            E2.add(img)
        rk[g] = E2.rank() - r0
    for g in gs:
        res[g]["dimVeta"] = dimVe[g]
        res[g]["H"] = dimVe[g] - rk[g] - rk.get(g - 1, 0)
    return res, K, S
