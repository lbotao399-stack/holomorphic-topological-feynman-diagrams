"""Budzik-Gaiotto meson complex C_.(gl'_{k|k}[A^!] | b_{k|k}; Ber^N) at k = 1 (arXiv:2611.18539,
(93)-(95), (102), (117)-(120)).

Chains: graded polynomials in traces T_w = Tr W_w and mesons
    M_w = psibar W_w psi,  P_w = psibar W_w u,  MM_w = ubar W_w psi,  R_w = ubar W_w u,
psi, psibar odd (gh 0), u, ubar even (gh -1; u antifield of psibar, ubar of psi).
d = Q + Delta_BV, Delta_BV = Tr d_psibar d_u + Tr d_ubar d_psi (left derivatives, Tr 1 = N).
Relative to b_{1|1} = span(P_0, MM_0, R_0):  drop P_0, MM_0, R_0; b_0-coinvariants with Ber^N fix
#M - #R = N; b_1-coinvariants quotient by D_r = {R_0, .}.
    R_N = C'_N / D_r(C'_{N+1}),   H_g(R_N),  g = ghost number (d raises g by 1).
"""
from . import rank_mod
from .uN_q0 import sort_sign

T, M, P, MM, R = 0, 1, 2, 3, 4
NAMES = {T: "T", M: "M", P: "P", MM: "MM", R: "R"}
START_PB = {M, P}          # starts with psibar (odd)
START_UB = {MM, R}         # starts with ubar (even)
END_PSI = {M, MM}          # ends with psi (odd)
END_U = {P, R}             # ends with u (even)


def mtype(start_pb, end_psi):
    if start_pb:
        return M if end_psi else P
    return MM if end_psi else R


class BG1:
    def __init__(self, alg):
        self.alg = alg
        self._par = {}

    # ---------------------------------------------------------------- generator data
    def par(self, g):
        p = self._par.get(g)
        if p is None:
            t, w = g
            pw = self.alg.pword(w)
            p = pw if t in (T, M, R) else 1 - pw
            self._par[g] = p
        return p

    def gh(self, g):
        t, w = g
        x = self.alg.ghword(w)
        return x - {T: 0, M: 0, P: 1, MM: 1, R: 2}[t]

    def ber(self, g):
        t, w = g
        return 1 if t == M else (-1 if t == R else 0)

    def is_zero_mode_b(self, g):
        return g[0] in (P, MM, R) and len(g[1]) == 0

    def canon(self, gens):
        return sort_sign(gens, self.par)

    # ---------------------------------------------------------------- Q
    def Q_gen(self, g):
        t, w = g
        out = []
        for nw, c in self.alg.Qword(w).items():
            if t == T:
                s, cw = self.alg.canon_cyclic(nw)
                if s:
                    out.append(((T, cw), c * s))
            elif t in START_PB:
                out.append(((t, nw), -c))
            else:
                out.append(((t, nw), c))
        return out

    def Q(self, gens):
        out = {}
        pre = 0
        for pos, g in enumerate(gens):
            s0 = -1 if pre else 1
            for ng, c in self.Q_gen(g):
                new = list(gens)
                new[pos] = ng
                sg, key = self.canon(new)
                if sg:
                    out[key] = out.get(key, 0) + s0 * c * sg
            pre ^= self.par(g)
        return out

    # ---------------------------------------------------------------- Delta_BV
    def Delta(self, gens, N, only=None):
        """Delta_BV on the ordered product gens.  only=i keeps contractions involving position i."""
        out = {}
        n = len(gens)
        pars = [self.par(g) for g in gens]
        prefix = [0] * (n + 1)
        for i in range(n):
            prefix[i + 1] = prefix[i] + pars[i]

        def emit(newgens, coef):
            sg, key = self.canon(newgens)
            if sg:
                out[key] = out.get(key, 0) + coef * sg

        # type 1: d_psibar d_u : psibar at start of gens[a], u at end of gens[b]
        for b in range(n):
            tb, wb = gens[b]
            if tb not in END_U:
                continue
            for a in range(n):
                ta, wa = gens[a]
                if ta not in START_PB:
                    continue
                if only is not None and only not in (a, b):
                    continue
                s = -1 if prefix[a] & 1 else 1
                if a == b:
                    rest = gens[:a] + gens[a + 1:]
                    if len(wa) == 0:
                        if N:
                            emit(rest, s * N)
                    else:
                        sc, cw = self.alg.canon_cyclic(wa)
                        if sc:
                            emit(gens[:a] + [(T, cw)] + gens[a + 1:], s * sc)
                    continue
                pA = 1 - pars[a]
                pB = pars[b]
                newg = (mtype(tb in START_PB, ta in END_PSI), wb + wa)
                if a < b:
                    mid = prefix[b] - prefix[a + 1]
                    sm = -1 if (pA * (mid + pB)) & 1 else 1
                else:
                    mid = prefix[a] - prefix[b + 1]
                    sm = -1 if (pA * mid) & 1 else 1
                new = list(gens)
                new[b] = newg
                del new[a]
                emit(new, s * sm)
        # type 2: d_ubar d_psi : psi at end of gens[a], ubar at start of gens[b]
        for a in range(n):
            ta, wa = gens[a]
            if ta not in END_PSI:
                continue
            s = -1 if (prefix[a] + pars[a] - 1) & 1 else 1
            for b in range(n):
                tb, wb = gens[b]
                if tb not in START_UB:
                    continue
                if only is not None and only not in (a, b):
                    continue
                if a == b:
                    rest = gens[:a] + gens[a + 1:]
                    if len(wa) == 0:
                        if N:
                            emit(rest, s * N)
                    else:
                        sc, cw = self.alg.canon_cyclic(wa)
                        if sc:
                            emit(gens[:a] + [(T, cw)] + gens[a + 1:], s * sc)
                    continue
                pA = 1 - pars[a]
                pB = pars[b]
                newg = (mtype(ta in START_PB, tb in END_PSI), wa + wb)
                if b > a:
                    mid = prefix[b] - prefix[a + 1]
                    sm = -1 if (pB * mid) & 1 else 1
                else:
                    mid = prefix[a] - prefix[b + 1]
                    sm = -1 if (pB * (mid + pA)) & 1 else 1
                new = list(gens)
                new[a] = newg
                del new[b]
                emit(new, s * sm)
        return {k: v for k, v in out.items() if v}

    def d(self, mono, N):
        gens = list(mono)
        out = self.Q(gens)
        for k, v in self.Delta(gens, N).items():
            out[k] = out.get(k, 0) + v
        return {k: v for k, v in out.items() if v}

    def Dr(self, mono):
        """D_r = {R_0, .}: cross contractions of R_0 (placed first) with mono, zero modes dropped."""
        gens = [(R, ())] + list(mono)
        raw = self.Delta(gens, 0, only=0)
        return {k: v for k, v in raw.items() if not any(self.is_zero_mode_b(g) for g in k)}

    # ---------------------------------------------------------------- sector enumeration
    def generators(self, lam):
        alg = self.alg
        words = alg.words_upto(lam)
        gens = []
        cyc = set()
        for w in words:
            s, cw = alg.canon_cyclic(w)
            if s:
                cyc.add(cw)
        gens += [(T, w) for w in sorted(cyc)]
        for t in (M, P, MM, R):
            gens += [(t, w) for w in words]
        return gens

    def chains(self, lam, ber):
        """Monomials of C'_ber with weight lam: dict gh -> list of monomials (canonical tuples)."""
        alg = self.alg
        lam = tuple(lam)
        gens = self.generators(lam)
        W = [alg.wword(g[1]) for g in gens]
        Pp = [self.par(g) for g in gens]
        B = [self.ber(g) for g in gens]
        out = {}
        M0 = (M, ())

        def rec(start, cur, rem, bsum):
            if all(r == 0 for r in rem):
                a = ber - bsum
                if a >= 0:
                    mono = list(cur) + [M0] * a
                    sg, key = self.canon(mono)
                    if sg:
                        gh = sum(self.gh(g) for g in key)
                        out.setdefault(gh, []).append(key)
                return
            for idx in range(start, len(gens)):
                wi = W[idx]
                if all(x <= y for x, y in zip(wi, rem)):
                    nrem = tuple(y - x for x, y in zip(wi, rem))
                    cur.append(gens[idx])
                    rec(idx + 1 if Pp[idx] else idx, cur, nrem, bsum + B[idx])
                    cur.pop()
        rec(0, [], lam, 0)
        for gh in out:
            out[gh] = sorted(set(out[gh]))
        return out


def homology(bg, lam, N, field=rank_mod.P1, verbose=False):
    """dim H_g of R_N = C'_N / D_r C'_{N+1} in sector lam.  Returns {g: (dim C'_N, dim S, rank d, dim H)}."""
    CN = bg.chains(lam, N)
    CN1 = bg.chains(lam, N + 1)
    ghs = sorted(set(CN) | set(g - 1 for g in CN1))
    # S_g = D_r(C'_{N+1, g+1}) inside C'_{N,g}
    Sspan = {}
    for g in ghs:
        sp = rank_mod.make_span(field)
        for m in CN1.get(g + 1, []):
            v = bg.Dr(m)
            if v:
                sp.add(v)
        Sspan[g] = sp
    dimS = {g: Sspan[g].dim() for g in ghs}
    rankd = {}
    for g in ghs:
        # rank of d'_g : C'_{N,g}/S_g -> C'_{N,g+1}/S_{g+1}
        base = rank_mod.make_span(field)
        if g + 1 in Sspan:
            for piv in Sspan[g + 1].piv.values():
                base.add(piv)
        r0 = base.dim()
        for m in CN.get(g, []):
            v = bg.d(m, N)
            if v:
                base.add(v)
        rankd[g] = base.dim() - r0
    out = {}
    for g in ghs:
        dimC = len(CN.get(g, []))
        dimR = dimC - dimS[g]
        h = dimR - rankd.get(g, 0) - rankd.get(g - 1, 0)
        out[g] = (dimC, dimS[g], rankd.get(g, 0), h)
    return out


def check_d2(bg, lam, N, ber, field=rank_mod.P1, absolute=False, maxn=200):
    """d^2 = 0 on chains of sector lam (relative C' unless absolute)."""
    ch = bg.chains(lam, ber)
    bad = 0
    tested = 0
    for g, monos in ch.items():
        for m in monos[:maxn]:
            v = bg.d(m, N)
            acc = {}
            for k, c in v.items():
                if not absolute and any(bg.is_zero_mode_b(x) for x in k):
                    continue
                for k2, c2 in bg.d(k, N).items():
                    acc[k2] = acc.get(k2, 0) + c * c2
            acc = {k: x for k, x in acc.items() if x % field}
            tested += 1
            if acc:
                bad += 1
    return tested, bad
