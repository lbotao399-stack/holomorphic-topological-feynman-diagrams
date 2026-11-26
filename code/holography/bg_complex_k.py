"""Budzik-Gaiotto meson complex C_.(gl'_{k|k}[A^!] | b_{k|k}; Ber^N) for general k (here k = 1, 2).

Generators (t, w, a, b): T_w = (0, w, 0, 0); mesons with flavour a of the left field and b of the right:
    M = psibar_a W psi_b, P = psibar_a W u_b, MM = ubar_a W psi_b, R = ubar_a W u_b.
Delta_BV = sum_a Tr d_{psibar_a} d_{u_a} + sum_a Tr d_{ubar_a} d_{psi_a}.
Charges: Lc[x] = #psibar_x - #u_x, Rc[x] = #psi_x - #ubar_x.
b_{k|k} = span(P_0, MM_0, R_0): the zero modes are dropped and the complex is replaced by the
b-coinvariants of the Ber^N-twisted action:
    R_N = V_T / ( sum_{x!=y} rho(P0_xy) V_{T-ch} + sum_{x!=y} rho(MM0_xy) V_{T-ch} + sum_{x,y} rho(R0_xy) V_{T-ch} ),
    T = (Lc, Rc) = (N,...,N; N,...,N);  the diagonal generators act by N - Lc[x], N - Rc[x] (zero at T).
"""
from itertools import combinations_with_replacement
from . import rank_mod
from .uN_q0 import sort_sign

T, M, P, MM, R = 0, 1, 2, 3, 4
START_PB = {M, P}
START_UB = {MM, R}
END_PSI = {M, MM}
END_U = {P, R}


def mtype(start_pb, end_psi):
    if start_pb:
        return M if end_psi else P
    return MM if end_psi else R


class BGk:
    def __init__(self, alg, k):
        self.alg = alg
        self.k = k
        self._par = {}

    def par(self, g):
        p = self._par.get(g)
        if p is None:
            t, w = g[0], g[1]
            pw = self.alg.pword(w)
            p = pw if t in (T, M, R) else 1 - pw
            self._par[g] = p
        return p

    def gh(self, g):
        return self.alg.ghword(g[1]) - (0, 0, 1, 1, 2)[g[0]]

    def charge(self, g):
        k = self.k
        Lc = [0] * k
        Rc = [0] * k
        t, w, a, b = g
        if t == M:
            Lc[a] += 1; Rc[b] += 1
        elif t == P:
            Lc[a] += 1; Lc[b] -= 1
        elif t == MM:
            Rc[a] -= 1; Rc[b] += 1
        elif t == R:
            Rc[a] -= 1; Lc[b] -= 1
        return tuple(Lc) + tuple(Rc)

    def is_b_mode(self, g):
        return g[0] in (P, MM, R) and len(g[1]) == 0

    def canon(self, gens):
        return sort_sign(gens, self.par)

    # ---------------------------------------------------------------- Q
    def Q_gen(self, g):
        t, w, a, b = g
        out = []
        for nw, c in self.alg.Qword(w).items():
            if t == T:
                s, cw = self.alg.canon_cyclic(nw)
                if s:
                    out.append(((T, cw, 0, 0), c * s))
            elif t in START_PB:
                out.append(((t, nw, a, b), -c))
            else:
                out.append(((t, nw, a, b), c))
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

        for b in range(n):
            tb, wb, fab, fbb = gens[b]
            if tb not in END_U:
                continue
            for a in range(n):
                ta, wa, faa, fba = gens[a]
                if ta not in START_PB or faa != fbb:
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
                            emit(gens[:a] + [(T, cw, 0, 0)] + gens[a + 1:], s * sc)
                    continue
                pA = 1 - pars[a]
                pB = pars[b]
                newg = (mtype(tb in START_PB, ta in END_PSI), wb + wa, fab, fba)
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
        for a in range(n):
            ta, wa, faa, fba = gens[a]
            if ta not in END_PSI:
                continue
            s = -1 if (prefix[a] + pars[a] - 1) & 1 else 1
            for b in range(n):
                tb, wb, fab, fbb = gens[b]
                if tb not in START_UB or fba != fab:
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
                            emit(gens[:a] + [(T, cw, 0, 0)] + gens[a + 1:], s * sc)
                    continue
                pA = 1 - pars[a]
                pB = pars[b]
                newg = (mtype(ta in START_PB, tb in END_PSI), wa + wb, faa, fbb)
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
        for kk, v in self.Delta(gens, N).items():
            out[kk] = out.get(kk, 0) + v
        return {kk: v for kk, v in out.items() if v}

    def rho(self, beta, mono):
        """{beta, mono}: cross contractions with beta placed first; b-modes dropped from the result."""
        raw = self.Delta([beta] + list(mono), 0, only=0)
        return {kk: v for kk, v in raw.items() if not any(self.is_b_mode(g) for g in kk)}

    # ---------------------------------------------------------------- enumeration
    def generators(self, lam):
        alg = self.alg
        k = self.k
        words = alg.words_upto(lam)
        cyc = set()
        for w in words:
            s, cw = alg.canon_cyclic(w)
            if s:
                cyc.add(cw)
        gens = [(T, w, 0, 0) for w in sorted(cyc)]
        for t in (M, P, MM, R):
            for w in words:
                for a in range(k):
                    for b in range(k):
                        gens.append((t, w, a, b))
        return gens

    def chains(self, lam, bers):
        """Monomials of weight lam with #M - #R in bers: dict (charge, gh) -> sorted list."""
        alg = self.alg
        k = self.k
        lam = tuple(lam)
        gens = self.generators(lam)
        W = [alg.wword(g[1]) for g in gens]
        Pp = [self.par(g) for g in gens]
        B = [1 if g[0] == M else (-1 if g[0] == R else 0) for g in gens]
        M0s = [(M, (), a, b) for a in range(k) for b in range(k)]
        out = {}

        def finish(cur, bsum):
            for ber in bers:
                a = ber - bsum
                if a < 0:
                    continue
                for combo in combinations_with_replacement(M0s, a):
                    mono = list(cur) + list(combo)
                    sg, key = self.canon(mono)
                    if not sg:
                        continue
                    ch = [0] * (2 * k)
                    gh = 0
                    for g in key:
                        c = self.charge(g)
                        for i in range(2 * k):
                            ch[i] += c[i]
                        gh += self.gh(g)
                    out.setdefault((tuple(ch), gh), set()).add(key)

        def rec(start, cur, rem, bsum):
            if all(r == 0 for r in rem):
                finish(cur, bsum)
                return
            for idx in range(start, len(gens)):
                wi = W[idx]
                if all(x <= y for x, y in zip(wi, rem)):
                    nrem = tuple(y - x for x, y in zip(wi, rem))
                    cur.append(gens[idx])
                    rec(idx + 1 if Pp[idx] else idx, cur, nrem, bsum + B[idx])
                    cur.pop()
        rec(0, [], lam, 0)
        return {key: sorted(v) for key, v in out.items()}

    def b_ops(self):
        """List of (beta, charge shift, gh shift) for the off-diagonal even and all odd b-generators."""
        k = self.k
        ops = []
        for x in range(k):
            for y in range(k):
                if x != y:
                    for t in (P, MM):
                        beta = (t, (), x, y)
                        ops.append((beta, self.charge(beta), 0))
                beta = (R, (), x, y)
                ops.append((beta, self.charge(beta), -1))
        return ops


def homology(bg, lam, N, field=rank_mod.P1):
    k = bg.k
    Tch = tuple([N] * k + [N] * k)
    ch = bg.chains(lam, (k * N, k * N + 1))
    ghs = sorted(set(g for (c, g) in ch if c == Tch))
    if not ghs:
        return {}
    ops = bg.b_ops()
    allg = sorted(set(g for (c, g) in ch))
    S = {}
    for g in set(ghs) | set(x + 1 for x in ghs):
        sp = rank_mod.make_span(field)
        for beta, cb, dg in ops:
            src_c = tuple(t - c for t, c in zip(Tch, cb))
            src_g = g - dg          # rho(R0) lowers gh by 1: source at g + 1
            for m in ch.get((src_c, src_g), []):
                v = bg.rho(beta, m)
                if v:
                    sp.add(v)
        S[g] = sp
    out = {}
    rankd = {}
    for g in ghs + [ghs[0] - 1]:
        base = rank_mod.make_span(field)
        if g + 1 in S:
            for piv in S[g + 1].piv.values():
                base.add(piv)
        r0 = base.dim()
        for m in ch.get((Tch, g), []):
            v = bg.d(m, N)
            if v:
                base.add(v)
        rankd[g] = base.dim() - r0
    for g in ghs:
        dimC = len(ch.get((Tch, g), []))
        dimS = S[g].dim()
        h = dimC - dimS - rankd.get(g, 0) - rankd.get(g - 1, 0)
        out[g] = (dimC, dimS, rankd.get(g, 0), h)
    return out
