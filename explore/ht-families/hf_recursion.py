"""Super-Gaussian recursion for rho_n on a series-parallel network (exact, at rational t).

rho_n(M) = pi^{-n/2} int d^n theta d^n a exp(-sum_e [t_e A_e^2 + dt_e A_e Theta_e]),
A = C^T a, Theta = C^T theta  (C = C_rec, see cycle_basis).

Edge factor exp(-t A^2 - dt A Theta) times an insertion I_e(A, Theta) (polynomial in A,
linear in Theta, coefficients forms in the dt of the original edges).  Chain: common (A, Theta).
Parallel merge of k chains (R_i, w_i = t_red / R_i):
   A_i = w_i J + j_i,  Theta_i = w_i Theta + vartheta_i,  sum_i j_i = sum_i vartheta_i = 0,
   I_new(J, Theta) = (1 + dt_red J Theta) * < int d^{k-1} vartheta  prod_i (1 - dR_i A_i Theta_i) P_i(A_i, Theta_i) >_j,
   < >_j: normalized Gaussian, < j_a j_b > = (1/2) (M_loc^{-1})_ab, (M_loc)_ab = R_a delta_ab + R_k  (a,b < k).
   rho_n = U^{-1/2} I_final(0, 0),  U = prod_m det M_loc(m).

Grassmann elements: dict {(odd tuple sorted, even exponent tuple): Fraction}.
Odd generator ids: 0..N0-1 = dt_e of the original edges; THETA = 10**6 (placeholder of an
insertion); VT0 + a = vartheta_a.  Even variables: index 0 = J (or A), 1..k-1 = j_a.
"""
from fractions import Fraction as Fr
import itertools

from hf_brute import merge_sign, fdet, finv

THETA = 10 ** 6
VT0 = 2 * 10 ** 6


def g_mul(X, Y):
    out = {}
    for (oa, ea), ca in X.items():
        for (ob, eb), cb in Y.items():
            s = merge_sign(oa, ob)
            if s == 0:
                continue
            key = (tuple(sorted(oa + ob)), tuple(x + y for x, y in zip(ea, eb)))
            out[key] = out.get(key, 0) + s * ca * cb
    return {k: v for k, v in out.items() if v != 0}


def g_add(X, Y, c=1):
    out = dict(X)
    for k, v in Y.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def g_const(c, nev):
    return {((), (0,) * nev): Fr(c)}


def g_oneform(f, nev):
    """1-form dict {(e,): c} -> element"""
    return {((e,), (0,) * nev): Fr(c) for (e,), c in f.items() if c != 0}


def g_even(lin, nev):
    """linear even element sum_i lin[i] * var_i (lin: dict var index -> coeff)"""
    out = {}
    for i, c in lin.items():
        ex = [0] * nev
        ex[i] = 1
        out[((), tuple(ex))] = Fr(c)
    return out


def g_pow(X, p, nev):
    out = g_const(1, nev)
    for _ in range(p):
        out = g_mul(out, X)
    return out


def substitute(ins, A_elem, TH_elem, nev):
    """ins: insertion {(odd tuple possibly containing THETA, (p,)): c} as a polynomial in (A, Theta).
    Returns ins(A_elem, TH_elem) as an element with nev even variables.  THETA is the largest
    odd id among the form generators, so it sits at the end of the sorted tuple."""
    out = {}
    cache = {}
    for (odd, (p,)), c in ins.items():
        if p not in cache:
            cache[p] = g_pow(A_elem, p, nev)
        Ap = cache[p]
        if odd and odd[-1] == THETA:
            base = {(odd[:-1], (0,) * nev): c}
            term = g_mul(g_mul(base, TH_elem), Ap)
        else:
            term = g_mul({(odd, (0,) * nev): c}, Ap)
        out = g_add(out, term)
    return out


def berezin_right(X, ids):
    """coefficient of prod_{a in ids} vartheta_a (ids sorted, they are the largest generators)
    standing at the right end."""
    ids = tuple(ids)
    out = {}
    m = len(ids)
    for (odd, ev), c in X.items():
        if m == 0:
            out[(odd, ev)] = out.get((odd, ev), 0) + c
            continue
        if odd[-m:] == ids and not (set(odd[:-m]) & set(ids)):
            key = (odd[:-m], ev)
            out[key] = out.get(key, 0) + c
    return {k: v for k, v in out.items() if v != 0}


def gauss_moment(alpha, S, memo):
    """E[prod_a j_a^alpha_a] for a centered Gaussian with covariance S (exact)"""
    alpha = tuple(alpha)
    if alpha in memo:
        return memo[alpha]
    tot = sum(alpha)
    if tot == 0:
        return Fr(1)
    if tot % 2:
        return Fr(0)
    a = next(i for i, x in enumerate(alpha) if x)
    rest = list(alpha)
    rest[a] -= 1
    val = Fr(0)
    for b, x in enumerate(rest):
        if x and S[a][b] != 0:
            r2 = list(rest)
            r2[b] -= 1
            val += S[a][b] * x * gauss_moment(r2, S, memo)
    memo[alpha] = val
    return val


def gauss_average(X, S):
    """X with even vars (J, j_1..j_{k-1}): average over j -> element with even var J only"""
    memo = {}
    out = {}
    for (odd, ev), c in X.items():
        m = gauss_moment(ev[1:], S, memo)
        if m == 0:
            continue
        key = (odd, (ev[0],))
        out[key] = out.get(key, 0) + c * m
    return {k: v for k, v in out.items() if v != 0}


class Recursion:
    """plan: list of (p, q, chains, signs) with chains = lists of edge ids (original 0..N0-1,
    new N0, N0+1, ... in the order of the moves), signs[i][j] = +1 if edge j of chain i is
    traversed tail -> head when going from p to q."""

    def __init__(self, nv, edges, plan):
        self.nv, self.edges, self.plan = nv, list(edges), plan
        self.N0 = len(edges)

    def run(self, t):
        N0 = self.N0
        tt = {e: Fr(t[e]) for e in range(N0)}
        dt = {e: {(e,): Fr(1)} for e in range(N0)}
        ins = {e: {((), (0,)): Fr(1)} for e in range(N0)}   # insertions, even var = A
        U = Fr(1)
        self.local = []
        self.sign = 1
        alive, par = [], {}
        for mi, (p, q, chains, signs) in enumerate(self.plan):
            k = len(chains)
            R = [sum(tt[e] for e in ch) for ch in chains]
            dR = []
            for ch in chains:
                f = {}
                for e in ch:
                    for key, c in dt[e].items():
                        f[key] = f.get(key, 0) + c
                dR.append({kk: v for kk, v in f.items() if v != 0})
            tred = 1 / sum(1 / r for r in R)
            w = [tred / r for r in R]
            dtred = {}
            for i in range(k):
                for key, c in dR[i].items():
                    dtred[key] = dtred.get(key, 0) + w[i] ** 2 * c
            dtred = {kk: v for kk, v in dtred.items() if v != 0}
            nev = k                       # J, j_1..j_{k-1}
            # consumed insertions, multiplied in creation order (alive list), with the sign of
            # moving them to the right of the unconsumed ones
            where = {}
            for i, (ch, sg) in enumerate(zip(chains, signs)):
                for e, s in zip(ch, sg):
                    where[e] = (i, s)
            consumed = [e for e in alive if e in where]
            uncons = [e for e in alive if e not in where]
            sgn = 1
            for x, e in enumerate(alive):
                if e in where:
                    for u in alive[x + 1:]:
                        if u not in where and par[e] and par[u]:
                            sgn = -sgn
            self.sign *= sgn
            Z = g_const(1, nev)
            Avars, THvars = [], []
            for i in range(k):
                if i < k - 1:
                    Ai = g_add(g_even({0: w[i]}, nev), g_even({i + 1: 1}, nev))
                    THi = g_add({((THETA,), (0,) * nev): w[i]}, {((VT0 + i,), (0,) * nev): Fr(1)})
                else:
                    Ai = g_add(g_even({0: w[i]}, nev), g_even({a + 1: -1 for a in range(k - 1)}, nev))
                    THi = {((THETA,), (0,) * nev): w[i]}
                    for a in range(k - 1):
                        THi[((VT0 + a,), (0,) * nev)] = Fr(-1)
                Avars.append(Ai)
                THvars.append(THi)
                odd_part = g_add(g_const(1, nev), g_mul(g_mul(g_oneform(dR[i], nev), Ai), THi), -1)
                Z = g_mul(Z, odd_part)
            for e in consumed:
                i, s = where[e]
                Ie = ins[e]
                if s < 0:   # (A, Theta) -> (-A, -Theta)
                    Ie = {(odd, (pp,)): c * (-1) ** pp * (-1 if (odd and odd[-1] == THETA) else 1)
                          for (odd, (pp,)), c in Ie.items()}
                Z = g_mul(Z, substitute(Ie, Avars[i], THvars[i], nev))
            Z = berezin_right(Z, [VT0 + a for a in range(k - 1)])
            Mloc = [[(R[a] if a == b else 0) + R[k - 1] for b in range(k - 1)] for a in range(k - 1)]
            S = [[x / 2 for x in row] for row in finv(Mloc)] if k > 1 else []
            U *= fdet(Mloc) if k > 1 else 1
            Z = gauss_average(Z, S) if k > 1 else {(odd, (ev[0],)): c for (odd, ev), c in Z.items()}
            # multiply by (1 + dt_red J Theta)
            corr = g_add({((), (0,)): Fr(1)}, g_mul(g_mul(g_oneform(dtred, 1), {((), (1,)): Fr(1)}),
                                                     {((THETA,), (0,)): Fr(1)}))
            Inew = g_mul(corr, Z)
            enew = N0 + mi
            tt[enew], dt[enew], ins[enew] = tred, dtred, Inew
            alive = uncons + [enew]
            par[enew] = (sum(par[e] for e in consumed) + k - 1) % 2
            self.local.append(dict(R=R, w=w, tred=tred))
        self.U = U
        self.ins, self.tt, self.dt = ins, tt, dt
        final = ins[N0 + len(self.plan) - 1]
        rho = {}
        for (odd, (pp,)), c in final.items():
            if pp == 0 and not (odd and odd[-1] == THETA):
                rho[odd] = rho.get(odd, 0) + c
        return {kk: self.sign * v for kk, v in rho.items() if v != 0}, U

    @staticmethod
    def _mul_ins(X, Y):
        return g_mul(X, Y)

    def cycle_basis(self):
        """Integer cycle basis C_rec implied by the recursion: for move m, loop a (a < k):
        +chain a - chain k; the unit flow p -> q of a new edge is routed through its chain k.
        Row order: moves in reverse order, a ascending within a move."""
        N0 = self.N0
        flow = {e: {e: 1} for e in range(N0)}       # unit flow tail->head of each edge id

        def chain_flow(ch, sg):
            f = {}
            for e, s in zip(ch, sg):
                for ee, c in flow[e].items():
                    f[ee] = f.get(ee, 0) + s * c
            return f

        rows_by_move = []
        for mi, (p, q, chains, signs) in enumerate(self.plan):
            k = len(chains)
            fk = chain_flow(chains[-1], signs[-1])
            rows = []
            for a in range(k - 1):
                fa = chain_flow(chains[a], signs[a])
                row = [0] * N0
                for ee, c in fa.items():
                    row[ee] += c
                for ee, c in fk.items():
                    row[ee] -= c
                rows.append(row)
            rows_by_move.append(rows)
            flow[N0 + mi] = fk
        C = []
        for rows in reversed(rows_by_move):
            C.extend(rows)
        return C
