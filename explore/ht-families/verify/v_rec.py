"""Own implementation of the recursion of claim a3, written from the statement of a3 only, compared with
Pf(-1/2 dM M^-1 dM) for C_rec (own vcore) and U = det M(C_rec); also a4 (plain merge closed form) and a6.

Elements: dict {(odd bitmask, even exponent tuple): Fraction}.  Odd bits: 0..N0-1 = dt_e of the original
edges, N0 = Theta, N0+1.. = vartheta_1..vartheta_{k-1}.  Even variables: (J, j_1, ..., j_{k-1}).
Insertions: dict {(odd bitmask over dt bits and the Theta bit, (p,)): c} = polynomial in A (power p), linear in Theta.
Rule (a3): A_i = w_i J + j_i, Theta_i = w_i Theta + vartheta_i, j_k = -sum j_a, vartheta_k = -sum vartheta_a,
 I_m = (1 + dt_m J Theta) < int d^{k-1}vartheta prod_i (1 - dR_i A_i Theta_i) prod_{consumed e} I_e(s_e A_i, s_e Theta_i) >_j,
 <j_a j_b> = (w_a delta_ab - w_a w_b)/(2 t_m); int d^{k-1}vartheta = coefficient of vartheta_1..vartheta_{k-1}
 standing at the right end; consumed insertions in creation order, sign of moving them to the right of the
 unconsumed odd ones; rho(C_rec) = U^{-1/2} I_final(0,0), U = prod_m t_m^{k-1} / prod_i w_i.
"""
import sys, random, math, itertools, time
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
import vcore as vc
from vcore import _sign


def mul(X, Y):
    out = {}
    for (a, ea), x in X.items():
        for (b, eb), y in Y.items():
            if a & b:
                continue
            k = (a | b, tuple(p + q for p, q in zip(ea, eb)))
            out[k] = out.get(k, 0) + _sign(a, b) * x * y
    return {k: v for k, v in out.items() if v != 0}


def add(X, Y, c=1):
    out = dict(X)
    for k, v in Y.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def const(c, nev):
    return {(0, (0,) * nev): Fr(c)}


def lift_form(f, nev):
    """form {mask: c} -> element"""
    return {(m, (0,) * nev): Fr(c) for m, c in f.items() if c != 0}


def even_lin(coeffs, nev):
    out = {}
    for i, c in coeffs.items():
        if c != 0:
            e = [0] * nev
            e[i] = 1
            out[(0, tuple(e))] = Fr(c)
    return out


def odd_lin(coeffs, nev):
    return {(1 << b, (0,) * nev): Fr(c) for b, c in coeffs.items() if c != 0}


def power(X, p, nev):
    out = const(1, nev)
    for _ in range(p):
        out = mul(out, X)
    return out


def parity(X):
    ps = set(bin(m).count("1") % 2 for (m, e) in X)
    assert len(ps) <= 1, "inhomogeneous parity"
    return ps.pop() if ps else 0


def wick_moment(alpha, S, memo):
    alpha = tuple(alpha)
    if alpha in memo:
        return memo[alpha]
    if sum(alpha) == 0:
        return Fr(1)
    if sum(alpha) % 2:
        return Fr(0)
    a = next(i for i, x in enumerate(alpha) if x)
    rest = list(alpha)
    rest[a] -= 1
    v = Fr(0)
    for b, x in enumerate(rest):
        if x and S[a][b] != 0:
            r2 = list(rest)
            r2[b] -= 1
            v += x * S[a][b] * wick_moment(r2, S, memo)
    memo[alpha] = v
    return v


def run_recursion(N0, plan, t):
    TH = N0
    tt = {e: Fr(t[e]) for e in range(N0)}
    dt = {e: {1 << e: Fr(1)} for e in range(N0)}
    ins = {e: {(0, (0,)): Fr(1)} for e in range(N0)}
    alive, par = [], {}
    sign = 1
    U = Fr(1)
    created = []
    for mi, (p, q, chains, signs) in enumerate(plan):
        k = len(chains)
        nev = k
        R = [sum(tt[e] for e in ch) for ch in chains]
        dR = []
        for ch in chains:
            f = {}
            for e in ch:
                f = vc.fadd(f, dt[e])
            dR.append(f)
        tm = 1 / sum(1 / r for r in R)
        w = [tm / r for r in R]
        dtm = {}
        for i in range(k):
            dtm = vc.fadd(dtm, dR[i], w[i] ** 2)
        U *= tm ** (k - 1) / math.prod(w)
        # chain variables
        A, TV = [], []
        for i in range(k):
            if i < k - 1:
                A.append(even_lin({0: w[i], i + 1: 1}, nev))
                TV.append(odd_lin({TH: w[i], TH + 1 + i: 1}, nev))
            else:
                A.append(even_lin(dict([(0, w[i])] + [(a + 1, -1) for a in range(k - 1)]), nev))
                TV.append(odd_lin(dict([(TH, w[i])] + [(TH + 1 + a, -1) for a in range(k - 1)]), nev))
        Z = const(1, nev)
        for i in range(k):
            Z = mul(Z, add(const(1, nev), mul(mul(lift_form(dR[i], nev), A[i]), TV[i]), -1))
        where = {}
        for i, (ch, sg) in enumerate(zip(chains, signs)):
            for e, s in zip(ch, sg):
                where[e] = (i, s)
        consumed = [e for e in alive if e in where]
        uncons = [e for e in alive if e not in where]
        for x, e in enumerate(alive):
            if e in where:
                for u in alive[x + 1:]:
                    if u not in where and par[e] and par[u]:
                        sign = -sign
        for e in consumed:
            i, s = where[e]
            sub = {}
            for (m, (pp,)), c in ins[e].items():
                base = {(m & ~(1 << TH), (0,) * nev): c}
                term = mul(base, power(add({}, A[i], s), pp, nev))
                if m >> TH & 1:
                    term = mul(mul(base, add({}, TV[i], s)), power(add({}, A[i], s), pp, nev))
                sub = add(sub, term)
            Z = mul(Z, sub)
        # Berezin over vartheta_1..vartheta_{k-1} at the right end
        vmask = sum(1 << (TH + 1 + a) for a in range(k - 1))
        Zb = {}
        for (m, ev), c in Z.items():
            if m & vmask == vmask:
                rest = m & ~vmask
                # term = c * (gens of rest) (vartheta block); vartheta bits are the largest -> already at the right end
                assert rest < (1 << (TH + 1))
                key = (rest, ev)
                Zb[key] = Zb.get(key, 0) + c
        # Gaussian average over j
        S = [[(w[a] * (a == b) - w[a] * w[b]) / (2 * tm) for b in range(k - 1)] for a in range(k - 1)]
        memo = {}
        Zj = {}
        for (m, ev), c in Zb.items():
            mom = wick_moment(ev[1:], S, memo)
            if mom:
                key = (m, (ev[0],))
                Zj[key] = Zj.get(key, 0) + c * mom
        Zj = {kk: v for kk, v in Zj.items() if v != 0}
        corr = add(const(1, 1), mul(mul(lift_form(dtm, 1), {(0, (1,)): Fr(1)}), {(1 << TH, (0,)): Fr(1)}))
        Inew = mul(corr, Zj)
        en = N0 + mi
        tt[en], dt[en], ins[en] = tm, dtm, Inew
        par[en] = parity(Inew)
        alive = uncons + [en]
        created.append(dict(move=mi, R=R, w=w, tm=tm, dtm=dtm, ins=Inew, chains=chains))
    final = ins[N0 + len(plan) - 1]
    rho = {}
    for (m, (pp,)), c in final.items():
        if pp == 0 and not (m >> TH & 1):
            rho[m] = rho.get(m, 0) + c
    return {kk: sign * v for kk, v in rho.items() if v != 0}, U, created


def crec(N0, plan):
    flow = {e: {e: 1} for e in range(N0)}
    rows_by = []
    for mi, (p, q, chs, sgs) in enumerate(plan):
        def cf(ch, sg):
            f = {}
            for e, s in zip(ch, sg):
                for ee, c in flow[e].items():
                    f[ee] = f.get(ee, 0) + s * c
            return f
        fk = cf(chs[-1], sgs[-1])
        rr = []
        for a in range(len(chs) - 1):
            fa = cf(chs[a], sgs[a])
            rr.append([fa.get(e, 0) - fk.get(e, 0) for e in range(N0)])
        rows_by.append(rr)
        flow[N0 + mi] = fk
    return [r for rr in reversed(rows_by) for r in rr]


def is_poly_in_jhat(I, tm, dtm, TH):
    """I(J, Theta) == sum_p c_p (J + dtm/(2 tm) Theta)^p with c_p the Theta-free coefficients?"""
    pred = {}
    for (m, (pp,)), c in I.items():
        if m >> TH & 1:
            continue
        pred[(m, (pp,))] = pred.get((m, (pp,)), 0) + c
        if pp >= 1:
            # c J^pp -> + c pp J^{pp-1} (dtm/2tm) Theta
            term = mul(mul({(m, (pp - 1,)): c * pp}, lift_form(vc.fscale(dtm, 1 / (2 * tm)), 1)), {(1 << TH, (0,)): Fr(1)})
            for kk, v in term.items():
                pred[kk] = pred.get(kk, 0) + v
    pred = {kk: v for kk, v in pred.items() if v != 0}
    return pred == I


if __name__ == "__main__":
    sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
    import htk
    from htk import spred
    from v_b1 import GR
    from v_b2 import banana, nested
    rng = random.Random(123)
    cases = []
    for k in range(5):
        cases.append(("G%d" % k, 6, GR[k], 2))
    cases.append(("theta(1,2,2)", 4, [(0, 1), (1, 3), (0, 2), (2, 3), (0, 3)], 2))
    for lens in [(1, 3, 3, 3, 3), (2, 2, 3, 3, 3)]:
        nv, e, _ = banana(lens)
        cases.append(("banana%s" % (lens,), nv, e, 3))
    for inner in [(1, 3, 3), (2, 2, 3)]:
        nv, e = nested(inner)
        cases.append(("nested%s" % (inner,), nv, e, 3))
    merges = []
    for name, nv, edges, D in cases:
        G = htk.Graph(nv, edges)
        for prefer in ("large", "small"):
            red = spred.sp_reduce(G, D=D, prefer=prefer)
            plan = [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]
            Cr = crec(len(edges), plan)
            net = vc.Net(nv, edges)
            B = net.Bmat()
            assert all(sum(B[e][v] * r[e] for e in range(len(edges))) == 0 for r in Cr for v in range(nv - 1))
            assert len(Cr) == net.n
            res = []
            t0 = time.time()
            for _ in range(2):
                t = [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(len(edges))]
                rho, U, created = run_recursion(len(edges), plan, t)
                pf, Upf, _ = vc.rho_forms(Cr, t, [{1 << e: Fr(1)} for e in range(len(edges))])
                res.append((rho == pf and len(pf) > 0, U == Upf))
            # a6 classification on the last point
            for cinfo in created:
                chs = cinfo["chains"]
                kinds = []
                for ch in chs:
                    red_e = [e for e in ch if e >= len(edges)]
                    kinds.append("plain" if not red_e else ("single" if len(ch) == 1 else "mixed"))
                merges.append((name, prefer, cinfo["move"], kinds, is_poly_in_jhat(cinfo["ins"], cinfo["tm"], cinfo["dtm"], len(edges))))
            print("A3 %-22s plan %-5s moves %d: own recursion == Pf(C_rec), U == det M(C_rec) at 2 points: %s  (%.1fs)" % (
                name, prefer, len(plan), res, time.time() - t0))
    # a6 summary
    tot = len(merges)
    rule_ok = all(poly == all(kd in ("plain", "single") for kd in kinds) for (_, _, _, kinds, poly) in merges)
    npoly = sum(1 for m in merges if m[4])
    print("A6 merges tested: %d ; polynomial in own Jhat: %d ; rule 'polynomial iff every chain plain or a single edge' holds for all: %s" % (
        tot, npoly, rule_ok))
    for m in merges:
        if m[4] != all(kd in ("plain", "single") for kd in m[3]):
            print("   exception:", m)
    sub = [m for m in merges if m[0].startswith("G")]
    print("A6 restricted to G0..G4 (both plans): %d merges, %d polynomial in Jhat, rule holds: %s" % (
        len(sub), sum(1 for m in sub if m[4]), all(m[4] == all(kd in ("plain", "single") for kd in m[3]) for m in sub)))

    # a4: plain merge of k single edges, closed form
    for k in range(2, 7):
        oks = []
        for _ in range(2):
            N0 = k
            plan = [(0, 1, [[i] for i in range(k)], [[1] for _ in range(k)])]
            t = [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(k)]
            rho, U, created = run_recursion(N0, plan, t)
            c = created[0]
            tm, w, dtm = c["tm"], c["w"], c["dtm"]
            TH = N0
            dw = []
            for i in range(k):
                R = Fr(t[i])
                dw.append(vc.fadd(vc.fscale(dtm, 1 / R), {1 << i: Fr(1)}, -tm / R ** 2))
            pref = {0: Fr((-1) ** ((k - 1) * (k - 2) // 2)) * tm ** (k - 1) / math.prod(w)}
            for i in range(k - 1):
                pref = vc.wedge(pref, dw[i])
            p = k - 1
            herm = {}
            for mm in range(p // 2 + 1):
                cm = Fr((-1) ** mm * math.factorial(p), math.factorial(mm) * math.factorial(p - 2 * mm)) / (4 * tm) ** mm
                q = p - 2 * mm
                herm = add(herm, {(0, (q,)): cm})
                if q >= 1:
                    herm = add(herm, mul(mul({(0, (q - 1,)): cm * q}, lift_form(vc.fscale(dtm, 1 / (2 * tm)), 1)), {(1 << TH, (0,)): Fr(1)}))
            closed = mul(lift_form(pref, 1), herm)
            oks.append(closed == c["ins"])
        print("A4 k=%d plain merge: own recursion insertion == s_k t^{k-1}/prod w dw_1..dw_{k-1} :Jhat^{k-1}:_t (2 points): %s" % (k, oks))
