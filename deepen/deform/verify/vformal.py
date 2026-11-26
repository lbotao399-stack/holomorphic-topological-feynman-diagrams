"""Formal-N (trace basis, faithful for N >= arity) verification of deform claims d4-d8 (formal parts)."""
import time, sys
from fractions import Fraction as Fr
from vce import *
from vce import _polymul
from vla import Echelon, kernel, rank_of
from vw import W133, W223

t0 = time.time()


def cbasis(k, a, b, sl):
    E = Echelon()
    basis = []
    ps = [p for p in ALLPERMS[k] if not (sl and has_fixed(p))]
    for m in monos(k, a, b):
        for s in ps:
            w = alt({(s, m, 0): Fr(1)}, k)
            if sl:
                w = restrict_sl(w)
            if not w:
                continue
            ok, _ = E.add(w)
            if ok:
                basis.append(w)
    return basis


def pq(k, a, b):
    return bracket_pq(k, a, b)


def J_of(Tlist, sl):
    """[p0p1] * T, T = list of (coef, sigma)"""
    w = from_poly_tau(pq(3, 0, 1), Tlist)
    return restrict_sl(w) if sl else w


T_D = [(1, (1, 2, 0)), (1, (2, 0, 1))]                 # tr(012) + tr(021)
T_S1 = [(1, (0, 2, 1)), (1, (2, 1, 0)), (1, (1, 0, 2))]  # tr0 tr12 + tr1 tr02 + tr2 tr01
T_S0 = [(1, (0, 1, 2))]                                  # tr0 tr1 tr2


def is_alt(w):
    k = arity(w)
    return all(lin((1, act(pi, w)), (-perm_sign(pi), w)) == {} for pi in ALLPERMS[k])


def F_of(W, sl):
    P = polyprod(pq(4, 0, 1), pq(4, 2, 3))
    out = {}
    for mm, c in P.items():
        for (sig, n), v in W.items():
            add(out, (sig, mm, n), Fr(c) * v)
    out = alt(out, 4)
    return restrict_sl(out) if sl else out


def strip(w):
    return {(s, m, 0): c for (s, m, n), c in w.items()}


def pieces(w):
    out = {}
    for (s, m, n), c in w.items():
        out.setdefault(n, {})[(s, m, 0)] = c
    return out


def fmtL(dct):
    """Laurent poly dict npow->Fraction to string"""
    if not dct:
        return "0"
    return " + ".join(f"({v})N^{p}" for p, v in sorted(dct.items()))


for sl in (False, True):
    tag = "sl_N" if sl else "gl_N"
    print(f"================ formal {tag} ================")
    # C^2
    for (a, b) in ((2, 0), (1, 1), (0, 2)):
        print(f"  dim C^2_({a},{b}) = {len(cbasis(2, a, b, sl))}")
    # H^3 at total weight 2
    H3 = {}
    for (a, b) in ((1, 1), (2, 0), (0, 2)):
        C3 = cbasis(3, a, b, sl)
        dC3 = [d(w, sl) for w in C3]
        r = rank_of(dC3)
        ker = kernel(dC3)
        Z = [lin(*[(c, C3[i]) for i, c in rel.items()]) for rel in ker]
        print(f"  ({a},{b}): dim C^3 = {len(C3)}, rank d3 = {r}, dim Z^3 = dim H^3 = {len(C3) - r}")
        H3[(a, b)] = (C3, Z)
    C3, Z3 = H3[(1, 1)]
    Js = {"D": J_of(T_D, sl)}
    if not sl:
        Js["S1"] = J_of(T_S1, sl)
        Js["S0"] = J_of(T_S0, sl)
    EZ = Echelon()
    for z in Z3:
        EZ.add(z)
    for nm, J in Js.items():
        # alt([p1p2]T) = 6 [p1p2] Sym T  and J alternating
        Tl = {"D": T_D, "S1": T_S1, "S0": T_S0}[nm]
        raw = from_poly_tau(pq(3, 0, 1), [(c, s) for c, s in Tl])
        raw = restrict_sl(raw) if sl else raw
        print(f"  {nm}: alternating {is_alt(J)}, closed {d(J, sl) == {}}, in Z^3 {EZ.coords(J) is not None}")
    print(f"  rank span(J's) = {rank_of(list(Js.values()))}  (= dim Z^3_(1,1) = {len(Z3)})")
    # alt([p1p2] tau) = [p1p2] sum_s s.tau for every tau
    okall = True
    for s in ALLPERMS[3]:
        if sl and has_fixed(s):
            continue
        raw = from_poly_tau(pq(3, 0, 1), [(1, s)])
        A = alt(raw, 3)
        symT = {}
        for pi in ALLPERMS[3]:
            # (pi.tau_s)(X) = tau_s(X_pi) ; with momentum part [p0p1] unchanged
            sig2 = [0] * 3
            for cc in range(3):
                sig2[pi[cc]] = pi[s[cc]]
            add(symT, tuple(sig2), Fr(1))
        B = from_poly_tau(pq(3, 0, 1), [(c, sg) for sg, c in symT.items()])
        if sl:
            A = restrict_sl(A); B = restrict_sl(B)
        okall &= (lin((1, A), (-1, B)) == {})
    print(f"  alt([p1p2] tau) = [p1p2] sum_pi pi.tau for all taus: {okall}")
    if not sl:
        # E_f, E_mix
        Ef = {}
        for i in range(3):
            e = [0] * 6; e[2 * i] = 1; e[2 * i + 1] = 1
            for mm, c in subst(tuple(e), tuple(((s, 1),) for s in range(3)), 3).items():
                for ct, sg in [(2, (1, 2, 0)), (-2, (2, 0, 1))]:
                    add(Ef, (sg, mm, 0), Fr(c * ct))
        e1 = [0] * 6; e1[0] = 1; e1[3] = 1
        e2 = [0] * 6; e2[1] = 1; e2[2] = 1
        raw = {}
        for e in (e1, e2):
            for mm, c in subst(tuple(e), tuple(((s, 1),) for s in range(3)), 3).items():
                add(raw, ((0, 2, 1), mm, 0), Fr(c))
        Emix = alt(raw, 3)
        print(f"  E_f alternating {is_alt(Ef)}; rank span(D,S1,S0,E_f,E_mix) = "
              f"{rank_of(list(Js.values()) + [Ef, Emix])}; rank(dE_f, dE_mix) = {rank_of([d(Ef), d(Emix)])}")
    else:
        Ef = {}
        for i in range(3):
            e = [0] * 6; e[2 * i] = 1; e[2 * i + 1] = 1
            for mm, c in subst(tuple(e), tuple(((s, 1),) for s in range(3)), 3).items():
                for ct, sg in [(2, (1, 2, 0)), (-2, (2, 0, 1))]:
                    add(Ef, (sg, mm, 0), Fr(c * ct))
        print(f"  rank span(D,E_f) = {rank_of([Js['D'], Ef])}; d E_f != 0: {d(Ef, True) != {}}")

    # H^4_(2,2)
    C3b = cbasis(3, 2, 2, sl)
    C4 = cbasis(4, 2, 2, sl)
    dC3 = [d(w, sl) for w in C3b]
    dC4 = [d(w, sl) for w in C4]
    rB = rank_of(dC3)
    r4 = rank_of(dC4)
    ker4 = kernel(dC4)
    Z4 = [lin(*[(c, C4[i]) for i, c in rel.items()]) for rel in ker4]
    print(f"  (2,2): dim C^3 = {len(C3b)}, dim C^4 = {len(C4)}, rank d3 = {rB}, rank d4 = {r4}, "
          f"dim Z^4 = {len(Z4)}, dim B^4 = {rB}, dim H^4 = {len(Z4) - rB}")
    print(f"  Z^4 vectors closed exactly: {all(d(z, sl) == {} for z in Z4)};  "
          f"C^5 coordinates touched: {len(set(k for v in dC4 for k in v))}")
    # B^4 inside Z^4 check
    EZ4 = Echelon()
    for z in Z4:
        EZ4.add(z)
    print(f"  B^4 subset Z^4: {all(EZ4.coords(b) is not None for b in dC3)}")

    F133 = F_of(W133(sl), sl)
    F223 = F_of(W223(sl), sl)
    print(f"  F_133 N-powers {sorted(pieces(F133))}, F_223 N-powers {sorted(pieces(F223))}; "
          f"dF_133 = 0 {d(F133, sl) == {}}, dF_223 = 0 {d(F223, sl) == {}}")
    # H^4 representatives: F_223 pieces first, then fill from Z4
    EB = Echelon()
    for i, b in enumerate(dC3):
        EB.add(b, ('B', i))
    reps = []
    P223 = pieces(F223)
    for n in sorted(P223, reverse=True):
        ok, _ = EB.add(P223[n], ('h', f"F223[N^{n}]"))
        if ok:
            reps.append(f"F223[N^{n}]")
    for i, z in enumerate(Z4):
        ok, _ = EB.add(z, ('h', f"z{i}"))
        if ok:
            reps.append(f"z{i}")
    print(f"  H^4 representatives: {reps}  (count {len(reps)})")

    def cls(w):
        """class of a cocycle (Laurent-in-N) -> dict rep -> Laurent dict"""
        out = {}
        for n, pc in pieces(w).items():
            co = EB.coords(pc)
            assert co is not None, "not in Z^4"
            for (t, nm), c in co.items():
                if t == 'h':
                    out.setdefault(nm, {})
                    out[nm][n] = out[nm].get(n, 0) + c
        return {k: {p: v for p, v in x.items() if v != 0} for k, x in out.items() if any(v != 0 for v in x.values())}

    print(f"  [F_223] = {cls(F223)}")
    print(f"  [F_133] = {cls(F133)}")
    print(f"  [F_133 + F_223] = {cls(lin((1, F133), (1, F223)))}")
    print(f"  [-3F_133 + 2F_223] = {cls(lin((-3, F133), (2, F223)))}   [3F_133 + 2F_223] = {cls(lin((3, F133), (2, F223)))}"
          f"   [-3F_133 - 2F_223] = {cls(lin((-3, F133), (-2, F223)))}")
    # bracket table
    names = list(Js)
    tab = {}
    for i, x in enumerate(names):
        for y in names[i:]:
            br = bracket(Js[x], Js[y], sl)
            closed = d(br, sl) == {}
            bd = sorted({bideg(m) for (_, m, _) in br})
            tab[(x, y)] = cls(br) if br else {}
            print(f"  {{{x},{y}}}: closed {closed}, bidegrees {bd}, class {tab[(x, y)]}")
    sys.stdout.flush()
    if not sl:
        # Moyal
        jM = scale(Js["D"], Fr(1, 2))
        P2 = polyprod(pq(3, 0, 1), pq(3, 0, 1))
        iM = scale(from_poly_tau(P2, [(1, (1, 2, 0)), (-1, (2, 0, 1))]), Fr(1, 8))
        print(f"  Moyal: i_M alternating {is_alt(iM)}; d j_M = 0 {d(jM) == {}}; "
              f"(1/2){{j_M,j_M}} == d i_M {lin((Fr(1, 2), bracket(jM, jM)), (-1, d(iM))) == {}}; "
              f"{{j_M,j_M}} != 0 {bracket(jM, jM) != {}}")
        # independent derivation of s_e expansion: alt(tr(012) e^{(e/2)[p0p1]})/3 order by order
        for order in (1, 2):
            Pk = pq(3, 0, 1)
            Pn = {tuple([0] * 6): 1}
            for _ in range(order):
                Pn = _polymul(Pn, Pk)
            coef = Fr(1, 2 ** order * factorial(order))
            raw = scale(from_poly_tau(Pn, [(1, (1, 2, 0))]), coef)
            se = scale(alt(raw, 3), Fr(1, 3))
            ref = jM if order == 1 else iM
            print(f"  s_e order e^{order}: alt(tr(012) e^((e/2)[p0p1]))/3 coefficient == {'j_M' if order == 1 else 'i_M'}: "
                  f"{lin((1, se), (-1, ref)) == {}}")
        s0 = from_poly_tau({tuple([0] * 6): 1}, [(1, (1, 2, 0)), (-1, (2, 0, 1))])
        print(f"  s_e order e^0 == s_0: {lin((1, scale(alt(from_poly_tau({tuple([0] * 6): 1}, [(1, (1, 2, 0))]), 3), Fr(1, 3))), (-1, s0)) == {}}")
        # basis change pi^*D_sl, M1, U
        print("  (basis change checked in vfixed.py at fixed N)")
    else:
        # sl: coboundary claim and G
        DD = bracket(Js["D"], Js["D"], True)
        P2 = polyprod(pq(3, 0, 1), pq(3, 0, 1))
        base = from_poly_tau(P2, [(1, (1, 2, 0)), (-1, (2, 0, 1))])
        eta = lin((Fr(1, 2), base), (Fr(2), {(s, m, n - 2): c for (s, m, n), c in base.items()}))
        F223s = {(s, m, n - 2): c / 12 for (s, m, n), c in F223.items()}
        lhs = scale(DD, Fr(1, 2))
        rhs = lin((1, F223s), (1, d(eta, True)))
        print(f"  (1/2){{D,D}}_sl == F_223/(12N^2) + d eta_N, eta_N = ((N^2+4)/(2N^2))[p1p2]^2 tr(1[23]): "
              f"{lin((1, lhs), (-1, rhs)) == {}}")
        G = restrict_sl(F_of({((1, 0, 3, 2), 0): Fr(1)}, True))
        print(f"  G = alt([p1p2][p3p4] tr(12)tr(34)): dG = 0 {d(G, True) == {}}, [G] = {cls(G)}")
    print(f"  time {time.time() - t0:.1f}s")
    sys.stdout.flush()
