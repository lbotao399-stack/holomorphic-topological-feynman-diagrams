"""Computational check of Theorems A, B (REPORT.txt) on the jet model.

For theory th, weight w = (w_tau, w_zb1, w_zb2), holomorphic bidegree (a,b), ghost g:
  computed  : dim H^g(V)_{w,(a,b)}                      (V = jets mod total derivatives, full BRST differential)
  predicted : dim H^{g+|S|}(Q)_{(a,b)} if w = -1_{S^c} (w in {0,-1}^3, S = {v : w_v = 0}),  0 otherwise,
              Q = C*_red(g[[z1,z2]]) / (D_z1, D_z2), computed by the same code with holomorphic ghost jets only.
Representatives: for every predicted nonzero entry, a basis of Z^n(Q) is mapped by
  psi_S(phi) = coefficient of dx_S in Phi(A)   (int d^2z dx_{S^c} Phi(A))
and it is checked that psi_S(Z) is closed in V and spans H^g(V) modulo B^g(V).
eta-invariants: V^eta = {x in V : eta_v x = 0 in V, v = tau, zb1, zb2}; dim H^g(V^eta)_{w} vs dim H^{g+3}(Q) (w = 0), 0 (w != 0).
usage: python3 check_h1.py <case>   (case in abelian1, abelian2, sl2a, sl2b, eta)
"""
import sys, time
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import jets as J
import cohom as C
import sla

T0 = time.time()


def log(*a):
    print(*a, flush=True)


_Qcache = {}


def Qdims(th, ab, ns, charge):
    key = (th.name, ab, charge)
    if key not in _Qcache:
        _Qcache[key] = {}
    need = [n for n in ns if n not in _Qcache[key] and n >= 1]
    if need:
        pc = C.Piece(th, None, ab, holo=True, charge=charge)
        res = C.cohomology(pc, need)
        for n in need:
            r = res[n]
            assert r["H"][0] == r["H"][1] == r["H"][2], ("modp/exact disagreement", r)
            _Qcache[key][n] = r["H"][2]
    return {n: (_Qcache[key][n] if n >= 1 else 0) for n in ns}


def S_of(w):
    if all(x in (0, -1) for x in w):
        return sum(b for x, b in zip(w, (1, 2, 4)) if x == 0)
    return None


def sf_poly_of_mono(m):
    """Phi(A): substitute each holomorphic ghost generator (a,0,al) by the superfield sum_T (a,T,al) dx_T."""
    X = {((), 0): 1}
    for g in m:
        a, T, al = g
        assert T == 0
        X = J.R_mul(X, J.superfield(a, al))
    return X


def psi(S, poly):
    out = {}
    for m, c in poly.items():
        X = sf_poly_of_mono(m)
        for mm, cc in J.R_coeff(X, S).items():
            J.padd(out, mm, c * cc)
    return out


def Zbasis_Q(th, ab, n, charge):
    pc = C.Piece(th, None, ab, holo=True, charge=charge)
    Mn = pc.Mg(n)
    imgs = [th.Q({m: 1}) for m in Mn]
    ker = sla.kernel_exact(imgs, pc.DM(n + 1), len(Mn))
    return [{Mn[i]: c for i, c in v.items()} for v in ker]


def check_case(th, w, ab, gs, charge=None, reps=True):
    t = time.time()
    pc = C.Piece(th, w, ab, charge=charge)
    res = C.cohomology(pc, gs)
    S = S_of(w)
    s = bin(S).count("1") if S is not None else None
    pred = {}
    if S is not None:
        qd = Qdims(th, ab, [g + s for g in gs], charge)
        pred = {g: qd[g + s] for g in gs}
    else:
        pred = {g: 0 for g in gs}
    allok = True
    for g in gs:
        r = res[g]
        Hm1, Hm2, He = r["H"]
        ok = (Hm1 == Hm2 == He == pred[g])
        allok &= ok
        line = "  %-8s w=%-12s (a,b)=%s g=%2d  dim M=%6d rkD=%s rkQbar=%s  H(V)=%d [modp %d %d]  pred H^%s(Q)=%d  %s" % (
            th.name, w, ab, g, r["dimM"], r["rkD"][2], r["rkQbar"][2], He, Hm1, Hm2,
            (g + s) if s is not None else "-", pred[g], "OK" if ok else "MISMATCH")
        log(line)
        if reps and S is not None and pred[g] > 0:
            n = g + s
            Z = Zbasis_Q(th, ab, n, charge)
            imgs = [psi(S, z) for z in Z]
            # closedness in V
            Bnext = pc.DM(g + 1)
            rB = sla.rank_exact(Bnext) if Bnext else 0
            closed = all(sla.rank_exact(Bnext + [th.Q(x)]) == rB for x in imgs if th.Q(x))
            Bg = [v for v in pc.QM(g - 1) if v] + pc.DM(g)
            r0 = sla.rank_exact(Bg) if Bg else 0
            r1 = sla.rank_exact(Bg + [x for x in imgs if x])
            ok2 = closed and (r1 - r0 == pred[g])
            allok &= ok2
            log("      psi_S, S=%d: dim Z^%d(Q)=%d, psi(Z) closed in V: %s, rank increase over B^%d(V) = %d (= H): %s" % (
                S, n, len(Z), closed, g, r1 - r0, ok2))
    log("    [%.1fs]" % (time.time() - t))
    return allok


def eta_case(th, w, ab, gs, charge=None, bits=(1, 2, 4)):
    """cohomology of V^eta at weight w (eta_v for v in bits; bits 1 = tau, 2 = zb1, 4 = zb2)."""
    t = time.time()
    pc = C.Piece(th, w, ab, charge=charge)
    ups = []
    for vb, bi in [x for x in ((1, 0), (2, 1), (4, 2)) if x[0] in bits]:
        w2 = list(w); w2[bi] += 1
        ups.append((vb, C.Piece(th, tuple(w2), ab, charge=charge)))
    K = {}
    for g in sorted(set(gs) | {g - 1 for g in gs} | {g + 1 for g in gs}):
        Mg = pc.Mg(g)
        imgs = []
        subs = []
        for vb, pu in ups:
            pass
        for i, m in enumerate(Mg):
            im = {}
            for vb, pu in ups:
                for k, c in J.eta({m: 1}, vb).items():
                    im[(vb, k)] = c
            imgs.append(im)
        for vb, pu in ups:
            for v in pu.DM(g - 1):
                subs.append({(vb, k): c for k, c in v.items()})
        ker = sla.kernel_exact(imgs, subs, len(Mg)) if Mg else []
        K[g] = [{Mg[i]: c for i, c in v.items()} for v in ker]
    allok = True
    S = S_of(w)
    pred = {}
    for g in gs:
        if tuple(bits) == (1, 2, 4):
            pred[g] = Qdims(th, ab, [g + 3], charge)[g + 3] if w == (0, 0, 0) else 0
        elif tuple(bits) == (2, 4):
            if w == (0, 0, 0):
                pred[g] = Qdims(th, ab, [g + 3], charge)[g + 3]
            elif w == (-1, 0, 0):
                pred[g] = Qdims(th, ab, [g + 2], charge)[g + 2]
            else:
                pred[g] = 0
    for g in gs:
        D0 = pc.DM(g)
        rD = sla.rank_exact(D0) if D0 else 0
        dimK = len(K[g]) - rD
        def rkQ(gg):
            src = [th.Q(x) for x in K[gg]]
            src = [x for x in src if x]
            Dn = pc.DM(gg + 1)
            a = sla.rank_exact(Dn) if Dn else 0
            b = sla.rank_exact(Dn + src) if (Dn or src) else 0
            return b - a
        H = dimK - rkQ(g) - rkQ(g - 1)
        ok = H == pred[g]
        extra = ""
        if w == (0, 0, 0) and tuple(bits) == (1, 2, 4):
            n = g + 3
            if n >= 1:
                ph = C.Piece(th, None, ab, holo=True, charge=charge)
                dq = len(ph.Mg(n)) - (sla.rank_exact(ph.DM(n)) if ph.DM(n) else 0)
            else:
                dq = 0
            extra = "  dim Q^%d=%d (cochains; V^eta_0 = psi(Q): %s)" % (n, dq, dq == dimK)
            ok = ok and (dq == dimK)
        allok &= ok
        log("  %-8s V^eta%s w=%-12s (a,b)=%s g=%2d  dim V^eta=%5d  H(V^eta)=%d  pred=%d%s  %s" % (
            th.name, "" if tuple(bits) == (1, 2, 4) else "(zb1,zb2 only)", w, ab, g, dimK, H, pred[g], extra,
            "OK" if ok else "MISMATCH"))
    log("    [%.1fs]" % (time.time() - t))
    return allok


ALLW = [(0, 0, 0), (-1, -1, -1), (0, -1, -1), (-1, 0, -1), (-1, -1, 0), (0, 0, -1), (0, -1, 0), (-1, 0, 0)]
POSW = [(1, 0, 0), (0, 1, 0), (1, -1, 0), (0, 1, -1), (1, 1, -1), (-1, 0, 1)]


def main(case):
    ok = True
    if case == "abelian1":
        th = J.abelian(1)
        for ab in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 1), (2, 2)]:
            for w in ALLW:
                ok &= check_case(th, w, ab, list(range(-3, 4)))
        for ab in [(0, 0), (1, 0), (1, 1)]:
            for w in POSW:
                ok &= check_case(th, w, ab, list(range(-3, 3)))
    elif case == "abelian2":
        th = J.abelian(2)
        for ab in [(0, 0), (1, 0), (1, 1)]:
            for w in ALLW:
                ok &= check_case(th, w, ab, list(range(-3, 3)))
        for ab in [(0, 0), (1, 0)]:
            for w in POSW[:3]:
                ok &= check_case(th, w, ab, list(range(-3, 3)))
    elif case == "sl2a":
        th = J.sl2()
        for ab in [(0, 0), (1, 0)]:
            for w in ALLW:
                ok &= check_case(th, w, ab, list(range(-3, 2)), charge=0)
        for w in [(0, 0, 0), (-1, -1, -1), (0, 0, -1), (-1, 0, 0)]:
            ok &= check_case(th, w, (1, 1), list(range(-2, 2)), charge=0)
        for w in [(1, 0, 0), (0, 1, -1)]:
            ok &= check_case(th, w, (0, 0), list(range(-2, 2)), charge=0)
            ok &= check_case(th, w, (1, 0), list(range(-2, 2)), charge=0)
    elif case == "sl2nc":
        # sl2 without the charge-0 restriction, smallest weights
        th = J.sl2()
        for w in [(0, 0, 0), (0, 0, -1), (-1, -1, -1)]:
            ok &= check_case(th, w, (0, 0), list(range(-3, 2)))
            ok &= check_case(th, w, (1, 0), list(range(-3, 2)))
    elif case == "sl2b":
        th = J.sl2()
        ok &= check_case(th, (0, 0, 0), (2, 1), [0, 1], charge=0)
    elif case == "sl2c":
        th = J.sl2()
        ok &= check_case(th, (0, 0, 0), (2, 2), [0], charge=0)
    elif case == "sl2d":
        th = J.sl2()
        ok &= check_case(th, (0, 0, -1), (1, 1), [2, 3], charge=0)
        ok &= check_case(th, (-1, -1, -1), (1, 1), [4, 5], charge=0)
        ok &= check_case(th, (0, 0, 0), (1, 1), [2], charge=0)
    elif case == "gl2":
        th = J.gl2()
        for w in ALLW:
            ok &= check_case(th, w, (0, 0), list(range(-3, 2)), charge=0)
        for w in [(0, 0, 0), (0, 0, -1), (-1, -1, -1)]:
            ok &= check_case(th, w, (1, 0), list(range(-3, 2)), charge=0)
        ok &= check_case(th, (0, 0, 0), (1, 1), [-1, 0], charge=0)
    elif case == "abelian2c":
        th = J.abelian(2)
        for w in ALLW[1:]:
            ok &= check_case(th, w, (2, 2), list(range(-1, 4)))
        for w in [(1, 0, 0), (0, 1, -1)]:
            ok &= check_case(th, w, (2, 1), list(range(-1, 2)))
    elif case == "gl2b":
        th = J.gl2()
        ok &= check_case(th, (0, 0, -1), (1, 1), [0, 1], charge=0)
        ok &= check_case(th, (-1, 0, 0), (1, 1), [0, 1], charge=0)
    elif case == "abelian2b":
        th = J.abelian(2)
        ok &= check_case(th, (0, 0, 0), (2, 2), [0, 1])
        ok &= check_case(th, (0, 0, 0), (2, 1), [-1, 0, 1])
    elif case == "eta":
        th = J.abelian(1)
        for ab in [(0, 0), (1, 0), (1, 1), (2, 1)]:
            ok &= eta_case(th, (0, 0, 0), ab, list(range(-3, 2)))
        for ab in [(0, 0), (1, 0)]:
            ok &= eta_case(th, (1, 0, 0), ab, list(range(-3, 2)))
            ok &= eta_case(th, (0, 1, 1), ab, list(range(-3, 2)))
        th = J.sl2()
        for ab in [(0, 0), (1, 0)]:
            ok &= eta_case(th, (0, 0, 0), ab, list(range(-2, 2)), charge=0)
        ok &= eta_case(th, (1, 0, 0), (0, 0), list(range(-2, 2)), charge=0)
    if case == "etaneg":
        th = J.abelian(1)
        for ab in [(0, 0), (1, 0), (1, 1)]:
            for w in [(-1, 0, 0), (0, -1, 0), (0, 0, -1), (-1, -1, -1)]:
                ok &= eta_case(th, w, ab, list(range(-3, 4)))
    if case == "eta12":
        th = J.abelian(1)
        for ab in [(0, 0), (1, 0), (1, 1)]:
            for w in [(0, 0, 0), (-1, 0, 0), (1, 0, 0)]:
                ok &= eta_case(th, w, ab, list(range(-3, 3)), bits=(2, 4))
        th = J.sl2()
        for ab in [(0, 0), (1, 0)]:
            for w in [(0, 0, 0), (-1, 0, 0)]:
                ok &= eta_case(th, w, ab, list(range(-2, 2)), charge=0, bits=(2, 4))
    log("case %s: all checks passed: %s   total %.1fs" % (case, ok, time.time() - T0))


if __name__ == "__main__":
    main(sys.argv[1])
