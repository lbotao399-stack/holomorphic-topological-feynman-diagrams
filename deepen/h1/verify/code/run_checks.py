"""Compare H^g(V_t) with the prediction of Theorem A; psi_S representative test; V^eta checks.
usage: python3 run_checks.py <alg> <mode> [charge0] [exact]
  mode 'A'   : Theorem A table on the case list for <alg>
  mode 'eta' : V^eta, V^{eta_1b,eta_2b} checks
"""
import sys, time
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from vjet import Jets, Echelon, P1, P2, add, popc, shuffle_sign
from vcohom import cohom_V, kernel_eta, cohom_Veta, Slice, Cache, e_, sub, addt

alg = sys.argv[1]
mode = sys.argv[2]
charge0 = "charge0" in sys.argv[3:]
exact = "exact" in sys.argv[3:]
J = Jets(alg)
CE = {}


def ce(a, b):
    if (a, b) not in CE:
        out, _ = cohom_V(J, (0, 0, 0, a, b), charge0, exact=exact)
        CE[(a, b)] = {g: o[P1] for g, o in out.items()}
        if exact:
            for g, o in out.items():
                assert o[P1] == o[P2] == o[None], ("CE field mismatch", a, b, g, o)
    return CE[(a, b)]


def pred(t, g):
    n = t[:3]
    if any(x > 1 for x in n):
        return 0
    s = sum(n)
    return ce(t[3], t[4]).get(g + s, 0)


# ---------------------------------------------------------------- psi_S
def superfield(x):
    """generator x = (a, 0, al) -> R-element sum_T c^a_{T,al} dx_T."""
    a, T, al = J.ginfo[x]
    assert T == 0
    return {((J.g(a, TT, al),), TT): 1 for TT in range(8)}


def rmul(X, Y):
    Z = {}
    for (m1, T1), c1 in X.items():
        for (m2, T2), c2 in Y.items():
            if T1 & T2:
                continue
            s = shuffle_sign(T1, T2) * (-1 if (popc(T1) * J.mono_par(m2)) % 2 else 1)
            sg, m = J.norm(m1 + m2)
            if sg:
                add(Z, (m, T1 | T2), s * sg * c1 * c2)
    return Z


def psi(P, S):
    """P: polynomial in holomorphic ghost jets (dict mono->c); returns [Phi(A)]_{dx_S}."""
    out = {}
    for m, c in P.items():
        X = {((), 0): 1}
        for x in m:
            X = rmul(X, superfield(x))
        for (mm, T), cc in X.items():
            if T == S:
                add(out, mm, c * cc)
    return out


def cocycles_C(a, b, n, p):
    """basis (mod p) of Z^n(C_g^0)_{(a,b)}: x with Q x in D_z M."""
    cache = Cache(J)
    S = Slice(J, (0, 0, 0, a, b), charge0, cache)
    E1 = Echelon(p)
    for r in S.Irows.get(n + 1, []):
        E1.add(r)
    L = S.M.get(n, [])
    E = Echelon(p)
    base = 10 ** 8
    kers = []
    for i, m in enumerate(L):
        w = E1.nf(S.Qrow(n, m)) if L else {}
        vec = dict(w)
        vec[base + i] = 1
        r, c = E.reduce(vec)
        if c is not None and c >= base:
            kers.append({k - base: x for k, x in r.items()})
        else:
            E.add(vec)
    return [{L[i]: x for i, x in kv.items()} for kv in kers]


def rep_test(t, g, p=P1):
    """psi_S of a Z^n basis: closed in V_t, rank modulo B^g(V_t)."""
    S = sum(1 << v for v in range(3) if t[v] == 1)
    n = g + popc(S)
    Z = cocycles_C(t[3], t[4], n, p)
    cache = Cache(J)
    Sl = Slice(J, t, charge0, cache)
    EI = Echelon(p)
    for r in Sl.Irows.get(g + 1, []):
        EI.add(r)
    closed = True
    reps = []
    for z in Z:
        y = psi(z, S)
        y = {m: c % p for m, c in y.items() if c % p}
        reps.append(y)
        qy = {}
        for m, c in y.items():
            for mm, cc in J.Qmono(m).items():
                k = Sl.idx[g + 1][mm]
                qy[k] = (qy.get(k, 0) + c * cc) % p
        if EI.reduce(qy)[1] is not None:
            closed = False
    EB = Echelon(p)
    for r in Sl.Irows.get(g, []):
        EB.add(r)
    for m in Sl.M.get(g - 1, []):
        EB.add(Sl.vec(g, J.Qmono(m)))
    r0 = EB.rank()
    for y in reps:
        EB.add(Sl.vec(g, y))
    return closed, EB.rank() - r0, len(Z)


# ---------------------------------------------------------------- case lists
def box(ab_list):
    out = []
    for ab in ab_list:
        for n0 in (1, 0):
            for n1 in (1, 0):
                for n2 in (1, 0):
                    out.append((n0, n1, n2) + ab)
    return out


CASES = {
    "ab1": box([(0, 0), (1, 0), (0, 1), (1, 1), (2, 1), (2, 2), (3, 1)])
    + [n + ab for ab in [(0, 0), (1, 0), (1, 1), (2, 1)] for n in [(2, 0, 0), (0, 2, 0), (0, 0, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2), (2, 2, 0), (2, 2, 2), (3, 0, 0), (3, 1, 1)]],
    "ab2": box([(0, 0), (1, 0), (1, 1), (2, 1)]) + [n + ab for ab in [(0, 0), (1, 0), (1, 1)] for n in [(2, 0, 0), (2, 1, 1), (1, 2, 0)]],
    "aff1": box([(0, 0), (1, 0), (0, 1), (1, 1), (2, 1)]) + [n + ab for ab in [(0, 0), (1, 0)] for n in [(2, 0, 0), (2, 1, 1), (0, 2, 1)]],
    "heis": box([(0, 0), (1, 0)]) + [(2, 1, 1, 0, 0), (2, 0, 0, 1, 0)],
    "sl2": box([(0, 0), (1, 0), (0, 1), (1, 1)]) + [(2, 1, 1, 0, 0), (2, 0, 0, 1, 0), (1, 2, 0, 0, 0), (0, 2, 1, 1, 0), (1, 1, 2, 0, 0), (2, 1, 1, 1, 0)],
    "gl2": box([(0, 0), (1, 0)]) + [(1, 1, 1, 1, 1), (0, 1, 1, 1, 1), (1, 1, 0, 1, 1)] + [(2, 1, 1, 0, 0)],
}
if alg == "sl2" and not charge0:
    CASES["sl2"] = box([(0, 0), (1, 0)])
if alg == "sl2" and "big" in sys.argv[3:]:
    CASES["sl2"] = [(1, 1, 1, 2, 1), (1, 1, 1, 2, 0), (0, 1, 1, 2, 1)]
for a in sys.argv[3:]:
    if a.startswith("cases="):
        CASES[alg] = [tuple(int(x) for x in c.split(",")) for c in a[6:].split(";")]
GHOSTS = None
if "ghosts" in " ".join(sys.argv[3:]):
    for a in sys.argv[3:]:
        if a.startswith("ghosts="):
            GHOSTS = set(int(x) for x in a[7:].split(","))

if mode == "A":
    tot = 0
    bad = 0
    for t in CASES[alg]:
        t0 = time.time()
        out, S = cohom_V(J, t, charge0, exact=exact, ghosts=GHOSTS)
        for g in sorted(out):
            o = out[g]
            if o[P1] is None:
                continue
            pr = pred(t, g)
            fields = [o[P1], o[P2]] + ([o[None]] if exact else [])
            st = "OK" if all(x == pr for x in fields) else "MISMATCH"
            tot += 1
            bad += st != "OK"
            line = f"{alg} t={t} w={tuple(x - 1 for x in t[:3])} (a,b)={t[3:]} g={g} dimM={o['dimM']} dimV={o.get('dimV')} H(p1,p2{',Q' if exact else ''})={fields} pred={pr} {st}"
            if pr and st == "OK":
                cl, rk, nz = rep_test(t, g)
                line += f" | psi_S rep: closed={cl} rank mod B={rk} (dim Z^n={nz}) {'REP-OK' if cl and rk == pr else 'REP-FAIL'}"
            print(line, flush=True)
        print(f"# {t} done in {time.time() - t0:.1f}s", flush=True)
    print(f"TOTAL {alg} charge0={charge0} exact={exact}: compared {tot} mismatches {bad}")

if mode == "eta":
    tot = bad = 0
    tl = [(1, 1, 1, 0, 0), (1, 1, 1, 1, 0), (1, 1, 1, 1, 1), (1, 1, 1, 2, 1), (2, 1, 1, 0, 0), (1, 2, 1, 0, 0), (2, 2, 2, 0, 0), (2, 1, 1, 1, 0),
          (0, 1, 1, 0, 0), (1, 0, 1, 1, 0), (0, 0, 0, 1, 1), (1, 1, 0, 0, 0), (0, 1, 1, 1, 0)]
    if alg in ("sl2", "gl2"):
        tl = [(1, 1, 1, 0, 0), (1, 1, 1, 1, 0), (2, 1, 1, 0, 0), (0, 1, 1, 0, 0), (1, 1, 0, 1, 0), (0, 1, 1, 1, 0)]
    for vs, name in [((0, 1, 2), "Veta"), ((1, 2), "Veta12")]:
        for t in tl:
            t0 = time.time()
            res, K, S = cohom_Veta(J, t, vs, charge0, p=(P2 if 'p2' in sys.argv[3:] else P1))
            for g in sorted(res):
                r = res[g]
                n = t[:3]
                if vs == (0, 1, 2):
                    pr = ce(t[3], t[4]).get(g + 3, 0) if n == (1, 1, 1) else 0
                else:
                    if n == (1, 1, 1):
                        pr = ce(t[3], t[4]).get(g + 3, 0)
                    elif n == (0, 1, 1):
                        pr = ce(t[3], t[4]).get(g + 2, 0)
                    else:
                        pr = 0
                extra = ""
                if vs == (0, 1, 2) and n == (1, 1, 1):
                    extra = f" dimC^{g + 3}_0={len(Slice(J, (0, 0, 0, t[3], t[4]), charge0, Cache(J)).M.get(g + 3, [])) - 0}"
                if vs == (0, 1, 2):
                    # B1: V^eta = eta_U V  (dimension of eta_U image mod I)
                    tm = sub(t, (1, 1, 1, 0, 0))
                    if min(tm) >= 0:
                        Mm = J.slice(tm, charge0).get(g + 3, [])
                        E = Echelon(P1)
                        for rr in S.Irows.get(g, []):
                            E.add(rr)
                        r0 = E.rank()
                        for m in Mm:
                            P = {m: 1}
                            for v in (2, 1, 0):
                                P = J.etapoly(P, v)
                            if P:
                                E.add(S.vec(g, P))
                        extra += f" dim(eta_U V)={E.rank() - r0}"
                    else:
                        extra += " dim(eta_U V)=0"
                st = "OK" if r["H"] == pr else "MISMATCH"
                tot += 1
                bad += st != "OK"
                print(f"{alg} {name} t={t} w={tuple(x - 1 for x in t[:3])} g={g} dimVeta={r['dimVeta']} I_in_K={r['I_in_K']} H={r['H']} pred={pr} {st}{extra}", flush=True)
            print(f"# {name} {t} {time.time() - t0:.1f}s", flush=True)
    print(f"TOTAL {alg} eta charge0={charge0} prime={P2 if 'p2' in sys.argv[3:] else P1}: compared {tot} mismatches {bad}")
