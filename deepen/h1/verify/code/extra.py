"""Explicit functionals: S_free, I, G_eps components, T_lam; eta_v I, eta_v S_free; exact over Q and mod P1."""
import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from vjet import Jets, Echelon, P1, add, popc, shuffle_sign, eps_vT
from vcohom import Slice, Cache, e_, addt

Z5 = (0, 0, 0, 0, 0)


def kappa(alg):
    if alg == "sl2":
        return {(0, 1): 1, (1, 0): 1, (2, 2): 2}
    if alg == "gl2":
        return {(0, 1): 1, (1, 0): 1, (2, 2): 2, (3, 3): 2}
    if alg == "ab1":
        return {(0, 0): 1}
    if alg == "ab2":
        return {(0, 0): 1, (1, 1): 1}


def rmul(J, X, Y):
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


def radd(X, Y, c=1):
    Z = dict(X)
    for k, v in Y.items():
        add(Z, k, c * v)
    return Z


def A_(J, a, al=Z5):
    return {((J.g(a, T, al),), T): 1 for T in range(8)}


def dA_(J, a):
    """jet at 0 of dA^a = sum_v dx_v d_v A^a :  d(c y) = (-1)^{|c|} c dy."""
    Z = {}
    for T in range(8):
        for v in range(3):
            if (T >> v) & 1:
                continue
            al = [0] * 5
            al[v] = 1
            s = (-1 if (1 - popc(T)) % 2 else 1) * eps_vT(v, T)
            add(Z, ((J.g(a, T, tuple(al)),), T | (1 << v)), s)
    return Z


def comp(X, T):
    return {m: c for (m, TT), c in X.items() if TT == T}


def cls(J, t, g, polys, charge0=False, p=None):
    """for polynomials in the slice V_t^g: closed?, rank modulo B^g, kernel description."""
    cache = Cache(J)
    S = Slice(J, t, charge0, cache)
    EI1 = Echelon(p)
    for r in S.Irows.get(g + 1, []):
        EI1.add(r)
    closed = []
    for P in polys:
        q = S.vec(g + 1, J.Qpoly(P)) if J.Qpoly(P) else {}
        closed.append(EI1.reduce(q)[1] is None)
    EB = Echelon(p)
    for r in S.Irows.get(g, []):
        EB.add(r)
    for m in S.M.get(g - 1, []):
        EB.add(S.vec(g, J.Qmono(m)))
    r0 = EB.rank()
    ranks = []
    for P in polys:
        EB.add(S.vec(g, P))
        ranks.append(EB.rank() - r0)
    return closed, ranks, S


def in_DM(J, t, g, P, p=None):
    S = Slice(J, t, False, Cache(J))
    E = Echelon(p)
    for r in S.Irows.get(g, []):
        E.add(r)
    if not P:
        return True
    return E.reduce(S.vec(g, P))[1] is None


U = 7
for alg in ["ab1", "sl2", "gl2"]:
    J = Jets(alg)
    k = kappa(alg)
    dim = J.dim
    A = [A_(J, a) for a in range(dim)]
    dA = [dA_(J, a) for a in range(dim)]
    Sf = {}
    for (a, b), kv in k.items():
        Sf = radd(Sf, rmul(J, A[a], dA[b]), 3 * kv)      # 6 * (1/2) kappa(A, dA)
    Sfree6 = comp(Sf, U)
    Iv = {}
    for (a, d), kv in k.items():
        for (b, c), fd in J.f.items():
            if fd.get(d):
                Iv = radd(Iv, rmul(J, rmul(J, A[a], A[b]), A[c]), kv * fd[d])   # 6 * (1/6) kappa(A,[A,A])
    I6 = comp(Iv, U)
    t = (1, 1, 1, 0, 0)
    for p in (None, P1):
        closed, ranks, S = cls(J, t, 0, [Sfree6, I6], p=p)
        print(f"{alg} field={'Q' if p is None else p}: 6S_free, 6I closed in V: {closed}; cumulative rank mod B^0(V): {ranks}")
    # find x,y with x*6S_free + y*6I in B^0
    cache = Cache(J)
    S = Slice(J, t, False, cache)
    EB = Echelon(None)
    for r in S.Irows.get(0, []):
        EB.add(r)
    for m in S.M.get(-1, []):
        EB.add(S.vec(0, J.Qmono(m)))
    a1 = EB.nf(S.vec(0, Sfree6))
    a2 = EB.nf(S.vec(0, I6)) if I6 else {}
    # proportionality of normal forms
    ratio = None
    if a2:
        ks = set(a1) | set(a2)
        rs = set()
        for kk in ks:
            x, y = a1.get(kk, 0), a2.get(kk, 0)
            rs.add((x / y) if y else ("inf" if x else None))
        ratio = rs
    print(f"{alg}: NF(6S_free)/NF(6I) over all coordinates = {ratio}   (=> (6S_free) - r (6I) in B^0)")
    for v in range(3):
        tv = addt(t, e_(v))
        eI = J.etapoly(I6, v)
        eS = J.etapoly(Sfree6, v)
        print(f"{alg} v={v}: eta_v(6I) in D M: {in_DM(J, tv, -1, eI)};  eta_v(6S_free) in D M: {in_DM(J, tv, -1, eS)}")
    # G_eps components: [kappa(d_zj A, A)]_{T}, T = {tau, 2b} for j = z2 (u=4), T = {tau, 1b} for j = z1 (u=3)
    for (u, T, tt) in [(4, 1 | 4, (1, 0, 1, 0, 1)), (3, 1 | 2, (1, 1, 0, 1, 0))]:
        G = {}
        for (a, b), kv in k.items():
            al = [0] * 5
            al[u] = 1
            G = radd(G, rmul(J, A_(J, a, tuple(al)), A[b]), kv)
        Gc = comp(G, T)
        for p in (None, P1):
            closed, ranks, S = cls(J, tt, 0, [Gc], p=p)
            print(f"{alg} G-component d_z{u - 2} T={T} t={tt} ghost 0: closed {closed}, rank mod B^0 {ranks}  (field {'Q' if p is None else p})")
    if alg == "gl2":
        Tl = comp(A[3], 1)   # lam = u^*, [lam(A)]_{dtau}
        for p in (None, P1):
            closed, ranks, S = cls(J, (1, 0, 0, 0, 0), 0, [Tl], p=p)
            print(f"gl2 T_lam = [u(A)]_dtau t=(1,0,0,0,0) ghost 0: closed {closed}, rank mod B^0 {ranks} (field {'Q' if p is None else p})")
