"""Extra exact checks.
(1) eta_S (pi^* phi) = +- psi_S(phi) as polynomials in M (phi random holomorphic ghost polynomial), all S.
(2) sl2, kappa = tr (kappa(e,f) = 1, kappa(h,h) = 2):  S_free = (1/2) kappa(A, dA), I = (1/6) kappa(A,[A,A]),
    N = (1/2) kappa(A, A) (top-form coefficients).  Q S_free, Q I in D M;  N (vanishes identically); x S_free + y I in B^0(V) solved;
    I not in B^0(V);  eta_v I in D M;  eta_v S_free not in D M.
"""
import sys, random
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import jets as J
import cohom as C
import sla
from check_h1 import psi

random.seed(7)


def eta_S(P, S):
    for vb in (1, 2, 4):
        if S & vb:
            P = J.eta(P, vb)
    return P


def check1(th, nmax=4):
    gens = [g for g in J.gens_within(th.dim, (0, 0, 0, 2, 2), holo_only=True)]
    allok = True
    signs = {}
    for k in range(1, nmax + 1):
        for trial in range(4):
            phi = {}
            for _ in range(3):
                m = {(): 1}
                for g in random.sample(gens, k):
                    m = J.pmul(m, {(g,): 1})
                for mm, c in m.items():
                    J.padd(phi, mm, random.randint(-3, 3) * c)
            if not phi:
                continue
            for S in range(8):
                a = eta_S(dict(phi), S)
                b = psi(S, phi)
                diffp = {}
                for kk, v in a.items(): J.padd(diffp, kk, v)
                for kk, v in b.items(): J.padd(diffp, kk, -v)
                summ = {}
                for kk, v in a.items(): J.padd(summ, kk, v)
                for kk, v in b.items(): J.padd(summ, kk, v)
                if not diffp:
                    s = +1
                elif not summ:
                    s = -1
                else:
                    s = 0
                    allok = False
                signs.setdefault((k, S), set()).add(s)
    print(th.name, "eta_S pi^*phi = sign(k,S) psi_S(phi), arity k = 1..%d, all S:" % nmax, allok)
    print("   signs (k,S):", {key: sorted(v) for key, v in sorted(signs.items())})
    return allok and all(len(v) == 1 for v in signs.values())


def top(X):
    return J.R_coeff(X, 7)


def check2():
    th = J.sl2()
    kap = {(0, 1): 1, (1, 0): 1, (2, 2): 2}
    A = [J.superfield(a) for a in range(3)]
    dA = [J.R_d(A[a]) for a in range(3)]
    # 6 S_free = 3 sum kappa_ab A^a dA^b ; 2N = sum kappa_ab A^a A^b ; 6 I = sum kappa_ad f^d_bc A^a A^b A^c
    def lin(terms):
        out = {}
        for c, X in terms:
            for k, v in X.items():
                nv = out.get(k, 0) + c * v
                if nv: out[k] = nv
                else: out.pop(k, None)
        return out
    S6 = top(lin([(3 * c, J.R_mul(A[a], dA[b])) for (a, b), c in kap.items()]))
    N2 = top(lin([(c, J.R_mul(A[a], A[b])) for (a, b), c in kap.items()]))
    terms = []
    for (a, d), c in kap.items():
        for b in range(3):
            for cc in range(3):
                fd = th.f[b][cc].get(d, 0)
                if fd:
                    terms.append((c * fd, J.R_mul(J.R_mul(A[a], A[b]), A[cc])))
    I6 = top(lin(terms))
    pc = C.Piece(th, (0, 0, 0), (0, 0), charge=0)
    D0, D1, Dm1 = pc.DM(0), pc.DM(1), pc.DM(-1)
    rD1 = sla.rank_exact(D1)
    closedS = sla.rank_exact(D1 + [th.Q(S6)]) == rD1 if th.Q(S6) else True
    closedI = sla.rank_exact(D1 + [th.Q(I6)]) == rD1 if th.Q(I6) else True
    ker = sla.kernel_exact([th.Q(N2), S6, I6], D0, 3)
    B0 = [v for v in pc.QM(-1) if v] + D0
    rB0 = sla.rank_exact(B0)
    I_exact = sla.rank_exact(B0 + [I6]) == rB0
    S_exact = sla.rank_exact(B0 + [S6]) == rB0
    print("sl2: Q(6 S_free) in DM:", closedS, "  Q(6 I) in DM:", closedI)
    print("sl2: kernel of (Q(2N), 6S_free, 6I) mod D M^0:", ker)
    print("sl2: 6I exact in V:", I_exact, "  6S_free exact in V:", S_exact)
    print("sl2: N = (1/2) kappa(A,A) top coefficient vanishes identically:", not N2)
    kerSI = sla.kernel_exact([S6, I6], B0, 2)
    print("sl2: kernel of (6S_free, 6I) modulo B^0(V) = Q M^-1 + D M^0:", kerSI)
    # eta-invariance
    for vb, name in ((1, "tau"), (2, "zb1"), (4, "zb2")):
        bi = {1: 0, 2: 1, 4: 2}[vb]
        w2 = [0, 0, 0]; w2[bi] = 1
        pu = C.Piece(th, tuple(w2), (0, 0), charge=0)
        Du = pu.DM(-1)
        r = sla.rank_exact(Du) if Du else 0
        eI = J.eta(I6, vb); eS = J.eta(S6, vb)
        inI = (not eI) or sla.rank_exact(Du + [eI]) == r
        inS = (not eS) or sla.rank_exact(Du + [eS]) == r
        print("sl2: eta_%s(6I) in DM: %s   eta_%s(6S_free) in DM: %s" % (name, inI, name, inS))


if __name__ == "__main__":
    ok = check1(J.abelian(2)) and check1(J.sl2())
    print("check1 all:", ok)
    check2()
