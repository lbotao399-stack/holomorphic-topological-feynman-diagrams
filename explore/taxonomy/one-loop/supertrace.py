"""Loop supertrace S_H(x) of every T = 0 twist and comparison with the physical gauge anomaly.

S_H(x) = sum_{T*[-1] pairs (V on Omega^{0,*}, V^vee[H-1] on Omega^{H,*})} Str_V rho_V(x)^{H+1}
       + 1/2 sum_{self-paired V} Str_V rho_V(x)^{H+1},        parity of V = statistics of its Omega^{0,0} part.

Checks (exact, random rational x):
 (1) every pair evaluated from both sides: Str_V(rho^k) = Str_{V^vee[H-1]}((-rho^T)^k);
 (2) self-paired V: invariant form omega (rho^T omega + omega rho = 0) of the stated symmetry type;
 (3) S_H(x) = -I_phys(x),  I_phys = 1/2 sum_{fermion fields of ESW, Pi S_chi (x) E} chi Tr_E x^{H+1};
 (4) formal expression of S_H in Tr_adj, Tr_R, Tr_F, Tr_U.
Field contents: ESW = arXiv:2002.10517 Part II, local files 09_logs/taxonomy_notion_import/generated/*.md.
"""
from fractions import Fraction as F
import random

from common import Log
from lie import (basis_sl, adjoint, sym2, dual, block_diag, mmul, madd, transpose, mscale, tr_pow,
                 killing_form, mat)

log = Log("supertrace.log")
rng = random.Random(7)


def rand_elem(B):
    X = mat(len(B[0]))
    for b in B:
        X = madd(X, b, F(rng.randint(-9, 9), rng.randint(1, 5)))
    return X


def eps_blocks(k):
    """symplectic form on (C^2)^k"""
    W = mat(2 * k)
    for a in range(k):
        W[2 * a][2 * a + 1] = F(1)
        W[2 * a + 1][2 * a] = F(-1)
    return W


def setup_sl3():
    B = basis_sl(3)
    x = rand_elem(B)
    ad = adjoint(B, x)
    R = block_diag(x, x, sym2(x))                   # R = 3 + 3 + 6
    Fm = block_diag(dual(x), dual(x), x)            # Fermi F = 3b + 3b + 3
    U = block_diag(R, dual(R))                      # U = T*R
    nR = len(R)
    wU = mat(2 * nR)
    for i in range(nR):
        wU[i][nR + i] = F(1)
        wU[nR + i][i] = F(-1)
    kil = killing_form(B)
    return dict(name="sl(3), R = 3+3+6, F = 3b+3b+3, U = T*R", adj=ad, coadj=dual(ad), R=R, Rd=dual(R),
                F=Fm, Fd=dual(Fm), U=U, omegaU=wU, kil=kil, g=x)


def setup_sl2():
    B = basis_sl(2)
    x = rand_elem(B)
    ad = adjoint(B, x)
    U = block_diag(x, x, x)                         # three half-hypers in the pseudo-real 2
    R = block_diag(x, x)                            # for BF forms: R = 2+2 (U_R = T*R)
    Fm = block_diag(x)
    kil = killing_form(B)
    return dict(name="sl(2), U = 2+2+2 (pseudo-real), R = 2+2, F = 2", adj=ad, coadj=dual(ad), R=R, Rd=dual(R),
                F=Fm, Fd=dual(Fm), U=U, omegaU=eps_blocks(3), kil=kil, g=x)


def Str(M, parity, k, mult=1):
    return mult * (-1) ** parity * tr_pow(M, k)


# twisted theories: ('pair', label, rep, parity, mult, src) | ('self', label, rep, parity, mult, form, symmetry, src)
# physical fermions: (label, chirality, rep, mult, src)
SRC = "09_logs/taxonomy_notion_import/generated/"
THEORIES = [
    dict(name="2d N=(2,2), chirals R", H=1, esw="Thm 12.1: T*[-1]Map(C, T[1](R/g))",
         src=SRC + "12_dimension_2.md l.39 (Thm 12.1), field list l.56-60",
         twist=[("pair", "c+A01 | B", "adj", 1, 1), ("pair", "T[1]: a+chi01 (g[2])", "adj", 0, 1),
                ("pair", "gamma0+gamma01 (R)", "R", 0, 1), ("pair", "T[1]: psi~0+.. (R[1])", "R", 1, 1)],
         phys=[("lambda_+ (S+ x W+ x g)", +1, "adj", 2), ("lambda_- (S- x W- x g)", -1, "adj", 2),
               ("psi (S+ x R)", +1, "R", 1), ("psi (S- x R)", -1, "R", 1),
               ("psi (S+ x R*)", +1, "Rd", 1), ("psi (S- x R*)", -1, "Rd", 1)],
         phys_src=SRC + "12_dimension_2.md l.20-21"),
    dict(name="2d N=(0,2) [ESW (2,0)], chirals R, Fermis F", H=1,
         esw="Thm 12.16: T*[-1]Map(C, R/g); Fermi multiplet = bc system (not in ESW)",
         src=SRC + "12_dimension_2.md l.233-244",
         twist=[("pair", "c+A01 | B", "adj", 1, 1), ("pair", "gamma0+gamma01 | beta", "R", 0, 1),
                ("pair", "Fermi eta | etabar [drv]", "F", 1, 1)],
         phys=[("lambda (S+ x W+ x g), dim W+ = 2", +1, "adj", 2), ("psi (S- x R*)", -1, "Rd", 1),
               ("psibar (S- x R)", -1, "R", 1), ("Fermi eta (S+ x F) [drv]", +1, "F", 1),
               ("Fermi etabar (S+ x F*) [drv]", +1, "Fd", 1)],
         phys_src=SRC + "12_dimension_2.md l.193-197"),
    dict(name="2d N=(4,0), hypers U", H=1, esw="Thm 12.17: T*[-1]Sect(C, (U x K^1/2)//g)",
         src=SRC + "12_dimension_2.md l.291-302",
         twist=[("pair", "c+A01 | A*", "adj", 1, 1), ("pair", "lambda+ | lambda-_10 (g^vee[-1])", "coadj", 1, 1),
                ("pair", "gamma+psi+* | psi_+ + gamma*", "U", 0, 1)],
         phys=[("lambda (S+ x (W~+W~) x g)", +1, "adj", 4), ("psi_pm (S- x (U+U))", -1, "U", 2)],
         phys_src=SRC + "12_dimension_2.md l.252-255"),
    dict(name="2d N=(4,4), hypers U", H=1, esw="Thm 12.7: T*[-1]Sect(C, T[1]((U x K^1/2)//g))",
         src=SRC + "12_dimension_2.md l.133",
         twist=[("pair", "c+A01", "adj", 1, 1), ("pair", "T[1] of g[1]", "adj", 0, 1),
                ("pair", "mu: g^vee[-1]", "coadj", 1, 1), ("pair", "T[1] of g^vee[-1]", "coadj", 0, 1),
                ("pair", "U", "U", 0, 1), ("pair", "T[1] of U", "U", 1, 1)],
         phys=[("(4,0) gauginos", +1, "adj", 4), ("(0,4) gauginos", -1, "adj", 4),
               ("(4,0) hyperini", -1, "U", 2), ("(0,4) hyperini", +1, "U", 2)],
         phys_src="non-chiral (reduction of 3d N=4) [drv]"),
    dict(name="2d N=(8,8)", H=1, esw="Thm (12.3.1): T*[-1]Map(C_Dol, T[1] g/g)",
         src=SRC + "12_dimension_2.md l.177",
         twist=[("pair", "eps^0: g[1]", "adj", 1, 1), ("pair", "eps^0: g (adjoint scalar of g/g)", "adj", 0, 1),
                ("pair", "eps^0 T[1]: g[2]", "adj", 0, 1), ("pair", "eps^0 T[1]: g[1]", "adj", 1, 1),
                ("pair", "dz (C_Dol): g[0]", "adj", 0, 1), ("pair", "dz: g[-1]", "adj", 1, 1),
                ("pair", "dz T[1]: g[1]", "adj", 1, 1), ("pair", "dz T[1]: g[0]", "adj", 0, 1)],
         phys=[("S+ x W+ x g, dim W+ = 8", +1, "adj", 8), ("S- x W- x g, dim W- = 8", -1, "adj", 8)],
         phys_src="non-chiral [drv]"),
    dict(name="2d N=(N+,0), N+ = 2..7", H=1, esw="Thm 12.18: T*[-1]Sect(C, (g^{N+-2} x K^1/2)/g)",
         src=SRC + "12_dimension_2.md l.333-341",
         family="Nplus",
         twist=[("pair", "A0+A01 | A*", "adj", 1, 1),
                ("self", "lambda~ in Omega^{0,*}(K^1/2) x g x C^{N+-2}", "adj", 1, "N-2", "kil", "sym")],
         phys=[("lambda (S+ x W+ x g), dim W+ = N+", +1, "adj", "N")],
         phys_src=SRC + "12_dimension_2.md l.306-309"),
    dict(name="4d N=1, chirals R", H=2, esw="Thm 10.1: T*[-1]Map(C^2, R/g)",
         src=SRC + "10_dimension_4.md l.68-86",
         twist=[("pair", "A0+A01+A02 | B", "adj", 1, 1), ("pair", "gamma0+gamma01+gamma02 | beta", "R", 0, 1)],
         phys=[("lambda_+ (S+ x g)", +1, "adj", 1), ("lambda_- (S- x g)", -1, "adj", 1),
               ("psi_- (S+ x R*)", +1, "Rd", 1), ("psi_+ (S- x R)", -1, "R", 1)],
         phys_src=SRC + "10_dimension_4.md l.15-17"),
    dict(name="4d N=2, hypers U", H=2, esw="Thm 10.2: T*[-1]Sect(C^2, (U x K^1/2[2n])//g)",
         src=SRC + "10_dimension_4.md l.142-171",
         twist=[("pair", "c+A01+B", "adj", 1, 1), ("pair", "a+b+xi* (g[2-4n])", "adj", 0, 1),
                ("pair", "phi+psi+vsigma* (U[2n])", "U", 0, 1)],
         phys=[("lambda_+ (S+ x W x g)", +1, "adj", 2), ("lambda_- (S- x W* x g)", -1, "adj", 2),
               ("psi_- (S+ x U)", +1, "U", 1), ("psi_+ (S- x U)", -1, "U", 1)],
         phys_src=SRC + "10_dimension_4.md l.94-96"),
    dict(name="4d N=4", H=2, esw="Thm 10.9: T*[-1]Map(C^2_Dol, Bg)",
         src=SRC + "10_dimension_4.md l.315-334",
         twist=[("pair", "Omega^0: c+A01+B (g[1])", "adj", 1, 1),
                ("pair", "Omega^1: phi+psi+vs* (C^2 x g)", "adj", 0, 2),
                ("pair", "Omega^2: xi+b*+a* (g[-1]), Map(C^2_Dol,Bg) side", "adj", 1, 1)],
         phys=[("S+ x W x g, dim W = 4", +1, "adj", 4), ("S- x W* x g", -1, "adj", 4)],
         phys_src="reduction of 10d S+ [drv]"),
    dict(name="6d N=(1,0), hypers U, form U//g", H=3, esw="Thm 8.1: Sect(C^3, (U x K^1/2)//g)",
         src=SRC + "8_dimension_6.md l.63-74",
         twist=[("pair", "c+A01+B | xi+..", "adj", 1, 1),
                ("self", "phi+psi in Omega^{0,*}(U x K^1/2), psi^dbar psi", "U", 0, 1, "omegaU", "antisym")],
         phys=[("gauge fermions (S+ x W+ x g), dim W+ = 2", +1, "adj", 2), ("matter fermions (S- x U)", -1, "U", 1)],
         phys_src=SRC + "8_dimension_6.md l.27, l.40-42"),
    dict(name="6d N=(1,0), U = T*R, form BF R/g", H=3, esw="Thm 8.2: T*[-1]Map(C^3, R/g)",
         src=SRC + "8_dimension_6.md l.82-92",
         twist=[("pair", "c+A01+B", "adj", 1, 1), ("pair", "phi_R+psi_R (R) | (R* x K)[2]", "R", 0, 1)],
         phys=[("gauge fermions (S+ x W+ x g)", +1, "adj", 2), ("matter fermions (S- x (R+R*))", -1, "R", 1),
               ("", -1, "Rd", 1)],
         phys_src=SRC + "8_dimension_6.md l.40-42"),
    dict(name="6d N=(1,1)", H=3, esw="Thm 8.3: T*[-1]Map(C^3, g/g)",
         src=SRC + "8_dimension_6.md l.115",
         twist=[("pair", "c+A01+B", "adj", 1, 1), ("pair", "adjoint scalar of g/g (R = g)", "adj", 0, 1)],
         phys=[("S+ x W+ x g", +1, "adj", 2), ("S- x W- x g", -1, "adj", 2)],
         phys_src=SRC + "8_dimension_6.md l.95-100 [(1,1) = (1,0) + U = T*g]"),
    dict(name="8d N=1", H=4, esw="Thm 6.1: T*[-1]Map(C^4, Bg)",
         src=SRC + "6_dimension_8.md l.38-62",
         twist=[("pair", "c+A01+B+C*+a* | a+C+B*+A01*+c*", "adj", 1, 1)],
         phys=[("chi,B,chi~ (Omega^{0,even} = S+)", +1, "adj", 1), ("rho,C (S-)", -1, "adj", 1)],
         phys_src=SRC + "6_dimension_8.md l.36"),
    dict(name="10d N=(1,0)", H=5, esw="Thm 4.4: Map(C^5, Bg) (hCS)",
         src=SRC + "4_dimension_10.md l.109 (Thm 4.4), S_0 l.123",
         twist=[("self", "A = c+A01+B+..", "adj", 1, 1, "kil", "sym")],
         phys=[("rho,B,chi (S+ x g)", +1, "adj", 1)],
         phys_src=SRC + "4_dimension_10.md l.60"),
]


def check_form(M, W, sym):
    """rho^T W + W rho = 0 and W^T = +-W"""
    inv = madd(mmul(transpose(M), W), mmul(W, M))
    ok = all(v == 0 for r in inv for v in r)
    WT = transpose(W)
    if sym == "sym":
        ok &= WT == W
    else:
        ok &= WT == mscale(W, -1)
    return ok


def evaluate(th, S, N=None):
    H = th["H"]
    k = H + 1
    tot = F(0)
    side_ok = True
    form_ok = True
    for ent in th["twist"]:
        kind, label, rep, par, mult = ent[:5]
        if mult == "N-2":
            mult = N - 2
        M = S[rep]
        if kind == "pair":
            a = Str(M, par, k, mult)
            b = mult * (-1) ** ((par + H - 1) % 2) * tr_pow(dual(M), k)
            side_ok &= (a == b)
            tot += a
        else:
            W = S[ent[5]]
            form_ok &= check_form(M, W, ent[6])
            tot += F(1, 2) * Str(M, par, k, mult)
    phys = F(0)
    for (lab, chi, rep, mult) in th["phys"]:
        if mult == "N":
            mult = N
        phys += F(1, 2) * chi * mult * tr_pow(S[rep], k)
    return tot, phys, side_ok, form_ok


SETUPS = [setup_sl3(), setup_sl2()]
for S in SETUPS:
    log("=== algebra:", S["name"])
    # invariance of the Killing form under ad (used for self-paired adjoint fields)
    log("  Killing form ad-invariant and symmetric:", check_form(S["adj"], S["kil"], "sym"),
        "; omega_U invariant and antisymmetric:", check_form(S["U"], S["omegaU"], "antisym"))
    for th in THEORIES:
        Ns = range(2, 8) if th.get("family") == "Nplus" else [None]
        for N in Ns:
            tot, phys, side_ok, form_ok = evaluate(th, S, N)
            k = th["H"] + 1
            tag = th["name"] + ("" if N is None else " [N+=%d]" % N)
            log("  %-48s H=%d  S_H = %-14s I_phys = %-14s S_H = -I_phys: %s   sides: %s  forms: %s"
                % (tag, th["H"], tot, phys, tot == -phys, side_ok, form_ok))

# ---------------------------------------------------------------- formal expressions
log("=== formal S_H (k = H+1; Tr_X := Tr_X x^k; Tr_{R*} = (-1)^k Tr_R; Tr_coadj = (-1)^k Tr_adj)")
SYM = {"adj": "Tr_adj", "coadj": "Tr_adj", "R": "Tr_R", "Rd": "Tr_R", "F": "Tr_F", "Fd": "Tr_F", "U": "Tr_U"}


def formal(th, N=None):
    H = th["H"]
    k = H + 1
    co = {}
    for ent in th["twist"]:
        kind, label, rep, par, mult = ent[:5]
        if mult == "N-2":
            mult = "(N+-2)"
        c = (-1) ** par * (1 if rep not in ("coadj", "Rd", "Fd") else (-1) ** k)
        c = F(c) * (F(1, 2) if kind == "self" else 1)
        key = SYM[rep] + ("" if not isinstance(mult, str) else "*" + mult)
        m = mult if isinstance(mult, int) else 1
        co[key] = co.get(key, 0) + c * m
    return " + ".join("%s %s" % (v, kk) for kk, v in co.items() if v != 0) or "0"


for th in THEORIES:
    log("  %-48s S_%d = %s" % (th["name"], th["H"], formal(th)))
log("  remark: for k odd (H even) Tr_adj x^k = 0 for every g with invariant form; Tr_U x^k = 0 for symplectic U")
log.close()
