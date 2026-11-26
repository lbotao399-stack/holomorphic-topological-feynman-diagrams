"""D_Gamma in bracket notation and the displacement rule X_{P_k} = sigma^P_k z_P + (1 - sigma^P_k) Zbar,
Zbar = sum_P w_P z_P, sigma^P_k = sum_{j>k} x^P_j, z_P = displacement on the spoke of P (oriented from O),
checked symbolically for all Q_2 graphs (5d, 7d, 3d)."""
import itertools
import sympy as sp
import htk_path  # noqa
from theta_kernels import ThetaKernel
from kernels_check import CASES

for name, G, paths, H, zspec in CASES:
    TK = ThetaKernel(G, paths, H, name)
    TK.check_pointwise(npts=2, seed=4)
    zP = [[sp.Symbol('Z%s_%d' % ("ABC"[P], i)) for i in range(H)] for P in range(3)]
    zedge = {TK.paths[P][0]: [TK.signs[P][0] * zP[P][i] for i in range(H)] for P in range(3)}
    X = TK.displacement(zedge)
    ok = True
    for P in range(3):
        for k_, v in enumerate(TK.interior[P]):
            sig = sum(TK.x[P][j] for j in range(k_ + 1, TK.nP[P]))
            for i in range(H):
                Zb = sum(TK.w[Q] * zP[Q][i] for Q in range(3))
                pred = sig * zP[P][i] + (1 - sig) * Zb
                ok &= sp.expand(X[v][i] - pred.subs(TK.w[0], 1 - TK.w[1] - TK.w[2])
                                .subs({TK.x[Q][0]: 1 - sum(TK.x[Q][1:]) for Q in range(3)})
                                - 0) == 0 or sp.expand((X[v][i] - pred).subs(TK.w[0], 1 - TK.w[1] - TK.w[2]).subs({TK.x[Q][0]: 1 - sum(TK.x[Q][1:]) for Q in range(3)})) == 0
    # bracket form of D: Laplace expansion of det[[Lam_A, Lam_B, 0], [0, Lam_B, Lam_C]]
    lab = {}
    for P in range(3):
        for k_, v in enumerate(TK.interior[P]):
            lab[v] = "%s%d" % ("ABC"[P], k_ + 1)
    A, B, C = TK.interior
    terms = []
    D = sp.expand(TK.Dtheta())
    cand = []
    nb = H - len(A)
    for Bp in itertools.combinations(B, nb):
        rest = [b for b in B if b not in Bp]
        left = list(A) + list(Bp)
        right = rest + list(C)
        def br(vs):
            return sp.Matrix([[TK.lam[v][i] for v in vs] for i in range(H)]).det() if H else 1
        cand.append((left, right, sp.expand(br(left) * br(right))))
    # fit coefficients +-1
    cs = sp.symbols('c0:%d' % len(cand))
    eq = sp.expand(D - sum(c * t for c, (_, _, t) in zip(cs, cand)))
    sol = sp.solve(sp.Poly(eq, *[TK.lam[v][i] for v in range(G.nv - 1) for i in range(H)]).coeffs(), cs, dict=True)
    sol = sol[0] if sol else {}
    s_ = " + ".join("(%s)[%s][%s]" % (sol.get(c, c), " ".join(lab[v] for v in l), " ".join(lab[v] for v in r))
                    for c, (l, r, _) in zip(cs, cand))
    print("%s: X rule: %s;  D = %s" % (name, ok, s_))
