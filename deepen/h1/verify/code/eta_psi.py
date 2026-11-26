"""eta_S phi = +- psi_S(phi) for monomials phi in holomorphic ghost jets (exact integers); sign table by (arity, S)."""
import sys, random, itertools
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from vjet import Jets, add, popc, shuffle_sign
random.seed(3)

for alg in ["ab2", "sl2", "aff1"]:
    J = Jets(alg)
    def superfield(x):
        a, T, al = J.ginfo[x]
        return {((J.g(a, TT, al),), TT): 1 for TT in range(8)}
    def rmul(X, Y):
        Z = {}
        for (m1, T1), c1 in X.items():
            for (m2, T2), c2 in Y.items():
                if T1 & T2: continue
                s = shuffle_sign(T1, T2) * (-1 if (popc(T1) * J.mono_par(m2)) % 2 else 1)
                sg, m = J.norm(m1 + m2)
                if sg: add(Z, (m, T1 | T2), s * sg * c1 * c2)
        return Z
    def psi(m, S):
        X = {((), 0): 1}
        for x in m: X = rmul(X, superfield(x))
        return {mm: c for (mm, T), c in X.items() if T == S}
    signs = {}
    ok = True
    for k in range(1, 5):
        for trial in range(15):
            gens = [J.g(random.randrange(J.dim), 0, (0, 0, 0, random.randrange(3), random.randrange(3))) for _ in range(k)]
            sg, m = J.norm(gens)
            if not sg: continue
            for S in range(8):
                P = {m: 1}
                for v in (2, 1, 0):          # eta_S := eta_0 eta_1 eta_2 restricted to S (apply rightmost first)
                    if (S >> v) & 1: P = J.etapoly(P, v)
                Y = psi(m, S)
                if not Y and not P: continue
                ratio = None
                for mm in set(P) | set(Y):
                    a, b = P.get(mm, 0), Y.get(mm, 0)
                    if b == 0 or a == 0 or (ratio is not None and a * 1.0 / b != ratio) or abs(a) != abs(b):
                        ok = False
                    ratio = a / b if b else None
                signs.setdefault((k % 2, S), set()).add(ratio)
    print(alg, "eta_S phi = +- psi_S(phi) for all tested (k<=4, all S):", ok, "; sign sets by (k mod 2, S):", {key: sorted(v) for key, v in sorted(signs.items())})
