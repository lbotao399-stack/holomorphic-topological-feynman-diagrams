"""Q^2 = 0, [Q,D_u] = 0, Q eta_v + eta_v Q = -2 D_v (Q = 2 x CE), eta_v eta_w + eta_w eta_v = 0, on slices."""
import sys, random
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from vjet import Jets, add

random.seed(1)
for name in ["ab2", "sl2", "aff1", "heis", "gl2"]:
    J = Jets(name)
    monos = []
    for t in [(1, 0, 0, 1, 0), (0, 1, 1, 0, 1), (1, 1, 1, 0, 0), (1, 1, 0, 1, 1), (0, 0, 0, 2, 1), (2, 0, 1, 0, 1)]:
        sl = J.slice(t)
        allm = [m for g in sl for m in sl[g]]
        random.shuffle(allm)
        monos += allm[:25]
    ok = {"Q2": True, "QD": True, "Qeta": True, "etaeta": True}
    for m in monos:
        P = {m: 1}
        if J.Qpoly(J.Qpoly(P)):
            ok["Q2"] = False
        for u in range(5):
            if J.Qpoly(J.Dpoly(P, u)) != J.Dpoly(J.Qpoly(P), u):
                ok["QD"] = False
        for v in range(3):
            L = J.Qpoly(J.etapoly(P, v))
            for k, c in J.etapoly(J.Qpoly(P), v).items():
                add(L, k, c)
            for k, c in J.Dpoly(P, v).items():
                add(L, k, 2 * c)
            if L:
                ok["Qeta"] = False
            for w in range(3):
                L = J.etapoly(J.etapoly(P, v), w)
                for k, c in J.etapoly(J.etapoly(P, w), v).items():
                    add(L, k, c)
                if L:
                    ok["etaeta"] = False
    print(name, "tested monomials", len(monos), ok)
