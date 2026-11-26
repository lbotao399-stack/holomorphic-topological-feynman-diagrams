"""Structural tests of the jet model: Q^2 = 0, [Q, D_u] = 0, {Q, eta_v} = -D_v, eta_v eta_w + eta_w eta_v = 0,
on all generators with derivative orders <= 1 and on random products of generators."""
import random, sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import jets as J

random.seed(1)


def neg(P):
    return {k: -v for k, v in P.items()}


def add(*Ps):
    R = {}
    for P in Ps:
        for k, v in P.items():
            J.padd(R, k, v)
    return R


def run(th):
    gens = J.gens_within(th.dim, (1, 1, 1, 1, 1))
    polys = [{(g,): 1} for g in gens]
    for _ in range(60):
        k = random.randint(2, 3)
        m = {(): 1}
        for g in random.sample(gens, k):
            m = J.pmul(m, {(g,): 1})
        if m:
            polys.append(m)
    ok_QQ = all(not th.Q(th.Q(P)) for P in polys)
    ok_QD = all(not add(th.Q(J.D(P, u)), neg(J.D(th.Q(P), u))) for P in polys for u in range(5))
    ok_eta = True
    ok_etaeta = True
    for vb, u in ((1, 0), (2, 3), (4, 4)):
        for P in polys:
            if add(th.Q(J.eta(P, vb)), J.eta(th.Q(P), vb), J.D(P, u)):
                ok_eta = False
            for wb in (1, 2, 4):
                if add(J.eta(J.eta(P, vb), wb), J.eta(J.eta(P, wb), vb)):
                    ok_etaeta = False
    print(th.name, "polys tested:", len(polys))
    print("  Q^2 = 0:", ok_QQ)
    print("  [Q, D_u] = 0 (u = tau,z1,z2,zb1,zb2):", ok_QD)
    print("  Q eta_v + eta_v Q = -D_v (v = tau,zb1,zb2):", ok_eta)
    print("  eta_v eta_w + eta_w eta_v = 0:", ok_etaeta)
    return ok_QQ and ok_QD and ok_eta and ok_etaeta


if __name__ == "__main__":
    r = all(run(th) for th in (J.abelian(2), J.sl2(), J.gl2()))
    print("all:", r)
