"""(a) Graphs at order hbar^3 other than Theta(J_1), (H,T) = (2,1):
 (1) one-loop graphs (b_1 = 1) with two J_1 vertices or one I_2 vertex: the only tight one-loop graph with
     |I| = Dn+1 = 4 is polygon(4); the T = 1 form rho_1 = U^{-1/2}<psi^tau_1> vanishes identically, so every
     holomorphic current insertion gives 0:  spred density (symbolic sources) == 0; brute-force top form at random
     rational t empty; kernels with insertions [y_e y_e'] (J_1 at an internal pair), [y y]^2 (I_2 type),
     [lam y][y y] (two J_1 vertices) == 0.
 (2) three-loop graphs with S_0 vertices (b_1 = 3): all 16 tight graphs of htk.catalog.tight_graphs(3, 3):
     brute-force top form rho_3 at random rational t empty (rho_odd = 0); spred kernels 0 where series-parallel.
Non-tight graphs: K = 0 by form degree (insertions are 0-forms and do not change it)."""
import random
from common import *
from htk import catalog, brute

rng = random.Random(11)
G = htk.polygon(4)
lamS = [[sp.Symbol('l%d_%d' % (v, c)) for c in range(2)] for v in range(G.nv - 1)]
red = htk.sp_reduce(G, prefer="large")
dens = red.density(lamS, T=1, check=2)
print("polygon(4) edges %s: T=1 density polynomial factors %s ; num == 0: %s" % (G.edges, dens.factors, dens.num == 0))
for t_ in range(3):
    t = [Fr(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(G.ne)]
    lam = [[Fr(rng.randint(-9, 9)) for _ in range(2)] for _ in range(G.nv - 1)]
    d = brute.density(G, t, lam, T=1)
    print("  brute density at t = %s: coeff %s" % ([str(x) for x in t], d.coeff))


def brk(a, b):
    return a[0] * b[1] - a[1] * b[0]


ins_list = {
    "[y0 y1]": lambda y: brk(y[0], y[1]),
    "[y0 y1][y2 y3]": lambda y: brk(y[0], y[1]) * brk(y[2], y[3]),
    "[y0 y1]^2": lambda y: brk(y[0], y[1]) ** 2,
    "[y0 y3][y1 y2]": lambda y: brk(y[0], y[3]) * brk(y[1], y[2]),
}
for nm, ins in ins_list.items():
    K = red.integrate(lamS, 1, insertion=ins, dens=dens)
    print("  polygon(4) insertion %s: K = %s" % (nm, K))

Gs = catalog.tight_graphs(3, 3)
print("three-loop tight graphs at D = 3: %d" % len(Gs))
for G in Gs:
    deg = sorted(sum(1 for e in G.edges if v in e) for v in range(G.nv))
    zero = True
    for t_ in range(2):
        t = [Fr(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(G.ne)]
        lam = [[Fr(rng.randint(-9, 9)) for _ in range(2)] for _ in range(G.nv - 1)]
        d = brute.density(G, t, lam, T=1)
        zero &= (d.coeff == 0)
    try:
        lam = [[Fr(rng.randint(-9, 9)) for _ in range(2)] for _ in range(G.nv - 1)]
        K = htk.kernel_exact(G, lam, T=1)
    except Exception as ex:
        K = "not series-parallel"
    print("  degrees %s: brute top form zero at 2 random points: %s ; spred K = %s" % (deg, zero, K))
