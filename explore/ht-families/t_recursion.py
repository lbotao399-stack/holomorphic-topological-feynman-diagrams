"""Check: rho_n from the super-Gaussian recursion == U^{-1/2} Pf(-1/2 dM M^-1 dM) for C_rec (exact)."""
import sys, random, time
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
import enumerate_graphs as eg
import hf_brute as hb
from hf_recursion import Recursion

def plan_of(red):
    return [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]

def check(name, nv, edges, plan, npts=2, seed=0):
    rec = Recursion(nv, edges, plan)
    C = rec.cycle_basis()
    Gm = hb.G(nv, edges)
    # C must be a basis of the cycle space with C B = 0 and unimodular w.r.t. the reference basis
    B = Gm.incidence()
    assert all(sum(row[e] * B[e][v] for e in range(len(edges))) == 0 for row in C for v in range(nv - 1))
    Cref = Gm.cycles()
    Mt = [[B[e][v] for e in range(len(edges))] for v in range(nv - 1)]
    sgn = hb.fdet(Mt + C)
    rng = random.Random(seed)
    ok = True
    for _ in range(npts):
        t = [Fr(rng.randint(1, 30), rng.randint(1, 12)) for _ in edges]
        t0 = time.time()
        rho_r, U_r = rec.run(t)
        t1 = time.time()
        pf, U = hb.rho_form(Gm, t, C)
        t2 = time.time()
        same = (rho_r == pf) and (U_r == U)
        ok &= same
        print("  %s: recursion == Pf(C_rec): %s, U equal: %s, #components %d (rec %.2fs, Pf %.2fs)" % (
            name, rho_r == pf, U_r == U, len(pf), t1 - t0, t2 - t1))
    print("  det[B^T; C_rec] =", sgn)
    return ok

if __name__ == "__main__":
    gs = [Gx for Gx in eg.enumerate_tight(2, 4) if eg.is_sp(Gx)]
    for k, Gx in enumerate(gs):
        G, _ = htk.from_networkx(Gx, None)
        for prefer in ("large", "small"):
            red = spred.sp_reduce(G, prefer=prefer)
            print("graph %d (%s plan): moves %s" % (k, prefer, [(p, q, ch) for p, q, ch, s in plan_of(red)]))
            check("graph%d" % k, G.nv, G.edges, plan_of(red))
    # theta graphs and bananas
    for (p, q, r) in [(1, 2, 2), (2, 2, 3), (1, 3, 3)]:
        G, interior = htk.theta(p, q, r)
        red = spred.sp_reduce(G, prefer="large")
        print("theta", (p, q, r))
        check("theta", G.nv, G.edges, plan_of(red))
