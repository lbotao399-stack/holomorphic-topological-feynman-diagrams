"""(b) n = 4 at (H,T) = (2,1) (D = 3, |V| = 10, |I| = 13): the two 5-chain bananas (1,3,3,3,3), (2,2,3,3,3)
and two nested graphs (theta(1,3,3) with the middle edge of a 3-chain replaced by theta(1,3,3) resp.
theta(2,2,3)).  Tightness and the subgraph bound checked with catalog.sparse_ok.
Sources: random integer 2-vectors at all non-root vertices (symbolic sources too slow at |V| = 10).
 E1 engine kernel (exact, PiNumber), engine self-check
 E2 my exact pointwise check of the engine density (own brute force, own t-map)
 E3 QMC integral of the float brute-force density, Dirichlet(alpha+1) and Dirichlet(1/2) sampling
 E4 banana closed form K = s (3/4) D_Gamma prod_P Gamma(H+3/2-a_P)/(a_P-1)! / Gamma(sum_P (H+3/2-a_P))
"""
import sys, random, time
from fractions import Fraction as Fr
import numpy as np
import networkx as nx
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred, catalog
import hf_brute as hb
import hf_numeric as hn
import c_pi


def gen_lam(rng, nvert, H):
    """random integer H-vectors with all 2x2 minors between distinct vertices nonzero (H = 2)"""
    while True:
        lam = [[Fr(rng.randint(-12, 12)) for _ in range(H)] for _ in range(nvert)]
        if H == 2 and any(lam[a][0] * lam[b][1] - lam[a][1] * lam[b][0] == 0 for a in range(nvert) for b in range(a + 1, nvert)):
            continue
        if any(x == 0 for l in lam for x in l):
            continue
        return lam


def banana(lens):
    nv = 2 + sum(a - 1 for a in lens)
    p, q = 0, nv - 1
    edges, v = [], 1
    interiors = []
    for a in lens:
        prev, ints = p, []
        for j in range(a):
            if j == a - 1:
                nxt = q
            else:
                nxt = v
                ints.append(v)
                v += 1
            edges.append((prev, nxt))
            prev = nxt
        interiors.append(ints)
    return nv, edges, interiors


def replace_edge(nv, edges, e, lens):
    """replace edge e = (u, v) by a banana with chains lens between u and v"""
    u, v = edges[e]
    out = [x for i, x in enumerate(edges) if i != e]
    nxt = nv
    for a in lens:
        prev = u
        for j in range(a):
            if j == a - 1:
                w = v
            else:
                w = nxt
                nxt += 1
            out.append((prev, w))
            prev = w
    return nxt, out


def normalize(nv, edges):
    """relabel so that vertices are 0..nv-1 with the last label as root (keep), sort edges"""
    E = sorted(tuple(sorted(x)) for x in edges)
    return nv, E


def tight_ok(nv, edges, D):
    Gx = nx.Graph()
    Gx.add_nodes_from(range(nv))
    Gx.add_edges_from(edges)
    return catalog.tight_ok(Gx, D) and Gx.number_of_edges() == len(edges)


def run(name, nv, edges, H=2, D=3, interiors=None, lens=None, seed=0):
    G = htk.Graph(nv, edges)
    print("=" * 100)
    print("%s: |V| = %d, |I| = %d, n = %d, tight + subgraph bound (D=3): %s" % (name, nv, len(edges), G.loops, tight_ok(nv, edges, D)))
    red = spred.sp_reduce(G, D=D, prefer="large")
    plan = [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]
    L, p, expU = c_pi.loops_and_p(G.ne, plan)
    print("  plan:", [(pp, qq, ch) for pp, qq, ch, s in plan], " p =", p)
    rng = random.Random(seed)
    Gm = hb.G(nv, edges)
    C = Gm.cycles()
    rm = hn.RatioMap(G.ne, plan)
    for trial in range(2):
        lam = gen_lam(rng, nv - 1, H)
        t0 = time.time()
        dens = red.density(lam, T=1, check=1)
        K = dens.integrate()
        print("  lambda = %s" % [[int(x) for x in l] for l in lam])
        print("    E1  K = %s   (engine, %.1fs; exponents %s)" % (K, time.time() - t0, dens.exponents()))
        ok = True
        for _ in range(3):
            pt = red.random_point(rng)
            t_own = [Fr(1)] * G.ne
            for e in range(G.ne):
                for b, ex in rm.texp[e].items():
                    t_own[e] *= Fr(pt[b]) ** ex
            rows = [t_own]
            for c in rm.coords:
                blk = next(bl for bl in rm.blocks if c in bl)
                rows.append([t_own[e] * (Fr(rm.texp[e].get(c, 0)) / Fr(pt[c]) - Fr(rm.texp[e].get(blk[0], 0)) / Fr(pt[blk[0]]))
                             for e in range(G.ne)])
            J = hb.fdet(rows)
            coeff, U = hb.density(Gm, t_own, lam, T=1, C=C)
            lhs = coeff * abs(J)
            rat, rad = dens.evaluate(tuple(Fr(pt[c]) for c in red.coords))
            ok &= ((lhs > 0) - (lhs < 0) == (rat > 0) - (rat < 0)) and lhs ** 2 / U == rat ** 2 * rad
        print("    E2  pointwise g == D_brute |J| at 3 rational points: %s" % ok)
        exact = float(K)
        lamf = np.array([[float(x) for x in l] for l in lam])
        ex = dens.exponents()
        alpha = {b: max(0.5, float(ex.get(b, 0)) + 1) for b in rm.names}
        t1 = time.time()
        mean, err, _ = hn.qmc_integral(nv, edges, C, plan, lamf, T=1, m=13, nrep=16, seed=7000 + 10 * trial + seed, alpha=alpha)
        print("    E3  QMC Dirichlet(alpha+1), 2^13 x 16: %.10g +- %.2g ; exact %.12g ; rel diff %.2g  (%.1fs)" % (
            mean, err, exact, abs(mean - exact) / abs(exact) if exact else float("nan"), time.time() - t1))
        t1 = time.time()
        mean, err, _ = hn.qmc_integral(nv, edges, C, plan, lamf, T=1, m=13, nrep=16, seed=8000 + 10 * trial + seed, alpha=None)
        print("    E3' QMC Dirichlet(1/2), 2^13 x 16:       %.10g +- %.2g ; rel diff %.2g  (%.1fs)" % (
            mean, err, abs(mean - exact) / abs(exact) if exact else float("nan"), time.time() - t1))
        if lens is not None:
            # E4: closed form with D_Gamma = det of the 2*4 x 2*4 matrix [[Lam_1 -Lam_5 ...]] (fibre forms), see report
            cs = [Fr(0)]
            num = Fr(1)
            g = sp.Integer(1)
            for a in lens:
                g *= sp.gamma(sp.Rational(2 * H + 3 - 2 * a, 2)) / sp.factorial(a - 1)
            g /= sp.gamma(sum(sp.Rational(2 * H + 3 - 2 * a, 2) for a in lens))
            Dg = D_banana(lam, interiors, H)
            val = sp.Rational(3, 4) * Dg * g
            print("    E4  (3/4) D_Gamma prod Gamma(H+3/2-a)/(a-1)! / Gamma(sum) = %s ;  K / that = %s" % (
                sp.nsimplify(val), sp.nsimplify(K.sympy() / val) if val != 0 else None))


def D_banana(lam, interiors, H):
    """D_Gamma = det of the (H*4) x (sum (a_P - 1)) matrix of fibre forms:
    rows (i, a), a = 1..4 (loop a = chain a - chain 5); block columns per chain P (interior sources):
    entry for chain P in row (i,a): Lambda^i_P if P = a, -Lambda^i_P if P = 5, 0 otherwise;
    Lambda_P = (lambda at interior vertices of P, in path order)."""
    k = len(interiors)
    cols = [(P, v) for P in range(k) for v in interiors[P]]
    rows = []
    for i in range(H):
        for a in range(k - 1):
            row = []
            for (P, v) in cols:
                if P == a:
                    row.append(lam[v - 1 + 1 - 1][i] if False else lam_at(lam, v, i))
                elif P == k - 1:
                    row.append(-lam_at(lam, v, i))
                else:
                    row.append(Fr(0))
            rows.append(row)
    return sp.Rational(hb.fdet(rows))


def lam_at(lam, v, i):
    return Fr(lam[v][i])


if __name__ == "__main__":
    which = sys.argv[1].split(",") if len(sys.argv) > 1 else ["b1", "b2", "n1", "n2"]
    for w in which:
        if w == "b1":
            nv, edges, ints = banana((1, 3, 3, 3, 3))
            run("banana(1,3,3,3,3)", nv, edges, interiors=ints, lens=(1, 3, 3, 3, 3), seed=1)
        if w == "b2":
            nv, edges, ints = banana((2, 2, 3, 3, 3))
            run("banana(2,2,3,3,3)", nv, edges, interiors=ints, lens=(2, 2, 3, 3, 3), seed=2)
        if w in ("n1", "n2"):
            nv, edges, ints = banana((1, 3, 3))
            # middle edge of the second chain (chain B = edges 1,2,3): edge index 2
            inner = (1, 3, 3) if w == "n1" else (2, 2, 3)
            nv2, edges2 = replace_edge(nv, edges, 2, inner)
            # root = last vertex: make the outer terminal q the root (relabel q <-> nv2-1)
            q = nv - 1
            perm = list(range(nv2))
            perm[q], perm[nv2 - 1] = nv2 - 1, q
            edges3 = [(perm[a], perm[b]) for a, b in edges2]
            run("theta(1,3,3) with B_2 -> theta%s" % (inner,), nv2, edges3, seed=3 if w == "n1" else 4)

