"""t2_moves.py -- pointwise exact checks of the moves on the topological factor with sources kappa.

 chain moves (series k=1, parallel a_i=1, general chains): every pair (p,q) of every admissible graph
   D <= 4, n <= 3 with a chain bundle (root not interior), all chains and single chains; bananas and
   multigraphs; random rational t and kappa at all non-root vertices.
 star-triangle (+ two parallel merges, LLC App. E): bitriangle and four larger graphs, random t, kappa.
"""
import random, time
from fractions import Fraction as Fr
import tb
from tb import Graph
import graphs as gr
from moves import find_chains, chain_move_check, chain_interior, star_triangle_check

rng = random.Random(314159)


def rand_kappa(G, zero=False):
    if zero:
        return [Fr(0)] * (G.nv - 1)
    return [Fr(rng.randint(-6, 6), rng.randint(1, 5)) for _ in range(G.nv - 1)]


def chain_cases(G):
    out = []
    for p in range(G.nv):
        for q in range(G.nv):
            if p == q:
                continue
            chs = find_chains(G, p, q)
            if not chs:
                continue
            ints = [chain_interior(G, ch, p) for ch in chs]
            if any(G.root in vs for vs in ints):
                continue
            if any(len(vs) == 0 for vs in ints) and len(chs) == 1:
                continue        # a single direct edge: nothing to reduce
            if len(chs) >= 1 and p < q:
                out.append((p, q, chs))
                if len(chs) >= 2:
                    for ch in chs:
                        if len(ch) >= 2:
                            out.append((p, q, [ch]))
    return out


def main():
    t0 = time.time()
    print(__doc__)
    tot = {}
    fails = 0
    graphsets = []
    for D in (1, 2, 3, 4):
        for n in (1, 2, 3):
            graphsets += [((D, n), G) for G in gr.admissible(D, n)]
    # multigraphs: theta with a doubled edge, banana with subdivided strands, two-level bundles
    extra = [
        Graph(4, [(0, 1), (0, 1), (1, 3), (0, 2), (2, 3), (0, 3)], "multi1"),
        Graph(5, [(0, 1), (1, 4), (0, 2), (2, 4), (0, 3), (3, 4), (0, 4), (0, 4)], "banana_sub"),
        Graph(6, [(0, 1), (1, 2), (2, 5), (0, 3), (3, 5), (0, 4), (4, 5), (1, 3), (2, 4)], "mixed"),
    ]
    graphsets += [(("x", 0), G) for G in extra]
    for key, G in graphsets:
        cases = chain_cases(G)
        for (p, q, chs) in cases:
            for zero in (False, True):
                t = tb.rand_t(G.ne, rng)
                kap = rand_kappa(G, zero)
                res = chain_move_check(G, p, q, chs, t, kap)
                kk = (key, "k=%d" % len(chs), "kappa=0" if zero else "kappa")
                s = tot.setdefault(kk, [0, 0])
                s[0] += 1
                good = all(v for kname, v in res.items() if not kname.startswith("_"))
                s[1] += good
                if not good:
                    fails += 1
                    print("FAIL", G.name, p, q, [len(c) for c in chs], res)
    print("chain moves: (graph class (D,n), number of chains k, sources) -> passed/total")
    for kk in sorted(tot, key=str):
        print("   %-40s %d/%d" % (kk, tot[kk][1], tot[kk][0]))
    # star-triangle
    print("\nstar-triangle (+ parallel merges), edges a=vu, b=vw, c=vx, d=ux, f=wx:")
    st = []
    # v=0,u=1,w=2,x=3 root
    st.append((Graph(4, [(0, 1), (0, 2), (0, 3), (1, 3), (2, 3)], "bitriangle"), (0, 1, 2, 3, 0, 1, 2, 3, 4)))
    st.append((Graph(5, [(0, 1), (0, 2), (0, 3), (1, 3), (2, 3), (1, 4), (2, 4), (3, 4)], "K5minus"), (0, 1, 2, 3, 0, 1, 2, 3, 4)))
    st.append((Graph(5, [(0, 1), (0, 2), (0, 3), (1, 3), (2, 3), (1, 4), (2, 4)], "ST_n3"), (0, 1, 2, 3, 0, 1, 2, 3, 4)))
    st.append((Graph(6, [(0, 1), (0, 2), (0, 3), (1, 3), (2, 3), (1, 4), (4, 5), (2, 5), (3, 5)], "ST_n4"), (0, 1, 2, 3, 0, 1, 2, 3, 4)))
    st.append((Graph(6, [(0, 1), (0, 2), (0, 3), (1, 3), (2, 3), (1, 4), (2, 4), (3, 4), (4, 5), (1, 5)], "ST_n5"), (0, 1, 2, 3, 0, 1, 2, 3, 4)))
    stok = True
    for G, (v, u, w, x, a, b, c, d, f) in st:
        for zero in (False, True):
            for rep in range(2):
                t = tb.rand_t(G.ne, rng)
                kap = rand_kappa(G, zero)
                res = star_triangle_check(G, v, u, w, x, a, b, c, d, f, t, kap)
                good = all(val for kname, val in res.items() if not kname.startswith("_"))
                stok &= good
                print("   %-11s n=%d kappa%s: %s  (#components %d)" % (G.name, G.loops, "=0" if zero else "!=0",
                      ", ".join("%s: %s" % (k_, v_) for k_, v_ in res.items() if not k_.startswith("_")), res["_ncomp"]))
    print("\nsummary: chain-move failures %d; star-triangle all %s; time %.1fs" % (fails, stok, time.time() - t0))


if __name__ == "__main__":
    main()
