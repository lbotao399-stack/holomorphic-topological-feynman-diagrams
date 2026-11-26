"""(c) pi-power of K(lambda;0) from the block structure of a series-parallel plan.

p(plan) = sum over simplex blocks of floor(h/2),
  x-block of chain (m,i):  h = #{j : L(e_{m,i,j}) odd},
  w-block of move m:       h = #{i : L(chain_{m,i}) even},
L = loop number of the subnetwork an edge stands for (0 for original edges, sum_{moves inside}(k-1)).
Checks: (1) exponent of every ratio in U equals L (x) resp. -1-L (w)  (exact, from the plan's t-map
and U = det M at a random rational point: U(xi) = prod b^{e_b} on the nose);
(2) half-integer exponents of the engine density = the odd-e_b ratios;
(3) survey of p over all SP graphs (D = 2: n = 2, 4, 6; plans 'small' and 'large').
"""
import sys, random, json
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred, catalog
import enumerate_graphs as eg
import hf_brute as hb
import hf_numeric as hn


def loops_and_p(N0, plan):
    L = {e: 0 for e in range(N0)}
    p = 0
    expU = {}
    for mi, (pp, qq, chains, signs) in enumerate(plan):
        k = len(chains)
        Lch = [sum(L[e] for e in ch) for ch in chains]
        if k >= 2:
            h = sum(1 for x in Lch if x % 2 == 0)
            p += h // 2
            for i in range(k):
                expU["w%d_%d" % (mi + 1, i + 1)] = -1 - Lch[i]
        for i, ch in enumerate(chains):
            if len(ch) >= 2:
                h = sum(1 for e in ch if L[e] % 2 == 1)
                p += h // 2
                for j, e in enumerate(ch):
                    expU["x%d_%d_%d" % (mi + 1, i + 1, j + 1)] = L[e]
        L[N0 + mi] = sum(Lch) + k - 1
    return L, p, expU


def check_U(G, plan, rng):
    rm = hn.RatioMap(G.ne, plan)
    L, p, expU = loops_and_p(G.ne, plan)
    Gm = hb.G(G.nv, G.edges)
    for _ in range(2):
        pt = {}
        for blk in rm.blocks:
            ns = [rng.randint(1, 30) for _ in blk]
            for b, x in zip(blk, ns):
                pt[b] = Fr(x, sum(ns))
        t = [Fr(1)] * G.ne
        for e in range(G.ne):
            for b, ex in rm.texp[e].items():
                t[e] *= pt[b] ** ex
        U = hb.fdet(hb.period(Gm, t, Gm.cycles()))
        mono = Fr(1)
        for b, ex in expU.items():
            mono *= pt[b] ** ex
        if U != mono:
            return False, p, expU
    return True, p, expU


def survey(D, n, rng, verbose=False):
    gs = [Gx for Gx in eg.enumerate_tight(D, n) if eg.is_sp(Gx)]
    rows = []
    for k, Gx in enumerate(gs):
        G, _ = htk.from_networkx(Gx, None)
        ps = {}
        for prefer in ("small", "large"):
            red = spred.sp_reduce(G, D=D, prefer=prefer)
            plan = [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]
            okU, p, expU = check_U(G, plan, rng)
            ps[prefer] = (p, okU)
        degs = sorted(G.degree(v) for v in range(G.nv))
        rows.append((k, degs, G.edges, ps))
    return rows


if __name__ == "__main__":
    rng = random.Random(99)
    out = {}
    for (D, n) in [(2, 2), (2, 4), (2, 6), (3, 2), (4, 2)]:
        rows = survey(D, n, rng)
        dist = {}
        allU = True
        indep = True
        for k, degs, edges, ps in rows:
            allU &= all(v[1] for v in ps.values())
            pv = set(v[0] for v in ps.values())
            indep &= len(pv) == 1
            for v in pv:
                dist[v] = dist.get(v, 0) + 1
        print("D=%d n=%d: %d SP graphs; U = prod b^{L or -1-L} exactly (2 points, both plans): %s; p same for both plans: %s; distribution of p: %s" % (
            D, n, len(rows), allU, indep, dict(sorted(dist.items()))))
        if n == 4 or n == 6:
            for k, degs, edges, ps in rows:
                if any(v[0] != n // 2 for v in ps.values()):
                    print("    p != n/2: graph %d degrees %s edges %s p(small,large) = %s" % (
                        k, degs, edges, (ps["small"][0], ps["large"][0])))
        out["%d,%d" % (D, n)] = [(k, degs, edges, {kk: v[0] for kk, v in ps.items()}) for k, degs, edges, ps in rows]
    json.dump(out, open("c_pi_survey.json", "w"))
