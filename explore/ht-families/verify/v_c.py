"""Verification of c1, c3 (and the counting part of c2) with own SP decompositions (vsp).
Graph lists: code/enumerate_graphs.py (enumerate_tight + is_sp), counts compared with results/graph_counts.txt.
 C1  U(t(xi)) == prod_b b^{Uexp_b}, Uexp(x of child c of an S-node) = L(c), Uexp(w of child c of a P-node) = -1 - L(c),
     own cycle basis, exact, 2 rational points, for 1 deterministic + 4 random decompositions (2 flattened,
     2 non-flattened) of every SP graph at (D, n) = (2,2), (2,4), (2,6), (3,2), (3,3), (4,2).
 C3  p := sum over blocks floor(h/2) (h = #children with L odd for S-nodes, #children with L even for P-nodes)
     == ceil(n/2) + 2E, E = sum over S-nodes floor(#odd-L children / 2), every decomposition; spread of E over
     the 5 decompositions + 20 further random ones; E distribution.
 C2  exponent structure on the b1 graphs and the b2 graphs: probed leading exponent of every ratio
     (own float density) == -Uexp/2 mod 1, and >= -1/2.
"""
import sys, random, math, time
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import vcore as vc
import vsp
import enumerate_graphs as eg


def graphs(d, n):
    out = []
    for Gx in eg.enumerate_tight(d, n):
        if not eg.is_sp(Gx):
            continue
        nodes = sorted(Gx.nodes())
        out.append((len(nodes), sorted(tuple(sorted((nodes.index(a), nodes.index(b)))) for a, b in Gx.edges())))
    return out


def counts(rc):
    p, E = 0, 0
    for nd, nm in rc.blocks:
        if nd.kind == 'S':
            h = sum(1 for c in nd.children if c.L % 2 == 1)
            E += h // 2
        else:
            h = sum(1 for c in nd.children if c.L % 2 == 0)
        p += h // 2
    return p, E


if __name__ == "__main__":
    rng = random.Random(8)
    for (d, n) in [(2, 2), (2, 4), (2, 6), (3, 2), (3, 3), (4, 2)]:
        t0 = time.time()
        gs = graphs(d, n)
        okU, okp = True, True
        Edist = {}
        Espread = []
        for gi, (nv, edges) in enumerate(gs):
            net = vc.Net(nv, edges)
            C = net.cycle_basis()
            decs = [vsp.sp_decompose(nv, edges)] + [vsp.sp_decompose(nv, edges, rng=random.Random(1000 * gi + s), flatten=(s < 2)) for s in range(4)]
            Es = set()
            for root in decs:
                rc = vsp.RatioCoords(root, net.N)
                assert root.L == n
                for _ in range(2):
                    full = rc.complete(rc.random_point(rng))
                    t = rc.t_of(full)
                    M = [[sum(C[r][e] * C[s][e] * t[e] for e in range(net.N)) for s in range(n)] for r in range(n)]
                    U = vc.det_fr(M)
                    mono = math.prod(full[b] ** ex for b, ex in rc.Uexp.items())
                    okU &= U == mono
                p, E = counts(rc)
                okp &= p == -(-n // 2) + 2 * E
                Es.add(E)
            for s in range(20):
                rc = vsp.RatioCoords(vsp.sp_decompose(nv, edges, rng=random.Random(77 + 1000 * gi + s), flatten=(s % 2 == 0)), net.N)
                p, E = counts(rc)
                okp &= p == -(-n // 2) + 2 * E
                Es.add(E)
            if len(Es) > 1:
                Espread.append((gi, sorted(Es)))
            Edist[min(Es)] = Edist.get(min(Es), 0) + 1
        print("(D,n)=(%d,%d): %d SP graphs ; C1 U == monomial (5 decompositions x 2 points): %s ; C3 p == ceil(n/2)+2E (25 decompositions): %s ; "
              "E distribution (min over decompositions) %s ; graphs with decomposition-dependent E: %s  (%.1fs)" % (
                  d, n, len(gs), okU, okp, dict(sorted(Edist.items())), Espread, time.time() - t0))

    # C2: exponents on the b1 and b2 graphs
    from v_b1 import GR
    from v_b2 import banana, nested
    cases = [("G%d" % k, 6, GR[k], 1) for k in range(5)]
    for lens in [(1, 3, 3, 3, 3), (2, 2, 3, 3, 3)]:
        nv, e, _ = banana(lens)
        cases.append(("banana%s" % (lens,), nv, e, 2))
    for inner in [(1, 3, 3), (2, 2, 3)]:
        nv, e = nested(inner)
        cases.append(("nested%s" % (inner,), nv, e, 2))
    for name, nv, edges, H in cases:
        net = vc.Net(nv, edges)
        C = net.cycle_basis()
        rc = vsp.RatioCoords(vsp.sp_decompose(nv, edges), net.N)
        lamf = np.array([[rng.uniform(-3, 3) for _ in range(H)] for _ in range(nv - 1)])
        f = lambda t, J, full: vc.density_batch(nv, edges, C, t, lamf) * J
        ex = vsp.probe_exponents(rc, f, random.Random(3), npts=4)
        par = all((Fr(ex[b]).limit_denominator(2) + Fr(rc.Uexp[b], 2)).denominator == 1 for b in rc.names)
        p, E = counts(rc)
        print("C2 %-18s exponent == -Uexp/2 mod 1 for all %d ratios: %s ; min exponent %s ; p = %d, E = %d" % (
            name, len(rc.names), par, min(ex.values()), p, E))
