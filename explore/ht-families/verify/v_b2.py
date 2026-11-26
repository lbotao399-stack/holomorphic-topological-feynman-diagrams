"""Verification of b2 and of the banana closed form of a5 at (H,T) = (2,1), n = 4 (|V| = 10, |I| = 13).
Graphs rebuilt from the construction stated in b21.py (banana(lens): p = 0, q = root = nv-1, chains in order,
interior vertices numbered consecutively; nested: theta(1,3,3) with its edge 2 replaced by a banana between
its end points, then q <-> last label).
 W1  exact: engine density at the explorer's printed source sets (plan 'large'), own pointwise check
     (3 rational points, sign and square), own Dirichlet integration; compared with the claimed values.
 W2  banana closed form ((n-1)!!/2^{n/2}) D_Gamma prod_P Gamma(H+3/2-a_P)/(a_P-1)! / Gamma(H+(n+1)/2) at the
     printed and at 2 fresh source sets (exact K from W1-type computation).
 W3  own QMC (own SP decomposition, own float density), one source set per graph.
 W4  tightness: |I| = 3n+1, |V| = 2n+2, and 2|I_S| <= 3|V_S| - 4 for every vertex subset (own check).
"""
import sys, random, time, math, itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import vcore as vc
import vsp
import vdir
from v_b1 import engine_exact, own_qmc


def banana(lens):
    nv = 2 + sum(a - 1 for a in lens)
    p, q = 0, nv - 1
    edges, v, ints = [], 1, []
    for a in lens:
        prev, iv = p, []
        for j in range(a):
            if j == a - 1:
                nxt = q
            else:
                nxt = v
                iv.append(v)
                v += 1
            edges.append((prev, nxt))
            prev = nxt
        ints.append(iv)
    return nv, edges, ints


def nested(inner):
    nv, edges, _ = banana((1, 3, 3))
    u, v = edges[2]
    out = [x for i, x in enumerate(edges) if i != 2]
    nxt = nv
    for a in inner:
        prev = u
        for j in range(a):
            if j == a - 1:
                w = v
            else:
                w = nxt
                nxt += 1
            out.append((prev, w))
            prev = w
    q = nv - 1
    perm = list(range(nxt))
    perm[q], perm[nxt - 1] = nxt - 1, q
    return nxt, [(perm[a], perm[b]) for a, b in out]


def tight_check(nv, edges, D=3):
    N = len(edges)
    n = N - nv + 1
    ok = N == D * n + 1 and nv == (D - 1) * n + 2
    for k in range(2, nv + 1):
        for Sset in itertools.combinations(range(nv), k):
            s = set(Sset)
            e = sum(1 for a, b in edges if a in s and b in s)
            if e and (D - 1) * e > D * k - D - 1:
                return False
    return ok


def D_banana(lam, ints, H):
    k = len(ints)
    cols = [(P, v) for P in range(k) for v in ints[P]]
    rows = []
    for i in range(H):
        for a in range(k - 1):
            rows.append([lam[v][i] if P == a else (-lam[v][i] if P == k - 1 else Fr(0)) for (P, v) in cols])
    return vc.det_fr(rows)


def banana_formula(lam, ints, lens, H, n):
    pref = sp.Rational(sp.factorial2(n - 1), 2 ** (n // 2))
    g = sp.Integer(1)
    for a in lens:
        g *= sp.gamma(sp.Rational(2 * H + 3 - 2 * a, 2)) / sp.factorial(a - 1)
    g /= sp.gamma(sp.Rational(2 * H + n + 1, 2))
    Dg = D_banana(lam, ints, H)
    return sp.nsimplify(pref * sp.Rational(Dg.numerator, Dg.denominator) * g)


PRINTED = {
    "banana(1,3,3,3,3)": [
        ([[7, 12], [12, -12], [10, 2], [-4, 11], [-5, 6], [-9, -2], [-12, -12], [-12, 8], [5, -12]], sp.Rational(1044576, 5) * sp.pi ** 2),
        ([[-2, 2], [7, -12], [-5, 8], [-7, 5], [6, -7], [-10, 5], [-4, -11], [9, -10], [-10, -12]], sp.Rational(77064, 7) * sp.pi ** 2)],
    "banana(2,2,3,3,3)": [
        ([[-11, -10], [-10, -1], [-7, 11], [9, -3], [-4, 7], [-6, 7], [-11, 6], [9, -7], [1, 8]], sp.Rational(19327113, 280) * sp.pi ** 2),
        ([[9, 2], [2, -1], [6, 11], [5, 11], [2, 3], [9, -5], [-2, 10], [-7, 7], [-4, 12]], -3136 * sp.pi ** 2)],
    "nested theta(1,3,3)": [
        ([[6, 1], [6, -5], [-2, 9], [-12, -4], [7, 9], [10, -7], [10, -2], [5, 6], [6, -9]], 139072 * sp.pi ** 2),
        ([[-9, 7], [4, -4], [1, 8], [11, 10], [-5, -3], [1, -4], [4, -3], [5, -2], [-12, 1]], -sp.Rational(58327, 200) * sp.pi ** 2)],
    "nested theta(2,2,3)": [
        ([[-1, -4], [12, -7], [-9, -4], [-6, -12], [8, -4], [-4, -6], [-7, -3], [-3, 8], [11, -1]], -sp.Rational(353210, 3) * sp.pi ** 2),
        ([[-10, 1], [-6, 8], [8, 2], [-4, -7], [-1, 1], [11, 6], [-2, 8], [5, -6], [-2, -9]], -sp.Rational(587279, 300) * sp.pi ** 2)],
}


def kexpr(K):
    return sum(sp.Rational(v.numerator, v.denominator) * sp.pi ** p for p, v in K.items())


if __name__ == "__main__":
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    rng = random.Random(77)
    cases = []
    for lens in [(1, 3, 3, 3, 3), (2, 2, 3, 3, 3)]:
        nv, edges, ints = banana(lens)
        cases.append(("banana(%s)" % ",".join(map(str, lens)), nv, edges, ints, lens))
    for inner in [(1, 3, 3), (2, 2, 3)]:
        nv, edges = nested(inner)
        cases.append(("nested theta(%s)" % ",".join(map(str, inner)), nv, edges, None, None))
    for name, nv, edges, ints, lens in cases:
        print("=" * 90)
        print("%s: nv %d, edges %s" % (name, nv, edges))
        print("  W4 tight + subgraph bound (D=3): %s" % tight_check(nv, edges))
        for idx, (lam_i, claim) in enumerate(PRINTED[name]):
            lam = [[Fr(x) for x in l] for l in lam_i]
            t0 = time.time()
            K, ok, dens = engine_exact(nv, edges, lam, rng, npts=3)
            Ke = kexpr(K)
            print("  printed set %d: W1 pointwise %s ; own Dirichlet K = %s ; claim %s ; equal %s  (%.1fs)" % (
                idx, ok, Ke, claim, sp.simplify(Ke - claim) == 0, time.time() - t0))
            if lens is not None:
                f = banana_formula(lam, ints, lens, 2, 4)
                print("     W2 closed form = %s ; K / closed form = %s" % (f, sp.simplify(Ke / f)))
            if idx == 0:
                t0 = time.time()
                lamf = np.array([[float(x) for x in l] for l in lam])
                mean, err, nex = own_qmc(nv, edges, lamf, m, 16, 300 + len(name))
                ex = float(claim)
                print("     W3 own QMC 2^%d x 16: %.10g +- %.2g ; claim %.10g ; diff/err %.2f ; rel diff %.1e ; %d points exact  (%.1fs)" % (
                    m, mean, err, ex, (mean - ex) / err, abs(mean - ex) / abs(ex), nex, time.time() - t0))
        if lens is not None:
            for trial in range(2):
                lam = [[Fr(rng.randint(-9, 9)) for _ in range(2)] for _ in range(nv - 1)]
                K, ok, dens = engine_exact(nv, edges, lam, rng, npts=3)
                Ke = kexpr(K)
                f = banana_formula(lam, ints, lens, 2, 4)
                print("  fresh set %d %s: W1 pointwise %s ; K = %s ; W2 K / closed form = %s" % (
                    trial, [[int(x) for x in l] for l in lam], ok, Ke, sp.simplify(Ke / f) if f != 0 else ("K=%s, form=0" % Ke)))
