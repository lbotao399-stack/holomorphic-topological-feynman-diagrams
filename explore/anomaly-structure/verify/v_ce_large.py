"""VERIFY supplement: the zero-tests of v_lie_ce.py repeated with large random ranges
(coordinates and momenta in [-R, R]), so that a nonzero polynomial of total degree <= 9 vanishes at
one sample with probability <= 9/(2R+1) (Schwartz-Zippel); 6 samples per test.
(b4), (c5) closedness; exactness solutions (c6 sl2, gl2; c7 relation) re-verified; (c8) cocycles."""
import sys, time, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/verify')
import lie_eval as le
import ce_eval as ce
import v_lie_ce_defs as V

rng = random.Random(161803)
R = 30


def big_args(alg, npts, H):
    return ce.random_args(alg, rng, npts, H, crange=R, prange=R)


def closed(alg, w, H, ns=6):
    return all(ce.d_eval(alg, w, big_args(alg, 5, H)) == 0 for _ in range(ns))


print("R = %d, Schwartz-Zippel bound per sample <= 9/%d = %.3f, for 6 samples <= %.1e" % (R, 2 * R + 1, 9 / (2 * R + 1), (9 / (2 * R + 1)) ** 6))
for g in [le.sl(2), le.gl(2), le.sl(3), le.gl(3)]:
    t0 = time.time()
    W1, W2, Wb = V.colour(g)
    print(" %s: d omega_box = 0 (H=3): %s; d omega_133 = 0: %s; d omega_223 = 0: %s (%.0fs)" % (
        g.name, closed(g, ce.TensorCochain(V.PBOX, Wb), 3), closed(g, ce.TensorCochain(V.P133, W1), 2),
        closed(g, ce.TensorCochain(V.P223, W2), 2), time.time() - t0), flush=True)

mono22 = ce.monomials(2, (2, 2))
for g, which in [(le.sl(2), ["omega", "F_223"]), (le.gl(2), ["omega", "F_223"]), (le.sl(3), ["rel"]), (le.gl(3), ["rel"]), (le.gl(4), ["rel"])]:
    t0 = time.time()
    W1, W2, Wb = V.colour(g)
    w133 = ce.TensorCochain(V.P133, W1)
    w223 = ce.TensorCochain(V.P223, W2)
    F223 = ce.TensorCochain(V.PF, W2)
    tg = {"omega": lambda a: w133(a) + w223(a), "F_223": F223, "rel": lambda a: 2 * w133(a) - 3 * w223(a)}
    etas = [ce.TraceCochain(g, mu, cyc, 3) for mu in mono22 for cyc in ce.trace_types(3)]
    rows, rhs = [], {k: [] for k in which}
    for _ in range(80):
        a = ce.random_args(g, rng, 4, 2)
        rows.append([ce.d_eval(g, f, a) for f in etas])
        for k in which:
            rhs[k].append(tg[k](a))
    for k in which:
        x = ce.solve_Q(rows, rhs[k])
        ok = x is not None
        if ok:
            for _ in range(6):
                a = big_args(g, 4, 2)
                ok &= (sum(xj * ce.d_eval(g, f, a) for xj, f in zip(x, etas) if xj != 0) == tg[k](a))
        print(" %s: %s = d eta, eta solved on 80 small samples, identity at 6 large samples: %s (%.0fs)" % (
            g.name, "omega_133 - (3/2) omega_223 (x2)" if k == "rel" else k, ok, time.time() - t0), flush=True)

# (c8): the kernel of d: C^4 -> C^5 found on samples consists of cocycles
t0 = time.time()
g = le.gl(5)
mono4 = ce.monomials(3, (2, 2))
funcs = [ce.TraceCochain(g, mu, cyc, 4) for mu in mono4 for cyc in ce.trace_types(4)]
rows = [[f(a) for f in funcs] for a in [ce.random_args(g, rng, 4, 2) for _ in range(60)]]
cols = ce.pivots_Q(rows)
basis = [funcs[j] for j in cols]
S5 = [ce.random_args(g, rng, 5, 2) for _ in range(50)]
drows = [[ce.d_eval(g, f, a) for f in basis] for a in S5]
# null space of drows
n, m = len(drows), len(basis)
M = [[Fr(x) for x in r] for r in drows]
piv, r = [], 0
for c in range(m):
    p = next((i for i in range(r, n) if M[i][c] != 0), None)
    if p is None:
        continue
    M[r], M[p] = M[p], M[r]
    pv = M[r][c]
    M[r] = [x / pv for x in M[r]]
    for i in range(n):
        if i != r and M[i][c] != 0:
            f = M[i][c]
            M[i] = [x - f * y for x, y in zip(M[i], M[r])]
    piv.append(c)
    r += 1
free = [c for c in range(m) if c not in piv]
kern = []
for fc in free:
    v = [Fr(0)] * m
    v[fc] = Fr(1)
    for i, pc in enumerate(piv):
        v[pc] = -M[i][fc]
    kern.append(v)
ok = True
for v in kern:
    for _ in range(4):
        a = ce.random_args(g, rng, 5, 2, crange=12, prange=12)
        ok &= (sum(vj * ce.d_eval(g, f, a) for vj, f in zip(v, basis) if vj != 0) == 0)
# independence of the kernel vectors as cochains (values on 4-point samples)
vals = [[sum(vj * f(a) for vj, f in zip(v, basis) if vj != 0) for v in kern] for a in [ce.random_args(g, rng, 4, 2) for _ in range(25)]]
print(" gl5 (c8): dim C^4 basis %d, kernel of d on samples: %d vectors; each is a cocycle at 4 fresh samples (range 12): %s; independent as cochains: rank %d (%.0fs)" % (
    len(basis), len(kern), ok, ce.rank_Q(vals), time.time() - t0), flush=True)
