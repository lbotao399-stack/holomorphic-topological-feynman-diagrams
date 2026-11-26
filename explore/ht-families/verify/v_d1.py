"""Verification of d1: G2 = triangle with an ear on each side, (H,T) = (1,1), n = 4.
Own decomposition built by hand (small-plan structure): ear A = P((0,1), S((0,2),(1,2))),
ear B = P((0,4), S((0,5),(4,5))), ear C = P((1,4), S((1,3),(3,4))), centre = P(A, S(B, C)).
 Z1  exact (sympy, symbolic sources l0..l4, ear ratios fixed at random rationals): with the central coordinates
     X (x of ear C in S(B,C)) and W (w of S(B,C) in the centre), every current y_e satisfies
     d(y_e, l)/d(X, W) = 0, l = current from vertex 0 into ear A; dl != 0; sqrt(Pi) on the central square,
     Pi = t_A t_B t_C / (t_A + t_B + t_C)^3 with the gauge t_centre = 1.
 Z2  central-square integral at fixed ear point (float): int g(ears, X, W) e^{-sum y.z} F(y) dX dW with
     X = sin^2 phi, W = sin^2 chi and 80 x 80 Gauss-Legendre nodes, compared with int |.|; control insertion
     X (not a function of l), which must not vanish.  4 ear points x 3 (z, F) choices.
 Z3  full 8-dim QMC of K(lambda; z) at z = 0 and at random z (own coordinates, own density).
"""
import sys, random, math, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
import vcore as vc
import vsp

edges = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (3, 4), (4, 5)]
nv, N = 6, 9
E = lambda e: vsp.Node('E', edge=e)
earA = vsp.Node('P', [E(0), vsp.Node('S', [E(1), E(4)])])
earB = vsp.Node('P', [E(2), vsp.Node('S', [E(3), E(8)])])
earC = vsp.Node('P', [E(6), vsp.Node('S', [E(5), E(7)])])
sBC = vsp.Node('S', [earB, earC])
root = vsp.Node('P', [earA, sBC])
rc = vsp.RatioCoords(root, N)
blk = {id(nd): nm for nd, nm in rc.blocks}
Xn = blk[id(sBC)][1]     # x of ear C
Wn = blk[id(root)][1]    # w of S(B, C)
central = {Xn, Wn, blk[id(sBC)][0], blk[id(root)][0]}
earnames = [b for b in rc.names if b not in central]
net = vc.Net(nv, edges)
C = net.cycle_basis()
rng = random.Random(5)
print("blocks:", [(nd.kind, nm) for nd, nm in rc.blocks], " central X =", Xn, " W =", Wn)


def ear_point(rng):
    pt = rc.random_point(rng)
    return {b: v for b, v in pt.items() if b not in (Xn, Wn)}


if __name__ == "__main__":
    # ---------------- Z1 (exact, forward-mode derivatives with dual numbers over Fractions)
    class Dual:
        __slots__ = ("v", "d")

        def __init__(self, v, d=(0, 0)):
            self.v, self.d = Fr(v), tuple(Fr(x) for x in d)

        def _c(self, o):
            return o if isinstance(o, Dual) else Dual(o)

        def __add__(self, o):
            o = self._c(o)
            return Dual(self.v + o.v, (self.d[0] + o.d[0], self.d[1] + o.d[1]))
        __radd__ = __add__

        def __sub__(self, o):
            o = self._c(o)
            return Dual(self.v - o.v, (self.d[0] - o.d[0], self.d[1] - o.d[1]))

        def __rsub__(self, o):
            return self._c(o) - self

        def __mul__(self, o):
            o = self._c(o)
            return Dual(self.v * o.v, (self.d[0] * o.v + self.v * o.d[0], self.d[1] * o.v + self.v * o.d[1]))
        __rmul__ = __mul__

        def __truediv__(self, o):
            o = self._c(o)
            return Dual(self.v / o.v, ((self.d[0] * o.v - self.v * o.d[0]) / o.v ** 2, (self.d[1] * o.v - self.v * o.d[1]) / o.v ** 2))

        def __rtruediv__(self, o):
            return self._c(o) / self

        def __neg__(self):
            return Dual(-self.v, (-self.d[0], -self.d[1]))

        def __pow__(self, k):
            out = Dual(1)
            base = self if k >= 0 else Dual(1) / self
            for _ in range(abs(k)):
                out = out * base
            return out

        def __eq__(self, o):
            o = self._c(o)
            return self.v == o.v and self.d == o.d

        def __ne__(self, o):
            return not self.__eq__(o)

    def dsolve(A, b):
        n = len(A)
        Mx = [list(A[i]) + [b[i]] for i in range(n)]
        for c in range(n):
            p = next(r for r in range(c, n) if Mx[r][c].v != 0)
            Mx[c], Mx[p] = Mx[p], Mx[c]
            piv = Mx[c][c]
            Mx[c] = [x / piv for x in Mx[c]]
            for r in range(n):
                if r != c and (Mx[r][c].v != 0 or Mx[r][c].d != (0, 0)):
                    f = Mx[r][c]
                    Mx[r] = [a - f * bb for a, bb in zip(Mx[r], Mx[c])]
        return [Mx[i][n] for i in range(n)]

    t0 = time.time()
    B = net.Bmat()
    okall, dlall, ndep = True, True, set()
    for trial in range(3):
        ep = ear_point(rng)
        for cp in range(3):
            Xv0 = Fr(rng.randint(1, 96), 97)
            Wv0 = Fr(rng.randint(1, 96), 97)
            vals = {b: Dual(v) for b, v in ep.items()}
            vals[Xn], vals[Wn] = Dual(Xv0, (1, 0)), Dual(Wv0, (0, 1))
            for nd, nm in rc.blocks:
                vals[nm[0]] = 1 - sum((vals[b] for b in nm[1:]), Dual(0))
            t = []
            for e in range(N):
                v = Dual(1)
                for b, ex in rc.texp[e].items():
                    v = v * vals[b] ** ex
                t.append(v)
            Lam = [[sum((Dual(B[e][a] * B[e][b]) / t[e] for e in range(N)), Dual(0)) for b in range(5)] for a in range(5)]
            # response to unit sources at each vertex v: y[v][e]
            Y = []
            for v in range(5):
                phi = dsolve(Lam, [Dual(int(u == v)) for u in range(5)])
                Y.append([sum((B[e][u] * phi[u] for u in range(5)), Dual(0)) / t[e] for e in range(N)])
            Lr = [Y[v][0] + Y[v][1] for v in range(5)]    # current out of vertex 0 into ear A
            for e in range(N):
                for v in range(5):
                    for v2 in range(5):
                        s = (Y[v][e].d[0] * Lr[v2].d[1] - Y[v][e].d[1] * Lr[v2].d[0]) + (Y[v2][e].d[0] * Lr[v].d[1] - Y[v2][e].d[1] * Lr[v].d[0])
                        okall &= s == 0
                for v in range(5):
                    if Y[v][e].d[1] != 0:
                        ndep.add(e)
            dlall &= any(Lr[v].d != (0, 0) for v in range(5))
            # Pi on the central square
            tA = 1 / (1 / t[0].v + 1 / (t[1].v + t[4].v))
            tB = 1 / (1 / t[2].v + 1 / (t[3].v + t[8].v))
            tC = 1 / (1 / t[6].v + 1 / (t[5].v + t[7].v))
            Pi = tA * tB * tC / (tA + tB + tC) ** 3
            okall &= Pi == Xv0 * (1 - Xv0) * Wv0 * (1 - Wv0) ** 2
    print("Z1 3 ear points x 3 central points: d(y_e, l)/d(X, W) == 0 as a quadratic form in the sources for all 9 edges, and "
          "Pi == X(1-X) W (1-W)^2: %s ; dl != 0: %s ; edges with dy_e/dW != 0: %s  (%.1fs)" % (okall, dlall, sorted(ndep), time.time() - t0))

    # ---------------- Z2
    nq = 80
    xg, wg = np.polynomial.legendre.leggauss(nq)
    ph = (xg + 1) * np.pi / 4
    wph = wg * np.pi / 4
    PH, CH = np.meshgrid(ph, ph, indexing="ij")
    WW = np.outer(wph, wph).ravel()
    Xv = np.sin(PH).ravel() ** 2
    Wv = np.sin(CH).ravel() ** 2
    jac_sub = (2 * np.sin(PH) * np.cos(PH) * 2 * np.sin(CH) * np.cos(CH)).ravel()
    lamf = np.array([[3.], [-2.], [5.], [1.], [-4.]])
    for trial in range(4):
        ep = ear_point(rng)
        full = {b: np.full(Xv.shape, float(v)) for b, v in ep.items()}
        full[Xn], full[Wn] = Xv, Wv
        for nd, nm in rc.blocks:
            full[nm[0]] = 1 - sum(full[b] for b in nm[1:])
        t, J = rc.t_and_jac(full)
        D, y = vc.density_batch(nv, edges, C, t, lamf, return_y=True)
        g = D * np.abs(J) * jac_sub * WW
        y = y[:, :, 0]
        zs = np.array([rng.uniform(-0.3, 0.3) for _ in range(N)])
        cases = [("z=0, F=1", np.ones_like(g)),
                 ("z random, F=1", np.exp(-(y * zs[None, :]).sum(axis=1))),
                 ("z random, F=y0^2 y3", np.exp(-(y * zs[None, :]).sum(axis=1)) * y[:, 0] ** 2 * y[:, 3]),
                 ("control F = X", Xv)]
        out = []
        for name, F in cases:
            I = (g * F).sum()
            Ia = np.abs(g * F).sum()
            out.append("%s: %.3e / %.3e" % (name, I, Ia))
        print("Z2 ear point %d: int_central g F / int_central |g F| -> %s" % (trial, " ; ".join(out)))
    # convergence of the control with nq
    for nq2 in (40, 80, 120):
        xg, wg = np.polynomial.legendre.leggauss(nq2)
        ph = (xg + 1) * np.pi / 4
        wph = wg * np.pi / 4
        PH, CH = np.meshgrid(ph, ph, indexing="ij")
        WW2 = np.outer(wph, wph).ravel()
        Xv2, Wv2 = np.sin(PH).ravel() ** 2, np.sin(CH).ravel() ** 2
        js2 = (2 * np.sin(PH) * np.cos(PH) * 2 * np.sin(CH) * np.cos(CH)).ravel()
        full = {b: np.full(Xv2.shape, float(v)) for b, v in ep.items()}
        full[Xn], full[Wn] = Xv2, Wv2
        for nd, nm in rc.blocks:
            full[nm[0]] = 1 - sum(full[b] for b in nm[1:])
        t, J = rc.t_and_jac(full)
        D, y = vc.density_batch(nv, edges, C, t, lamf, return_y=True)
        g = D * np.abs(J) * js2 * WW2
        print("   nq = %3d: int g = %.3e ; int g X = %.12e" % (nq2, g.sum(), (g * Xv2).sum()))

    # ---------------- Z3
    f = lambda t, J, full: vc.density_batch(nv, edges, C, t, lamf) * J
    ex = vsp.probe_exponents(rc, f, random.Random(1))
    alpha = {b: max(0.5, ex[b] + 1.25) for b in rc.names}
    zs = np.array([rng.uniform(-0.3, 0.3) for _ in range(N)])
    lamF = [[Fr(float(x)) for x in row] for row in lamf]

    def exact(xi, mode):
        full = rc.complete(xi)
        tt = rc.t_of(full)
        Jx = vc.det_fr([tt] + rc.jac_rows(full, tt))
        coef, U, y = vc.density_slice(net, C, tt, lamF)
        g = float(coef * abs(Jx)) / math.sqrt(float(U))
        if mode == 1:
            g *= math.exp(-sum(float(y[e][0]) * zs[e] for e in range(N)))
        if mode == 2:
            g = abs(g)
        return g

    for mode, name, func in [(0, "z = 0", lambda t, J, full: vc.density_batch(nv, edges, C, t, lamf) * J),
                             (1, "z random", lambda t, J, full: (lambda D, y: D * J * np.exp(-(y[:, :, 0] * zs[None, :]).sum(axis=1)))(
                                 *vc.density_batch(nv, edges, C, t, lamf, return_y=True))),
                             (2, "|g|", lambda t, J, full: np.abs(vc.density_batch(nv, edges, C, t, lamf) * J))]:
        t0 = time.time()
        mean, err, _, nex = vsp.qmc_integrate(rc, func, alpha, m=15, nrep=16, seed=40,
                                              exact_cb=lambda xi, mode=mode: exact(xi, mode))
        print("Z3 %-9s own QMC 2^15 x 16: %.5g +- %.2g ; %d points exact  (%.1fs)" % (name, mean, err, nex, time.time() - t0))
