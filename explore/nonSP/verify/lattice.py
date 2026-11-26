"""Trapezoidal rule in log-resistance space for K4 subdivisions (direct, no reduction formulas).

K = int_{R^5} F(s) d^5 s,  F(s) = prod_E R_E * int dx D(t(R,x)) prod_E R_E^{a_E-1},  R = (e^{s_0},...,e^{s_4}, 1).
Lattice h Z^5 + shift; flood fill from s = 0 through points with max_S |F_S| > tol * running max.
x-integral exact (Gauss product rule, m = 2 verified exact against m = 3, 4, 6 by xdeg.py).
"""
import time
import numpy as np
from hepp import K4Int


def lattice(I, h, tol, shift=0.0, chunk=3000, maxpts=20_000_000):
    seen = set()
    origin = (0,) * 5
    front = [origin]; seen.add(origin)
    total = 0.0; vmax = 0.0; n = 0; below = 0.0
    while front:
        nxt = []
        for i in range(0, len(front), chunk):
            blk = front[i:i + chunk]
            s = np.array(blk, dtype=float) * h + shift
            R = np.ones((len(blk), 6)); R[:, :5] = np.exp(s)
            F = I.f_R(R) * np.prod(R, axis=1)[:, None]          # (P, S)
            n += len(blk)
            total = total + F.sum(0)
            a = np.abs(F).max(1)
            vmax = max(vmax, a.max())
            for k, av in zip(blk, a):
                if av > tol * vmax:
                    for d in range(5):
                        for sg in (-1, 1):
                            kk = list(k); kk[d] += sg; kk = tuple(kk)
                            if kk not in seen:
                                seen.add(kk); nxt.append(kk)
                else:
                    below += av
        front = nxt
        if n > maxpts:
            raise RuntimeError("too many points")
    return total * h ** 5, below * h ** 5, n


def run(nv, E, lams, h, tol, label=""):
    I = K4Int(nv, E, lams, m=2)
    t0 = time.time()
    v0, b0, n0 = lattice(I, h, tol, 0.0)
    v1, b1, n1 = lattice(I, h, tol, h / 2)
    return v0, v1, b0, n0 + n1, time.time() - t0
