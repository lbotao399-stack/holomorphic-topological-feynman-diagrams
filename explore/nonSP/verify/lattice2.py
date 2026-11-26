"""Vectorized trapezoidal rule in log-space (same integrand as lattice.py; flood fill done with numpy keys).

K = int_{R^d} F(s) d^d s, F(s) = Fcall(exp(s)); lattice h Z^d + shift; BFS from 0 through points with
max_S |F_S| > tol * max (running max over evaluated layers).  Returns sum over all evaluated points (incl. the
sub-threshold boundary layer), the boundary-layer mass, and the number of points."""
import time
import numpy as np

OFF, BASE = 256, 512


def enc(c):
    k = np.zeros(len(c), dtype=np.int64)
    for d in range(c.shape[1] - 1, -1, -1):
        k = k * BASE + (c[:, d] + OFF)
    return k


def dec(k, dim):
    c = np.zeros((len(k), dim), dtype=np.int64)
    for d in range(dim):
        c[:, d] = k % BASE - OFF
        k = k // BASE
    return c


def lattice(Fcall, dim, h, tol, shift=0.0, chunk=4000, maxpts=30_000_000, verbose=False):
    front = np.zeros((1, dim), dtype=np.int64)
    visited = enc(front)
    total = 0.0; vmax = 0.0; below = 0.0; n = 0; layer = 0
    t0 = time.time()
    while len(front):
        vals = []
        for i in range(0, len(front), chunk):
            s = front[i:i + chunk].astype(float) * h + shift
            vals.append(Fcall(np.exp(s)))
        F = np.concatenate(vals, 0)                 # (P, S)
        n += len(front)
        total = total + F.sum(0)
        a = np.abs(F).max(1)
        vmax = max(vmax, a.max())
        act = a > tol * vmax
        below += a[~act].sum()
        A = front[act]
        if len(A) == 0:
            break
        nb = np.concatenate([A + np.eye(dim, dtype=np.int64)[d] * sg for d in range(dim) for sg in (-1, 1)], 0)
        kn = np.unique(enc(nb))
        kn = kn[~np.isin(kn, visited, assume_unique=True)]
        visited = np.union1d(visited, kn)
        front = dec(kn, dim)
        layer += 1
        if verbose:
            print("  layer %d front %d n %d  %.0fs" % (layer, len(front), n, time.time() - t0), flush=True)
        if n > maxpts:
            raise RuntimeError("too many points")
    return total * h ** dim, below * h ** dim, n
