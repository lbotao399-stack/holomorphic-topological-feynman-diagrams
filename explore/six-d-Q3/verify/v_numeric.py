"""C5 numerical cross-check: float integral of the brute-force density of S2 eq-kernel-def
(own vectorized nodal Kirchhoff solver, own chart t(r) with analytic forward derivatives),
integrated over D_Gamma with scrambled Sobol points (sorted coordinates on the simplex blocks),
compared with the exact I_Gamma(lambda;0) of v_integrals.py at the explorer's sources."""
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from scipy.stats import qmc

import vk

NPTS = int(sys.argv[1]) if len(sys.argv) > 1 else 2 ** 16
NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 8
EXACT = {'G1': Fr(52924313, 311040), 'G2': Fr(-85183, 25920), 'G3': Fr(58733, 25920),
         'G4': Fr(7757, 1440), 'G5': Fr(-18409, 51840)}
LAM = {
    'G1': [[0, 2, -1], [-1, 3, 3], [3, 4, 4], [3, 2, 0], [-4, 3, -4]],
    'G2': [[3, 2, -1], [1, 4, 4], [3, 1, -4], [-4, -2, 4], [2, 0, 3]],
    'G3': [[-1, -1, 4], [-1, -4, -3], [-4, -4, 2], [3, 3, -1], [-4, 3, 3]],
    'G4': [[3, -4, 1], [0, -1, -3], [-1, -3, -1], [0, -1, -3], [1, -3, -2]],
    'G5': [[-1, -2, -2], [1, 3, -1], [-4, -1, 0], [1, -1, -4], [0, 4, 2]],
}


class VD:
    """vectorized dual: v (N,), g (N, k)"""

    def __init__(self, v, g):
        self.v, self.g = v, g

    def __add__(self, o):
        o = self._l(o)
        return VD(self.v + o.v, self.g + o.g)
    __radd__ = __add__

    def __sub__(self, o):
        o = self._l(o)
        return VD(self.v - o.v, self.g - o.g)

    def __rsub__(self, o):
        return self._l(o) - self

    def __mul__(self, o):
        o = self._l(o)
        return VD(self.v * o.v, self.v[:, None] * o.g + o.v[:, None] * self.g)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._l(o)
        return VD(self.v / o.v, (self.g * o.v[:, None] - self.v[:, None] * o.g) / (o.v ** 2)[:, None])

    def _l(self, x):
        if isinstance(x, VD):
            return x
        return VD(np.full_like(self.v, float(x)), np.zeros_like(self.g))


def chart_vec(name, r):
    N, k = r.shape
    R = [VD(r[:, i].copy(), np.tile(np.eye(k)[i], (N, 1))) for i in range(k)]
    recs = vk.build_moves(name)
    tval = {}
    kind, mv, roles, _ = recs[-1]
    b = 3 * (len(recs) - 1)
    a1, a2, a3 = R[b], R[b + 1], R[b + 2]
    tval[roles['t1']] = 1 - a1
    tval[roles['t2']] = a1 - a2
    tval[roles['t3']] = a2 - a3
    tval[roles['t4']] = a3
    for mi in range(len(recs) - 2, -1, -1):
        kind, mv, roles, mk = recs[mi]
        a1, a2, a3 = R[3 * mi], R[3 * mi + 1], R[3 * mi + 2]
        tr = tval.pop(mk)
        if kind == 'chain':
            T = tr / (1 - a1)
            tval[roles['t0']] = tr / a1
            tval[roles['t1']] = T * (1 - a2)
            tval[roles['t2']] = T * (a2 - a3)
            tval[roles['t3']] = T * a3
        else:
            T2 = tr / a1
            T1 = tr / (1 - a1)
            tval[roles['vq']] = a2 * T1
            tval[roles['pv']] = (1 - a2) * T1
            tval[roles['v2q']] = a3 * T2
            tval[roles['pv2']] = (1 - a3) * T2
    m = len(vk.GRAPHS[name]['edges'])
    ts = [tval[('e', i)] for i in range(m)]
    return np.stack([x.v for x in ts], 1), np.stack([x.g for x in ts], 1)   # (N,m), (N,m,k)


def density_vec(name, r, lam, tree):
    g = vk.GRAPHS[name]
    net = vk.Network(g['edges'], g['verts'])
    B = np.array(net.B, dtype=float)       # m x nv
    t, dt = chart_vec(name, r)
    N, m = t.shape
    nv = B.shape[1]
    L = np.einsum('ei,ne,ej->nij', B, 1.0 / t, B)
    lamv = np.array([[lam.get(v, (0, 0, 0))[a] for a in range(3)] for v in g['verts']], dtype=float)
    rhs = np.concatenate([np.broadcast_to(lamv, (N, nv, 3)), np.broadcast_to(B.T, (N, nv, m))], axis=2)
    sol = np.linalg.solve(L, rhs)
    phi, X = sol[:, :, :3], sol[:, :, 3:]
    y = np.einsum('ei,nia->nea', B, phi) / t[:, :, None]
    bphi = np.einsum('ji,nia->nja', B, phi)                        # (N, m_j, 3)
    # D[n,e,a,j] = (1/t_e) sum_i B[e,i] X[i,j] bphi[j,a] / t_j^2 - delta_ej y_ea / t_e
    BX = np.einsum('ei,nij->nej', B, X)
    D = BX[:, :, None, :] * (bphi / (t ** 2)[:, :, None]).transpose(0, 2, 1)[:, None, :, :] / t[:, :, None, None]
    idx = np.arange(m)
    D[:, idx, :, idx] -= (y / t[:, :, None]).transpose(1, 0, 2)
    dydr = np.einsum('neaj,njk->neak', D, dt)
    sT, cot = net.tree_sign(tree)
    rows = dydr[:, cot, :, :].reshape(N, 3 * len(cot), -1)
    jac = np.linalg.det(rows)
    Mo = np.concatenate([t[:, None, :], dt.transpose(0, 2, 1)], axis=1)
    orient = np.sign(np.linalg.det(Mo))
    return float(sT) * jac * orient


def map_points(name, u):
    cols = []
    vol = 1.0
    for m, mv in enumerate(vk.GRAPHS[name]['moves']):
        a = u[:, 3 * m:3 * m + 3]
        if mv[0] == 'chain':
            s = np.sort(a[:, 1:], axis=1)[:, ::-1]
            cols += [a[:, :1], s]
            vol *= 0.5
        elif mv[0] == 'bubble':
            cols += [a]
        else:
            cols += [np.sort(a, axis=1)[:, ::-1]]
            vol /= 6.0
    return np.concatenate(cols, axis=1), vol


if __name__ == '__main__':
    names = sys.argv[3:] or ['G1', 'G2', 'G3', 'G4', 'G5']
    for name in names:
        t0 = time.time()
        g = vk.GRAPHS[name]
        lam = {v: tuple(x) for v, x in zip(g['outputs'], LAM[name])}
        net = vk.Network(g['edges'], g['verts'])
        tree = net.random_tree(__import__('random').Random(5))
        ests = []
        for rep in range(NREP):
            u = qmc.Sobol(9, scramble=True, seed=1000 + rep).random(NPTS)
            u = np.clip(u, 1e-12, 1 - 1e-12)
            r, vol = map_points(name, u)
            vals = []
            for i0 in range(0, NPTS, 8192):
                vals.append(density_vec(name, r[i0:i0 + 8192], lam, tree))
            ests.append(vol * np.concatenate(vals).mean())
        ests = np.array(ests)
        mean, err = ests.mean(), ests.std(ddof=1) / np.sqrt(NREP)
        ex = float(EXACT[name])
        rel = abs(mean - ex) / abs(ex)
        print('%s: exact I(lambda;0) = %.10g; Sobol %d x %d of the brute-force density (cotree %s): %.10g +- %.2g; '
              'rel. dev. %.2e (%.1f digits), deviation/sigma %.2f  [%.0fs]'
              % (name, ex, NREP, NPTS, net.tree_sign(tree)[1], mean, err, rel, -np.log10(rel), abs(mean - ex) / err,
                 time.time() - t0), flush=True)
