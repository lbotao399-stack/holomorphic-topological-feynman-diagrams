"""Part (c): two-loop obstruction of topological-holomorphic Chern-Simons on R x C^2 as a
Chevalley-Eilenberg 4-cochain of g[[z_1, z_2]] modulo translations, its closedness and exactness.

omega = sum_Gamma |Aut Gamma|^{-1} alt( K_Gamma(p_1..p_4) W_Gamma ),  Gamma = theta(1,3,3), theta(2,2,3),
K_Gamma from the engine (checked in c_kernels_T1.py), W_Gamma by the orientation rule of cochain.py
(same vertex labels and edge orientations as the engine graph).  Common factors dropped:
K_133 = -(pi/20) [p1 p2][p3 p4],  K_223 = (pi/30) [p1 p2][p3 p4],  |Aut| = 4 for both, so
    omega  propto  omega_133 + omega_223,   omega_G := alt((60/pi) K_G W_G).
Closedness: d omega = 0 (canonical form modulo conservation, exact integers).
Exactness: omega in d(C^3) at weight (2,2): spanning set alt(T (x) mu) of C^3, T over invariant
3-tensors (trace monomials, f), mu over the 9 canonical monomials of degree (2,2); ranks exact.
For sl2 also with all 27 elementary 3-tensors (no invariance assumed).
"""
import sys, time, itertools
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import htk
from htk import theta
import lie, cochain as cc

H = 2
mom = cc.Momenta(6, H)
P = mom.p


def kernel_poly(pqr):
    G, interior = theta(*pqr)
    ext = [v for Pp in interior for v in Pp]
    lam = {v: (list(P[ext.index(v)]) if v in ext else [0, 0]) for v in range(G.nv - 1)}
    K = sp.expand(htk.kernel_exact(G, lam, T=1) * 60 / sp.pi)
    return G, ext, K


def trace_tensors(g, k):
    """all products of traces tr(X_{a_s(1)} ...) over permutations s of S_k (cycle decomposition)."""
    mats = g.mats
    d = g.dim
    out = []
    for s in itertools.permutations(range(k)):
        # cycles of s
        seen, cycles = set(), []
        for i in range(k):
            if i in seen:
                continue
            c, j = [], i
            while j not in seen:
                seen.add(j); c.append(j); j = s[j]
            cycles.append(c)
        T = np.zeros((d,) * k, dtype=np.int64)
        for idx in itertools.product(range(d), repeat=k):
            val = 1
            for c in cycles:
                m = np.eye(mats[0].shape[0], dtype=np.int64)
                for j in c:
                    m = m @ mats[idx[j]]
                val *= int(np.trace(m))
                if val == 0:
                    break
            T[idx] = val
        out.append((s, T))
    return out


def run(g, full_basis=False):
    t0 = time.time()
    res = {}
    omegas = {}
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        G, ext, K = kernel_poly(pqr)
        W, sc = cc.colour_weight(g, G.nv, G.edges, ext)
        omegas[pqr] = cc.alternate(mom, K, W, 4)
        res['W%d%d%d_scale' % pqr] = sc
    w133, w223 = omegas[(1, 3, 3)], omegas[(2, 2, 3)]
    omega = {m: w133.get(m, 0) + w223.get(m, 0) for m in set(w133) | set(w223)}
    omega = {m: T for m, T in omega.items() if np.any(T)}
    res['omega133_zero'] = cc.is_zero(w133)
    res['omega223_zero'] = cc.is_zero(w223)
    res['omega_zero'] = cc.is_zero(omega)
    d133 = cc.d_CE(mom, g, w133, 4)
    d223 = cc.d_CE(mom, g, w223, 4)
    res['d_omega133_zero'] = cc.is_zero(d133)
    res['d_omega223_zero'] = cc.is_zero(d223)
    dsum = {m: d133.get(m, 0) + d223.get(m, 0) for m in set(d133) | set(d223)}
    res['d_omega_zero'] = cc.is_zero(dsum)
    monos5 = sorted(set(d133) | set(d223))
    if monos5:
        A = np.stack([cc.flatten(d133, monos5, g.dim, 5), cc.flatten(d223, monos5, g.dim, 5)], axis=1)
        res['rank[d omega133, d omega223]'] = cc.exact_rank(A)
    # exactness
    monos4 = cc.all_monos(mom, 4, (2, 2))
    monos3 = cc.all_monos(mom, 3, (2, 2))
    if full_basis:
        tens = []
        for idx in itertools.product(range(g.dim), repeat=3):
            T = np.zeros((g.dim,) * 3, dtype=np.int64); T[idx] = 1
            tens.append(T)
    else:
        tens = [T for s, T in trace_tensors(g, 3)] + [g.fl]
    cols = []
    for T in tens:
        for mu in monos3:
            e = mom.mono_expr(mu, 3, [P[0], P[1]])
            eta = cc.alternate(mom, e, T, 3)
            if not eta:
                cols.append(np.zeros(len(monos4) * g.dim ** 4, dtype=np.int64))
                continue
            cols.append(cc.flatten(cc.d_CE(mom, g, eta, 3), monos4, g.dim, 4))
    A = np.stack(cols, axis=1)
    om = cc.flatten(omega, monos4, g.dim, 4)
    r1 = cc.exact_rank(A)
    r2 = cc.exact_rank(np.concatenate([A, om[:, None]], axis=1))
    res['n_eta_spanning'] = len(cols)
    res['rank d(C3)'] = r1
    res['rank [d(C3) | omega]'] = r2
    res['exact'] = (r1 == r2)
    # same for omega133, omega223 separately
    for name, w in (('133', w133), ('223', w223)):
        v = cc.flatten(w, monos4, g.dim, 4)
        res['exact_' + name] = cc.exact_rank(np.concatenate([A, v[:, None]], axis=1)) == r1
    res['time'] = round(time.time() - t0, 1)
    return res, omegas


if __name__ == "__main__":
    for g, full in [(lie.sl(2), False), (lie.sl(2), True), (lie.gl(2), False), (lie.sl(3), False), (lie.gl(3), False)]:
        res, _ = run(g, full)
        print(g.name, "(all 27 elementary 3-tensors)" if full else "(invariant 3-tensors)", res, flush=True)
