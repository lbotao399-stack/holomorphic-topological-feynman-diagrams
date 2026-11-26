"""Float Monte Carlo (independent oracle, not exact) for T = 1, H = 1 kernels in a composite chart:
a chain of star-triangle moves (htk.star.SingleMove / PairedMove, each acting on the reduced graph
of the previous one) followed by the series-parallel chart (htk.spred) of the last reduced graph.

    K = int D_Gamma(t(chi)) |det[t; dt/dchi]| dchi,  |det| = |J_red| prod_moves |J_loc|,
    single: |J_loc| = T_uw T_ux^2 T_wx^2 / (r1^2 r2^2 (1-r1)^3 (1-r2)^3 P^3),
    paired: |J_loc| = T_uA T_uB T_AB^2 / (r1^2 (1-r1)^2 r2^2 (1-r2)^2 r3^2 r4^2 (1-r3-r4)^2) prod_i R^i_uA R^i_uB R^i_AB / P_i^3,
sampling f ~ Beta(1/2,1/2), (r3, r4, 1-r3-r4) ~ Dirichlet(1/2,1/2,1/2), Dirichlet(1/2) on spred simplices.
D_Gamma: htk.oracle.density_batch (float brute force, no reduction)."""
import math
import sys
import time

import numpy as np

from htk import star, catalog, sp_reduce, oracle


def find_chain(G, depth=3):
    """Shortest chain of moves ending in a series-parallel graph: list of move objects."""
    if catalog.is_series_parallel(G):
        return []
    if depth == 0:
        return None
    best = None
    cands = []
    for (v, u, w, x) in star.find_single(G):
        cands.append(star.SingleMove(G, v))
    for (v1, v2, u, A, B) in star.find_paired(G):
        cands.append(star.PairedMove(G, v1, v2))
    for mv in cands:
        rest = find_chain(mv.red, depth - 1)
        if rest is not None and (best is None or len(rest) + 1 < len(best)):
            best = [mv] + rest
            if len(best) == 1:
                return best
    return best


def _single_map(mv, f, T):
    r1, r2 = f[:, 0], f[:, 1]
    Tuw, Tux, Twx = T[:, 0], T[:, 1], T[:, 2]
    Rux, Rwx = Tux / (1 - r1), Twx / (1 - r2)
    P = Tuw + Rux + Rwx
    t = np.zeros((T.shape[0], mv.G.ne))
    t[:, mv.e["vu"]] = Tuw * Rux / P
    t[:, mv.e["vw"]] = Tuw * Rwx / P
    t[:, mv.e["vx"]] = Rux * Rwx / P
    t[:, mv.e["ux"]] = Tux / r1
    t[:, mv.e["wx"]] = Twx / r2
    for j, e in enumerate(mv.others):
        t[:, e] = T[:, 3 + j]
    J = Tuw * Tux ** 2 * Twx ** 2 / (r1 ** 2 * r2 ** 2 * (1 - r1) ** 3 * (1 - r2) ** 3 * P ** 3)
    return t, J


def _paired_map(mv, f, T):
    r1, r2, r3, r4 = f[:, 0], f[:, 1], f[:, 2], f[:, 3]
    TuA, TuB, TAB = T[:, 0], T[:, 1], T[:, 2]
    s = 1 - r3 - r4
    Rs = [(TuA / r1, TuB / r2, TAB / r3), (TuA / (1 - r1), TuB / (1 - r2), TAB / r4)]
    t = np.zeros((T.shape[0], mv.G.ne))
    J = TuA * TuB * TAB ** 2 / (r1 ** 2 * (1 - r1) ** 2 * r2 ** 2 * (1 - r2) ** 2 * r3 ** 2 * r4 ** 2 * s ** 2)
    for i, R in ((1, Rs[0]), (2, Rs[1])):
        P = R[0] + R[1] + R[2]
        t[:, mv.e["%du" % i]] = R[0] * R[1] / P
        t[:, mv.e["%dA" % i]] = R[0] * R[2] / P
        t[:, mv.e["%dB" % i]] = R[1] * R[2] / P
        J = J * R[0] * R[1] * R[2] / P ** 3
    t[:, mv.e["AB"]] = TAB / s
    for j, e in enumerate(mv.others):
        t[:, e] = T[:, 3 + j]
    return t, J


def _sample_fractions(mv, rng, M):
    if isinstance(mv, star.SingleMove):
        f = rng.beta(0.5, 0.5, size=(M, 2))
        logq = np.sum(-0.5 * np.log(f) - 0.5 * np.log1p(-f) - math.log(math.pi), axis=1)
        return f, logq
    a = rng.beta(0.5, 0.5, size=(M, 2))
    d = rng.dirichlet([0.5] * 3, size=M)
    f = np.concatenate([a, d[:, :2]], axis=1)
    logq = np.sum(-0.5 * np.log(a) - 0.5 * np.log1p(-a) - math.log(math.pi), axis=1)
    logq += math.lgamma(1.5) - 3 * math.lgamma(0.5) + np.sum(-0.5 * np.log(d), axis=1)
    return f, logq


def _sample_spred(spr, rng, M):
    logq = np.zeros(M)
    bary = {}
    for blk in spr.blocks:
        m = len(blk.bary)
        x = rng.dirichlet([0.5] * m, size=M)
        for j, b in enumerate(blk.bary):
            bary[b] = x[:, j]
        logq += math.lgamma(0.5 * m) - m * math.lgamma(0.5) + np.sum(-0.5 * np.log(x), axis=1)
    return bary, logq


def _spred_times(spr, red, bary, M):
    Nr = red.ne
    T = np.ones((M, Nr))
    for e in range(Nr):
        for b, k in spr.texp[e].items():
            T[:, e] *= bary[b] ** k
    cols = []
    for c in spr.coords:
        blk = next(bb for bb in spr.blocks if c in bb.bary)
        b0 = blk.bary[0]
        col = np.zeros((M, Nr))
        for e in range(Nr):
            col[:, e] = spr.texp[e].get(c, 0) / bary[c] - spr.texp[e].get(b0, 0) / bary[b0]
        cols.append(col)
    L = np.stack([np.ones((M, Nr))] + cols, axis=1)
    return T, np.abs(np.prod(T, axis=1) * np.linalg.det(L))


def _safe_density(G, t, lam, powers, chunk=5000):
    """oracle.density_batch in chunks; a chunk with a numerically singular loop matrix is evaluated
    row by row and the singular rows (measure-zero corners) are set to 0."""
    out = np.zeros(t.shape[0])
    with np.errstate(all="ignore"):
        rng_ok = np.log(np.max(t, axis=1) / np.min(t, axis=1)) < math.log(1e12)
    for a in range(0, t.shape[0], chunk):
        sl = slice(a, a + chunk)
        try:
            out[sl] = oracle.density_batch(G, t[sl], lam, T=1, powers=powers)
        except np.linalg.LinAlgError:
            for i in range(a, min(a + chunk, t.shape[0])):
                try:
                    out[i] = oracle.density_batch(G, t[i:i + 1], lam, T=1, powers=powers)[0]
                except np.linalg.LinAlgError:
                    out[i] = 0.0
    # samples with time ratios beyond 1e12 (float cancellation in the oracle density) are dropped
    return np.where(rng_ok, out, 0.0)


def estimate(G, lam, Mper, batches, seed=0, chain=None, powers=None):
    chain = find_chain(G) if chain is None else chain
    if chain is None:
        raise ValueError("no chart found")
    last = chain[-1].red if chain else G
    spr = sp_reduce(last, prefer="large")
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(batches):
        bary, logq = _sample_spred(spr, rng, Mper)
        T, J = _spred_times(spr, last, bary, Mper)
        for mv in reversed(chain):
            f, lq = _sample_fractions(mv, rng, Mper)
            logq = logq + lq
            T, Jl = (_single_map if isinstance(mv, star.SingleMove) else _paired_map)(mv, f, T)
            J = J * Jl
        D = _safe_density(G, T, lam, powers)
        w = D * J * np.exp(-logq)
        w = w[np.isfinite(w)]
        vals.append(np.mean(w))
    vals = np.array(vals)
    return float(vals.mean()), float(vals.std(ddof=1) / math.sqrt(len(vals))), chain


def describe(chain):
    out = []
    for mv in chain:
        if isinstance(mv, star.SingleMove):
            out.append("Y(%d;%d,%d,%d)" % (mv.v, mv.u, mv.w, mv.x))
        else:
            out.append("YY(%d,%d;%d,%d,%d)" % (mv.v1, mv.v2, mv.u, mv.A, mv.B))
    return " ".join(out) if out else "SP"
