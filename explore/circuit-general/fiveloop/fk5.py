"""Bare T=0 kernels via rank-constrained sparse modular Newton interpolation with the
five-loop period function per5.period (2-sums, K5-e = dual prism, W_n, prism/K33/five-loop tables).
  K = int P / U^q,  P = U^(n+H) D  (degree n^2-1, <= n-1 per time),  q = n + H.
Monomials pushed forward along series paths to the series skeleton (p unchanged), contracted with
exact periods reduced mod p, CRT + rational reconstruction per label; stops when two successive
reconstructions agree (>= min_primes primes)."""
import os
import sys
import time
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CORES = os.path.abspath(os.path.join(HERE, '..', 'cores'))
sys.path.insert(0, CORES)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..', 'code')))

import per5                                     # noqa: E402
import fastker                                  # noqa: E402  (FastKernel, grid_points, interpolate_mod, PRIMES)
import modsolve                                 # noqa: E402
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', 'assembly')))
from numerator import BareNumerator             # noqa: E402


def kernel(nv, edges, lam, primes=fastker.PRIMES, verbose=True, min_primes=2, batch=20000):
    fk = fastker.FastKernel(nv, edges, lam)
    n, N, H = fk.n, fk.N, fk.H
    q = n + H
    t0 = time.time()
    numerator = BareNumerator(fk.G, fk.lam_nr)
    groups, core_edges = fk.paths
    alpha = numerator.exponents() + 1
    if not len(alpha):
        return {}
    core_alpha = np.stack([alpha[:, g].sum(axis=1) for g in groups], axis=1).astype(np.int16)
    if verbose:
        print('  rank-constrained grid %d points, core edges %d  %.1fs' % (
            numerator.lowerset.M, len(core_edges), time.time() - t0), flush=True)
    results = []
    rec = None
    pcache = {}
    for ip, p in enumerate(primes):
        tp = time.time()
        coef = numerator.coefficients_mod(p, chunk=batch)
        ts = time.time()
        fact = [1]
        for i in range(1, int(alpha.max()) * N + 2):
            fact.append(fact[-1] * i % p)
        factn = np.array(fact, dtype=np.int64)
        sf = np.ones(len(alpha), dtype=np.int64)
        for g in groups:
            if len(g) > 1:
                for e in g:
                    sf = sf * factn[alpha[:, e].astype(np.int64) - 1] % p
                ap = alpha[:, g].sum(axis=1).astype(np.int64)
                sf = sf * fastker.modinv_vec(factn[ap - 1], p) % p
        w = coef * sf % p
        nzm = w != 0
        # Exact row keys avoid overflow of the old base**number_of_paths
        # encoding, which can silently identify different higher-loop weights.
        uk, inv = np.unique(core_alpha[nzm], axis=0, return_inverse=True)
        sums = np.zeros(len(uk), dtype=np.int64)
        np.add.at(sums, inv, w[nzm])
        sums %= p
        tot = {}
        for j in range(len(uk)):
            if not sums[j]:
                continue
            ca = tuple(int(x) for x in uk[j])
            val = pcache.get(ca)
            if val is None:
                val = per5.period(core_edges, list(ca), q)
                pcache[ca] = val
            for lab, x in val.items():
                tot[lab] = (tot.get(lab, 0) + int(sums[j]) * (x.numerator % p) * pow(x.denominator % p, p - 2, p)) % p
        results.append((p, tot))
        if verbose:
            print('  prime %d: numerator %.1fs, %d monomials, %d core vectors (%d periods cached), total %.1fs'
                  % (ip + 1, ts - tp, int(nzm.sum()), len(uk), len(pcache), time.time() - tp), flush=True)
        labels = sorted({l for _, t in results for l in t})
        M = 1
        accv = {l: 0 for l in labels}
        for pp, t in results:
            for l in labels:
                a = accv[l]
                x = t.get(l, 0)
                accv[l] = a + M * ((x - a) * pow(M, -1, pp) % pp)
            M *= pp
        new = {}
        okr = True
        for l in labels:
            r = modsolve.ratrec(accv[l], M)
            if r is None:
                okr = False
                break
            if r:
                new[l] = r
        if okr and rec == new and ip + 1 >= min_primes:
            if verbose:
                print('  stable after %d primes, total %.1fs' % (ip + 1, time.time() - t0), flush=True)
            return new
        rec = new if okr else None
    raise ArithmeticError('reconstruction did not stabilize')


def core_numerator(nv, edges, lam, primes=fastker.PRIMES, verbose=True, min_primes=2, batch=20000):
    """exact series-skeleton numerator: K = sum_alpha c_alpha I_q(core_edges; alpha), returns
    (core_edges, q, {alpha: Fraction}); c_alpha reconstructed by CRT until two successive primes agree."""
    fk = fastker.FastKernel(nv, edges, lam)
    n, N, H = fk.n, fk.N, fk.H
    q = n + H
    t0 = time.time()
    numerator = BareNumerator(fk.G, fk.lam_nr)
    groups, core_edges = fk.paths
    alpha = numerator.exponents() + 1
    if not len(alpha):
        return core_edges, q, {}
    core_alpha = np.stack([alpha[:, g].sum(axis=1) for g in groups], axis=1).astype(np.int16)
    uk, inv = np.unique(core_alpha, axis=0, return_inverse=True)
    acc, M, prev = None, 1, None
    for ip, p in enumerate(primes):
        tp = time.time()
        coef = numerator.coefficients_mod(p, chunk=batch)
        fact = [1]
        for i in range(1, int(alpha.max()) * N + 2):
            fact.append(fact[-1] * i % p)
        factn = np.array(fact, dtype=np.int64)
        sf = np.ones(len(alpha), dtype=np.int64)
        for g in groups:
            if len(g) > 1:
                for e in g:
                    sf = sf * factn[alpha[:, e].astype(np.int64) - 1] % p
                ap = alpha[:, g].sum(axis=1).astype(np.int64)
                sf = sf * fastker.modinv_vec(factn[ap - 1], p) % p
        w = coef * sf % p
        sums = np.zeros(len(uk), dtype=np.int64)
        np.add.at(sums, inv, w)
        sums %= p
        if acc is None:
            vecs = [tuple(int(x) for x in row) for row in uk]
            acc = [int(x) for x in sums]
            M = p
        else:
            minv = pow(M, -1, p)
            acc = [a + M * ((int(s) - a) * minv % p) for a, s in zip(acc, sums)]
            M *= p
        rec = {}
        ok = True
        for v, a in zip(vecs, acc):
            r = modsolve.ratrec(a, M)
            if r is None:
                ok = False
                break
            if r:
                rec[v] = r
        if verbose:
            print('  prime %d: %d core vectors, %d nonzero, reconstructed %s  %.1fs' % (
                ip + 1, len(uk), len(rec) if ok else -1, ok, time.time() - tp), flush=True)
        if ok and prev == rec and ip + 1 >= min_primes:
            return core_edges, q, rec
        prev = rec if ok else None
    raise ArithmeticError('core numerator did not stabilize')


def kernel_from_core(core_edges, q, num, period=None):
    period = period or per5.period
    tot = {}
    for a, c in num.items():
        tot = per5.vadd(tot, period(core_edges, list(a), q), c)
    return tot
