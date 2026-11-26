"""Independent whole-graph route for one insertion at D = 3, n = 4 (no chain decomposition, no skeleton
current polynomials):  P_e^i = U^(n+H) D * (U y^i_e), degree n^2 - 1 + n, degree <= n per variable,
interpolated mod p on its lower set ([0,4]^12, |e| <= 19: 4.4e7 points), pushed forward along the
series chains (fastker), contracted with I_{n+H+1}(skeleton) (q = 8 tables);  compared mod p with the
assembled exact moment.   usage: python3 check_direct_n4.py IDX E I [E I ...]   (bulk3, standard sources)"""
import sys
import time

import numpy as np

import sets
import assembly
import fastsample
from numerator import BareNumerator, LowerSet
from corelookup import CoreLookup
from modlin import fmod
from run_batch import std_sources

P = 2147483629


def main(argv):
    idx = int(argv[0])
    pairs = [(int(argv[k]), int(argv[k + 1])) for k in range(1, len(argv), 2)]
    G, O, S0 = dict(sets.items('bulk3'))[idx]
    Gh, lab, lam, full = sets.instance('bulk3', idx, G, O, S0, 3)
    legs, fstd = std_sources(Gh, lab, S0, O, 3)
    net = assembly.Network(Gh.nv, Gh.edges, fstd)
    t0 = time.time()
    res = assembly.Assembly(net, order=1, verbose=False).run()
    print('assembled moments (exact) %.0fs' % (time.time() - t0), flush=True)
    num = BareNumerator(net.G, net.lam)
    n, N, H = net.n, net.N, net.H
    total = n * n - 1 + n
    t0 = time.time()
    ls = LowerSet(N - 1, n, total)
    print('lower set %d points (%.0fs)' % (ls.M, time.time() - t0), flush=True)
    groups = [c['edges'] for c in net.chains]
    cl = CoreLookup(net.skel_edges, n + H + 1)
    for (e, i) in pairs:
        t0 = time.time()
        vals = np.empty(ls.M, dtype=np.int64)
        ch = 1000000
        for s in range(0, ls.M, ch):
            d = ls.dig[s:s + ch].astype(np.int64)
            tp = np.concatenate([np.ones((len(d), 1), dtype=np.int64), 1 + d], axis=1)
            vals[s:s + ch] = fastsample.samples_ins(tp, num.C, num.y0, n, H, P, np.array([e]), np.array([i]))
        t1 = time.time()
        coef = ls.newton_to_monomial(vals, P)
        del vals
        t2 = time.time()
        # push-forward (chunks): core key -> sum coef * prod Gamma(a_e) / Gamma(a_c)
        fact = [1]
        for k in range(1, 4 * total + 8):
            fact.append(fact[-1] * k % P)
        fa = np.array(fact, dtype=np.int64)
        ifa = np.array([pow(int(x), P - 2, P) for x in fa], dtype=np.int64)
        base = 4 * total + 4
        keys_all, vals_all = [], []
        for s in range(0, ls.M, ch):
            d = ls.dig[s:s + ch].astype(np.int64)
            a = np.concatenate([(total - d.sum(axis=1))[:, None], d], axis=1) + 1
            w = coef[s:s + ch] % P
            core = np.stack([a[:, g].sum(axis=1) for g in groups], axis=1)
            for ci, g in enumerate(groups):
                if len(g) > 1:
                    for ee in g:
                        w = w * fa[a[:, ee] - 1] % P
                    w = w * ifa[core[:, ci] - 1] % P
            key = np.zeros(len(core), dtype=np.int64)
            for j in range(core.shape[1]):
                key = key * base + core[:, j]
            nz = w != 0
            uk, inv = np.unique(key[nz], return_inverse=True)
            acc = np.zeros(len(uk), dtype=np.int64)
            np.add.at(acc, inv.ravel(), w[nz])
            keys_all.append(uk)
            vals_all.append(acc % P)
        uk, inv = np.unique(np.concatenate(keys_all), return_inverse=True)
        acc = np.zeros(len(uk), dtype=np.int64)
        np.add.at(acc, inv.ravel(), np.concatenate(vals_all))
        acc %= P
        core = np.zeros((len(uk), len(groups)), dtype=np.int64)
        x = uk.copy()
        for j in range(len(groups) - 1, -1, -1):
            core[:, j] = x % base
            x //= base
        prep = cl.prepare(core)
        V, ok = cl.evaluate(core, prep, P)
        if (acc[~ok] != 0).any():
            print('  divergent core period with nonzero coefficient', flush=True)
        direct = {lb: int((acc * v % P).sum() % P) for lb, v in V.items()}
        ref = {lb: fmod(x, P) for lb, x in res.get(('y', e, i), {}).items()}
        labs = set(direct) | set(ref)
        same = all(direct.get(lb, 0) == ref.get(lb, 0) for lb in labs)
        print('bulk3 #%d (std sources) m[%d][%d]: direct whole-graph route mod p %s, assembled %s: %s  '
              '(samples %.0fs, Newton %.0fs, push-forward + periods %.0fs; %d core monomials)' % (
                  idx, e, i, direct, ref, 'EQUAL' if same else 'DIFFERENT', t1 - t0, t2 - t1, time.time() - t2, len(uk)), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
