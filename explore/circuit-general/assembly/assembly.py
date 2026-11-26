"""Kernels with current insertions assembled through the series reduction Gamma -> skeleton.

Chain c (skeleton edge p_c -> q_c), edges e_{c,1..L}, signs s_{c,j} (+1 if oriented p -> q),
interior vertices u_{c,1..L-1}; t_{c,j} = x_{c,j} R_c, sigma_{c,k} = sum_{j>k} x_{c,j}.
  y_{e_{c,j}} = s_{c,j} [ Y_c(R; lambda^s(sigma)) - mhat_c + Q_{c,j} ],
  mhat_c = sum_k sigma_{c,k} lambda_{u_{c,k}},  Q_{c,j} = sum_{k<j} lambda_{u_{c,k}},
  lambda^s_w = lambda_w + sum_{c,k} lambda_{u_{c,k}} ([p_c = w] sigma_{c,k} + [q_c = w](1 - sigma_{c,k})),
  U_s Y_c = sum_w Phi_{c,w}(R) lambda^s_w,   Phi_{c,w} = sum_T flow_T(c; w -> z) prod_{f not in T} R_f
[src htk CONVENTIONS §5 (chain move), atlas REPORT §1.4; chk Network.check_point, check_pointwise.py].
Integrals: K[w(x) F(R)/U_s^k] = sum_alpha c_alpha  prod_c Dir_c(alpha, w)  I_{q0+k}(skeleton; alpha_c + gamma),
P = U^(n+H) D = sum c_alpha t^(alpha-1),  Dir_c = prod_{e in c} Gamma(alpha_e + beta_e) / Gamma(alpha_c + |beta_c|)
for w = x^beta (series push-forward, cores/fastker.py), q0 = n + H.
"""
import os
import sys
import time
from fractions import Fraction as F
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..', 'code')))

from htk import brute                       # noqa: E402
from htk.graph import Graph                 # noqa: E402

import modlin                               # noqa: E402
from modlin import CRT, PRIMES              # noqa: E402
from numerator import BareNumerator         # noqa: E402
from corelookup import CoreLookup           # noqa: E402


class Network:
    """htk graph with integer sources; chains and series skeleton."""

    def __init__(self, nv, edges, lam):
        self.G = G = Graph(nv, [tuple(e) for e in edges])
        lam_nr = brute.sources(G, lam)
        self.H = H = len(lam_nr[0])
        for v in lam_nr:
            for x in v:
                if F(x).denominator != 1:
                    raise ValueError('integer sources required')
        lam_full = [[int(x) for x in v] for v in lam_nr]
        lam_full.append([-sum(v[i] for v in lam_full) for i in range(H)])
        self.lam = lam_full
        self.n, self.N = G.loops, G.ne
        deg = [0] * nv
        inc = [[] for _ in range(nv)]
        for e, (a, b) in enumerate(G.edges):
            deg[a] += 1
            deg[b] += 1
            inc[a].append(e)
            inc[b].append(e)
        self.deg = deg
        skel = [v for v in range(nv) if deg[v] != 2]
        if len(skel) < 2:
            raise ValueError('no skeleton')
        used = set()
        chains = []
        for s in skel:
            for e0 in inc[s]:
                if e0 in used:
                    continue
                es, sg, interior = [], [], []
                prev, e = s, e0
                while True:
                    a, b = G.edges[e]
                    nxt = b if a == prev else a
                    es.append(e)
                    sg.append(1 if (a, b) == (prev, nxt) else -1)
                    used.add(e)
                    if deg[nxt] != 2:
                        break
                    interior.append(nxt)
                    e = [f for f in inc[nxt] if f != e][0]
                    prev = nxt
                if nxt == s:
                    raise ValueError('chain from %d returns to itself' % s)
                chains.append(dict(p=s, q=nxt, edges=es, signs=sg, interior=interior))
        self.chains = chains
        self.skel = skel
        self.sidx = {v: i for i, v in enumerate(skel)}
        self.skel_edges = [(self.sidx[c['p']], self.sidx[c['q']]) for c in chains]
        self.Ns = len(chains)
        self.Vs = len(skel)
        # position of every original edge: (chain, j) 0-based j
        self.where = {}
        for ci, c in enumerate(chains):
            for j, e in enumerate(c['edges']):
                self.where[e] = (ci, j)
        # x-weights: index 0 = constant 1; one per (chain, j >= 1) with j 0-based  (x_{c,j+1})
        self.xw = []
        for ci, c in enumerate(chains):
            for j in range(1, len(c['edges'])):
                self.xw.append((ci, j))
        self.xidx = {cj: 1 + k for k, cj in enumerate(self.xw)}
        self.nw = 1 + len(self.xw)
        self._trees()

    # ---------------------------------------------------------------- skeleton trees
    def _trees(self):
        Vs, Es = self.Vs, self.skel_edges
        z = Vs - 1                      # sink of the skeleton flows
        trees = []
        for T in combinations(range(len(Es)), Vs - 1):
            parent = list(range(Vs))

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            ok = True
            for e in T:
                a, b = find(Es[e][0]), find(Es[e][1])
                if a == b:
                    ok = False
                    break
                parent[a] = b
            if ok:
                trees.append(T)
        self.trees = trees
        nT = len(trees)
        # flow[T][c][w]: flow on skeleton edge c (p->q) for unit source at w, sink z
        flow = np.zeros((nT, len(Es), Vs), dtype=np.int64)
        cot = np.ones((nT, len(Es)), dtype=np.int64)
        for ti, T in enumerate(trees):
            adj = {v: [] for v in range(Vs)}
            for e in T:
                a, b = Es[e]
                adj[a].append((b, e))
                adj[b].append((a, e))
                cot[ti, e] = 0
            par = {z: (None, None)}
            order = [z]
            for v in order:
                for w, e in adj[v]:
                    if w not in par:
                        par[w] = (v, e)
                        order.append(w)
            for w in range(Vs):
                if w == z:
                    continue
                v = w
                while v != z:
                    pv, e = par[v]
                    flow[ti, e, w] += 1 if Es[e][0] == v else -1
                    v = pv
        self.flow = flow          # T x c x w
        self.cot = cot            # T x c  (0/1 cotree indicator)

    # ---------------------------------------------------------------- source structure
    def sigma_vec(self, ci, k):
        """sigma_{c,k} (k 1-based interior index) as a vector over x-weights"""
        v = np.zeros(self.nw, dtype=object)
        L = len(self.chains[ci]['edges'])
        for j in range(k, L):          # 0-based edges j >= k  <->  x_{c,j+1}, j+1 > k
            v[self.xidx[(ci, j)]] += 1
        return v

    def skeleton_source_weights(self, i):
        """Lam[w, beta]: coefficient of x-weight beta in lambda^{s,i}_w(sigma)  (Vs x nw)"""
        Lam = np.zeros((self.Vs, self.nw), dtype=object)
        for v in self.skel:
            Lam[self.sidx[v], 0] += self.lam[v][i]
        for ci, c in enumerate(self.chains):
            pw, qw = self.sidx[c['p']], self.sidx[c['q']]
            for k, u in enumerate(c['interior'], start=1):
                lu = self.lam[u][i]
                if not lu:
                    continue
                sv = self.sigma_vec(ci, k)
                Lam[pw] += lu * sv
                Lam[qw, 0] += lu
                Lam[qw] -= lu * sv
        return Lam

    def insertion_coefficients(self, e, i):
        """y^i_e = s [ sum_{T,beta} A1[T,beta] (x^beta Phi-term of tree T)/U_s
                       + sum_beta A0[beta] x^beta ]  ->  (s, A1 (nT x nw), A0 (nw))"""
        ci, j = self.where[e]
        c = self.chains[ci]
        s = c['signs'][j]
        Lam = self.skeleton_source_weights(i)                 # Vs x nw
        A1 = self.flow[:, ci, :].astype(object) @ Lam         # nT x nw
        A0 = np.zeros(self.nw, dtype=object)
        for k, u in enumerate(c['interior'], start=1):
            lu = self.lam[u][i]
            if not lu:
                continue
            A0 -= lu * self.sigma_vec(ci, k)                  # - mhat
            if k <= j:                                        # Q_{c,j}: interior vertices before edge j
                A0[0] += lu
        return s, A1, A0

    def check_point(self, seed=1):
        """exact check of the chain decomposition of the currents at a random rational point"""
        import random
        rng = random.Random(seed)
        G = self.G
        t = [F(rng.randint(1, 30), rng.randint(1, 7)) for _ in range(self.N)]
        y = brute.currents(G, t, self.lam[:-1])
        R = [sum(t[e] for e in c['edges']) for c in self.chains]
        x = {}
        for ci, c in enumerate(self.chains):
            for j, e in enumerate(c['edges']):
                x[(ci, j)] = t[e] / R[ci]
        xv = [F(1)] + [x[cj] for cj in self.xw]
        cotm = [np.prod([R[f] for f in range(self.Ns) if self.cot[ti, f]]) for ti in range(len(self.trees))]
        U = sum(cotm)
        for e in range(self.N):
            for i in range(self.H):
                s, A1, A0 = self.insertion_coefficients(e, i)
                val = sum(A1[ti, b] * xv[b] * cotm[ti] for ti in range(len(self.trees)) for b in range(self.nw)) / U
                val += sum(A0[b] * xv[b] for b in range(self.nw))
                if s * val != y[e][i]:
                    raise AssertionError('current decomposition fails at edge %d comp %d' % (e, i))
        return True


class Assembly:
    """K(lambda; 0), first moments int D y^i_e and (order 2) second moments int D y^i_e y^j_f,
    modulo primes, assembled on the series skeleton."""

    def __init__(self, net, order=1, verbose=True, pairs=None):
        self.net = net
        self.order = order
        self.verbose = verbose
        self.num = BareNumerator(net.G, net.lam)
        self.q0 = net.n + net.H
        t0 = time.time()
        ex = self.num.exponents()
        alpha = ex.astype(np.int16) + 1                               # t^(alpha-1)
        self.alpha = alpha
        groups = [c['edges'] for c in net.chains]
        core = np.stack([alpha[:, g].sum(axis=1) for g in groups], axis=1).astype(np.int32)
        base = int(core.max()) + 2
        key = np.zeros(len(core), dtype=np.int64)
        for j in range(core.shape[1]):
            key = key * base + core[:, j]
        uk, first, inv = np.unique(key, return_index=True, return_inverse=True)
        del key, uk
        self.ucore = core[first].astype(np.int64)                     # Mc x Ns
        inv = inv.ravel()
        self.order_idx = np.argsort(inv, kind='stable')
        sinv = inv[self.order_idx]
        self.seg_starts = np.concatenate([[0], np.nonzero(np.diff(sinv))[0] + 1])
        del inv, sinv
        self.core_rows = core                                         # M x Ns (int32)
        if verbose:
            print('  numerator monomials %d, core monomials %d (%.1fs)' % (len(alpha), len(self.ucore), time.time() - t0), flush=True)
        # quadratic weights (b1 <= b2) over linear weights 0..nw-1
        nw = net.nw
        self.qpairs = [(b1, b2) for b1 in range(nw) for b2 in range(b1, nw)] if order >= 2 else [(0, b) for b in range(nw)]
        self.qidx = {pr: k for k, pr in enumerate(self.qpairs)}
        self.G2 = np.array([[self.qidx[(min(a, b), max(a, b))] if (min(a, b), max(a, b)) in self.qidx else -1
                             for b in range(nw)] for a in range(nw)])
        # period lookups (prime independent)
        t0 = time.time()
        self.cl = [CoreLookup(net.skel_edges, self.q0 + k) for k in range(min(order, 2) + 1)]
        self.prep0 = self.cl[0].prepare(self.ucore)
        self.prep1 = []
        if order >= 1:
            for ti in range(len(net.trees)):
                self.prep1.append(self.cl[1].prepare(self.ucore + net.cot[ti][None, :]))
        self.ksums, self.KS, self.prep2 = [], None, []
        if order >= 2:
            nT = len(net.trees)
            kmap = {}
            KS = np.zeros((nT, nT), dtype=np.int64)
            for a in range(nT):
                for b in range(nT):
                    v = tuple(int(x) for x in net.cot[a] + net.cot[b])
                    if v not in kmap:
                        kmap[v] = len(self.ksums)
                        self.ksums.append(np.array(v, dtype=np.int64))
                    KS[a, b] = kmap[v]
            self.KS = KS
            cache2 = len(self.ksums) * len(self.ucore) <= 300000000
            for v in self.ksums:
                self.prep2.append(self.cl[2].prepare(self.ucore + v[None, :]) if cache2 else None)
        if verbose:
            print('  core %s: level %d %d/%d convergent; level %d: %d trees; level %d: %d tree sums (%.1fs)' % (
                self.cl[0].kind, self.q0, int(self.prep0['ok'].sum()), len(self.ucore), self.q0 + 1, len(self.prep1),
                self.q0 + 2, len(self.ksums), time.time() - t0), flush=True)
        self.pairs = pairs

    # --------------------------------------------------------------- push-forward
    def pushforward(self, coef, p):
        """PF[g] (len(qpairs) x Mc) residues: sum over numerator monomials of coef * Dirichlet(alpha, x^g)"""
        net, alpha = self.net, self.alpha
        top = int(self.core_rows.max()) + 4
        fact = [1]
        for i in range(1, top):
            fact.append(fact[-1] * i % p)
        fa = np.array(fact, dtype=np.int64)
        ifa = modlin.modinv_vec(fa.copy(), p)
        inv_int = np.array([0] + [pow(i, p - 2, p) for i in range(1, top)], dtype=np.int64)
        base = coef % p
        for ci, c in enumerate(net.chains):
            if len(c['edges']) < 2:
                continue
            for e in c['edges']:
                base = base * fa[alpha[:, e].astype(np.int64) - 1] % p
            base = base * ifa[self.core_rows[:, ci].astype(np.int64) - 1] % p
        out = np.zeros((len(self.qpairs), len(self.ucore)), dtype=np.int64)
        for g, (b1, b2) in enumerate(self.qpairs):
            w = base
            if b2 == 0:
                out[g] = self._segsum(w, p)
                continue
            if b1 == 0:
                ci, j = net.xw[b2 - 1]
                e = net.chains[ci]['edges'][j]
                ac = self.core_rows[:, ci].astype(np.int64)
                w = w * alpha[:, e].astype(np.int64) % p * inv_int[ac] % p
            else:
                c1, j1 = net.xw[b1 - 1]
                c2, j2 = net.xw[b2 - 1]
                e1, e2 = net.chains[c1]['edges'][j1], net.chains[c2]['edges'][j2]
                a1 = alpha[:, e1].astype(np.int64)
                a2 = alpha[:, e2].astype(np.int64)
                ac1 = self.core_rows[:, c1].astype(np.int64)
                if c1 != c2:
                    ac2 = self.core_rows[:, c2].astype(np.int64)
                    w = w * a1 % p * inv_int[ac1] % p * a2 % p * inv_int[ac2] % p
                elif e1 != e2:
                    w = w * a1 % p * a2 % p * inv_int[ac1] % p * inv_int[ac1 + 1] % p
                else:
                    w = w * a1 % p * (a1 + 1) % p * inv_int[ac1] % p * inv_int[ac1 + 1] % p
            out[g] = self._segsum(w, p)
        return out

    def _segsum(self, w, p):
        s = np.add.reduceat(w[self.order_idx], self.seg_starts)
        return s % p

    def _contract(self, PF, A, cl, prep, p, what):
        """dict label -> int64 vector over weights g: sum_core PF[g] * I(core)"""
        if prep is None:
            prep = cl.prepare(A)
        V, ok = cl.evaluate(A, prep, p)
        if (PF[:, ~ok] != 0).any():
            raise ArithmeticError('divergent period with nonzero coefficient (%s)' % what)
        labs = sorted(V)
        if not labs:
            return {}
        R = modlin.dotmod(PF, np.stack([V[l] for l in labs], axis=1), p)
        return {lab: R[:, k] for k, lab in enumerate(labs)}

    # --------------------------------------------------------------- one prime
    def run_prime(self, p):
        net = self.net
        t0 = time.time()
        coef = self.num.coefficients_mod(p)
        t1 = time.time()
        PF = self.pushforward(coef, p)
        del coef
        nT, nw, nq = len(net.trees), net.nw, len(self.qpairs)
        T0 = self._contract(PF, self.ucore, self.cl[0], self.prep0, p, 'level q0')
        labels = set(T0)
        res = {('bare', lab): int(T0[lab][0]) for lab in T0}
        if self.order >= 1:
            S1 = {}
            for ti in range(nT):
                r = self._contract(PF, self.ucore + net.cot[ti][None, :], self.cl[1], self.prep1[ti], p, 'level q0+1')
                for lab, v in r.items():
                    S1.setdefault(lab, np.zeros((nT, nq), dtype=np.int64))[ti] = v
            labels |= set(S1)
            lin = [self.qidx[(0, b)] for b in range(nw)]
            coeffs = {}
            for e in range(net.N):
                for i in range(net.H):
                    s, A1, A0 = net.insertion_coefficients(e, i)
                    A1 = A1.astype(np.int64)
                    A0 = A0.astype(np.int64)
                    coeffs[(e, i)] = (s, A1, A0)
                    for lab in labels:
                        val = 0
                        if lab in S1:
                            val += int((A1 * S1[lab][:, lin]).sum())
                        if lab in T0:
                            val += int((A0 * T0[lab][lin]).sum())
                        val = s * val % p
                        if val:
                            res[(('y', e, i), lab)] = val
        if self.order >= 2:
            S2 = {}
            for k, v in enumerate(self.ksums):
                r = self._contract(PF, self.ucore + v[None, :], self.cl[2], self.prep2[k], p, 'level q0+2')
                for lab, x in r.items():
                    S2.setdefault(lab, np.zeros((len(self.ksums), nq), dtype=np.int64))[k] = x
            labels |= set(S2)
            G2 = self.G2
            mats = {}
            for lab in labels:
                T0m = T0[lab][G2] if lab in T0 else np.zeros((nw, nw), dtype=np.int64)
                S1m = S1[lab][:, G2] if lab in S1 else np.zeros((nT, nw, nw), dtype=np.int64)
                if lab in S2:
                    S2m = S2[lab][self.KS[:, None, :, None], G2[None, :, None, :]]   # T b T' b'
                else:
                    S2m = np.zeros((nT, nw, nT, nw), dtype=np.int64)
                mats[lab] = (T0m, S1m, S2m)
            keys = sorted(coeffs)
            pairs = self.pairs if self.pairs is not None else [(a, b) for ia, a in enumerate(keys) for b in keys[ia:]]
            for a, b in pairs:
                s, A1, A0 = coeffs[a]
                s2, B1, B0 = coeffs[b]
                for lab in labels:
                    T0m, S1m, S2m = mats[lab]
                    val = int(A0 @ (T0m % p) @ B0)
                    val += int(np.einsum('tb,tbc,c->', A1, S1m, B0) + np.einsum('b,tbc,tc->', A0, S1m, B1))
                    val += int(np.einsum('tb,tbuc,uc->', A1, S2m, B1, optimize=True))
                    val = s * s2 * val % p
                    if val:
                        res[(('yy',) + a + b, lab)] = val
        if self.verbose:
            print('  prime %d: numerator %.1fs, assembly %.1fs' % (p, t1 - t0, time.time() - t1), flush=True)
        return res

    def run(self, primes=PRIMES, min_primes=3):
        crt = CRT()
        t0 = time.time()
        for p in primes:
            crt.add(p, self.run_prime(p))
            rec = crt.stable(min_primes)
            if rec is not None:
                if self.verbose:
                    print('  stable after %d primes (%.1fs)' % (crt.nprimes, time.time() - t0), flush=True)
                return unpack(rec)
        raise ArithmeticError('reconstruction did not stabilize')


def unpack(rec):
    out = {'bare': {}}
    for (k, lab), v in rec.items():
        if k == 'bare':
            out['bare'][lab] = v
        else:
            out.setdefault(k, {})[lab] = v
    return out


def kernel_moments(nv, edges, lam, order=1, verbose=True, check=True, primes=PRIMES, min_primes=3, pairs=None):
    """{'bare': K(lambda;0), ('y',e,i): int D y^i_e, ('yy',e,i,f,j): int D y^i_e y^j_f}, values {label: Fraction}"""
    net = Network(nv, edges, lam)
    if check:
        net.check_point()
    A = Assembly(net, order=order, verbose=verbose, pairs=pairs)
    return A.run(primes, min_primes)
