"""Order-n nilpotency sum_{a+b=n} Q_a Q_b O = 0 with pooled-colour plane-wave functionals (engine: ../qn-holomorphic).

Change w.r.t. fmode.evaluate: functional i of set f has colour pool[type_i][r_{f,i}] (R pool vectors per block) and
momentum base_f[m_{f,i}] (nmom classes); (r, m) distinct within a block.  Per history the colour network is contracted
once with the pool (tensor R^N, independent of F) and gathered at the N! slot permutations, instead of the dense
(F, n_0..n_{N-1}) tensor.  Kernels, Koszul signs, norms: fmode unchanged.
usage: poolnilp.py THEORY GAUGE MATTER n O_blocks func_blocks F seed [--R r] [--nmom m] [--only ka|kb,..] [--save f.json]
       [--kappa f.json] [--jets k:a.b;..] [--unicyclic] [--prime p]
"""
import argparse
import itertools
import json
import os
import random
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ENG = os.path.join(HERE, '..', 'qn-holomorphic')
sys.path.insert(0, ENG)

ap = argparse.ArgumentParser()
ap.add_argument('theory')
ap.add_argument('gauge')
ap.add_argument('matter')
ap.add_argument('n', type=int)
ap.add_argument('O')
ap.add_argument('ftypes')
ap.add_argument('F', type=int)
ap.add_argument('seed', type=int)
ap.add_argument('--R', type=int, default=3)
ap.add_argument('--nmom', type=int, default=4)
ap.add_argument('--jets', default='')
ap.add_argument('--only', default='')
ap.add_argument('--save', default='')
ap.add_argument('--prime', type=int, default=0)
ap.add_argument('--chunk', type=int, default=2048)
ap.add_argument('--unicyclic', action='store_true')
ap.add_argument('--list', action='store_true')
ap.add_argument('--fsub', default='')
ap.add_argument('--flav', type=int, default=0, help='6d10 fund*K: u colours only on flavour 1 (local u positions 0,1,2K,2K+1)')
args = ap.parse_args()

import kern_h
if args.prime:
    kern_h.set_prime(args.prime)
import polyops as po
po.set_prime(kern_h.P)
import b_fastdet  # noqa: F401  (patch)
import fmode
fmode.P = kern_h.P
import theories
from wick import Algebra
from explicit import graph_ops
from fmode import FOp, Context, apply_types, colour_values, key_kernel

P = kern_h.P


def blocks_for(sl):
    N = sl.N
    names = sorted(sl.blocks, key=lambda b: (len(sl.blocks[b]), b))
    cover = {}
    for b in names:
        ids = sl.blocks[b]
        if len({int(sl.pi[i]) for i in ids}) != 1:
            continue
        if any(i in cover for i in ids):
            continue
        for i in ids:
            cover[i] = b
    out = {}
    for i in range(N):
        b = cover.get(i, 'par%d' % int(sl.pi[i]))
        out.setdefault(b, []).append(i)
    return out


class PoolFunctionals:
    def __init__(self, ctx, types, F, rng, R, nmom, amp_col=3, amp_mom=5):
        self.types = list(types)
        self.N = N = len(types)
        self.F = F
        H = ctx.alg.H
        dim = ctx.alg.N
        self.utypes = sorted(set(types))
        self.R = R
        self.pool = {}
        for t in self.utypes:
            M = np.zeros((R, dim), dtype=np.int64)
            for r in range(R):
                for a in allowed(t):
                    M[r, a] = rng.randint(-amp_col, amp_col)
            self.pool[t] = M
        self.ridx = np.zeros((F, N), dtype=np.int64)
        self.u = np.zeros((F, N, H), dtype=np.int64)
        self.mcls = np.zeros((F, N), dtype=np.int64)
        for f in range(F):
            base = [[rng.randint(-amp_mom, amp_mom) for _ in range(H)] for _ in range(nmom)]
            for t in self.utypes:
                ii = [i for i in range(N) if types[i] == t]
                grid = [(r, m) for r in range(R) for m in range(nmom)]
                assert len(grid) >= len(ii), 'R*nmom too small for block %s' % t
                pick = rng.sample(grid, len(ii))
                for i, (r, m) in zip(ii, pick):
                    self.ridx[f, i] = r
                    self.u[f, i] = base[m]
            seen = {}
            for i in range(N):
                self.mcls[f, i] = seen.setdefault(tuple(self.u[f, i]), len(seen))
        perms = np.array(list(itertools.permutations(range(N))), dtype=np.int64)
        self.perms = perms
        self.pat_of = {}
        for r, s in enumerate(perms):
            self.pat_of.setdefault(tuple(types[i] for i in s), []).append(r)
        self.pat_of = {k: np.array(v, dtype=np.int64) for k, v in self.pat_of.items()}
        odd = np.array([ctx.bpar[t] for t in types])
        sg = []
        for s in perms:
            o = [i for i in s if odd[i]]
            inv = sum(1 for a in range(len(o)) for b in range(a + 1, len(o)) if o[a] > o[b])
            sg.append(-1 if inv & 1 else 1)
        self.sg = np.array(sg, dtype=np.int64)
        self.nrows = len(perms)
        # pseudo functional set for the pooled contraction: F = 1, rows = pool vectors by block
        rows, self.prow = [], {}
        for t in self.utypes:
            self.prow[t] = np.arange(len(rows), len(rows) + R, dtype=np.int64)
            rows += list(self.pool[t])
        self.pseudo = type('PF', (), {})()
        self.pseudo.F = 1
        self.pseudo.Af = np.array(rows, dtype=np.float64)[None, :, :]


def evaluate_pool(ctx, op1, op2, O, O_cols, func, jets=None, chunk=2048, stats=None):
    alg = ctx.alg
    H = alg.H
    jets = {t: tuple(v) for t, v in (jets or {}).items() if v}
    t0 = time.time()
    O = tuple(O)
    Otags = [t for t, _ in O]
    Gsum, infos = {}, {}
    nhist = 0
    F, N = func.F, func.N
    norm = pow(alg.Gden, op1.nG + op2.nG, P)
    norm = pow(norm, P - 2, P)
    if op1.gop is None:
        norm = norm * pow(2, P - 2, P) % P
    if op2.gop is None:
        norm = norm * pow(2, P - 2, P) % P
    for bind1 in itertools.permutations(Otags, op1.nin):
        terms1 = apply_types(op1, ctx, O, bind1)
        for s1, w1, combo1, q1 in terms1:
            tags1 = [t for t, _ in w1]
            for bind2 in itertools.permutations(tags1, op2.nin):
                terms2 = apply_types(op2, ctx, w1, bind2)
                for s2, w2, combo2, q2 in terms2:
                    slots = [t for t, _ in w2]
                    pat = tuple(b for _, b in w2)
                    if pat not in func.pat_of:
                        continue
                    nhist += 1
                    wtypes = dict(O)
                    wtypes.update(dict(zip(op1.outtags(), combo1)))
                    wtypes.update(dict(zip(op2.outtags(), combo2)))
                    jetloc = []
                    for t in Otags:
                        if t in jets:
                            if t in bind1:
                                jetloc.append((t, ('i', bind1.index(t))))
                            elif t in bind2:
                                jetloc.append((t, ('o', bind2.index(t))))
                            else:
                                jetloc.append((t, ('s', slots.index(t))))
                    jetloc = tuple(jetloc)
                    smap = tuple(b[1] if b[0] == op1.prefix else None for b in bind2)
                    if op1.gop is None:
                        key = ('Q0in', op2.kind, smap if jetloc else None, jetloc)
                    elif op2.gop is None:
                        key = ('Q0out', op1.kind, smap, jetloc)
                    else:
                        key = ('GG', op1.kind, op2.kind, smap, jetloc)
                    infos[key] = dict(ops=(op1, op2), smap=smap, jets=jets, jetloc=jetloc)
                    slot_funcs = [func.prow[pat[s]] for s in range(N)]
                    Gv = colour_values(ctx, (op1, op2), (bind1, bind2), wtypes, O_cols, slots, slot_funcs,
                                       func.pseudo, qs=(q1, q2))[0]
                    if not Gv.any():
                        continue
                    rows = func.pat_of[pat]
                    acc = Gsum.setdefault(key, {}).setdefault(pat, np.zeros((F, len(rows)), dtype=np.int64))
                    S = func.perms[rows]
                    idx = tuple(func.ridx[:, S[:, s]] for s in range(N))
                    vals = Gv[idx] % P
                    acc += vals * (s1 * s2 % P) % P
                    acc %= P
    tcol = time.time() - t0
    total = np.zeros(F, dtype=np.int64)
    nkeys, nker = 0, 0
    for key, pacc in Gsum.items():
        FI, RI, VV = [], [], []
        for pat, acc in pacc.items():
            rows = func.pat_of[pat]
            acc = acc * (func.sg[rows][None, :] % P) % P
            fi, li = np.nonzero(acc)
            FI.append(fi)
            RI.append(rows[li])
            VV.append(acc[fi, li])
        fi, ri, vv = np.concatenate(FI), np.concatenate(RI), np.concatenate(VV)
        if len(fi) == 0:
            continue
        nkeys += 1
        S_all = func.perms[ri]
        code = fi.copy()
        for s in range(N):
            code = code * N + func.mcls[fi, S_all[:, s]]
        ucode, first, inv = np.unique(code, return_index=True, return_inverse=True)
        gs = np.zeros(len(ucode), dtype=np.int64)
        np.add.at(gs, inv, vv)
        gs %= P
        nz = gs != 0
        uf, ur, gs = fi[first][nz], ri[first][nz], gs[nz]
        nker += len(uf)
        for c0 in range(0, len(uf), chunk):
            fc, rc = uf[c0:c0 + chunk], ur[c0:c0 + chunk]
            L = len(fc)
            S = func.perms[rc]
            mom = [np.array([func.u[fc, S[:, s], a] % P for a in range(H)], dtype=np.int64) for s in range(N)]
            K = key_kernel(infos[key], mom, L, H)
            contrib = K * gs[c0:c0 + chunk] % P
            np.add.at(total, fc, contrib)
            total %= P
    total = total * norm % P
    if stats is not None:
        stats.update(hist=nhist, keys=nkeys, nker=nker, tcol=tcol, ttot=time.time() - t0)
    return total % P


matter = None if args.matter == '-' else args.matter
H, sl = theories.build(args.theory, args.gauge, matter)
alg = Algebra(sl, H)
blocks = blocks_for(sl)
ctx = Context(alg, blocks)
n = args.n
rng = random.Random(args.seed)
Otypes = args.O.split(',')
ftypes = args.ftypes.split(',')
for t in Otypes + ftypes:
    assert t in blocks, 'unknown block %s; blocks: %s' % (t, sorted(blocks))
jets = {}
if args.jets:
    for item in args.jets.split(';'):
        k, comps = item.split(':')
        jets[('O', int(k))] = tuple(int(c) for c in comps.split('.'))
O = tuple((('O', k), t) for k, t in enumerate(Otypes))
def allowed(t):
    if args.flav and t == 'u':
        K = args.flav
        return [blocks[t][i] for i in (0, 1, 2 * K, 2 * K + 1)]
    return blocks[t]


O_cols = {}
for k, t in enumerate(Otypes):
    v = np.zeros(alg.N)
    for a in allowed(t):
        v[a] = rng.randint(-2, 2)
    O_cols[('O', k)] = v
func = PoolFunctionals(ctx, ftypes, args.F, rng, args.R, args.nmom)
if args.fsub:
    a, b = map(int, args.fsub.split(':'))
    func.F = b - a
    func.ridx, func.u, func.mcls = func.ridx[a:b], func.u[a:b], func.mcls[a:b]
    print('functional subset %d:%d of %d' % (a, b, args.F), flush=True)
print('theory %s gauge %s matter %s H=%d N=%d: order %d, O = %s, jets %s, functionals %s, F = %d, seed %d, R = %d, '
      'nmom = %d, P = %d' % (args.theory, args.gauge, matter, H, alg.N, n, ','.join(Otypes), jets, ','.join(ftypes),
                            args.F, args.seed, args.R, args.nmom, P), flush=True)
print('blocks:', {b: len(ids) for b, ids in blocks.items()}, flush=True)
G = {m: graph_ops(alg, m, unicyclic=args.unicyclic, drop_zero=False) for m in range(1, n + 1)}
for m in G:
    print('  order %d graphs: %s' % (m, [(g.name, g.aut) for g in G[m]]), flush=True)
if args.list:
    sys.exit(0)
import b_theta333
for m in G:
    for op in G[m]:
        if not any(len(mv.chains) >= 3 for mv in op.gk.red.moves):
            continue
        gk_sym = op.gk
        if gk_sym.red.tight_plan:
            gk_bt, sw = b_theta333.make_gk(H, dict(nb=op.nb, edges=op.edges))
        else:
            gk_bt, sw = gk_sym, True
        gk_bt.density_factors()
        if m < n and gk_sym.red.tight_plan:
            gk_bt.kernel_symbolic = gk_sym.kernel_symbolic
            gk_bt.red.tight_plan = True
        op.gk = gk_bt
        print('  [%s: block-triangular numeric density%s]' % (op.name, ', symbolic kernel from the tight plan'
                                                              if m < n else ''), flush=True)
T0 = time.time()
E = {}
only = set(args.only.split(',')) if args.only else None
for a in range(0, n + 1):
    b = n - a
    outer = [None] if a == 0 else G[a]
    inner = [None] if b == 0 else G[b]
    for ga in outer:
        for gb in inner:
            ka = 'Q0' if ga is None else ga.name
            kb = 'Q0' if gb is None else gb.name
            if only is not None and '%s|%s' % (ka, kb) not in only:
                continue
            op1 = FOp(alg, gb, 'i')
            op2 = FOp(alg, ga, 'o')
            st = {}
            val = evaluate_pool(ctx, op1, op2, O, O_cols, func, jets=jets, chunk=args.chunk, stats=st)
            E['%s|%s' % (ka, kb)] = [int(x) for x in val]
            print('  <hatQ[%s] hatQ[%s] O>: %s   (histories %d, keys %d, kernel pts %d, colour %.0fs, total %.0fs)'
                  % (ka, kb, [int(x) for x in val], st['hist'], st['keys'], st['nker'], st['tcol'], st['ttot']),
                  flush=True)
            if args.save:
                json.dump(dict(P=P, F=func.F, fsub=args.fsub, E=E), open(args.save, 'w'), indent=0)
print('total time %.0fs' % (time.time() - T0), flush=True)
