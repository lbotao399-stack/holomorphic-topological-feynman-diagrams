"""[fork B copy of run_nilp_f.py with --custom ALG and --H overrides: ALG in {berezin:G:k, glmn:m:n}]
Order-n nilpotency  sum_{a+b=n} Q_a Q_b O = 0  with plane-wave functionals (fmode), modulo P.

Q_n = c^n sum_Gamma kappa_Gamma hatQ_Gamma (c drops out).  kappa_{<n} read from --kappa (json), kappa_n solved from the
F functional sets (rank printed) and compared with --expect if given.  Writes the extended json to --out.
usage: python3 run_nilp_f.py THEORY GAUGE MATTER n O_blocks func_blocks F seed [--jets k:a.b;..] [--kappa f.json]
       [--out f.json] [--prime p] [--unicyclic]
"""
import argparse
import json
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('theory')
ap.add_argument('gauge')
ap.add_argument('matter')
ap.add_argument('n', type=int)
ap.add_argument('O')
ap.add_argument('ftypes')
ap.add_argument('F', type=int)
ap.add_argument('seed', type=int)
ap.add_argument('--jets', default='')
ap.add_argument('--kappa', default='')
ap.add_argument('--out', default='')
ap.add_argument('--prime', type=int, default=0)
ap.add_argument('--chunk', type=int, default=1024)
ap.add_argument('--mode', default='symbolic')
ap.add_argument('--only', default='')
ap.add_argument('--blocks', default='', help='block partition: name of a scheme in blocks_for()')
ap.add_argument('--unicyclic', action='store_true')
ap.add_argument('--custom', default='')
ap.add_argument('--H', type=int, default=0)
args = ap.parse_args()

import kern_h
if args.prime:
    kern_h.set_prime(args.prime)
import polyops as po
po.set_prime(kern_h.P)
import fmode
fmode.P = kern_h.P
import theories
from wick import Algebra
from explicit import graph_ops
from fmode import FOp, Context, Functionals, evaluate

P = kern_h.P


def blocks_for(sl):
    """partition of the basis into constant-parity blocks: finest available from sl.blocks, else parity."""
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


matter = None if args.matter == '-' else args.matter
if args.custom:
    import superalg as _sa
    kind = args.custom.split(':')
    if kind[0] == 'berezin':
        sl = _sa.berezin_tensor(theories.gauge(kind[1]), int(kind[2]))
    elif kind[0] == 'glmn':
        import b_superlie
        sl = b_superlie.glmn(int(kind[1]), int(kind[2]))
    H = args.H
else:
    H, sl = theories.build(args.theory, args.gauge, matter)
    if args.H:
        H = args.H
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
O_cols = {}
for k, t in enumerate(Otypes):
    v = np.zeros(alg.N)
    for a in blocks[t]:
        v[a] = rng.randint(-2, 2)
    O_cols[('O', k)] = v
func = Functionals(ctx, ftypes, args.F, rng)
print('theory %s gauge %s matter %s H=%d N=%d: order %d, O = %s, jets %s, functionals %s, F = %d, seed %d, P = %d'
      % (args.theory, args.gauge, matter, H, alg.N, n, ','.join(Otypes), jets, ','.join(ftypes), args.F, args.seed, P),
      flush=True)
print('blocks:', {b: len(ids) for b, ids in blocks.items()}, flush=True)
G = {m: graph_ops(alg, m, unicyclic=args.unicyclic, drop_zero=False) for m in range(1, n + 1)}
import b_theta333
for m in G:
    for op in G[m]:
        gk, sw = b_theta333.make_gk(H, dict(nb=op.nb, edges=op.edges))
        if sw:
            op.gk = gk
            print('  [%s: block-triangular numeric density]' % op.name, flush=True)
for m in G:
    print('  order %d graphs: %s' % (m, [(g.name, g.aut, g.edges) for g in G[m]]), flush=True)
kap = {}
if args.kappa:
    kap = {k: Fr(v) for k, v in json.load(open(args.kappa)).items()}
T0 = time.time()
E = {}
for a in range(0, n + 1):
    b = n - a
    outer = [None] if a == 0 else G[a]
    inner = [None] if b == 0 else G[b]
    for ga in outer:
        for gb in inner:
            op1 = FOp(alg, gb, 'i')
            op2 = FOp(alg, ga, 'o')
            ka = 'Q0' if ga is None else ga.name
            kb = 'Q0' if gb is None else gb.name
            if args.only and '%s|%s' % (ka, kb) not in args.only.split(','):
                continue
            st = {}
            val = evaluate(ctx, op1, op2, O, O_cols, func, jets=jets, mode=args.mode, chunk=args.chunk, stats=st)
            E[(ka, kb)] = val
            print('  <hatQ[%s] hatQ[%s] O>: %s   (histories %d, keys %d, colour %.0fs, total %.0fs)'
                  % (ka, kb, [int(x) for x in val], st['hist'], st['keys'], st['tcol'], st['ttot']), flush=True)
if args.only:
    sys.exit(0)
kp = {'Q0': 1}
for m in G:
    if m == n:
        continue
    for g in G[m]:
        if m == 1 and g.name not in kap:
            kap[g.name] = Fr(1, g.aut)
        kp[g.name] = po.modq(kap[g.name])
names = [g.name for g in G[n]]
rows, rhs = [], []
for fi in range(args.F):
    rows.append([int((E[('Q0', nm)][fi] + E[(nm, 'Q0')][fi]) % P) for nm in names])
    r = 0
    for (ka, kb), v in E.items():
        if ka != 'Q0' and kb != 'Q0':
            r = (r - kp[ka] * kp[kb] % P * int(v[fi])) % P
    rhs.append(int(r))
from nilp_explicit import rank_solve
rk, aug, sol = rank_solve(rows, rhs, P)
print('kappa_%d system: %d equations, rank %d, augmented rank %d (unknowns %d)' % (n, args.F, rk, aug, len(names)),
      flush=True)
if n == 1:
    print('order 1 residuals {Q0,hatQ_G}:', {nm: [r[i] for r in rows] for i, nm in enumerate(names)}, flush=True)
if sol is not None:
    for nm, s, g in zip(names, sol, G[n]):
        kap[nm] = po.ratrec(s)
        print('  kappa[%s] = %s   |Aut| = %d   kappa*|Aut| = %s' % (nm, kap[nm], g.aut, kap[nm] * g.aut), flush=True)
    # residual with the solution
    res = []
    for fi in range(args.F):
        rr = (sum(int(sol[i]) * rows[fi][i] for i in range(len(names))) - rhs[fi]) % P
        res.append(rr)
    print('residual with solved kappa_%d:' % n, res, flush=True)
elif rk == aug:
    print('consistent but underdetermined (rank %d < %d)' % (rk, len(names)), flush=True)
else:
    print('INCONSISTENT', flush=True)
if args.out:
    json.dump({k: str(v) for k, v in kap.items()}, open(args.out, 'w'), indent=1)
print('total time %.0fs' % (time.time() - T0), flush=True)
