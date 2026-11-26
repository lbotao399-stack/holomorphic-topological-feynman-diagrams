"""Order-n nilpotency  sum_{a+b=n} Q_a Q_b = 0  on a test operator, modulo a prime.

  Q_0 = i hatQ_0,  Q_n = i/(8 pi^3)^n sum_Gamma kappa_Gamma hatQ_Gamma,  kappa_Gamma = wick_sign / |Aut|
  =>  sum_{a+b=n} sum_{Gamma_a, Gamma_b} kappa kappa' <Phi, hatQ_{Gamma_a} hatQ_{Gamma_b} O> = 0   (kappa_{Q_0} = 1).
Usage: python3 run_nilp.py n algebra O_types func_types F seed [--prime p] [--jets k:a.b] [--mode joint] [--nmom m]
  algebra: su2 | so5 | so6;  O_types e.g. b,b,b,b,b;  func_types e.g. c,c,c,c,c,g1,g2
"""
import sys
import time
import random
import argparse

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('n', type=int)
ap.add_argument('alg')
ap.add_argument('O')
ap.add_argument('ftypes')
ap.add_argument('F', type=int)
ap.add_argument('seed', type=int)
ap.add_argument('--prime', type=int, default=0)
ap.add_argument('--jets', default='', help='k:a.b;k2:c  (letter k of O gets d_{a+1} d_{b+1} ...)')
ap.add_argument('--mode', default='symbolic', help='symbolic | joint (kernels of two graph operators)')
ap.add_argument('--chunk', type=int, default=1024, help='kernel rows per evaluation')
ap.add_argument('--nmom', type=int, default=0, help='number of distinct functional momenta per set')
ap.add_argument('--only', default='', help='comma list of composition labels outer|inner to evaluate')
args = ap.parse_args()

import kern
if args.prime:
    kern.set_prime(args.prime)
import comp
from comp import Op, Functionals, evaluate, so_n, su2, check_lie
from kern import all_graphs, modq, ratrec

n = args.n
alg = args.alg
Otypes = args.O.split(',')
ftypes = args.ftypes.split(',')
F = args.F
seed = args.seed
jets = {}
if args.jets:
    for item in args.jets.split(';'):
        k, comps = item.split(':')
        jets[('O', int(k))] = tuple(int(c) for c in comps.split('.'))
P = kern.P
rng = random.Random(seed)
dim, f = su2() if alg == 'su2' else so_n(int(alg[2:]))
check_lie(f)
assert len(Otypes) == n + 1 and len(ftypes) == n + 3
O = tuple((('O', k), t) for k, t in enumerate(Otypes))
O_cols = {('O', k): np.array([rng.randint(-2, 2) for _ in range(dim)], dtype=np.int64) for k in range(len(Otypes))}
func = Functionals(ftypes, dim, F, rng, nmom=args.nmom)
print('order %d, %s, O = %s, jets %s, functionals %s, F = %d, seed %d, P = %d, mode %s, nmom %d'
      % (n, alg, ''.join(Otypes), jets, ''.join(ftypes), F, seed, P, args.mode, args.nmom), flush=True)
G = {m: all_graphs(m) for m in range(1, n + 1)}
T0 = time.time()
E = {}
for a in range(0, n + 1):
    b = n - a
    outer = [None] if a == 0 else G[a]
    inner = [None] if b == 0 else G[b]
    for ga in outer:
        for gb in inner:
            op1 = Op(gb, 'i')
            op2 = Op(ga, 'o')
            st = {}
            ka = 'Q0' if ga is None else ga.name
            kb = 'Q0' if gb is None else gb.name
            if args.only and '%s|%s' % (ka, kb) not in args.only.split(','):
                continue
            val = evaluate(op1, op2, O, O_cols, func, f, jets=jets, mode=args.mode, chunk=args.chunk, stats=st)
            E[(ka, kb)] = val
            print('  <hatQ[%s] hatQ[%s] O>: %s   (histories %d, keys %d, colour %.0fs, total %.0fs)'
                  % (ka, kb, [int(x) for x in val], st['hist'], st['keys'], st['tcol'], st['ttot']), flush=True)
if args.only:
    print('total time %.0fs' % (time.time() - T0))
    sys.exit(0)
kap = {'Q0': 1}
for m in G:
    for g in G[m]:
        kap[g.name] = modq(g.kappa)
tot = np.zeros(F, dtype=np.int64)
for (ka, kb), v in E.items():
    tot = (tot + kap[ka] * kap[kb] % P * v) % P
print('sum kappa kappa E (must vanish mod P):', [int(x) for x in tot], flush=True)
# kappa_n from the linear system, lower orders fixed
names = [g.name for g in G[n]]
rows, rhs = [], []
for fi in range(F):
    rows.append([int((E[('Q0', nm)][fi] + E[(nm, 'Q0')][fi]) % P) for nm in names])
    r = 0
    for (ka, kb), v in E.items():
        if ka != 'Q0' and kb != 'Q0':
            r = (r - kap[ka] * kap[kb] % P * int(v[fi])) % P
    rhs.append(int(r))


def rank_mod(M):
    M = [list(r) for r in M]
    rk = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % P), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], P - 2, P)
        M[rk] = [x * inv % P for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % P:
                fac = M[i][c]
                M[i] = [(x - fac * y) % P for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk, M


rk, _ = rank_mod(rows)
rka, red = rank_mod([r + [b] for r, b in zip(rows, rhs)])
print('kappa_%d system: %d equations, rank %d, augmented rank %d (unknowns %d)' % (n, F, rk, rka, len(names)))
if rk == len(names) == rka:
    sol = {}
    for i, nm in enumerate(names):
        sol[nm] = ratrec(red[i][-1])
    print('solution:', {k: str(v) for k, v in sol.items()})
    print('wick/Aut:', {g.name: str(g.kappa) for g in G[n]})
print('total time %.0fs' % (time.time() - T0))
