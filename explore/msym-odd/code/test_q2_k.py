"""{Q_0, sum_Gamma k_Gamma Qhat_Gamma} = 0 with the generic code (salgk, opsk) for
  5d: H = 2, k = 2, graphs (A), (B) of S2 (cross-check of ops5 / test_q2_5d);
  7d: H = 3, k = 1, theta(1,4,4), theta(2,3,4), theta(3,3,3) with O = N;
  3d: H = 1, k = 3, theta(1,2,2) with O = N.
Usage: python3 test_q2_k.py 5d|7d|3d so3|so5 [nrandom] [seed] [--ops op1,op2,...]
Operators: tokens  <type>_<colour>[^j1j2j3], type = mask of the theta component (0 = c ...)."""
import random
import sys
import time

import sympy as sp

import htk_path  # noqa: F401
import htk
from salgk import SAlg, structure_constants
from opsk import TheoryK, ThetaOp, TwoLoopK
from kernels_check import theta_rootN

case = sys.argv[1]
gname = sys.argv[2]
fixed = None
argv = list(sys.argv)
if '--ops' in argv:
    i = argv.index('--ops')
    fixed = argv[i + 1].split(';')
    argv = argv[:i] + argv[i + 2:]
nrand = int(argv[3]) if len(argv) > 3 else 30
seed = int(argv[4]) if len(argv) > 4 else 1

t0 = time.time()
if case == '5d':
    H, k = 2, 2
    graphs = {'A': (htk.graph_A(), [[3], [0, 1, 2], [6, 5, 4]]), 'B': (htk.graph_B(), [[0, 1], [2, 3], [6, 5, 4]])}
    aut = {'A': 2, 'B': 2}
elif case == '5dN':
    H, k = 2, 2
    graphs, aut = {}, {}
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        graphs[str(pqr)] = theta_rootN(*pqr)
        aut[str(pqr)] = 2
elif case == '7d':
    H, k = 3, 1
    graphs, aut = {}, {}
    for pqr, a in [((1, 4, 4), 2), ((2, 3, 4), 1), ((3, 3, 3), 6)]:
        graphs[str(pqr)] = theta_rootN(*pqr)
        aut[str(pqr)] = a
elif case == '3d':
    H, k = 1, 3
    graphs = {'(1, 2, 2)': theta_rootN(1, 2, 2)}
    aut = {'(1, 2, 2)': 2}
else:
    raise SystemExit("case must be 5d, 7d or 3d")
alg = SAlg(k)
n, f = structure_constants(gname)
th = TheoryK(alg, n, f)
ops = {nm: ThetaOp(G, paths, H, nm) for nm, (G, paths) in graphs.items()}
T2 = TwoLoopK(th, ops)
kinds = list(graphs)
print('%s (H=%d, k=%d), %s: graphs %s, |Aut| = %s, s = %s  (%.1fs)'
      % (case, H, k, gname, kinds, [aut[g] for g in kinds], [ops[g].s for g in kinds], time.time() - t0))
for g in kinds:
    w0 = ops[g].weights(((0, 0, 0),) * 3)
    print('   %s: D = %s;  K(lam;0)/pi = %s * D' % (g, sp.factor(ops[g].D), sp.cancel(ops[g].TK.K0()[1] / ops[g].TK.Dtheta())))


def letter(tok):
    if '^' in tok:
        a, j = tok.split('^')
        jets = tuple(int(c) for c in j) + (0,) * (3 - len(j))
    else:
        a, jets = tok, (0, 0, 0)
    mask, A = a.split('_')
    return (1 + int(mask), int(A)) + jets


def op_of(tokens):
    p = alg.const(1)
    for tk in tokens:
        p = alg.mul(p, alg.gen(letter(tk)))
    return p


# Q_0^2 = 0
rng = random.Random(seed)
masks = list(range(1 << k))
tests0 = []
for _ in range(8):
    toks = ['%d_%d^%s' % (rng.choice(masks), rng.randrange(n), ''.join(str(rng.randint(0, 1)) for _ in range(H)))
            for _ in range(rng.randint(1, 3))]
    tests0.append((toks, op_of(toks)))
t1 = time.time()
ok0 = all(th.Q0(th.Q0(p)) == {} for _, p in tests0 if p)
print('Q0^2 = 0 on %d random operators: %s (%.1fs)' % (len([1 for _, p in tests0 if p]), ok0, time.time() - t1), flush=True)

if fixed:
    tests = [(s, op_of(s.split())) for s in fixed]
else:
    tests, seen = [], set()
    jets1 = [tuple(1 if i == a else 0 for i in range(H)) for a in range(H)]
    while len(tests) < nrand:
        nl = rng.choice([3, 3, 3, 4])
        toks = []
        nder = rng.choice([0, 1, 1, 2])
        for i in range(nl):
            j = ''
            if i < nder:
                j = '^' + ''.join(map(str, rng.choice(jets1)))
            toks.append('%d_%d%s' % (rng.choice(masks), rng.randrange(n), j))
        key = tuple(sorted(toks))
        if key in seen:
            continue
        seen.add(key)
        p = op_of(toks)
        if p:
            tests.append((' '.join(toks), p))
rows = []
nz = 0
nq = 0
for nm, p in tests:
    t1 = time.time()
    vals = []
    qs = [len(T2.Q2(g, p)) for g in kinds]
    nq += any(qs)
    for g in kinds:
        a = alg.add(th.Q0(T2.Q2(g, p)), T2.Q2(g, th.Q0(p)))
        vals.append(a)
    mons = set().union(*[set(v) for v in vals])
    eqs = [[sp.Rational(v.get(m, 0)) for v in vals] for m in mons]
    rows += eqs
    nz += bool(eqs)
    if eqs or fixed:
        rk = sp.Matrix(eqs).rank() if eqs else 0
        cum = [list(v) for v in sp.Matrix(rows).nullspace()] if rows else 'all'
        print('%-30s |Q_G O|: %s, |{Q0,Q_G}O|: %s, equations %d, rank %d; cumulative: %d equations, null space %s  (%.1fs)'
              % (nm, qs, [len(v) for v in vals], len(eqs), rk, len(rows), cum, time.time() - t1), flush=True)
if rows:
    Mx = sp.Matrix(rows)
    ns = Mx.nullspace()
    print('%d operators, %d with Q_G O != 0, %d with nonzero anticommutators; total equations %d, rank %d, null space %s'
          % (len(tests), nq, nz, len(rows), Mx.rank(), [list(v) for v in ns]))
    pred = [sp.Rational(1, aut[g]) for g in kinds]
    print('   1/|Aut| per graph:', pred)
else:
    print('%d operators, %d with Q_G O != 0: all anticommutators vanish' % (len(tests), nq))
print('total %.1fs' % (time.time() - t0))
