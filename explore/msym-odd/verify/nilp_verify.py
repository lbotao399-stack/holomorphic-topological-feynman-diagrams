"""Closure test {Q0hat, sum_G k_G Qhat_G} O = 0 with the independent implementation (valg, vops, vkern).
Usage: python3 nilp_verify.py case alg nrandom seed [--ops "tok tok ..;tok .."] [--mutate double_z]
case: 5d (S2 graphs A, B with S2 eq-res-TA/TB), 5dN (theta_rootN (1,3,3), (2,2,3), arm tensor),
      7d (theta_rootN (1,4,4), (2,3,4), (3,3,3), arm), 3d (theta_rootN (1,2,2), arm).
token: <S>_<A>[^<jet digits>], S = digits of the theta indices (0-based) or 'e' for the empty set."""
import random
import sys
import time
from fractions import Fraction as Fr

from valg import Alg, Q0, so_structure_constants
from vops import KernelOp, TwoLoop, rank_null

args = list(sys.argv[1:])
ops_fixed, mutate = None, None
if '--ops' in args:
    i = args.index('--ops'); ops_fixed = args[i + 1].split(';'); del args[i:i + 2]
if '--mutate' in args:
    i = args.index('--mutate'); mutate = args[i + 1]; del args[i:i + 2]
case, alg = args[0], args[1]
nrand = int(args[2]) if len(args) > 2 else 20
seed = int(args[3]) if len(args) > 3 else 0

T0 = time.time()
if case == '5d':
    H, k = 2, 2
    names = ['A', 'B']; kops = {'A': KernelOp('S2A', 2, mutate=mutate), 'B': KernelOp('S2B', 2, mutate=mutate)}
    tensor = {'A': 'S2A', 'B': 'S2B'}
elif case == '5dN':
    H, k = 2, 2
    names = ['133', '223']
    kops = {'133': KernelOp('N', 2, (1, 3, 3), mutate), '223': KernelOp('N', 2, (2, 2, 3), mutate)}
    tensor = {n: 'arm' for n in names}
elif case == '7d':
    H, k = 3, 1
    names = ['144', '234', '333']
    kops = {n: KernelOp('N', 3, tuple(int(c) for c in n), mutate) for n in names}
    tensor = {n: 'arm' for n in names}
elif case == '3d':
    H, k = 1, 3
    names = ['122']
    kops = {'122': KernelOp('N', 1, (1, 2, 2), mutate)}
    tensor = {'122': 'arm'}
else:
    raise SystemExit('case?')
N = int(alg[2:])
n, f = so_structure_constants(N)
S = Alg(k, H)
q0 = Q0(S, n, f)
TL = TwoLoop(S, q0, kops, tensor)
print('%s, %s (dim %d), H=%d, k=%d, graphs %s, signs s = %s, mutate=%s' % (case, alg, n, H, k, names, [kops[x].s for x in names], mutate))
for x in names:
    print('   %s: K(lam;0)/pi weights for undifferentiated letters: %s' % (x, {m: str(v) for m, v in kops[x].weights(((0,) * H,) * 3).items()}))


def tok_letter(tk):
    if '^' in tk:
        a, j = tk.split('^'); jet = tuple(int(c) for c in j)
    else:
        a, jet = tk, (0,) * H
    Ss, A = a.split('_')
    Sset = () if Ss == 'e' else tuple(sorted(int(c) for c in Ss))
    return S.letter(Sset, int(A), jet)


def op_of(toks):
    p = {((), ()): Fr(1)}
    for tk in toks:
        p = S.mul(p, S.gen(tok_letter(tk)))
    return p


rng = random.Random(seed)
subsets = S.subsets


def rand_tok(nder):
    Sset = rng.choice(subsets)
    Ss = 'e' if not Sset else ''.join(map(str, Sset))
    jet = [0] * H
    for _ in range(nder):
        jet[rng.randrange(H)] += 1
    return '%s_%d^%s' % (Ss, rng.randrange(n), ''.join(map(str, jet)))


# Q0^2 = 0
t1 = time.time()
okq = True
cnt = 0
for _ in range(8):
    toks = [rand_tok(rng.choice([0, 1])) for _ in range(rng.randint(1, 3))]
    p = op_of(toks)
    if p:
        cnt += 1
        okq &= q0(q0(p)) == {}
print('Q0hat^2 = 0 on %d random operators: %s (%.1fs)' % (cnt, okq, time.time() - t1), flush=True)

if ops_fixed:
    tests = [(s_, op_of(s_.split())) for s_ in ops_fixed]
else:
    tests, seen = [], set()
    while len(tests) < nrand:
        nl = rng.choice([3, 3, 4])
        tot = rng.choice([1, 1, 2, 2, 3])
        ders = [0] * nl
        for _ in range(tot):
            ders[rng.randrange(nl)] += 1
        toks = [rand_tok(d_) for d_ in ders]
        key = tuple(sorted(toks))
        if key in seen:
            continue
        seen.add(key)
        p = op_of(toks)
        if p:
            tests.append((' '.join(toks), p))
rows = []
nq = nz = 0
for nm, p in tests:
    t1 = time.time()
    qs = [TL.Q(x, p) for x in names]
    an = [TL.anti(x, p) for x in names]
    nq += any(bool(q) for q in qs)
    mons = set().union(*[set(a) for a in an])
    eqs = [[a.get(mn, Fr(0)) for a in an] for mn in mons]
    rows += eqs
    nz += bool(eqs)
    if eqs or ops_fixed:
        rk, ns = rank_null(eqs, len(names))
        rkc, nsc = rank_null(rows, len(names))
        print('%-34s |Q_G O| %s, |{Q0,Q_G}O| %s, eqs %d, rank %d, null %s; cumulative %d eqs, rank %d, null %s (%.1fs)'
              % (nm, [len(q) for q in qs], [len(a) for a in an], len(eqs), rk, [[str(c) for c in v] for v in ns],
                 len(rows), rkc, [[str(c) for c in v] for v in nsc], time.time() - t1), flush=True)
if rows:
    rk, ns = rank_null(rows, len(names))
    print('TOTAL: %d operators, %d with some Q_G O != 0, %d with nonzero anticommutators; %d equations, rank %d, null space %s'
          % (len(tests), nq, nz, len(rows), rk, [[str(c) for c in v] for v in ns]))
else:
    print('TOTAL: %d operators, %d with some Q_G O != 0: all anticommutators vanish' % (len(tests), nq))
print('total %.1fs' % (time.time() - T0))
