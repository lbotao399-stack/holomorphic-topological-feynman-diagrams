import sys, time
from fractions import Fraction as Fr
from salg import structure_constants, add, mul, scale
from ops import Theory, L

n, f = structure_constants('su2')
T = Theory(n, f)

def show(name, p):
    print(name, 'terms:', len(p), flush=True)

# Q0^2 = 0
tests0 = [L('c', 0), L('g1', 1), L('b', 2, (1, 0, 0)), mul(L('g1', 0), L('g2', 1)), mul(L('b', 0), L('c', 1, (0, 1, 0)))]
print('Q0^2 = 0:', all(T.Q0(T.Q0(p)) == {} for p in tests0))

# sample value: Q1hat(gamma^1_A gamma^2_B)
p = mul(L('g1', 0), L('g2', 1))
r = T.Q1(p)
print('Q1hat(g1^0 g2^1) =', {k: str(v) for k, v in r.items()})

tests = [
    mul(L('g1', 0), L('g2', 1)),
    mul(L('g1', 0), L('g2', 0)),
    mul(L('b', 0), L('b', 1)),
    mul(L('b', 0), L('g1', 1)),
    mul(L('b', 0), L('c', 1, (1, 0, 0))),
    mul(L('g1', 0, (1, 0, 0)), L('g2', 1)),
    mul(L('g1', 0, (1, 0, 0)), L('g2', 1, (0, 1, 0))),
    mul(L('b', 2, (0, 0, 1)), L('b', 1, (1, 0, 0))),
    mul(L('b', 0, (1, 1, 0)), L('g2', 1)),
    mul(mul(L('g1', 0), L('g2', 1)), L('b', 2)),
    mul(mul(L('b', 0), L('b', 1)), L('g1', 2, (0, 1, 0))),
]
ok = True
for p in tests:
    t0 = time.time()
    a = T.Q0(T.Q1(p)); b = T.Q1(T.Q0(p))
    s = add(a, b)
    print('  {Q0,Q1} on', list(p)[0], '->', 'ZERO' if s == {} else f'NONZERO ({len(s)} terms)', '| Q0Q1 terms', len(a), f'| {time.time()-t0:.1f}s', flush=True)
    ok &= (s == {})
print('{Q0hat, Q1hat} = 0 on all tests:', ok)
