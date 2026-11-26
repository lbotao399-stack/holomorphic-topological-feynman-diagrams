"""Q_0^2 = 0 and {hat Q_0, hat Q_1} = 0 for holomorphic CS on C^{4|1} (adapted from S2 code/test_q0q1.py)."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
from ops_c41 import structure_constants, add, mul, Theory, L

name = sys.argv[1] if len(sys.argv) > 1 else 'su2'
n, f = structure_constants(name)
T = Theory(n, f)
print('algebra', name, 'dim', n)
e = lambda a: tuple(1 if i == a else 0 for i in range(4))

tests0 = [L('c', 0), L('g', 1), L('g', 2, e(0)), mul(L('g', 0), L('g', 1)), mul(L('g', 0), L('c', 1, e(1)))]
print('Q0^2 = 0:', all(T.Q0(T.Q0(p)) == {} for p in tests0), flush=True)

p = mul(L('g', 0), L('g', 1))
r = T.Q1(p)
print('hatQ1(g_0 g_1): %d terms' % len(r))
for k, v in list(r.items())[:6]:
    print('   ', k, v)

tests = [
    ('g0 g1', mul(L('g', 0), L('g', 1))),
    ('g0 g0', mul(L('g', 0), L('g', 0))),
    ('g0 c1', mul(L('g', 0), L('c', 1))),
    ('c0 c1', mul(L('c', 0), L('c', 1))),
    ('d1g0 g1', mul(L('g', 0, e(0)), L('g', 1))),
    ('d1g0 d2g1', mul(L('g', 0, e(0)), L('g', 1, e(1)))),
    ('d1d2g0 g1', mul(L('g', 0, (1, 1, 0, 0)), L('g', 1))),
    ('d4g2 d1c1', mul(L('g', 2, e(3)), L('c', 1, e(0)))),
    ('d1d1c0 g1', mul(L('c', 0, (2, 0, 0, 0)), L('g', 1))),
    ('g0 g1 g2', mul(mul(L('g', 0), L('g', 1)), L('g', 2))),
    ('g0 g1 c2', mul(mul(L('g', 0), L('g', 1)), L('c', 2))),
    ('g0 c1 d3c2', mul(mul(L('g', 0), L('c', 1)), L('c', 2, e(2)))),
]
ok = True
for nm, p in tests:
    t0 = time.time()
    a = T.Q0(T.Q1(p)); b = T.Q1(T.Q0(p))
    s = add(a, b)
    print('  {Q0,Q1} on %-12s -> %s | terms Q1(p) %d, Q0Q1(p) %d, Q1Q0(p) %d | %.1fs'
          % (nm, 'ZERO' if s == {} else 'NONZERO (%d terms)' % len(s), len(T.Q1(p)), len(a), len(b), time.time() - t0), flush=True)
    ok &= (s == {})
print('{hatQ0, hatQ1} = 0 on all tests:', ok)
