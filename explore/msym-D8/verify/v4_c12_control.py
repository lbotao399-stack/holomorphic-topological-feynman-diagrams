"""Control for C12: the same system with deliberately wrong first-order theta weights
(inputs z_1 and z_2 exchanged in the first-order weights of theta(2,3,4)); expected: inconsistent."""
import sys, time, pickle
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify')
import v4_ops as V
from v4_ops import Theory, letter as Lt, op, e, padd, Fr, log
from v4_c12 import solve, THETAS

V.LOG = open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/v4_c12_control.log', 'w')
tw = pickle.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/theta_weights.pkl', 'rb'))
W1 = tw[(2, 3, 4)]['W1']
for b in range(4):
    W1[(0, b)], W1[(1, b)] = W1[(1, b)], W1[(0, b)]
n, f = V.su2()
T = Theory(n, f, tw)
g = lambda A, k=(0, 0, 0, 0): Lt('g', A, k)
c = lambda A, k=(0, 0, 0, 0): Lt('c', A, k)
p = op(c(0, e(1)), g(1), g(2))
t0 = time.time()
lhs = T.Q1(T.Q1(p))
q0p = T.Q0(p)
anti = {G: padd(T.Q0(T.Qtheta(G, p)), T.Qtheta(G, q0p)) for G in THETAS}
mons = set(lhs)
for a in anti.values():
    mons |= set(a)
eqs = [[Fr(anti[G].get(m, 0)) for G in THETAS] + [-Fr(lhs.get(m, 0))] for m in mons]
r_, s_, c_ = solve(eqs, 3)
log('CONTROL (z_1 <-> z_2 in the first-order weights of theta(2,3,4)), operator d2c0 g1 g2: %d equations, rank %s, consistent %s, solution %s (%.1fs)'
    % (len(eqs), r_, c_, s_, time.time() - t0))
