"""(d), (e) Cochains (formal gl_N, kappa = tr; sl_N by restriction):
    o'        = sum_G |Aut G|^{-1} alt(K_G/pi W_G)                         (o_2 = pi nu o')
    theta'(T) = sum_G |Aut G|^{-1} sum_v alt(K_G^(v)/pi W_G^(v,T))           (Theta(J_1) = pi nu x theta'(J), J_1 = sum x_T J_T)
W^(v,T): vertex v carries T in {d, s1, s0} instead of f; |Aut| = 4, 4.
Checks: o' = (1/240)(-3 F_133 + 2 F_223) (deform spaces.obstruction_cochains);
        (e) d theta'(T) = {J_T, o'} exactly (formal N), T = d, s1, s0 (gl_N) and T = d (sl_N);
        controls: d theta'(T) != -{J_T, o'};  {J_T, o'} != 0.
Saves data/theta.pkl."""
import time
from common import *

t0 = time.time()
KD = load('kernels.pkl')
AUT = {(1, 3, 3): 4, (2, 2, 3): 4}
op = {}
th = {T: {} for T in ('d', 's1', 's0')}
for pqr, dat in KD.items():
    G, ext, lam = graph_data(pqr)
    W = colour_weight(G.nv, G.edges, ext)
    op = ce.lin((1, op), (Fr(1, AUT[pqr]), alt_cochain(dat['K0'], W)))
    for T in th:
        for v in range(G.nv):
            Wv = colour_weight(G.nv, G.edges, ext, {v: T})
            th[T] = ce.lin((1, th[T]), (Fr(1, AUT[pqr]), alt_cochain(dat['Kv'][v], Wv)))
    print("theta%s assembled  t=%.1f" % (pqr, time.time() - t0), flush=True)

obs = spaces.obstruction_cochains()
ref = ce.lin((Fr(-3, 240), obs[(1, 3, 3)]['F']), (Fr(2, 240), obs[(2, 2, 3)]['F']))
print("o' == (1/240)(-3 F_133 + 2 F_223) (formal gl_N):", is_zero(ce.lin((1, op), (-1, ref))))
print("o' alternating:", ce.is_alternating(op, 4), "; d o' = 0:", is_zero(ce.d_ce(op, 4)))
for T in th:
    print("theta'(%s): %d coordinates, N-powers %s, alternating %s, weights %s" % (
        T, len(th[T]), sorted({p for d in th[T].values() for p in d}), ce.is_alternating(th[T], 4), ce.weight_of(th[T])), flush=True)

res = {}
for sl in (False, True):
    nm = 'sl_N' if sl else 'gl_N'
    o = ce.restrict_sl(op) if sl else op
    for T in (('d',) if sl else ('d', 's1', 's0')):
        J = J_cochain(T, sl)
        thT = ce.restrict_sl(th[T]) if sl else th[T]
        dth = ce.d_ce(thT, 4, sl)
        bro = ce.bracket(J, 3, o, 4, sl=sl)
        ok = is_zero(ce.lin((1, dth), (-1, bro)))
        ctrl = is_zero(ce.lin((1, dth), (1, bro)))
        print("%s, J = %s: d theta' == {J, o'}: %s ; control d theta' == -{J, o'}: %s ; {J,o'} = 0: %s ; d theta' = 0: %s"
              % (nm, T, ok, ctrl, is_zero(bro), is_zero(dth)), flush=True)
        res[(nm, T)] = ok
save(dict(op=op, th=th, checks=res), 'theta.pkl')
print("time %.1fs" % (time.time() - t0))
