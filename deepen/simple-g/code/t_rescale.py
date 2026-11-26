"""kappa -> c kappa (c = 2, 3): f -> c f, kappa^{-1} -> kappa^{-1}/c, W -> W/c (6 vertices, 7 edges),
F_G -> F_G/c, beta -> c beta; exact comparison of the b_theta evaluations."""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import numpy as np
from fractions import Fraction as Fr
import algebras as al, core
for key in sys.argv[1:]:
    L = al.build(key); T = core.Tensors(L); th = L.theta
    base = {G: core.F_eval(T, G, th) for G in ('133', '223')}
    b0 = core.beta_eval(L, th)
    for c in (2, 3):
        Lc = al.LA(L.name + ' x%d' % c, L.fu, c * L.kappa, '%d * form' % c, cartan=L.cartan)
        Tc = core.Tensors(Lc)
        ok = []
        for G in ('133', '223'):
            Fc = core.F_eval(Tc, G, th)
            # F_true = F_int / s^7 ; check F_true(c kappa) = F_true(kappa) / c  <=>  Fc * s^7 * c = F * sc^7
            ok.append(np.array_equal(Fc.astype(object) * (T.s7 * c), base[G].astype(object) * Tc.s7))
        okb = np.array_equal(core.beta_eval(Lc, th), c * b0)
        print("%-12s kappa -> %d kappa: W_133, W_223 (via F) scale by 1/%d: %s ; beta scales by %d: %s" % (L.name, c, c, ok, c, okb))
