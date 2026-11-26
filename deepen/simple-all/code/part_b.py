"""Part (b): probe functionals, vacuum graphs and ladder moments.

Probe tensors (pairings pi of the slots 0..3):  kk_pi(x) = kappa(x_pi0,x_pi1) kappa(x_pi2,x_pi3),
ff_pi(x) = kappa([x_pi0,x_pi1],[x_pi2,x_pi3]),  pi in P3 = ((0,1,2,3), (0,2,1,3), (0,3,1,2)).
For a G x SL2-invariant 4-cochain X = sum_mu mu(p) X^mu (mu over QB):  v(X) = (<X^mu, kk_pi>, <X^mu, ff_pi>)_{mu,pi}.
Ladder moments (operators on g (x) g, channel (01)->(23)):
  m_k = tr(Om^k) = <Om^k, kk_(0213)>,  n_k = tr(P Om^k) = <Om^k, kk_(0312)>,  e_k = <Om^k, kk_(0123)> = d w1^k.
Killing normalization: a quantity with scaling exponent e (value(c kappa) = c^e value(kappa)) is multiplied by r^e,
r = Killing/kappa.  Exponents: <W,kk> -3, <W,ff> -4, m_k, n_k, e_k: -k, <X^mu,kk>: e_X - 2, <X^mu,ff>: e_X - 3,
e_F = -1, e_beta = 1, e_Gkk = 2.
Output: one JSON record per algebra (logs/b_data.json) and a printed table.
"""
import sys, time, json, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/code')
import numpy as np
from fractions import Fraction as Fr
import lie, t4

P3 = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]


def F2s(x):
    return str(x)


def build(key):
    if key.startswith('so') and key.endswith('tr'):
        return lie.so_tr(int(key[2:-2]))
    if key.startswith('sp') and key.endswith('tr'):
        return lie.sp_tr(int(key[2:-2]) // 2)
    if key.startswith('sl') and key.endswith('tr'):
        return lie.sl_tr(int(key[2:-2]))
    return lie.build(key)


def run(key):
    t0 = time.time()
    L = build(key)
    r, okr = L.killing_ratio()
    E = t4.Engine(L)
    W133, W223 = E.weights()
    Om, Om2, Om3 = E.Om, E.Om2, E.Om3
    Om4 = E.compose(Om2, Om2)
    kk_up = E.raise_all(E.kk)
    ff_up = E.raise_all(E.ff)
    KK = [t4.perm(kk_up, p) for p in P3]
    FF = [t4.perm(ff_up, p) for p in P3]
    KKl = [t4.perm(E.kk, p) for p in P3]
    FFl = [t4.perm(E.ff, p) for p in P3]
    rec = {'key': key, 'name': L.name, 'dim': L.dim, 'r': F2s(r), 's': L.s}
    K = lambda x, e: x * r ** e            # Killing normalization
    # ladder moments
    mom = {}
    for k, T in ((1, Om), (2, Om2), (3, Om3), (4, Om4)):
        mom['e%d' % k] = K(E.pair(T, KK[0]), -k)
        mom['m%d' % k] = K(E.pair(T, KK[1]), -k)
        mom['n%d' % k] = K(E.pair(T, KK[2]), -k)
    rec['moments'] = {k: F2s(v) for k, v in mom.items()}
    rec['n4_form'] = F2s(E.pair(Om4, KK[2]))          # n4 in the form used (no rescaling)
    # Gram
    gram = {}
    for i in range(3):
        for j in range(3):
            gram['kk%d.kk%d' % (i, j)] = K(E.pair(KKl[i], KK[j]), 0)
            gram['kk%d.ff%d' % (i, j)] = K(E.pair(KKl[i], FF[j]), -1)
            gram['ff%d.ff%d' % (i, j)] = K(E.pair(FFl[i], FF[j]), -2)
    rec['gram'] = {k: F2s(v) for k, v in gram.items()}
    # W numbers
    wn = {}
    for g, W in (('133', W133), ('223', W223)):
        for i in range(3):
            wn['W%s.kk%d' % (g, i)] = K(E.pair(W, KK[i]), -3)
            wn['W%s.ff%d' % (g, i)] = K(E.pair(W, FF[i]), -4)
    rec['W'] = {k: F2s(v) for k, v in wn.items()}
    # cochain probe vectors
    F133, F223 = t4.alt_bb(W133), t4.alt_bb(W223)
    B = t4.beta_t(E.ff)
    Gkk = t4.alt_bb(E.kk)
    vec = {}
    for name, X, eX in (('beta', B, 1), ('F133', F133, -1), ('F223', F223, -1), ('Gkk', Gkk, 2)):
        v = []
        for mu in range(6):
            for i in range(3):
                v.append(K(E.pair(X[mu], KK[i]), eX - 2))
            for i in range(3):
                v.append(K(E.pair(X[mu], FF[i]), eX - 3))
        vec[name] = v
    rec['v'] = {k: [F2s(x) for x in v] for k, v in vec.items()}
    rec['time'] = time.time() - t0
    return rec


if __name__ == "__main__":
    out = sys.argv[1]
    recs = []
    for key in sys.argv[2:]:
        rec = run(key)
        recs.append(rec)
        m = rec['moments']
        print("%-16s d %3d r %s : e1..e4 %s ; m1..m4 %s ; n1..n4 %s (%.1fs)" % (
            rec['name'], rec['dim'], rec['r'], [m['e%d' % k] for k in range(1, 5)], [m['m%d' % k] for k in range(1, 5)],
            [m['n%d' % k] for k in range(1, 5)], rec['time']), flush=True)
        print("     W: " + " ".join("%s=%s" % (k, v) for k, v in rec['W'].items()), flush=True)
        json.dump(recs, open(out, 'w'), indent=0)
