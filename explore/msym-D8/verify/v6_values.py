"""The explorer's quoted QMC reference values (C5, C7, C8, C9) recomputed from the own bracket coefficients
(theta_weights.pkl, v3_thetas.py) and the own pentagon formula K = <0123>/24, at the explorer's QMC sources."""
import sys, pickle
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify')
from fractions import Fraction as Fr
from common import bracket
from v3_thetas import closed_form

lam = {1: [2, 3, 2, 2], 2: [3, 4, -2, -3], 3: [3, 2, 5, 4], 4: [-3, -4, 2, -1], 5: [-3, -4, 3, 5], 6: [-5, 4, 1, 2]}
lam_p = {0: [2, 3, 2, 2], 1: [3, 4, -2, -3], 2: [3, 2, 5, 4], 3: [-3, -4, 2, -1]}
print('pentagon: <0123>/24 = %s (explorer exact 9/2)' % (Fr(bracket(lam_p, (0, 1, 2, 3)), 24)))
tw = pickle.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/theta_weights.pkl', 'rb'))
quoted = {(1, 4, 4): Fr(-287899, 1440), (2, 3, 4): Fr(-16549, 80), (3, 3, 3): Fr(-2243657, 8640)}
for G, q in quoted.items():
    own = sum(c * bracket(lam, V) * bracket(lam, Vp) for (V, Vp), c in tw[G]['C0'].items())
    print('theta%s: own bracket coefficients K = %s, claimed closed form = %s, explorer exact = %s, all equal: %s'
          % (G, own, closed_form(G, lam), q, own == closed_form(G, lam) == q))

# exact polynomial identity: own W0 (monomials of K(lam;0)) == expansion of the claimed bracket closed forms
import itertools
from common import perm_sign
CF = {
    (1, 4, 4): [(Fr(-1, 11520) * a * b, S, T) for a, S in ((1, (1, 2, 3, 4)), (2, (1, 2, 3, 5)), (3, (1, 2, 3, 6)))
                for b, T in ((1, (1, 4, 5, 6)), (2, (2, 4, 5, 6)), (3, (3, 4, 5, 6)))],
    (2, 3, 4): [(Fr(-1, 17280) * a * b, S, T) for a, S in ((1, (1, 2, 3, 4)), (2, (1, 2, 3, 5)), (3, (1, 2, 3, 6)))
                for b, T in ((3, (1, 4, 5, 6)), (2, (2, 4, 5, 6)), (4, (3, 4, 5, 6)))],
    (3, 3, 3): [(Fr(c, 8640), S, T) for c, S, T in [
        (-1, (1, 2, 3, 4), (1, 3, 5, 6)), (-2, (1, 2, 3, 4), (1, 4, 5, 6)), (-2, (1, 2, 3, 4), (2, 3, 5, 6)),
        (-4, (1, 2, 3, 4), (2, 4, 5, 6)), (-1, (1, 2, 3, 5), (1, 4, 5, 6)), (-2, (1, 2, 3, 5), (2, 4, 5, 6)),
        (-1, (1, 2, 3, 5), (3, 4, 5, 6)), (-2, (1, 2, 3, 6), (1, 4, 5, 6)), (-4, (1, 2, 3, 6), (2, 4, 5, 6)),
        (-2, (1, 2, 3, 6), (3, 4, 5, 6)), (1, (1, 2, 4, 5), (1, 3, 5, 6)), (2, (1, 2, 4, 5), (2, 3, 5, 6)),
        (-2, (1, 2, 4, 5), (3, 4, 5, 6)), (2, (1, 2, 4, 6), (1, 3, 5, 6)), (4, (1, 2, 4, 6), (2, 3, 5, 6)),
        (-4, (1, 2, 4, 6), (3, 4, 5, 6))]],
}


def expand(terms, outputs):
    idx = {v: i for i, v in enumerate(outputs)}
    out = {}
    perms = [(p, perm_sign(p)) for p in itertools.permutations(range(4))]
    for cf, S, T in terms:
        for p1, s1 in perms:
            for p2, s2 in perms:
                m = [[0] * 4 for _ in outputs]
                for b in range(4):
                    m[idx[S[b]]][p1[b]] += 1
                    m[idx[T[b]]][p2[b]] += 1
                key = tuple(tuple(r) for r in m)
                out[key] = out.get(key, 0) + cf * s1 * s2
    return {k: v for k, v in out.items() if v}


for G in CF:
    W0 = tw[G]['W0']
    E = expand(CF[G], tw[G]['outputs'])
    print('theta%s: claimed closed form expanded in monomials (%d) == own K(lam;0) monomials (%d): %s' % (G, len(E), len(W0), E == W0))
