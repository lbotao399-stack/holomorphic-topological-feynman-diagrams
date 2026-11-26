"""(e) {J_T, F_G} for T = d, s1, s0, G = theta(1,3,3), theta(2,2,3) (formal gl_N and sl_N), and a control
{D, r} with r a random combination of the C^4_(2,2) basis."""
import random
from common import *

obs = spaces.obstruction_cochains()
for sl in (False, True):
    for T in (('d', 's1', 's0') if not sl else ('d',)):
        J = J_cochain(T, sl)
        for G in [(1, 3, 3), (2, 2, 3)]:
            F = ce.restrict_sl(obs[G]['F']) if sl else obs[G]['F']
            print("%s: {J_%s, F_%d%d%d} == 0: %s" % ('sl_N' if sl else 'gl_N', T, *G, is_zero(ce.bracket(J, 3, F, 4, sl=sl))))
C4, k4 = spaces.basis(4, (2, 2), False)
rng = random.Random(1)
r = {}
for c in C4[:10]:
    r = ce.lin((1, r), (Fr(rng.randint(-5, 5)), spaces.to_np(c)))
print("control: {D, r} (r random in C^4_(2,2), gl_N) number of nonzero coordinates:", len(ce.clean(ce.bracket(J_cochain('d'), 3, r, 4))))
