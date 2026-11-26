"""(g) d([p1p2]^3 s0) = 0 (formal gl_N)."""
from common import *
u = ce.bracket_poly(0, 1, 3); u3 = ce.poly_mul(ce.poly_mul(u, u), u)
s0 = {((0,), (1,), (2,)): {0: Fr(1)}}
c = ce.canon_raw(ce.raw_from(u3, s0), 3)
print('d([p1p2]^3 s0) == 0:', is_zero(ce.d_ce(c, 3)), '; [p1p2]^3 s0 != 0:', not is_zero(c))
