"""(f) SL2 structure of H^4 at total momentum degree 6 (formal N, trace basis, certified ranks).
(i)  dim H^3, H^4 at bidegrees (3,3), (4,2), (5,1), (6,0):  m_0 = dim H^4_(3,3) - dim H^4_(4,2) = dim (H^4_(3,3))^{SL2}.
(ii) SL2-invariant subcomplex at (3,3): C^k_inv = span alt(B tau), B = products of three brackets [p_i p_j] of the
     k momenta (sum p = 0), tau trace monomials; Z^4_inv, B^4_inv = d C^3_inv, H^4_inv.
Usage: python3 o7_sl2_structure.py sl|gl"""
import sys, time, itertools
from common import *

sl = sys.argv[1] == 'sl'
nm = 'sl_N' if sl else 'gl_N'
t0 = time.time()


def rank_of(cochs, keys):
    cochs = [c for c in cochs if c]
    if not cochs:
        return 0, None
    _, M = la.matrix_from(cochs, keys)
    c = la.certified_rank(M)
    return c['rank'], c


def H_at(deg):
    C2, _ = spaces.basis(2, deg, sl)
    C3, k3 = spaces.basis(3, deg, sl)
    C4, k4 = spaces.basis(4, deg, sl)
    k5 = spaces.keys_of(5, deg, sl)
    r2, _ = rank_of([spaces.const_part(ce.d_ce(spaces.to_np(c), 2, sl)) for c in C2], k3)
    r3, _ = rank_of([spaces.const_part(ce.d_ce(spaces.to_np(c), 3, sl)) for c in C3], k4)
    r4, c4 = rank_of([spaces.const_part(ce.d_ce(spaces.to_np(c), 4, sl)) for c in C4], k5)
    print("  %s %s: dim C^2,C^3,C^4 = %d,%d,%d; rank d2,d3,d4 = %d,%d,%d; H^3 = %d, H^4 = %d   (t=%.0f)"
          % (nm, deg, len(C2), len(C3), len(C4), r2, r3, r4, len(C3) - r3 - r2, len(C4) - r4 - r3, time.time() - t0), flush=True)
    return len(C4) - r4 - r3


print("(i) formal %s, total degree 6" % nm)
dims = {}
for deg in [(3, 3), (4, 2), (5, 1), (6, 0)]:
    dims[deg] = H_at(deg)
print("  m_0 = dim H^4_(3,3) - dim H^4_(4,2) = %d" % (dims[(3, 3)] - dims[(4, 2)]))


def inv_span(k):
    """alt(B tau), B = product of three brackets over the k momenta (raw, slot-major)."""
    pairs = list(itertools.combinations(range(k), 2))
    out = []
    for trip in itertools.combinations_with_replacement(pairs, 3):
        poly = {tuple([0] * (2 * k)): 1}
        for (i, j) in trip:
            poly = ce.poly_mul(poly, ce.bracket_poly(i, j, k))
        for tau in spaces.taus(k, sl):
            c = ce.alt_raw({(tau, m): {0: Fr(v)} for m, v in poly.items()}, k, sl=sl)
            c = ce.clean(c)
            if c:
                out.append(spaces.const_part(c))
    return out


print("(ii) SL2-invariant subcomplex at (3,3)")
deg = (3, 3)
k3, k4, k5 = spaces.keys_of(3, deg, sl), spaces.keys_of(4, deg, sl), spaces.keys_of(5, deg, sl)
S3 = inv_span(3)
S4 = inv_span(4)
r3s, c3s = rank_of(S3, k3)
r4s, c4s = rank_of(S4, k4)
C3i = [S3[j] for j in c3s['pivcols']]
C4i = [S4[j] for j in c4s['pivcols']]
rB, cBi = rank_of([spaces.const_part(ce.d_ce(spaces.to_np(c), 3, sl)) for c in C3i], k4)
dC4 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 4, sl)) for c in C4i]
rZ, cZi = rank_of(dC4, k5)
print("  dim C^3_inv = %d, dim C^4_inv = %d, rank d3|inv = %d, rank d4|inv = %d: dim Z^4_inv = %d, dim H^4_inv = %d   (t=%.0f)"
      % (len(C3i), len(C4i), rB, rZ, len(C4i) - rZ, len(C4i) - rZ - rB, time.time() - t0))
save(dict(C3i=C3i, C4i=C4i, dims=dims), 'sl2_structure_%s.pkl' % ('sl' if sl else 'gl'))
