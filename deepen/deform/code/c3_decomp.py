"""C^3(1,1) = [p1p2] (x) (Sym^3 g*)^g  +  E_f  +  E_mix (gl_N); d E_f, d E_mix independent and nonzero.
E_f = alt(p1^a p1^b tr(1[23])) (S_3-invariant momentum factor times f), E_mix = alt((p1^a p2^b + p1^b p2^a) tr(1) tr(23))."""
import sys
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
import ce, la, spaces
u = ce.bracket_poly(0, 1, 3)
for sl in (False, True):
    T = {"D": {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}},
         "S1": {((0,), (1, 2)): {0: Fr(1)}, ((0, 2), (1,)): {0: Fr(1)}, ((0, 1), (2,)): {0: Fr(1)}},
         "S0": {((0,), (1,), (2,)): {0: Fr(1)}}}
    if sl:
        T = {"D": T["D"]}
    J = {k: ce.canon_raw(ce.raw_from(u, t), 3, sl=sl) for k, t in T.items()}
    mono = {ce.monomial_full(3, {(0, 0): 1, (0, 1): 1}): 1}      # p1^a p1^b
    Ef = ce.alt_raw(ce.raw_from(mono, {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(-1)}}), 3, sl=sl)
    extra = {"E_f": Ef}
    if not sl:
        extra["E_mix"] = ce.alt_raw(ce.raw_from({ce.monomial_full(3, {(0, 0): 1, (1, 1): 1}): 1, ce.monomial_full(3, {(0, 1): 1, (1, 0): 1}): 1}, {((0,), (1, 2)): {0: Fr(1)}}), 3)
    allc = list(J.values()) + list(extra.values())
    keys, M = la.matrix_from([spaces.const_part(c) for c in allc])
    r = la.certified_rank(M)['rank']
    C3, _ = spaces.basis(3, (1, 1), sl)
    dE = [spaces.const_part(ce.d_ce(c, 3, sl)) for c in extra.values()]
    keys4, M4 = la.matrix_from(dE)
    r4 = la.certified_rank(M4)['rank']
    print("%s: rank{%s} = %d = dim C^3(1,1) = %d: %s;  rank{d %s} = %d" % (
        "sl_N" if sl else "gl_N", ", ".join(list(J) + list(extra)), r, len(C3), r == len(C3), ", d ".join(extra), r4))
    # E_f momentum factor: alt(p1^a p1^b f) = (sum_i p_i^a p_i^b) * 2 f ?
    sym = {}
    for i in range(3):
        sym[ce.monomial_full(3, {(i, 0): 1, (i, 1): 1})] = 1
    chk = ce.canon_raw(ce.raw_from(sym, {((0, 1, 2),): {0: Fr(2)}, ((0, 2, 1),): {0: Fr(-2)}}), 3, sl=sl)
    print("   E_f = 2 (p1^a p1^b + p2^a p2^b + p3^a p3^b) tr(1[23]):", ce.clean(ce.lin((1, Ef), (-1, chk))) == {})
