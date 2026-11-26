"""zeta(3) cancellation at n = 3: the zeta(3)-carrying graphs are #4, #6, #7, #8 (K4 subdivisions).
Formal gl(t) colour identities among their Sym parts (and full multilinear tensors), and
sum_Gamma [c_Gamma]_{zeta3} |Aut|^{-1} Sym W_Gamma = 0."""
import json
from fractions import Fraction as Fr
import colour as C
from assemble_n3 import aut
recs = {r['idx']: r for r in json.load(open('kernels_n3.json'))}
def W(i, sym=True):
    r = recs[i]
    rib = C.ribbon(r['nv'], [tuple(e) for e in r['edges']], r['legs'] + [r['nv'] - 1])
    return C.ribbon_sym(rib) if sym else rib
S = {i: W(i) for i in (4, 6, 7, 8)}
print('Sym W_6 == Sym W_4:', S[6] == S[4], '; Sym W_8 == Sym W_4:', S[8] == S[4], '; Sym W_7 == 2 Sym W_4:', S[7] == {k: 2 * v for k, v in S[4].items()})
z = {i: Fr(recs[i]['bare'].get('z3', '0')) / aut(recs[i]['nv'], [tuple(e) for e in recs[i]['edges']]) for i in (4, 6, 7, 8)}
print('zeta3 coefficients c/|Aut|:', {i: str(v) for i, v in z.items()}, '; weighted sum (units of Sym W_4):', z[4] + z[6] + 2 * z[7] + z[8])
