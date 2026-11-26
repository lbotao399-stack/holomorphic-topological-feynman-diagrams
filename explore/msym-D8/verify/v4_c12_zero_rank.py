"""Rank of the C12 system restricted to the zero-derivative test operators (rows from v4_c12_rows.pkl)."""
import sys, pickle
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify')
from fractions import Fraction as Fr
from v4_c12 import solve, THETAS
rows = pickle.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/v4_c12_rows.pkl', 'rb'))
zero = [nm for nm in rows if 'd' not in nm.split('[')[0]]
for sel in ([n for n in zero if 'explorer' in n], [n for n in zero if 'explorer' not in n], zero):
    eqs = []
    for nm in sel:
        lhs, anti = rows[nm]
        mons = set(lhs)
        for a in anti.values():
            mons |= set(a)
        eqs += [[Fr(anti[G].get(m, 0)) for G in THETAS] + [-Fr(lhs.get(m, 0))] for m in mons]
    r, s, c = solve(eqs, 3)
    print('zero-derivative operators %s: %d equations, rank %d, consistent %s, solution %s' % (sel, len(eqs), r, c, None if s is None else [str(x) for x in s]))
