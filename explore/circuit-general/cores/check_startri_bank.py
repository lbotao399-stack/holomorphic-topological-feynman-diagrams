import sys, json
sys.path.insert(0,'.')
from fractions import Fraction as F
import startri, vac
PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]
data = json.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/nonSP/verify/prism_forcer/periods.json'))
ok=bad=n=0
for ex, co in zip(data['exponents'], data['coefficients']):
    pw = [k+1 for k in ex]
    us = startri.unique_vertices(PRISM, pw, 6)
    if not us: continue
    n+=1
    ref = {k:F(v) for k,v in zip(('1','z3','z5'), co) if F(v)}
    ne, npw = startri.star_to_triangle(PRISM, pw, us[0], 6)
    try:
        P = vac.full_period(ne, npw, 6, mark=None)
    except vac.Stuck as e:
        print('stuck', e); continue
    if P == ref: ok+=1
    else:
        bad+=1
        if bad < 4: print('MISMATCH', pw, P, ref)
print('unique-vertex prism entries', n, 'ok', ok, 'bad', bad)
