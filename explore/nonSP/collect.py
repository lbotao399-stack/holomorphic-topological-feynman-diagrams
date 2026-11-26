"""Collect exact results c_Gamma = s * int P/Q^10 (s from the pointwise comparison g_ref/g_sym)."""
import json, re, os, glob
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 50

def sign_of(base):
    ref = json.load(open(base + '_ref.json'))
    log = open(base + '.log').read()
    gs = {int(m.group(1)): Fr(m.group(2).replace(' ', '')) for m in re.finditer(r'gsym(\d+) InputForm\[([-0-9/ ]+)\]', log)}
    ratios = [Fr(g) / gs[k] for k, g in enumerate(ref['gref'])]
    assert len(set(ratios)) == 1 and abs(ratios[0]) == 1, ratios
    return int(ratios[0]), ref

def result_of(base):
    txt = open(base + '_int.log').read()
    m = re.search(r'RESULT c = InputForm\[(.*?)\]  N', txt, re.S)
    return m.group(1).replace('\n', '') if m else None

def parse(expr):
    """a + b*Zeta[3] forms -> (a, b) rationals"""
    e = expr.replace(' ', '')
    if e == '0':
        return Fr(0), Fr(0)
    m = re.fullmatch(r'\((-?\d+)([+-]\d+)\*Zeta\[3\]\)/(\d+)', e)
    if m:
        a, b, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
        return Fr(a, d), Fr(b, d)
    m = re.fullmatch(r'(-?\d+)/(\d+)', e)
    if m:
        return Fr(int(m.group(1)), int(m.group(2))), Fr(0)
    m = re.fullmatch(r'(-?\d+)', e)
    if m:
        return Fr(int(m.group(1))), Fr(0)
    return None

if __name__ == "__main__":
    rows = {}
    for w in range(6):
        for d in (0, 1):
            base = 'out/case%d_D%d' % (w, d)
            if not (os.path.exists(base + '.log') and os.path.exists(base + '_int.log')):
                continue
            r = result_of(base)
            if r is None:
                continue
            s, ref = sign_of(base)
            ab = parse(r)
            a, b = (s * ab[0], s * ab[1]) if ab else (None, None)
            val = mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator * mp.zeta(3) if ab else None
            rows.setdefault(ref['name'], []).append({'Dc': ref['Dc'], 'tri': ref['tri'], 'sign': s, 'integral': r,
                                                    'c_rational_part': str(a), 'c_zeta3_part': str(b),
                                                    'c_numeric': mp.nstr(val, 40) if val is not None else None})
    for k, v in rows.items():
        print(k)
        for x in v:
            print('   Dc', x['Dc'], 'tri', x['tri'], 's', x['sign'], ' c =', x['c_rational_part'], '+', x['c_zeta3_part'], '* zeta(3) =', x['c_numeric'])
    json.dump(rows, open('results_exact_C3.json', 'w'), indent=1)
