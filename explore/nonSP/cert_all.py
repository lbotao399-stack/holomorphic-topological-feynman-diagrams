"""Certified high-precision values (two quadratures of the exact 1D xi-integrand) and PSLQ for the C^3 K4 cases."""
import json, os, sys
import mpmath as mp
import xi_cert
from collect import sign_of, result_of, parse
from pslq_test import tests

bases = sys.argv[1:] if len(sys.argv) > 1 else ['case0_D0', 'case0_D1', 'case1_D1', 'case2_D0', 'case2_D1', 'case3_D0', 'case3_D1', 'case5_D0', 'case5_D1']
out = {}
for b in bases:
    base = 'out/' + b
    if not os.path.exists(base + '_HLexport.txt'):
        continue
    s, ref = sign_of(base)
    v1, v2, neg, tail, rem = xi_cert.run(base, s)
    mp.mp.dps = 80
    diff = abs(v1 - v2)
    a, z = parse(result_of(base))
    exact = s * (mp.mpf(a.numerator) / a.denominator + mp.mpf(z.numerator) / z.denominator * mp.zeta(3))
    agree = -int(mp.floor(mp.log10(abs(diff / v1)))) if diff != 0 else 99
    rel = tests(v1, 80, min(agree - 5, 60))
    out[b] = {'name': ref['name'], 'Dc': ref['Dc'], 'sign': s, 'value_tanh_sinh': mp.nstr(v1, 65), 'value_gauss_legendre': mp.nstr(v2, 65),
              'abs_difference': mp.nstr(diff, 3), 'agreeing_significant_digits': agree,
              'closed_form_minus_numeric': mp.nstr(exact - v1, 3), 'negative_powers_left': str(neg),
              'pslq_tol_digits': min(agree - 5, 60), 'pslq': {k: (list(map(int, v)) if v else None) for k, v in rel.items()}}
    print(b, ref['name'], 'agree', agree, 'digits; closed form - numeric =', mp.nstr(exact - v1, 3), '; PSLQ', out[b]['pslq'])
json.dump(out, open('results_cert_C3.json', 'w'), indent=1)
