"""Per-monomial certification for symbolic-source cases (blocks 'mon j k' in base_HLexport_src.txt):
value of [SA^j SB^k] I by series[0,1/2] + tanh-sinh / Gauss-Legendre(geometric) on [1/2,1]; PSLQ on each."""
import sys, re, json, os
from fractions import Fraction as Fr
import mpmath as mp
import xi_cert
from pslq_test import tests
from collect import sign_of

def blocks(path):
    lines = open(path).read().strip().split('\n')
    out = {}; cur = None; i = 0
    while i < len(lines):
        m = re.match(r'mon (\d+) (\d+)', lines[i])
        if m:
            cur = (int(m.group(1)), int(m.group(2))); out[cur] = []; i += 1; continue
        m = re.match(r'k (\d+) p (\d+)', lines[i])
        k, p = int(m.group(1)), int(m.group(2))
        cs = [Fr(c) for c in lines[i + 1].split()] if lines[i + 1].strip() else []
        out[cur].append((k, p, cs)); i += 2
    return out

def run_block(terms, sign, M=300):
    tmp = 'out/_tmp_block_HLexport.txt'
    with open(tmp, 'w') as f:
        for k, p, cs in terms:
            f.write('k %d p %d\n%s\n' % (k, p, ' '.join(str(c) for c in cs)))
    return xi_cert.run('out/_tmp_block', sign, M=M)

if __name__ == "__main__":
    base = 'out/' + sys.argv[1]
    s, ref = sign_of(base)
    res = {}
    for (j, k), terms in sorted(blocks(base + '_HLexport_src.txt').items()):
        if all(not cs or all(c == 0 for c in cs) for _, _, cs in terms):
            continue
        v1, v2, neg, tail, rem = run_block(terms, s)
        mp.mp.dps = 80
        diff = abs(v1 - v2)
        agree = -int(mp.floor(mp.log10(abs(diff / v1)))) if diff != 0 and v1 != 0 else 99
        rel = tests(v1, 80, min(agree - 5, 60)) if v1 != 0 else {}
        res['SA^%d SB^%d' % (j, k)] = {'value': mp.nstr(v1, 50), 'agree_digits': agree, 'neg_powers': str(neg),
                                      'pslq': {kk: (list(map(int, vv)) if vv else None) for kk, vv in rel.items()}}
        print(sys.argv[1], 'SA^%d SB^%d' % (j, k), mp.nstr(v1, 30), 'agree', agree, 'PSLQ', res['SA^%d SB^%d' % (j, k)]['pslq'])
    json.dump(res, open('results_cert_%s.json' % sys.argv[1], 'w'), indent=1)
