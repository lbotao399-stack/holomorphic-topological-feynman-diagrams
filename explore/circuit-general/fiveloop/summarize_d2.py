"""Markdown table of run_d2.jsonl (D=2, n=5 kernels)."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
cen = json.load(open(os.path.join(HERE, '..', 'atlas', 'census.json')))
gs = cen['data']['2,5']['graphs']
rows = {}
for fn in ('run_d2.jsonl', 'run_d2_w2.jsonl'):
    p = os.path.join(HERE, fn)
    if not os.path.exists(p):
        continue
    for l in open(p):
        d = json.loads(l)
        if d['index'] not in rows or rows[d['index']].get('K') is None:
            rows[d['index']] = d


def fmtK(K):
    if K is None:
        return '[open]'
    if not K:
        return '$0$'
    order = ['1', 'z3', 'z5', 'z7', 'z3*z3']
    parts = []
    for k in sorted(K, key=lambda s: order.index(s) if s in order else 99):
        v = K[k]
        lab = {'1': '', 'z3': r'\zeta_3', 'z5': r'\zeta_5', 'z7': r'\zeta_7', 'z3*z3': r'\zeta_3^2'}.get(k, k)
        num, _, den = v.partition('/')
        s = (r'\frac{%s}{%s}' % (num.lstrip('-'), den)) if den else num.lstrip('-')
        sign = '-' if v.startswith('-') else '+'
        parts.append((sign, s + lab))
    out = ''.join(('%s%s' % (sg if i or sg == '-' else '', t)) for i, (sg, t) in enumerate(parts))
    return '$' + out + '$'


print('| # | core | gen class | $K(\\lambda;0)$ | core vectors | s |')
print('|---|---|---|---|---|---|')
for i in sorted(rows):
    d = rows[i]
    g = gs[i]
    cls = 'SP' if g['sp'] else ('B' if g['B'] else 'R%s' % g['k_gen'])
    print('| %d | %s | %s | %s | %s | %s |' % (i, d['core'], cls, fmtK(d['K']), d.get('core_vectors', '—'), d['seconds']))
missing = [i for i in range(len(gs)) if i not in rows]
print('\nnot yet run:', missing, file=sys.stderr)
