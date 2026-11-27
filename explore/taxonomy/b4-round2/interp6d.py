"""6d (1,0) sl2 + K fund, external u colours on flavour 1 only: each composition value is a polynomial in K of degree
<= number of closed flavour loops (<= 2 at order 2).  Interpolate from K = 1,2,3 to K = 16 (mod P); K = 4 checks the
degree bound.  usage: interp6d.py tag  (reads o2_6d10/E_f{K}_{tag}.json, writes o2_6d10/E_f16interp_{tag}.json)"""
import json
import sys

tag = sys.argv[1]
D = {K: json.load(open('o2_6d10/E_f%d_%s.json' % (K, tag))) for K in (1, 2, 3, 4)}
P = D[1]['P']
F = D[1]['F']


def lag(xs, x):
    out = []
    for j, xj in enumerate(xs):
        num, den = 1, 1
        for i, xi in enumerate(xs):
            if i != j:
                num *= (x - xi)
                den *= (xj - xi)
        out.append(num * pow(den % P, P - 2, P) % P)
    return out


c4, c16 = lag([1, 2, 3], 4), lag([1, 2, 3], 16)
E16 = {}
ok = True
for key in D[1]['E']:
    vals = [D[K]['E'][key] for K in (1, 2, 3)]
    pred4 = [sum(c * v[f] for c, v in zip(c4, vals)) % P for f in range(F)]
    good = (pred4 == D[4]['E'][key]) if key in D[4]['E'] else 'K=4 n/a'
    ok &= (good is True) if key in D[4]['E'] else True
    E16[key] = [sum(c * v[f] for c, v in zip(c16, vals)) % P for f in range(F)]
    nz = [K for K in (1, 2, 3, 4) if any(D[K]['E'].get(key, []))]
    print('%-12s K=4 predicted from K=1,2,3: %s   nonzero at K in %s' % (key, good, nz))
print('degree <= 2 in K on all compositions with K=4 data: %s (%d/%d compositions checked)' % (ok, sum(k in D[4]['E'] for k in D[1]['E']), len(D[1]['E'])))
json.dump(dict(P=P, F=F, E=E16), open('o2_6d10/E_f16interp_%s.json' % tag, 'w'), indent=0)
