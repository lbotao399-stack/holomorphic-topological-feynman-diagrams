"""Fork B: 8d Q_2 theta kernels (no input jets) as products of 4x4 determinants.  Log: b_theta_forms.log"""
import itertools
import sympy as sp
import kern_h
import b_fastdet
import kernels_closed as kc
import rgraphs
from kern_h import GraphKernel

H = 4
for idx, g in enumerate(rgraphs.rooted_trees(4, 2)):
    gk = GraphKernel(H, g['nb'], g['edges'])
    e = sp.expand(kc.poly_of(gk))
    outs = gk.outputs
    L = {v: sp.Matrix([sp.Symbol('l%d_%d' % (v, a)) for a in range(H)]) for v in outs}
    const, facs = sp.factor_list(e)
    print('#%d %s |Aut| %d outputs %s: I = %s * prod of %d factors' % (idx + 1, g['edges'], g['aut'], outs, const, len(facs)), flush=True)
    # try each factor = det(L_a, L_b, L_c, sum_v x_v L_v) with {a,b,c} = 3 outputs, x integer in [-3,3]
    for f, mult in facs:
        found = None
        for trip in itertools.combinations(outs, 3):
            rest = [v for v in outs if v not in trip]
            for xs in itertools.product(range(-3, 4), repeat=len(rest)):
                if not any(xs):
                    continue
                M = sp.Matrix.hstack(L[trip[0]], L[trip[1]], L[trip[2]], sum((x * L[v] for x, v in zip(xs, rest)), sp.zeros(4, 1)))
                d = sp.expand(M.det())
                if d == 0:
                    continue
                r = sp.simplify(f / d)
                if r.is_number:
                    found = (r, trip, dict(zip(rest, xs)))
                    break
            if found:
                break
        print('    factor^%d = %s' % (mult, ('%s * det(l%d, l%d, l%d, %s)' % (found[0], found[1][0], found[1][1], found[1][2],
              ' + '.join('%d l%d' % (x, v) for v, x in found[2].items() if x))) if found else str(f)[:200]), flush=True)
