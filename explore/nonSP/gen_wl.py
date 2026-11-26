"""Generate the Mathematica script for one K4 subdivision: symbolic density of the 9-form
wedge_r d^3 y_{c_r} (c_r = first edges of the chains on the triangle K4 edges), bordered-determinant form
   det(d y/d xi) = det[[Q, dQ],[N, dN]] / Q^10,   y = N/Q,
evaluated at Python's rational test points, then integrated exactly over the chain variables sigma."""
import sys, random, json
from fractions import Fraction as Fr
import sympy as sp
from k4exact import K4Case
from k4formula import formula_currents

def mname(s):
    return str(s).replace('_', 'x')

def build(which, Dc_index=0, npts=3, seed=7):
    c = K4Case(which, Dc_index=Dc_index)
    y, Q = formula_currents(c)
    # cotree: first edges of the chains on the triangle K4 edges
    cot = []
    for ci, (path, es) in enumerate(c.chains):
        if c.Dc not in (path[0], path[-1]):
            cot.append(es[0])
    cot.sort()
    ren = {s: sp.Symbol(mname(s)) for s in c.xi}
    ys = [y[e][i].subs(ren) for e in cot for i in range(3)]
    xi = [ren[s] for s in c.xi]
    rng = random.Random(seed)
    pts, gref = [], []
    for _ in range(npts):
        pt = c.random_point(rng)
        g, t = c.g_ref(pt)
        pts.append(pt); gref.append(g)
    # sigma integration domains: per long chain, ordered simplex 1 > s1 > s2 > ... > 0
    doms = []
    for ci, (path, es) in enumerate(c.chains):
        if len(es) >= 2:
            doms.append([mname(c.sig[(ci, k)]) for k in range(1, len(es))])
    return c, cot, ys, Q.subs(ren), xi, pts, gref, doms

def wl_script(which, Dc_index, outbase):
    c, cot, ys, Q, xi, pts, gref, doms = build(which, Dc_index)
    L = []
    L.append('ys = {%s};' % ', '.join(sp.mathematica_code(e) for e in ys))
    L.append('Qp = %s;' % sp.mathematica_code(Q))
    L.append('vars = {%s};' % ', '.join(str(v) for v in xi))
    L.append('t0 = AbsoluteTime[];')
    L.append('Ns = Together[ys*Qp];')
    L.append('If[Union[Denominator /@ Ns] =!= {1}, Print["NOT POLYNOMIAL"]; Exit[1]];')
    L.append('Ns = Expand[Ns];')
    L.append('bord = Prepend[Table[Prepend[D[Ns[[j]], #] & /@ vars, Ns[[j]]], {j, 9}], Prepend[D[Qp, #] & /@ vars, Qp]];')
    L.append('detB = Expand[Det[bord]];')
    L.append('Print["detB terms ", Length[detB], " time ", AbsoluteTime[] - t0];')
    for k, pt in enumerate(pts):
        rule = ', '.join('%s -> %s' % (v, '%d/%d' % (p.numerator, p.denominator)) for v, p in zip(xi, pt))
        L.append('Print["gsym%d ", InputForm[(detB/Qp^10) /. {%s}]];' % (k, rule))
    # sigma integration of the polynomial detB
    L.append('sigInt[poly_, dom_List] := Module[{p = poly, d = dom, v},')
    L.append('  Do[p = Integrate[p, {d[[k]], 0, If[k == 1, 1, d[[k - 1]]]}], {k, Length[d], 1, -1}]; p];')
    L.append('P = detB;')
    for d in doms:
        L.append('P = Expand[sigInt[P, {%s}]];' % ', '.join(d))
    L.append('Print["P terms ", Length[P], " time ", AbsoluteTime[] - t0];')
    L.append('Print["degrees gA gB rA rB rC: ", Exponent[P, #] & /@ {gA, gB, rA, rB, rC}];')
    L.append('Put[P, "%s_P.m"];' % outbase)
    L.append('Put[Qp, "%s_Q.m"];' % outbase)
    return '\n'.join(L), c, pts, gref

if __name__ == "__main__":
    which = int(sys.argv[1]); dci = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    base = 'out/case%d_D%d' % (which, dci)
    s, c, pts, gref = wl_script(which, dci, base)
    open('wl/case%d_D%d.wl' % (which, dci), 'w').write(s)
    json.dump({'name': c.name, 'Dc': c.Dc, 'tri': c.tri, 'pts': [[str(p) for p in pt] for pt in pts],
               'gref': [str(g) for g in gref]}, open(base + '_ref.json', 'w'))
    print(c.name, 'Dc', c.Dc, 'tri', c.tri, 'written')
