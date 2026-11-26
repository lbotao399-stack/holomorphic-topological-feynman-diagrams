"""Mathematica script for the C^2 prism/K33 cases: density det(dy/dxi) of wedge_r d^2 y_{c_r} (Kruskal cotree),
bordered form det[[Q,dQ],[N,dN]]/Q^9, check at rational points vs brute, integrate the fractions over (0,1)."""
import sys, random, json
from fractions import Fraction as Fr
import sympy as sp
from c2case import prism_case, k33_case, prism3_case, k33_3_case

def mname(s):
    return str(s).replace('_', 'x')

def write(mk, dci, base, npts=3, seed=11):
    c = mk(dci)
    cot = c.G.cotree()
    ren = {s: sp.Symbol(mname(s)) for s in c.xi}
    ys = [c.y_expr[e][i].subs(ren) for e in cot for i in range(2)]
    xi = [ren[s] for s in c.xi]
    rng = random.Random(seed)
    pts, gref, srcs = [], [], []
    for _ in range(npts):
        pt = c.random_point(rng)
        SA, SB = Fr(rng.randint(-9, 9), 7), Fr(rng.randint(1, 9), 5)
        pts.append(pt); srcs.append((SA, SB)); gref.append(c.g_ref(pt, SA, SB))
    L = []
    L.append('ys = {%s};' % ', '.join(sp.mathematica_code(e) for e in ys))
    L.append('Qp = %s;' % sp.mathematica_code(c.Q.subs(ren)))
    L.append('vars = {%s};' % ', '.join(str(v) for v in xi))
    L.append('t0 = AbsoluteTime[];')
    L.append('Ns = Together[ys*Qp];')
    L.append('If[Union[Denominator /@ Ns] =!= {1}, Print["NOT POLYNOMIAL ", Union[Denominator /@ Ns]]; Exit[1]];')
    L.append('Ns = Expand[Ns];')
    L.append('bord = Prepend[Table[Prepend[D[Ns[[j]], #] & /@ vars, Ns[[j]]], {j, 8}], Prepend[D[Qp, #] & /@ vars, Qp]];')
    L.append('detB = Expand[Det[bord]];')
    L.append('Print["detB terms ", Length[detB], " time ", AbsoluteTime[] - t0];')
    for k, (pt, (SA, SB)) in enumerate(zip(pts, srcs)):
        rule = ', '.join('%s -> %d/%d' % (v, p.numerator, p.denominator) for v, p in zip(xi, pt))
        rule += ', SA -> %d/%d, SB -> %d/%d' % (SA.numerator, SA.denominator, SB.numerator, SB.denominator)
        L.append('Print["gsym%d ", InputForm[(detB/Qp^9) /. {%s}]];' % (k, rule))
    L.append('P = detB;')
    for f in xi[:len(c.fr)]:
        L.append('P = Expand[Integrate[P, {%s, 0, 1}]];' % f)
    L.append('Print["P terms ", Length[P], " time ", AbsoluteTime[] - t0];')
    L.append('Put[P, "%s_P.m"]; Put[Qp, "%s_Q.m"];' % (base, base))
    open('wl/' + base.split('/')[-1] + '.wl', 'w').write('\n'.join(L))
    json.dump({'name': c.name, 'Dc': c.Dc, 'tri': c.tri, 'cotree': cot, 'gref': [str(g) for g in gref]},
              open(base + '_ref.json', 'w'))
    print(c.name, 'Dc', c.Dc, 'cotree', cot, 'written')

if __name__ == "__main__":
    which = sys.argv[1]; dci = int(sys.argv[2])
    mk = {'prism': prism_case, 'K33': k33_case, 'prism3': prism3_case, 'K33_3': k33_3_case}[which]
    write(mk, dci, 'out/%s_D%d' % (which, dci))
