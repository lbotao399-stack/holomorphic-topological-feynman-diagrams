"""Exact rational-source current moments of the four-loop prism and K3,3."""
from collections import defaultdict
from fractions import Fraction as F
from math import factorial,lcm,prod
from pathlib import Path
import sys

from flint import fmpz_mpoly_ctx
import sympy as sp

from periods import CUTS, PERMS, PRISM, edge_order, form, read_vectors, tagged
from shift import K33, period_evaluator, shifted_polynomial
from check import AUTOMORPHISMS
from vk import G

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'k4_periods'))
from kernel import numerator, tree_flows


def current_polynomials(edges,source):
    if edges!=PRISM and edges!=K33:
        raise ValueError('Expected the project prism or K3,3 edge order.')
    source = [tuple(map(F,row)) for row in source]
    if len(source)!=6 or any(len(row)!=2 for row in source) or any(map(sum,zip(*source))):
        raise ValueError('Expected six conserved two-component momenta.')
    denominator = lcm(*(x.denominator for row in source for x in row))
    integer = [[int(x*denominator) for x in row] for row in source]
    graph = G(6,edges)
    trees = tree_flows(graph,integer)
    raw = numerator(graph,integer,trees)
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(9)))
    currents = [[ctx.from_dict({tuple(int(j in cut) for j in range(9)):v[e][i]
                                for cut,v in trees.items() if v[e][i]}) for i in range(2)] for e in range(9)]
    return ctx.from_dict(raw), currents,denominator


def evaluate(folder,executable,edges,source,insertions):
    """Evaluate a batch of monomial insertions of the same order. Empty batch is allowed."""
    if not insertions:
        return []
    order=len(insertions[0])
    if any(len(row)!=order for row in insertions):
        raise ValueError('Group insertions of equal total degree.')
    density,currents,denominator = current_polynomials(edges,source)
    perms = PERMS if edges==PRISM else AUTOMORPHISMS
    lookup = {}
    rows = []
    for row in insertions:
        polynomial = density
        for edge,component in row:
            if not 0<=edge<9 or component not in (0,1):
                raise ValueError('Unknown edge or current component.')
            polynomial = polynomial*currents[edge][component]
        out=defaultdict(int)
        for index,c in polynomial.terms():
            ex=tuple(map(int,index))
            if ex not in lookup:
                orbit={tuple(ex[j] for j in p) for p in perms}
                lookup.update(dict.fromkeys(orbit,min(orbit)))
            out[lookup[ex]]+=int(c)
        rows.append({ex:c for ex,c in out.items() if c})
    if not any(rows):
        return [sp.Integer(0)]*len(rows)
    norm=factorial(5+order)*denominator**(8+order)
    cuts = []
    for strategy in CUTS:
        cut = strategy if edges==PRISM else 0
        choose = min if edges==K33 and strategy==1 else max
        tagged_rows=[]
        for idx,row in enumerate(rows):
            if not row:
                continue
            terms=[]
            for ex,c in row.items():
                pattern=choose(tuple(ex[p[e]] for e in edge_order(cut)) for p in perms)
                den='*'.join(f'1/(p{j}.p{j})^{k+1}' for j,k in enumerate(pattern) if j)
                terms.append(f'+({F(c*prod(map(factorial,ex)),norm)})*'+den)
            tagged_rows.append((idx,''.join(terms)))
        text = form(folder,executable,f'insert{strategy}',tagged(tagged_rows),cut,edges=edges)
        cuts.append({} if text.strip()=='0' else read_vectors(text))
    shifts=4+order
    def valuation(s):
        return next(i for i,c in enumerate(shifted_polynomial(s,shifts)) if c)
    last=max([0]+[valuation(b)-valuation(a) for cut in cuts for row in cut.values() for a,b in row.values()])
    integral=period_evaluator(folder,executable,last,shifts,edges)
    values=[[integral(cut[i]) if i in cut else sp.Integer(0) for i in range(len(rows))] for cut in cuts]
    assert values[0]==values[1]
    return values[0]


if __name__=='__main__':
    from replay import run
    source=[(1,0),(0,1),(1,0),(-1,0),(-1,1),(0,-2)]
    def fixture(folder,executable):
        # Swapping 0 and 2 reverses the ordered graph and fixes the momentum assignment.
        assert evaluate(folder,executable,K33,source,[()])==[0]
        result=evaluate(folder,executable,K33,source,[[(e,j)] for e in (0,1) for j in (0,1)])
        coefficients = [(F(35039171,65318400),F(316373,226800),-F(32,15)),
                        (F(1820521,10886400),F(60881,113400),-F(493,630)),
                        (F(16083553,65318400),F(141,224),-F(407,420)),
                        (-F(4695731,13063680),-F(153887,151200),F(107,70))]
        z3,z5 = sp.symbols('z3 z5')
        assert result==[a+b*z3+c*z5 for a,b,c in coefficients]
        assert sp.expand(-3*result[1]+4*result[2]+result[3])==F(2695393,21772800)-F(16711,151200)*z3
        print('K3,3 first moments on 01 and 03',result)
    run(fixture)
