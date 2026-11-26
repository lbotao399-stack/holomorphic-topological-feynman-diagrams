"""All convergent primitive four-loop monomials at Kirchhoff pole order three."""
from argparse import ArgumentParser
from fractions import Fraction as F
from itertools import combinations
from math import factorial, prod
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import sympy as sp

from quartic_projective_ibp import AUTOMORPHISMS, K33, PERMS, PRISM, fold
from quartic_planar_periods import forest_data
from four_loop import WHEEL, edge_order, form, period_evaluator, read_vectors, shifted_polynomial, tagged
from wheel import WHEEL_PERMS


def cut_periods(name,rows,q,folder,executable):
    graph,group=(K33,AUTOMORPHISMS) if name=='K33' else (PRISM,PERMS)
    reductions=[]
    for strategy in (0,1):
        cut=strategy if name=='prism' else 0
        ordering=edge_order(cut)
        integrands=[]
        for i,p in enumerate(rows):
            lifted=p[:6]+(0,)+p[6:] if name=='wheel' else p
            images=[tuple(lifted[g[j]] for j in ordering) for g in group]
            levels=sorted({a[0] for a in images if a[0]},reverse=True)
            level=levels[min(strategy,len(levels)-1) if name=='K33' else 0]
            choice=max(a for a in images if a[0]==level)
            integrands.append((i,'*'.join(f'1/(p{j}.p{j})^{a}' for j,a in enumerate(choice[1:],1) if a)))
        vector=read_vectors(form(folder,executable,f'{name}_terminal_{strategy}',tagged(integrands),cut,edges=graph))
        assert set(vector)==set(range(len(rows)))
        reductions.append(vector)
    def valuation(p):
        return next(i for i,c in enumerate(shifted_polynomial(p,q-2)) if c)
    last=max([0]+[valuation(b)-valuation(a) for red in reductions for cols in red.values()
                  for a,b in cols.values()])
    evaluate=period_evaluator(folder,executable,last,q-2,graph)
    z3,z5=sp.symbols('z3 z5')
    result=[]
    for i,p in enumerate(rows):
        values=[evaluate(red[i]) for red in reductions]
        assert values[0]==values[1]
        norm=F(prod(factorial(a-1) for a in p),factorial(q-1))
        value=(values[0]*norm).expand()
        terms=(value.subs({z3:0,z5:0}),value.coeff(z3),value.coeff(z5))
        assert all(t.is_Rational for t in terms) and value==terms[0]+terms[1]*z3+terms[2]*z5
        result.append({'powers':p,'1,zeta3,zeta5':tuple(map(str,terms))})
    return result


def terminal(name,folder,executable):
    graph,group={'prism':(PRISM,PERMS),'K33':(K33,AUTOMORPHISMS),'wheel':(WHEEL,WHEEL_PERMS)}[name]
    count=len(graph)
    u=forest_data(tuple(graph),max(max(e) for e in graph)+1)[1]
    cotrees=[tuple(map(int,ex)) for ex,c in u.terms()]
    subsets=[(ids,min(sum(t[i] for i in ids) for t in cotrees)) for n in range(1,count)
             for ids in combinations(range(count),n)]
    rows=set()
    for doubled in combinations(range(count),12-count):
        p=tuple(1+int(i in doubled) for i in range(count))
        if all(sum(p[i] for i in ids)>3*h for ids,h in subsets):
            rows.add(fold(p,group))
    result=cut_periods(name,sorted(rows),3,folder,executable)
    print(name,len(rows),'complete convergent q=3 orbits',flush=True)
    return result


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('form',type=Path)
    parser.add_argument('forcer',type=Path)
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    with TemporaryDirectory(prefix='bcov-terminal-') as temporary:
        folder=Path(temporary)
        for item in ('forcer','forcer.h'):
            (folder/item).symlink_to((args.forcer/item).resolve(strict=True))
        certificate={name:terminal(name,folder,args.form.resolve(strict=True)) for name in ('prism','K33','wheel')}
    path=Path(__file__).with_name('quartic_terminal_periods.json')
    content=json.dumps(certificate,indent=2)+'\n'
    if args.save: path.write_text(content)
    else: assert path.read_text()==content
