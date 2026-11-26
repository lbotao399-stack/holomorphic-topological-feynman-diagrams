"""Verify the complete traceless four-input QME at order four before rank specialization."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations
import json
import sys

from four_point import HERE,columns,kernel
sys.path.insert(0,str(HERE.parent))
from order_four import ce,cubic,marked_thetas,ribbon
from theta_moyal import alternating,serialize,to_ce


def bare_terms(graphs):
    output = [defaultdict(F),defaultdict(F)]
    for index,graph in enumerate(graphs):
        edges=[edge for path in graph['chains'] for edge in zip(path,path[1:])]
        colours=ribbon.colour_weight(10,edges,list(range(6,10)))
        for trace,ranks in colours.items():
            if ce.has_singleton(trace):
                continue
            for rank,c in ranks.items():
                for monomial,pair in columns(index).items():
                    for j,coefficient in enumerate(pair):
                        output[j][trace,rank,monomial]+=coefficient*c/graph['aut']
    return [to_ce(serialize(alternating(raw))) for raw in output]


def marked_terms():
    output=[defaultdict(F),defaultdict(F)]
    for shape in ((1,3,3),(2,2,3)):
        for colours,integrated,power,coefficient,slots in marked_thetas(shape):
            for trace,ranks in colours.items():
                if ce.has_singleton(trace):
                    continue
                for k,c in ranks.items():
                    for monomial,v in integrated.items():
                        output[slots[0]][trace,k+power,monomial]+=coefficient*c*v
    return [to_ce(serialize(alternating(raw))) for raw in output]


def verify():
    graphs=json.loads((HERE.parent/'four-loop-graphs.json').read_text())
    bare0,bare2=bare_terms(graphs)
    single,double=marked_terms()
    df=ce.d_ce(cubic(4,'f'),3,sl=True)
    tree13=ce.bracket(cubic(1,'d'),3,cubic(3,'d'),3,sl=True)
    tree22=ce.bracket(cubic(2,'f'),3,cubic(2,'f'),3,sl=True)
    a2={2:F(-1,4)}
    k2={0:F(-1,2),2:F(-3,20)}
    k3={0:F(-1,6),2:F(-1,15)}
    total=ce.lin((ce.np_mul(a2,k3),tree13),(ce.np_scale(ce.np_mul(k2,k2),F(1,2)),tree22),
                 (k2,single),(a2,double),(1,bare0))
    coeff0={0:F(1,24),2:F(1102559,5443200),4:F(4216,212625)}
    coeff2={2:F(-35701,1935360),4:F(-130727,92897280)}
    for lhs,coeff in ((total,coeff0),(bare2,coeff2)):
        assert not ce.clean(ce.lin((1,lhs),(ce.np_scale(coeff,-1),df)))
    print('complete formal-rank QME coefficients',len(total),len(bare2),'CE coordinates')


if __name__=='__main__':
    verify()
