"""Original nine-edge C² prism to its rational + zeta(3) source polynomial.

No exported density or downstream algebraic result is used in the evaluation.
Invoke "python prism_series.py four 0" for four sources and "three 0" for
the distinct three-source choice. Source variables are formal in four-source
mode. Three independent parallel-conductance splitting fractions each have
uniform cube measure in the Jacobian change of variables.
"""
import json
from ast import literal_eval
import sys
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import comb, prod
from pathlib import Path
from random import Random

import networkx as nx
import sympy as sp
from flint import fmpz_mpoly_ctx
from graph_k4_series import wedge_determinant
from reduced_k4_series import beta_finite_parts, convergence, one_dimension

sys.dont_write_bytecode=True
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base))
from c2case import C2Case, PRISM, prism_case,prism3_case


def nonedge_case(choice):
    remaining=[0,1,2,4]
    center=remaining[choice]
    return C2Case('prism-nonedge',PRISM,(0,1,2,4),[(3,(0,2,4)),(5,(1,2,4))],
                  center,[v for v in remaining if v!=center])


def run(kind,center):
    case={'four':prism_case,'three':prism3_case,'nonedge':nonedge_case}[kind](center)
    xi=case.xi
    symbols=xi+[case.SA,case.SB]
    ctx=fmpz_mpoly_ctx.get(tuple(map(str,symbols)))
    gen=dict(zip(symbols,ctx.gens()))
    one=ctx.constant(1)
    @cache
    def rational(expr):
        expr=sp.sympify(expr)
        if expr.is_Integer:return ctx.constant(int(expr)),one
        if expr.is_Symbol:return gen[expr],one
        if expr.is_Pow:
            x,y=rational(expr.args[0])
            n=int(expr.args[1])
            return (x**n,y**n) if n>=0 else (y**(-n),x**(-n))
        n,d=ctx.constant(0 if expr.is_Add else 1),one
        for arg in expr.args:
            a,b=rational(arg)
            if expr.is_Add:n,d=n*b+a*d,d*b
            elif expr.is_Mul:n,d=n*a,d*b
            else:raise ValueError(expr)
        g=n.gcd(d)
        return n//g,d//g
    def poly(f):
        n,d=rational(f)
        return n/d
    q=poly(case.Q)
    num=[[poly(case.Q*v) for v in row] for row in case.y_expr]
    graph=nx.Graph(case.edges)
    edge_id={edge:i for i,edge in enumerate(case.edges)}
    for v in range(6):
        for i in range(2):
            assert sum((int(a==v)-int(b==v))*row[i] for (a,b),row in zip(case.edges,num))==q*poly(case.lamsym[v][i])
    ratios=[rational(t) for t in case.t_expr]
    cycles=nx.cycle_basis(graph)
    assert len(cycles)==4
    for cycle in cycles:
        ids=[]
        for a,b in zip(cycle,cycle[1:]+cycle[:1]):
            ids.append((edge_id.get((a,b),edge_id.get((b,a))),1 if (a,b) in edge_id else -1))
        common=ctx.constant(1)
        for e,_ in ids:
            denominator=ratios[e][1]
            common*=denominator//common.gcd(denominator)
        for j in range(2):
            assert sum(s*num[e][j]*ratios[e][0]*(common//ratios[e][1]) for e,s in ids)==0
    cotree=case.G.cotree()
    rows=[q]+[v for e in cotree for v in num[e]]
    columns=[rows]+[[r.derivative(i) for r in rows] for i in range(8)]
    det=wedge_determinant(columns,ctx.constant(0),ctx.constant(1))
    # Schur complement gives the current Jacobian det/Q^9. The chart is the
    # composition of bijective positive star/triangle and parallel split moves.
    sample=case.random_point(Random(619))
    source=[F(2),F(1)]
    point=sample+source
    powers=[[v**i for i in range(int(d)+1)] for v,d in zip(point,det.degrees())]
    numerators=defaultdict(lambda:defaultdict(F))
    value=F(0)
    for k,v in det.terms():
        k=tuple(map(int,k));v=int(v)
        value+=v*prod(p[i] for i,p in zip(k,powers))
        # Integral over three independent parallel split fractions.
        numerators[k[-2:]][k[3:8]]+=F(v,prod(i+1 for i in k[:3]))
    sign=case.g_ref(sample,*source)*F(case.Q.subs(dict(zip(xi,sample))))**9/value
    assert sign in (1,-1)
    coefs={}
    for ex,terms in sorted(numerators.items()):
        terms={k:v for k,v in terms.items() if v}
        if not terms:continue
        D=max(7,max(k[0]+k[1] for k in terms))  # ensure an auxiliary beta integral
        homog=defaultdict(F)
        for (g,h,*r),v in terms.items():
            n=D-g-h
            for i in range(n+1):
                for j in range(n-i+1):
                    homog[(g+i,h+j,n-i-j,*r)]+=v*comb(n,i)*comb(n-i,j)
        homog={k:v for k,v in homog.items() if v}
        margin=convergence(homog,D,9)
        val={k:str(sign*v) for k,v in beta_finite_parts(one_dimension(D,homog,9)).items()}
        assert sum(ex)<=4 and set(val)<={0,3}
        coefs[str(ex)]=dict(margin=margin,terms=val)
    return dict(mode=kind,center=case.Dc,cotree=cotree,determinant_monomials=len(det),
                orientation=int(sign),coefficients=coefs)


def compare(result):
    """Compare after the integration, using the original bracket-polynomial table."""
    values={literal_eval(k):{int(n):F(v) for n,v in data["terms"].items()}
            for k,data in result["coefficients"].items()}
    if result["mode"]=="three":
        assert values=={(0,0):{0:F(1,34560)}}
        return
    if result["mode"]=="nonedge":
        # Ward prediction for the signed-standard coefficient on (01)|(24).
        from prism_rung import bracket, projection
        projected=defaultdict(F)
        for sources,factor in projection(((0,1),(2,4))):
            a,b,c=[bracket(sources[i],sources[j]) for i,j in ((0,1),(0,2),(1,2))]
            for (sa,sb),pair in values.items():
                for n,v in pair.items():
                    projected[n]+=factor*v*a**(4-sa-sb)*b**sb*(-c)**sa
        assert dict(projected)=={0:-F(1,28800),3:F(1,9600)}
        return
    reference=json.loads((base/'results_nonSP_C2.json').read_text())['prism_4src']['table']
    for name,pair in reference.items():
        i,j,k=map(int,name[-3:])
        # <lambda0,lambda1>=1, <lambda0,lambda2>=SB, <lambda1,lambda2>=-SA.
        expected={weight:(-1)**k*F(v) for weight,v in zip((0,3),pair) if F(v)}
        assert values.get((k,j),{})==expected,(name,values.get((k,j)),expected)
    assert len(reference)==15


if __name__=='__main__':
    result=run(sys.argv[1],int(sys.argv[2]))
    compare(result)
    print(json.dumps(result,sort_keys=True))
