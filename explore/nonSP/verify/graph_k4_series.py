"""Exact graph-to-period replay. Requires SymPy and python-flint.

For the project's star-triangle chart verify *all* full-graph Kirchhoff equations
as polynomials. They uniquely fix each current in the positive domain. C d(t y)=0
and d y=C^T d ell in cotree coordinates give M d ell=-C diag(y) dt. Thus the
geometric density equals the current Jacobian up to a constant orientation sign.
One regular rational point fixes that sign from the original graph density.
An exterior product forms the bordered Jacobian anew over Z in nine variables.
Ordered-simplex monomial integration eliminates the four internal chain fractions.
The remaining five variables are integrated by reduced_k4_series.py directly.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import reduce
from itertools import product
from math import comb
from pathlib import Path
from random import Random
import sys

import networkx as nx
import sympy as sp
from flint import fmpz_mpoly_ctx
from reduced_k4_series import one_dimension, convergence, beta_finite_parts, read_homogeneous

sys.dont_write_bytecode=True
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base))
from k4exact import K4Case
from k4formula import formula_currents


def wedge_determinant(columns,zero,one):
    terms={0:one}
    for column in columns:
        next_terms={}
        for mask,minor in terms.items():
            for i,entry in enumerate(column):
                bit=1<<i
                if mask&bit or not entry:continue
                # Append row i to increasing row indices.
                term=minor*entry
                if (mask>>(i+1)).bit_count()%2:term=-term
                dst=mask|bit
                next_terms[dst]=next_terms.get(dst,zero)+term
        terms={key:v for key,v in next_terms.items() if v}
    return terms[(1<<len(columns))-1]


def run(which,center):
    case=K4Case(which,Dc_index=center)
    ys,Q=formula_currents(case)
    xi=case.xi
    ctx=fmpz_mpoly_ctx.get(tuple(map(str,xi)))
    gens=ctx.gens()
    def poly(f):
        return ctx.from_dict({m:int(v) for m,v in sp.Poly(f,*xi,domain=sp.ZZ).terms()})
    q=poly(Q)
    numerator=[[poly(sp.cancel(Q*v)) for v in row] for row in ys]
    graph=nx.Graph(case.edges)
    index={edge:i for i,edge in enumerate(case.edges)}
    # KCL and KVL identically in the chart, independently of formula_currents.
    for v in graph:
        for j in range(3):
            assert sum((int(a==v)-int(b==v))*row[j] for (a,b),row in zip(case.edges,numerator))==q*case.lam.get(v,[0]*3)[j]
    numden=[tuple(map(poly,sp.cancel(t).as_numer_denom())) for t in case.t_expr]
    cycles=nx.cycle_basis(graph)
    assert len(cycles)==3
    for cycle in cycles:
        signed=[]
        for a,b in zip(cycle,cycle[1:]+cycle[:1]):
            signed.append((index.get((a,b),index.get((b,a))),1 if (a,b) in index else -1))
        for j in range(3):
            common=ctx.constant(1)
            for e,_ in signed:common*=numden[e][1]
            assert sum(sign*numerator[e][j]*numden[e][0]*(common//numden[e][1]) for e,sign in signed)==0
    cotree=sorted(es[0] for path,es in case.chains if case.Dc not in (path[0],path[-1]))
    rows=[q]+[v for e in cotree for v in numerator[e]]
    columns=[rows]+[[row.derivative(i) for row in rows] for i in range(9)]
    # Multilinearity: remove all chain fractions from column zero, as N is affine there.
    columns[0]=[v-sum(gens[i]*columns[i+1][j] for i in range(4)) for j,v in enumerate(rows)]
    det=wedge_determinant(columns,ctx.constant(0),ctx.constant(1))
    chains=[len(es)-1 for path,es in case.chains if len(es)>1]
    unhomogeneous=defaultdict(F)
    sample=case.random_point(Random(718),den=31)
    # Evaluate at a fixed regular rational point termwise, without numerical operations.
    powers=[[v**n for n in range(d+1)] for v,d in zip(sample,det.degrees())]
    value=F(0)
    for ex,coef in det.terms():
        ex=tuple(map(int,ex))
        at=0;den=1
        for size in chains:
            accum=0
            for j,k in enumerate(reversed(ex[at:at+size]),1):
                accum+=k;den*=j+accum
            at+=size
        unhomogeneous[ex[4:]]+=F(int(coef),den)
        value+=int(coef)*reduce(lambda a,b:a*b,(p[n] for p,n in zip(powers,ex)),F(1))
    reference,_=case.g_ref(sample)
    qr=F(Q.subs(dict(zip(xi,sample))))
    sign=reference*qr**10/value
    assert sign in (1,-1)
    terms=defaultdict(F)
    D=max(k[0]+k[1] for k,v in unhomogeneous.items() if v)
    for (g,h,*r),c in unhomogeneous.items():
        if not c:continue
        n=D-g-h
        for i in range(n+1):
            for j in range(n-i+1):
                terms[(g+i,h+j,n-i-j,*r)]+=c*comb(n,i)*comb(n-i,j)
    terms={k:v for k,v in terms.items() if v}
    # Agreement is an additional check; the graph-generated polynomial is the input below.
    assert (D,terms)==read_homogeneous(base/'out'/f'case{which}_D{center}_P.m')
    margin=convergence(terms,D)
    answer={k:v*sign for k,v in beta_finite_parts(one_dimension(D,terms)).items()}
    print(case.name,'star center',case.Dc,'exact bordered terms',len(det),'chain density terms',len(terms),
          'orientation',sign,'convergence margin',margin,'coefficients 1,zeta(k)',answer,flush=True)


if __name__=='__main__':
    for name in sys.argv[1:] or ['0:0']:
        run(*map(int,name.split(':')))
