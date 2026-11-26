"""Integrate both identical hard charges and both soft jumps on a fixed retained metric graph."""
from argparse import ArgumentParser
from fractions import Fraction as F
from itertools import combinations,product

from flint import fmpq,fmpq_mat,fmpz_mat
from fourfold_mixed_soft_ward import V,hard_bipartite,integrated_soft_placements
from mixed_bare_ranks import WAVES
from point_graph import solve_graph


def integral(poly):
    return sum(c/F(i+1)/F(j+1) for (i,j),c in poly.items())


def quadratic(f):
    o,x,y,xy,m,n=f
    return {(0,0):o,(1,0):4*m-3*o-x,(0,1):4*n-3*o-y,
            (2,0):2*(x+o-2*m),(0,2):2*(y+o-2*n),(1,1):xy-x-y+o}


class HardPair:
    def __init__(self,edges,cycles,times,source,plus,minus):
        self.edges,self.cycles,self.times=edges,cycles,times
        self.source,self.plus,self.minus=source,plus,minus
        count=len(edges)
        assert count==13 and len(cycles)==4 and all(t>0 for t in times)
        self.columns=list(zip(*cycles))
        self.soft=list(combinations(range(4),2))
        flow,potential=solve_graph(edges,source,times,0)
        self.constant=[p[3] for p in flow]
        self.voltage=potential[plus][3]-potential[minus][3]
        # The remaining two unit charges balance the prescribed sources.
        assert tuple(map(sum,zip(*source)))==(0,0,0,-2)
        loads=[[int(i==j) for j in range(len(source))] for i in range(len(source))]
        self.response,phi=solve_graph(edges,loads,times,0)
        self.test=[x-y for x,y in zip(phi[plus],phi[minus])]
        matrix=[[sum(times[k]*a[k]*b[k] for k in range(count)) for b in cycles] for a in cycles]
        self.gram=[[2*(matrix[i][k]*matrix[j][l]-matrix[i][l]*matrix[j][k])/sum(times)**2
                    for k,l in self.soft] for i,j in self.soft]
        constant=[[1]*count]+[[fmpq(str(flow[e][mu]*r[e])) for e in range(count)]
                               for mu in range(3) for r in cycles]
        self.minors={}
        for i,j in combinations(range(count),2):
            fixed=[[v for h,v in enumerate(row) if h!=i and h!=j] for row in constant]
            ds=[]
            for a,b in self.soft:
                hard=[[*fixed[0],0,0]]
                for mu in range(3):
                    for k in range(4):
                        hard.append([*fixed[1+4*mu+k],V[mu]*int(a==k),V[mu]*int(b==k)])
                ds.append((-1)**(1+i+j)*F(str(fmpq_mat(hard).det())))
            if any(ds):self.minors[i,j]=ds

    def terms(self,e,f):
        cols=self.columns
        out=[]
        for (i,j),ds in self.minors.items():
            det=int(fmpz_mat([cols[k] for k in (i,j,e,f)]).det())
            if det:out.append((i,j,[det*v for v in ds]))
        return out

    def coefficients(self,e,f):
        terms=self.terms(e,f)
        if not terms:return {}
        u,v=self.edges[e];w,x=self.edges[f]
        tests=((0,0),(1,0),(0,1),(1,1),(F(1,2),0),(0,F(1,2)))
        evaluations=[]
        for s,t in tests:
            flow=[c+r[u]+r[w]+s*(r[v]-r[u])+t*(r[x]-r[w]) for c,r in zip(self.constant,self.response)]
            evaluations.append([sum(ds[a]*flow[i]*flow[j] for i,j,ds in terms) for a in range(6)])
        minors=[quadratic(column) for column in zip(*evaluations)]
        voltage={(0,0):self.voltage+self.test[u]+self.test[w],
                 (1,0):self.test[v]-self.test[u],(0,1):self.test[x]-self.test[w]}
        gram={k:[sum(g[a]*p[k] for a,p in enumerate(minors)) for g in self.gram] for k in minors[0]}
        square={}
        for (i,j),(l,m) in product(minors[0],repeat=2):
            key=i+l,j+m
            square[key]=square.get(key,0)+sum(p[i,j]*v for p,v in zip(minors,gram[l,m]))
        out={}
        for (i,j),c in square.items():
            for (l,m),v in voltage.items():
                key=i+l,j+m
                out[key]=out.get(key,0)+c*v
        return {k:v for k,v in out.items() if v}

    def ordinary(self,e,f,s,t):
        edges=[];times=[];parents=[];source=list(self.source)
        for h,((u,v),length) in enumerate(zip(self.edges,self.times)):
            if h in (e,f):
                x=s if h==e else t
                assert 0<x<1
                k=len(source);source.append(WAVES[3])
                edges.extend(((u,k),(k,v)));times.extend((length*x,length*(1-x)));parents.extend((h,h))
            else:
                edges.append((u,v));times.append(length);parents.append(h)
        flow,phi=solve_graph(edges,source,times,0)
        cycles=[[row[p] for p in parents] for row in self.cycles]
        return integrated_soft_placements(flow,cycles,times,V)*(phi[self.plus][3]-phi[self.minus][3])


def fixture(extended=False):
    edges,cycles,times=hard_bipartite((0,1,4,8))
    source=[(0,)*4]*6+[WAVES[k] for k in (0,1,2,5)]
    graph=HardPair(edges,cycles,times,source,6,7)
    total=F();nonzero=0;direct=0
    for e,f in combinations(range(len(edges)),2):
        polynomial=graph.coefficients(e,f)
        if not polynomial:
            if extended:assert graph.ordinary(e,f,F(2,5),F(3,7))==0
            continue
        assert max(i+j for i,j in polynomial)<=5
        nonzero+=1
        if direct<3 or extended:
            for s,t in ((F(2,5),F(3,7)),(F(4,9),F(5,11))):
                assert graph.ordinary(e,f,s,t)==sum(c*s**i*t**j for (i,j),c in polynomial.items())
            direct+=1
        if direct==1 or extended:
            # Six-point open Newton-Cotes on each interval; verified monomial by monomial.
            points=[(F(i+1,7),F(w,1440)) for i,w in enumerate((611,-453,562,562,-453,611))]
            assert all(sum(w*x**i for x,w in points)==F(1,i+1) for i in range(6))
            ordinary=sum(a*b*graph.ordinary(e,f,s,t) for (s,a),(t,b) in product(points,repeat=2))
            assert ordinary==integral(polynomial)
        total+=2*times[e]*times[f]*integral(polynomial)
    assert nonzero==74
    assert total==-F(629164240596633792744612123035776339502508524544,
                     1667494697978527475332008799247164349889438453144844885)
    print('Signed two-hard two-soft integral, nonzero unordered segment pairs:',total,nonzero,flush=True)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--all',action='store_true',help='check every surviving pair by direct segmented quadrature')
    fixture(parser.parse_args().all)
