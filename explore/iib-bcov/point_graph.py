"""Direct heat and polyvector minors for marked BCOV graphs in five dimensions."""
from fractions import Fraction as F
from math import lcm,prod
from pathlib import Path
import sys
from flint import fmpq_mat,fmpz_mat
from marked_theta import MOMENTA,partitions,square

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'twoloop-universal'))
import twoloop


def solve_graph(edges,source,t,root):
    retained=[v for v in range(len(source)) if v!=root]
    incidence=fmpq_mat([[int(i==v)-int(j==v) for v in retained] for i,j in edges])
    diag=fmpq_mat([[str(F(1,x)) if i==j else 0 for j in range(len(t))]
                   for i,x in enumerate(t)])
    potential=(incidence.transpose()*diag*incidence).solve(fmpq_mat([source[v] for v in retained]))
    flow=diag*incidence*potential
    assert incidence.transpose()*flow==fmpq_mat([source[v] for v in retained])
    phi=dict(zip(retained,([F(str(x)) for x in row] for row in potential.tolist())))
    phi[root]=[F()]*len(source[root])
    return [[F(str(x)) for x in row] for row in flow.tolist()],phi


def raw(edges,y,root,at,u):
    vs=sorted(set(x for edge in edges for x in edge)-{root})
    def component(i,j,k):
        return [(int(i==v)-int(j==v))*int(k==h) for v in vs for h in range(5)]
    rows=[]
    denom=1
    for (i,j),current in zip(edges,y):
        scale=lcm(*(v.denominator for v in current))
        yy=[int(v*scale) for v in current]
        r=next(k for k in range(5) if yy[k])
        denom*=(-1)**r*scale*yy[r]**3
        base=component(i,j,r)
        for k in range(5):
            if k!=r:
                rows.append([yy[r]*a-yy[k]*b for a,b in zip(component(i,j,k),base)])
    def coefficient(row):
        return F(int(fmpz_mat([*rows,row]).det()),denom)
    odd=coefficient([a*int(v==at) for v in vs for a in u])
    facets=[]
    for (i,j),yy in zip(edges,y):
        r=next(k for k in range(5) if yy[k])
        facets.append(coefficient(component(i,j,r))/yy[r])
    return odd,facets


def theta(shape):
    nv,edges=twoloop.theta_graph(*shape)
    th=twoloop.Theta(nv,edges)
    times=[2+(i*i+3*i)%13 for i in range(len(edges))]
    sizes=[sum(times[j] for j,_ in th.paths[p]) for p in 'ABC']
    kirchhoff=sum(sizes[i]*sizes[(i+1)%3] for i in range(3))
    for groups in partitions(shape):
        if not square(groups):
            continue
        labels=sum(groups,())
        src=[(0,)*5,*[MOMENTA[i] for i in labels],(0,)*5]
        root=1+labels.index(7)
        at=1+labels.index(6)
        current,phi=solve_graph(edges,src,times,root)
        if not all(any(row) for row in current):
            continue
        q=twoloop.Q_form(th,src,twoloop.L_of_t(th,times,src))
        predicted=F(prod(sizes),kirchhoff**2)*q
        if not predicted:
            continue
        odd,facets=raw(edges,current,root,at,(1,1,1,1,0))
        for t,facet in zip(times,facets):
            assert abs(facet/t)==abs(predicted)
            assert odd*t==facet*sum(phi[at][:4])
        assert odd and sum(phi[at][:4])>0
        return
    raise AssertionError(shape)


def chain_graph(n):
    edges=[(i,(i+1)%6) for i in range(6)]
    left,right=2,3
    size=6
    for _ in range(1,n):
        path=[left,*range(size,size+4),right]
        edges+=list(zip(path,path[1:]))
        left,right=size+1,size+2
        size+=4
    return edges,size


def ladder(n):
    edges,size=chain_graph(n)
    degree=[sum(i in e for e in edges) for i in range(size)]
    src=[(0,)*5 for _ in degree]
    for j,i in enumerate((i for i in range(2,size) if degree[i]==2),1):
        src[i]=(1,j,j*j,j**3,j**4)
    src[1]=(0,1,2,3,4)
    src[0]=tuple(-sum(col) for col in zip(*src))
    t=[1+(i*i+7*i)%23 for i in range(len(edges))]
    flow,phi=solve_graph(edges,src,t,0)
    # Cycle basis from the directed boundary of each polygon.
    cycles=[[1 if i<6 else 0 for i in range(len(edges))]]
    for step in range(1,n):
        start=6+5*(step-1)
        old=2 if step==1 else start-3
        cycles.append([int(start<=i<start+5)-int(i==old) for i in range(len(edges))])
    for c in cycles:
        assert all(sum(x*(int(i==v)-int(j==v)) for x,(i,j) in zip(c,edges))==0
                   for v in range(size))
    minor=fmpq_mat([[str(c[j]*flow[j][i]) for j in range(1,len(edges))]
                    for i in range(5) for c in cycles])
    assert minor.det() and phi[1][0]>0
    return n


if __name__=='__main__':
    for shape in ((1,5,5),(2,4,5),(3,3,5),(3,4,4)):
        theta(shape)
    print('44 direct heat facets and four independent closed-index contractions checked.')
    print('Nonzero marked path-graph minors and strictly positive vector potentials at loops',
          [ladder(n) for n in range(1,9)])
