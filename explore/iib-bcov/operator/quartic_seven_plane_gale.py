"""All retained-label permutations from one canonical current nullspace."""
from itertools import combinations,combinations_with_replacement,permutations
from fractions import Fraction as F
from math import prod

from random import Random
from flint import fmpq_mat,fmpz_mat,nmod_mat
from quartic_seven_plane_lift import NORMAL,fixture,cubic_cores,G,solve_graph
from square_three import tree_flows
from current_minor import numerator


def cubic_nodes():
    u=tuple(tuple(int(i==j) for i in range(4)) for j in range(4))
    return [*u,*[tuple(x+sign*y for x,y in zip(u[i],u[j])) for i,j in combinations(range(4),2) for sign in (1,-1)],
            *[tuple(map(sum,zip(*(u[i] for i in inds)))) for inds in combinations(range(4),3)]]


def monomials(x):
    return [prod(x[i] for i in ids) for ids in combinations_with_replacement(range(4),3)]


def gale(q):
    return [(-1)**i*int(fmpz_mat(q[:i]+q[i+1:4]).det()) for i in range(4)]


def canonical(edges,t):
    cycles=G(11,edges).C
    unit=[tuple(int(i==j)-int(i==4) for j in range(4)) for i in range(5)]
    flow,phi=solve_graph(edges,[(0,0,0,0)]*6+unit,t,0)
    gram=fmpz_mat([[sum(a*b*h for a,b,h in zip(x,y,t)) for x in cycles] for y in cycles])
    u=int(gram.det())
    y=[[int(a*u) for a in line] for line in flow]
    # Row order is cycle, canonical source. Kirchhoff permits deleting the first column.
    a=[[cs[e]*y[e][j] for e in range(1,14)] for cs in cycles for j in range(4)]
    r,n=fmpq_mat(a).transpose().rref()
    if n<13:
        return None
    piv=[next(i for i in range(16) if r[j,i]) for j in range(13)]
    omit=[i for i in range(16) if i not in piv]
    base=(-1)**sum(omit)*int(fmpz_mat([a[i] for i in piv]).det())
    null=[[F() for _ in omit] for i in range(16)]
    for j,k in enumerate(omit):
        null[k][j]=F(1)
        for h,i in enumerate(piv):
            null[i][j]=-F(str(r[h,k]))
    return base,null,u,cycles,y,phi


def cubic(null,m):
    x=[[sum(m[a]*null[4*r+a][j] for a in range(4)) for j in range(3)] for r in range(4)]
    def det(a,b,c):
        return a[0]*(b[1]*c[2]-c[1]*b[2])-a[1]*(b[0]*c[2]-c[0]*b[2])+a[2]*(b[0]*c[1]-c[0]*b[1])
    return [(-1)**r*det(*[p for j,p in enumerate(x) if j!=r]) for r in range(4)]


def replay(edges,t,normal=NORMAL):
    data=canonical(edges,t)
    assert data
    base,null,u,cycles,y,phi=data
    z=(1,-1,0,0,0)
    nodes=cubic_nodes()
    fit=fmpq_mat([monomials(x) for x in nodes]).inv()*fmpq_mat([[str(base*a) for a in cubic(null,x)] for x in nodes])
    m=gale(normal)+[0]
    determinant=sum(a*b for a,b in zip(m,z))
    for order in permutations(range(5)):
        q=[normal[i] for i in order]
        d=gale(q); sz=[z[i] for i in order]
        expected=[int(base*sum(a*b for a,b in zip(d,sz))*v) for v in cubic(null,d)]
        transformed=[m[i]-m[order[-1]] for i in order[:4]]
        fast=fmpq_mat([monomials(transformed)])*fit*determinant
        assert expected==[int(x) for x in fast.tolist()[0]]
        cy=[[sum(a*b for a,b in zip(row,x)) for row in y] for x in zip(*q[:4],strict=True)]
        zz=[sum(a*b for a,b in zip(row,sz[:4])) for row in y]
        norms=[[cs[e]*c[e] for e in range(1,14)] for c in cy for cs in cycles]
        direct=[int(fmpz_mat([*norms,[cs[e]*zz[e] for e in range(1,14)]]).det()) for cs in cycles]
        assert expected==direct
    return True


def canonical_sources(m):
    # Choose unit coordinate at the pivot to give integral normal waves and lift.
    pivot=next(i for i,v in enumerate(m) if abs(v)==1)
    rest=[i for i in range(4) if i!=pivot]
    q=[tuple(int(i==j) for j in rest) if i!=pivot else tuple(-m[j]//m[pivot] for j in rest) for i in range(4)]
    d=gale(q)
    ratio=d[pivot]//m[pivot]
    assert d==[ratio*x for x in m] and abs(ratio)==1
    p=[(int(i==pivot),*v) for i,v in enumerate(q)]
    p.append(tuple(-sum(col) for col in zip(*p)))
    return [(0,)*4]*6+p,d[pivot]*ratio**3


def symbolic_support(paths):
    g=G(11,fixture(cubic_cores()[0],paths))
    nodes=cubic_nodes()
    inv=fmpq_mat([monomials(x) for x in nodes]).inv().tolist()
    vals=[]
    for x,y,c,d in nodes:
        source,scale=canonical_sources((x+c,y+c,c,d))
        flow=next(iter(tree_flows(g,source).values()))
        vals.append([scale*(-1)**r*numerator(g,source,flow,[(0,j) for j in range(4) if j!=r]) for r in range(4)])
    ctx=vals[0][0].context()
    support=[]
    words=list(combinations_with_replacement("xycd",3))
    for r in range(4):
        coeff=[sum((vals[j][r]*int(6*v) for j,v in enumerate(row) if v),ctx.constant(0)) for row in inv]
        support.append(["".join(k) for k,v in zip(words,coeff) if v])
    return support


def resistance_table(q,allowed,target=(0,1),axis=(0,0,1)):
    m=gale(q)+[0]
    pairs=list(combinations(range(5),2))
    tab=[[0]*len(allowed) for _ in pairs]
    for order in permutations(range(5)):
        s=[m[order[0]]-m[order[2]],m[order[1]]-m[order[2]],m[order[2]]-m[order[4]],m[order[3]]-m[order[4]]]
        mono=[prod(s[i] for i in k) for k in allowed]
        charge=[sum(x*y for x,y in zip(q[j],axis)) for j in order]
        label=[int(j==target[0])-int(j==target[1]) for j in order]
        for h,(i,j) in enumerate(pairs):
            v=label[i]*charge[j]+label[j]*charge[i]
            if v:
                tab[h]=[a+v*b for a,b in zip(tab[h],mono)]
    for i in range(5):
        assert all(v==0 for v in map(sum,zip(*(row for row,pair in zip(tab,pairs) if i in pair))))
    return fmpz_mat(tab)


def primitive_rank(core):
    rng=Random(1)
    prime=1000003
    rows=[]
    edges=fixture(core,(0,2,3,6,8))
    for _ in range(5):
        data=canonical(edges,[rng.randrange(1,30) for _ in range(14)])
        columns=[[v.numerator*pow(v.denominator,-1,prime)%prime for v in cubic(data[1],x)] for x in cubic_nodes()]
        rows.extend(zip(*columns))
    return int(nmod_mat(rows,prime).det())


if __name__=='__main__':
    assert fmpz_mat([monomials(x) for x in cubic_nodes()]).rank()==20
    m=gale(NORMAL)+[0]
    assert fmpz_mat([monomials([m[o[i]]-m[o[4]] for i in range(4)]) for o in permutations(range(5))]).rank()==20
    core=cubic_cores()[0]
    for paths in ((0,0,1,7,8),(0,1,4,7,8),(0,0,1,7,7)):
        assert replay(fixture(core,paths),[2+i*i%17 for i in range(14)])
    assert [primitive_rank(c) for c in cubic_cores()[-2:]]==[990338,896731]
    all_words=list(combinations_with_replacement(range(4),6))
    for cut,expected in ((0,3),(2,5),(6,5)):
        words=[w for w in all_words if sum(j>=2 for j in w)<=cut]
        rays=resistance_table(NORMAL,words)
        assert rays.rank()==expected
        print(len(words),'Gale monomials, resistance rank',expected)
    for paths in ((0,0,1,7,8),(0,0,0,7,8),(0,0,1,7,7),(0,0,0,7,7)):
        channels=symbolic_support(paths)
        second=['xxc','xxd','xyc','xyd','yyc','yyd'] if paths[1]!=paths[2] else []
        assert channels==[['xxx','xxy','xyy','yyy'],second,[],[]]
        print(paths,'four polynomial supports',channels)
