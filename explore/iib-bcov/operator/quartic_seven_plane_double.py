"""Exact seven-planar quadratic lift on the four-loop two-double core.

Each normal is at a labeled bivalent vertex. All fourteen segment times and
all seven ordered planar insertions are integrated. Values strip the planar
Cauchy–Binet tensor, which is 144 t_2^2 on 3e1,3e2,t with an e1 lift.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache, reduce
from itertools import combinations_with_replacement, permutations, product
from math import factorial, gcd, lcm, prod
from pathlib import Path
import sys

from flint import fmpq_mat, fmpz_mpoly_ctx

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from bulk_three import forest_data
from square_three import G, Graph, plan_generic, sp_monomial
from current_minor import numerator, polynomial_det
from quartic_seven_plane_lift import NORMAL, fixture, density, solve_graph
from quartic_seven_plane_bare import cofactor_square

CORE=((0,1),(0,1),(0,2),(1,3),(2,3),(2,4),(3,5),(4,5),(4,5))
CYCLES=((1,-1,0,0,0,0,0,0,0),(-1,0,1,-1,1,0,0,0,0),
        (0,0,0,0,-1,1,-1,1,0),(0,0,0,0,0,0,0,-1,1))
Q=13
SCALE=factorial(23)**3*factorial(13)**3


def graph_data(marks):
    edges=fixture(CORE,marks)
    paths=[]
    for i in range(9):
        start=sum(map(len,paths))
        paths.append(tuple(range(start,start+marks.count(i)+1)))
    graph=G(11,edges)
    graph.C=[[row[i] for i,path in enumerate(paths) for _ in path] for row in CYCLES]
    ctx=fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(14)))
    times=ctx.gens()
    gram=[[sum(x*y*t for x,y,t in zip(row,col,times)) for col in graph.C] for row in graph.C]
    adj=[[(-1)**(i+j)*polynomial_det([row[:i]+row[i+1:] for k,row in enumerate(gram) if k!=j],ctx)
          for j in range(4)] for i in range(4)]
    _,u,dist=forest_data(edges,11)
    u=ctx.from_dict({k:int(v) for k,v in u.terms()})
    assert u==polynomial_det(gram,ctx)
    groups=(sum(paths[:4],()),paths[4],sum([paths[i] for i in (7,8,5,6)],()))
    splits=((len(paths[0]),len(paths[1])),(len(paths[7]),len(paths[8])))
    return graph,ctx,u,adj,dist,groups,splits


def tree_flow(graph,sources):
    _,parent=graph._tree()
    flows=[[0]*4 for _ in graph.edges]
    for v,source in enumerate(sources):
        while parent[v] is not None:
            w,e=parent[v]
            sign=1 if graph.edges[e]==(v,w) else -1
            for i,x in enumerate(source):
                flows[e][i]+=sign*x
            v=w
    return flows


def polynomials(data,order,original=False,*,normal=NORMAL,slopes=(1,-1,0,0,0),pair=(0,1),polarization=(0,0,1)):
    """Current and voltage polynomials at integral momenta, lifts and polarization."""
    graph,ctx,u,adj,dist,_,_=data
    if any(map(sum,zip(*normal))) or sum(slopes):
        raise ValueError('Normal sources and tangential variations must conserve momentum.')
    src=[(0,0,0,0)]*6+[(slopes[i],*normal[i]) for i in order]
    plus,minus=(6+order.index(i) for i in pair)
    charges=[sum(x*y for x,y in zip(row,polarization)) for row in normal]
    v=-sum(charges[k]*(dist[plus][6+j]-dist[minus][6+j]) for j,k in enumerate(order))/2
    assert all(c.denominator==1 for _,c in v.terms())
    voltage=ctx.from_dict({m:int(c) for m,c in v.terms()})
    flow=tree_flow(graph,src)
    ds=[(-1)**r*numerator(graph,src,flow,[(0,j) for j in range(4) if j!=r]) for r in range(4)]
    times=[2+(i*i+3*i)%17 for i in range(14)]
    value=lambda poly: int(poly(*times))
    y,phi=solve_graph(graph.edges,src,times,10)
    potential=sum(z*(a-b) for z,a,b in zip(polarization,phi[plus][1:],phi[minus][1:]))
    assert value(voltage)==value(u)*potential
    for r,d in enumerate(ds):
        rows=[[graph.C[s][e]*y[e][i] for e in range(1,14)] for i in (1,2,3) for s in range(4)]
        rows.append([graph.C[r][e]*y[e][0] for e in range(1,14)])
        direct=F(str(fmpq_mat([[str(x) for x in row] for row in rows]).det()))*value(u)**4/times[0]
        assert value(d)==direct
    signal=sum(value(ds[r])*value(adj[r][s])*value(ds[s]) for r in range(4) for s in range(4))
    if original:
        plane=(2,3)
        lifts=[(0,0)]*10+[(-3-plane[0],-3-plane[1])]
        values=[]
        for lam in (-1,0,1):
            tangential=[(x+lam*q[0],y) for (x,y),q in zip(lifts,src)]
            values.append(cofactor_square(graph.edges,graph.C,times,[q[1:] for q in src],plane,10,
                                         tangential=tangential))
        assert (values[0]+values[2]-2*values[1])*value(u)**7==2*144*plane[1]**2*signal
        # The unintegrated angular weight has scaling degree +2. The physical
        # projective density includes four further inverse Kirchhoff factors.
        default=(tuple(normal)==NORMAL and tuple(slopes)==(1,-1,0,0,0) and tuple(pair)==(0,1)
                 and tuple(polarization)==(0,0,1))
        if default:
            assert density(graph.edges,times,orders=[order])[1]/value(u)**4==F(signal*value(voltage),value(u)**13)
        if default and order==tuple(range(5)):
            sheared=[(2*q[0]+q[1]-q[3],-q[0]+q[2]+q[3]) for q in src]
            weights=(-4,-8,-14,7,0)
            wx,wy=(F(sum(a*b[j] for a,b in zip(weights,sheared[6:])),4) for j in range(2))
            assert (wx,wy)==(2,-1)
            mixed=[]
            for lam in (-1,0,1):
                tangential=[(a+lam*x,b+lam*y) for (a,b),(x,y) in zip(lifts,sheared)]
                mixed.append(cofactor_square(graph.edges,graph.C,times,[q[1:] for q in src],plane,10,
                                             tangential=tangential))
            x,y=plane
            planar=36*(4*(y*y*wx*wx+x*x*wy*wy)-2*x*y*wx*wy)
            assert (mixed[0]+mixed[2]-2*mixed[1])*value(u)**7==2*planar*signal
    return ds,voltage


def primitive(terms):
    denominator=lcm(*(v.denominator for _,v in terms))
    integer=[(m,int(v*denominator)) for m,v in terms]
    content=reduce(gcd,(v for _,v in integer))
    if integer[0][1]<0:
        content=-content
    return tuple((m,v//content) for m,v in integer),F(content,denominator)


def decompose(poly,groups):
    # Slice by the middle exponents and by the two local degrees. The remaining
    # coefficient matrix is exactly its pivot columns times its nonzero RREF rows.
    slices=defaultdict(dict)
    for mon,c in poly.terms():
        left,middle,right=(tuple(int(mon[i]) for i in group) for group in groups)
        slices[sum(left),middle,sum(right)][left,right]=int(c)
    for (_,middle,_),entries in slices.items():
        left=sorted({a for a,b in entries})
        right=sorted({b for a,b in entries})
        mat=fmpq_mat([[entries.get((a,b),0) for b in right] for a in left])
        rows,rank=mat.rref()
        for i in range(rank):
            pivot=next(j for j in range(len(right)) if rows[i,j])
            x,c=primitive([(m,F(str(mat[k,pivot]))) for k,m in enumerate(left) if mat[k,pivot]])
            y,d=primitive([(m,F(str(rows[i,k]))) for k,m in enumerate(right) if rows[i,k]])
            yield x,y,middle,c*d


@cache
def local_moment(split,exponents):
    # Two inner paths in parallel, in series with every remaining segment.
    k,l=split
    a=sum(exponents[:k])+k
    b=sum(exponents[k:k+l])+l
    rest=exponents[k+l:]
    power=sum(rest)+len(rest)+a+b-Q
    if min(Q-a,Q-b,a+b-Q)<=0:
        raise ArithmeticError(('nonintegrable local factor',split,exponents))
    num=prod(map(factorial,exponents))*prod(map(factorial,(Q-1-a,Q-1-b,a+b-Q-1)))
    den=prod(map(factorial,(a-1,b-1,2*Q-1-a-b,power-1)))
    v,remainder=divmod(SCALE*num,den)
    assert not remainder
    return v


class Arm:
    def __init__(self,group,split):
        self.group=group
        self.split=split
        self.ctx=fmpz_mpoly_ctx.get(tuple(f'x{i}' for i in range(len(group))))
        self.ids={}
        self.polys=[]
        self.degrees=[]
        self.products={}
        self.moments={}

    def register(self,terms):
        if terms not in self.ids:
            self.ids[terms]=len(self.polys)
            self.polys.append(self.ctx.from_dict(dict(terms)))
            self.degrees.append(sum(terms[0][0]))
        return self.ids[terms]

    def integrate(self,i,j,k):
        if (i,j,k) not in self.moments:
            if (i,j) not in self.products:
                self.products[i,j]=self.polys[i]*self.polys[j]
            poly=self.products[i,j]*self.polys[k]
            self.moments[i,j,k]=sum(int(c)*local_moment(self.split,tuple(map(int,m))) for m,c in poly.terms())
        return self.moments[i,j,k]


def period(data,order,original=False,**kinematics):
    graph,_,_,adj,_,groups,splits=data
    ds,voltage=polynomials(data,order,original,**kinematics)
    left,right=(Arm(group,split) for group,split in zip((groups[0],groups[2]),splits))
    def parts(poly):
        return [(left.register(l),right.register(r),m,a) for l,r,m,a in decompose(poly,groups)]
    dparts=list(map(parts,ds))
    total=F()
    tested=False
    for i,j in combinations_with_replacement(range(4),2):
        if not ds[i] or not ds[j] or not adj[i][j]:
            continue
        terms=defaultdict(F)
        for x,y,z in product(dparts[i],dparts[j],parts(voltage*adj[i][j])):
            m=tuple(sum(col) for col in zip(x[2],y[2],z[2]))
            key=(*sorted([x[0],y[0]]),z[0],*sorted([x[1],y[1]]),z[1],m)
            terms[key]+=x[3]*y[3]*z[3]
        for (a,b,c,d,e,f,m),coef in terms.items():
            if not coef:
                continue
            sL=sum(left.degrees[h] for h in (a,b,c))+len(left.group)-Q
            sR=sum(right.degrees[h] for h in (d,e,f))+len(right.group)-Q
            sM=sum(m)+len(m)
            assert sL+sR+sM==2*Q and max(sL,sR,sM)<Q
            outer=F(prod(factorial(Q-1-s) for s in (sL,sM,sR)),factorial(Q-1))
            middle=F(prod(map(factorial,m)),factorial(sM-1))
            if not tested:
                ex=[0]*14
                factored=middle*outer
                for arm,ids in ((left,(a,b,c)),(right,(d,e,f))):
                    exp=tuple(sum(col) for col in zip(*(next(iter(arm.polys[h].monoms())) for h in ids)))
                    for index,power in zip(arm.group,exp):
                        ex[index]=power
                    factored*=F(local_moment(arm.split,exp),SCALE)
                for index,power in zip(groups[1],m):
                    ex[index]=power
                assert factored==sp_monomial(ex,Q,plan_generic(Graph(11,tuple(graph.edges))))
                tested=True
            total+=(1+(i!=j))*coef*left.integrate(a,b,c)*right.integrate(d,e,f)*middle*outer/SCALE**2
    return total


@cache
def symmetries():
    # A path image carries orientation, so interior vertices on reversed paths reverse order.
    actions=[]
    for flip,exchange,pair1,pair2 in product(range(2),repeat=4):
        vertex=list(range(6))
        if flip:
            vertex=[1,0,3,2,5,4]
        if exchange:
            vertex=[[4,5,2,3,0,1][v] for v in vertex]
        action=[]
        for i,(a,b) in enumerate(CORE):
            image=(vertex[a],vertex[b])
            target=tuple(sorted(image))
            choices=[j for j,e in enumerate(CORE) if e==target]
            j=choices[0]
            if len(choices)==2:
                j=choices[(i in (1,8))^(pair1 if j==0 else pair2)]
            action.append((j,image!=target))
        actions.append(tuple(action))
    assert len(set(actions))==16
    return tuple(actions)


def image(configuration,action):
    result=[() for _ in range(9)]
    for values,(j,reverse) in zip(configuration,action):
        result[j]=values[::-1] if reverse else values
    return tuple(result)


def shapes():
    candidates={tuple(sorted((*l,m,*r))) for l in combinations_with_replacement((0,1),2)
                for r in combinations_with_replacement((7,8),2) for m in range(2,7)}
    while candidates:
        marks=min(candidates)
        orbit={tuple(sorted(action[h][0] for h in marks)) for action in symmetries()}
        candidates-=orbit
        yield marks,len(orbit)


def orders(marks):
    stabilizer=[a for a in symmetries() if sorted(a[h][0] for h in marks)==list(marks)]
    unseen=set(permutations(range(5)))
    while unseen:
        order=min(unseen)
        config=tuple(tuple(order[i] for i,h in enumerate(marks) if h==j) for j in range(9))
        orbit={sum(image(config,action),()) for action in stabilizer}
        unseen-=orbit
        yield order,len(orbit)


if __name__=='__main__':
    shapes_all=list(shapes())
    assert sum(size for _,size in shapes_all)==45
    if '--list' in sys.argv:
        for index,(marks,size) in enumerate(shapes_all):
            print(index,marks,size,len(list(orders(marks))))
    else:
        by_shape=(F(63498894979,3300620400000000), F(14713395553,825155100000000),
                  F(623468081,141455160000000), F(26791687,3857868000000),
                  F(7418973031,275051700000000), F(4310186269,168781725000000),
                  F(53129669,4973032968750))
        assert sum(s*w/16 for s,(_,w) in zip(by_shape,shapes_all))==F(367378814983,7129340064000000)
        selected=range(7) if '--all' in sys.argv else ([int(x) for x in sys.argv[1:]] or [6])
        for i in selected:
            marks,size=shapes_all[i]
            data=graph_data(marks)
            total=0
            for order,weight in orders(marks):
                v=period(data,order,original=order==tuple(range(5)))
                if marks==(0,1,4,7,8) and order==tuple(range(5)):
                    assert v==F(968497,4059618750000)
                total+=weight*v
                print('order',i,order,weight,v,flush=True)
            assert total==by_shape[i]
            local_moment.cache_clear()
            print('shape',i,marks,size,total,total*size/16,flush=True)
