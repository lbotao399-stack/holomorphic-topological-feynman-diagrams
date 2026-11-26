"""Four-loop kernels for the three series-parallel cubic graphs and their cyclic trace."""
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import json
import sys

import sympy as sp

here=Path(__file__).resolve().parent
root=here.parents[3]
sys.path.insert(0,str(root/'explore/ht-families/verify'))
import vcore
sys.path.insert(0,str(root/'deepen/o3/code'))
import ribbon

graphs=json.loads((here.parent/'four-loop-graphs.json').read_text())[:3]
probe=[(1,0),(0,1),(1,1),(-2,-2)]


def cross(p,q):
    return p[0]*q[1]-p[1]*q[0]


def chart(chains,coordinates):
    a,b,u,v,x,y,s,t,ma,na,mb,nb=coordinates
    c=1-a-b
    times={}
    for path,sides,w,share,frac,right,m,n in [
        ([0,1,2,3],[(p,False) for p in chains if p[0]==1 and p[-1]==2],a,u,x,s,ma,na),
        ([0,5,4,3],[(p,True) for p in chains if p[0]==4 and p[-1]==5],b,v,y,t,mb,nb)]:
        times[frozenset(path[:2])]=(1-frac-right)/w
        times[frozenset(path[2:])]=right/w
        splits=[[1],[1-m,m-n,n]] if len(sides[0][0])==2 else [[1-m,m],[1-n,n]]
        for (vertices,reverse),fractions,part in zip(sides,splits,[share,1-share]):
            vs=vertices[::-1] if reverse else vertices
            for i,r in enumerate(fractions):
                times[frozenset(vs[i:i+2])]=r*frac/(w*part)
    times[frozenset((0,3))]=1/c
    edges=[tuple(pair) for chain in chains for pair in zip(chain,chain[1:])]
    return edges,[times[frozenset(e)] for e in edges]


def kernel_density_squared(shapes,at,p):
    a,b,u,v,x,y,s,t,ma,na,mb,nb=at
    c=1-a-b
    ua=u*u if shapes[0]==(1,3) else -u*(1-u)
    ub=v*v if shapes[1]==(1,3) else -v*(1-v)
    amp=F(3,4)*ua*ub*(ma-na)*(mb-nb)*cross(*p[:2])**2*cross(*p[2:])**2
    return amp, amp**2*c**3*x*y/(u*(1-u)*v*(1-v))


def check_shape(index,fixtures):
    chains=graphs[index]['chains']
    shapes=[tuple(len(ch)-1 for ch in chains if ch[0]==i and ch[-1]==j) for i,j in [(1,2),(4,5)]]
    coordinates=sp.symbols('a b u v x y s t m n M N')
    edges,res=chart(chains,coordinates)
    matrix=sp.Matrix(res)
    jac=matrix.row_join(matrix.jacobian(coordinates))
    net=vcore.Net(10,edges)
    cycles=net.cycle_basis()
    sign=None
    for at,p in fixtures:
        vals=dict(zip(coordinates,at))
        times=[F(r.subs(vals)) for r in res]
        volume=abs(F(jac.subs(vals).det(method='domain-ge')))
        lam=[(0,0)]*6+list(p[:3])
        coef,U,_=vcore.density_slice(net,cycles,times,lam)
        amp,square=kernel_density_squared(shapes,at,p)
        assert (coef*volume)**2==U*square, ((coef*volume)**2/(U*square) if square else (coef,square))
        if amp:
            orientation=1 if amp*coef>0 else -1
            if sign is None:
                sign=orientation
            assert sign==orientation
    return sign


def cyclic(indices,kernels):
    result={}
    for index in indices:
        edges=[tuple(pair) for chain in graphs[index]['chains'] for pair in zip(chain,chain[1:])]
        tensor=ribbon.colour_weight(10,edges,list(range(6,10)))
        by_rank={}
        for s in permutations(range(4)):
            p=[probe[i] for i in s]
            coefficient=ribbon.perm_sign(s)*kernels[index]*cross(*p[:2])**2*cross(*p[2:])**2
            for cycles,ranks in tensor.items():
                ordered=ribbon.canon_tau([tuple(s[i] for i in cycle) for cycle in cycles])
                if ordered==((0,1,2,3),):
                    for power,weight in ranks.items():
                        by_rank[power]=by_rank.get(power,F(0))+coefficient*weight
        result[index]={degree:val/graphs[index]['aut'] for degree,val in by_rank.items() if val}
    return result


def exterior_identity():
    a,b,c,x,y,K,J=sp.symbols('a b c x y K J')
    total=x*K+y*J
    first,second=a*total-x*K,b*total-y*J
    variables=(a,b,x,y)
    euler={(2,3):c,(0,2):y*b,(1,2):y*(b+c),(0,3):-x*(a+c),
           (1,3):-x*a,(0,1):3*x*y}
    top=0
    for (i,j),value in euler.items():
        k,l=[s for s in range(4) if s not in (i,j)]
        minor=first.diff(variables[k])*second.diff(variables[l])-first.diff(variables[l])*second.diff(variables[k])
        top+=ribbon.perm_sign((i,j,k,l))*value*minor
    assert sp.expand((top-3*c*x*y*K*J).subs(c,1-a-b))==0


def matrix_unit_sum(index,kernel):
    # Independently glue fundamental row/column indices, pinning the four matrix units.
    chains=graphs[index]['chains']
    edges=[pair for ch in chains for pair in zip(ch,ch[1:])]
    half=[[] for _ in range(10)]
    for v in range(6,10):
        half[v].append(v-6)
    for e,(v,w) in enumerate(edges):
        half[v].append(4+2*e)
        half[w].append(5+2*e)
    sign=ribbon.perm_sign([h for slots in half for h in slots])
    answer={}
    for chosen in product((0,1),repeat=10):
        parents=list(range(60))
        def representative(i):
            while i!=parents[i]:
                parents[i]=parents[parents[i]]
                i=parents[i]
            return i
        def glue(i,j):
            parents[representative(i)]=representative(j)
        for i in range(13):
            h=4+2*i
            glue(2*h,2*(h+1)+1)
            glue(2*h+1,2*(h+1))
        for v,slots in enumerate(half):
            h=slots[::-1] if chosen[v] else slots
            for i,t in enumerate(h):
                glue(2*t+1,2*h[(i+1)%3])
        free=len({representative(h) for h in range(60)}-{representative(h) for h in range(8)})
        vertex_sign=sign*(-1)**sum(chosen)
        for perm in permutations(range(4)):
            labels={}
            for h,colour in enumerate(perm):
                for t,fixed in [(2*h,colour),(2*h+1,(colour+1)%4)]:
                    labels.setdefault(representative(t),set()).add(fixed)
            if any(len(values)>1 for values in labels.values()):
                continue
            p=[probe[i] for i in perm]
            answer[free]=answer.get(free,F(0))+vertex_sign*ribbon.perm_sign(perm)*kernel*cross(*p[:2])**2*cross(*p[2:])**2
    return {power:v/graphs[index]['aut'] for power,v in answer.items() if v}


def main():
    exterior_identity()
    at1=tuple(F(a,b) for a,b in [(1,4),(1,5),(1,3),(2,5),(2,7),(3,8),(1,9),(1,11),(3,4),(1,4),(4,5),(2,5)])
    at2=tuple(F(a,b) for a,b in [(2,5),(1,3),(2,7),(3,5),(1,5),(3,10),(2,9),(1,7),(2,3),(1,3),(5,6),(2,7)])
    fixtures=[(at1,[probe[i] for i in [0,2,1,3]]),(at2,[(1,2),(2,-1),(0,1),(-3,-2)])]
    signs=[check_shape(i,fixtures) for i in range(3)]
    print('oriented pointwise density signs',signs,flush=True)
    assert signs==[-1,1,1]
    radial=F(3,4)*F(4,35)*F(4,15)**2
    kernel=radial*F(3,8)**2*F(1,6)**2
    assert kernel==F(1,42000)
    result=cyclic(range(3),[signs[0]*kernel,0,0])
    assert result=={0:{4:F(1,10500)},1:{},2:{}}
    assert matrix_unit_sum(0,signs[0]*kernel)==result[0]
    print('cyclic contributions / pi^2',result)

if __name__=='__main__':
    main()
