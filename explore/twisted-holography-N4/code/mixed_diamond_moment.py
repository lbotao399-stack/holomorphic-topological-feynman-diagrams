"""Full-support Newton moment of the three-path diamond with Euclidean vector charges."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial,prod
from flint import fmpq, fmpq_mpoly_ctx


@lru_cache(maxsize=4096)
def cell_moment(powers):
    a,b,*rest=powers
    return fmpq(factorial(a)*factorial(b),factorial(a+b+2)*prod(i+1 for i in rest))


def determinant(rows):
    if len(rows)==1: return rows[0][0]
    return sum((-1)**i*entry*determinant([r[:i]+r[i+1:] for r in rows[1:]])
               for i,entry in enumerate(rows[0]) if entry)


@lru_cache(maxsize=4096)
def diamond(gram,powers=(0,0,0,0,0)):
    """Gram of source, first and second intermediate charges. Path order 0,1L,1R,2L,2R."""
    total=fmpq(0)
    for first,second in product(('L','R','LR'),repeat=2):
        states=(first,second)
        names=('a','b',*('s'+str(i) for i,x in enumerate(states) if x=='LR'))
        ctx=fmpq_mpoly_ctx.get(names)
        variables=ctx.gens()
        s=[]
        for i,state in enumerate(states):
            s.append(variables[names.index('s'+str(i))] if state=='LR' else ctx.constant(int(state=='R')))
        direction=(ctx.constant(1),*s)
        shares=(1-variables[0]-variables[1],variables[0],variables[1])
        current=[tuple(shares[0]*z for z in direction)]
        for i in range(2):
            left=tuple(shares[i+1]*z-(s[i] if j==i+1 else 0) for j,z in enumerate(direction))
            current.extend((left,tuple(z+int(j==i+1) for j,z in enumerate(left))))
        def energy(row):
            return sum((fmpq(str(gram[i][j]))*(1 if i==j else 2)*row[i]*row[j]
                        for i in range(3) for j in range(i,3)),ctx.constant(0))/2
        heights=tuple(map(energy,current))
        active=[1+2*i+j for i,x in enumerate(states) for j,c in enumerate('LR') if c in x]
        # These interlaced variables have positive orientation relative to edge lengths/root length.
        coords=[j for i,x in enumerate(states) for j in
                ((i,names.index('s'+str(i))) if x=='LR' else (i,))]
        jac=determinant([[heights[e].derivative(j) for j in coords] for e in active])
        if not jac: continue
        density=jac
        for e,(height,k) in enumerate(zip(heights,powers)):
            if e in active: density*=height**k
            else: density*=height**(k+1)/(k+1)
        # The divergence, inactive intervals and the concave height Hessian give net sign (-1)^5.
        total-=sum((c*cell_moment(k) for k,c in density.terms()),fmpq(0))
    return F(str(total))/prod(map(factorial,powers))


def graph_value(edges,charges,powers):
    degree=tuple(sum(i in e for e in edges) for i in range(len(charges)))
    if sorted(degree)!=[2,2,3,3] or len(set(edges))!=5:
        raise ValueError('Expected the five distinct edges of a diamond.')
    source,sink=[i for i,d in enumerate(degree) if d==3]
    u,v=[i for i,d in enumerate(degree) if d==2]
    order=((source,sink),(source,u),(u,sink),(source,v),(v,sink))
    by_edge={tuple(sorted(e)):k for e,k in zip(edges,powers)}
    gram=tuple(tuple(sum(a*b for a,b in zip(charges[i],charges[j])) for j in (source,u,v))
               for i in (source,u,v))
    return diamond(gram,tuple(by_edge[tuple(sorted(e))] for e in order))


def check():
    from scalar_flow_recursion import coefficient
    from betti_first_variation import normal_derivative
    edges=((0,3),(0,1),(1,3),(0,2),(2,3))
    for p in ((1,-1,-1,1),(2,-3,1,0),(1,2,3,-6)):
        for k in ((0,)*5,(0,2,0,0,1)):
            assert graph_value(edges,tuple((x,) for x in p),k)==coefficient(edges,p,k)
    p=(1,-1,0,0)
    q=(0,1,0,-1)
    for k in ((0,)*5,(1,0,0,0,0)):
        d=5+sum(k)
        weights=[sum(F(-1,j*j) for j in range(1,d+1))]
        weights.extend(F(1,j*j)*prod(F(-s*s,j*j-s*s) for s in range(1,d+1) if s!=j)
                       for j in range(1,d+1))
        values=[graph_value(edges,tuple(zip(p,(j*x for x in q))),k) for j in range(d+1)]
        assert sum(w*x for w,x in zip(weights,values))==normal_derivative(edges,p,q,k)
    charges=((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))
    perm=(2,3,0,1)
    expected=graph_value(edges,charges,(0,)*5)
    assert expected==graph_value(tuple((perm[a],perm[b]) for a,b in edges),
                                 tuple(charges[perm.index(i)] for i in range(4)),(0,)*5)
    print('Unit three-axis diamond:',expected,flush=True)
    print('Independent scalar graph and polarized facet controls passed.',flush=True)


if __name__=='__main__':
    check()
