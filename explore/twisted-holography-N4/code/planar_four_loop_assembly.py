"""Global leg/normal expansion and first-KdV check of the genus-four planar stable-graph stencil."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import comb, factorial, prod

import networkx as nx
from planar_four_loop import LEVELS, graphs, stencil


@cache
def cotangent(g, ds):
    """First KdV, string and dilaton; eliminate the largest degree first."""
    n=len(ds)
    if g<0 or 2*g+n<=2 or sum(ds)!=3*g-3+n:
        return F()
    if not g:
        return F(factorial(n-3),prod(map(factorial,ds)))
    def child(h, degrees):
        return cotangent(h,tuple(sorted(degrees,reverse=True)))
    if ds[-1]==0:
        return sum((child(g,(*ds[:i],d-1,*ds[i+1:-1])) for i,d in enumerate(ds[:-1]) if d),F())
    if ds[-1]==1 and n>1:
        return (2*g+n-3)*cotangent(g,ds[:-1])
    # Apply the first KdV relation to tau_1 tau_0 with the first exponent raised.
    raised=(ds[0]+1,*ds[1:])
    rhs=child(g-1,(*raised,0,0,0,0))/12
    for mask in range(1,2**n-1):
        left=tuple(d for i,d in enumerate(raised) if mask>>i&1)
        right=tuple(d for i,d in enumerate(raised) if not (mask>>i&1))
        genus,mod=divmod(sum(left)+1-len(left),3)
        if not mod:
            rhs+=child(genus,(*left,0,0))*child(g-genus,(*right,0,0))/2
    reduced=rhs/(2*g+n-1)
    return reduced-sum((child(g,(raised[0],*ds[1:i],d-1,*ds[i+1:]))
                        for i,d in enumerate(ds[1:],1)),F())


def compose(total, slots):
    if not slots:
        if not total: yield ()
        return
    for bars in combinations(range(total+slots-1),slots-1):
        cuts=(-1,*bars,total+slots-1)
        yield tuple(b-a-1 for a,b in zip(cuts,cuts[1:]))


@cache
def expansions(size,degree=8):
    result=[]
    for line_degree in range(degree+1-size):
        for legs in compose(line_degree,3):
            leg_weight=F(2**legs[2],2**line_degree*prod(map(factorial,legs)))
            for normal in compose(degree-size-line_degree,size):
                for left in product(*(range(k+1) for k in normal)):
                    right=tuple(k-a for k,a in zip(normal,left))
                    factor=leg_weight*prod(comb(k,a) for k,a in zip(normal,left))
                    result.append((legs,normal,left,right,factor))
    return result


def independent_marks(edges,n):
    graph=nx.Graph()
    graph.add_nodes_from(range(n))
    multiplicity=Counter(edges)
    for (a,b),count in multiplicity.items():
        graph.add_edge(a,b,weight=count)
    matcher=nx.isomorphism.GraphMatcher(graph,graph,edge_match=lambda x,y:x['weight']==y['weight'])
    symmetry=list(matcher.isomorphisms_iter())
    parallel=prod(factorial(w) for w in multiplicity.values())
    pending=set(product(range(n),repeat=3))
    while pending:
        placement=min(pending)
        images={tuple(f[v] for v in placement) for f in symmetry}
        pending.difference_update(images)
        yield placement,F(len(images),parallel*len(symmetry))


@cache
def local_degree(values):
    g,mod=divmod(sum(values)+3-len(values),3)
    if mod or g<0: return 0
    return cotangent(g,tuple(sorted(values,reverse=True)))


def independent(size,genus=4,level_sets=LEVELS,*,shape_subset=None):
    result=Counter()
    for edges in graphs(size) if shape_subset is None else shape_subset:
        n=1+max((v for edge in edges for v in edge),default=0)
        if size-n+1>genus: continue
        flags=[[] for _ in range(n)]
        for i,(a,b) in enumerate(edges):
            flags[a].append((0,i)); flags[b].append((1,i))
        for marks,scale in independent_marks(edges,n):
            legs=tuple(tuple(i for i,v in enumerate(marks) if v==h) for h in range(n))
            charge=tuple((int(0 in ids)-int(2 in ids),int(1 in ids)-int(2 in ids)) for ids in legs)
            if size and (not any(map(any,charge)) or any(len(f)==1 and not any(q) for f,q in zip(flags,charge))):
                continue
            # Before distributing the eight edge/leg powers, discard any channel whose
            # minimum nonnegative local power budgets already exceed the global budget.
            allowed=[]
            for c,levels in enumerate(level_sets):
                needed=0
                for incident,ids in zip(flags,legs):
                    dim=len(incident)+len(ids); current=sum(levels[i] for i in ids)
                    least=max(0,(current-dim+5)//3,(4-dim)//2)
                    needed+=3*least-3+dim-current
                if needed<=2*genus-size: allowed.append((c,levels))
            if not allowed:continue
            for extra,normal,left,right,factor in expansions(size,2*genus):
                halves=(left,right)
                branch=tuple(tuple(halves[s][j] for s,j in incident) for incident in flags)
                for channel,levels in allowed:
                    value=1
                    for h,ids in enumerate(legs):
                        ds=(*branch[h],*(levels[i]+extra[i] for i in ids))
                        value*=local_degree(ds)
                        if not value:break
                    if value:result[edges,charge,normal,channel]+=scale*factor*value
    return result


def check(sizes):
    for g in range(1,5):
        assert cotangent(g,(3*g-2,))==F(1,24**g*factorial(g))
    assert cotangent(2,(3,2))==F(29,5760)
    for size in sizes:
        actual=independent(size)
        reference=stencil(4,size,LEVELS)
        diff={k:actual[k]-reference[k] for k in actual.keys()|reference.keys() if actual[k]!=reference[k]}
        assert not diff,(size,len(diff),next(iter(diff.items())))
        print(size,len(actual),tuple(sum(k[-1]==i for k in actual) for i in range(4)),flush=True)


def star_coordinates(chart,resistance):
    """Recover all cube coordinates from positive original edge resistances."""
    kind=chart[0]
    if kind=='terminal':
        def descend(node):
            if node[0]=='edge':return resistance[node[1]],()
            left,x=descend(node[2]); right,y=descend(node[3])
            total=left+right
            if node[0]=='parallel':
                equivalent=left*right/total; t=right/total
                assert (equivalent/t,equivalent/(1-t))==(left,right)
                return equivalent,(t,*x,*y)
            t=left/total
            assert (total*t,total*(1-t))==(left,right)
            return total,(t,*x,*y)
        return descend(chart[2])[1]
    r=dict(resistance)
    if kind=='parallel':
        _,i,j,_,fresh,child=chart
        a,b=r.pop(i),r.pop(j)
        r[fresh]=a*b/(a+b)
        t=b/(a+b)
        assert (r[fresh]/t,r[fresh]/(1-t))==(a,b)
        return (t,*star_coordinates(child,r))
    if kind=='series':
        _,_,_,_,i,_,j,_,fresh,child=chart
        a,b=r.pop(i),r.pop(j)
        r[fresh]=a+b
        t=b/(a+b)
        assert (r[fresh]*(1-t),r[fresh]*t)==(a,b)
        return (t,*star_coordinates(child,r))
    _,center,star,triangle,child=chart
    leg={(v if u==center else u):r.pop(i) for i,u,v in star}
    reciprocal=sum(1/x for x in leg.values())
    for i,u,v in triangle:r[i]=leg[u]*leg[v]*reciprocal
    # The inverse is homogeneous and maps the whole positive triangle cone back to the star.
    perimeter=sum(r[i] for i,_,_ in triangle)
    for v,x in leg.items():
        assert prod(r[i] for i,a,b in triangle if v in (a,b))/perimeter==x
    return star_coordinates(child,r)


def star_check():
    from betti_star_charts import plan,solve
    from flint import fmpq_mat
    from planar_four_loop import sp
    checked=0
    for size in (6,7,8):
        fixtures=0
        for ends in graphs(size):
            if sp(ends) or not nx.is_biconnected(nx.Graph(ends)):continue
            n=1+max(v for p in ends for v in p)
            if n<4 or size-n+1>4:continue
            edges=tuple((i,a,b) for i,(a,b) in enumerate(ends))
            selections=(*combinations(range(n),2),*combinations(range(n),3))
            for active in selections:
                empty=frozenset(set(range(n))-set(active))
                chart=plan(edges,empty)
                if chart is None:continue
                q={i:[F(0),F(0)] for i in range(n)}
                charges=((1,2),(-1,-2)) if len(active)==2 else ((1,2),(3,-1),(-4,-1))
                for i,row in zip(active,charges):q[i]=list(map(F,row))
                r={i:F(3+i+((i*i+fixtures)%7),2+i%3) for i in range(size)}
                cube=star_coordinates(chart,r)
                assert len(cube)==size-1 and all(0<t<1 for t in cube)
                actual=solve(chart,q,cube,F())
                inc=fmpq_mat([[int(i==a)-int(i==b) for i in range(1,n)] for _,a,b in edges])
                conductance=fmpq_mat([[str(1/r[i]) if i==j else 0 for j in range(size)] for i in range(size)])
                currents=conductance*inc*(inc.transpose()*conductance*inc).solve(fmpq_mat([[str(x) for x in q[i]] for i in range(1,n)]))
                assert all(actual[e]==[F(str(v)) for v in row] for e,row in enumerate(currents.tolist()))
                fixtures+=1
        checked+=fixtures
    print('Star-chart current and inverse-resistance controls:',checked,flush=True)


if __name__=='__main__':
    import sys
    if sys.argv[1:]!=['--stars-only']:
        check(list(map(int,sys.argv[1:])) or range(9))
    star_check()
