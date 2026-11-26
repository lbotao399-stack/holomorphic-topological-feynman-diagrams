"""Eight Betti factors on the genus-four three-mark planar simplex."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations,permutations,product
from math import comb,factorial,prod

from betti_flow_polytope import weak
from betti_power_graphs import vector_coefficient
from betti_series_parallel import coefficient
from soft_quartic_psi_bound import intersection


def labelings(edges):
    n=max(max(e) for e in edges)+1
    adjacency=[[] for _ in range(n)]
    for a,b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    colors=tuple(map(len,adjacency))
    while True:
        signatures=[(c,tuple(sorted(colors[j] for j in adjacent))) for c,adjacent in zip(colors,adjacency)]
        palette={key:i for i,key in enumerate(sorted(set(signatures)))}
        refined=tuple(palette[key] for key in signatures)
        if len(set(refined))==len(set(colors)): break
        colors=refined
    blocks=[tuple(i for i,c in enumerate(refined) if c==value) for value in sorted(set(refined))]
    for pieces in product(*(permutations(block) for block in blocks)):
        perm=[0]*n
        for j,v in enumerate(v for group in pieces for v in group): perm[v]=j
        yield tuple(perm)


def image(edges,p):
    return tuple(sorted((min(p[i],p[j]),max(p[i],p[j])) for i,j in edges))


@cache
def graphs(size):
    if size==0: return ((),)
    if size==1: return (((0,1),),)
    raw=set()
    for edges in graphs(size-1):
        n=max(max(e) for e in edges)+1
        for e in combinations(range(n+1),2):
            raw.add(tuple(sorted((*edges,e))))
    return tuple(sorted({min(image(edges,p) for p in labelings(edges)) for edges in raw}))


@cache
def marked_shapes(edges):
    if not edges: return ((1,(0,0,0),1),)
    n=max(max(e) for e in edges)+1
    autom=[p for p in labelings(edges) if image(edges,p)==edges]
    assert autom
    pending=set(product(range(n),repeat=3))
    result=[]
    parallel=prod(factorial(m) for m in Counter(edges).values())
    while pending:
        marks=min(pending)
        orbit={tuple(p[i] for i in marks) for p in autom}
        pending-=orbit
        result.append((n,marks,len(autom)*parallel//len(orbit)))
    return tuple(result)


LEVELS=((4,0,0),(3,1,0),(2,2,0),(2,1,1))
LEVELS_BY_GENUS={3:((3,0,0),(2,1,0),(1,1,1)),4:LEVELS}
EXPECTED={
    3:(Fraction(7,10368000), Fraction(17,7257600), Fraction(11,2177280)),
    4:(Fraction(731,1024192512000), Fraction(1493,438939648000),
       Fraction(157,29262643200), Fraction(433,37623398400)),
}


def stencil(genus,size,level_sets,*,shape_subset=None):
    """Impose component dimension before requesting a graph moment."""
    answer=Counter()
    loads=((1,0),(0,1),(-1,-1))
    for edges in graphs(size) if shape_subset is None else shape_subset:
        for n,marks,aut in marked_shapes(edges):
            remainder=genus-(size-n+1)
            if remainder<0: continue
            ends=[[i for i,e in enumerate(edges) if j in e] for j in range(n)]
            sites=[[i for i,m in enumerate(marks) if m==j] for j in range(n)]
            slots=[len(ss)+len(ee) for ss,ee in zip(sites,ends)]
            charges=tuple(tuple(sum(loads[i][k] for i in ss) for k in (0,1)) for ss in sites)
            if size and not any(map(any,charges)): continue
            # A branch carrying zero total current vanishes in its full-support moment.
            if any(len(e)==1 and not any(q) for e,q in zip(ends,charges)): continue
            for h in weak(remainder,n):
                if any(2*t-2+b<=0 for t,b in zip(h,slots)): continue
                for channel,levels in enumerate(level_sets):
                    missing=[3*t-3+f-sum(levels[i] for i in ss) for t,f,ss in zip(h,slots,sites)]
                    if min(missing)<0: continue
                    choices=[]
                    for t,ss,ee,amount in zip(h,sites,ends,missing):
                        cell=[]
                        for powers in weak(amount,len(ss)+len(ee)):
                            leg=powers[:len(ss)]
                            branch=powers[len(ss):]
                            value=intersection(t,(*[levels[i]+a for i,a in zip(ss,leg)],*branch))
                            for i,a in zip(ss,leg):
                                value*=Fraction(sum(x*x for x in loads[i])**a,2**a*factorial(a))
                            if value: cell.append((dict(zip(ee,branch)),value))
                        choices.append(cell)
                    for selection in product(*choices):
                        weight=prod(w for dd,w in selection)/aut
                        normal=[]
                        for i,(a,b) in enumerate(edges):
                            x,y=selection[a][0][i],selection[b][0][i]
                            normal.append(x+y)
                            weight*=comb(x+y,x)
                        answer[edges,charges,tuple(normal),channel]+=weight
    return answer


def sp(edges):
    neighbors={}
    for u,v in edges:
        neighbors.setdefault(u,set()).add(v)
        neighbors.setdefault(v,set()).add(u)
    while neighbors:
        leaf=next((u for u,V in neighbors.items() if len(V)<=2),None)
        if leaf is None: return False
        V=neighbors.pop(leaf)
        for u in V: neighbors[u].remove(leaf)
        if len(V)==2:
            u,v=V
            neighbors[u].add(v)
            neighbors[v].add(u)
    return True


@cache
def projection_edges(edges):
    return tuple(i for i in range(len(edges)) if sp(edges[:i]+edges[i+1:]))


def moment(edges,charges,powers):
    if not edges: return Fraction(1)
    if sp(edges): return vector_coefficient(edges,charges,powers)
    roots=projection_edges(edges)
    return coefficient(edges,charges,powers,roots[0] if roots else 0)


def calculate(genus,size,check_roots=False):
    levels=LEVELS_BY_GENUS[genus]
    weights=stencil(genus,size,levels)
    grouped={}
    for (edges,charges,powers,i),w in weights.items():
        grouped.setdefault((edges,charges,powers),[Fraction() for _ in levels])[i]+=w
    totals=[Fraction() for _ in levels]
    checked=0
    for (edges,q,powers),row in grouped.items():
        value=moment(edges,q,powers)
        if check_roots and not sp(edges):
            roots=projection_edges(edges)
            if len(roots)<2:
                raise ValueError(f'Two distinct projection edges are required for {edges}')
            assert value==coefficient(edges,q,powers,roots[-1]),(edges,q,powers,roots)
            checked+=1
        for i,w in enumerate(row): totals[i]+=w*value
    print(genus,size,len(graphs(size)),len(grouped),*(str(t) for t in totals),flush=True)
    if check_roots:
        print('Non-series-parallel moments verified with two projections:',checked,flush=True)
    return totals


if __name__=='__main__':
    from argparse import ArgumentParser
    p=ArgumentParser(description=__doc__)
    p.add_argument('--genus',choices=(3,4),type=int,default=4)
    p.add_argument('--edges',type=int,nargs='*')
    p.add_argument('--check-roots',action='store_true',
                   help='Reintegrate every non-series-parallel input from a second projected edge.')
    args=p.parse_args()
    complete=list(range(2*args.genus+1))
    sizes=complete if args.edges is None else args.edges
    if not sizes or any(e not in complete for e in sizes) or len(set(sizes))!=len(sizes):
        p.error('Select distinct edge counts between zero and twice the genus.')
    rows=[calculate(args.genus,e,args.check_roots) for e in sizes]
    total=tuple(map(sum,zip(*rows)))
    print('Sum of selected edge counts:',*total,flush=True)
    if sorted(sizes)==complete:
        assert total==EXPECTED[args.genus],(args.genus,total)
