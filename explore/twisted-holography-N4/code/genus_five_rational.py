"""Rational genus-five three-mark stable sums with capacity reductions before integration."""
from collections import Counter
from fractions import Fraction as F
from functools import cache, lru_cache
from math import comb, factorial, prod

from genus_five_graph_census import LEVELS, safe_roots
from planar_four_loop import coefficient, graphs, image, labelings, marked_shapes, moment, sp, stencil


COMPLETE = {
    (5,0,0): F(179,442451165184000),
    (4,1,0): F(14951,6083703521280000),
    (3,2,0): F(463,86910050304000),
    (3,1,1): F(5783,506975293440000),
    (2,2,1): F(51787,2897001676800000),
}


def reduced_input(edges,charges,powers):
    packed = sorted(zip(edges,powers))
    return tuple(e for e,k in packed),tuple(charges),tuple(k for e,k in packed)


@cache
def block_plan(edges):
    import networkx as nx
    net=nx.Graph(edges)
    if nx.is_biconnected(net): return ()
    result=[]
    for block in nx.biconnected_components(net):
        verts=sorted(block)
        local=[(e,(verts.index(i),verts.index(j))) for e,(i,j) in enumerate(edges) if i in block and j in block]
        rest=net.copy()
        rest.remove_edges_from(edges[e] for e,uv in local)
        sides=[tuple(sorted(nx.node_connected_component(rest,v))) for v in verts]
        assert sum(map(len,sides))==len(net)
        result.append((local,sides))
    return tuple(result)


@cache
def primitive_labelings(edges):
    """Precompute the edge transports of the canonical simple graph."""
    canonical=min(image(edges,p) for p in labelings(edges))
    transports=[]
    for p in labelings(edges):
        mapped=tuple(tuple(sorted((p[u],p[v]))) for u,v in edges)
        if tuple(sorted(mapped))!=canonical:continue
        transports.append((tuple(p.index(v) for v in range(len(p))),
                           tuple(mapped.index(e) for e in canonical)))
    return canonical,tuple(transports)


def canonical_primitive(edges,charges,powers):
    canonical,transports=primitive_labelings(edges)
    # The scalar product is unchanged under the eight signed coordinate permutations.
    load_orbit={tuple((a*s,b*t) for a,b in (reversed(q) if swap else q for q in charges))
                for s in (-1,1) for t in (-1,1) for swap in (False,True)}
    inputs=min((tuple(load[i] for i in vertices),tuple(powers[j] for j in order))
               for vertices,order in transports for load in load_orbit)
    return canonical,*inputs


@lru_cache(maxsize=32768)
def capacity_moment(edges,charges,powers):
    if not edges: return F(1)
    if not any(map(any,charges)): return F()
    if len(edges)==1:
        return -F(sum(x*x for x in charges[0])**(powers[0]+1),
                  2**(powers[0]+1)*factorial(powers[0]+1))
    parallel=Counter(edges)
    if len(parallel)<len(edges):
        factor=F(1)
        simple=tuple(sorted(parallel))
        merged=[]
        for side in simple:
            ks=[k for e,k in zip(edges,powers) if e==side]
            degree=len(ks)+sum(ks)-1
            merged.append(degree)
            # Ball capacities add. Integrate their positive shares at fixed sum.
            factor*=F(2**(len(ks)-1)*factorial(degree)*prod(factorial(2*k+1) for k in ks),
                      factorial(2*degree+1)*prod(map(factorial,ks)))
        return factor*capacity_moment(simple,charges,tuple(merged))
    for v,q in enumerate(charges):
        inc=[(j,w if u==v else u) for j,(u,w) in enumerate(edges) if v in (u,w)]
        if len(inc)==1 and not any(q): return F()
        if len(inc)!=2 or any(q): continue
        (i,u),(j,w)=inc
        new=[(e,k) for n,(e,k) in enumerate(zip(edges,powers)) if n not in (i,j)]
        degree=powers[i]+powers[j]+1
        new.append((tuple(sorted((u,w))),degree))
        renumber=[tuple(k-(k>v) for k in pair) for pair,power in new]
        smaller=reduced_input(renumber,charges[:v]+charges[v+1:],[k for pair,k in new])
        # The effective squared capacity is the minimum of the squared capacities.
        # Independent exponential variables have additive rates.
        return -comb(degree+1,powers[i]+1)*capacity_moment(*smaller)
    blocks=block_plan(edges)
    if blocks:
        result=F(1)
        for local,sides in blocks:
            loads=tuple(tuple(map(sum,zip(*(charges[v] for v in side)))) for side in sides)
            result*=capacity_moment(tuple(e for j,e in local),loads,tuple(powers[j] for j,e in local))
            if not result: break
        return result
    return primitive_moment(*canonical_primitive(edges,charges,powers))


@lru_cache(maxsize=8192)
def primitive_moment(edges,charges,powers):
    if sp(edges):
        return moment(edges,charges,powers)
    roots=safe_roots(edges) or tuple(range(len(edges)))
    degrees=Counter(v for pair in edges for v in pair)
    # A charged high-valence root keeps more large inactive blocks out of its
    # projection. The ten-edge prism with a diagonal has its shortest checked
    # projection on the edge joining the two four-valent charged vertices.
    def priority(e):
        u,v=edges[e]
        return (any(charges[u])+any(charges[v]),degrees[u]+degrees[v],-e)
    for root in sorted(roots,key=priority,reverse=True):
        try:
            return coefficient(edges,charges,powers,root=root)
        except ValueError as error:
            if 'No polynomial electrical reduction' not in str(error): raise
    raise ArithmeticError((edges,charges,powers))


def eligible_shape(edges,budget):
    if not edges:return True
    degrees=Counter(v for pair in edges for v in pair)
    # At most three vertices are marked. Every other vertex consumes at
    # least this much of the remaining leg/flag degree.
    costs=sorted(((budget+1 if d==1 else 2 if d==2 else d-3) for d in degrees.values()),reverse=True)
    return sum(costs[3:])<=budget


@cache
def exceptional_faces(size):
    """Return projected roots and exact face omissions for the two cancellation channels."""
    if size not in (8,9,10): return {}
    from genus_five_nine_faces import active_face, loads
    from genus_five_graph_census import census
    rows=census(size)[1]
    faces={}
    for ch in (1,3):
        reserved=[(g,m) for g,m,aut,channels,bad in rows if ch in channels]
        roots={8:(7,),9:(3,3,2),10:(3,3,2,3,2,3)}[size]
        assert len(reserved)==len(roots)
        for (g,m),root in zip(reserved,roots):
            key=(g,loads(m,max(max(e) for e in g)+1),ch)
            if size==8:
                # The independent bank includes this entire eight-edge marked graph.
                faces[key]=None
            elif size==9:
                faces[key]=(root,(tuple(range(9)),))
            else:
                _,missing,_,_=active_face(g,m,root)
                faces[key]=(root,(tuple(e for e in range(10) if e!=missing),))
    return faces


def subtotals(size,channels=(4,),check=False,check_weights=False,*,shape_subset=None,swap_marks=False):
    levels=tuple(LEVELS[i] for i in channels)
    if swap_marks:
        if any(i in (1,3) for i in channels):
            raise ValueError('The exceptional face bank uses the original ordered momentum triangle.')
        levels=tuple(level[::-1] for level in levels)
    omitted=exceptional_faces(size) if any(i in (1,3) for i in channels) else {}
    totals=[F() for i in levels]
    counts=[0 for i in levels]
    for edges in graphs(size) if shape_subset is None else shape_subset:
        if not eligible_shape(edges,10-size): continue
        table=stencil(5,size,levels,shape_subset=(edges,))
        if check_weights:
            from planar_four_loop_assembly import independent
            assert table==independent(size,5,levels,shape_subset=(edges,)),edges
        for (g,q,p,i),w in table.items():
            exclusion=(g,q,channels[i])
            if exclusion in omitted:
                entry=omitted[exclusion]
                if entry is None:val=F()
                else:
                    root,skip=entry
                    val=coefficient(g,q,p,root=root,omit_active=skip)
            else:
                val=capacity_moment(g,q,p)
                if check:
                    assert val==moment(g,q,p),(g,q,p,val)
            totals[i]+=w*val
            counts[i]+=1
        # The program streams unlabeled shapes; mark orbits from old shapes are disposable.
        marked_shapes.cache_clear()
    for level,count,total in zip(levels,counts,totals):
        print(size,level,count,total,"polynomial charts" if omitted else "",flush=True)
    return tuple(totals)


if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--edges',nargs='+',type=int,default=range(11))
    parser.add_argument('--channels',nargs='+',type=int,choices=range(5),default=(4,),
                        help='indices in (5,0,0), (4,1,0), (3,2,0), (3,1,1), (2,2,1)')
    parser.add_argument('--swap-marks',action='store_true',help='exchange the first and last cotangent levels in polynomial channels')
    parser.add_argument('--check',action='store_true',help='also evaluate every term using unreduced moments')
    parser.add_argument('--check-weights',action='store_true',help='check every weight by a global half-edge expansion and first KdV')
    args=parser.parse_args()
    result=[F() for i in args.channels]
    for size in args.edges:
        result=[a+b for a,b in zip(result,subtotals(size,args.channels,args.check,args.check_weights,swap_marks=args.swap_marks))]
    if any(i in (1,3) for i in args.channels) and set(args.edges).intersection((8,9,10)):
        if set((8,9,10))<=set(args.edges):
            from genus_five_nine_star import summed_obstructions
            for j,i in enumerate(args.channels):
                if i in (1,3):result[j]+=summed_obstructions(i)
        else:
            print('The selected polynomial charts omit the grouped eight- through ten-edge terms.')
    if sorted(args.edges)==list(range(11)):
        for channel,value in zip(args.channels,result):
            level=LEVELS[channel]
            assert value==COMPLETE[level],(level,value,COMPLETE[level])
    print('Selected edge sum:',*result)
