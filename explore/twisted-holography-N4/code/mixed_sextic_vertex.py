"""Six mixed Betti factors, including the full tetrahedron period in two channels."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb,factorial,prod
from operator import index
from betti_power_graphs import graph_terms,weighted_vertex
from betti_series_parallel import coefficient as face_integral
from mixed_k4_interior import interior


@lru_cache(maxsize=512)
def tetrahedron(charges,root=0):
    edges=tuple(combinations(range(4),2))
    proper=face_integral(edges,charges,(0,)*6,root=root,proper=True)
    # A single current direction gives a length Hessian of rank at most three.
    pivot=next((row for row in charges if any(row)),None)
    if pivot is None: return F(),F()
    j=next(i for i,v in enumerate(pivot) if v)
    if all(row[i]*pivot[j]==row[j]*c for row in charges for i,c in enumerate(pivot)): return proper,F()
    rational,zeta=interior(charges,root)
    return proper+rational,zeta


def partitions(count):
    """Place labels at up to four unordered vertices, keeping only nonempty groups."""
    if not count:
        yield ()
        return
    for groups in partitions(count-1):
        for i,group in enumerate(groups):
            yield (*groups[:i],(*group,count-1),*groups[i+1:])
        if len(groups)<4: yield (*groups,(count-1,))


def sextic(genus,charges,levels):
    """Pair Theta(charges)^6/6! against standard mark cotangents on the compact moduli space."""
    try:
        genus=index(genus)
        levels=tuple(map(index,levels))
    except TypeError as exc:
        raise ValueError('The genus and cotangent degrees must be integers.') from exc
    charges=tuple(tuple(row) for row in charges)
    n=len(charges)
    if (genus<0 or n!=len(levels) or not n or len(set(map(len,charges)))!=1 or any(map(sum,zip(*charges)))
        or any(k<0 for k in levels) or sum(levels)!=3*genus-3+n-6):
        raise ValueError('Check conserved charge vectors and complementary cotangent degree.')
    # Conservation removes one independent source; each target contributes at most genus Betti pairs.
    if genus*min(len(charges[0]),n-1)<6:
        return F(),F()
    a=b=F()
    for groups,flags,extra,weight in graph_terms(6,tuple(charges),genus,omit_tetrahedron=True):
        shifted=tuple(x+y for x,y in zip(extra,levels))
        for group,branch in zip(groups,flags):
            weight*=weighted_vertex(group,branch,shifted)
            if not weight: break
        a+=weight
    if genus<3: return a,b
    for occupied in partitions(n):
        empty=4-len(occupied)
        groups=(*occupied,*((),)*empty)
        weight=prod(weighted_vertex(g,(0,0,0),levels) for g in groups)/factorial(empty)
        if not weight: continue
        q=tuple(tuple(sum(charges[i][j] for i in g) for j in range(len(charges[0]))) for g in groups)
        r,z=tetrahedron(tuple(sorted(q)))
        a+=weight*r
        b+=weight*z
    return a,b


def labeled_sum(genus,charges,levels):
    """Independent all-label sum, with explicit vertex genera and half-edge degrees."""
    from betti_flow_polytope import weak
    from betti_power_graphs import shapes,vector_coefficient
    from soft_quartic_psi_bound import intersection
    n=len(charges)
    total=[F(),F()]
    for size in range(7):
        for vertices,edges,aut,edge_aut in shapes(size):
            left=genus-(size-vertices+1)
            if left<0: continue
            ends=[[e for e,pair in enumerate(edges) if v in pair] for v in range(vertices)]
            for assignment in product(range(vertices),repeat=n):
                groups=[[i for i,v in enumerate(assignment) if v==j] for j in range(vertices)]
                slots=[len(g)+len(ee) for g,ee in zip(groups,ends)]
                staple=tuple(tuple(sum(charges[i][j] for i in g) for j in range(len(charges[0]))) for g in groups)
                for genera in weak(left,vertices):
                    missing=[3*h-3+k-sum(levels[i] for i in g) for h,k,g in zip(genera,slots,groups)]
                    if min(missing)<0 or any(2*h-2+k<=0 for h,k in zip(genera,slots)): continue
                    expansions=[]
                    for g,flags,h,d in zip(groups,ends,genera,missing):
                        choices=[]
                        for deg in weak(d,len(g)+len(flags)):
                            legs=deg[:len(g)]
                            branches=deg[len(g):]
                            t=intersection(h,(*[levels[i]+k for i,k in zip(g,legs)],*branches),last=True)
                            for i,k in zip(g,legs):
                                t*=F(sum(v*v for v in charges[i])**k,2**k*factorial(k))
                            if t: choices.append((t,dict(zip(flags,branches))))
                        expansions.append(choices)
                    for terms in product(*expansions):
                        weight=prod(v[0] for v in terms)/F(len(aut)*edge_aut)
                        normal=[]
                        for e,(i,j) in enumerate(edges):
                            u,v=terms[i][1][e],terms[j][1][e]
                            normal.append(u+v)
                            weight*=comb(u+v,u)
                        if size==6 and vertices==4 and len(set(edges))==6:
                            coeff=tetrahedron(staple,root=5)
                        else:
                            coeff=vector_coefficient(edges,staple,tuple(normal)),F()
                        for ch,c in enumerate(coeff): total[ch]+=weight*c
    return tuple(total)


def check(explicit=False):
    q=((1,0),(0,1),(-1,-1),(0,0))
    values=[tetrahedron(q,r) for r in range(6)]
    assert len(set(values))==1 and values[0]==(-F(4657,217728000),0)
    print('Rank-two tetrahedron:',values[0],flush=True)
    q=q[:3]
    control=((2,1),(1,2),(-3,-3))
    fixtures=(((3,0,0),F(7,10368000)),((2,1,0),F(17,7257600)),((1,1,1),F(11,2177280)))
    for levels,wanted in fixtures:
        value=sextic(3,q,levels)
        assert value==(wanted,0)
        assert sextic(3,control,levels)==(3**6*wanted,0)
        if explicit: assert labeled_sum(3,q,levels)==value
        print('Genus-three two-target vertex',levels,*value,flush=True)

    assert sextic(2,q,(0,0,0))==(0,0)
    if explicit: assert labeled_sum(2,q,(0,0,0))==(0,0)
    spectator=(*q,(0,0))
    for levels,wanted in fixtures:
        assert sextic(3,spectator,(*levels,1))==(7*wanted,0)
    for levels in ((4,0,0),(3,1,0),(2,1,1)):
        wanted=F()
        for i,a in enumerate(levels):
            if a:
                lower=tuple(k-int(j==i) for j,k in enumerate(levels))
                wanted+=dict(fixtures)[tuple(sorted(lower,reverse=True))]
        assert sextic(3,spectator,(*levels,0))==(wanted,0)
    four=((1,0),(0,1),(1,2),(-2,-3))
    value=sextic(3,four,(1,1,1,1))
    assert value==(F(545359,33530112),F(2,1617))
    print('Genus-three four-mark positive zeta(3) control:',*value,flush=True)


if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--labeled',action='store_true',help='Replay all three vertices with explicit genera and labeled placements')
    check(parser.parse_args().labeled)
