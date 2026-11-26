"""Three-dimensional theta integrals remaining after four saturated planar jumps."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations, product
from math import factorial, lcm, prod
from pathlib import Path
import sys

from flint import fmpz_mat, fmpz_mpoly_ctx
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bulk_eight import BULK
from marked_boundaries import quartics
from planar_trace import matchings


CTX = fmpz_mpoly_ctx.get(tuple('s'+str(i) for i in range(4))+tuple('w'+str(i) for i in range(3)))
S = CTX.gens()[:4]
W = CTX.gens()[4:]
DEN = factorial(18)**4 * factorial(11)
NORMAL = ((1,2,0),(3,-1,0),(-2,1,1),(0,2,2),(-2,-4,-3))


def determinant(u,v,w):
    return (u[0]*(v[1]*w[2]-v[2]*w[1])
           +u[1]*(v[2]*w[0]-v[0]*w[2])
           +u[2]*(v[0]*w[1]-v[1]*w[0]))


def chambers():
    """Quotient the three paths by permutation; retain orders of the labeled positions."""
    for path in product(range(3),repeat=4):
        if tuple(sorted(set(path),key=path.index))!=tuple(range(max(path)+1)):
            continue
        groups=tuple(tuple(i for i in range(4) if path[i]==k) for k in range(3))
        if max(map(len,groups))>2:
            continue
        for ordered in product(*(permutations(g) for g in groups)):
            yield path,ordered


def cofactor(p,path):
    zero=(0,0,0)
    columns=[(q if h!=2 else zero)+(q if h!=0 else zero) for q,h in zip(p,path)]
    coef=CTX.constant(0)
    for i,j in permutations(range(4),2):
        if path[i]!=path[j]:
            m=int(fmpz_mat([p[i]+zero,zero+p[j],*columns]).det())
            coef+=m*S[i]*S[j]
    return coef


@cache
def moments(ordered,k,shift):
    # The pullback to conductance fractions has density product w_h^(3-occupancy).
    # Ordered position monomials and the conductance simplex each give a factorial ratio.
    ss=k[:4]
    ww=k[4:]
    den=1
    for group in ordered:
        total=0
        for i in group:
            total+=ss[i]+1
            den*=total
    deg=[3-len(g)+ww[h]-int(h==shift) for h,g in enumerate(ordered)]
    den*=factorial(2+sum(deg))
    num=DEN*prod(factorial(a) for a in deg)
    assert num%den==0
    return int(num//den)


def integrate_chamber(poly,ordered,shift=-1):
    return sum(int(c)*moments(ordered,tuple(map(int,k)),shift) for k,c in poly.terms())


def vertex_determinant(p,mark,kind,ordered,path,l):
    if kind=='branch':
        # Incident currents on paths 0,1 at the grounded endpoint.
        h=[[sum(((1-S[i])*p[i][j] for i in g),CTX.constant(0))
            for j in range(3)] for g in ordered[:2]]
        vert=(determinant(mark,h[0],h[1])+W[0]*determinant(mark,l,h[1])
              +W[1]*determinant(mark,h[0],l))
    else:
        # Event zero is the bivalent trace. Its jump is mark[0]+mark[1].
        hp=path[0]
        indices=ordered[hp]
        mixed=[determinant(*mark,q) for q in p]
        vert=W[hp]*sum(S[i]*x for i,x in enumerate(mixed))
        for pos,i in enumerate(indices):
            vert+=(int(pos>indices.index(0))-S[i])*mixed[i]
    return vert


def terms(p,mark,external_sign,kind,direction):
    # B = ground (s=0), A = other endpoint (s=1), sum w=1.
    l=[sum(S[i]*q[j] for i,q in enumerate(p)) for j in range(3)]
    gamma=[sum(x*y for x,y in zip(q,direction)) for q in p]
    global_v=sum(S[i]*x for i,x in enumerate(external_sign))*sum(S[i]*x for i,x in enumerate(gamma))
    base=0
    for path,ordered in chambers():
        quad=cofactor(p,path)
        if not quad:
            continue
        vert=vertex_determinant(p,mark,kind,ordered,path,l)
        core=quad**2*vert**4
        if not core:
            continue
        base+=integrate_chamber(core*global_v,ordered)
        for h,g in enumerate(ordered):
            local=CTX.constant(0)
            for ai,i in enumerate(g):
                for aj,j in enumerate(g):
                    local+=external_sign[i]*gamma[j]*(S[g[min(ai,aj)]]-S[i]*S[j])
            if local:
                base+=integrate_chamber(core*local,ordered,h)
    return Fraction(base,DEN)


def normal_component(normal,kind='both',target=(0,1),direction=(0,0,1)):
    """Difference of two normal potentials, polarized in a normal direction."""
    if (len(normal)!=5 or any(len(q)!=3 for q in normal) or any(map(sum,zip(*normal)))
            or len(direction)!=3 or len(target)!=2 or any(i not in range(5) for i in target)
            or kind not in ('both','branch','interior')):
        raise ValueError('Expected five conserved three-component waves, a target pair, and a direction.')
    scale=lcm(*(Fraction(v).denominator for q in normal for v in q))
    normal=tuple(tuple(int(v*scale) for v in q) for q in normal)
    vscale=lcm(*(Fraction(v).denominator for v in direction))
    direction=tuple(int(v*vscale) for v in direction)
    first,second=target
    total=Fraction()
    if kind in ('both','branch'):
        for n,r in enumerate(normal):
            ids=tuple(i for i in range(5) if i!=n)
            signs=[int(i==first)-int(i==second) for i in ids]
            total+=terms([normal[i] for i in ids],r,signs,'branch',direction)
    if kind in ('both','interior'):
        for i,j in combinations(range(5),2):
            merged=tuple(a+b for a,b in zip(normal[i],normal[j]))
            ids=tuple(k for k in range(5) if k not in (i,j))
            signs=[int(first in (i,j))-int(second in (i,j))]+[int(k==first)-int(k==second) for k in ids]
            total+=terms([merged,*[normal[k] for k in ids]],(normal[i],normal[j]),signs,'interior',direction)/2
    return total/(vscale*scale**25)


def planar_factor(plane):
    axes=[(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)]
    value=Fraction()
    for group in combinations(range(7),3):
        rows=[plane[i]+(0,0,0) for i in group]+axes
        rows.append(tuple(-sum(col) for col in zip(*rows)))
        weight=-sum(x*y for x,y in zip(BULK,quartics(rows)))
        # matchings already counts the two orders of the pair factors.
        value+=weight*matchings([q for i,q in enumerate(plane) if i not in group])
    return value


def box_normal(normal):
    from fusion_transfer_box import shift
    totals=Counter()
    edges=[(i,j) for i in (0,1) for j in (2,3,4)]
    for mask in range(64):
        point=normal
        for k,(i,j) in enumerate(edges):
            if mask>>k&1:
                point=shift(point,i,j,(1,0,0))
        totals[tuple(point)]+=(-1)**mask.bit_count()
    return tuple(sum(sign*normal_component(point,kind) for point,sign in totals.items() if sign)
                 for kind in ('interior','branch'))


def evaluate(poly,x):
    return sum(int(c)*prod(z**int(n) for z,n in zip(x,k)) for k,c in poly.terms())


def check_planar_placement():
    paths=((1,0),(-1,-1),(0,1))
    lengths=(2,3,5)
    kirchhoff=sum(x*y for x,y in combinations(lengths,2))
    for waves,expected in ((((1,0),(0,1),(2,-1),(1,2)),80724),
                           (((1,0),(1,0),(0,1),(0,1)),3844)):
        weight=0
        for assignment in product(range(3),repeat=4):
            cols=[[c*x for c in paths[e] for x in q] for e,q in zip(assignment,waves)]
            weight+=prod(lengths[e] for e in assignment)*int(fmpz_mat(cols).det())**2
        assert weight==kirchhoff**2*matchings(waves)==expected


def check_network(kind):
    from point_graph import solve_graph
    from scalar_fan import cofactor as graph_cofactor
    # At the trace use the oriented currents of a grounded network with unequal lengths.
    path=(0,1,0,2)
    ordered=((0,2),(1,),(3,))
    pos=[Fraction(i,7) for i in (2,3,5,4)]
    w=[Fraction(i,10) for i in (2,3,5)]
    if kind=='branch':
        events=list(NORMAL[:4])
        mark=NORMAL[4]
        source=[mark,(0,0,0),*events]
        at=0
    else:
        mark=(NORMAL[0],NORMAL[3])
        events=[tuple(map(sum,zip(*mark))),NORMAL[1],NORMAL[2],NORMAL[4]]
        source=[(0,0,0),(0,0,0),*events]
        at=2
    edges,times,parents=[],[],[]
    for h,g in enumerate(ordered):
        vs=[0,*[2+i for i in g],1]
        s=[0,*[pos[i] for i in g],1]
        for j,(a,b) in enumerate(zip(vs,vs[1:])):
            edges.append((a,b))
            times.append((s[j+1]-s[j])/w[h])
            parents.append(h)
    currents,phi=solve_graph(edges,source,times,0)
    cycles=[[int(h==k)-int(h==1) for h in parents] for k in (0,2)]
    l=[sum(S[i]*q[j] for i,q in enumerate(events)) for j in range(3)]
    q=evaluate(cofactor(events,path),pos+w)
    u=sum(1/(w[i]*w[j]) for i,j in combinations(range(3),2))
    assert graph_cofactor(currents,cycles,times)==abs(q)/(prod(w)*u**2) and q
    y=[tuple((int(a==at)-int(b==at))*v for v in row) for (a,b),row in zip(edges,currents) if at in (a,b)]
    test=determinant(mark,y[0],y[1]) if kind=='branch' else determinant(*mark,y[0])
    computed=evaluate(vertex_determinant(events,mark,kind,ordered,path,l),pos+w)
    assert computed**4==test**4 and test
    for i in range(4):
        for k in range(3):
            expected=sum(events[j][k]*(pos[i]*pos[j]+
                (min(pos[i],pos[j])-pos[i]*pos[j])/w[path[i]]*(path[i]==path[j])) for j in range(4))
            assert phi[2+i][k]==expected
    # Restore all four planar ordinary waves on distinct pieces of the sampled network.
    planar=((1,0),(1,0),(0,1),(0,1))
    grouped=((1,0),(0,1),(-3,-3))
    def jump_determinant(assignment):
        rows=[[c[e]*x for c in cycles for x in q] for e,q in zip(assignment,planar)]
        return abs(int(fmpz_mat(rows).det()))
    labels=next(es for es in permutations(range(len(edges)),4) if jump_determinant(es))
    deleted=jump_determinant(labels)
    src=[(0,0,*z) for z in source]
    src[at]=(*map(sum,zip(*grouped)),*source[at])
    src.extend((*p,0,0,0) for p in planar)
    full_edges,full_times,which=[],[],[]
    for e,(a,b) in enumerate(edges):
        if e in labels:
            node=6+labels.index(e)
            full_edges.extend(((a,node),(node,b)))
            full_times.extend((times[e]/3,2*times[e]/3))
            which.extend((e,e))
        else:
            full_edges.append((a,b))
            full_times.append(times[e])
            which.append(e)
    flow,full_phi=solve_graph(full_edges,src,full_times,0)
    assert graph_cofactor(flow,[[c[e] for e in which] for c in cycles],full_times)==deleted*graph_cofactor(currents,cycles,times)
    assert all(full_phi[i][2:]==phi[i] for i in range(6))
    inside=[tuple((int(b==at)-int(a==at))*v for v in row) for (a,b),row in zip(full_edges,flow) if at in (a,b)]
    marked_waves=[(*v,0,0,0) for v in grouped]+[(0,0,*x) for x in ([mark] if kind=='branch' else mark)]+inside
    scale=prod(z.denominator for row in inside for z in row)
    marked_waves=[tuple(int(scale*v) for v in q) for q in marked_waves]
    reference=[(*v,0,0,0) for v in grouped]+[(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)]
    reference.append(tuple(-sum(col) for col in zip(*reference)))
    assert quartics(marked_waves)==tuple(v*computed**4*scale**20 for v in quartics(reference))


if __name__=='__main__':
    assert len(list(chambers()))==24
    check_network('branch')
    check_network('interior')
    check_planar_placement()
    plane=((1,0),)*3+((0,1),)*3+((-3,-3),)
    assert planar_factor(plane)==Fraction(17,160)
    raw=tuple(normal_component(NORMAL,k) for k in ('interior','branch'))
    assert raw==(Fraction(1864769,81000),Fraction(4213502549,170100000))
    assert sum(raw)*planar_factor(plane)==Fraction(138201796633,27216000000)
    print('physical theta',sum(raw)*planar_factor(plane),flush=True)
    if '--box' in sys.argv:
        diff=box_normal(NORMAL)
        assert diff==(Fraction(-175476439,220500),Fraction(825565581019,2546775000))
        assert sum(diff)*planar_factor(plane)==Fraction(-20420183920327,407484000000)
        print('interior and endpoint normal differences',diff,'physical theta difference',sum(diff)*planar_factor(plane))
