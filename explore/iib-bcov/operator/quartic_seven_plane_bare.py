"""Bare four-loop seven-planar-wave integrand after summing and integrating all seven placements."""
from fractions import Fraction as F
from collections import Counter
from functools import cache
from itertools import combinations, product
from math import lcm, prod
from pathlib import Path
import sys

from flint import fmpq_mat, fmpz_mat
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from point_graph import solve_graph
from scalar_fan import cofactor
from planar_trace import matchings


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def minors(rows):
    """Alternating maximal minors of a cycle-column matrix with one extra column."""
    rank = len(rows[0])
    assert len(rows) == rank-1
    return tuple((-1)**j*int(fmpz_mat([r[:j]+r[j+1:] for r in rows]).det()) for j in range(rank))


def compressed_network(edges,cycles,lengths,normal):
    # J maps unit conserved endpoint sources (sink at vertex 0) to edge currents.
    n=len(normal)
    basis=[tuple(int(i==j)-int(i==0) for j in range(1,n)) for i in range(n)]
    eflow,_=solve_graph(edges,basis,lengths,0)
    currents,phi=solve_graph(edges,[list(map(str,row)) for row in normal],lengths,0)
    cols=list(zip(*cycles))
    rank = len(cycles)
    reduced=[[1]*len(edges)]+[
        [c[r]*y[k] for c,y in zip(cols,currents)] for k in range(3) for r in range(rank)]
    assert len(reduced)==len(edges)-1
    def signed_minor(j):
        return (-1)**j*F(str(fmpq_mat([[str(v) for v in row[:j]+row[j+1:]] for row in reduced]).det()))
    kernel=tuple(map(signed_minor,range(len(edges))))
    assert all(dot(row,kernel)==0 for row in reduced)
    potentials=[(F(),)*rank]+[
        tuple(sum(w*c[r]*flow[j] for w,c,flow in zip(kernel,cols,eflow)) for r in range(rank))
        for j in range(n-1)]
    return cols,potentials,phi


def triples(edges,lengths,cols,potentials,indices):
    # For each distinct triple, integrate its three positions and include all six label orders.
    H=[[F()]*4 for _ in range(4)]
    g=[F()]*4
    c=F()
    for chosen in combinations(indices,3):
        a=minors([cols[i] for i in chosen])
        if not any(a):
            continue
        mean=var=F()
        for i in chosen:
            left,right=edges[i]
            x,y=dot(a,potentials[left]),dot(a,potentials[right])
            mean+=(x+y)/2
            var+=(x-y)**2/12
        weight=6*prod(lengths[i] for i in chosen)
        for j in range(4):
            for k in range(4):
                H[j][k]+=weight*a[j]*a[k]
            g[j]+=weight*a[j]*mean
        c+=weight*(mean*mean+var)
    return H,g,c


def planar_sources(count,t,last,tangential):
    seventh=tuple(t) if isinstance(t,(tuple,list)) else (t,t)
    if tangential is None:
        tangential=[(0,0)]*count
        tangential[last]=tuple(-3-x for x in seventh)
    if len(tangential)!=count or tuple(sum(row[i] for row in tangential)+3+seventh[i] for i in range(2))!=(0,0):
        raise ValueError('Expected conserved tangential sources on the reduced vertices.')
    return seventh,tangential


def cofactor_square(edges,cycles,lengths,normal,t,last,eligible=None,tangential=None):
    """All seven planar positions; existing five normal vertices and times are fixed.

    The planar waves are three e1, three e2 and one covector t. A scalar t
    denotes (t,t). By default the vertex last absorbs the tangential recoil.
    'eligible' restricts the placement support for independent tests.
    """
    indices=range(len(edges)) if eligible is None else eligible
    cols,potentials,_=compressed_network(edges,cycles,lengths,normal)
    H,g,c=triples(edges,lengths,cols,potentials,indices)
    seventh,planar=planar_sources(len(normal),t,last,tangential)
    offset=[[sum(p[k]*y[j] for p,y in zip(planar,potentials)) for j in range(4)] for k in range(2)]
    def act(v):
        return [dot(row,v) for row in H]
    value=F()
    for i in indices:
        a,b=edges[i]
        x=cols[i]
        he=act(x)
        mass=dot(x,he)
        d=[y-w for y,w in zip(potentials[a],potentials[b])]
        var=2*prod(seventh)**2*(mass*dot(d,act(d))-dot(he,d)**2)/12
        z=[[s*(y+w)/2+p for y,w,p in zip(potentials[a],potentials[b],h)] for s,h in zip(seventh,offset)]
        second=[c+2*dot(g,k)+dot(k,act(k)) for k in z]
        mixed=[dot(x,g)+dot(he,k) for k in z]
        term=mass*(seventh[0]**2*second[1]+seventh[1]**2*second[0])-2*prod(seventh)*prod(mixed)
        assert term>=0 and var>=0
        value+=lengths[i]*(term+var)
    return value/sum(lengths)**2


@cache
def wedge_entries(degree):
    """Nonzero entries of left exterior multiplication in the ordered basis."""
    source=list(combinations(range(8),degree))
    target={s:i for i,s in enumerate(combinations(range(8),degree+1))}
    entries=[]
    for j,slots in enumerate(source):
        for h in range(8):
            if h not in slots:
                i=target[tuple(sorted((*slots,h)))]
                entries.append((i,j,h,(-1)**sum(k<h for k in slots)))
    return len(target),len(source),entries


def wedge(degree,vector):
    rows,cols,entries=wedge_entries(degree)
    result=fmpz_mat(rows,cols)
    for i,j,h,sign in entries:
        if vector[h]:
            result[i,j]=sign*vector[h]
    return result


def covariance_square(edges,cycles,lengths,normal,plane,tangential,eligible=None):
    """Integrate seven arbitrary planar labels and the given conserved tangential lifts."""
    if (len(plane)!=7 or len(tangential)!=len(normal) or any(map(sum,zip(*[*plane,*tangential])))):
        raise ValueError('Expected seven planar momenta and conserved tangential lifts.')
    indices=range(len(edges)) if eligible is None else eligible
    cols,potentials,_=compressed_network(edges,cycles,lengths,normal)
    if not any(any(row) for row in potentials):
        return F()
    planar_den=lcm(*(F(x).denominator for row in [*plane,*tangential] for x in row))
    plane=[tuple(int(x*planar_den) for x in row) for row in plane]
    lifts=[tuple(int(x*planar_den) for x in row) for row in tangential]
    normal_den=2*lcm(*(x.denominator for row in potentials for x in row))
    integer_k=[[int(x*normal_den) for x in row] for row in potentials]
    length_den=lcm(*(F(x).denominator for x in lengths))
    # W is the exterior product of the jumps, U=W exterior (total residual planar current).
    # On each edge both are propagated by a block matrix affine in the new position.
    initial=fmpz_mat([[1],*[[sum(q[h]*k[j] for q,k in zip(lifts,integer_k))]
                          for h in range(2) for j in range(4)]])
    gram=initial*initial.transpose()
    for k,q in enumerate(plane):
        dimw,dimu=wedge_entries(k)[0],wedge_entries(k+1)[0]
        updated=fmpz_mat(dimw+dimu,dimw+dimu)
        for e in indices:
            left,right=(integer_k[v] for v in edges[e])
            jump=[z*h for z in q for h in cols[e]]
            midpoint=[z*((a+b)//2) for z in q for a,b in zip(left,right)]
            slope=[z*(b-a) for z in q for a,b in zip(left,right)]
            w,u=wedge(k,jump),wedge(k+1,jump)
            mixed=(wedge(k+1,midpoint)*w)*((-1)**(k+1))
            change=(wedge(k+1,slope)*w)*((-1)**(k+1))
            transfer=fmpz_mat([*[row+[0]*u.ncols() for row in w.tolist()],
                                *[row+tail for row,tail in zip(mixed.tolist(),u.tolist())]])
            diff=fmpz_mat([*[[0]*transfer.ncols() for _ in range(dimw)],
                           *[row+[0]*u.ncols() for row in change.tolist()]])
            updated+=(12*(transfer*gram*transfer.transpose())+diff*gram*diff.transpose())*int(lengths[e]*length_den)
        gram=updated
    denominator=planar_den**16*normal_den**2*(12*length_den)**7*sum(lengths)**2
    return F(int(gram[gram.nrows()-1,gram.nrows()-1]))/denominator


def quadratic_lift_square(edges,cycles,lengths,normal,plane,lifts):
    """Quadratic coefficient on varying the tangential lifts along a conserved direction."""
    if any(map(sum,zip(*lifts))):
        raise ValueError('The variation of the tangential lifts must be conserved.')
    _,potentials,_=compressed_network(edges,cycles,lengths,normal)
    current=[[sum(q[h]*k[j] for q,k in zip(lifts,potentials)) for j in range(4)] for h in range(2)]
    matrix=fmpq_mat([[str(sum(z*a*b for z,a,b in zip(lengths,c,d))) for d in cycles] for c in cycles])
    gram=matrix.det()**2*matrix.inv()
    total=F()
    for i,q in enumerate(plane):
        projection=fmpq_mat([[str(q[1]*a-q[0]*b)] for a,b in zip(*current)])
        remaining=[u for j,u in enumerate(plane) if j!=i]
        total+=matchings(remaining)*F(str((projection.transpose()*gram*projection)[0,0]))
    return total/sum(lengths)**2


def placement(edges,cycles,lengths,normal,t,last,assigned,positions,tangential=None):
    """Restore all seven bivalent vertices and compare with the reduced determinant."""
    seventh,planar=planar_sources(len(normal),t,last,tangential)
    waves=[(1,0)]*3+[(0,1)]*3+[seventh]
    cols,potentials,_=compressed_network(edges,cycles,lengths,normal)
    v=[sum(p[k]*y[j] for p,y in zip(planar,potentials)) for k in range(2) for j in range(4)]
    data=[]
    for wave,i,s in zip(waves,assigned,positions):
        a,b=edges[i]
        ks=[(1-s)*x+s*y for x,y in zip(potentials[a],potentials[b])]
        data.append([h*r for h in wave for r in cols[i]])
        v=[u+d for u,d in zip(v,[h*k for h in wave for k in ks])]
    value=abs(F(str(fmpq_mat([[str(x) for x in row] for row in [*data,v]]).det())))/sum(lengths)
    full_edges,full_times,parent=[],[],[]
    sources=[(*p,*q) for p,q in zip(planar,normal)]
    sources.extend((*u,0,0,0) for u in waves)
    for i,(a,b) in enumerate(edges):
        slots=sorted((s,len(normal)+j) for j,(e,s) in enumerate(zip(assigned,positions)) if e==i)
        path=[a,*[j for _,j in slots],b]
        cuts=[F(),*[s for s,_ in slots],F(1)]
        for k,(u,w) in enumerate(zip(path,path[1:])):
            full_edges.append((u,w))
            full_times.append(lengths[i]*(cuts[k+1]-cuts[k]))
            parent.append(i)
    yy,_=solve_graph(full_edges,[list(map(str,row)) for row in sources],full_times,0)
    direct=cofactor(yy,[[c[i] for i in parent] for c in cycles],full_times)
    assert value==direct
    return direct


def literal_placements(edges,cycles,lengths,normal,plane,lifts,assignments):
    # The eighth row is affine in all seven uniform positions. Evaluate its
    # mean and each independent variance directly for each individual assignment.
    cols,potentials,_=compressed_network(edges,cycles,lengths,normal)
    offset=[sum(p[k]*y[j] for p,y in zip(lifts,potentials)) for k in range(2) for j in range(4)]
    total=F()
    for assigned,orderings in assignments:
        if max(Counter(assigned).values())>2:
            continue
        seven=[[h*r for h in wave for r in cols[i]] for wave,i in zip(plane,assigned)]
        k=[(-1)**j*int(fmpz_mat([r[:j]+r[j+1:] for r in seven]).det()) for j in range(8)]
        if not any(k):
            continue
        mean=dot(k,offset)
        variance=F()
        for i,wave in zip(assigned,plane):
            a,b=edges[i]
            ka=[x*wave[0]+y*wave[1] for x,y in zip(k[:4],k[4:])]
            left,right=dot(ka,potentials[a]),dot(ka,potentials[b])
            mean+=(left+right)/2
            variance+=(right-left)**2/12
        total+=orderings*prod(lengths[i] for i in assigned)*(mean*mean+variance)
    return total/sum(lengths)**2


def literal_moments(edges,cycles,lengths,normal,t,last,eligible,tangential=None):
    seventh,lifts=planar_sources(len(normal),t,last,tangential)
    groups=list(combinations(eligible,3))
    assignments=((first+second+(final,),36) for first in groups for second in groups for final in eligible)
    return literal_placements(edges,cycles,lengths,normal,[(1,0)]*3+[(0,1)]*3+[seventh],lifts,assignments)


def fixture():
    normal=((1,2,0),(3,-1,0),(-2,1,1),(0,2,2),(-2,-4,-3))
    basis=((1,0),(0,1),(-1,-1))
    edges,cols=[],[]
    occupied=(0,1,4,6,8)
    for i in range(3):
        for j in range(3):
            h=3*i+j
            path=[i,6+occupied.index(h),3+j] if h in occupied else [i,3+j]
            for edge in zip(path,path[1:]):
                edges.append(edge)
                cols.append([a*b for a in basis[i] for b in basis[j]])
    return edges,list(zip(*cols)),[1+(7*i*i+2*i)%17 for i in range(14)],[(0,0,0)]*6+list(normal)


def fixed_density(edges,cycles,lengths,normal,t=-3,pair=(6,7),polarization=(0,0,1),last=10,tangential=None,plane=None):
    # The polarization annihilates the two-plane.
    tree=F(str(fmpq_mat([[str(sum(z*x*y for z,x,y in zip(lengths,c,d))) for c in cycles] for d in cycles]).det()))
    _,phi=solve_graph(edges,[list(map(str,row)) for row in normal],lengths,0)
    voltage=dot(polarization,phi[pair[0]])-dot(polarization,phi[pair[1]])
    if plane is None:
        square=cofactor_square(edges,cycles,lengths,normal,t,last,tangential=tangential)
    else:
        if tangential is None:
            tangential=[(0,0)]*len(normal)
            tangential[last]=tuple(-sum(col) for col in zip(*plane))
        square=covariance_square(edges,cycles,lengths,normal,plane,tangential)
    return voltage*square/tree**5


def fixed_transfer(edges,cycles,lengths,normal):
    from fusion_transfer_box import shift
    points=Counter({tuple(normal[6:]):1})
    for i in (0,1):
        for j in (2,3,4):
            new=Counter()
            for point,weight in points.items():
                new[shift(point,i,j,(1,0,0))]+=weight
                new[point]-=weight
            points=new
    return sum(weight*fixed_density(edges,cycles,lengths,normal[:6]+list(point))
               for point,weight in points.items() if weight)


if __name__=='__main__':
    edges,cycles,lengths,normal=fixture()
    positions=tuple(F(i,9) for i in range(1,8))
    for assigned in ((0,2,8,1,5,7,3),(0,2,5,2,5,7,7)):
        for t in (-3,1):
            assert placement(edges,cycles,lengths,normal,t,10,assigned,positions)
    slots=(0,2,5,7,9)
    lift=[(0,0)]*6+[(1,0),(-2,3),(0,-1),(1,-2),(-4,-1)]
    assert placement(edges,cycles,lengths,normal,(1,-2),10,(0,2,5,2,5,7,7),positions,lift)
    mixed=cofactor_square(edges,cycles,lengths,normal,(1,-2),10,slots,lift)
    assert mixed==literal_moments(edges,cycles,lengths,normal,(1,-2),10,slots,lift)>0
    general=[(1,0),(0,1),(1,1),(2,-1),(1,2),(-1,3),(2,3)]
    arbitrary_lifts=[(0,0)]*6+[(1,0),(-2,3),(0,-1),(1,-2),(-6,-9)]
    check_edges=(0,2,5,7)
    brute=literal_placements(edges,cycles,lengths,normal,general,arbitrary_lifts,
                             ((e,1) for e in product(check_edges,repeat=7)))
    assert covariance_square(edges,cycles,lengths,normal,general,arbitrary_lifts,check_edges)==brute>0
    variation=[(0,0)]*6+[(1,-1),(0,2),(-2,2),(3,1),(-2,-4)]
    res=[covariance_square(edges,cycles,lengths,normal,general,
                           [tuple(x+k*y for x,y in zip(a,b)) for a,b in zip(arbitrary_lifts,variation)]) for k in (-1,0,1)]
    assert quadratic_lift_square(edges,cycles,lengths,normal,general,variation)==(res[0]-2*res[1]+res[2])/2>0
    for t,expected in (
        (-3,F(361526804441142177573036508129317414785560855040,1035559322835992133562676685423141457713)),
        (1,F(172006718575322401100505019811358704033216614912,9320033905523929202064090168808273119417)),
    ):
        assert cofactor_square(edges,cycles,lengths,normal,t,10,slots)==literal_moments(edges,cycles,lengths,normal,t,10,slots)>0
        x=cofactor_square(edges,cycles,lengths,normal,t,10)
        assert x==expected
        seventh,ls=planar_sources(len(normal),t,10,None)
        assert covariance_square(edges,cycles,lengths,normal,[(1,0)]*3+[(0,1)]*3+[seventh],ls)==x
        print(t,'all seven placements, fixed-time cofactor square',x,flush=True)
    density=fixed_density(edges,cycles,lengths,normal)
    assert density == F(-643589294212512422231164445698794603166425863020497920,
                        271268107591917906878240705413980381876127749731946347799257191725718629099)
    assert fixed_density(edges,cycles,[F(z,2) for z in lengths],normal)==2**14*density
    transfer=fixed_transfer(edges,cycles,lengths,normal)
    assert transfer == F(396915990365908696602971048851299327965554795336359510016,
                         2441412968327261161904166348725823436885149747587517130193314725531467661891)
    print('fixed-time bare sixth transfer density',transfer)
