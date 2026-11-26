"""Compact cotangent tests of exponential powers of a scalar Abel-Jacobi Betti form."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, factorial, prod

from betti_flow_polytope import charge_cuts, connected, normal_coefficient, weak
from soft_quartic_psi_bound import intersection
from soft_quartic_two_betti import LINKS, descendant_vertex


@lru_cache(maxsize=None)
def shapes(k):
    if not k: return ((1,(),((0,),),1),)
    result=[]
    for v in range(2,k+2):
        labels=tuple(range(v))
        perms=tuple(permutations(labels))
        seen=set()
        for edges in combinations_with_replacement(combinations(labels,2),k):
            if edges in seen or not connected(set(labels),edges): continue
            autos=[]
            for perm in perms:
                image=tuple(sorted(tuple(sorted((perm[a],perm[b]))) for a,b in edges))
                seen.add(image)
                if image==edges: autos.append(perm)
            result.append((v,edges,tuple(autos),prod(factorial(m) for m in Counter(edges).values())))
    return tuple(result)


def norm2(p):
    return sum(v*v for v in p) if isinstance(p,tuple) else p*p


@lru_cache(maxsize=16384)
def vector_coefficient(edges,charges,powers):
    """Mixed Betti moment on trees of parallel bundles, unicyclic simplifications and the diamond."""
    if not edges: return F(1)
    bundles = sorted(set(edges))
    if len(bundles)==len(charges)-1:
        weight=F(1)
        for edge in bundles:
            side={edge[0]}
            for _ in range(len(charges)-1):
                for a,b in bundles:
                    if (a,b)==edge: continue
                    if a in side: side.add(b)
                    if b in side: side.add(a)
            q=tuple(map(sum,zip(*(charges[i] for i in side))))
            m=[powers[i] for i,e in enumerate(edges) if e==edge]
            d=len(m)+sum(m)
            weight *= -F(norm2(q)**d*prod(factorial(2*k+1) for k in m),
                         factorial(2*d)*prod(2**k*factorial(k) for k in m))
        return weight
    if len(edges)==3 and not any(powers) and len(charges)==3:
        a,b,c=map(norm2,charges)
        return -F((a+b+c)**3,960)+F((a+b+c)*(a*b+a*c+b*c),120)-F(7*a*b*c,240)
    if len(bundles)==len(edges)==5 and len(charges)==4:
        from mixed_diamond_moment import graph_value
        return graph_value(edges,charges,powers)
    if len(bundles)!=len(charges):
        from betti_series_parallel import coefficient
        return coefficient(edges,charges,powers)
    from cycle_height_moments import cycle_value, merge_parallel
    active = {}
    weight=F(1,prod(factorial(k) for k in powers))
    for edge in bundles:
        k,factor = merge_parallel(tuple(p for e,p in zip(edges,powers) if e==edge))
        active[edge]=k
        weight*=F(factor)
    neighbors={v:set() for v in range(len(charges))}
    for u,v in active:
        neighbors[u].add(v)
        neighbors[v].add(u)
    loads=list(charges)
    while any(len(n)==1 for n in neighbors.values()):
        leaf=next(v for v,n in neighbors.items() if len(n)==1)
        parent,=neighbors.pop(leaf)
        neighbors[parent].remove(leaf)
        k=active.pop(tuple(sorted((leaf,parent))))
        weight*=-F(norm2(loads[leaf])**(k+1),(2**(k+1))*(k+1))
        if not weight: return weight
        loads[parent]=tuple(x+y for x,y in zip(loads[parent],loads[leaf]))
    start=min(neighbors)
    v=start
    prev=None
    offset=(0,)*len(charges[0])
    offsets=[]
    cycle_powers=[]
    while True:
        nxt=min(neighbors[v]-{prev})
        offset=tuple(x+y for x,y in zip(offset,loads[v]))
        offsets.append(offset)
        cycle_powers.append(active[tuple(sorted((v,nxt)))])
        prev,v=v,nxt
        if v==start: break
    return weight*cycle_value(offsets,tuple(cycle_powers))


@lru_cache(maxsize=6)
def graph_terms(d,charges,max_genus=None,variation=None,omit_tetrahedron=False):
    """Build one stable-graph integration stencil. Vertex genera follow from dimension."""
    n=len(charges)
    terms=Counter()
    if variation is not None:
        from betti_first_variation import normal_derivative
    for e in range(d+1):
        for v,edges,autos,edge_aut in shapes(e):
            if max_genus is not None and len(edges)-v+1>max_genus: continue
            if omit_tetrahedron and e==6 and v==4 and len(set(edges))==6: continue
            valences=tuple(sum(i in e for e in edges) for i in range(v))
            for labels in product(range(v),repeat=n):
                images={tuple(perm[j] for j in labels) for perm in autos}
                if min(images)!=labels: continue
                sym=F(len(images),len(autos)*edge_aut)
                groups=tuple(tuple(i for i,j in enumerate(labels) if j==vertex) for vertex in range(v))
                if max_genus is not None:
                    minimum=sum(max(0,(2-valence-len(group))//2+1)
                                for valence,group in zip(valences,groups))
                    if minimum+e-v+1>max_genus: continue
                if variation is not None:
                    base=tuple(sum(charges[i] for i in group) for group in groups)
                    tangent=tuple(sum(variation[i] for i in group) for group in groups)
                    local=lambda ex:normal_derivative(edges,base,tangent,ex)
                elif isinstance(charges[0],tuple):
                    vertex_charge=tuple(tuple(sum(charges[i][j] for i in group) for j in range(len(charges[0])))
                                        for group in groups)
                    local=lambda ex:vector_coefficient(edges,vertex_charge,ex)
                else:
                    vertex_charge=tuple(sum(charges[i] for i in group) for group in groups)
                    if e and not any(vertex_charge): continue
                    if d>=5:
                        from scalar_flow_recursion import coefficient
                        local=lambda ex:coefficient(edges,vertex_charge,ex)
                    else:
                        cuts=charge_cuts(edges,vertex_charge)
                        if e and (not cuts or len(set().union(*(set(c) for c,_ in cuts)))!=e): continue
                        local=lambda ex:normal_coefficient(cuts,ex)
                for branch_degree in range(d-e+1):
                    for normal in weak(branch_degree,e):
                        weight=sym*local(normal)
                        if not weight: continue
                        flag_choices=Counter()
                        for selected in product(*(range(b+1) for b in normal)):
                            flags=[[] for _ in groups]
                            multiplier=1
                            for (a,b),k,l in zip(edges,normal,selected):
                                flags[a].append(l)
                                flags[b].append(k-l)
                                multiplier*=comb(k,l)
                            flag_choices[tuple(tuple(sorted(f)) for f in flags)]+=multiplier
                        leg_degree=d-e-branch_degree
                        for extra in weak(leg_degree,n):
                            leg_factor=F(prod(norm2(c)**k for c,k in zip(charges,extra)),
                                         2**leg_degree * prod(factorial(k) for k in extra))
                            if leg_factor:
                                for flags,mult in flag_choices.items():
                                    terms[groups,flags,extra]+=mult*leg_factor*weight
    return tuple((g,f,k,w) for (g,f,k),w in terms.items() if w)


@lru_cache(maxsize=262144)
def weighted_vertex(group, flags, powers):
    degrees=(*flags,*(powers[i] for i in group))
    g,mod=divmod(sum(degrees)+3-len(degrees),3)
    return F() if mod or g<0 else intersection(g,degrees)


def scalar_power(genus,d,charges,degrees):
    """Pair Theta(charges)^d/d! with the indicated compact cotangent classes."""
    unbalanced = any(map(sum,zip(*charges))) if isinstance(charges[0],tuple) else sum(charges)
    if (unbalanced or len(charges)!=len(degrees) or any(k<0 for k in degrees)
        or sum(degrees)+d!=3*genus-3+len(charges)):
        raise ValueError('Check the momenta and the complementary cotangent degree.')
    result=F()
    for groups,flags,extra,weight in graph_terms(d,tuple(charges),genus):
        powers=tuple(x+y for x,y in zip(extra,degrees))
        for group,branch in zip(groups,flags):
            weight*=weighted_vertex(group,branch,powers)
            if not weight: break
        result+=weight
    return result


def retained_bound(pair, d):
    """Keep d copies of the same Betti factor; expand the other cotangent averages."""
    degrees=Counter({(0,)*5:1})
    for i,j,m in LINKS:
        count=m-(d if (i,j)==pair else 0)
        if count<0: raise ValueError('Not enough occurrences of this Betti form')
        updated=Counter()
        for powers,weight in degrees.items():
            for a in range(count+1):
                v=list(powers)
                v[i]+=a
                v[j]+=count-a
                updated[tuple(v)]+=comb(count,a)*weight
        degrees=updated
    charges=tuple(int(k==pair[0])-int(k==pair[1]) for k in range(5))
    total=sum(weight*scalar_power(4,d,charges,powers) for powers,weight in degrees.items())
    return 2*factorial(d)*total/(2**(14-d)*prod(factorial(m) for i,j,m in LINKS))


def check():
    from one_target_dr import homogeneous_dr
    for charges in ((1,-1),(2,-3,1),(2,-3,0,1)):
        n=len(charges)
        for powers in weak(n+1,n):
            actual=scalar_power(2,2,charges,powers)
            assert actual==descendant_vertex(2,tuple((p,) for p in charges),powers)
        for genus in (3,4):
            if genus==4 and n>3: continue
            dr=homogeneous_dr(genus,charges).to_dict()
            for powers in weak(2*genus-3+n,n):
                actual=scalar_power(genus,genus,charges,powers)
                assert F(str(dr.get(powers,0)))==actual,(genus,charges,powers,dr.get(powers,0),actual)
        print('All cotangent distributions verified on',charges,flush=True)
    old=retained_bound((0,1),2)
    assert old==F(19566009209,2191721811148800)
    print('Validated two-factor upper bound',old,flush=True)
    for pair in ((0,1),(2,4)):
        for d in (3,4):
            print(pair,d,retained_bound(pair,d),flush=True)


if __name__=='__main__':
    check()
