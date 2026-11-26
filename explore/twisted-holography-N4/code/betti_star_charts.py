"""Polynomial unit-cube charts for graphs admitting elimination of empty electrical stars."""
from functools import lru_cache
from itertools import combinations
from math import prod

from betti_series_parallel import currents, network, vertices


@lru_cache(maxsize=8192)
def plan(edges, empty):
    """Use oriented original edges and retain every independent resistance ratio."""
    _,source,sink=edges[0]
    try:
        return ('terminal',source,network(edges,source,sink))
    except ValueError:
        pass
    fresh=max(e for e,u,v in edges)+1
    for first,second in combinations(edges,2):
        i,u,v=first
        j,a,b=second
        if {u,v}=={a,b}:
            remaining=tuple(e for e in edges if e not in (first,second))+((fresh,u,v),)
            child=plan(remaining,empty)
            if child:
                return ('parallel',i,j,1 if a==u else -1,fresh,child)
            break
    for center in sorted(empty & vertices(edges)):
        star=tuple(e for e in edges if center in e[1:])
        neighbors=tuple(v if center==u else u for i,u,v in star)
        if len(star)!=3 or len(set(neighbors))!=3:
            continue
        triangle=tuple((fresh+i,a,b) for i,(a,b) in enumerate(combinations(neighbors,2)))
        remaining=tuple(e for e in edges if e not in star)+triangle
        child=plan(remaining,empty-{center})
        if child:
            return ('star',center,star,triangle,child)
    for center in sorted(vertices(edges)):
        series=tuple(e for e in edges if center in e[1:])
        if len(series)!=2: continue
        (i,u,v),(j,w,z)=series
        a,b=(v if u==center else u),(z if w==center else w)
        if a==b: continue
        remaining=tuple(e for e in edges if e not in series)+((fresh,a,b),)
        changed=empty-{center} if center in empty else empty-{center,a,b}
        child=plan(remaining,changed)
        if child:
            return ('series',center,a,b,i,1 if u==a else -1,j,1 if w==center else -1,fresh,child)
    return None


def solve(chart,charges,parameter,zero):
    """Recover every original edge current in the independent cube variables."""
    kind=chart[0]
    if kind=='terminal':
        _,source,tree=chart
        return currents(tree,charges,parameter,charges[source],zero)
    if kind=='parallel':
        _,i,j,sign,fresh,child=chart
        result=solve(child,charges,parameter[1:],zero)
        net=result.pop(fresh)
        t=parameter[0]
        result[i]=[t*x for x in net]
        result[j]=[sign*(1-t)*x for x in net]
    elif kind=='series':
        _,center,a,b,i,si,j,sj,fresh,child=chart
        t=parameter[0]
        load=charges[center]
        effective=dict(charges)
        effective[a]=[x+t*y for x,y in zip(charges[a],load)]
        effective[b]=[x+(1-t)*y for x,y in zip(charges[b],load)]
        result=solve(child,effective,parameter[1:],zero)
        net=result.pop(fresh)
        result[i]=[si*(x-t*y) for x,y in zip(net,load)]
        result[j]=[sj*(x+(1-t)*y) for x,y in zip(net,load)]
    else:
        _,center,star,triangle,child=chart
        result=solve(child,charges,parameter,zero)
        for i,u,v in star:
            neighbor=v if u==center else u
            direction=1 if neighbor==u else -1
            net=[zero]*len(charges[center])
            for k,a,b in triangle:
                if neighbor in (a,b):
                    outward=1 if neighbor==a else -1
                    net=[x+outward*y for x,y in zip(net,result[k])]
            result[i]=[direction*x for x in net]
        for k,a,b in triangle: del result[k]
    return result


def face_sign(jacobian):
    """The projected energy Jacobian has a fixed sign on the positive cube."""
    from flint import fmpq
    dimension=jacobian.context().nvars()
    middle=jacobian(*([fmpq(1,2)]*dimension))
    if not middle:
        middle=sum(c/prod(i+1 for i in exp) for exp,c in jacobian.terms())
    assert middle
    return -1 if middle>0 else 1
