"""Seven-gauge single-trace Ward equation with a moving divergence-free trivector."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from external_closed import stable_colour,relabel,parity


@cache
def alternant(entries):
    result=Counter()
    for order in product(*entries):
        if len(set(order))==len(order):
            result[tuple(sorted(order))]+=parity(order)
    return {s:c for s,c in result.items() if c}


def add(out,left,right,coefficient):
    for a,x in alternant(left).items():
        for b,y in alternant(right).items():
            out[a,b]+=coefficient*x*y


def wheel(rank,target=(tuple(range(7)),)):
    result=Counter()
    for (power,words),weight in stable_colour(3).items():
        for s,sign in relabel(words,target):
            p=[(i,) for i in s]
            for i in range(5):
                add(result,tuple(p[:5]),(p[5],p[6],p[i]),(5-i)*weight*rank**power*sign)
    return result


def cyclic_primitive(left,right):
    out=Counter()
    for j in range(7):
        order=[(j+i)%7 for i in range(7)]
        momenta=[tuple(order[:2])]+[(i,) for i in order[2:]]
        for direction in [-1,1]:
            for first in range(6):
                s=[(first+i*direction)%6 for i in range(6)]
                weight=-parity(order)*parity(s)
                add(out,tuple(momenta[s[i]] for i in left),tuple(momenta[s[i]] for i in right),weight)
    return out


@cache
def coordinates(slots):
    # p0,...,p4 are the five coordinate covectors; p5,p6 have independent entries.
    out=Counter()
    for components in product(*(range(5) if i>=5 else (i,) for i in slots)):
        if len(set(components))!=len(components):
            continue
        variables=tuple(sorted((i,j) for i,j in zip(slots,components) if i>=5))
        out[tuple(sorted(components)),variables]+=parity(components)
    return {key:value for key,value in out.items() if value}


@cache
def actual_column(a,b):
    # Every trivector transverse to q is contraction of a four-vector with q for q nonzero.
    result=Counter()
    for (_,left),x in coordinates(a).items():
        for i in range(7):
            if i in b:
                continue
            for (rows,right),y in coordinates((i,*b)).items():
                result[rows,tuple(sorted(left+right))]-=x*y
    return {key:value for key,value in result.items() if value}


def actual(symbol):
    result=Counter()
    for (a,b),weight in symbol.items():
        for key,value in actual_column(a,b).items():
            result[key]+=value*weight
    return {key:value for key,value in result.items() if value}


def fixed_lift(rank):
    """Differential of the previously stored (2,6) primitive on a transverse two-form q wedge lambda."""
    import json
    from pathlib import Path
    rows=json.loads(Path(__file__).with_name('multitrace_preimages.json').read_text())['closed_two_form_preimages']
    primitive=Counter()
    for power,terms in rows.items():
        for left,right,value in terms:
            primitive[tuple(left),tuple(right)]+=rank**int(power)*F(value)
    # Use label 7 for lambda and 8 for q=-sum_{i=0}^6 p_i.
    formal=Counter()
    for j in range(7):
        s=[(j+i)%7 for i in range(7)]
        momenta=[tuple(s[:2])]+[(i,) for i in s[2:]]
        for direction in (-1,1):
            for start in range(6):
                order=[(start+direction*i)%6 for i in range(6)]
                entries={i+2:momenta[index] for i,index in enumerate(order[:5])}
                entries.update({7:(8,),8:(8,),9:(7,)})
                for (left,right),weight in primitive.items():
                    add(formal,tuple(entries[i] for i in left),tuple(entries[i] for i in right),
                        -F(parity(order),12)*weight)
    expanded=Counter()
    for (u,v),value in formal.items():
        if value:
            add(expanded,tuple(tuple(range(7)) if i==8 else (i,) for i in u),
                tuple(tuple(range(7)) if i==8 else (i,) for i in v),
                value*(-1)**(u.count(8)+v.count(8)))
    result=Counter()
    for (a,b),value in expanded.items():
        if not value:continue
        for (_,left),x in coordinates(a).items():
            for (_,right),y in coordinates(b).items():
                variables=tuple(sorted(left+right))
                i=variables[-1][1]
                assert variables[-1][0]==7
                # mu(p,r,s)=det(q,lambda,p,r,s): match the four-vector coordinates in actual().
                rows=tuple(j for j in range(5) if j!=i)
                result[rows,variables[:-1]]+=(-1)**(5-i)*value*x*y
    return {key:value for key,value in result.items() if value}


def second_disk_exchange():
    from primary_disk import solve as disk
    _,_,coefficients=disk(3,3,(-2,-1))
    result=Counter()
    for (external,internal),value in coefficients.items():
        # The gauge pair 0,1 induces a trivector dual to p0 wedge p1.
        # Internal field momentum is p0+p1 and external field momentum is -sum(p).
        det_slots=[(0,),(1,)]+[tuple(range(7)) if i==-2 else (i+2,) for i in internal]
        polyvector_slots=[(0,1) if i==-1 else (i+2,) for i in external]
        add(result,tuple(det_slots),tuple(polyvector_slots),-value if -2 in internal else value)
    return actual(result)


def project(columns,targets):
    pivots={}
    for name,values in columns:
        reduced=Counter(values)
        combination=Counter({name:F(1)})
        for pivot,(stored,coefficients) in pivots.items():
            factor=reduced[pivot]
            if factor:
                for key,value in stored.items():reduced[key]-=factor*value
                for key,value in coefficients.items():combination[key]-=factor*value
        pivot=next((key for key,value in sorted(reduced.items()) if value),None)
        if pivot is not None:
            scale=F(1,reduced[pivot])
            pivots[pivot]=({key:value*scale for key,value in reduced.items() if value},
                           {key:value*scale for key,value in combination.items() if value})
    output=[]
    for target in targets:
        residual=Counter(target)
        rep=Counter()
        for pivot,(column,combination) in pivots.items():
            factor=residual[pivot]
            for key,value in column.items():residual[key]-=factor*value
            for key,value in combination.items():rep[key]+=factor*value
        output.append(({key:value for key,value in rep.items() if value},
                       {key:value for key,value in residual.items() if value}))
    return len(pivots),output


def solve():
    target=actual(wheel(32))
    candidates=[(a,b) for a in combinations(range(6),5) for b in combinations(range(6),3)]
    columns=[actual(cyclic_primitive(a,b)) for a,b in candidates]
    rank,[(rep,residual)]=project(zip(candidates,columns),[target])
    assert not residual and rank==4
    lift=fixed_lift(32)
    assert {key:-5760*value for key,value in lift.items()}==target
    print('single-trace differential rank',rank,'verified monomial coefficients',len(target))
    print('existing (2,6) lift multiplier',-5760)
    exchange=second_disk_exchange()
    for rank in (7,32):  # Rank enters the wheel at degree at most one.
        assert actual(wheel(rank,((0,1),tuple(range(2,7)))))=={key:80*c for key,c in exchange.items()}
        assert not actual(wheel(rank,((0,1,2),tuple(range(3,7)))))
    print('closed disk exchange multiplier',80,'verified monomial coefficients',len(exchange))
    print({key:value for key,value in rep.items() if value})


if __name__=='__main__':
    solve()
