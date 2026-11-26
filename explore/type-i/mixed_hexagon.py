"""Exact geometric factor for each derivative assignment on the one-closed-edge wheel."""
from fractions import Fraction as F
from itertools import permutations
from math import factorial
from random import Random
from hexagon import det


def current_insertions(momenta,closed_edge,slots,t=None):
    """slots contain four oriented entries ('external',covector) or ('edge',(number,sign))."""
    if (len(momenta)!=6 or any(len(p)!=5 for p in momenta) or any(map(sum,zip(*momenta)))
            or closed_edge not in range(6) or len(slots)!=4):
        raise ValueError('expected six conserved five-component vertex momenta, one cut, and four disk insertions')
    t=[1]*6 if t is None else t
    partial=[[sum(momenta[s][k] for s in range(e+1)) for k in range(5)] for e in range(6)]
    mean=[sum(F(w,sum(t))*p[k] for w,p in zip(t,partial)) for k in range(5)]
    currents=[tuple(p[k]-mean[k] for k in range(5)) for p in partial]
    rows=[currents[closed_edge]]
    for kind,v in slots:
        if kind=='external':
            rows.append(v)
        elif kind=='edge':
            e,sgn=v
            if sgn not in (-1,1) or e not in range(6):
                raise ValueError('an edge insertion needs an edge number and endpoint sign')
            rows.append(tuple(sgn*a for a in currents[e]))
        else:
            raise ValueError('unknown derivative insertion')
    return det(rows)


def geometric(momenta,closed_edge,slots):
    return -F(1,factorial(5))*det(momenta[:5])*current_insertions(momenta,closed_edge,slots)


def check():
    rng=Random(25415)
    seen=set()
    for closed in range(6):
        p=[[rng.randrange(-4,5) for _ in range(5)] for _ in range(5)]
        p.append([-sum(c) for c in zip(*p)])
        for incidence in range(4):
            slots=[('edge',((closed-1)%6,-1)) if incidence&1 else ('external',p[0]),
                   ('edge',((closed+1)%6,1)) if incidence&2 else ('external',p[1]),
                   ('external',p[2]),('external',p[3])]
            corners=[current_insertions(p,closed,slots,[int(i==j) for i in range(6)]) for j in range(6)]
            mean=sum(corners,F(0))/6
            assert geometric(p,closed,slots)==-det(p[:5])*mean/120
            t=[rng.randrange(1,10) for _ in range(6)]
            assert current_insertions(p,closed,slots,t)==sum(F(w,sum(t))*a for w,a in zip(t,corners))
            if mean:seen.add(incidence)
    assert seen==set(range(4))
    print('all six closed-edge positions; all four internal/external derivative cases; affine Dirichlet integral')


if __name__=='__main__':
    check()


def single_trace():
    """Polarized disk tensor with the displayed determinant and iterated commutator."""
    from collections import Counter
    from open_theta import relabelings,reduce_minor,straight
    words=Counter({(0,):1})
    for leg in range(1,5):
        commutator=Counter()
        for word,c in words.items():
            commutator[(leg,)+word]+=c
            commutator[word+(leg,)]-=c
        words=commutator
    result=Counter()
    for word,c in words.items():
        for perm,sign in relabelings((word+(5,6,7),)):
            first=reduce_minor(tuple(perm[k] for k in range(5)))
            for j,weight in ((1,3),(2,2),(3,1)):
                second=reduce_minor(tuple(perm[k] for k in (0,5,6,7,j)))
                for a,x in first:
                    for b,y in second:
                        for monomial,sign2 in straight(a,b):
                            result[monomial]+=c*sign*weight*x*y*sign2
    return {k:v for k,v in result.items() if v}


def disk_signs():
    from hexagon import parity
    rng=Random(961440)
    def three_form(q,d1,d3):
        return sum(parity(s)*det([q,d1,*[d3[j] for j in s[1:]]]) for s in permutations(range(4)))
    # Color inputs of the rooted tree are split (two,four). Its spatial heat edge points into the integrated vertex.
    for _ in range(4):
        p=[[rng.randrange(-3,4) for _ in range(5)] for _ in range(5)]
        p.append([-sum(x) for x in zip(*p)])
        q=[-a-b for a,b in zip(*p[:2])]
        d1=[a-b for a,b in zip(*p[:2])]
        assert three_form(q,d1,p[2:])==-48*det(p[:5])
    # In the mixed hexagon the last three waves occur at the same PV³ disk.
    for _ in range(4):
        p=[[rng.randrange(-3,4) for _ in range(5)] for _ in range(7)]
        p.append([-sum(x) for x in zip(*p)])
        partial=[tuple(sum(p[j][i] for j in range(k+1)) for i in range(5)) for k in range(5)]+[(0,)*5]
        q=[-sum(F(row[i],6) for row in partial) for i in range(5)]
        d1=[2*a+b for a,b in zip(p[0],q)]
        y4=[a+b for a,b in zip(partial[4],q)]
        last=[3*a+2*b+c for a,b,c in zip(*p[1:4])]
        assert three_form(q,d1,[*p[5:],y4])==8*det([p[0],*p[5:],last])
    # The open slots at the disks are ordered (external,internal), (three external,internal).
    global_slots=[('x',i) for i in range(8)]+[(e,s) for e in range(5) for s in range(2)]
    local=[('x',0),(0,0)]
    for i in range(1,5):
        local += [('x',i),(i-1,1),(i,0)]
    local += [('x',5),('x',6),('x',7),(4,1)]
    assert parity([global_slots.index(s) for s in local])==-1
    tree_heat=-1
    assert tree_heat*(-48)==48
    assert -F(8,96*factorial(3)*120)==-F(1,8640)
    print('rooted tree sign and four cubic Wick signs; mixed subtraction -1/8640')


if __name__=='__main__':
    disk_signs()
