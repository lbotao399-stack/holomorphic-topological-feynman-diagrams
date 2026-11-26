"""Contraction at determinant weight -3 on block-diagonal S_2 + S_3 vector fields."""
from collections import Counter
from itertools import combinations, product
from math import comb, factorial
from fractions import Fraction


def decompositions(total,minimum=1):
    if total==0:
        yield ()
    for part in range(minimum,total+1):
        for rest in decompositions(total-part,part):
            yield (part,*rest)


def planar_dimensions(grade,output_weight):
    # R_h in two dimensions is the binary-form representation Sym^(h+2).
    by_arity=Counter()
    for profile in decompositions(grade):
        weights=Counter({0:1})
        for h,multiplicity in Counter(profile).items():
            individual=range(-h-2,h+3,2)
            exterior=Counter(map(sum,combinations(individual,multiplicity)))
            nxt=Counter()
            for a,x in weights.items():
                for b,y in exterior.items():nxt[a+b]+=x*y
            weights=nxt
        by_arity[len(profile)]+=weights[output_weight]-weights[output_weight+2]
    return {k:v for k,v in sorted(by_arity.items()) if v}


def falling(n,r):
    return factorial(n)//factorial(n-r) if n>=r else 0


def transvectant(f,g,order):
    result=Counter()
    for (a,b),x in f.items():
        for (c,d),y in g.items():
            coefficient=sum((-1)**i*comb(order,i)*falling(a,order-i)*falling(b,i)
                            *falling(c,i)*falling(d,order-i) for i in range(order+1))
            if coefficient:result[a+c-order,b+d-order]+=coefficient*x*y
    return {key:value for key,value in result.items() if value}


def parity(order):
    return (-1)**sum(i>j for pos,i in enumerate(order) for j in order[pos+1:])


def cochain(name,inputs):
    if any(not x for x in inputs):return {}
    degrees=[sum(next(iter(x))) for x in inputs]
    order=sorted(range(len(inputs)),key=degrees.__getitem__)
    profile=tuple(degrees[i] for i in order)
    args=[inputs[i] for i in order]
    if name=='S' and profile==(5,5):out=transvectant(*args,5)
    elif name=='P' and profile==(4,5):out=transvectant(*args,4)
    elif name=='A' and profile==(3,3,5):out=transvectant(transvectant(*args[:2],1),args[2],4)
    elif name=='B' and profile==(3,4,4):out=transvectant(args[0],transvectant(*args[1:],3),2)
    else:return {}
    return {key:parity(order)*value for key,value in out.items()}


def differential(name,args):
    out=Counter()
    for i,j in combinations(range(len(args)),2):
        pair=transvectant(args[i],args[j],1)
        rest=[v for h,v in enumerate(args) if h not in (i,j)]
        for key,val in cochain(name,[pair,*rest]).items():out[key]+=(-1)**(i+j)*val
    return {key:value for key,value in out.items() if value}


def monomials(degree):
    return [{(i,degree-i):1} for i in range(degree+1)]


def main():
    assert planar_dimensions(6,0)=={2:1,3:1}
    assert planar_dimensions(5,1)=={2:1,3:2,4:1}
    assert differential('S',[{(3,0):1},{(0,4):1},{(3,2):1}])=={(0,0):17280}
    coefficient=None
    for profile in [(3,3,5),(3,4,4)]:
        for args in product(*(monomials(n) for n in profile)):
            raw=Counter(differential('P',args))
            raw.update(cochain('A',args))
            residual={k:v for k,v in raw.items() if v}
            other=cochain('B',args)
            if other:
                key=next(iter(other))
                ratio=Fraction(residual[key],other[key])
                if coefficient is None:coefficient=ratio
                assert coefficient==ratio
            expected={k:coefficient*v for k,v in other.items()} if other else {}
            assert residual==expected
    witness=None
    for args in product(monomials(3),monomials(3),monomials(3),monomials(4)):
        a=differential('A',args);b=differential('B',args)
        assert a=={key:coefficient*val for key,val in b.items()}
        if b and witness is None:witness=([next(iter(f)) for f in args],b)
    assert coefficient==7
    assert witness==([(0,3),(1,2),(2,1),(3,1)],{(0,1):-2592})
    print('planar scalar dimensions',planar_dimensions(6,0),'scalar differential probe',17280)
    print('planar current dimensions',planar_dimensions(5,1),'dP=-A+cB, c=',coefficient)
    print('dB witness',witness)


if __name__=='__main__':
    main()
