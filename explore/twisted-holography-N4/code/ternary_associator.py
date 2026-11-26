"""An explicit quantum primitive for the SU(2) cyclic three-graviton associator.

The coefficient units in the marked product are 1, i/L, 1/L^2, i/L^3,
where L=32*pi^2. In its primitive the units are 1/L, -i/L^2, 1/L^3.
All entries omit the common three fundamental trace factors, 1/8.
The zero and top components in the superfield are ig*c and b/(ig).
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
import ternary_product as p
from ternary_product_cubic import scalar


def invariant(letters, colors, edges):
    """Each edge contracts a pair of derivative slots with the spinor epsilon."""
    out=defaultdict(int)
    for indices,color_factor in colors:
        for spins in product(p.SPIN,repeat=len(edges)):
            jets=[[] for _ in letters]
            coefficient=color_factor
            for (i,j),(a,b,eps) in zip(edges,spins):
                jets[i].append(a)
                jets[j].append(b)
                coefficient*=eps
            if any(kind=='ghost' and not jet for (kind,_),jet in zip(letters,jets)):
                continue
            word=tuple(p.letter(kind,flavor,color,jet)
                       for (kind,flavor),color,jet in zip(letters,indices,jets))
            p.add(out,word,coefficient)
    return {word:value for word,value in out.items() if value}


def quartic_candidates():
    for i,j in combinations(range(3),2):
        letters=(('beta',i),('beta',j),('gamma',i),('gamma',j))
        for mate in (1,2,3):
            colors=[(tuple(c if n in (0,mate) else d for n in range(4)),1)
                    for c,d in product(range(3),repeat=2)]
            for edge in combinations(range(4),2):
                yield (i,j,mate,edge),invariant(letters,colors,[edge])


def cubic_candidates():
    colors=[((a,b,c),sign) for a,b,c,sign in p.COLORS]
    patterns=[(('gamma',0),('gamma',1),('gamma',2)),(('top',0),('ghost',0),('ghost',0))]
    patterns += [(('ghost',0),('beta',f),('gamma',f)) for f in range(3)]
    for letters in patterns:
        for edges in combinations_with_replacement(combinations(range(3),2),2):
            yield (letters,edges),invariant(letters,colors,edges)


def reduce_column(vector,coefficients,basis):
    while vector:
        pivot=min(vector)
        if pivot not in basis:
            break
        image,source=basis[pivot]
        factor=-vector[pivot]
        vector=p.plus((vector,1),(image,factor))
        coefficients=p.plus((coefficients,1),(source,factor))
    return vector,coefficients


def preimage(target,candidates):
    basis={}
    sources={}
    for label,source in candidates:
        sources[label]=source
        image,coefficients=reduce_column(p.classical(source),{label:1},basis)
        if image:
            pivot=min(image)
            scale=Fraction(1,image[pivot])
            basis[pivot]=p.plus((image,scale)),p.plus((coefficients,scale))
    remainder,coefficients=reduce_column(target,{},basis)
    assert not remainder
    solution={label:-value for label,value in coefficients.items()}
    primitive=p.plus(*[(sources[label],value) for label,value in solution.items()])
    assert not p.plus((p.classical(primitive),1),(target,-1))
    return primitive,solution


def associator():
    left,right=p.first_order()
    difference=p.plus((left,1),(right,-1))
    one,first_coefficients=preimage(difference,quartic_candidates())

    left,right=p.second_order()
    difference=p.plus((left,1),(right,-1),(p.first_fan_displaced(one),Fraction(1,2)))
    two,second_coefficients=preimage(difference,cubic_candidates())

    cubic,quadratic=scalar()
    difference=p.plus((cubic,Fraction(1,24)-Fraction(1,72)),
                      (p.first_fan_displaced(two),-Fraction(1,2)))
    last_coefficient=Fraction(1,72)
    three=p.plus((quadratic,last_coefficient))
    assert not p.plus((p.classical(three),1),(difference,-1))
    assert not p.first_fan_displaced(three)
    return (one,two,three),(first_coefficients,second_coefficients,last_coefficient)


if __name__=='__main__':
    primitives,coefficients=associator()
    print('Component monomials in the three primitive coefficients:',*[len(x) for x in primitives])
    print('First-order invariant contractions:',coefficients[0])
    print('Second-order invariant contractions:',coefficients[1])
    print('Coefficient of contracted third ghost jets in the primitive:',coefficients[2])
