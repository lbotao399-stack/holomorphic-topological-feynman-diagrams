"""Genus-one BCOV vertices at arbitrary descendant levels and conserved rational momenta."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import factorial, prod

from flint import fmpq_mat


@lru_cache(maxsize=4096)
def elliptic_psi(degrees):
    """Top psi intersection on Mbar_1,n, by the string and dilaton equations."""
    if not degrees or sum(degrees) != len(degrees):
        return Fraction(0)
    if degrees[0]:
        return Fraction(factorial(len(degrees)-1),24)
    rest = degrees[1:]
    return sum(elliptic_psi(tuple(sorted((*rest[:i],k-1,*rest[i+1:]))))
               for i,k in enumerate(rest) if k)


def det_square(rows):
    return Fraction(str(fmpq_mat([[str(x) for x in row] for row in rows]).det()))**2


def blocks(marks, count, levels, excess):
    """Keep only partitions with the required number of genus-one blocks."""
    if count == 1:
        if sum(levels[i] for i in marks)==len(marks)-1+excess:
            yield (marks,)
        return
    first,*others = marks
    for size in range(len(marks)-count+1):
        for extra in combinations(others,size):
            block = (first,*extra)
            spent = sum(levels[i] for i in block)-size
            if not 0<=spent<=excess:
                continue
            remaining = tuple(j for j in others if j not in extra)
            for rest in blocks(remaining,count-1,levels,excess-spent):
                yield (block,*rest)


def vertex(momenta,levels):
    if not momenta or not momenta[0] or len(momenta)!=len(levels):
        raise ValueError('Supply a momentum and descendant level for each mark.')
    dimension = len(momenta[0])
    if any(len(p)!=dimension for p in momenta) or any(map(sum,zip(*momenta))):
        raise ValueError('The equal-dimensional momenta must be conserved.')
    if any(type(k) is not int or k<0 for k in levels):
        raise ValueError('Descendant levels must be nonnegative integers.')
    if sum(levels)!=len(levels)-dimension or len(levels)<=dimension:
        return Fraction(0)
    keep=[i for i,(p,k) in enumerate(zip(momenta,levels)) if k!=1 or any(p)]
    if len(keep)<len(levels):
        if len(keep)<=dimension:
            return Fraction(0)
        return (prod(range(len(keep),len(levels))) *
                vertex([momenta[i] for i in keep],[levels[i] for i in keep]))
    total=Fraction()
    for partition in blocks(tuple(range(len(levels))),dimension+1,levels,1):
        ordinary=[]
        weight=Fraction(1)
        for group in partition:
            degrees=tuple(sorted(levels[i] for i in group))
            degree=sum(degrees)
            if degree==len(group)-1:
                ordinary.append(tuple(map(sum,zip(*(momenta[i] for i in group)))))
                weight*=Fraction(factorial(degree),prod(map(factorial,degrees)))
            else:
                weight*=elliptic_psi(degrees)
        total+=weight*det_square(ordinary)
    return total


def elementary_one_loop(scalars,levels,vector,polarization):
    """Return (2 pi i)^5 Q_1 phi, for scalar A_k waves and one transverse primary vector."""
    if len(scalars)!=len(levels) or len(vector)!=5 or len(polarization)!=5:
        raise ValueError('Expected fivefold momenta with one level per scalar.')
    if any(len(p)!=5 for p in scalars) or any(type(k) is not int or not 0<=k<=4 for k in levels):
        raise ValueError('Minimal fivefold scalar levels range from zero to four.')
    if sum(p*v for p,v in zip(vector,polarization)):
        raise ValueError('The primary vector must be transverse.')
    waves=[*scalars,vector]
    root=tuple(-sum(column) for column in zip(*waves))
    return sum(p*v for p,v in zip(root,polarization))*vertex([*waves,root],[*levels,0,0])


def cut_recursion(momenta,levels,fixed=False,last=False):
    """Independent replay by reducing one cotangent line at a time on rational tails."""
    dim,n=len(momenta[0]),len(levels)
    if not any(levels):
        return det_square(momenta[:-1]) if fixed and n==dim+1 else Fraction(0)
    indices=[i for i,k in enumerate(levels) if k]
    i=indices[-1] if last else indices[0]
    lowered=list(levels)
    lowered[i]-=1
    total=Fraction(0) if fixed else cut_recursion(momenta,lowered,True,last)/24
    others=[j for j in range(n) if j!=i]
    for size in range(1,len(others)+1):
        for extra in combinations(others,size):
            group=(i,*extra)
            if sum(lowered[j] for j in group)!=len(group)-2:
                continue
            outside=[j for j in others if j not in extra]
            degree=Fraction(factorial(len(group)-2),prod(factorial(lowered[j]) for j in group))
            fused=tuple(map(sum,zip(*(momenta[j] for j in group))))
            total+=degree*cut_recursion([momenta[j] for j in outside]+[fused],
                                       [levels[j] for j in outside]+[0],fixed,last)
    return total


def check():
    from itertools import combinations_with_replacement
    from random import Random
    rng=Random(440024)
    for dim in (1,2,3,5):
        plane=[tuple(int(i==j) for i in range(dim)) for j in range(dim)]
        simplex=plane+[(-1,)*dim]
        assert vertex(simplex,[1]+[0]*dim)==Fraction(1,24)
        assert (vertex([*simplex,*([(0,)*dim]*64)],[1,*([0]*dim),*([1]*64)])
                ==Fraction(prod(range(dim+1,dim+65)),24))
        for extra in (1,2,3):
            n=dim+extra
            points=[tuple(rng.randrange(-2,3) for _ in range(dim)) for _ in range(n-1)]
            points.append(tuple(-sum(c) for c in zip(*points)))
            for slots in combinations_with_replacement(range(n),extra):
                levels=tuple(slots.count(i) for i in range(n))
                result=vertex(points,levels)
                assert result==cut_recursion(points,levels)==cut_recursion(points,levels,last=True)
    probe=tuple((*row,0) for row in ((1,0),(0,1),(1,2),(0,-1),(2,-3),(-4,1)))
    assert vertex(probe,(0,0,0,1,1,1))==0
    unit=[tuple(int(i==j) for i in range(5)) for j in range(5)]
    scalars=[*unit[2:],(1,1,0,0,0),(-2,-2,-1,-1,-1)]
    assert elementary_one_loop(scalars,(0,0,0,2,0),unit[1],unit[0])==Fraction(1,2)
    assert elementary_one_loop(scalars,(0,0,0,1,1),unit[1],unit[0])==Fraction(7,12)
    background=[*unit[2:],(-1,)*5,(0,)*5]
    for m in range(9):
        assert (elementary_one_loop([*background,*([(0,)*5]*m)],
                                    [0,0,0,0,2,*([1]*m)],unit[1],unit[0])
                == Fraction(prod(range(7,7+m)),4))
    print('Partition sum and both rational-tail recursions agree in dimensions 1, 2, 3, 5.')


if __name__=='__main__':
    check()
