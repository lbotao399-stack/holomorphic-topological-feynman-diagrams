"""Shuffle the odd fibre contractions past de Rham kernels, then count the root orbits."""
from fractions import Fraction
from itertools import permutations
from math import factorial, comb


def par(a):
    return (-1)**sum(x>y for j,x in enumerate(a) for y in a[j+1:])


for cut in range(13):
    # Target order is all form factors left, then all contraction operators.
    # P forms and contraction operators are odd; the heat form is even.
    ordered=[]
    for e in range(13):
        if cut!=e:
            ordered.append(e)
        ordered.append(13+e)
    assert par(ordered)==(-1)**cut
    assert (comb(13,2)-cut)%2==cut%2
print("All 13 cut positions reproduce the single boundary sign.")


def automorphisms(lengths):
    vertices=["N","R"]+[(path,i) for path,n in enumerate(lengths) for i in range(1,n)]
    maps=[]
    for pi in permutations(range(3)):
        if any(lengths[i]!=lengths[j] for i,j in enumerate(pi)):
            continue
        for reverse in (False,True):
            h={"N":"R" if reverse else "N","R":"N" if reverse else "R"}
            for path,n in enumerate(lengths):
                for i in range(1,n):h[path,i]=(pi[path],n-i if reverse else i)
            maps.append(h)
    return vertices,maps

for shape in ((1,6,6),(2,5,6),(3,4,6),(3,5,5),(4,4,5)):
    vertices,group=automorphisms(shape)
    orbits={frozenset(h[v] for h in group) for v in vertices}
    markedsum=0
    for orbit in orbits:
        v=next(iter(orbit))
        stabilizer=sum(h[v]==v for h in group)
        assert len(orbit)*stabilizer==len(group)
        markedsum+=factorial(11)//stabilizer
    raw=factorial(12)//len(group)
    assert markedsum==raw and Fraction(raw,factorial(11))==Fraction(12,len(group))
    print(shape, 'root orbit multiplicities',sorted(map(len,orbits)), 'raw labelings',raw)

bulk = Fraction(-73,8648640)*Fraction(-1,4096)
assert bulk == Fraction(73,35424829440)
assert 12*bulk == Fraction(73,2952069120)
assert factorial(12)*bulk == Fraction(3285,3328)
print("In units of 1/pi^10 for the coefficient of g: 1/12! cut",bulk,
      "Q2 I",12*bulk,"raw twelve-bracket",factorial(12)*bulk)

