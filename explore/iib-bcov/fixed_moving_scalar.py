"""One fixed flat IIB scalar and vector profile measures the threefold period at every genus."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json

from collinear_abel_moments import genus_three
from elementary_scalar import elementary_three_loop
from elliptic_descendants import vertex as elliptic
from genus_one_graph_ward import plus
from two_loop_moving_descendant import BASIS, descendant_vertices


def coefficient(genus, a, b, simplex, fourfold):
    """Return (-1)^g (2 pi i)^(5g) Q_g in the period basis of the supplied vertices."""
    return tuple((a*b)**genus*d + (a*b)**(genus-1)*(a+b)**2*w/factorial(genus-1)
                 for d,w in zip(simplex,fourfold))


def smaller_genus(genus, a, b):
    total = F()
    for m in range(2*genus+1):
        l = 2*genus-m
        regular = (*BASIS[:3], *(BASIS[3],)*m, *(BASIS[4],)*l)
        waves = (*regular,tuple(-x for x in plus(regular)))
        value = (elliptic(waves,[int(i==2) for i in range(len(waves))])
                 if genus==1 else descendant_vertices(waves)[2])
        assert bool(value) == (m==1 if genus==1 else abs(m-genus)<=1)
        total += value*a**m*b**l/(factorial(m)*factorial(l))
    return total


def positive_primary_profile(genus):
    """Sum all three primary scalar multiplicities, with one vector and one A_1 wave."""
    if genus not in (1, 2):
        raise ValueError('The independent moving-descendant graph banks cover genera one and two.')
    total = F()
    for first in range(1, 2*genus+2):
        for fourth in range(2*genus+2-first):
            fifth = 2*genus+1-first-fourth
            regular = (BASIS[1], BASIS[2], *(BASIS[0],)*first,
                       *(BASIS[3],)*fourth, *(BASIS[4],)*fifth)
            waves = (*regular,tuple(-x for x in plus(regular)))
            value = (elliptic(waves,[int(i==1) for i in range(len(waves))])
                     if genus==1 else descendant_vertices(waves)[1])
            assert value >= 0
            assert bool(value) == bool(fourth and fifth)
            total += first*value/(factorial(first)*factorial(fourth)*factorial(fifth))
    return total


def third_genus():
    # Compute every possible one-tail primary using the independent nine-field trace.
    expected = genus_three(1,1)
    totals = []
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        source = Path(__file__).resolve().parent/'operator'/'nine_sextics.cpp'
        run(['c++','-std=c++17','-O3',str(source),'-o',str(exe)],check=True)
        for m in range(7):
            l = 6-m
            value = [F(),F()]
            for axis,count in ((3,m),(4,l)):
                if not count:
                    continue
                fused = plus((BASIS[2],BASIS[axis]))
                scalars = [BASIS[0],fused]+[BASIS[3]]*(m-int(axis==3))+[BASIS[4]]*(l-int(axis==4))
                # The heat-stripped primary trace is -p_root(v)*V, with p_root(v)=-1.
                row = elementary_three_loop(scalars,BASIS[1],BASIS[0],exe=exe)
                value = [x+count*y for x,y in zip(value,row)]
            multiplicity = {2:24,3:36,4:24}.get(m,0)
            assert value==[multiplicity*w for w in expected]
            totals.append(tuple(value))
    assert all(sum((-1)**(6-m)*row[i]/(factorial(m)*factorial(6-m))
                   for m,row in enumerate(totals))==0 for i in (0,1))
    print('All seven genus-three multiplicities and both primary period columns passed.')


def check():
    for genus,simplex,fourfold in ((1,F(1,24),F()),(2,F(1,9720),F(1,62208))):
        for a,b in ((1,-1),(1,1),(2,-3),(F(2,3),F(3,5)),(0,1)):
            assert smaller_genus(genus,a,b)==coefficient(genus,a,b,(simplex,),(fourfold,))[0]
        assert positive_primary_profile(genus)==(simplex if genus==1 else 5*simplex+15*fourfold)
    assert 5*F(1,9720)+15*F(1,62208)==F(47,62208)
    third_genus()
    bank = Path(__file__).resolve().parents[1]/'kodaira-spencer'
    d4 = json.loads((bank/'four_loop_descendant.json').read_text())
    w4 = json.loads((bank/'four_loop_arakelov_periods.json').read_text())
    rows = [row for group in d4['graphs'].values() for row in group]
    assert len(rows)==w4['checked_graph_count']==50
    simplex = tuple(sum(map(F,col)) for col in zip(*rows))
    assert simplex==tuple(map(F,d4['sum']))==coefficient(4,1,-1,simplex,tuple(map(F,w4['equal_strength_W'])))
    print('Fixed-profile heat coefficients at genera one and two: 1/24, 1/9720.')
    print('Genus four, in the 1, zeta(3), zeta(5) basis:',*simplex)


if __name__=='__main__':
    check()
