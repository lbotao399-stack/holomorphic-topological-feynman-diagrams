"""Planar primary loops on z²/2 + y²(c+exp x)/2 along its punctured critical line."""
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
from math import factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from planar_cubic_primary import HARD, ROOT, UNIT, labeled_vertex, labeled_planar_rows
from worldsheet_primary import fourfold_primary


def close_fourfold(points):
    points = [p[:4] for p in HARD]+[(0,0,*p) for p in points]
    return points+[tuple(-sum(col) for col in zip(*points))]


def critical_line(genus,exe):
    coefficients = []
    for count in range(genus+1):
        fixtures,weights = [],[]
        for slopes in product((-1,0,1),repeat=genus):
            normals = tuple((int(i<count),a) for i,a in enumerate(slopes))
            fixtures.append(normals)
            weights.append(prod(F(1,2) if a else F(-1) for a in slopes)/
                           (factorial(count)*factorial(genus-count)))
        five = [tuple((*v,0) for v in normal)+(UNIT[2],)*genus for normal in fixtures]
        values = labeled_vertex(genus,five,exe)
        coefficient=[F(0),F(0)]
        for pts,w,full in zip(fixtures,weights,values):
            lower=fourfold_primary(genus,close_fourfold(pts))
            assert tuple(a/factorial(genus) for a in full)==lower
            for j,column in enumerate(lower):
                coefficient[j]+=w*column
        coefficients.append(tuple(coefficient))
    assert coefficients[:2]==[(0,0),(0,0)]
    return coefficients


def soft_simplex(genus, direction=(-1,1)):
    # At two exponential slots the center-of-mass slope is a normal shear.
    # Interpolate the two-loop constant or the coefficient of t² at genus 3.
    import sympy as sp
    t = sp.Symbol("t")
    degree = 2*genus if genus>2 else 0
    if len(direction)==2:
        direction=(0,0,*direction)
    simplex=[tuple(int(i==j) for i in range(4)) for j in range(4)]
    values=[]
    for k in range(degree+1):
        points=simplex+[tuple(k*v for v in direction)]*(genus-2)
        values.append(fourfold_primary(genus,[*points,tuple(-sum(c) for c in zip(*points))]))
    return tuple(F(sp.interpolate(list(enumerate(col)),t).expand().coeff(t,2*(genus-2)))
                 for col in zip(*values))


def frequency_jet(frequencies):
    """Reduced labeled symbol after extracting y²/2 at each of g slots."""
    genus=len(frequencies)
    coefficients=[F(0),F(0)]
    for slopes in product((-1,0,1),repeat=genus):
        weight=prod(F(1,2) if a else F(-1) for a in slopes)
        vertex=fourfold_primary(genus,close_fourfold(zip(frequencies,slopes)))
        for i,x in enumerate(vertex):
            coefficients[i]+=weight*x
    return tuple(coefficients)


def polarized_pair(genus):
    for p,q in ((1,-1),(2,3),(-1,-1)):
        measured=frequency_jet((p,q)+(0,)*(genus-2))
        target=soft_simplex(genus,(-q,p))
        assert measured == tuple(6*(p*q)**2*x for x in target)
        if (p,q)==(1,-1):
            print('opposite exponential coefficient',genus,
                  tuple(x/factorial(genus-2) for x in measured),flush=True)


def tangent_deformations(genus,exe):
    for p,q,L,M in ((2,-3,(1,-1),(2,3)),(1,1,(0,0),(-2,-3))):
        rows,weights=[],[]
        for slopes in product((-1,0,1),repeat=genus):
            slots=[(*L,p,slopes[0],0),(*M,q,slopes[1],0)]
            soft=[(0,0,0,t,0) for t in slopes[2:]]
            rows.append([*slots,*soft,*([(0,0,0,0,1)]*genus),*HARD])
            weights.append(prod(F(1,2) if s else F(-1) for s in slopes))
        sums=tuple(sum(weight*x for weight,x in zip(weights,col))/
                   (factorial(genus)*factorial(genus-2))
                   for col in zip(*labeled_planar_rows(genus,rows,exe)))
        direction=(q*L[0]-p*M[0],q*L[1]-p*M[1],-q,p)
        assert len(set((*direction,0)))==5
        factor=F(6*(p*q)**2*(1+L[0]+M[0]),factorial(genus-2))
        assert sums==tuple(factor*x for x in soft_simplex(genus,direction))
        print('fully distinct simplex coordinates',genus,direction,sums,flush=True)


def cubic_jet_polynomial():
    """Interpolate the degree-six symmetric polynomial on the three x frequencies."""
    import sympy as sp
    partitions=[k for k in combinations_with_replacement(range(7),3) if sum(k)==6]
    words=[tuple(set(permutations(k))) for k in partitions]
    def row(points):
        return [sum(prod(x**e for x,e in zip(points,word)) for word in group)
                for group in words]
    points=list(combinations_with_replacement(range(5),3))
    pivot_indices=sp.Matrix([row(p) for p in points]).T.rref()[1]
    nodes=[points[k] for k in pivot_indices]
    assert len(nodes)==len(partitions)
    coefficients=sp.Matrix([row(p) for p in nodes]).inv()*sp.Matrix([frequency_jet(p) for p in nodes])
    for p in ((-1,2,3),(2,-3,4),(5,3,-1)):
        assert tuple(sp.Matrix([row(p)])*coefficients) == frequency_jet(p)
    # Equal background fields give 1/3!, with one copy of each distinct permutation.
    return {p:tuple(F(x)*len(group)/6 for x in coefficients.row(i))
            for i,(p,group) in enumerate(zip(partitions,words)) if any(coefficients.row(i))}


def normal_residues(coefficients):
    import sympy as sp
    x,y=sp.symbols('x y')
    profile=(1-x**2)**2
    root=1-x**2
    two=profile.diff(x,2)**2/20736
    three=[sum(value[i]*prod(profile.diff(x,j) for j in partition)
               for partition,value in coefficients.items()) for i in range(2)]
    for point,sign in ((-1,1),(1,-1)):
        assert sp.residue(two/root,x,point)==F(sign,648)
        for i,term in enumerate(three):
            assert sp.residue(term/root,x,point)==sign*256*coefficients[2,2,2][i]
    # Mix the normal coordinates after a nonlinear line reparametrization.
    f=x+x**2
    Y=y/f.diff(x)
    X=f+Y
    assert sp.cancel(sp.Matrix([X,Y]).jacobian((x,y)).det()) == 1
    assert sp.cancel((Y**2*X**2).diff(y,2).subs(y,0)/2-f**2/f.diff(x)**2)==0
    double_profile=sp.cancel(f**2/f.diff(x)**2)
    assert double_profile.diff(x,2).subs(x,0)==2


if __name__ == "__main__":
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        for g in (2,3):
            columns = critical_line(g,exe)
            assert columns[2] == tuple(3*x/factorial(g-2) for x in soft_simplex(g))
            print('genus',g,'polynomial coefficients in increasing exponential multiplicity',columns,flush=True)
            polarized_pair(g)
            tangent_deformations(g,exe)
        coefficients=cubic_jet_polynomial()
        normal_residues(coefficients)
        print('three-loop differential polynomial',coefficients,flush=True)
