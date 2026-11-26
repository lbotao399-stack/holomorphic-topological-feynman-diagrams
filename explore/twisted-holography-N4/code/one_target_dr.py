"""BSSZ (arXiv:1211.5273, Theorems 1-2) and the one-mark two-Betti intersection."""
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, prod
from pathlib import Path
import sys

from flint import fmpq, fmpq_mpoly_ctx


def exact(x):
    return fmpq(str(x))


def homogeneous_dr(genus, momenta):
    """Return the symmetric descendant polynomial at degree 2g in conserved weights.

    The scalar graph deletion and contraction identity identifies its coefficients
    with the one-target BCOV descendant vertices at every positive genus.
    """
    n = len(momenta)
    if n<1 or sum(momenta) or genus<1:
        raise ValueError('Supply at least one conserved momentum, at positive genus.')
    ctx = fmpq_mpoly_ctx.get(tuple('z'+str(i) for i in range(n)))
    z = ctx.gens()
    zero = ctx.constant(0)
    if not any(momenta):
        return zero
    if n==2:
        return (sum(z)**(2*genus-1)*exact(momenta[0])**(2*genus)
                / (4**genus*factorial(2*genus+1)))
    if sum(a==0 for a in momenta)>1:
        # A pair of zero weights makes an identically vanishing denominator.
        # Specialize after cancellation by polynomial interpolation in a conserved perturbation.
        direction = (*range(1,n),-n*(n-1)//2)
        points = []
        t = 0
        while len(points)<=2*genus:
            t += 1
            shifted = tuple(a+t*b for a,b in zip(momenta,direction))
            if all(shifted):
                points.append((t,homogeneous_dr(genus,shifted)))
        return sum((poly*exact(prod(Fraction(-s,t-s) for s,_ in points if s!=t))
                    for t,poly in points),zero)
    a = tuple(map(exact,momenta))
    pairs = {(i,j):a[i]*z[j]-a[j]*z[i] for i,j in combinations(range(1,n),2)}
    numerator = zero
    for tail in permutations(range(1,n)):
        order = (0,*tail)
        prefix_a, prefix_z = a[0],z[0]
        cross = []
        for j in tail:
            cross.append(prefix_a*z[j]-a[j]*prefix_z)
            prefix_a, prefix_z = prefix_a+a[j], prefix_z+z[j]
        # First determinant cancels its own denominator. The last determinant is -a_last * sum(z).
        factor = -a[order[-1]]*prod(z[i] for i in order[1:-1])*prod(cross[1:-1])
        arcs = {tuple(sorted((i,j))) for i,j in zip(tail,tail[1:])}
        factor *= (-1)**sum(i>j for i,j in zip(tail,tail[1:]))
        factor *= prod(val for edge,val in pairs.items() if edge not in arcs)
        coeff = [ctx.constant(1)]+[zero]*genus
        for value in cross:
            series = [value**(2*k)/(4**k*factorial(2*k+1)) for k in range(genus+1)]
            coeff = [sum((coeff[i-k]*series[k] for k in range(i+1)),zero) for i in range(genus+1)]
        numerator += factor*coeff[genus]
    return numerator / prod(pairs.values())


def degrees(total, count):
    if count==1:
        yield (total,)
        return
    for first in range(total+1):
        for tail in degrees(total-first,count-1):
            yield (first,*tail)


def single_mark_pair(genus, a, b, mark):
    """Two actual Betti forms times psi_mark^(3g-5+n), on every stable case with g>=1."""
    if genus<1 or len(a)!=len(b) or sum(a) or sum(b) or not 0<=mark<len(a):
        raise ValueError('Supply conserved charge lists, a valid mark, and positive genus.')
    if 3*genus-5+len(a)<0:
        raise ValueError('Negative descendant degree.')
    labels = [i for i in range(len(a)) if i!=mark]
    value = exact(Fraction(3*genus*(genus-1),5))*sum(a[i]**2*b[i]**2 for i in labels)
    for i,j in combinations(labels,2):
        value += genus**2*(a[i]**2*b[j]**2+a[j]**2*b[i]**2)-2*genus*a[i]*a[j]*b[i]*b[j]
    return value/(24**genus*factorial(genus))


def check_single_mark():
    """The remaining coefficients have degree <=2 in g after division by <tau_(3g-2)>.

    The D D term has at most total cotangent exponent two away from the marked slot.
    String and dilaton leave only <tau_2 tau_(3g-3)> / <tau_(3g-2)>,
    equal to (12g^2-12g+5)/5 by the first Virasoro recursion. Every non-root
    component of a separating tree has at least two flags and at most one cotangent
    power in total, so has genus zero. On the two-parallel-edge graph the root has
    genus g-1 with all other powers zero, giving <tau_(3g-5)>/<tau_(3g-2)>=24g.
    This proves the degree bound for g>=2. Unused markings lower the marked power
    by the string equation. Since each pair of edges uses at most four other wave
    labels, four independent weights on each polarization test every monomial.
    """
    from soft_quartic_two_betti import four_times_betti_pair
    ctx = fmpq_mpoly_ctx.get(tuple('u'+str(i) for i in range(4))+tuple('v'+str(i) for i in range(4)))
    variables = ctx.gens()
    a = (-sum(variables[:4]),*variables[:4])
    b = (-sum(variables[4:]),*variables[4:])
    pairs = tuple(combinations(range(5),2))
    for g in (1,2,3,4,5):
        # g=2,3,4 determines the quadratic uniquely. g=1 is checked separately.
        degrees = (3*g,0,0,0,0)
        table = {(p,q):four_times_betti_pair(g,degrees,p,q,False)/4 for p in pairs for q in pairs}
        raw = sum(a[i]*a[j]*b[k]*b[l]*exact(val) for ((i,j),(k,l)),val in table.items())
        assert raw==single_mark_pair(g,a,b,0)
    print('Independent scalar weights on four labels verify every incidence type at genera 1-5.')


def check():
    from soft_quartic_two_betti import descendant_vertex
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'iib-bcov'))
    from elliptic_descendants import vertex
    for weights in ((1,-1),(2,-3,1),(1,2,-4,1),(1,2,0,-3),(1,0,0,-1),(2,-1,3,-5,1)):
        n = len(weights)
        p = tuple((a,) for a in weights)
        for g in (1,2):
            table = homogeneous_dr(g,weights).to_dict()
            signatures = list(degrees(n+2*g-3,n))
            # Full three- and four-point checks, then several nonsymmetric five-point rows.
            if n==5:
                signatures = signatures[::13]+signatures[-1:]
            for k in signatures:
                expected = vertex(p,k) if g==1 else descendant_vertex(g,p,k)
                assert Fraction(str(table.get(k,0)))==expected, (weights,k,g)
            if g==2:
                for i in range(n):
                    rest = [a for j,a in enumerate(weights) if j!=i]
                    single = sum(Fraction(a**4,1920) for a in rest)
                    single += sum(Fraction(a*a*b*b,576) for a,b in combinations(rest,2))
                    assert Fraction(str(table.get(tuple((n+1)*int(i==j) for j in range(n)),0)))==single
        print(weights,'independent genus-one and genus-two intersection tests passed.',flush=True)
    check_single_mark()
    weights = (2,-3,1)
    table = homogeneous_dr(2,weights).to_dict()
    for k in ((4,0,0),(3,1,0),(2,1,1),(0,2,2)):
        print(weights,k,table[k])


if __name__=='__main__':
    check()
