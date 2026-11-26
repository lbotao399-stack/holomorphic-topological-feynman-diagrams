"""Five-edge diamond scalar charge comparison by two freely summed residues."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, prod

from betti_flow_polytope import charge_cuts, normal_coefficient


EDGES=((0,3),(0,1),(1,3),(0,2),(2,3))


def residue_sum(modulus,charges,powers):
    # Net outgoing convention. The remaining edge weights are fixed by the
    # direct-edge source and the two independent currents into the middle sites.
    source,left,right,_=charges
    f=[j*(modulus-j) for j in range(modulus)]
    upper=[c**powers[1]*f[(j+left)%modulus]**powers[2] for j,c in enumerate(f)]
    lower=[c**powers[3]*f[(j+right)%modulus]**powers[4] for j,c in enumerate(f)]
    direct=[v**powers[0] for v in f]
    return sum(up*sum(down*direct[(source-u-v)%modulus] for v,down in enumerate(lower))
               for u,up in enumerate(upper))


def modulus_constant(charges,powers):
    degree=2*sum(powers)+2
    start=1+sum(q for q in charges if q>0)
    nodes=range(start,start+degree+1)
    values=[residue_sum(r,charges,powers) for r in nodes]
    # Coefficient of r^2 of the unnormalized sum by Lagrange interpolation.
    result=F()
    for i,r in enumerate(nodes):
        other=[s for s in nodes if s!=r]
        second=sum(F(1,a*b) for a,b in combinations(other,2))
        result+=F(prod(other)*values[i],(-1)**(degree-i)*factorial(i)*factorial(degree-i))*second
    assert residue_sum(nodes[-1]+1,charges,powers)==sum(
        (-1)**(degree-i)*comb(degree+1,i)*y for i,y in enumerate(values))
    return result


def compare(charges,powers=(1,1,1,1,1)):
    if len(charges)!=4 or sum(charges):
        raise ValueError('Supply four conserved scalar vertex charges.')
    d=sum(powers)
    rows=[modulus_constant(tuple(t*p for p in charges),powers) for t in range(2*d+1)]
    discrete=sum((-1)**t*comb(2*d,t)*row for t,row in enumerate(rows))
    discrete*=F((-1)**(d-5),2**d*factorial(2*d)*prod(factorial(m) for m in powers))
    continuous=normal_coefficient(charge_cuts(EDGES,charges),tuple(m-1 for m in powers))
    assert discrete==continuous,(charges,powers,discrete,continuous)
    print(charges,powers,discrete,flush=True)
    return discrete,continuous


if __name__=='__main__':
    for charges in ((1,-1,-1,1),(1,1,-1,-1),(2,-3,1,0),(2,-1,2,-3),(1,2,3,-6)):
        compare(charges)
    for powers in ((2,1,1,1,1),(1,2,1,1,1),(1,1,1,2,2)):
        compare((1,-1,-1,1),powers)
