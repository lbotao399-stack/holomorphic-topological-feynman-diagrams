"""Integrated splitting ratios and tetrahedron conductance cuts on the genus-five full face.

The output is positive projected area with the divided branch power on 34.
The lower-face contribution has the opposite sign. For the reduced current
j on 34 and source p=(1,0), the top row after integrating the parallel share
is [p,j]^2 dt/56 + [p,j][j,dj]/252. Its other five rows are dE_e.
"""
from fractions import Fraction
from itertools import combinations
from math import comb

import sympy as sp
from sympy.polys.rings import ring


EDGES = tuple(combinations(range(4),2))
TRANSVERSE = (1,2,3,1)


def incidence(edges,n):
    B = sp.zeros(n-1,len(edges))
    for i,(a,b) in enumerate(edges):
        if a<n-1: B[a,i]=1
        if b<n-1: B[b,i]=-1
    return B


def determinants(transverse):
    # Delete and contract 12, fixing the effective 34 conductance to one.
    # t divides the marked load and z=c12*R/(1+c12*R).
    base,t,z = ring('t,z',sp.QQ)
    B = incidence(EDGES,4)
    g = (0,*transverse,1)
    C = sp.diag(*g)
    S = (B*C*B.T).inv()
    b = B[:,0]
    R = (b.T*S*b)[0]
    d = S*b
    transfer = [C*B.T*S,C*B.T*(S-d*d.T/R)]
    transfer[1][0,:] = b.T*S/R
    def constant(x):
        return base.ground_new(sp.QQ.convert(x))
    # Reduced charges of original nodes 1,2,3 are (0,1),(-1,-1),(t,0).
    endpoints = [[[-constant(M[i,1])+t*constant(M[i,2]),constant(M[i,0]-M[i,1])]
                  for i in range(6)] for M in transfer]
    def lerp(x,y):
        return x+z*(y-x)
    currents = [[lerp(i,j) for i,j in zip(*pair)] for pair in zip(*endpoints)]
    response = [sp.eye(6)-M*B for M in transfer]
    derivatives = []
    for i,jj in enumerate(currents):
        partial = [[j.diff(v) for j in jj] for v in (t,z)]
        for l in range(1,5):
            left,right = [[constant(matrix[i,l]/g[l])*x for x in ends[l]]
                          for matrix,ends in zip(response,endpoints)]
            partial.append([lerp(x,y) for x,y in zip(left,right)])
        derivatives.append(partial)
    energies = [[sum(x*y for x,y in zip(jj,dj)) for dj in partial]
                for jj,partial in zip(currents,derivatives)]
    jx,jy = currents[-1]
    # Density with the root divided branch power two included.
    multiplier = ((jx*jx+jy*jy)/2)**3/6
    D0 = [jy*jy,*([base.zero]*5)]
    D1 = [jy*(jx*dy-jy*dx) for dx,dy in derivatives[-1]]
    def exterior_det(rows):
        partial = {0:base.one}
        for row in rows:
            nxt = {}
            for used,value in partial.items():
                for col,x in enumerate(row):
                    if not x or used & (1<<col): continue
                    new = used | (1<<col)
                    sign = (-1)**((used>>(col+1)).bit_count())
                    nxt[new] = nxt.get(new,base.zero) + sign*value*x
            partial = nxt
        return partial.get(63,base.zero)
    return tuple(multiplier*exterior_det([first,*energies[:5]]) for first in (D0,D1))


def original_density(t,s,z,transverse,charges=((1,0),(0,1),(-1,-1),(0,0))):
    """Direct derivative of all seven original (nonroot) energies."""
    B = incidence(EDGES,4)
    CB = sp.diag(0,*transverse,1)
    S = (B*CB*B.T).inv()
    b = B[:,0]
    R = (b.T*S*b)[0]
    x = z/(R*(1-z))
    edges = tuple(e for e in combinations(range(5),2) if e not in ((0,1),(0,2)))
    ground = incidence(edges,5)
    c = sp.diag(s/(1-t),s/t,x,*transverse,1-s)
    du = sp.zeros(8,7)
    # Columns are t,s,z followed by the four fixed transverse conductances.
    du[0,0],du[0,1] = s/(1-t)**2,1/(1-t)
    du[1,0],du[1,1] = -s/t**2,1/t
    du[-1,1],du[2,2] = -1,1/(R*(1-z)**2)
    for l in range(1,5):
        du[2,l+2] = x*(B[:,l].T*S*b)[0]**2/R
        du[l+2,l+2] = 1
    transfer = c*ground.T*(ground*c*ground.T).inv()
    current = transfer*sp.Matrix(charges)
    gram = current*current.T
    derivative = (sp.eye(8)-transfer*ground).multiply_elementwise(gram)*c.inv()*du
    return (gram[-1,-1]/2)**3*derivative[:7,:].det()/6


def evaluate(poly,t,z):
    return sum(Fraction(str(c))*t**i*z**j for (i,j),c in poly.items())


def integrate_ratios(first,second):
    integrated = first/56+second/252
    return sum(Fraction(str(c))/((i+1)*(j+1)) for (i,j),c in integrated.items())


def two_edge_cut():
    """Integrate c13 too. Its only finite poles are in deletion/contraction tree polynomials."""
    x=sp.Symbol('x')
    B=incidence(EDGES,4)
    deleted=B*sp.diag(0,x,2,3,1,1)*B.T
    assert deleted.det()==15*x+14
    assert (B[:,0].T*deleted.adjugate()*B[:,0])[0]==4*x+15
    def denominator(i):
        return ((15*i+14)*(4*i+15))**23
    def value(i):
        first,second=determinants((i,2,3,1))
        return sp.Rational(integrate_ratios(first,second))
    # The two tree polynomials give a common denominator for each endpoint
    # current. Four transverse derivatives can square their denominators.
    # The six coordinate derivatives include dc13, gaining two powers at infinity.
    values=[value(i)*denominator(i) for i in range(1,46)]
    numerator=sp.Poly(0,x,domain=sp.QQ)
    binomial=sp.Poly(1,x,domain=sp.QQ)
    for j in range(45):
        numerator += binomial.mul_ground(values[0])
        values=[b-a for a,b in zip(values,values[1:])]
        binomial=(binomial*sp.Poly(x-j-1,x)).mul_ground(sp.Rational(1,j+1))
    assert all(numerator.eval(i)/denominator(i)==value(i) for i in (sp.Rational(1,3),49))
    monic_numerator=numerator.mul_ground(sp.Rational(1,60**23))
    poles=(sp.Rational(14,15),sp.Rational(15,4))
    def principal(a,b):
        shifted=monic_numerator.shift(-a)
        return {23-j:sum(shifted.nth(j-r)*(-1)**r*comb(22+r,r)/(b-a)**(23+r)
                         for r in range(j+1)) for j in range(23)}
    parts=[principal(*poles),principal(*reversed(poles))]
    for probe in (0,1,sp.Rational(7,2)):
        assert sum(c/(probe+a)**j for a,part in zip(poles,parts) for j,c in part.items()) == (
            numerator.eval(probe)/denominator(probe))
    assert parts[0][1]+parts[1][1]==0
    rational=sum(c*a**(1-j)/(j-1) for a,part in zip(poles,parts) for j,c in part.items() if j>1)
    logarithm=parts[0][1]
    assert rational == sp.Rational(-71161938221472988962540094215797989902554359889,
                                    13598660191230190000557188561633690815707322368000000000000000)
    assert logarithm == sp.Rational(23062673138156617081,5980085648358972073353442597384800)
    print('Second edge integrated; coefficients of 1 and log(225/56):',rational,logarithm)
    return rational,logarithm


def direct_quadrature():
    """Integrate the original eight-edge determinant without reduced-current sensitivities."""
    # Equally spaced interior nodes, exact through degree 18 in each variable.
    nodes = [sp.Rational(i,20) for i in range(1,20)]
    vandermonde = sp.Matrix([[node**i for node in nodes] for i in range(19)])
    weights = vandermonde.inv()*sp.Matrix([sp.Rational(1,i+1) for i in range(19)])
    # s(1-s)^6 times an affine function has its exact weighted average at 2/9.
    s = sp.Rational(2,9)
    total = 0
    for t,wt in zip(nodes,weights):
        for z,wz in zip(nodes,weights):
            total += wt*wz*original_density(t,s,z,TRANSVERSE)
    return total/(56*s*(1-s)**6)


def check():
    first,second = determinants(TRANSVERSE)
    for it,ir in ((1,2),(2,4),(3,1),(4,2)):
        for ins in (1,3):
            t,s,z = (sp.Rational(a,b) for a,b in ((it,5),(ins,4),(ir,5)))
            direct = original_density(t,s,z,TRANSVERSE)
            a,b = (evaluate(poly,t,z) for poly in (first,second))
            assert direct == (1-s)**6*(s*a+s*s*b)
    integrated = first/56+second/252
    cut = integrate_ratios(first,second)
    assert cut==Fraction(465489528164329449978660994907,22686383745706107326428233439728280322801141376)
    print('Signed coordinate density integrated over s,t,z:',cut)
    print('Degrees in t,z:',integrated.degree(0),integrated.degree(1))
    print('Transverse conductances and effective root:',TRANSVERSE,1)
    return cut


if __name__=='__main__':
    import argparse
    options = argparse.ArgumentParser(description=__doc__)
    options.add_argument('--direct',action='store_true',help='repeat the cut with 361 original eight-edge determinants')
    options.add_argument('--two-edges',action='store_true',help='integrate the second tetrahedral conductance as well')
    args = options.parse_args()
    cut = check()
    if args.direct:
        assert direct_quadrature()==cut
        print('361 independent original eight-edge Jacobians integrate to the same cut.')
    if args.two_edges:
        two_edge_cut()
