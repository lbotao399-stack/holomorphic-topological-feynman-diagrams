"""Two-loop null three-field projector. Coefficients strip w^3/(4*pi)^4."""
from itertools import permutations, product
from collections import Counter
import sympy as s

Q = s.Rational


def triangle_masters():
    c, v, y, t = s.symbols('c v y t')
    a, b = (1-c)*(1-v), (1-c)*v
    # Triangle k^2, (k+P)^2, (k-q)^2. Mean k=c*q-b*P,
    # Delta=c*(1-c)*(q+v*P)^2. P^2=0.
    # For three inner cubics only terms with one contraction in the triangle
    # can have degree six after the second integration.
    left, right = c-Q(1,2)+a, Q(1,2)-a
    first = left*(1-v**2) + right*v*(2-v)
    master_210 = -s.integrate((1-c)*first, (v,0,1), (c,0,1))/96

    # One quartic, three cubics. Its directed color contraction is X1-X2.
    k, q, p, h = s.symbols('k q p h')
    r = -p-h
    x1 = (-q+h-r)*(q-2*k-p)
    x2 = (k+p-r)*(-2*q+h+k)
    polynomial = s.expand(x1-x2)
    # Scalar letters here encode metric products, so all products of
    # external vectors p=P, h=S vanish, including p*h.
    def metric_projection(expr):
        return s.expand(expr).subs({p**2:0, h**2:0, p*h:0})
    shifted = metric_projection(polynomial.subs(k,c*q-b*p))
    aq = -c**2 + 4*c - 1
    al, am = 2*b*c-4*b-4*c+6, 4-6*c
    assert s.expand(shifted-(aq*q**2+al*p*q+am*q*h)) == 0
    linear = s.diff(polynomial,k).subs(k,c*q-b*p).expand()
    assert s.expand(linear - ((4-2*c)*q+(2*b-4)*p-6*h)) == 0

    # A contraction in the outer q integral uses variance 1/(2 T_q).
    # A contraction in the inner k integral uses variance 1/(2 T_k).
    # The remaining universal radial Gamma(2 epsilon) contributes 1/2.
    outer = v*(10*aq*y**2*t*v + 2*am*y*t*v - al*y**2)
    inner = (4*(c+b)-8*b*c)*y**2*v + 6*c*y*t*v**2 - 6*b*y*v
    def moment(poly):
        simplex = s.integrate(poly,(y,0,1-t),(t,0,1))
        return s.integrate((1-c)*simplex,(v,0,1),(c,0,1))/4
    pieces = [moment(outer), moment(inner)]
    assert pieces == [-Q(7,864), -Q(11,864)]
    assert master_210 == -Q(1,576)
    return master_210, sum(pieces)


def colors():
    f = {}
    for p in permutations((0,1,2)):
        f[p] = (-1)**sum(p[i] > p[j] for i,j in ((0,1),(0,2),(1,2)))
    def eps(i,j,k):
        return f.get((i,j,k),0)
    # External P,S,R have colors 0,1,2 in the order used below.
    total_300 = total_210 = total_200 = 0
    total_110 = [0,0,0]
    for a,b,c,d,i,j in product(range(3),repeat=6):
        total_300 += eps(a,b,c)*eps(d,b,c)*eps(0,a,i)*eps(1,i,j)*eps(2,j,d)
        total_210 += eps(d,b,c)*eps(a,i,c)*eps(0,b,i)*eps(1,a,j)*eps(2,j,d)
        # C1 at the quartic in the double-edge graph.
        total_200 += eps(0,a,j)*eps(b,c,j)*eps(d,b,c)*eps(1,a,i)*eps(2,i,d)
        # C1, C2, C3 in the remaining graph. The quartic colors are
        # (R,a,i,c); the three cubic tensors are common to all channels.
        cubic = eps(d,b,c)*eps(0,b,i)*eps(1,a,d)
        total_110[0] += eps(2,a,j)*eps(i,c,j)*cubic
        total_110[1] += eps(2,i,j)*eps(c,a,j)*cubic
        total_110[2] += eps(2,c,j)*eps(a,i,j)*cubic
    assert (total_300,total_210,total_200,total_110) == (2,1,-2,[1,-1,0])
    return total_300,total_210,total_200,total_110


def multiplicities():
    results = []
    for subdivisions, external_junction in (
        ((3,0,0),False), ((2,1,0),False), ((2,0,0),True), ((1,1,0),True)
    ):
        edges = Counter()
        next_vertex = 2
        for n in subdivisions:
            path = [0]+list(range(next_vertex,next_vertex+n))+[1]
            next_vertex += n
            for pair in zip(path,path[1:]):
                edges[tuple(sorted(pair))] += 1
        vertex_aut = 0
        for perm in permutations(range(next_vertex)):
            if external_junction and perm[0] != 0:
                continue
            if all(edges[tuple(sorted((perm[i],perm[j])))] == count for (i,j),count in edges.items()):
                vertex_aut += 1
        edge_aut = s.prod(s.factorial(count) for count in edges.values())
        results.append(s.factorial(3)/(vertex_aut*edge_aut))
    assert results == [Q(3,2),3,3,3]
    return results


def direct_quartic_gaussian():
    # Independent contractions in four real Cartesian coordinates. The two
    # complex external vectors span a null plane with omega(P,S)=2.
    P = s.Matrix([1,0,s.I,0])
    S = s.Matrix([0,1,0,-s.I])
    R = -P-S
    J = s.diag(*([s.Matrix([[0,1],[-1,0]])]*2))
    def omega(i,j):
        return (i.T*J*j)[0]
    def numerator(k,q):
        x1 = (-q+S-R).dot(q-2*k-P)
        x2 = (k+P-R).dot(-2*q+S+k)
        return omega(P,k)*omega(q,k)*omega(S,q)*(x1-x2)

    z = s.symbols('z')
    for c,v,y,t in ((Q(1,3),Q(2,3),Q(1,4),Q(1,2)),(Q(2,7),Q(1,5),Q(3,8),Q(1,3))):
        b = (1-c)*v
        q = y*S-t*v*P
        k = c*q-b*P
        # The graph has a conditional k fluctuation and an independent q
        # fluctuation. One covariance is the sum of the z^2 coefficients
        # along the four orthogonal fluctuation directions.
        inner = outer = 0
        for unit in s.eye(4).columnspace():
            inner += s.expand(numerator(k+z*unit,q)).coeff(z,2)
            outer += s.expand(numerator(k+c*z*unit,q+z*unit)).coeff(z,2)
        aq,al,am = -c**2+4*c-1,2*b*c-4*b-4*c+6,4-6*c
        assert inner/8 == (4*(c+b)-8*b*c)*y**2*t*v+6*c*y*t**2*v**2-6*b*y*t*v
        assert outer/8 == b*c*(10*aq*y**2*t*v+2*am*y*t*v-al*y**2)


def main():
    direct_quartic_gaussian()
    master_210, master_110 = triangle_masters()
    # For the doubled edge, X1=(-q-2P).(q-2k). The bubble with one
    # omega(q,k) gives q^2 B omega(P,q)/(d-1). Compare cubic_bubble.py.
    master_300 = Q(1,576)
    master_200 = -4*master_300
    assert master_110 == -Q(1,48)
    c300,c210,c200,c110 = colors()
    g300,g210,g200,g110 = multiplicities()
    five_cubic = g300*c300*master_300 + g210*c210*master_210
    # Both one-quartic graphs use omega(q,-k) at the internal cubic.
    # After the contractions C_2=C_3=-C_1/2 in the double-edge graph,
    # the channel sum is (3/2) C_1 X_1. Each quartic has coefficient 1/6.
    one_quartic = -g200*Q(1,4)*c200*master_200-g110*Q(1,6)*c110[0]*master_110
    assert five_cubic == one_quartic == 0
    print('Oriented residues (300,210,200,110):', master_300,master_210,master_200,master_110)
    print('Five cubic vertices:',five_cubic,'; one quartic and three cubics:',one_quartic)


if __name__ == '__main__':
    main()
