"""Integrate the decorated genus-five K4 face after the two path fractions.

The orientation is the positive projected area of genus_five_k4_cut.py.
The lower-face orientation and stable multiplier are applied separately.

For cycle matrix C and lengths L, write M=C diag(L) C^T; det(M)=U.
Put the particular t-dependent source on the root edge. With horizontal
and vertical loop currents as six coordinates, each nonroot energy has
gradient V_e=(j_ex C_e,j_ey C_e). Under a length variation the loop
coordinates vary by -diag(M,M)^(-1) V_e. At fixed lengths their t
derivative is -diag(M,M)^(-1) (L_root C_root,0). The explicit current
on the root has the extra differential (dt,0).

Let w be the six signed complementary minors of the five columns V_e.
Then the five-energy length determinant is
-w^T diag(M,M) w/U^2. The root one-form times 504 is
7 j_y^2 dt + 2 j_y (-j_y C_root,j_x C_root) applied to the loop
differentials. The Kirchhoff equation gives w^T V_root=0. Expanding
the remaining six-dimensional determinant and including the sign from
five resistance inversions yields a positive numerator. In the code
the w entries are represented by the integer polynomials U^3 w/L_root;
division by L_root U^2 in the five-current numerator is checked exactly.
The root scale is set to one only after forming the homogeneous density.
The root energy cube contributes (F_x^2+F_y^2)^3/(48 U^6).
The x-block reduced minors have t-degree at most two and the y-block
minors have degree at most three. The final numerator has t-degree at
most twelve. Its length degree is 42, giving the projective form U^-16.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from math import factorial, lcm, prod

from flint import fmpz_mat
from betti_series_parallel import coefficient
from mixed_k4_interior import CTX, PERMS, trees


# Fundamental cycles with chords 01,02,12 and tree 03,13,23.
CYCLES = ((1,0,-1,0,1,0), (0,1,-1,0,0,1), (0,0,0,1,-1,1))
TRIANGLES = ((0,1,3),(0,2,4),(1,2,5),(3,4,5))


@cache
def minor_patterns():
    """Laplace-expand the six complementary minors of the 6-by-5 current matrix."""
    patterns = []
    for omitted in range(6):
        first = [i for i in range(3) if i!=omitted]
        second = [i for i in range(3) if i+3!=omitted]
        block = []
        for cols in combinations(range(5),len(first)):
            rest = tuple(e for e in range(5) if e not in cols)
            a = int(fmpz_mat([[CYCLES[i][e] for e in cols] for i in first]).det())
            b = int(fmpz_mat([[CYCLES[i][e] for e in rest] for i in second]).det())
            sign = (-1)**(omitted+sum(cols)-sum(range(len(first))))
            if a*b:
                block.append((sign*a*b,cols,rest))
        patterns.append(block)
    return patterns


def poly_product(left,right):
    """Product in the parallel splitting fraction, with length-polynomial coefficients."""
    terms = [CTX.constant(0) for _ in range(len(left)+len(right)-1)]
    for i,x in enumerate(left):
        for j,y in enumerate(right):
            terms[i+j] += x*y
    return terms


def flag_weights(node,scale,u,x,y,p,alternate):
    """Primitive summed over a flag at each empty vertex, using squared currents."""
    norm = [a*a+b*b for a,b in zip(x,y)]
    # With load p at 0 the two outward flows are t*p-s*j and (1-t)*p+s*j.
    p2 = sum(i*i for i in p)
    dot = p[0]*x[5]+p[1]*y[5]
    left = [p2*(node*u)**2+norm[1]+norm[3],-2*node*u*dot,norm[5]]
    right = [p2*((scale-node)*u)**2+norm[2]+norm[4],2*(scale-node)*u*dot,norm[5]]
    root = [norm[5],-2*norm[5],norm[5]]
    if alternate:
        h = poly_product([6*norm[0],-6*norm[0]],poly_product(
            [a+r for a,r in zip(left,root)],[b+r for b,r in zip(right,root)]))
    else:
        # The root primitive produces 1/(1+number of flags on that root).
        ab = poly_product(left,right)
        one_root = poly_product(root,[a+b for a,b in zip(left,right)])
        two_roots = poly_product(root,root)
        h = poly_product(root,[6*z+3*i+2*j for z,i,j in zip(ab,one_root,two_roots)])
    # Integrate s(1-s) and s^2 against each coefficient of 6 times the primitive.
    metric = sum((2520//((i+2)*(i+3)))*v for i,v in enumerate(h))
    radial = sum((2520//(i+3))*v for i,v in enumerate(h))
    assert all(2520%((i+2)*(i+3))==0 for i in range(len(h)))
    return metric,radial


def fixed_numerator(node,scale,alternate=False,all_flags=False,permuted=False):
    """Integer numerator with both source coordinates scaled by scale, t=node/scale."""
    u,flow = trees()
    lengths = CTX.gens()
    x = [node*f[2]-scale*f[1] for f in flow]
    y = [scale*(f[0]-f[1]) for f in flow]
    p = (-1,-1) if permuted else (1,0)
    if permuted:
        x,y = [b-a for a,b in zip(x,y)],[-a for a in x]
    root_norm = x[5]**2+y[5]**2
    divisor = lengths[5]*u**2
    cofactor = []
    for block in minor_patterns():
        raw = sum((c*prod(x[e] for e in left)*prod(y[e] for e in right)
                   for c,left,right in block),CTX.constant(0))
        reduced,remainder = divmod(raw,divisor)
        assert not remainder
        cofactor.append(reduced)
    wx,wy = [[sum(CYCLES[i][e]*cofactor[shift+i] for i in range(3))
               for e in range(6)] for shift in (0,3)]
    assert x[5]*wx[5]+y[5]*wy[5]==0
    metric = sum(length*(a*a+b*b) for length,a,b in zip(lengths,wx,wy))
    area = (p[0]*y[5]-p[1]*x[5])**2
    direction = (p[0]*wx[5]+p[1]*wy[5])**2
    if all_flags:
        a,b = flag_weights(node,scale,u,x,y,p,alternate)
        return lengths[0 if alternate else 5]*(a*area*metric+b*lengths[5]*root_norm*direction)
    if alternate:
        # The stationary-current identity sum_e L_e dE_e=0 relates projected
        # areas by L_12/L_34. The divided root primitive changes by 3 E_12/E_34.
        # Before the parallel integral L_34=L_5/(1-s) and E_34=(1-s)^2 E_5.
        # Integration of s(1-s)^5 and s^2(1-s)^5 gives 1/42 and 1/168.
        return 9*lengths[0]*(x[0]**2+y[0]**2)*root_norm**2*(
            3*area*metric+lengths[5]*root_norm*direction)
    return lengths[5]*root_norm**3*(7*area*metric+2*lengths[5]*root_norm*direction)


def integrated_numerator(alternate=False,all_flags=False,permuted=False):
    """Degree <=12 in the splitting fraction. Integrate the interpolant exactly."""
    from sympy import Rational, integrate, prod as sprod, symbols
    t = symbols('t')
    weights = [F(integrate(sprod((12*t-j)*Rational(1,k-j) for j in range(13) if j!=k),
                           (t,0,1))) for k in range(13)]
    scale = lcm(*(w.denominator for w in weights))
    poly = CTX.constant(0)
    for k,w in enumerate(weights):
        poly += int(w*scale)*fixed_numerator(k,12,alternate,all_flags,permuted)
    return poly,scale*12**18*(8*6*2520 if all_flags else 48*504)


def check_original_cut(poly,denominator):
    """Integrate the resistance numerator over the earlier c12 cut in original coordinates."""
    from sympy import Matrix, Rational
    from genus_five_k4_cut import TRANSVERSE, EDGES, incidence
    B = incidence(EDGES,4)
    S = (B*Matrix.diag(0,*TRANSVERSE,1)*B.T).inv()
    R = F((B[:,0].T*S*B[:,0])[0])
    u,_ = trees()
    nodes = [Rational(i,20) for i in range(1,20)]
    weights = Matrix([[x**i for x in nodes] for i in range(19)]).inv()*Matrix(
        [Rational(1,i+1) for i in range(19)])
    total = F()
    for node,weight in zip(nodes,weights):
        z = F(node)
        times = (R*(1-z)/z,*(F(1,g) for g in TRANSVERSE),F(1))
        scale = lcm(*(i.denominator for i in times))
        point = tuple(int(scale*i) for i in times)
        density = F(scale**6*int(poly(*point)),denominator*int(u(*point))**16)
        dz_dc = R*(1-z)**2
        total += F(weight)*density*prod(i*i for i in times)/dz_dc
    assert total==F(465489528164329449978660994907,22686383745706107326428233439728280322801141376)


def check_grounded():
    from genus_five_k4_cut import EDGES, incidence, original_density
    import sympy as sp
    u,_ = trees()
    B = incidence(EDGES,4)
    for times,node,scale in (((2,3,5,7,11,1),2,5),((5,1,2,3,1,1),1,4),((1,2,3,2,1,1),3,4)):
        conductance = tuple(sp.Rational(1,k) for k in times)
        transverse = conductance[1:5]
        base = B*sp.diag(0,*transverse,1)*B.T
        R = (B[:,0].T*base.inv()*B[:,0])[0]
        z = conductance[0]*R/(1+conductance[0]*R)
        dz = R/(1+conductance[0]*R)**2
        t,s = sp.Rational(node,scale),sp.Rational(2,9)
        # One weighted sample is exact for the affine residual in the parallel share.
        direct = original_density(t,s,z,transverse)*dz/(56*s*(1-s)**6*prod(k*k for k in times))
        measured = F(int(fixed_numerator(node,scale)(*times)),
                     scale**18*48*504*int(u(*times))**16)
        assert measured==direct,(measured,direct)


def check_flag_jacobian(permuted=False):
    """Nine original-graph samples integrate the degree-eight parallel-share primitive."""
    import sympy as sp
    from genus_five_k4_cut import incidence, original_density
    original = tuple(e for e in combinations(range(5),2) if e not in ((0,1),(0,2)))
    original_inc = incidence(original,5)
    times = (2,3,5,7,11,1)
    conductance = tuple(sp.Rational(1,k) for k in times)
    transverse = conductance[1:5]
    B = incidence(tuple(combinations(range(4),2)),4)
    R = (B[:,0].T*(B*sp.diag(0,*transverse,1)*B.T).inv()*B[:,0])[0]
    z = conductance[0]*R/(1+conductance[0]*R)
    dz = R/(1+conductance[0]*R)**2
    nodes = [sp.Rational(i,10) for i in range(1,10)]
    weights = sp.Matrix([[x**i for x in nodes] for i in range(9)]).inv()*sp.Matrix(
        [sp.Rational(1,i+1) for i in range(9)])
    t = sp.Rational(2,5)
    charges = ((-1,-1),(1,0),(0,1),(0,0)) if permuted else ((1,0),(0,1),(-1,-1),(0,0))
    direct = [0,0]
    for s,w in zip(nodes,weights):
        c = sp.diag(s/(1-t),s/t,*conductance[:5],1-s)
        currents = c*original_inc.T*(original_inc*c*original_inc.T).inv()*sp.Matrix(charges)
        energy = [row.dot(row)/2 for row in (currents.row(i) for i in range(8))]
        jac = 6*original_density(t,s,z,transverse,charges)/energy[7]**3
        for left,pair in enumerate(original):
            if 3 not in pair: continue
            for right,pair in enumerate(original):
                if 4 not in pair: continue
                factors = energy[left]*energy[right]
                direct[0] += w*jac*energy[7]*factors/(1+(left==7)+(right==7))
                direct[1] += w*jac*energy[2]*factors*c[7,7]/c[2,2]
    u,_ = trees()
    for alt,value in enumerate(direct):
        measured = F(int(fixed_numerator(2,5,bool(alt),True,permuted)(*times)),
                     5**18*8*6*2520*int(u(*times))**16)
        assert measured==value*dz/prod(k*k for k in times)


def folded_weights(poly,q=16):
    """Fold convergent K4 periods by all 24 vertex permutations."""
    weights = defaultdict(int)
    orbit_cache = {}
    for raw,c in poly.terms():
        powers = tuple(int(i)+1 for i in raw)
        assert sum(powers)==3*q and min(powers)>=1 and max(powers)<q
        assert all(sum(powers[j] for j in tri)>q for tri in TRIANGLES)
        if powers not in orbit_cache:
            orbit = [tuple(powers[i] for i in p) for p in PERMS]
            representative = min(orbit)
            for key in orbit:
                orbit_cache[key]=representative
        weights[orbit_cache[powers]] += int(c)
    return {k:v for k,v in weights.items() if v}


@cache
def moment(powers):
    """Choose the tetrahedral cut with a short triangle recurrence."""
    from core import cut
    def cost(t):
        a,b,c,d,e = t[1:]
        return (min(a+b,c+d)+e,e,t)
    cut_powers = min((tuple(powers[j] for j in p) for p in PERMS),key=cost)
    q,rem = divmod(sum(powers),3)
    assert not rem
    factor = F(prod(factorial(k-1) for k in powers),factorial(q-1))
    return tuple(factor*c for c in cut(cut_powers))


def full_graph(alternate=False,all_flags=False,permuted=False):
    if permuted and not all_flags:
        raise ValueError("The permuted marking uses the full stable flag sum.")
    poly,denominator = integrated_numerator(alternate,all_flags,permuted)
    weights = folded_weights(poly)
    if not alternate and not all_flags:
        check_original_cut(poly,denominator)
        assert (len(poly),len(weights))==(40489,5648)
    full_face = [F(),F()]
    for p,c in weights.items():
        value = moment(p)
        # Reciprocal lengths exchange trees with cotrees, with dual edges opposite.
        dual = tuple(16-i for i in reversed(p))
        assert value==moment(dual)
        for j,x in enumerate(value):
            full_face[j] += c*x/denominator
    original = tuple(e for e in combinations(range(5),2) if e not in ((0,1),(0,2)))
    root = 2 if alternate else 7
    charges = (((-1,-1),(1,0),(0,1)) if permuted else ((1,0),(0,1),(-1,-1)))+((0,0),(0,0))
    if all_flags:
        proper = F()
        for left,pair in enumerate(original):
            if 3 not in pair: continue
            for right,pair in enumerate(original):
                if 4 not in pair: continue
                normal = tuple(int(e==left)+int(e==right) for e in range(8))
                proper += (1+(left==right))*coefficient(original,charges,normal,root=root,proper=True)
    else:
        proper = coefficient(original,charges,(0,)*7+(2,),root=root,proper=True)
        if not alternate:
            assert full_face==[F(99133977499,26125187194397491200000),-F(1,3692000632320)]
            assert proper==F(17914642896173,15342609934164344832000)
    graph = (proper-full_face[0],-full_face[1])
    if permuted:
        assert graph==(F(145775286798138013,542914046383572864000000),-F(9817,17998503082560))
    elif all_flags:
        assert graph==(F(3963229030865953661,147672620616331819008000000),-F(343723,719940123302400))
    else:
        assert graph==(F(3273677772720991,2812811821263463219200000),F(1,3692000632320))
    print('Omitted energy:',original[root],'convergent monomials and orbits:',len(poly),len(weights))
    print('Positive projected area; rational and zeta(3):',*full_face)
    print('Proper faces:',proper)
    print('All faces; rational and zeta(3):',*graph)
    if all_flags:
        print('All flags with stable multiplier 1/48:',*(v/48 for v in graph))
    return graph


if __name__=='__main__':
    check_grounded()
    print('Three original eight-edge Jacobians agree at unequal metrics and splitting fractions.',flush=True)
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all-flags',action='store_true',help='sum all 16 flag placements')
    parser.add_argument("--permuted",action="store_true",help="use the (4,1,0) marking and all its flags")
    args = parser.parse_args()
    flags = args.all_flags or args.permuted
    if flags:
        check_flag_jacobian(args.permuted)
    assert full_graph(False,flags,args.permuted)==full_graph(True,flags,args.permuted)
