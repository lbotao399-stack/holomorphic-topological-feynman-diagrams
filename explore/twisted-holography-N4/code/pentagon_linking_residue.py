"""Shrinking-link test after vertex integration, linear in the finite puncture block."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial, prod

import sympy as sp

from normal_flux_homotopy import holomorphic, linking_residue, R, W


def check_angular():
    # The covectors multiplying the normal forms U and b_i in a fixed flat frame.
    # U is represented by its normal vector via contraction with the volume.
    bar = sp.Matrix(sp.symbols('b1:4'))
    hol = sp.Matrix(sp.symbols('z1:4'))
    radius_squared = bar.dot(hol)
    for i in range(3):
        angular = bar.cross(sp.eye(3)[:,i])
        assert not bar.dot(angular)
        for j in range(3):
            assert all(sp.cancel(x)==0 for x in
                       ((angular/radius_squared**2).diff(hol[j])
                        +2*bar[j]*angular/radius_squared**3))
    for genus in range(2,7):
        for order in range(2*genus+1):
            derivative=(-1)**order*prod(range(3,3+order))*R[2]**order
            assert linking_residue(W[2]**(2*genus)*derivative)==(
                factorial(2*genus) if order==2*genus else 0)


def check_pentagons():
    # Vertex 0 is genus three; the other four vertices are cubic. The eleven
    # distinct markings are nine scalar slots, one explicit flux, and one test.
    # Orient one puncture line by its angular b end and its annular l end.
    histogram=Counter()
    for bulk_marks in combinations(range(11),7):
        for single in permutations(i for i in range(11) if i not in bulk_marks):
            vertices=(bulk_marks, *((i,) for i in single))
            flux=next(i for i,marks in enumerate(vertices) if 9 in marks)
            test=next(i for i,marks in enumerate(vertices) if 10 in marks)
            genus=3 if test==0 else 0
            for edge in range(5):
                for target,source in ((edge,(edge+1)%5),((edge+1)%5,edge)):
                    if test==source:
                        histogram['annular support misses link']+=1
                    elif flux==target:
                        histogram['horizontal wedge zero']+=1
                    else:
                        # Output divergence of the tangent flux bivector is tangential.
                        # Allow the puncture BV edge one normal derivative of b at test.
                        singularity=(5+2*genus if flux==test else
                                     3+2*genus+1 if test==target else 0)
                        power=8+5-singularity
                        assert power>=2
                        histogram[power]+=1
    assert sum(histogram.values())==79200
    print('Ordered ninth-pole assignments and lower bounds on radial powers:',dict(histogram))


def check_differentiated_angular_modes():
    # Normalize on the unit sphere by d^3w U. The projection of
    # b_i wedge dbar(w_j) divided by U is (delta_ij - w_i r_j)/2.
    tested = 0
    for degree in range(9):
        for word in combinations_with_replacement(range(3), degree):
            orders = tuple(word.count(j) for j in range(3))
            test = prod(w**n for w,n in zip(W,orders))
            angular = (-1)**degree*factorial(degree+1)*prod(r**n for r,n in zip(R,orders))
            for i in range(3):
                measured = linking_residue(test*angular*(1-W[i]*R[i])/2)
                predicted = sp.Rational((-1)**degree*prod(map(factorial,orders))*
                                        (degree+2-orders[i]), (degree+2)*(degree+3))
                assert measured == predicted
                tested += 1
    for degree in range(13):
        puncture = sp.S.One
        for _ in range(degree):
            puncture = holomorphic(puncture, 2, 2)
        assert puncture == (-1)**degree*factorial(degree+1)*R[2]**degree
        # Relative to the seed annular normalization, the dual mode
        # bar(w_1) w_3^k dbar(chi) pairs with w_1 d_3^j U.
        dual = W[2]**degree*W[0]*R[0]
        for j in range(13):
            flux = (-1)**j*sp.Rational(factorial(j+2),2)*R[2]**j
            norm = 3*linking_residue(dual*flux)
            want = (-1)**degree*sp.Rational(3*factorial(degree),degree+3) if j==degree else 0
            assert norm == want
    print('Differentiated angular modes:',tested,'exact quadratic link coefficients; 169 annular pairings.')


def check_angular_projectors():
    # At the north pole evaluate the reproducing kernel for harmonics of
    # bidegree (1,k). Unitary rotations cover every point of the sphere.
    north = {symbol: int(i==2) for variables in (W,R) for i,symbol in enumerate(variables)}
    tested = 0
    for k in range(8):
        top = R[2]*W[2]**k
        if k:
            top -= sp.Rational(k,k+2)*W[2]**(k-1)
        kernel = (k+1)*(k+2)*(k+3)*top/2
        assert kernel.subs(north) == (k+1)*(k+3)
        for word in combinations_with_replacement(range(3), k):
            monomial = sp.prod(R[j] for j in word)
            for i in range(3):
                harmonic = W[i]*monomial-sp.diff(monomial,R[i])/(k+2)
                assert linking_residue(kernel*harmonic) == harmonic.subs(north)
                tested += 1
    # With y on the third unit axis the link on w_3^m exp(p w_3)
    # against dbar(w_3) sums to 2 R^2 int_0^1 t^(m+1) exp(p t) dt.
    for m in (0,2,8):
        for n in range(13):
            k = m+n
            coefficient = (k+1)*(k+3)*linking_residue((W[2]*R[2])**k*(1-W[2]*R[2])/2)
            assert coefficient/factorial(n) == sp.Rational(2,(k+2)*factorial(n))
    # For H = w_1^a r_3^k pair its primitive with bar(w_1)^(a-1)
    # w_3^k dbar(w_1). The sphere radius supplies R^(2a).
    for a in range(1,8):
        for k in range(13):
            integrand = (W[0]*R[0])**(a-1)*(W[2]*R[2])**k*(1-W[0]*R[0])/(k+2)
            assert linking_residue(integrand) == sp.Rational(2*factorial(a-1)*factorial(k),factorial(a+k+2))
    print('Complete differentiated sector:',tested,'harmonic projector tests; 39 resummed coefficients; 91 higher-degree tests.')


def check_two_degree_selection():
    # Harmonic projection of w_3^a bar(w_3)^k. Its jth trace subtraction
    # multiplies |w|^(2j) w_3^(a-j) bar(w_3)^(k-j).
    cancelled = nonzero = 0
    for a in range(1,14):
        for k in range(14):
            trace = [Fraction(1)]
            for j in range(min(a,k)):
                trace.append(-trace[-1]*Fraction((a-j)*(k-j),(j+1)*(a+k+1-j)))
            for m in (0,2,8):
                for q in range(16):
                    p = q+k+1-a-m
                    if p<0:
                        continue
                    # After the central phase selection, each projected term
                    # integrates |w_3|^(2s) (1-|w_3|^2) on the unit sphere.
                    moment = Fraction(0)
                    for j, coefficient in enumerate(trace):
                        s = q+k-j
                        moment += coefficient*(a-j)*Fraction(4,(s+1)*(s+2)*(s+3))
                    if p+q < a+k-3-m:
                        assert moment == 0
                        cancelled += 1
                    nonzero += bool(moment)
    assert cancelled and nonzero
    print('Simultaneous degree selection:',cancelled,'central-phase-compatible cancellations;',nonzero,'nonzero controls.')


def check_integrated_angular_tail():
    # On cutting the special line, number the ordinary path from angular
    # vertex 0 to annular vertex 4. Prescribing an interior test decouples
    # its two arms, even with the radial Schwinger mass of one flux vertex.
    checks = 0
    for resistances in ((1,2,3,4),(7,3,5,2)):
        for test in (1,2,3):
            boundary = [0,test,4]
            internal = [j for j in range(5) if j not in boundary]
            for flux in (1,2,3,4):
                for mass in (0,sp.Rational(1,2),11):
                    laplacian = sp.zeros(5)
                    for j,length in enumerate(resistances):
                        conductance = sp.Rational(1,length)
                        laplacian[j,j] += conductance
                        laplacian[j+1,j+1] += conductance
                        laplacian[j,j+1] -= conductance
                        laplacian[j+1,j] -= conductance
                    laplacian[flux,flux] += mass
                    coupling = laplacian.extract(boundary,internal)
                    reduced = (laplacian.extract(boundary,boundary)
                               -coupling*laplacian.extract(internal,internal).inv()*coupling.T)
                    assert reduced[0,2] == 0
                    if 0<flux<test:
                        first = sp.Rational(1,sum(resistances[:flux]))
                        second = sp.Rational(1,sum(resistances[flux:test]))
                        precision = first*(second+mass)/(first+second+mass)
                        fraction = second/(second+mass)
                    else:
                        precision = sp.Rational(1,sum(resistances[:test]))
                        fraction = 1
                    assert reduced[0,0] == precision
                    assert -reduced[0,1]/precision == fraction
                    checks += 1
    moments = 0
    for a in range(1,15):
        for k in range(15):
            # int_C3 x1^(a-1) bar(x3)^(k+1)/|x|^(2k+4)
            # times exp(-c|x|^2+u x3+v bar(x1)), after angular integration.
            # Strip pi^3 u^(k+1) v^(a-1) c^(-a-1).
            sphere = linking_residue((W[0]*R[0])**(a-1)*(W[2]*R[2])**(k+1))
            angular = sphere*sp.Rational(factorial(a),2*factorial(a-1)*factorial(k+1))
            # Independent Schwinger variable t/(c+t) gives an integral of
            # z^(k+1)(1-z)^a. Integrate its finite binomial expansion.
            beta = sum((-1)**j*sp.binomial(a,j)/sp.Integer(k+2+j) for j in range(a+1))
            assert angular == beta/factorial(k+1) == sp.Rational(factorial(a),factorial(a+k+2))
            moments += 1
    print('Unrestricted angular vertex:',checks,'path reductions and',moments,'independent Gaussian moments.')


def check_angular_test_collision():
    # Unit-sphere convolution at x=(0,0,R). Compare the two-variable
    # Taylor coefficient of the heat kernel to the Poisson difference form.
    terms = 0
    for k in range(20):
        for j in range(13):
            n = k+j+3
            sphere = linking_residue(W[0]*R[0]*(W[2]*R[2])**(k+j))
            angular = (k+1)*(k+3)*sphere/(factorial(j)*factorial(k+j))
            poisson = sp.Rational(2*(n-j)*(n-j-2),factorial(n)*factorial(j))
            assert angular == poisson
            terms += 1
    L,theta,c,s = sp.symbols('L theta c s',positive=True)
    # Untruncated Poisson second moment and the eight-derivative link.
    expectation = (L-theta)**2+(L+theta)-2*(L-theta)
    assert sp.expand(expectation-(L**2-(2*theta+1)*L+theta**2+3*theta))==0
    radius = sp.Symbol('rho',positive=True)
    radial = (2*c**3*radius**2*expectation/L**3).subs({L:c*radius**2,theta:c*s})
    angular_factor = linking_residue((W[2]*R[2])**8*(1-W[2]*R[2])/2)
    assert angular_factor == sp.Rational(1,495)
    eighth = (angular_factor*s**11*sp.diff(radial/s**2,s,8)).subs(s,theta/c).expand()
    expected = (2*factorial(9)*theta*c-(2*factorial(9)*theta+4*factorial(8)*theta**2)/radius**2
                +6*factorial(8)*theta**2/(c*radius**4))/495
    assert (eighth-expected).expand()==0
    print('Angular-test collision:',terms,'independent sphere/Poisson coefficients; axial eighth derivative verified.')


def check_axial_endpoint_tensors():
    counts = Counter()
    survivors = Counter()
    saturated = Counter()
    b_required = 0
    for marks in combinations(range(11),7):
        for singles in permutations(i for i in range(11) if i not in marks):
            vertices = (marks,*((i,) for i in singles))
            test = next(i for i,group in enumerate(vertices) if 10 in group)
            flux = next(i for i,group in enumerate(vertices) if 9 in group)
            for sign in (-1,1):
                path = tuple((test+sign*j)%5 for j in range(5))
                if test==flux:
                    counts['horizontal'] += 1
                    continue
                def rank_two(A,B):
                    tangent = {6:(1,0),7:(0,1),8:(A,B)}
                    rows = [tuple(sum(tangent[i][j] for i in vertices[v] if i in tangent)
                                  for j in range(2)) for v in path[1:4]]
                    return any(a*d!=b*c for (a,b),(c,d) in combinations(rows,2))
                # The minors are affine in A,B since mark 8 occurs once.
                if not any(rank_two(A,B) for A,B in ((0,0),(1,0),(0,1))):
                    counts['tangent rank'] += 1
                elif sum(i>=6 for i in marks)<2:
                    counts['fivefold rank'] += 1
                elif flux==path[4] and 0 in path[1:4]:
                    counts['normal phase'] += 1
                else:
                    counts['retained'] += 1
                    survivors[path.index(flux)] += 1
                    b_required += not (rank_two(0,0) or rank_two(1,0))
                    hard_pure = sum(i<6 for i in marks)
                    single_pure_interior = any(len(vertices[v])==1 and vertices[v][0]<6
                                                for v in path[1:4])
                    if hard_pure==5:
                        if single_pure_interior:
                            saturated['zero jump',path.index(flux)] += 1
                        else:
                            assert len(vertices[path[4]])==1 and vertices[path[4]][0]<6
                            saturated['constant annulus',path.index(flux)] += 1
                    else:
                        assert hard_pure==4 and flux==0
                        saturated['four pure',path.index(flux)] += 1
    assert counts == {'horizontal':6048,'tangent rank':6744,'fivefold rank':240,
                      'normal phase':216,'retained':2592}
    assert survivors == {1:720,2:720,3:720,4:432} and b_required==360
    assert saturated == {**{('zero jump',j):144 for j in (1,2,3)},('zero jump',4):432,
                         **{('constant annulus',j):216 for j in (1,2,3)},
                         **{('four pure',j):360 for j in (1,2,3)}}
    # Flux at the annulus. Interior momenta are t dw2+dz2, dw2,
    # and u dw2+A dz1+B dz2. Use relative difference and three jumps
    # as current rows, stripping the invertible proper-time matrix.
    x = sp.Matrix(sp.symbols('x1:4'))
    bar = sp.Matrix(sp.symbols('b1:4'))
    delta = sp.Matrix(sp.symbols('d1:4'))
    A,B,t,u = sp.symbols('A B t u')
    rows = sp.Matrix([[0,0,*delta],[0,1,0,t,0],[0,0,0,1,0],[A,B,0,u,0]])
    cofactor = sp.Matrix([(-1)**j*rows[:,[i for i in range(5) if i!=j]].det() for j in range(5)])
    assert cofactor == sp.Matrix([0,0,-A*delta[2],0,A*delta[0]])
    root = (bar+delta).cross(bar)
    dolbeault = rows[:,:2].row_join(rows[:,2:]*root).row_join(rows[:,2:]*x).det()
    normal = sp.Matrix(cofactor[2:]).dot(bar)
    transverse = bar[0]*delta[2]-bar[2]*delta[0]
    assert (dolbeault*normal-A**2*x.dot(delta)*transverse**2).expand()==0
    print('Axial endpoint tensor census:',dict(counts),'; flux sites',dict(survivors),
          '; B needed for tangent rank',b_required)
    print('Pure wave saturation:',dict(saturated))


def check_normal_torus_and_minors():
    # Pull back the full Dolbeault wedge at each nonangular flux site.
    def exterior(x,y,bx,by,flux_vector,rows,site):
        def normal(v,slot):
            return sp.Matrix([[*v,*[sum(v[i]*rows[j,i+2] for i in range(3))
                                    if j<slot else 0 for j in range(4)]]])
        tangent = [sp.Matrix([[0,0,0,*rows[:,i]]]) for i in (0,1)]
        common = [*tangent,normal(by.cross(bx),0),normal(x,0),normal(y,4)]
        units = sp.eye(3)
        return sum((-1)**i*flux_vector[i]*sp.Matrix.vstack(
            *common,*[normal(units[:,j],site) for j in range(3) if j!=i]).det()
                   for i in range(3))
    x,y = sp.Matrix([3,4,1]),sp.Matrix([2,1,6])
    bx,by = sp.Matrix([2,5,3]),sp.Matrix([7,3,4])
    flux = sp.Matrix([1,3,2])
    rows = sp.Matrix([[2,1,0,2,3],[-1,2,5,1,0],[1,1,1,4,2],[-2,-4,-1,0,7]])
    for g in (sp.diag(2,3,5),sp.diag(7,1,3)):
        action = sp.diag(1,1,*g.inv().diagonal())
        for site in range(1,5):
            original = exterior(x,y,bx,by,flux,rows,site)
            pulled = exterior(g*x,g*y,g.inv()*bx,g.inv()*by,g.inv()*flux,rows*action,site)
            assert original and pulled*g.det()==original
        vertex = sp.Matrix([[2,1,3,2,5]])
        assert rows.col_join(vertex).det() == (rows.col_join(vertex)*action).det()*g.det()
    # Each of 126 five-minors changes by the same exceptional determinant
    # on all nine-slot vertices with five pure scalar marks and one mixed.
    basis = list(sp.eye(5).tolist())
    checked = 0
    for multiplicities in ((3,2),(2,3)):
        def frame(m,q,r):
            fixed = [basis[0]]*multiplicities[0]+[basis[1]]*multiplicities[1]+[m,q,r]
            return sp.Matrix([*fixed,[-sum(row[j] for row in fixed) for j in range(5)]])
        reference = frame(*basis[2:])
        for m,q,r in (([3,1,1,2,0],[4,3,-1,1,2],[5,7,2,-3,1]),
                      ([-2,9,1,2,0],[8,5,1,3,2],[-3,4,2,-3,1])):
            moved = frame(m,q,r)
            scale = sp.Matrix([*basis[:2],m,q,r]).det()
            assert scale
            for group in combinations(range(9),5):
                assert moved[list(group),:].det()==scale*reference[list(group),:].det()
                checked += 1
    print('Normal character on all four flux sites; independent five-minor comparisons:',checked)


if __name__=='__main__':
    check_angular()
    check_differentiated_angular_modes()
    check_angular_projectors()
    check_two_degree_selection()
    check_integrated_angular_tail()
    check_angular_test_collision()
    check_axial_endpoint_tensors()
    check_normal_torus_and_minors()
    check_pentagons()
