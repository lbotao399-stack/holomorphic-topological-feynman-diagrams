"""Cartan restriction and incoming one-loop images for the two-dz scalar modes."""
from fractions import Fraction
from itertools import product
from math import comb,factorial,prod


def polarizations(alpha):
    """Count the shuffles of dz1, dz2 and the three commuting dtheta's."""
    target=(1,1,*alpha)
    degree=sum(target)
    coefficients={(0,0,0,0,0):1}
    for _ in range(degree):
        expanded={}
        for powers,value in coefficients.items():
            for i,limit in enumerate(target):
                if powers[i]<limit:
                    shifted=powers[:i]+(powers[i]+1,)+powers[i+1:]
                    expanded[shifted]=expanded.get(shifted,0)+value
        coefficients=expanded
    return coefficients[target]


def check_hessian():
    """Trace J^{-1}(A_+^t J + J A_- + rho(F1)/rho(H) J) entrywise."""
    import sympy as sp
    a=sp.Matrix(2,2,sp.symbols('a:4'))
    b=sp.Matrix(2,2,sp.symbols('b:4'))
    j=sp.Matrix([[0,1],[-1,0]])
    r=sp.Symbol('r')
    equation=a.T*j+j*b+r*j
    assert sp.expand((j.inv()*equation).trace()-a.trace()-b.trace()-2*r)==0
    # F=(2X1,-X2+[X1,X2],-X3+[X1,X3]) preserves Tr X1[X2,X3].
    # The two full root traces are individually nonzero at a regular Cartan H.
    for rank in (2,3):
        cartan=((1,-1),(1,2,-3))[rank-2]
        for i,k in product(range(rank),repeat=2):
            root=cartan[i]-cartan[k]
            if not root:continue
            plus,minus=(-1+root)*sp.eye(2),(-1-root)*sp.eye(2)
            assert plus.T*j+j*minus+2*j==sp.zeros(2)
            assert (2+plus.trace(),2+minus.trace())==(2*root,-2*root)


def check_invariant_polynomials():
    """Separate color polynomials on traceless diagonals and count form components."""
    import sympy as sp
    for rank in range(2,6):
        eigenvalues = sp.symbols(f'h0:{rank-1}')
        diagonal = (*eigenvalues, -sum(eigenvalues))
        traces = [sum(h**p for h in diagonal) for p in range(2,rank+1)]
        jacobian = sp.Matrix(traces).jacobian(eigenvalues)
        # A regular traceless diagonal has distinct eigenvalues.
        assert jacobian.subs(dict(zip(eigenvalues,range(1,rank)))).det() != 0
        counts = [1]+[0]*10
        for power in range(2,rank+1):
            for degree in range(power,len(counts)):
                counts[degree] += counts[degree-power]
        for degree,count in enumerate(counts):
            components = sum(comb(2,d)*comb(degree-d+2,2) for d in range(min(degree,2)+1))
            assert components*count == (2*degree**2+2*degree+1)*count
        if rank==3:
            assert counts[6]==2
    values = [[sum(x*x for x in diagonal)**3,sum(x**3 for x in diagonal)**2]
              for diagonal in [(1,-1,0),(1,1,-2)]]
    assert sp.Matrix(values).det() == 288
    print('SU(3), degree six: 170 independent components, including 30 two-ghost components.')
    values = []
    for nonzero_pairs in (1,2,3):
        diagonal = [1,-1]*nonzero_pairs+[0]*(6-2*nonzero_pairs)
        traces = [sum(x**p for x in diagonal) for p in (2,4,6)]
        values.append([traces[2],traces[0]*traces[1],traces[0]**3])
    assert sp.Matrix(values).det() == 768
    print('SO(7) and Sp(3), degree six: 255 components each, including 45 two-ghost components.')


def check_one_coordinate_fan():
    """Evaluate both dotted endpoint tensors, including a transverse displacement."""
    checked = 0
    for n in range(1,10):
        for p,q in (((2,3),(4,6)), ((-3,2),(6,-4)), ((0,1),(0,7)),
                    ((2,3),(4,7)), ((-3,2),(5,-4))):
            for transverse in (0,1,2):
                left = (p[0],transverse*p[1])
                right = (transverse*q[1],-q[0])
                contracted = sum(prod(left[i]*right[i] for i in slots)
                                 for slots in product(range(2),repeat=n))
                assert contracted == (transverse*(p[0]*q[1]-p[1]*q[0]))**n
                checked += 1
    import sympy as sp
    # Generic SU(3) Cartan with an independently polarized ghost jet.
    h1,h2,p1,p2 = sp.symbols('h1 h2 p1 p2')
    scalar=(h1,h2,-h1-h2)
    ghost=(p1,p2,-p1-p2)
    trace=lambda power:sum(h**power for h in scalar)
    polarized=lambda power:sum(p*h**power for p,h in zip(ghost,scalar))
    assert sp.expand(polarized(8)-trace(2)**2*trace(3)*polarized(1)/4
                     -(trace(2)**3/8+trace(3)**2/9)*polarized(2))==0
    assert 9*polarized(8).subs({h1:1,h2:1,p1:1,p2:1})==-4590
    print('One-coordinate boundary fan:',checked,'endpoint contractions; ninth-pole Cartan coefficient -4590.')


def check_axial_pole():
    """Constant tangent jets of the mixed and axial genus-four contacts on the same Cartan test."""
    import json
    from pathlib import Path
    import sys
    root = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(root/'code'))
    from holography.derham import Phi
    exponents = (0, 0, -1, -1, -9)
    tangent_one = (2, exponents, (0,))
    tangent_two = (2, exponents, (1,))
    form = lambda i: ((0, 0), (i,), (), (0, 0, 8))
    assert Phi({tangent_one: -factorial(8)}) == {form(1): 1}
    assert Phi({tangent_two: factorial(8), tangent_one: -factorial(8)}) == {form(0): 1, form(1): 1}
    bank = json.loads((root/'explore/kodaira-spencer/four_loop_arakelov_periods.json').read_text())
    trace = 9*(1+1+(-2)**9)
    required = [-trace*Fraction(x) for x in bank['stripped_axial_pole_72C4']]
    assert required == [Fraction(1141633397798354170399423,564970494488400000000),
                        Fraction(141363079741091006137,26156041411500000),
                        -Fraction(38974345636001537,4744860120000)]
    print('Axial four-loop Cartan test, required compensation:', *required)


if __name__=='__main__':
    check_hessian()
    check_invariant_polynomials()
    check_one_coordinate_fan()
    check_axial_pole()
    total=0
    for m in range(1,15):
        for i in range(m+1):
            for j in range(m-i+1):
                alpha=(i,j,m-i-j)
                assert polarizations(alpha)==factorial(m+2)//prod(factorial(x) for x in alpha)
                total+=1
        power=m+2
        su2=Fraction(1,2)**power+Fraction(-1,2)**power
        assert su2==(0 if m%2 else Fraction(1,2**(m+1)))
        assert 2+(-2)**power
    print('Trace polarizations, SU(2)/SU(3) Cartan tests:',total)
    print('Opposite-root Hessian and trace cancellation verified.')
