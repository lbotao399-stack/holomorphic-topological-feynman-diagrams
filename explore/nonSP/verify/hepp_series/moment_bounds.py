"""Two-sided rational Hermite bounds for a positive measure on 1 <= U/max(t_I) <= M."""
from fractions import Fraction as F
from math import isqrt

import sympy as sp


def orthogonal_nodes(moments, degree, bits=16):
    """Choose rational squares near the nodes of the positive-moment Gauss rule."""
    if degree==0:
        return []
    matrix=sp.Matrix(degree,degree,lambda i,j:sp.Rational(moments[i+j]))
    rhs=sp.Matrix([-sp.Rational(moments[i+degree]) for i in range(degree)])
    coefficients=matrix.inv()*rhs
    x=sp.Symbol('x')
    polynomial=sp.Poly(x**degree+sum(c*x**j for j,c in enumerate(coefficients)),x)
    roots=polynomial.intervals(eps=sp.Rational(1,1<<(bits+2)))
    if len(roots)!=degree or any(m!=1 for (a,b),m in roots):
        raise ValueError('Moment measure must be positive with enough distinct support points')
    dyadic=1<<bits
    nodes=[]
    for (a,b),m in roots:
        midpoint=F((a+b)/2)
        lower=isqrt(midpoint.numerator*dyadic**2//midpoint.denominator)
        nodes.append(F(lower,dyadic)**2)
    return nodes


def hermite_integral(nodes, moments, alpha):
    """Single and double nodes give a polynomial lower or upper bound."""
    count=sum(mult for x,mult in nodes)
    rows=[]
    values=[]
    for x,mult in nodes:
        square_root=F(isqrt(x.numerator),isqrt(x.denominator))
        assert square_root**2==x
        value=square_root**(-int(2*alpha))
        rows.append([x**i for i in range(count)])
        values.append(value)
        if mult==2:
            rows.append([0]+[i*x**(i-1) for i in range(1,count)])
            values.append(-alpha*value/x)
    coefficients=sp.Matrix(rows).inv()*sp.Matrix(values)
    return sum(F(v)*m for v,m in zip(coefficients,moments))


def integer_enclosure(moments,alpha,tree_count):
    """Optimal Gauss and endpoint Lobatto values for integral denominator power."""
    moments = [sp.Rational(str(m)) for m in moments]
    r = len(moments)//2
    if not r or len(moments)!=2*r or int(alpha)!=alpha or alpha<=0:
        raise ValueError('Need an even number of moments and positive integral power.')
    if not moments[0] or tree_count==1:
        return F(moments[0]),F(moments[0])
    x = sp.Symbol('x')
    def monic(data,degree):
        if not degree:
            return sp.Integer(1)
        matrix = sp.Matrix(degree,degree,lambda i,j:data[i+j])
        assert all(matrix[:i,:i].det()>0 for i in range(1,degree+1))
        coeff = matrix.inv()*sp.Matrix([-data[i+degree] for i in range(degree)])
        return x**degree+sum(c*x**i for i,c in enumerate(coeff))
    weighted = [(1+tree_count)*moments[i+1]-moments[i+2]-tree_count*moments[i]
                for i in range(2*r-2)]
    def evaluate(polynomial):
        # On every node the polynomial remainder agrees with x**(-alpha).
        remainder = sp.Poly(sp.invert(x**int(alpha),polynomial,x),x)
        return F(sum(remainder.nth(i)*m for i,m in enumerate(moments)))
    lower = evaluate(monic(moments,r))
    upper = evaluate((x-1)*(x-tree_count)*monic(weighted,r-1))
    assert 0<lower<=upper
    return lower,upper


def enclose(moments, alpha, tree_count):
    """Use positive raw moments through an odd degree. Output rational endpoints."""
    moments=[F(str(m)) for m in moments]
    alpha=F(alpha)
    assert len(moments)%2==0 and alpha>0 and 2*alpha==int(2*alpha)
    if not moments[0]:
        return F(0),F(0)
    if tree_count==1:
        return moments[0],moments[0]
    r=len(moments)//2
    root=isqrt(tree_count*256)
    upper_x=F(root+int(root**2<tree_count*256),16)**2
    weighted=[-moments[i+2]+(1+upper_x)*moments[i+1]-upper_x*moments[i]
              for i in range(len(moments)-2)]
    for bits in (12,24,48,96):
        lower_nodes=orthogonal_nodes(moments,r,bits)
        upper_nodes=orthogonal_nodes(weighted,r-1,bits)
        if all(1<x<upper_x for x in lower_nodes+upper_nodes) and (
                len(set(lower_nodes))==r and len(set(upper_nodes))==r-1):
            break
    else:
        raise ArithmeticError('Unable to separate interpolation nodes')
    lower=hermite_integral([(x,2) for x in lower_nodes],moments,alpha)
    upper=hermite_integral([(F(1),1),(upper_x,1)]+[(x,2) for x in upper_nodes],moments,alpha)
    assert lower<=upper
    return max(lower,F(0)),upper
