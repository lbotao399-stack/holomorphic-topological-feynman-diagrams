"""Check the Poisson fan against grounded Laplacians and direct polynomial integration."""
from fractions import Fraction
from math import factorial, prod

import sympy as sp

from fan import bracket, dot, taylor, zero_displacement


def path_currents(r, p, q):
    n = len(r)//2
    result = [None]*n
    forward = [prod(r[:2*i:2]) for i in range(1, n+1)]
    result[-1] = [forward[-1]*a-r[-1]*b for a, b in zip(p, q)]
    for i in reversed(range(n-1)):
        v = r[2*i+1]
        result[i] = [forward[i]*v*a+(1-v)*b for a, b in zip(p, result[i+1])]
    return result


def vertex_product(path, p, q):
    left = [p, *path]
    right = [*path, [-b for b in q]]
    return prod(-bracket(a, b) for a, b in zip(left, right))


def density(r, p, q):
    n = len(r)//2
    weights = (r[2*i-2]**(2*(n-i)+1) * r[2*i-1]**2 * (1-r[2*i-1])**(2*i-1)
               for i in range(1, n))
    return bracket(p, q)**(2*n+1) * r[-2] * r[-1]**(2*n-1) * prod(weights)


def direct(n, p, q, seed):
    t = [Fraction((i+seed)**2 % 31 + 1, i+2) for i in range(2*n+1)]
    laplace = sp.diag(*(1/x for x in [*t[1:2*n:2],t[-1]]))
    for i in range(n):
        edge = sp.zeros(n+1, 1)
        edge[i], edge[i+1] = 1, -1
        laplace += edge*edge.T/t[2*i]
    potentials = laplace.inv()*sp.Matrix([p, *([(0,0)]*(n-1)), q])
    path = [[(potentials[i,j]-potentials[i+1,j])/t[2*i] for j in (0,1)] for i in range(n)]
    r = []
    spoke = t[1]
    for i in range(n-1):
        series = spoke+t[2*i]
        r.extend([spoke/series, series/(series+t[2*i+3])])
        spoke = series*t[2*i+3]/(series+t[2*i+3])
    total = t[-3]+spoke+t[-1]
    r.extend([spoke/total, t[-1]/total])
    assert path == path_currents(r, p, q)
    # Kirchhoff currents in this chart are multi-affine, so these two values give
    # each exact derivative even on a simplex coordinate.
    columns = []
    for j in range(2*n):
        one = path_currents(r[:j]+[1]+r[j+1:], p, q)
        zero = path_currents(r[:j]+[0]+r[j+1:], p, q)
        columns.append([x-y for hi, lo in zip(one, zero) for x, y in zip(hi, lo)])
    # Each of the first n-1 inverse series/parallel charts reverses orientation.
    jac = -sp.det(sp.Matrix(columns))
    assert jac * vertex_product(path, p, q) == density(r, p, q)


def integrate(poly, variables):
    # Intervals followed by the triangle with its induced Lebesgue measure.
    result = 0
    for powers, c in sp.Poly(poly, *variables).terms():
        a,b = powers[-2:]
        result += c*factorial(a)*factorial(b)/(factorial(a+b+2)*prod(i+1 for i in powers[:-2]))
    return result


def raw_moments(n, p, q, path_z, spoke_z, order):
    r = sp.symbols(f'r:{2*n}')
    path = path_currents(r, p, q)
    spokes = [[a-b for a,b in zip(p,path[0])]]
    for last, following in zip(path,path[1:]):
        spokes.append([a-b for a,b in zip(last,following)])
    spokes.append([a+b for a,b in zip(q,path[-1])])
    phase = sum(dot(y,z) for y,z in zip(path+spokes,path_z+spoke_z))
    jac = sp.Matrix(sum(path, [])).jacobian(r).det()
    weight = sp.factor(-jac*vertex_product(path,p,q))
    phase = sp.Poly(-phase, *r)
    poly = sp.Poly(weight, *r)
    result = []
    for m in range(order+1):
        result.append(integrate(poly, r)/factorial(m))
        poly *= phase
    return result


if __name__ == '__main__':
    for n in range(1,8):
        for seed, p, q in ((3, (2,-1), (1,3)), (7, (1,2), (2,-3))):
            direct(n,p,q,seed)
    for n, order in ((1,6),(2,4),(3,1)):
        p,q = (2,-1),(-3,2)
        path_z = [(i,2-i) for i in range(n)]
        spoke_z = [(i-1,3*i-2) for i in range(n+1)]
        coeff = taylor(n,p,q,path_z,spoke_z,order)
        assert coeff == raw_moments(n,p,q,path_z,spoke_z,order)
        assert coeff[0] == zero_displacement(n,p,q)
    print('grounded currents and Jacobians through 7 loops; all 14 direct Taylor coefficients agree')
