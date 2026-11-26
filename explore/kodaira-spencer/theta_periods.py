"""Integrated labeled KS theta contacts in the repository htk orientation."""
from fractions import Fraction as F
from math import factorial, prod
import sympy as sp

from check_bulk import complement, edge_rows, twoloop


P = [(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)]


def ordered_moment(poly, variables, groups):
    value = 0
    for powers,c in sp.Poly(poly,*variables).terms():
        divisor = 1
        for group in groups:
            tail = 0
            for k in reversed(group):
                tail += 1+powers[k]
                divisor *= tail
        value += c/divisor
    return value


def barycentric_moment(poly, variables, groups):
    # For ordered tail sums s_i use independent barycentric increments instead.
    replacements = {}
    for group in groups:
        for k,i in enumerate(group):
            replacements[variables[i]] = sum(variables[j] for j in group[k:])
    value = 0
    for powers,c in sp.Poly(poly.xreplace(replacements),*variables).terms():
        for group in groups:
            c *= sp.Rational(prod(factorial(powers[j]) for j in group),
                             factorial(len(group)+sum(powers[j] for j in group)))
        value += c
    return value


def simplex_moment(powers):
    return F(prod(factorial(k) for k in powers),factorial(2+sum(powers)))


def numerator(th, src, potential):
    # Q is quadratic in the vector L. Recover its six coefficients by polarization.
    qs = [twoloop.Q_form(th,src,p) for p in P[:3]]
    q = sum(c*x**2 for c,x in zip(qs,potential))
    for i in range(3):
        for j in range(i):
            unit_sum = [a+b for a,b in zip(P[i],P[j])]
            q += (twoloop.Q_form(th,src,unit_sum)-qs[i]-qs[j])*potential[i]*potential[j]
    return sp.factor(q)


def data(shape):
    th = twoloop.Theta(*twoloop.theta_graph(*shape))
    groups = [[v-1 for v in th.interior[path]] for path in 'ABC']
    src = [(0,0,0),*P,(0,0,0)]
    s = sp.symbols('s:4')
    potential = [sum(s[j]*P[j][k] for j in range(4)) for k in range(3)]
    q = numerator(th,src,potential)
    return th,groups,src,s,q


def blocks(shape, moment=ordered_moment):
    th,groups,src,s,q = data(shape)
    degrees = [4-k for k in shape]
    base = simplex_moment(degrees)
    ranks = [[base*s[i]*s[j] for j in range(4)] for i in range(4)]
    for k,group in enumerate(groups):
        reduced = degrees[:]
        reduced[k] -= 1
        mean = simplex_moment(reduced)
        for i in group:
            for j in group:
                ranks[i][j] += mean*s[max(i,j)]*(1-s[min(i,j)])
    return sp.Matrix(4,4,lambda i,j:moment(q*q*ranks[i][j],s,groups))


def direct(shape,seed):
    th,groups,src,s,q = data(shape)
    w = [F(seed,15+seed),F(3,10)]
    w.append(1-sum(w))
    x = {}
    sigma = []
    for strand,(path,k) in enumerate(zip('ABC',shape)):
        vec = [1+(i+seed+strand)**2 for i in range(k)]
        x[path] = [F(xi,sum(vec)) for xi in vec]
        sigma.extend(sum(x[path][i:]) for i in range(1,k))
    t = twoloop.chart_t(th,w,x)
    y = twoloop.currents(th.nv,th.edges,th.root,t,src)
    phi = []
    for i in range(4):
        same = next(k for k,group in enumerate(groups) if i in group)
        def response(j):
            r = sigma[i]*sigma[j]
            if j in groups[same]:
                r += sigma[max(i,j)]*(1-sigma[min(i,j)])/w[same]
            return r
        phi.append([sum(response(j)*P[j][k] for j in range(4)) for k in range(3)])
    field = [tuple(v) for v in ((0,1,0),(0,0,1),(1,0,0),(1,-1,0))]
    test = twoloop.det([[*u,sum(a*b for a,b in zip(u,z))] for u,z in zip(field,phi)])
    if not test:
        raise AssertionError('choose another nonsingular test')
    weight,rem = complement(edge_rows(th.nv,th.edges,y))
    ext = [[sum(a*b for a,b in zip(u,col)) for col in zip(*rem[3*v:3*v+3])]
           for v,u in enumerate(field,1)]
    odd = weight*twoloop.det(ext)
    reference = F(q.subs(dict(zip(s,sigma))))*prod(w)
    conversion = odd/(reference*test)
    assert abs(conversion) == 1
    scalar = th.order_sign*(-1)*twoloop.sgn_detBC(th.nv,th.edges,th.root,th.C3())**3
    chart_q = twoloop.closed_chart_density(th,w,x,src,3,0)[0]
    assert chart_q == -th.order_sign*F(q.subs(dict(zip(s,sigma))))*prod(a**(3-n) for a,n in zip(w,shape))
    return scalar*conversion


if __name__ == '__main__':
    for shape in ((1,3,3),(2,2,3)):
        orientations = {direct(shape,k) for k in (1,2,4)}
        assert len(orientations) == 1
        orient = orientations.pop()
        m = orient*blocks(shape)
        assert m == orient*blocks(shape,barycentric_moment)
        unit = sp.ilcm(*(x.q for x in m))
        print(shape, m*unit, '/',unit)
