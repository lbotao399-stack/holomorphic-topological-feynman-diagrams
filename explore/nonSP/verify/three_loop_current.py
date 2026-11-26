"""C², K4 with edge 12 subdivided, polynomial edge-current insertions."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb, prod
import sys
import sympy as s
from flint import fmpz_mpoly_ctx

from graph_k4_series import wedge_determinant
from reduced_k4_series import beta_finite_parts, convergence, one_dimension
from htk import brute
from htk.graph import Graph

EDGES = [(0,1),(0,2),(0,3),(1,4),(4,2),(1,3),(2,3)]
GRAPH = Graph(5,EDGES)
SOURCES = [(1,2),(2,-1),(-1,3),(-3,-4),(1,0)]
PAIRED_SOURCES = [(1,0),(0,0),(0,1),(-1,0),(0,-1)]


def data(variables, center, sources):
    split, first, second, *rho = variables
    outer = [i for i in range(4) if i != center]
    g = dict(zip(outer,(first,second,1-first-second)))
    r = dict(zip(outer,rho))
    p = {i:g[i]*(1-r[i]) for i in outer}
    denominator = sum(p.values())
    lam = [list(row) for row in sources]
    for i in range(2):
        lam[1][i] += (1-split)*lam[4][i]
        lam[2][i] += split*lam[4][i]
    mu = {i:[u+g[i]*v for u,v in zip(lam[i],lam[center])] for i in outer}
    links = {}
    for i,j in combinations(outer,2):
        k = next(iter(set(outer)-{i,j}))
        bare = [p[j]*x-p[i]*y for x,y in zip(mu[i],mu[j])]
        links[i,j] = [r[k]*x for x in bare]
        links[j,i] = [-x for x in links[i,j]]
    for i in outer:
        row = [g[i]*denominator*v for v in lam[center]]
        for j in outer:
            if i == j:
                continue
            k = next(iter(set(outer)-{i,j}))
            for d in range(2):
                row[d] -= (1-r[k])*(p[j]*mu[i][d]-p[i]*mu[j][d])
        links[center,i] = row
        links[i,center] = [-x for x in row]
    links[1,4] = [x-(1-split)*denominator*u for x,u in zip(links[1,2],lam[4])]
    links[4,2] = [x+split*denominator*u for x,u in zip(links[1,2],lam[4])]
    return denominator,[links[edge] for edge in EDGES]


def time_fractions(variables,center):
    split,first,second,*rho = variables
    outer = [i for i in range(4) if i != center]
    g = dict(zip(outer,(first,second,1-first-second)))
    r = dict(zip(outer,rho))
    link = {}
    for i in outer:
        link[tuple(sorted((i,center)))] = (1,g[i])
        j,k = sorted(set(outer)-{i})
        link[j,k] = (1-r[i],g[j]*g[k]*r[i])
    n,d = link[1,2]
    link[1,4],link[4,2] = (split*n,d),((1-split)*n,d)
    return [link[edge] for edge in EDGES]


def numerator(center,sources):
    symbols = 'split gA gB rA rB rC'.split()
    ctx = fmpz_mpoly_ctx.get(tuple(symbols))
    q,num = data(ctx.gens(),center,sources)
    for v in range(5):
        for d in range(2):
            assert sum((int(v==i)-int(v==j))*row[d] for (i,j),row in zip(EDGES,num)) == q*sources[v][d]
    fractions = time_fractions(ctx.gens(),center)
    common = prod(d for _,d in fractions)
    voltage_weights = [n*(common//d) for n,d in fractions]
    for cycle in GRAPH.cycle_matrix():
        for d in range(2):
            assert not sum(j*weight*row[d] for j,weight,row in zip(cycle,voltage_weights,num))
    rows = [q]+[x for e in GRAPH.cotree() for x in num[e]]
    columns = [rows]+[[row.derivative(i) for row in rows] for i in range(6)]
    try:
        det = wedge_determinant(columns,ctx.constant(0),ctx.constant(1))
    except KeyError as exc:
        if exc.args != (127,):
            raise
        det = ctx.constant(0)
    return q,num,det


def build(center, sources):
    q,num,det = numerator(center,sources)
    symbols = 'split gA gB rA rB rC'.split()
    xi = s.symbols(' '.join(symbols))
    ts = [n/d for n,d in time_fractions(xi,center)]
    jac = s.Matrix([ts]+[[s.diff(t,x) for t in ts] for x in xi])
    sign = None
    for point in ([F(1,3),F(1,5),F(1,3),F(1,2),F(1,4),F(3,5)],
                  [F(2,5),F(2,7),F(1,4),F(2,5),F(2,3),F(1,5)]):
        tt = [F(t.subs(dict(zip(xi,point)))) for t in ts]
        direct = brute.currents(GRAPH,tt,sources[:-1])
        qq,nn = data(point,center,sources)
        assert [[u/qq for u in row] for row in nn] == [list(row) for row in direct]
        raw = sum(F(int(co))*prod(x**int(a) for x,a in zip(point,ex)) for ex,co in det.terms())
        density = brute.density(GRAPH,tt,sources[:-1]).rational()
        orientation = density*abs(F(jac.subs(dict(zip(xi,point))).det()))*qq**7/raw
        assert orientation in (1,-1) and (sign is None or sign == orientation)
        sign = orientation
    return det,num,sign


def integrate(det,num,sign,insertions):
    poly = det
    for e,d in insertions:
        poly *= num[e][d]
    power = 7+len(insertions)
    answer,margin = integrate_numerator(poly,power,sign)
    print('insertions',insertions,'margin',margin,'coefficients',answer,flush=True)
    return answer


def integrate_numerator(poly,power,sign):
    if not poly:
        return {},None
    terms = defaultdict(F)
    for ex,coef in poly.terms():
        n,*other = map(int,ex)
        terms[tuple(other)] += sign*F(int(coef),n+1)
    D = max(power-2,max(g+h for (g,h,*r),v in terms.items() if v))
    homog = defaultdict(F)
    for (g,h,*r),v in terms.items():
        n = D-g-h
        for i in range(n+1):
            for j in range(n-i+1):
                homog[(g+i,h+j,n-i-j,*r)] += v*comb(n,i)*comb(n-i,j)
    homog = {k:v for k,v in homog.items() if v}
    margin = convergence(homog,D,power)
    answer = beta_finite_parts(one_dimension(D,homog,power))
    return answer,margin


if __name__ == '__main__':
    center = int(sys.argv[1]) if len(sys.argv)>1 else 0
    det,num,sign = build(center,SOURCES)
    print('center',center,'orientation',sign,'determinant monomials',len(det),flush=True)
    expected = [{0:-F(71,12)}, {0:F(222317,25200),3:-F(3021,700)},
                {0:-F(4003781,302400),3:F(3772,525)}]
    for powers,answer in enumerate(expected):
        assert integrate(det,num,sign,[(0,0)]*powers) == answer
    det,num,sign = build(center,PAIRED_SOURCES)
    assert integrate(det,num,sign,[]) == {}
    assert integrate(det,num,sign,[(0,0)]) == {0:F(17,9450),3:-F(1,700)}
    assert integrate(det,num,sign,[(0,1)]) == {0:-F(79,50400),3:F(1,700)}
