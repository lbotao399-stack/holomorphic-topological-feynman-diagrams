"""Independent symbolic contraction and integration of four labelled theta_166 flows."""

import sys
from fractions import Fraction
from math import prod
from pathlib import Path

from sympy import QQ
from sympy.polys.rings import ring

sys.dont_write_bytecode = True
verify = Path(__file__).resolve().parents[1] / "verify"
sys.path.insert(0, str(verify))
from vcommon import theta_graph, D_theta, kirchhoff
from v_integrand import flows_m1, grouped_integrand


def prefactor(b, c, case):
    common = (1+b)**3 * (1-c)**4
    if case == (1,0):
        return -common * (2*b*b + 4*b*c - b + 4*c*c - 3*c)
    if case == (11,0):
        return -common * (b+2*c) * (2*b+1)
    if case == (1,11):
        return common * (2*b*b + 4*b*c + b + 4*c*c - c)
    if case == (0,11):
        return common * (2*b+1) * (b+2*c-1)
    if case == "sum":
        return common * (2*c-1)
    raise ValueError(case)


def signed_numerator(b, c, case):
    minor = (b[2]-b[4])*(c[3]-c[4]) - (b[3]-b[4])*(c[2]-c[4])
    return prefactor(b[0],c[0],case) * minor


def fixture(order):
    n, edges, paths = theta_graph(1, 6, 6)
    oriented, types = next((fl, typ) for fl, typ in flows_m1(n, edges)
                           if typ[0] == "c" and all(h != 1 for t, h in fl))
    assert oriented == [(11,0),(1,0),(1,2),(2,3),(3,4),(4,5),(5,11),
                        (6,0),(7,6),(8,7),(9,8),(10,9),(11,10)]
    eye = [[int(i==j) for j in range(5)] for i in range(5)]
    lam = {0:[0]*5,11:[0]*5}
    lam.update({v+1:eye[v] for v in range(5)})
    lam.update({v+6:[-x for x in eye[i]] for v,i in enumerate(order)})
    assert all(sum(c) == 0 for c in zip(*lam.values()))
    ext = {v:dict(u=eye[0], vg=eye[1], lam=lam[v]) for v in range(n)}
    return n, edges, paths, oriented, types, lam, ext


def check_full_word_polynomial(order):
    n, edges, paths, _, _, lam, ext = fixture(order)
    _,*gen = ring("b0,b1,b2,b3,b4,c0,c1,c2,c3,c4", QQ)
    b,c = gen[:5],gen[5:]
    y = [[-x-z for x,z in zip(b,c)]]
    for start,path in zip((b,c),paths[1:]):
        current = start
        y.append(current)
        for v in path[0]:
            current = [old+injected for old,injected in zip(current,lam[v])]
            y.append(current)
    actual_sum = gen[0].ring.zero
    cases = []
    for oriented, types in flows_m1(n,edges):
        c_vertex = next(v for v in types if types[v] == "c")
        if not (c_vertex == 1 or
                (c_vertex in (0,11) and all(h != 1 for _,h in oriented))):
            continue
        case = c_vertex, oriented[0][0]
        cases.append(case)
        actual = grouped_integrand(n,edges,paths,oriented,types,y,ext)
        assert actual == signed_numerator(b,c,case)
        actual_sum += actual
        for r in (13,17,19):
            t = [Fraction(1+((e*e+3*e+5)%r),1+((5*e+1)%7)) for e in range(13)]
            physical = kirchhoff(n,edges,t,lam)
            assert (signed_numerator(physical[1],physical[7],case) ==
                    grouped_integrand(n,edges,paths,oriented,types,physical,ext))
    assert len(cases) == len(set(cases)) == 4 and actual_sum == signed_numerator(b,c,"sum")
    print(order, "four polynomial identities in ten current components and twelve Kirchhoff checks;", 
          "D_theta =", D_theta(lam,paths), flush=True)
    return cases, D_theta(lam,paths)


def rising(a,n):
    return prod((a+i for i in range(n)), start=Fraction(1))


def mean(order,case):
    # Given the maximum s on a path, the five descending coordinate means
    # are (5,4,3,2,1)*s/5. C orders only permuting the last three slots
    # leave s,t as the first components. The unscaled minor of these means:
    assert order[:2] == (0,1)
    B = [5,4,3,2,1]
    C = [B[order.index(i)] for i in range(5)]
    minor = (B[2]-B[4])*(C[3]-C[4]) - (B[3]-B[4])*(C[2]-C[4])
    if not minor:
        return Fraction(0)
    _,s,t,h,v = ring("s,t,h,v", QQ)
    # h=wB+wC, v=wB/h are independent Beta(1,11/2) and Beta(1/2,1/2).
    b = (h*v-1)*s-h*v*t
    c = h*(1-v)*s+(1-h*(1-v))*t
    polynomial = prefactor(b,c,case)
    result = Fraction(0)
    for (i,j,k,l), coef in polynomial.items():
        # -wA*s*t*minor/25, s,t independently Beta(5,1).
        path = Fraction(-minor,(6+i)*(6+j))
        hmoment = rising(1,k)*Fraction(11,2)/rising(Fraction(13,2),k+1)
        vmoment = rising(Fraction(1,2),l)/rising(1,l)
        result += Fraction(int(coef.numerator),int(coef.denominator))*path*hmoment*vmoment
    return result


if __name__ == "__main__":
    for order in ((0,1,3,2,4),(0,1,2,3,4)):
        cases,determinant = check_full_word_polynomial(order)
        expected = (
            {case:Fraction(k,44464437108000) for case,k in zip(
                [(1,0),(11,0),(1,11),(0,11)],[-65558839,-11982421,85509601,-27353656])}
            if order[2] == 3 else dict.fromkeys(cases,Fraction(0)))
        means = {case:mean(order,case) for case in cases}
        assert means == expected
        total = mean(order,"sum")
        assert total == sum(means.values())
        print("Signed oriented means:",means,"summed kernel/pi:",total*determinant/158400)
