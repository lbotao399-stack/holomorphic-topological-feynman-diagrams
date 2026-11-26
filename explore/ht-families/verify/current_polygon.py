"""Exact conditional current density on a two-loop theta, including all T6 examples."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from random import Random
import sys
import sympy as s

sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'code'))
from htk import brute
from htk.graph import theta


def geometry(lengths,y):
    p,q,r = lengths
    n = p+q+r
    paths = [range(p),range(p,p+q),range(p+q,n)]
    S = s.Matrix([[int(e in path) for e in range(n)] for path in paths])
    C = s.Matrix([list(S.row(1)-S.row(0)),list(S.row(2)-S.row(0))])
    psi = s.Matrix.vstack(*(C*s.diag(*y.col(i)) for i in range(y.cols)))
    # t=W R is the unique solution with the indicated path sums and zero cycle voltage.
    W = psi.col_join(S).inv()[:,-3:]
    assert S*W == s.eye(3) and psi*W == s.zeros(psi.rows,3)
    return S,C,psi,W


def polygon(W):
    # The normalized path times lie in R_A+R_B+R_C=1 and W R >= 0.
    candidates = set()
    for i,j in combinations(range(W.rows),2):
        v = W.row(i).cross(W.row(j)).T
        if not sum(v):
            continue
        v /= sum(v)
        if min(W*v)>=0:
            candidates.add((v[1],v[2]))
    points = sorted(candidates)
    def hull(chain):
        boundary = []
        for x,y in chain:
            while len(boundary)>=2:
                a,b = boundary[-2]
                c,d = boundary[-1]
                if (c-a)*(y-b)-(d-b)*(x-a)>0:
                    break
                boundary.pop()
            boundary.append((x,y))
        return boundary
    if len(points)<3:
        return []
    cyclic = hull(points)[:-1]+hull(reversed(points))[:-1]
    return [s.Matrix([1-x-y,x,y]) for x,y in cyclic]


def angle_cosines(vertices):
    def gram(u,v):
        a,b,c = u
        A,B,C = v
        return ((a+b)*(A+C)+(A+B)*(a+c))/2 - a*A
    cosines = []
    for i,v in enumerate(vertices):
        vv = gram(v,v)
        if not vv:
            cosines.append(1)
            continue
        u,w = vertices[i-1],vertices[(i+1)%len(vertices)]
        U = u-v*gram(u,v)/vv
        W = w-v*gram(w,v)/vv
        cosines.append(s.factor(-gram(U,W)/s.sqrt(gram(U,U)*gram(W,W))))
    return cosines


def check(lengths):
    h = (sum(lengths)-3)//2
    graph,interior = theta(*lengths)
    rng = Random(340+h+lengths[0])
    sources = [[rng.randint(-3,3) for _ in range(h)] for _ in range(graph.root)]
    while not brute.theta_D(sources,interior):
        sources = [[rng.randint(-3,3) for _ in range(h)] for _ in range(graph.root)]
    times = s.Matrix([rng.randint(1,9) for _ in graph.edges])
    y = s.Matrix(brute.currents(graph,list(map(F,times)),sources))
    S,C,psi,W = geometry(lengths,y)
    periods = S*times
    assert W*periods == times
    b,c = periods[1]/periods[0],periods[2]/periods[0]
    m = C*s.diag(*times)*C.T
    # d ell^i=-M^{-1} Psi^i. Differentiate the two path ratios as well.
    dy = s.Matrix.vstack(*[-m.inv()*psi[2*i:2*i+2,:] for i in range(h)])
    derivatives = s.Matrix.vstack((S.row(1)-b*S.row(0))/periods[0],
                                  (S.row(2)-c*S.row(0))/periods[0],dy)
    jac = derivatives[:,1:].det()/times[0]
    delta = (1+b)*(1+c)-1
    direct = brute.density(graph,list(map(F,times)),sources,T=1)
    assert F(jac**2/(4*delta**3)) == direct.square()
    source_sign = (-1)**((h+1)//2+(h+1)*(lengths[1]+1))*brute.theta_D(sources,interior)
    assert direct.sign() == (1 if source_sign>0 else -1)
    vertices = polygon(W)
    assert len(vertices)>=3 and all(min(W*v)>=0 for v in vertices)
    print(lengths,'current-fibre polygon sides',len(vertices),'original density matched')


def angular_example(h):
    lengths = (1,h+1,h+1)
    eye = s.eye(h)
    zero = s.zeros(1,h)
    path = [zero]+[2*h*eye.row(i) for i in range(h)]
    y = s.Matrix.vstack(s.ones(1,h),*path,*path)
    _,_,_,W = geometry(lengths,y)
    vertices = polygon(W)
    assert set(map(tuple,vertices)) == {(0,1,0),(0,0,1),(s.Rational(1,2),s.Rational(1,4),s.Rational(1,4))}
    assert sorted(angle_cosines(vertices)) == [-s.Rational(1,9),1,1]


if __name__=='__main__':
    for lengths in ((1,2,2),(1,3,3),(2,2,3),(1,6,6),(2,5,6),(3,4,6),(3,5,5),(4,4,5)):
        check(lengths)
    for h in (1,2,5):
        angular_example(h)
    print('Two ideal vertices and cos(angle)=-1/9: area = acos(1/9).')
