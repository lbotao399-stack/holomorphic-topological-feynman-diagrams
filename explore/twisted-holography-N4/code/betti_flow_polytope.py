"""Local scalar Betti height contributions on stable graphs from exact flow capacity bounds."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, factorial, prod

from flint import fmpq_mat


def fractions(matrix):
    return tuple(F(str(x)) for x in matrix.entries())


def connected(vertices, edges):
    if not vertices:
        return False
    reached = {next(iter(vertices))}
    for _ in range(len(vertices)-1):
        for a,b in edges:
            if a in reached and b in vertices: reached.add(b)
            if b in reached and a in vertices: reached.add(a)
    return reached==vertices


def charge_cuts(edges, charges):
    """Facet inequalities of the graphical zonotope of edge flows with fixed bounds."""
    vertices = set(range(len(charges)))
    cuts = []
    for mask in range(1,1<<len(charges),2):
        side = {i for i in vertices if mask & (1<<i)}
        charge = abs(sum(charges[i] for i in side))
        if charge and connected(side,edges) and connected(vertices-side,edges):
            cut = tuple(i for i,(a,b) in enumerate(edges) if (a in side)!=(b in side))
            cuts.append((cut,charge))
    return tuple(sorted(cuts))


@lru_cache(maxsize=1024)
def simplices(dim, constraints):
    """Exact pulling triangulation using the given rational facet inequalities."""
    halfspaces = [(tuple(-int(i==j) for j in range(dim)),0) for i in range(dim)]
    halfspaces.extend((tuple(int(i in cut) for i in range(dim)),bound) for cut,bound in constraints)
    vertices = set()
    for active in combinations(halfspaces,dim):
        lhs = fmpq_mat([row for row,_ in active])
        if not lhs.det(): continue
        pt = fractions(lhs.solve(fmpq_mat([[rhs] for _,rhs in active])))
        if all(sum(a*b for a,b in zip(row,pt))<=rhs for row,rhs in halfspaces):
            vertices.add(pt)
    vertices = tuple(sorted(vertices))
    faces = [frozenset(i for i,pt in enumerate(vertices) if sum(a*b for a,b in zip(row,pt))==rhs)
             for row,rhs in halfspaces]

    def dim_of(face):
        if not face: return -1
        if len(face)==1: return 0
        p,*rest = [vertices[i] for i in face]
        return fmpq_mat([[str(a-b) for a,b in zip(q,p)] for q in rest]).rank()

    def triangulate(face, degree):
        if len(face)==degree+1:
            yield tuple(sorted(face))
            return
        apex = min(face)
        for sub in sorted({face & facet for facet in faces if apex not in facet},key=sorted):
            if dim_of(sub)==degree-1:
                for cell in triangulate(sub,degree-1):
                    yield (apex,*cell)

    result=[]
    for cell in triangulate(frozenset(range(len(vertices))),dim):
        p,*rest = [vertices[i] for i in cell]
        jac = abs(F(str(fmpq_mat([[str(a-b) for a,b in zip(q,p)] for q in rest]).det())))
        assert jac
        result.append((jac,(p,*rest)))
    return tuple(result)


@lru_cache(maxsize=8192)
def polytope_moment(powers,constraints):
    # Disjoint sets of constrained variables factor as a Cartesian product.
    dim=len(powers)
    component={0}
    for _ in range(dim):
        for cut,_ in constraints:
            if component.intersection(cut): component.update(cut)
    if len(component)<dim:
        total=1
        for block in (sorted(component),sorted(set(range(dim))-component)):
            rows=tuple((tuple(block.index(i) for i in cut),v) for cut,v in constraints if cut[0] in block)
            total*=polytope_moment(tuple(powers[i] for i in block),rows)
        return total
    if len(constraints)==1:
        cut, charge=constraints[0]
        assert len(cut)==dim
        return F(charge**(dim+sum(powers))*prod(map(factorial,powers)),factorial(dim+sum(powers)))
    return sum(simplex_moment(powers,jac,cell) for jac,cell in simplices(dim,constraints))


def simplex_moment(powers,jac,cell):
    """Integral of one ambient monomial over a rational simplex of any dimension."""
    dimension=len(cell)-1
    poly={(0,)*len(cell):1}
    for j,count in enumerate(powers):
        for _ in range(count):
            next_poly={}
            for ex,coef in poly.items():
                for a,vertex in enumerate(cell):
                    if vertex[j]:
                        key=ex[:a]+(ex[a]+1,)+ex[a+1:]
                        next_poly[key]=next_poly.get(key,0)+vertex[j]*coef
            poly=next_poly
    return jac*sum(c*prod(map(factorial,k)) for k,c in poly.items())/factorial(dimension+sum(powers))


@lru_cache(maxsize=16384)
def normal_coefficient(cuts, powers):
    """After X_e = -psi_h-psi_h', powers refer to the remaining branch-psi factors."""
    if not powers: return F(1)
    edges=set(range(len(powers)))
    if not cuts or set().union(*(set(c) for c,_ in cuts))!=edges:
        return F()
    exponents=tuple(2*k+1 for k in powers)
    value=F()
    for size in range(1,len(cuts)+1):
        for bad in combinations(cuts,size):
            if set().union(*(set(c) for c,_ in bad))==edges:
                value+=(-1)**size*polytope_moment(exponents,bad)
    return value/prod(2**k*factorial(k) for k in powers)


def weak(total, count):
    if not count:
        if total==0: yield ()
        return
    for j in range(total+1):
        for tail in weak(total-j,count-1):
            yield (j,*tail)


def check():
    for n in range(1,6):
        edges=((0,1),)*n
        for charge in (1,2):
            cuts=charge_cuts(edges,(charge,-charge))
            for k in range(3):
                for powers in weak(k,n):
                    expected=-F(charge**(2*(n+k))*prod(factorial(2*m+1) for m in powers),
                                factorial(2*(n+k))*prod(2**m*factorial(m) for m in powers))
                    assert normal_coefficient(cuts,powers)==expected
    triangle=((0,1),(0,2),(1,2))
    for charges in ((1,0,-1),(1,2,-3),(4,-3,-1)):
        cuts=charge_cuts(triangle,charges)
        expected=F(sum(x*x for x in charges)**3,960)-F(7*prod(charges)**2,240)
        assert normal_coefficient(cuts,(0,0,0))==expected
    for powers in weak(3,3):
        # The uncoupled box integral cross-checks Jacobians and simplex factorials.
        bounds=(2,1,3)
        box=tuple(((j,),v) for j,v in enumerate(bounds))
        assert polytope_moment(powers,box)==prod(F(v**(p+1),p+1) for p,v in zip(powers,bounds))
    print('Bananas through five edges, triangle sextics and box moments passed.')


if __name__=='__main__':
    check()
