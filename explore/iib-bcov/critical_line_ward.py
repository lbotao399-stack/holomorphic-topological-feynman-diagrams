"""Normal reparametrization Ward graphs on the planar critical line."""
from fractions import Fraction as F
from itertools import combinations, product

from bulk_eight import BULK, SHAPES, potentials
from marked_boundaries import quartics
from planar_cubic_primary import HARD


def theta_variation(frequency, first=1, second=1):
    # The normal vector is [v^1]exp(ux+vy)(v d/dx - u d/dy).
    # The three y slopes have joint degree five. Central differences extract
    # precisely the two quadratic donor jets and the first vector jet.
    u=frequency
    tang=(0,1,0,0,0)
    z=(0,0,0,0,1)
    totals=[F(0) for _ in SHAPES]
    for s,t in product((-1,0,1),repeat=2):
        # Simultaneously reflecting all y slopes doubles the v=1 branch of
        # the odd difference. Two repeated z and two repeated mixed donors.
        # Return the symmetric polarization (one half of the mixed Taylor coefficient).
        weight=F(1,4)
        weight*=F(1,2) if s else -1
        weight*=F(1,2) if t else -1
        fields=[tang,(0,0,u,1,0),HARD[0],z,z,(0,0,first,s,0),(0,0,second,t,0)]
        total=tuple(map(sum,zip(*fields)))
        points=fields+[tuple(-a for a in total)]
        normal_total=total[2]-u*total[3]
        for k,shape in enumerate(SHAPES):
            phi=potentials(points,shapes=(shape,))
            root=phi[-1]
            tangent=phi[0][0]-root[0]
            normal=phi[1][2]-root[2]-u*(phi[1][3]-root[3])
            totals[k]+=weight*(total[0]*normal-normal_total*tangent)
    return totals


def exact_input_ward(p,q,L,M):
    """Fusions and original theta diagrams for one specified null-homotopy."""
    def dot(a,b):
        return sum(x*y for x,y in zip(a,b))
    def diff(a,b):
        return tuple(x-y for x,y in zip(a,b))
    channels={}
    shapes=[F(0) for _ in SHAPES]
    tangent=HARD[0]
    # First y derivative of (L1 d_y - s d_1)exp(P+sy)/(2 L1)
    # with the two second derivatives in the remaining y slopes.
    normal=(-1,0,0,L[0],0)
    for t,r in product((-1,0,1),repeat=2):
        weight=F(1,4*L[0])  # normalization of the primitive and the equal z donors
        weight*=F(1,2) if r else -1
        weight*=F(1,2) if t else -1
        z=(0,0,0,0,1)
        fields=[HARD[1],(*L,p,1,0),HARD[0],(*M,q,t,0),(0,0,0,r,0),z,z]
        root=tuple(-sum(col) for col in zip(*fields))
        points=[*fields,root]
        for k,shape in enumerate(SHAPES):
            phi=potentials(points,shapes=(shape,))
            shapes[k]+=weight*(dot(root,normal)*dot(diff(phi[0],phi[-1]),tangent)
                               -dot(root,tangent)*dot(diff(phi[1],phi[-1]),normal))
        for i,j in combinations(range(8),2):
            def insertion(slot,pol):
                if i==slot: return dot(points[j],pol)
                if j==slot: return dot(points[i],pol)
                return 0
            pairing=(-dot(root,normal)*(insertion(0,tangent)-insertion(7,tangent))
                     +dot(root,tangent)*(insertion(1,normal)-insertion(7,normal)))
            if pairing:
                fused=[a for n,a in enumerate(points) if n!=i and n!=j]
                fused.append(tuple(x+y for x,y in zip(points[i],points[j])))
                trace=sum(x*y for x,y in zip(BULK,quartics(fused)))
                channels[i,j]=channels.get((i,j),F(0))+weight*pairing*trace
    assert channels.get((1,4),0) == F(6*(p*q)**2*(1+L[0]+M[0]),62208)
    assert sum(channels.values())+sum(shapes)==0
    return {key:value for key,value in channels.items() if value}, shapes


def exact_input_hexagon(p,q,L,M):
    from bulk_three_hexagon import potentials as hexagon
    total=F(0)
    z=(0,0,0,0,1)
    for t,a,b in product((-1,0,1),repeat=3):
        weight=F(1,24*L[0])  # 2 L1, two equal y donors, three equal z donors
        for s in (t,a,b):
            weight*=F(1,2) if s else -1
        fields=[HARD[1],(*L,p,1,0),HARD[0],(*M,q,t,0),
                (0,0,0,a,0),(0,0,0,b,0),z,z,z]
        root=tuple(-sum(col) for col in zip(*fields))
        phi=hexagon(fields+[root])
        tangent=phi[0][0]-phi[-1][0]
        normal=L[0]*(phi[1][3]-phi[-1][3])-(phi[1][0]-phi[-1][0])
        total+=weight*((L[0]*root[3]-root[0])*tangent-root[0]*normal)
    return total


if __name__=='__main__':
    import sympy as sp
    u=sp.Symbol('u')
    rows=[theta_variation(k) for k in range(6)]
    polynomials=[sp.factor(sp.interpolate(list(enumerate(col)),u)) for col in zip(*rows)]
    expected=[u**2*(82*u+131)/2177280,u**2*(44*u+67)/544320,
              u**2*(2*u+3)/48384,u**2*(3*u+4)/90720]
    assert all(sp.expand(a-b)==0 for a,b in zip(polynomials,expected))
    # Homogeneity in (k,p,q) is five. At k=1 the integer simplex of degree
    # five is unisolvent; symmetry of the equal-background polarization
    # halves that grid. The last point is an independent signed holdout.
    nodes=[(1,p,q) for q in range(6) for p in range(q+1) if p+q<=5]
    for k,p,q in [*nodes,(-1,1,-1),(2,2,3)]:
        value=F(3*k*k*(p*p*q+q*q*p)+2*k**3*(p*p+q*q),20736)
        assert sum(theta_variation(k,p,q)) == value
    print('Original theta polynomials on the normal vector',polynomials,flush=True)
    print('total',sp.factor(sum(polynomials)),flush=True)
    for fixture in ((2,-3,(1,-1),(2,3)),(1,-1,(1,0),(-1,0))):
        channels,shapes=exact_input_ward(*fixture)
        print('Exact-input tree fusions',fixture,channels,flush=True)
        print('Exact-input theta integrals',shapes,flush=True)
        expected=channels[1,4]
        assert channels=={(1,4):expected}
        assert shapes==[-expected*x for x in (F(41,210),F(44,105),F(3,14),F(6,35))]
        assert exact_input_hexagon(*fixture)==0
        print('Exact-input three-loop hexagon: zero.',flush=True)
