"""Original five-position current and ordered-simplex checks for the two free prism stars."""
from fractions import Fraction as F
from functools import cache
from math import lcm, prod

from flint import fmpq, fmpq_mat
from quartic_seven_plane_lift_positions import CTX, data as electric, sampled
from quartic_seven_plane_lift import NORMAL,density,fixture
from quartic_seven_plane_double_tensor import frame
from quartic_seven_plane_prism import CORE,STARS,project


@cache
def chamber(exponents,marks):
    """Integrate the position monomial on one ordered simplex in each marked path."""
    d=1
    for edge in set(marks):
        n=0
        for degree,path in zip(exponents,marks):
            if path==edge:
                n+=degree+1
                d*=n
    return F(1,d)


def core_density(times,marks,order,normal):
    u,adj,cycles,intercept,slope,green=electric(CORE,times,marks)
    momenta=[normal[i] for i in order]
    tangential=[(1,-1,0,0,0)[i] for i in order]
    b=sampled(momenta,tangential,marks,cycles,intercept,slope)
    square=sum((adj[i][j]*(1+(i!=j))*b[i]*b[j]
                for i in range(4) for j in range(i,4)),CTX.constant(0))
    s=CTX.gens()
    voltage=CTX.constant(0)
    for i,sign in ((order.index(0),1),(order.index(1),-1)):
        a,b=CORE[marks[i]]
        for j,source in enumerate(momenta):
            charge=sign*source[2]
            if not charge:
                continue
            c,d=CORE[marks[j]]
            voltage+=charge*(int(green[a][c]*u)*(1-s[i])*(1-s[j])
                            +int(green[a][d]*u)*(1-s[i])*s[j]
                            +int(green[b][c]*u)*s[i]*(1-s[j])
                            +int(green[b][d]*u)*s[i]*s[j])
            if marks[i]==marks[j]:
                voltage+=charge*times[marks[i]]*u*(s[min(i,j)]-s[i]*s[j])
    polynomial=square*voltage
    denominator=6**10*u**21*sum(times)**2
    test=[F(1+marks[:i].count(e),1+marks.count(e)) for i,e in enumerate(marks)]
    lengths=[]
    for e,t in enumerate(times):
        fractions=[F(),*[test[i] for i,j in enumerate(marks) if j==e],F(1)]
        lengths.extend(t*(b-a) for a,b in zip(fractions,fractions[1:]))
    scale=lcm(*(x.denominator for x in lengths))
    actual=density(fixture(CORE,marks),[int(scale*x) for x in lengths],normal,orders=[order])[1]
    sample=sum((int(c)*prod(z**int(j) for z,j in zip(test,ex)) for ex,c in polynomial.terms()),F())
    assert sample/denominator==actual/(scale**2*u**4)
    integral=sum((int(c)*chamber(tuple(map(int,e)),marks) for e,c in polynomial.terms()),F())
    return integral*prod(times[i] for i in marks)/denominator,u


def dualize(times):
    z=list(map(F,times))
    jac=F(1)
    for ids in STARS:
        a,b,c=(times[i] for i in ids)
        delta=a*b+b*c+c*a
        new=[F(delta,times[i]) for i in ids]
        for i,l in zip(ids,new):z[i]=l
        jac*=prod(new)/sum(new)**3
    return z,jac


def kirchhoff(edges,times):
    verts=sorted({v for e in edges for v in e})
    b=fmpq_mat([[int(x==v)-int(y==v) for v in verts[1:]] for x,y in edges])
    d=fmpq_mat([[fmpq(1)/fmpq(str(t)) if i==j else 0 for j in range(len(times))] for i,t in enumerate(times)])
    return prod(times)*F(str((b.transpose()*d*b).det()))


def pair_graph():
    edges=list(CORE)
    for vertex,ids in ((0,STARS[0]),(4,STARS[1])):
        neighbors=[next(v for v in CORE[i] if v!=vertex) for i in ids]
        for j,index in enumerate(ids):edges[index]=tuple(neighbors[k] for k in range(3) if k!=j)
    return edges


def verify(marks,order,normal,timings,projector=project):
    transformed=projector(marks,[order],normal)
    edges=pair_graph()
    nonzero=False
    for times in timings:
        calculated,u=core_density(times,marks,order,normal)
        nonzero=nonzero or bool(calculated)
        z,jac=dualize(times)
        u_new=kirchhoff(edges,z)
        assert u_new==u*prod(sum(z[i] for i in ids) for ids in STARS)
        zz=tuple(fmpq(str(x)) for x in z)
        new=F(str(transformed(*zz)))/u_new**13
        assert new==calculated*jac
    assert nonzero


def triangle_order(marks):
    """Exact numerator order on the triangle (2,3,S-5), before parallel contraction."""
    from flint import fmpq_poly
    def scaled(s):
        times=[i*s for i in (2,3,5,7,11,13,17,19,23)]
        times[:3]=[3*(s-5),2*(s-5),6]
        density,u=core_density(times,marks,(0,2,1,3,4),NORMAL)
        return fmpq(str(6*(s-5)*density*u**13))
    # Original path numerator degree is 4*13-9=43. Multiplying the transformed
    # numerator by S^33 gives this path numerator evaluated at the
    # affine times multiplied by the triangle product 6*(S-5).
    values=[scaled(i) for i in range(6,51)]
    x=fmpq_poly([0,1])
    basis=fmpq_poly([1])
    result=fmpq_poly([])
    for k in range(45):
        result+=values[0]*basis
        values=[b-a for a,b in zip(values,values[1:])]
        basis=basis*(x-6-k)/(k+1)
    assert result(59)==scaled(59)
    lowest=next(k for k in range(len(result)) if result[k])
    return lowest-33,result[lowest]


if __name__=='__main__':
    from quartic_seven_plane_prism_paired import MARKS as paired_marks, project as paired_project
    for marks in ((3,3,3,7,7),(3,5,5,5,7),(3,3,3,5,7)):
        verify(marks,tuple(range(5)),NORMAL,[(2,3,5,7,11,13,17,19,23)])
        verify(marks,(0,2,1,3,4),frame(2),[(2,1,3,1,3,2,1,3,2),(2,3,5,7,11,13,17,19,23)])
    for marks in paired_marks:
        verify(marks,(0,2,1,3,4),frame(2),[(2,1,3,1,3,2,1,3,2)],paired_project)
    print('All five prism placements, two source frames, eleven direct chamber and star Jacobian checks.')

    import sys
    if '--occupied-star' in sys.argv:
        for marks,order in (((3,3,3,7,7),2),((0,3,3,7,7),-6),((0,1,3,7,7),-10)):
            measured,leading=triangle_order(marks)
            assert measured==order
            print(marks,'triangle order',measured,'leading coefficient',leading)
