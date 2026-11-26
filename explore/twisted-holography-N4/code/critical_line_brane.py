"""Regular D3 on a moving point of the critical line of (z²+y²h(x))/2."""
import sympy as sp

from quadratic_d3_factorization import exterior_matrices


def punctured_source():
    import projective_flux as pv
    x,y,z=pv.w
    background=(z*z+y*y)/2
    def differential(field):
        return pv.add(pv.differential(field),pv.divergence(pv.scale(field,background)))
    for H in ((1+x)*pv.z1**3*(1+pv.z2),pv.z1*x**4):
        source=H.diff(pv.z1)
        primitive=pv.add(pv.scale(pv.xi[3],y*source/2),pv.scale(pv.xi[0],-H/2))
        pv.check_zero(pv.divergence(primitive))
        lifted=differential(primitive)
        assert lifted[0]==y**2*source/2
        assert any(mask&7 for mask in lifted)
        pv.check_zero(pv.divergence(lifted))
        pv.check_zero(differential(lifted))
    print("Both punctured scalar sources are full flux boundaries on the quadratic background.")


def check():
    x,y,z,t = sp.symbols("x y z t")
    h=sp.Function("h")
    E,H=exterior_matrices()
    unit=sp.eye(8)
    D=(x-t)*H[0]+y*H[1]+z*H[2]+y*h(x)*E[1]/2+z*E[2]/2
    eta=(H[0],H[1]-h(x)*E[1]/2,H[2]-E[2]/2)
    def verify(matrix):
        assert matrix.applyfunc(sp.expand)==sp.zeros(8)
    verify(D**2-(z**2+y**2*h(x))*unit/2)
    assert D.subs({x:t,y:0,z:0})==sp.zeros(8)
    for i,coord in enumerate((x-t,y,z)):
        verify(D*E[i]+E[i]*D-coord*unit)
        verify(D*eta[i]+eta[i]*D)
        for j in range(3):
            square=(0,-h(x),-1)[i] if i==j else 0
            verify(eta[i]*eta[j]+eta[j]*eta[i]-square*unit)
    # Evaluate a non-polynomial analytic graph coefficient on the displaced unit.
    for response in (sp.exp(2*x),h(x).diff(x,2)**2):
        difference=response-response.subs(x,t)
        primitive=E[0]*difference/(x-t)
        verify(D*primitive+primitive*D-difference*unit)
    # Fixed x=t; use the two transverse Koszul factors in y,z.
    fibre=(0,2,4,6)
    parity=sp.diag(1,-1,-1,1)
    restriction=lambda a: a.subs(x,t).extract(fibre,fibre)
    dy,dz=(restriction(D.diff(coord)) for coord in (y,z))
    ordered=restriction(eta[1])*restriction(eta[2])
    numerator=sp.expand((parity*(dy*dz-dz*dy)*ordered).trace()/2)
    assert numerator==-h(t)
    assert sp.cancel(numerator/(h(t)*y*z)) == -1/(y*z)
    print("Moving D3 identity survives; odd squares are 0, -h(t)/2, -1/2.")
    print("The ordered transverse disk returns minus the scalar on each massive fibre.")


def sheared_transverse_disk():
    y,z,t,u,c = sp.symbols('y z t u c')
    E,H = exterior_matrices()
    fibre = (0,2,4,6)
    E,H = [[m.extract(fibre,fibre) for m in ops] for ops in (E,H)]
    Z = z+u*t*y
    D = y*(H[1]+c*t*t*E[1]/2)+Z*(H[2]+E[2]/2)
    unit = sp.eye(4)
    def zero(m):
        return m.applyfunc(sp.expand)==sp.zeros(4)
    assert zero(D**2-(Z**2+c*t*t*y*y)*unit/2)
    eta_z = H[2]-E[2]/2
    eta_y = H[1]-c*t*t*E[1]/2+u*t*eta_z
    assert zero(D*eta_y+eta_y*D) and zero(D*eta_z+eta_z*D)
    assert zero(2*eta_y**2+(c+u**2)*t*t*unit)
    assert zero(eta_y*eta_z+eta_z*eta_y+u*t*unit)
    dy,dz = D.diff(y),D.diff(z)
    parity = sp.diag(1,-1,-1,1)
    numerator = (parity*(dy*dz-dz*dy)*eta_y*eta_z).trace()/2
    jacobian = sp.hessian((Z**2+c*t*t*y*y)/2,(y,z)).det()
    assert jacobian.expand()==c*t*t and sp.cancel(numerator/jacobian)==-1
    print("On a sheared massive fibre the ordered disk remains minus the scalar.")


if __name__=="__main__":
    check()
    sheared_transverse_disk()
    punctured_source()
