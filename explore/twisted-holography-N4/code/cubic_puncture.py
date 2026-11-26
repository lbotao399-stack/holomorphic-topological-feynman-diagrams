"""Čech primitives of the planar loop scalar on the unfluxed punctured cubic."""

from itertools import combinations
import sympy as sp


def add(*terms):
    keys = set().union(*(term.keys() for term in terms))
    return {key: value for key in keys if
            (value := sp.factor(sum(term.get(key, 0) for term in terms))) != 0}


def times(coefficient, field):
    return {key: coefficient*value for key,value in field.items()}


def contract(field, apply):
    result = {}
    for indices, coefficient in field.items():
        for slot,index in enumerate(indices):
            key = indices[:slot]+indices[slot+1:]
            result[key] = result.get(key,0)+(-1)**slot*apply(index,coefficient)
    return add(result)


def gradient_charts(parameter):
    # All fields carry the same tangent wave exp(z1+z2). Only the first
    # tangent direction is needed in the lift; encode this direction by zero.
    w = sp.symbols("x y z")
    rho = sp.Symbol("rho")
    B = sum(c**3 for c in w)/6+parameter*sp.prod(w)
    gradient = [None, *(B.diff(c) for c in w)]

    def divergence(field):
        return contract(field, lambda i,c: c if i==0 else c.diff(w[i-1]))

    def differential(field):
        assert not divergence(field)
        return contract(field,lambda i,c: 0 if i==0 else rho*gradient[i]*c)

    def lift(indices):
        field = {indices: 1/sp.prod(gradient[i] for i in indices)}
        correction = {(0,*key):-c for key,c in divergence(field).items()}
        return add(field,correction)

    cochain = {}
    for n,coefficient in ((1,rho**2),(2,rho),(3,-1)):
        for indices in combinations(range(1,4),n):
            cochain[indices] = times(coefficient,lift(indices))
            faces = [times((-1)**k,cochain[indices[:k]+indices[k+1:]])
                     for k in range(n) if n>1]
            total = add(*faces,times((-1)**(n-1),differential(cochain[indices])))
            assert total == ({():rho**3} if n==1 else {})


def flux_primitive():
    # Check the identities as actual angular coefficient functions. The
    # antiholomorphic coordinates are independent of the holomorphic ones.
    x,y,z,b = sp.symbols("x y z b")
    normal = (x,y,z)
    bars = sp.symbols("X Y Z")
    radius = sum(n*q for n,q in zip(normal,bars))
    r = [q/radius for q in bars]
    angular = [[q.diff(c) for c in bars] for q in r]
    one_forms = [[r[(i+1)%3]*angular[(i+2)%3][j]-r[(i+2)%3]*angular[(i+1)%3][j]
                  for j in range(3)] for i in range(3)]
    primitive = [sp.cancel((x*bx).diff(x,2)/2-2*bx.diff(x)+b*(z*by).diff(x,2))
                 for bx,by in zip(one_forms[0],one_forms[1])]
    gradient = x*x/2+b*y*z
    for j,k in combinations(range(3),2):
        flux = sp.cancel(2*sum(r[i]*(angular[(i+1)%3][j]*angular[(i+2)%3][k]
                              -angular[(i+1)%3][k]*angular[(i+2)%3][j]) for i in range(3)))
        dolbeault = primitive[k].diff(bars[j])-primitive[j].diff(bars[k])
        assert sp.cancel(dolbeault-gradient*flux.diff(x,2)+flux) == 0


def full_flux_lift(quadratic_axis=False):
    import projective_flux as pv
    x,y,z = pv.w
    coefficient = pv.QQ(1,3)
    rho = 3
    axis = 2 if quadratic_axis else 0
    def normal_derivative(field):
        return pv.derivative(field,pv.w[axis])
    if quadratic_axis:
        # Only the z gradient enters on these polarizations. The x,y
        # potential and the current scalar graph coefficient can vary freely.
        background = rho*(z**2/2+x**4+y**3)
        T = pv.scale(normal_derivative(pv.U),-1)
        k = pv.scale(normal_derivative(pv.b[2]),-1)
        scalar_graph = 1+x**2+y**3
    else:
        background = rho*((x**3+y**3+z**3)/6+coefficient*x*y*z)
        T = normal_derivative(normal_derivative(pv.U))
        k = pv.add(pv.scale(normal_derivative(normal_derivative(pv.scale(pv.b[0],x))),pv.QQ(1,2)),
                   pv.scale(normal_derivative(pv.b[0]),-2),
                   pv.scale(normal_derivative(normal_derivative(pv.scale(pv.b[1],z))),coefficient))
        scalar_graph = 1
    def lift(f,vector):
        bivector = pv.add(pv.wedge(T,pv.wedge(pv.xi[2+axis],vector)),
                          pv.scale(pv.wedge(normal_derivative(T),pv.P),-f))
        return pv.add({0:rho*f},pv.scale(bivector,-pv.N),
                      pv.scale(pv.wedge(k,vector),rho*pv.N))

    for f in (pv.z1,pv.z1*pv.z2+pv.z2**3):
        field = lift(f,pv.hamiltonian(f))
        pv.check_zero(pv.divergence(field))
        pv.check_zero(pv.add(pv.differential(field),
                             pv.divergence(pv.scale(field,background))))

    # Now z1,z2 represent exp(z1), exp(z2) in the coefficient ring.
    # Normal derivatives are unchanged. Tangent derivatives act by Euler.
    def wave_divergence(field):
        return pv.add(*(pv.scale(pv.derivative(pv.contract(field,j+3),w),w if j<2 else 1)
                        for j,w in enumerate(pv.x)))

    first = lift(pv.z1,pv.scale(pv.xi[1],pv.z1))
    product = pv.wedge(first,pv.scale(pv.xi[0],pv.z2))
    output_wave = pv.z1*pv.z2
    output = pv.scale(lift(output_wave,wave_divergence(pv.scale(pv.P,output_wave))),scalar_graph)
    pv.check_zero(pv.wedge(product,pv.F))
    pv.check_zero(pv.add(wave_divergence(pv.scale(product,scalar_graph)),pv.scale(output,-1)))
    pv.check_zero(pv.add(pv.bar(output),wave_divergence(pv.scale(output,background)),
                         pv.scale(wave_divergence(pv.wedge(pv.F,output)),pv.N)))


def reduced_homotopy():
    # A local frame with dB=(1,0,0). Keep rho and N formal and all three
    # antiholomorphic directions, including a nonzero third-degree remainder.
    import projective_flux as pv
    rho = pv.z1
    e = pv.xi[1]
    H = pv.add(pv.xi[2],pv.scale(pv.xi[3],pv.v[2]),pv.scale(pv.xi[4],pv.v[1]))
    U = pv.wedge(pv.anti[0],pv.anti[1])
    k = pv.scale(pv.anti[2],pv.v[0])
    T = pv.add(U,pv.bar(k))
    lift = pv.add({0:rho},pv.scale(pv.wedge(T,pv.wedge(pv.xi[2],e)),-pv.N),
                  pv.scale(pv.wedge(k,e),rho*pv.N))

    def differential(field,flux):
        return pv.add(pv.bar(field),pv.scale(pv.contract(field,5),rho),
                      pv.scale(pv.wedge(pv.wedge(U,e),field),flux))

    def primitive(field,flux,degrees):
        result = {}
        residual = field
        for q in degrees:
            lowest = {m:c for m,c in residual.items() if (m&7).bit_count()==q}
            assert lowest
            term = {m:c.exquo(rho) for m,c in pv.wedge(H,lowest).items()}
            result = pv.add(result,term)
            residual = pv.add(residual,pv.scale(differential(term,flux),-1))
        assert not residual
        return result

    assert not differential(lift,pv.N)
    bare = primitive({0:rho**3},0,range(3))
    remaining = pv.add(pv.scale(lift,rho**2),pv.scale(differential(bare,pv.N),-1))
    correction = primitive(remaining,pv.N,range(1,4))
    assert differential(pv.add(bare,correction),pv.N) == pv.scale(lift,rho**2)


if __name__ == "__main__":
    gradient_charts(0)
    gradient_charts(sp.Rational(1,3))
    print("All seven Čech entries give rho^3 times the tangent wave on both smooth cubics.")
    flux_primitive()
    print("The Hesse flux primitive reproduces all three angular two-form coefficients.")
    full_flux_lift()
    full_flux_lift(quadratic_axis=True)
    print("Full cubic-plus-flux differential and free divergence vanish for both tangent tests.")
    reduced_homotopy()
    print("Formal flux homotopy removes rho^2 times the closed lift through form degree three.")
