"""Four-input tangent Hamiltonian tree on two Laurent states and two D3 shears."""
from math import factorial

import projective_flux as form


def normal_derivatives(value, order, coordinate=2):
    for _ in range(order):
        value = form.derivative(value,form.w[coordinate])
    return value


def lift(hamiltonian, order):
    # The two exterior sources are ordered before the Dolbeault forms.
    derivative = hamiltonian.diff(form.z1)
    angular = (form.U,form.b[0],form.b[1],{0:form.v[2]*form.u})
    tangent = (hamiltonian,derivative,derivative,derivative.diff(form.z1))
    return [form.scale(normal_derivatives(a,order),h) for h,a in zip(tangent,angular)]


def differential(inputs):
    base, first, second, pair = inputs
    def action(value,index):
        return form.scale(form.derivative(value,form.z1),form.w[index])
    return [form.bar(base),
            form.add(action(base,0),form.scale(form.bar(first),-1)),
            form.add(action(base,1),form.scale(form.bar(second),-1)),
            form.add(form.bar(pair),action(second,0),form.scale(action(first,1),-1))]


def poisson(left,right):
    return form.add(form.wedge(form.derivative(left,form.z1),form.derivative(right,form.z2)),
                    form.scale(form.wedge(form.derivative(left,form.z2),form.derivative(right,form.z1)),-1))


def bracket(left,right,left_degree=2):
    result = [{} for _ in range(4)]
    for i,a in enumerate(left):
        for j,b in enumerate(right):
            if i & j:
                continue
            inversions = int(bool(i & 2) and bool(j & 1))
            crossing = (left_degree-i.bit_count())*j.bit_count()
            result[i|j] = form.add(result[i|j],form.scale(poisson(a,b),(-1)**(inversions+crossing)))
    return result


def hamiltonian_coefficient(H,K,a,b):
    def pb(x,y):
        return x.diff(form.z1)*y.diff(form.z2)-x.diff(form.z2)*y.diff(form.z1)
    h,k = H.diff(form.z1),K.diff(form.z1)
    numerator = -factorial(a)*factorial(b+2)*pb(h.diff(form.z1),K)
    numerator -= factorial(a+2)*factorial(b)*pb(H,k.diff(form.z1))
    numerator += 2*factorial(a+1)*factorial(b+1)*pb(h,k)
    return numerator/factorial(a+b+3)


def primitive(hamiltonian, order):
    terms = { (i+1,*rest):c/(i+1) for (i,*rest),c in hamiltonian.items() }
    integral = form.R.from_dict(terms)
    assert integral.diff(form.z1)==hamiltonian
    first = normal_derivatives(normal_derivatives(form.U,order),1,0)
    pair = normal_derivatives(normal_derivatives(form.b[0],order),1,0)
    return [{},{},form.scale(first,-integral),form.scale(pair,hamiltonian)]


def check():
    H = form.z1**3/6
    K = form.z2**2/2
    for first,second in [(H,K),(H*(1+form.z2),K*(1+form.z1**2))]:
        for a,b in ((0,0),(1,1),(1,2),(3,1),(2,3)):
            A,B = lift(first,a),lift(second,b)
            for value in [*differential(A),*differential(B)]:
                form.check_zero(value)
            *lower,pair = bracket(A,B)
            for value in lower:
                form.check_zero(value)
            L = hamiltonian_coefficient(first,second,a,b)
            degree = a+b+1
            form.check_zero(form.add(pair,form.scale(normal_derivatives(form.U,degree),-L)))
            correction = primitive(L,degree)
            *lower,image = differential(correction)
            for value in [*lower,*bracket(correction,A,3),*bracket(correction,B,3),
                          *bracket(correction,correction,3)]:
                form.check_zero(value)
            form.check_zero(form.add(pair,form.scale(image,-1)))
    coefficient = bracket(lift(H,1),lift(K,1))[-1]
    third = normal_derivatives(form.U,3)
    form.check_zero(form.add(coefficient,form.scale(third,form.z2/20)))
    print("Two-shear, two-state tree: +1/20 times d_w3^3 U d/dz1.")


if __name__ == '__main__':
    check()
