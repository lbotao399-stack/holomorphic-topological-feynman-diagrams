"""Joint dilation, U(1)^2 characters and diagonal eta contraction at finite heat times."""
from itertools import combinations, product

import sympy as sp


def contraction(form, vector):
    result = {}
    for slots, coeff in form.items():
        for j,slot in enumerate(slots):
            shorter = slots[:j]+slots[j+1:]
            result[shorter] = result.get(shorter,0) + (-1)**j*vector[slot]*coeff
    return {slots:val for slots,val in result.items() if sp.simplify(val) != 0}


def check():
    r,t = sp.symbols("r t",positive=True)
    u,v = sp.symbols("u v",nonzero=True)
    tau,z1,z2,b1,b2 = sp.symbols("tau z1 z2 b1 b2")
    points = {tau:r*tau,z1:r*u*z1,z2:r*v*z2,b1:r*b1/u,b2:r*b2/v,t:r*r*t}
    differential = {0:r,1:r/u,2:r/v,3:r*r}
    kernel = sp.exp(-(tau*tau+z1*b1+z2*b2)/(4*t))/(4*(4*sp.pi*t)**sp.Rational(5,2))
    heat = {(0,1,2):kernel}
    adjoint = {slots:value/t for slots,value in contraction(heat,{0:tau/2,1:b1,2:b2}).items()}
    boundary = {(3,)+slots:-value for slots,value in adjoint.items()}
    character = 1/(r*r*u*v)
    for form,expected in [(heat,character),(adjoint,character/r**2),
                          ({**heat,**boundary},character)]:
        for slots,value in form.items():
            pullback = value.xreplace(points)*sp.prod(differential[i] for i in slots)
            assert sp.simplify(pullback/expected-value) == 0
    # Coordinate and time derivatives are checked on the actual scalar Gaussian.
    for variable,factor in ((tau,r),(z1,r*u),(z2,r*v),(b1,r/u),(b2,r/v),(t,r*r)):
        derivative = sp.diff(kernel,variable)
        assert sp.simplify(derivative.xreplace(points)*r**5*factor-derivative) == 0

    # Difference forms d(x_i-y_i) split into x and y slots. Simultaneous contractions cancel.
    for form in (heat,adjoint):
        expanded = {}
        for slots,value in form.items():
            for ends in product((0,1),repeat=len(slots)):
                labels = tuple(2*i+j for i,j in zip(slots,ends))
                expanded[labels] = (-1)**sum(ends)*value
        for axis in range(3):
            generator = {slot:int(slot//2 == axis) for slot in range(6)}
            assert not contraction(expanded,generator)
    print("Heat, gauge-fixed integrand, full Schwinger form, six derivatives and diagonal contractions pass")

    axes = (0,1,2)
    a,b = sp.symbols("a b",integer=True,nonnegative=True)
    for size in range(4):
        for selected in combinations(axes,size):
            w0,w1,w2 = (0 if i in selected else -1 for i in axes)
            dilation = w0+a+b-2+w1+w2
            spin = (a-1-w1)+(b-1-w2)
            assert sp.expand(dilation-spin) == w0+2*w1+2*w2
            assert (sp.expand(dilation-spin) == 0) == (selected == axes)
    print("All eight classical cohomology summands have the expected joint character")


if __name__ == "__main__":
    check()
