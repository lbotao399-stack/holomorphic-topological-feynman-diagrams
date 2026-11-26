"""Recover a local scalar primitive of a closed integrated one-vector jet cochain."""
from functools import cache
from itertools import combinations

import sympy as sp


DIMENSION = 5
ZERO = (0,)*DIMENSION
indices = {}


@cache
def jet(alpha):
    symbol = sp.Symbol('u_' + '_'.join(map(str,alpha)))
    indices[symbol] = alpha
    return symbol


def shift(alpha,axis,step):
    return tuple(value+step*(j==axis) for j,value in enumerate(alpha))


u = jet(ZERO)
du = tuple(jet(shift(ZERO,i,1)) for i in range(DIMENSION))


def derivative(poly,axis):
    return sp.expand(sum(poly.diff(x)*jet(shift(indices[x],axis,1))
                         for x in poly.free_symbols))


def gradient(poly):
    return [derivative(poly,i) for i in range(DIMENSION)]


def undifferentiate(form):
    """Integrate an exact horizontal one-form through its largest scalar jet."""
    first = form[0]
    primitive = sp.S.Zero
    while first:
        # A total derivative raises the largest jet of its primitive by e_0.
        high = max(first.free_symbols,key=lambda x:(sum(indices[x]),indices[x]))
        alpha = indices[high]
        if not alpha[0]:
            raise ValueError('The supplied one-form has no local horizontal primitive.')
        previous = jet(shift(alpha,0,-1))
        term = sp.integrate(first.diff(high),previous)
        primitive += term
        first = sp.expand(first-derivative(term,0))
    if any(sp.expand(a-b) for a,b in zip(form,gradient(primitive))):
        raise ValueError('The supplied one-form is not horizontally exact.')
    return primitive


def divide_curl_by_du(form):
    """Invert multiplication by du on the closed horizontal two-form dA."""
    curl = {(i,j):sp.expand(derivative(form[j],i)-derivative(form[i],j))
            for i,j in combinations(range(DIMENSION),2)}
    for i,j,k in combinations(range(DIMENSION),3):
        if sp.expand(du[i]*curl[j,k]-du[j]*curl[i,k]+du[k]*curl[i,j]):
            raise ValueError('The one-bivector Ward identity is violated.')
    lower = [sp.S.Zero]*DIMENSION
    for (i,j),coefficient in curl.items():
        for term in sp.Add.make_args(coefficient):
            powers = term.as_powers_dict()
            degree = sum(powers.get(x,0) for x in du)
            # {sum i_i d/du_i, du wedge} = first-jet degree + 5 - form degree.
            lower[j] += term.diff(du[i])/(degree+DIMENSION-2)
            lower[i] -= term.diff(du[j])/(degree+DIMENSION-2)
    return [sp.expand(a) for a in lower]


def variational_derivative(poly):
    total = sp.S.Zero
    for x in poly.free_symbols:
        part = poly.diff(x)
        alpha = indices[x]
        for axis,times in enumerate(alpha):
            for _ in range(times):
                part = derivative(part,axis)
        total += (-1)**sum(alpha)*part
    return sp.expand(total)


def scalar_primitive(form,field_count):
    """Return a polynomial density B and a gradient V with A=du*Euler(B)+dV."""
    if field_count<1:
        raise ValueError('Positive scalar field count is required.')
    chain = [[sp.sympify(v) for v in form]]
    for _ in range(field_count):
        chain.append(divide_curl_by_du(chain[-1]))
    if any(chain[-1]):
        raise ValueError('The zero-field term is outside the positive derivative sector.')
    previous = sp.S.Zero
    for current in reversed(chain[:-1]):
        corrected = [sp.expand(a+x*previous) for a,x in zip(current,du)]
        f = -previous
        previous = undifferentiate(corrected)
    density = sp.expand(u*f/field_count)
    if sp.expand(variational_derivative(density)-f):
        raise ValueError('The two-vector Ward identity is violated.')
    return density,previous


def hexagon_check():
    # For the integrated quartic, Euler(u det Hess(u)/24) = det Hess(u)/4.
    from genus_two_probes import det
    from quartic_unisolvence import multiplicities
    p=[tuple(int(i==j) for j in range(5)) for i in range(5)]+[(1,-2,2,-1,1)]
    k=tuple(-sum(row) for row in zip(*p))
    bare=[sum(det([*(p[i] for i in subset),k])**2
              *(2*k[j]+sum(p[i][j] for i in subset))
              for subset in combinations(range(6),4)) for j in range(5)]
    euler=[6*sum(row[j]*det(p[:i]+p[i+1:])**2 for i,row in enumerate(p))
           for j in range(5)]
    assert bare==[-180,72,-288,0,-180] and euler==[12,12,36,0,12]
    assert [a+b for a,b in zip(bare,euler)]==[84*x for x in k]
    target=multiplicities(20,(7,6,6,6,6))[1]-multiplicities(20,(6,6,6,6,6))[1]
    trace=multiplicities(19,(6,6,6,6,6))[0]
    assert (target,trace,target-trace)==(31268,637,30631)
    print('Quartic hexagon residual is 84 k / 24. Stable marked dimensions: 31268, 637, 30631.')


def check():
    # Four scalar fields and a one-vector cochain with five holomorphic derivatives.
    density = du[0]**2*du[1]**2
    initial_gradient = u**2*du[0]*jet((0,3,0,0,0))
    euler = variational_derivative(density)
    form = [sp.expand(x*euler+b) for x,b in zip(du,gradient(initial_gradient))]
    recovered,gradient_term = scalar_primitive(form,4)
    assert variational_derivative(recovered)==euler
    assert sp.expand(gradient_term-initial_gradient)==0
    bad = [x*du[0] for x in du]
    try:
        scalar_primitive(bad,2)
    except ValueError as e:
        assert 'two-vector' in str(e)
    else:
        raise AssertionError('A non-self-adjoint Frechet derivative must be detected.')
    print('Bivector descent, horizontal gradient, Euler primitive and Helmholtz failure verified.')


if __name__=='__main__':
    check()
    hexagon_check()
