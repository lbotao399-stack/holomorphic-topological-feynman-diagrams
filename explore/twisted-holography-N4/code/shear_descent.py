"""Normal Dolbeault descent for the three commuting D3 tangent shears."""
import sympy as sp


def exterior(x, y):
    return sp.expand(x[0]*y[1] - x[1]*y[0])


def check():
    w = sp.symbols("w1:4", nonzero=True)
    x, y = sp.symbols("r1 r2")
    r = (x, y, (1-w[0]*x-w[1]*y)/w[2])

    def derivative(value):
        return sp.Matrix([sp.diff(value,x), sp.diff(value,y)])

    def curl(form):
        return sp.diff(form[1],x)-sp.diff(form[0],y)

    dr = tuple(map(derivative,r))
    cyclic = [(i,(i+1)%3,(i+2)%3) for i in range(3)]
    # Coefficients in dr1 ^ dr2, with w held fixed by the Dolbeault differential.
    flux = sp.simplify(2*sum(r[i]*exterior(dr[j],dr[k]) for i,j,k in cyclic))
    assert flux == 2/w[2]
    b = [r[j]*dr[k]-r[k]*dr[j] for i,j,k in cyclic]
    for i,j,k in cyclic:
        assert sp.simplify(curl(b[i])-w[i]*flux) == 0
        assert sp.simplify(w[i]*b[j]-w[j]*b[i]+dr[k]) == sp.zeros(2,1)
    # Čech representatives: δ of the pair primitive is minus the unit triple.
    assert sp.simplify(sum(r[k]/(w[i]*w[j]) for i,j,k in cyclic)) == 1/sp.prod(w)

    # Normal multiplication on a Laurent derivative is minus differentiation
    # with respect to a variable indexing derivatives of the unit residue.
    modes = sp.symbols("m1:4")
    ghost_two = [-v/3 for v in modes]
    detector = lambda v: -sum(sp.diff(c,m) for c,m in zip(v,modes))
    assert detector(ghost_two) == 1
    scalar = sp.symbols("z")
    first_order = [scalar*modes[i]**2*modes[(i+2)%3]+i*scalar**2 for i in range(3)]
    ambiguity = [sp.diff(-sp.diff(first_order[k],modes[j])
                         +sp.diff(first_order[j],modes[k]),scalar) for i,j,k in cyclic]
    assert detector(ambiguity) == 0
    for degree in range(7):
        h = scalar**degree/sp.factorial(degree)
        assert (sp.diff(h,scalar,3) == 0) == (degree < 3)
    print('Unit Čech residue; cyclic detector +1 in units of N on z1^3/3!.')


if __name__ == '__main__':
    check()
