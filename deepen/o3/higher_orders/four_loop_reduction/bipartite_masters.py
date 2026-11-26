"""FORCER masters divided by G_0^2 G_1 at t=epsilon+7/2, through the orders of the reconstructed four-point cuts."""
from collections import defaultdict
import sympy as sp

from odd_cut import HALF, bubble, replay


# Columns in replay.MASTERS order. Each key is the exact power of t.
R = sp.Rational
NORMALIZED = (
    {0:R(-426489,25),1:(70337276781935407+32458975203154500*sp.pi**2)/sp.Integer(15292288761750)},
    {0:R(94208,1715),1:R(-258505869452,4990490505),
     2:(-138764493638084616019+31022124567494860800*sp.pi**2)/sp.Integer(121015527116108625)},
    {0:sp.S.One},
    {2:-22287*sp.pi**2/163840},
    {1:R(-2,15015),2:R(-21169,54108054),
     3:-(177249499589+43286443200*sp.pi**2)/sp.Integer(121864864621500)},
)
RESIDUE = sp.pi**2/6994593273507840


def period(row):
    """Integrate a convergent cut when its rational valuations fit the tabulated orders."""
    coefficients = defaultdict(lambda:sp.S.Zero)
    for i,(v,series) in enumerate(zip(row,NORMALIZED,strict=True)):
        laurent = replay.laurent(v,1-min(series))
        if i!=2:  # Column 2 is exactly one at all t.
            assert min(laurent,default=1)>=1-max(series),i
        for n,f in laurent.items():
            for p,g in series.items():
                if n+p<=1:
                    coefficients[n+p] += sp.Rational(f)*g
    assert all(sp.expand(v)==0 for n,v in coefficients.items() if n<1)
    return sp.expand(RESIDUE*coefficients[1])


def beta_series(exponents):
    n,leading,linear,quadratic = 0,sp.S.One,sp.S.Zero,sp.S.Zero
    for j,power in enumerate(exponents):
        v,f,g,h = bubble(j)
        n += power*v
        leading *= f**power
        linear += power*g
        quadratic += power*h
    return {n:leading,n+1:leading*linear,n+2:leading*(quadratic+linear**2/2)}


def verify(directory,executable):
    from bipartite_graphs import calibration_rows, eliminator
    rows = calibration_rows(directory,executable,range(4))
    basis = [beta_series(p) for p in ((0,0,0),(0,-1,1),(-1,0,1))]
    assert sp.expand(beta_series((2,1,0))[-1]/sp.gamma(HALF)-RESIDUE)==0
    for degrees,expected in zip(basis,NORMALIZED[2:]):
        for k,v in expected.items():
            assert sp.expand(degrees[k]-v)==0
    for aux in (2,3):
        eliminate = eliminator([rows[0],rows[1],rows[aux]])
        for k in (0,1):
            unit = [replay.FIELD.zero]*5
            unit[k] = replay.FIELD.one
            # At order t^k the weights of the zero differences are regular.
            # Dividing by the simple pole of G_0^2 G_1 fixes one extra order.
            reduced = eliminate(unit,k)
            last = k+1
            values = defaultdict(lambda:sp.S.Zero)
            for j,(gamma,coefficient) in enumerate(zip(basis,reduced)):
                rational = replay.laurent(coefficient,last-min(gamma))
                if j!=0:
                    assert min(rational,default=last)>=last-max(gamma)
                for n,v in rational.items():
                    for m,b in gamma.items():
                        if n+m<=last:
                            values[n+m] += sp.Rational(v)*b
            for n in set(values)|set(NORMALIZED[k]):
                assert sp.expand(values[n]-NORMALIZED[k].get(n,0))==0
    print('Two calibration systems give the same five normalized master expansions.')


if __name__=='__main__':
    replay.run(verify)
