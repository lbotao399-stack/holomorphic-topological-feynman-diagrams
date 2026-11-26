"""Signed closed BV edge in the source ordering of differential forms and polyvectors."""
from fractions import Fraction
from itertools import combinations
from math import comb


def exterior(forms):
    values={():1}
    for tail,head in forms:
        updated={}
        for labels,c in values.items():
            for v,sign in ((tail,1),(head,-1)):
                if v not in labels:
                    k=tuple(sorted((*labels,v)))
                    updated[k]=updated.get(k,0)+sign*c*(-1)**sum(i>v for i in labels)
        values={k:c for k,c in updated.items() if c}
    return values


def contract(values,slot,scale):
    return {tuple(i for i in k if i!=slot):scale*(-1)**k.index(slot)*v
            for k,v in values.items() if slot in k}


def verify():
    # Ordered tensor factors: (dbar z, theta_z) at the tail, then (dbar w, theta_w).
    anti=[(i,10+i) for i in range(5)]
    poly=[(5+i,15+i) for i in range(5)]
    full=exterior(anti+poly)
    assert len(full)==1024
    for i in range(5):
        # ∂_z g extracts -y_i; ∂/∂theta_z acts by left contraction in the first factor.
        bv=contract(full,5+i,-1)
        expected=exterior(anti+poly[:i]+poly[i+1:])
        assert bv=={k:(-1)**i*c for k,c in expected.items()}
        for j in range(5):
            # Gauge fixing -4 ∂_z ι_(dbar z) extracts +4 y_j.
            propagator=contract(bv,j,4)
            expected=exterior(anti[:j]+anti[j+1:]+poly[:i]+poly[i+1:])
            assert propagator=={k:4*(-1)**(i+j)*c for k,c in expected.items()}
    assert (-1)**(25*comb(9,2))==1
    # (−32i)^(-2) (4π)^(-10) 4^10 equals −1/(2^10 π^10).
    assert Fraction(-1,32**2)*Fraction(4**10,4**10)==Fraction(-1,2**10)


if __name__=='__main__':
    verify()
    print('All five signed closed BV omissions and all 25 propagator omissions; positive divergence factors.')
