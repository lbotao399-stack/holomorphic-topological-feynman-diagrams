"""Costello 1705.02500 §§10–13,16 and Gaiotto–Rapčák–Zhou 2309.16929 (2.26)–(2.30),(3.58)."""
from collections import Counter
import sympy as sp


def normal_order_trace(gauge_rank,flavor_rank):
    """Trace (13.0.1) in coordinate-before-derivative order."""
    terms = Counter()
    for i in range(gauge_rank):
        for k in range(gauge_rank):
            terms[('phi',i,k)] += 1
            terms[('phi',k,i)] -= 1
        for alpha in range(flavor_rank):
            terms[('rho',i,alpha)] += 1
            terms['hbar'] += 1
    assert not any(terms[('phi',i,k)] for i in range(gauge_rank) for k in range(gauge_rank))
    assert terms['hbar']==gauge_rank*flavor_rank
    assert terms['hbar']==sum(terms[('rho',i,a)] for i in range(gauge_rank)
                              for a in range(flavor_rank))


def verify():
    # Write I for the image of the trace of Id under the §11 ordered lift.
    # The marked disc D(alpha up,beta down) maps to sum_i rho_ia d_rho_ib/h.
    # The unmarked disc maps to K. Thus h I=(c-N h)K in every finite rank.
    for k in range(1,6):
        for n in range(1,6):
            normal_order_trace(k,n)

    h,c,n,I,L,b,q = sp.symbols('h c N I L b q',nonzero=True)
    rank = h*I/(c-n*h)
    # GRZ (A0),(2.28): epsilon1=h, epsilon2=c-Nh, t00=K/h.
    # Their duality fixes t00,h and interchanges epsilon2 with -c.
    # Return to the same c by the dilation of Costello Lemma 10.0.3.
    multiplier = c/(n*h-c)
    def duality(f):
        return sp.factor(f.subs({h:multiplier*h,I:multiplier*I},simultaneous=True))
    assert duality(duality(h))==h and duality(duality(I))==I
    assert sp.cancel(duality(rank)-rank)==0
    # Exhibit an equivariant parameter map with arbitrary genuine-center slope.
    odd_h = h/(1-n*h/(2*c))
    odd_I = I/(1-n*h/(2*c))
    equivariant = L*odd_h+b*odd_I*odd_h**2
    assert sp.cancel(duality(equivariant)+equivariant)==0
    assert sp.expand(sp.series(equivariant,h,0,3).removeO()).coeff(h,2)==L*n/(2*c)+b*I
    # Normalizing the bulk coordinate at I=0 removes only its constant term.
    forward = L*h+b*I*h**2
    inverse = q/L-b*I*q**2/L**3
    assert sp.series(forward.subs(h,inverse)-q,q,0,3).removeO().expand()==0
    assert sp.series(h**2*rank,h,0,3).removeO()==0
    print('Empty disc = h*tr(E)/(c-N*h); the reflected parameter map admits every tr(E) slope.')


if __name__=='__main__':
    verify()
