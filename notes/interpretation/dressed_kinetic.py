"""Exact symbol calculation on all eight partial-form components in the flagship gauge."""
from sympy import Matrix, SparseMatrix, symbols, zeros, eye


def wedge(i):
    entries = {}
    for mask in range(8):
        if not mask & (1<<i):
            sign = (-1)**((mask & ((1<<i)-1)).bit_count())
            entries[mask+(1<<i),mask] = sign
    return SparseMatrix(8,8,entries)


def zero(m):
    assert all(x.expand()==0 for x in m)


dt,b1,b2 = [wedge(i) for i in range(3)]
t,p1,p2,q1,q2 = symbols('t p1 p2 q1 q2')
D = dt*t+b1*q1+b2*q2
gauge = -dt.T*t-4*(b1.T*p1+b2.T*p2)
shift = b1*p2-b2*p1
heat,interval = symbols('heat interval')
P = gauge*interval
K = eye(8)*heat
zero(D*shift+shift*D)
zero(gauge*shift+shift*gauge)
zero(shift*shift)
zero(gauge*gauge)
zero(P*shift*P)
zero(K*shift*P+P*shift*K)
zero((D+shift)*gauge+gauge*(D+shift)+(t*t+4*(p1*q1+p2*q2))*eye(8))
print('All eight form components satisfy the regulated propagator, heat-edge and Laplacian identities.')
