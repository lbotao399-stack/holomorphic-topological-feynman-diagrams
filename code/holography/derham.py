"""Polynomial de Rham complex of the superspace C^{2|3}: z1, z2 even; theta1..3 odd; dz odd; dtheta even.

Form key (m, Sz, J, a):  the monomial  z^m dz^{Sz} theta^J dtheta^a  (in this order),
    m = (m1, m2) >= 0, Sz sorted subset of {0,1}, J sorted subset of {0,1,2}, a = (a1,a2,a3) >= 0.
Weight beta = (m1 + [0 in Sz], m2 + [1 in Sz], [I in J] + a_I);  degree p = |Sz| + |a|;
parity = |Sz| + |J| mod 2.
d (odd derivation, d z = dz, d theta = dtheta):
    d(z^m dz^S theta^J dtheta^a) = sum_{i not in S} m_i (-1)^{#(j in S, j<i)} z^{m-e_i} dz^{S+i} theta^J dtheta^a
                                  + (-1)^{|S|} sum_{I in J} (-1)^{pos(I,J)} z^m dz^S theta^{J-I} dtheta^{a+e_I}.
Bulk dictionary (A2 -> Omega, chain isomorphism, checked in tests):
    Phi(z^m w^(-a-1) xi_z^{Sz} xi_w^J) = eps(Sz) (-1)^{|a|} / a!  z^m dz^{{0,1}-Sz} theta^J dtheta^a,
    eps({0,1}) = 1, eps({1}) = 1, eps({0}) = -1, eps({}) = -1;     Phi(F) = 1,  Phi o del = d o Phi.
"""
from fractions import Fraction
from itertools import combinations
from math import factorial

from .linalg import rank, vadd

SZ_SETS = [(), (0,), (1,), (0, 1)]
J_SETS = [J for r in range(4) for J in combinations(range(3), r)]


def form_weight(key):
    m, Sz, J, a = key
    return (m[0] + (0 in Sz), m[1] + (1 in Sz), (0 in J) + a[0], (1 in J) + a[1], (2 in J) + a[2])


def form_degree(key):
    return len(key[1]) + sum(key[3])


def form_parity(key):
    return (len(key[1]) + len(key[2])) % 2


def forms(beta, p=None):
    out = []
    for Sz in SZ_SETS:
        m = (beta[0] - (0 in Sz), beta[1] - (1 in Sz))
        if min(m) < 0:
            continue
        for J in J_SETS:
            a = tuple(beta[2 + I] - (I in J) for I in range(3))
            if min(a) < 0:
                continue
            key = (m, Sz, J, a)
            if p is None or form_degree(key) == p:
                out.append(key)
    return out


def d_form(v):
    out = {}
    for (m, S, J, a), c in v.items():
        for i in range(2):
            if i in S or m[i] == 0:
                continue
            sgn = (-1) ** sum(1 for j in S if j < i)
            m2 = list(m)
            m2[i] -= 1
            vadd(out, {(tuple(m2), tuple(sorted(S + (i,))), J, a): Fraction(sgn * m[i])}, c)
        for pos, I in enumerate(J):
            sgn = (-1) ** len(S) * (-1) ** pos
            a2 = list(a)
            a2[I] += 1
            vadd(out, {(m, S, tuple(x for x in J if x != I), tuple(a2)): Fraction(sgn)}, c)
    return out


def B_dim(beta, p):
    """dim of exact p-forms at weight beta (rank of d: Omega^{p-1} -> Omega^p)."""
    if p <= 0:
        return 0
    return rank([d_form({k: 1}) for k in forms(beta, p - 1)])


def B_list(beta):
    pmax = sum(beta) + 1
    return [B_dim(beta, p) for p in range(pmax + 1)]


EPS = {(0, 1): 1, (1,): 1, (0,): -1, (): -1}


def Phi(v):
    """Bulk A2 element {(2, e, S): c} -> form."""
    out = {}
    for (typ, e, S), c in v.items():
        assert typ == 2
        m = e[:2]
        a = tuple(-x - 1 for x in e[2:])
        Sz = tuple(i for i in S if i < 2)
        J = tuple(i - 2 for i in S if i >= 2)
        Szc = tuple(i for i in (0, 1) if i not in Sz)
        coef = Fraction(EPS[Sz] * (-1) ** sum(a), factorial(a[0]) * factorial(a[1]) * factorial(a[2]))
        vadd(out, {(m, Szc, J, a): coef}, c)
    return out
