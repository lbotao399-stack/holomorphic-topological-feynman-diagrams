"""Concrete theories (real basis [t_a,t_b] = f_ab^c t_c)."""
from fractions import Fraction as Fr
from ht3d import Theory

def eps(a, b, c):
    if len({a, b, c}) < 3:
        return 0
    return 1 if (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] else -1

F_SU2 = {(a, b, c): eps(a, b, c) for a in range(3) for b in range(3) for c in range(3) if eps(a, b, c)}

def adj_rho():
    # rho(t_a)^b_c = f_ac^b  ([t_a, t_c] = f_ac^b t_b)
    return [{(b, c): eps(a, c, b) for b in range(3) for c in range(3) if eps(a, c, b)} for a in range(3)]

def fund_rho():
    # su(2) fundamental, real basis t_a = -(i/2) sigma_a is complex; use the real form sl(2) instead:
    # basis e,f,h with complex structure constants avoided -> use so(3)=su(2) via adjoint only.
    raise NotImplementedError

def su2_adj(nchi=1, W=None, kappa=0, **kw):
    mults = [(3, adj_rho(), 0) for _ in range(nchi)]
    return Theory(3, F_SU2, mults, W=W, kappa=kappa, name="SU(2)+%d adj" % nchi, **kw)

def u1(charges, W=None, kappa=0, **kw):
    mults = [(1, [{(0, 0): Fr(q)}], 0) for q in charges]
    return Theory(1, {}, mults, W=W, kappa=kappa, name="U(1) charges %s" % (charges,), **kw)

def W_N8():
    # W = eps_abc X^a Y^b Z^c  (= tr X[Y,Z] up to normalization), multiplets 0,1,2 = X,Y,Z
    return {((0, a), (1, b), (2, c)): eps(a, b, c) for a in range(3) for b in range(3) for c in range(3) if eps(a, b, c)}

def W_mass_adj(m=1):
    # (m/2) phi.phi for one adjoint
    return {((0, a), (0, a)): Fr(m, 2) for a in range(3)}

def gl_basis(N):
    return [(i, j) for i in range(N) for j in range(N)]

def F_gl(N):
    """[E_ij, E_kl] = delta_jk E_il - delta_li E_kj"""
    bas = gl_basis(N); ix = {b: n for n, b in enumerate(bas)}
    f = {}
    for a, (i, j) in enumerate(bas):
        for b, (k, l) in enumerate(bas):
            if j == k:
                c = ix[(i, l)]; f[(a, b, c)] = f.get((a, b, c), 0) + 1
            if l == i:
                c = ix[(k, j)]; f[(a, b, c)] = f.get((a, b, c), 0) - 1
    return {k: v for k, v in f.items() if v}

def fund_gl(N, anti=False):
    bas = gl_basis(N)
    rho = []
    for (i, j) in bas:
        # fundamental: E_ij e_j = e_i  -> rho^i_j = 1 ; antifundamental: -E^T: rho^j_i = -1
        rho.append({(j, i): -1} if anti else {(i, j): 1})
    return rho

def gl_theory(N, nf, na, W=None, **kw):
    mults = [(N, fund_gl(N), 0) for _ in range(nf)] + [(N, fund_gl(N, True), 0) for _ in range(na)]
    return Theory(N * N, F_gl(N), mults, W=W, kappa=0, name="GL(%d) %dF+%dAF" % (N, nf, na), **kw)
