"""Bulk symmetry algebra of twisted IIB on X = C^5 \\ C^2 with flux N F at t^0, and the C^{2|3} calculus.  Exact.

Conventions of bulk_x: A0 key (0, e, S), coordinates (z1,z2,w1,w2,w3) = 0..4, xi_x = d/dx, |a| = #xi mod 2.
g = ker phi_*|t^0 = {alpha in A0 : del alpha = 0, {F, alpha} = 0}  (bulk_x.ker_phi_t0).

Lie superalgebra on A0[1]:   [a, b] = (-1)^{|a|} {a, b},   shifted parity p(a) = |a| + 1.

Dictionary Psi : A0 -> PV(C^{2|3}) = C[z, theta, eta, zeta]   (z even, theta odd, eta_I = d/dtheta_I even, zeta_a = d/dz_a odd),
    the superalgebra isomorphism   z_a -> z_a,  w_I -> -eta_I,  xi_{z_a} -> zeta_a,  xi_{w_I} -> theta_I.
    Polyvector monomial (m, K, c, T) = z^m theta^K eta^c zeta^T (this order);  C^{2|3} polyvector degree k = |T| + |c|;
    C^{2|3} weight gamma = (W_z, -W_w).
    Psi(z^m w^c xi_z^T xi_w^K) = (-1)^{|c| + |T||K|} z^m theta^K eta^c zeta^T.
    Psi o del = Delta o Psi,  Delta = sum_a d_{z_a} d_{zeta_a} - sum_I d_{eta_I} d_{theta_I}: the Berezinian divergence
    (on even X = X^a d_{z_a} + X^I d_{theta_I}:  Delta X = sum_a d_a X^a - sum_I d_{theta_I} X^I).
Forms on C^{2|3}: derham keys (m, Sz, J, a) = z^m dz^Sz theta^J dtheta^a.
Contraction iota : PV(C^{2|3}) -> End Omega(C^{2|3}), superalgebra homomorphism:
    z, theta -> left multiplication,  eta_I -> d/d(dtheta_I),  zeta_a -> d/d(dz_a) (left);
    iota_{z^m theta^K eta^c zeta^T} = z^m theta^K iota_eta^c iota_zeta^T (operators applied right to left).
    Phi(alpha . zeta) = iota_{Psi alpha} Phi(zeta)                                  (test_symmetry)
    Phi({alpha, zeta}) = (L_{Psi alpha} - iota_{Psi del alpha}) Phi(zeta),  L_pi = d iota_pi - (-1)^{|pi|} iota_pi d.
psl(3|3): homogeneous coordinates (Y_0 | Y_1, Y_2 ; Y_3, Y_4, Y_5) = (1 | z1, z2 ; theta1, theta2, theta3) on P^{2|3},
    V_PQ = y_P d/dy_Q (Q != 0),  V_P0 = -y_P E,  E = sum_i y_i d/dy_i;   psi(E_PQ) = -V_PQ  is a Lie superalgebra map
    gl(3|3) -> A0[1] with kernel C.1 and del psi(A) = str(A)  (test_symmetry), so psi(sl(3|3)) = psl(3|3) inside g.
"""
from fractions import Fraction
from itertools import product, combinations

from .linalg import Echelon, kernel, rank, vadd, vscale, clean
from . import bulk_x as bx, derham as dr

ZERO5 = (0, 0, 0, 0, 0)
ONE5 = (1, 1, 1, 1, 1)
UNIT = {(0, ZERO5, ()): Fraction(1)}


# ----------------------------------------------------------------- Lie superalgebra on A0[1]
def apar(a):
    ps = {len(k[2]) % 2 for k in a}
    assert len(ps) <= 1, "inhomogeneous parity"
    return ps.pop() if ps else 0


def spar(a):
    return (apar(a) + 1) % 2


def lie(a, b):
    """[a,b] = (-1)^{|a|}{a,b}; super-antisymmetric and super-Jacobi for the shifted parity (test_symmetry)."""
    if not a or not b:
        return {}
    out = {}
    for ka, ca in a.items():
        s = (-1) ** (len(ka[2]) % 2)
        vadd(out, bx.bv_bracket({ka: ca}, b), s)
    return clean(out)


# ----------------------------------------------------------------- Psi and polyvector data
def psi_key(key):
    """A0 monomial -> (sign, (m, K, c, T)) with K = theta indices 0..2, T = zeta indices 0..1."""
    typ, e, S = key
    assert typ == 0
    T = tuple(x for x in S if x < 2)
    K = tuple(x - 2 for x in S if x >= 2)
    c = e[2:]
    sign = (-1) ** (sum(c) + len(T) * len(K))
    return sign, (e[:2], K, c, T)


def pv_degree(key):
    """C^{2|3} polyvector degree k = #zeta + #eta = |S_z| + w-degree."""
    typ, e, S = key
    return sum(1 for x in S if x < 2) + sum(e[2:])


def c23_weight(W):
    return (W[0], W[1], -W[2], -W[3], -W[4])


def label_of_regular(W):
    """Source label beta = W_r + (1^5) of the dictionary weight rule."""
    return tuple(x + 1 for x in W)


def gamma_of_label(beta):
    """C^{2|3} weight of the regular modes with source label beta."""
    return (beta[0] - 1, beta[1] - 1, 1 - beta[2], 1 - beta[3], 1 - beta[4])


def split_by_degree(a):
    out = {}
    for k, c in a.items():
        out.setdefault(pv_degree(k), {})[k] = c
    return out


def pretty_pv(a):
    """Human-readable C^{2|3} polyvector of a bulk A0 element."""
    names_z = ("z1", "z2")
    names_t = ("t1", "t2", "t3")
    terms = []
    for key, c in sorted(a.items()):
        s, (m, K, cc, T) = psi_key(key)
        c = c * s
        f = []
        for i in range(2):
            if m[i]:
                f.append(names_z[i] + ("^%d" % m[i] if m[i] > 1 else ""))
        for I in K:
            f.append(names_t[I])
        for I in range(3):
            if cc[I]:
                f.append("d" + names_t[I] + ("^%d" % cc[I] if cc[I] > 1 else ""))
        for a_ in T:
            f.append("d" + names_z[a_])
        mon = "*".join(f) if f else "1"
        terms.append("%s%s" % (("+" if c > 0 else "-") + (str(abs(c)) + "*" if abs(c) != 1 else ""), mon))
    return " ".join(terms) if terms else "0"


# ----------------------------------------------------------------- forms on C^{2|3}
def iota_zeta(v, a_):
    out = {}
    for (m, S, J, a), c in v.items():
        if a_ not in S:
            continue
        p = S.index(a_)
        vadd(out, {(m, tuple(x for x in S if x != a_), J, a): Fraction((-1) ** p)}, c)
    return out


def iota_eta(v, I):
    out = {}
    for (m, S, J, a), c in v.items():
        if a[I] == 0:
            continue
        a2 = list(a)
        a2[I] -= 1
        vadd(out, {(m, S, J, tuple(a2)): Fraction(a[I])}, c)
    return out


def mul_theta(v, I):
    out = {}
    for (m, S, J, a), c in v.items():
        if I in J:
            continue
        J2 = tuple(sorted(J + (I,)))
        s = (-1) ** len(S) * (-1) ** J2.index(I)
        vadd(out, {(m, S, J2, a): Fraction(s)}, c)
    return out


def mul_z(v, n):
    return {((m[0] + n[0], m[1] + n[1]), S, J, a): c for (m, S, J, a), c in v.items()}


def iota_pvmono(pvm, omega):
    m, K, c, T = pvm
    v = dict(omega)
    for a_ in reversed(T):
        v = iota_zeta(v, a_)
    for I in range(3):
        for _ in range(c[I]):
            v = iota_eta(v, I)
    for I in reversed(K):
        v = mul_theta(v, I)
    return mul_z(v, m)


def iota(alpha, omega):
    """iota_{Psi alpha} omega."""
    out = {}
    for key, c in alpha.items():
        s, pvm = psi_key(key)
        vadd(out, iota_pvmono(pvm, omega), c * s)
    return out


def phi_inv(omega):
    """Omega(C^{2|3}) -> A2, inverse of derham.Phi."""
    out = {}
    for (m, Sz, J, a), c in omega.items():
        Szc = tuple(i for i in (0, 1) if i not in Sz)
        coef = Fraction(dr.EPS[Szc] * (-1) ** sum(a) * _fact(a))
        e = (m[0], m[1]) + tuple(-x - 1 for x in a)
        S = Szc + tuple(2 + I for I in J)
        vadd(out, {(2, e, S): coef}, c)
    return out


def _fact(a):
    from math import factorial
    r = 1
    for x in a:
        r *= factorial(x)
    return r


def L_pv(alpha, omega):
    """Generalized Lie derivative L_{Psi alpha} = d iota - (-1)^{|alpha|} iota d (alpha homogeneous)."""
    s = (-1) ** apar(alpha)
    out = dr.d_form(iota(alpha, omega))
    vadd(out, iota(alpha, dr.d_form(omega)), -s)
    return out


def bulk_on_forms(alpha, omega):
    """Phi({alpha, Phi^{-1} omega}): the bulk bracket action of alpha in A0 on Laurent modes, in form language."""
    return dr.Phi(bx.bv_bracket(alpha, phi_inv(omega)))


# ----------------------------------------------------------------- forms: product and Lie derivative (independent)
def form_mul(u, v):
    """Product of polynomial forms on C^{2|3} (z, dtheta even; dz, theta odd)."""
    out = {}
    for (m1, S1, J1, a1), c1 in u.items():
        for (m2, S2, J2, a2), c2 in v.items():
            if set(S1) & set(S2) or set(J1) & set(J2):
                continue
            # z^m1 dz^S1 th^J1 dth^a1 . z^m2 dz^S2 th^J2 dth^a2: move dz^S2 left past th^J1
            sgn = (-1) ** (len(J1) * len(S2))
            sgn *= (-1) ** sum(1 for i in S1 for j in S2 if i > j)
            sgn *= (-1) ** sum(1 for i in J1 for j in J2 if i > j)
            key = ((m1[0] + m2[0], m1[1] + m2[1]), tuple(sorted(S1 + S2)), tuple(sorted(J1 + J2)),
                   tuple(x + y for x, y in zip(a1, a2)))
            vadd(out, {key: Fraction(sgn)}, c1 * c2)
    return out


def _gen_list(key):
    """Ordered generator list of a form monomial: ('z',a),('dz',a),('t',I),('dt',I) with multiplicity."""
    m, S, J, a = key
    g = []
    for i in range(2):
        g += [('z', i)] * m[i]
    g += [('dz', i) for i in S]
    g += [('t', I) for I in J]
    for I in range(3):
        g += [('dt', I)] * a[I]
    return g


_GPAR = {'z': 0, 'dz': 1, 't': 1, 'dt': 0}


def _gen_form(g):
    t, i = g
    if t == 'z':
        m = [0, 0]; m[i] = 1
        return {(tuple(m), (), (), (0, 0, 0)): Fraction(1)}
    if t == 'dz':
        return {((0, 0), (i,), (), (0, 0, 0)): Fraction(1)}
    if t == 't':
        return {((0, 0), (), (i,), (0, 0, 0)): Fraction(1)}
    a = [0, 0, 0]; a[i] = 1
    return {((0, 0), (), (), tuple(a)): Fraction(1)}


def lie_derivative(X, omega, parity):
    """Standard Lie derivative of a homogeneous vector field X on C^{2|3} (parity = its parity as a derivation).

    X = {coord: 0-form}, coord 0,1 = z1,z2 and 2,3,4 = theta1..3.  L_X y = X^y, L_X dy = (-1)^{|X|} d(X^y),
    extended as a derivation of parity |X| with Koszul signs over the ordered generators of each monomial.
    """
    out = {}
    for key, c in omega.items():
        gens = _gen_list(key)
        pre = 0
        for idx, g in enumerate(gens):
            t, i = g
            if t in ('z', 't'):
                img = X.get(i if t == 'z' else 2 + i, {})
            else:
                img = dr.d_form(X.get(i if t == 'dz' else 2 + i, {}))
                img = vscale(img, (-1) ** parity)
            if img:
                left = {((0, 0), (), (), (0, 0, 0)): Fraction(1)}
                for g2 in gens[:idx]:
                    left = form_mul(left, _gen_form(g2))
                right = {((0, 0), (), (), (0, 0, 0)): Fraction(1)}
                for g2 in gens[idx + 1:]:
                    right = form_mul(right, _gen_form(g2))
                term = form_mul(form_mul(left, img), right)
                vadd(out, term, c * (-1) ** (parity * pre))
            pre += _GPAR[t]
    return out


def vector_field_of(alpha):
    """For alpha in A0 with Psi(alpha) of polyvector degree 1: X = {coord: 0-form} with Psi(alpha) = sum X^a zeta_a + X^I eta_I,
    read as the vector field sum X^a d/dz_a + X^I d/dtheta_I."""
    X = {}
    for key, c in alpha.items():
        s, (m, K, cc, T) = psi_key(key)
        assert len(T) + sum(cc) == 1
        f = {(m, (), K, (0, 0, 0)): Fraction(c * s)}
        if T:
            # z^m theta^K zeta_a = (z^m theta^K) d/dz_a
            vadd(X.setdefault(T[0], {}), f)
        else:
            I = [i for i in range(3) if cc[i]][0]
            vadd(X.setdefault(2 + I, {}), f)
    return {k: v for k, v in X.items() if v}


# ----------------------------------------------------------------- psl(3|3)
PAR6 = (0, 0, 0, 1, 1, 1)


def _mono(e, S, c=1):
    return {(0, tuple(e), tuple(S)): Fraction(c)}


def y_elem(P):
    if P == 0:
        return dict(UNIT)
    if P in (1, 2):
        e = [0] * 5; e[P - 1] = 1
        return _mono(e, ())
    return _mono(ZERO5, (2 + P - 3,))


def dy_elem(Q):
    if Q in (1, 2):
        return _mono(ZERO5, (Q - 1,))
    e = [0] * 5; e[2 + Q - 3] = 1
    return _mono(e, (), -1)


def euler_elem():
    E = {}
    for i in range(1, 6):
        vadd(E, bx.mul(y_elem(i), dy_elem(i)))
    return E


def V_elem(P, Q):
    if Q != 0:
        return clean(bx.mul(y_elem(P), dy_elem(Q)))
    return clean(vscale(bx.mul(y_elem(P), euler_elem()), -1))


def psi_mat(M):
    """psi(A) = -sum A_PQ V_PQ for a supermatrix A = {(P,Q): coeff}."""
    out = {}
    for (P, Q), c in M.items():
        vadd(out, V_elem(P, Q), -c)
    return clean(out)


def supercomm(A, B):
    """[A,B] = AB - (-1)^{|A||B|} BA for homogeneous supermatrices (dicts)."""
    def par(M):
        ps = {(PAR6[P] + PAR6[Q]) % 2 for (P, Q) in M}
        assert len(ps) <= 1
        return ps.pop() if ps else 0
    out = {}
    for (P, Q), a in A.items():
        for (R, S), b in B.items():
            if Q == R:
                vadd(out, {(P, S): Fraction(1)}, a * b)
            if S == P:
                vadd(out, {(R, Q): Fraction(1)}, -(-1) ** (par(A) * par(B)) * a * b)
    return clean(out)


def psl33_basis():
    """[(name, supermatrix, bulk element)]: 30 root vectors E_PQ (P != Q) and 4 Cartan elements of psl(3|3)."""
    out = []
    for P in range(6):
        for Q in range(6):
            if P != Q:
                M = {(P, Q): Fraction(1)}
                out.append(("E%d%d" % (P, Q), M, psi_mat(M)))
    # sl(3|3) Cartan modulo the identity: H = E_00 - E_11, E_11 - E_22, E_22 + E_33, E_33 - E_44 (E_44 - E_55 is dependent mod 1)
    for (P, Q, s) in ((0, 1, -1), (1, 2, -1), (2, 3, 1), (3, 4, -1)):
        M = {(P, P): Fraction(1), (Q, Q): Fraction(s)}
        out.append(("H%d%d" % (P, Q), M, psi_mat(M)))
    return out


def coords_in(basis_vecs, v):
    """Coordinates of v in span(basis_vecs) (exact), or None."""
    E = Echelon(track=True)
    for i, b in enumerate(basis_vecs):
        E.add(b, label=i)
    row, comb = E.reduce(v)
    if row:
        return None
    return {i: -c for i, c in comb.items()}


def structure_constants(basis_vecs):
    """f[i][j] = coords of [b_i, b_j]."""
    E = Echelon(track=True)
    for i, b in enumerate(basis_vecs):
        indep, _ = E.add(b, label=i)
        assert indep, "basis dependent"
    f = {}
    for i, a in enumerate(basis_vecs):
        for j, b in enumerate(basis_vecs):
            v = lie(a, b)
            row, comb = E.reduce(v)
            assert not row, ("not closed", i, j)
            f[(i, j)] = {k: -c for k, c in comb.items()}
    return f


def killing_form(basis_vecs, f):
    """kappa(i,j) = str(ad_i ad_j) with shifted parities."""
    n = len(basis_vecs)
    par = [spar(b) for b in basis_vecs]
    kap = {}
    for i in range(n):
        for j in range(n):
            s = Fraction(0)
            for k in range(n):
                # (ad_i ad_j)(b_k) = sum_l f[j,k]_l [b_i, b_l] ; coefficient on b_k
                for l, c in f[(j, k)].items():
                    s += (-1) ** par[k] * c * f[(i, l)].get(k, 0)
            kap[(i, j)] = s
    return kap


# ----------------------------------------------------------------- g by weight, Z(A0), excess
def Z_A0(W):
    return bx.ker_del(0, W)[1]


def g_basis(W):
    return bx.ker_phi_t0(W)


def g_degree_dims(W):
    """{k: dim of the divergence-free part of C^{2|3}-polyvector degree k at bulk weight W} and dim g(W)."""
    B = bx.basis(0, W)
    byk = {}
    for b in B:
        byk.setdefault(pv_degree(b), []).append(b)
    dims = {}
    for k, keys in byk.items():
        imgs = [bx.delta({b: 1}) for b in keys]
        dims[k] = len(keys) - rank(imgs)
    return dims


def B1_dim(beta):
    """U(1) letters: dim B^1(beta) = [beta_z >= 0, beta_theta in {0,1}^3, beta != 0]."""
    if min(beta) < 0:
        return 0
    if any(x > 1 for x in beta[2:]):
        return 0
    return 0 if beta == ZERO5 else 1


def excess_formula(beta):
    """dim g(beta - 1^5) - dim coker(beta) = dim B^1(beta) - dim B^1(gamma(beta))  [derived]."""
    return B1_dim(beta) - B1_dim(gamma_of_label(beta))


def in_M1(b):
    return b != ZERO5 and all(x in (0, 1) for x in b[2:]) and min(b[:2]) == 0


# ----------------------------------------------------------------- action on coker phi_*
def coker_data(Wn):
    """(Z2 basis, B1 = im phi_* basis) at Laurent weight Wn."""
    return bx.coker_phi_t0(Wn)


def coker_complement(Wn):
    """Basis of a complement of B^1 in Z(A2)_Wn (representatives of coker phi_*)."""
    Z2, im = coker_data(Wn)
    E = Echelon()
    for v in im:
        E.add(v)
    comp = []
    for z in Z2:
        if E.add(z)[0]:
            comp.append(z)
    return comp, im


def rho_matrix_rank(alpha, Wa, Wn):
    """rank of rho(alpha): coker(Wn) -> coker(Wn + Wa), rho(alpha) = {alpha, -} mod B^1."""
    comp, _ = coker_complement(Wn)
    if not comp:
        return 0
    Wt = tuple(x + y for x, y in zip(Wn, Wa))
    _, imt = coker_data(Wt)
    E = Echelon()
    for v in imt:
        E.add(v)
    r0 = E.rank()
    for z in comp:
        E.add(bx.bv_bracket(alpha, z))
    return E.rank() - r0


# ----------------------------------------------------------------- flat (Koszul) source of the U(1) letter
def degree_split_Z2(Wn):
    """Closed A2 at Wn split by form degree of Phi: {p: basis of closed forms of degree p} (as A2 elements)."""
    B = bx.basis(2, Wn)
    byp = {}
    for b in B:
        f = next(iter(dr.Phi({b: 1})))
        byp.setdefault(dr.form_degree(f), []).append(b)
    out = {}
    for p, keys in byp.items():
        imgs = [bx.delta({b: 1}) for b in keys]
        out[p] = [{keys[i]: c for i, c in comb.items()} for comb in kernel(imgs)]
    return out


def u1_source(beta, P=bx.fischer_res):
    """For a B^1-weight beta with beta_z != (0,0): the line in Z(A0)(beta - 1^5) annihilated by P_A against
    the closed Laurent modes of form degree >= 2 (the flat source of the U(1) letter).  Returns basis list."""
    Wr = tuple(x - 1 for x in beta)
    Wn = (beta[0] - 1, beta[1] - 1) + tuple(-x - 1 for x in beta[2:])
    Z0 = Z_A0(Wr)
    byp = degree_split_Z2(Wn)
    high = [z for p, zs in byp.items() if p >= 2 for z in zs]
    rows = [{j: bx.pair(P, a, z) for j, z in enumerate(high)} for a in Z0]
    ker = kernel(rows)
    out = []
    for comb in ker:
        v = {}
        for i, c in comb.items():
            vadd(v, Z0[i], c)
        out.append(clean(v))
    return out


def alpha0_is_constant(alpha):
    """alpha|_{w=0, xi_z=0} in C  <=>  alpha in C.1 + (w, xi_z)."""
    for (typ, e, S), c in alpha.items():
        if any(x < 2 for x in S) or any(e[2:]):
            continue
        if e[:2] != (0, 0) or S:
            return False
    return True


def bulk_of_pv(m, K, c, T, coef=Fraction(1)):
    """Psi^{-1} of the C^{2|3} polyvector monomial z^m theta^K eta^c zeta^T (K, T sorted tuples)."""
    sign = (-1) ** (sum(c) + len(T) * len(K))
    e = (m[0], m[1]) + tuple(c)
    S = tuple(T) + tuple(2 + I for I in K)
    return {(0, e, S): Fraction(coef) * sign}


def hamiltonian_vf(h_exp, K):
    """theta^K H_h for h = z^{h_exp}, H_h = (d_2 h) d/dz_1 - (d_1 h) d/dz_2 (Hamiltonian for dz1^dz2), as bulk element.
    For h_exp = (0,0) returns theta^K d/dz_1 ^ d/dz_2 (the Poisson bivector)."""
    m1, m2 = h_exp
    out = {}
    if (m1, m2) == (0, 0):
        return bulk_of_pv((0, 0), K, (0, 0, 0), (0, 1))
    if m2:
        vadd(out, bulk_of_pv((m1, m2 - 1), K, (0, 0, 0), (0,), m2))
    if m1:
        vadd(out, bulk_of_pv((m1 - 1, m2), K, (0, 0, 0), (1,), -m1))
    return clean(out)


def u1_formula(beta):
    """Predicted flat source of the U(1) letter d(z^{beta_z} theta^{beta_theta}): theta^{1-beta_theta} H_{z^{beta_z}}."""
    K = tuple(I for I in range(3) if beta[2 + I] == 0)
    return hamiltonian_vf(beta[:2], K)


def proportional(u, v):
    if not u or not v or set(u) != set(v):
        return False
    r = {u[k] / v[k] for k in u}
    return len(r) == 1


def ks_source_z0(beta):
    """beta_z = (0,0), beta a B^1-weight: annihilator in Z(A0)(beta-1^5) of B^{>=2}(beta) under K = P_A o del^{-1}."""
    Wr = tuple(x - 1 for x in beta)
    Wn = (beta[0] - 1, beta[1] - 1) + tuple(-x - 1 for x in beta[2:])
    Z0 = Z_A0(Wr)
    byp = degree_split_Z2(Wn)
    high = [z for p, zs in byp.items() if p >= 2 for z in zs]
    B2 = bx.basis(2, Wn)
    prims = [bx.preimage(z, B2, bx.delta) for z in high]
    assert all(g is not None for g in prims)
    rows = [{j: bx.pair(bx.fischer_res, a, g) for j, g in enumerate(prims)} for a in Z0]
    out = []
    for comb in kernel(rows):
        v = {}
        for i, c in comb.items():
            vadd(v, Z0[i], c)
        out.append(clean(v))
    return out


def projected_pairing_data(beta):
    """Canonical D-compatible pairing P_D(alpha, zeta) = P_A(alpha, zeta_{p>=2}) (beta_z != (0,0)) or K on zeta_{p>=2}
    (beta_z = (0,0)), restricted to g(beta - 1^5) x coker(beta).  coker(beta) is represented by B^{>=2}(beta) (+ F at beta=0).
    Returns dims and rank."""
    Wr = tuple(x - 1 for x in beta)
    Wn = (beta[0] - 1, beta[1] - 1) + tuple(-x - 1 for x in beta[2:])
    g = g_basis(Wr)
    byp = degree_split_Z2(Wn)
    high = [z for p, zs in byp.items() if p >= 2 for z in zs]
    if beta == ZERO5:
        high = [z for zs in byp.values() for z in zs]
    if beta[:2] == (0, 0) and beta != ZERO5:
        B2 = bx.basis(2, Wn)
        high = [bx.preimage(z, B2, bx.delta) for z in high]
    rows = [{j: bx.pair(bx.fischer_res, a, z) for j, z in enumerate(high)} for a in g]
    r = rank(rows)
    return dict(dim_g=len(g), dim_coker=len(high), rank=r, left_rad=len(g) - r, right_rad=len(high) - r)


def residual_trivial(W, Wn_iter, stop_at=0):
    """Subspace of g(W) acting trivially on coker(Wn) for all Wn in Wn_iter (incremental; stops when dim <= stop_at).
    Returns (basis of the residual subspace, number of Wn used)."""
    cur = g_basis(W)
    used = 0
    for Wn in Wn_iter:
        if len(cur) <= stop_at:
            break
        comp, _ = coker_complement(Wn)
        if not comp:
            continue
        Wt = tuple(x + y for x, y in zip(Wn, W))
        if not bx.basis(2, Wt):
            continue
        used += 1
        _, imt = coker_data(Wt)
        E = Echelon()
        for v in imt:
            E.add(v)
        vecs = [dict() for _ in cur]
        for zi, z in enumerate(comp):
            for gi, a in enumerate(cur):
                row, _ = E.reduce(bx.bv_bracket(a, z))
                for k, c in row.items():
                    vecs[gi][(zi,) + k] = c
        ker = kernel(vecs)
        new = []
        for comb in ker:
            v = {}
            for i, c in comb.items():
                vadd(v, cur[i], c)
            new.append(clean(v))
        cur = new
    return cur, used
