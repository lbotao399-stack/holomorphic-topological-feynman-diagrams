"""Bulk complex of twisted IIB on X = C^5 \\ C^2 (Costello-Li 1606.00365 Sec. 9) at the level of
Dolbeault cohomology, with D = t*del + N*{F, -}.

Coordinates x = (z1, z2, w1, w2, w3) <-> index 0..4.  Odd generators xi_x = d/dx.
Monomial key  (typ, e, S):
    typ = 0 : A0 = H^0 = C[z,w] (x) Lambda[xi],          e = (m1, m2, b1, b2, b3), all >= 0
    typ = 2 : A2 = H^2 = C[z] (x) L (x) Lambda[xi],      e = (m1, m2, -a1-1, -a2-1, -a3-1), w-entries <= -1
    S = sorted tuple of xi indices;  the monomial is  z^m w^e  xi_{S[0]} xi_{S[1]} ...
Weight (torus T^5 scaling x):  W = e - ind(S).
Products: A0*A0 -> A0, A0*A2 -> A2 (the C[w]-module L: zero as soon as a w-exponent is >= 0),
A2*A2 = 0 (H^2.H^2 lands in H^4(C^3\\0) = 0).  Functions are even; parity = |S| mod 2 (Dolbeault degree
0 or 2 is even).  d/dx acts on exponents (on L: w^(-a-1) -> -(a+1) w^(-a-2)), d/dxi_x acts from the left.
del = sum_x d/dx d/dxi_x   (divergence for Omega = dz1 dz2 dw1 dw2 dw3).
BV bracket:  {a,b} = del(ab) - del(a) b - (-1)^{|a|} a del(b).
Explicit form (checked in tests):  {a,b} = sum_x [ (d_xi_x a)(d_x b) + (-1)^{|a|} (d_x a)(d_xi_x b) ].
F = (w1 w2 w3)^(-1) xi_z1 xi_z2  (key (2, (0,0,-1,-1,-1), (0,1))).

Full complex T = A0[[t]] (+) A2[[t]], key (k, typ, e, S) for t^k * monomial,
    D = t del + N {F, -}.
Grading preserved by D (finite pieces):
    What = W - (1,1,1,1,1) on A0,  What = W on A2;      h = k - 1 on A0,  h = k on A2;
    D maps T(What, h) -> T(What, h+1);   T(What, h) = t^h A2_What (+) t^(h+1) A0_(What+1).
"""
from fractions import Fraction
from itertools import combinations
from math import factorial

from .linalg import Echelon, kernel, rank, rel_rank, vadd, clean

ZERO5 = (0, 0, 0, 0, 0)
ONE5 = (1, 1, 1, 1, 1)
F_KEY = (2, (0, 0, -1, -1, -1), (0, 1))
F = {F_KEY: Fraction(1)}
XI_TOP = {(0, ZERO5, (0, 1, 2, 3, 4)): Fraction(1)}


# ----------------------------------------------------------------- monomials
def xi_mul(S, T):
    """xi^S xi^T = sign * xi^(S u T); None if they share an index."""
    if set(S) & set(T):
        return None
    inv = sum(1 for i in S for j in T if i > j)
    return (-1) ** inv, tuple(sorted(S + T))


def valid(typ, e):
    if typ == 0:
        return min(e) >= 0
    if typ == 2:
        return e[0] >= 0 and e[1] >= 0 and max(e[2:]) <= -1
    return False


def weight(key):
    typ, e, S = key
    return tuple(e[i] - (1 if i in S else 0) for i in range(5))


def parity(key):
    return len(key[2]) % 2


def mono_mul(k1, k2):
    t1, e1, S1 = k1
    t2, e2, S2 = k2
    typ = t1 + t2
    if typ > 2:
        return None
    e = tuple(a + b for a, b in zip(e1, e2))
    if not valid(typ, e):
        return None
    r = xi_mul(S1, S2)
    if r is None:
        return None
    return r[0], (typ, e, r[1])


def mul(a, b):
    out = {}
    for k1, c1 in a.items():
        for k2, c2 in b.items():
            r = mono_mul(k1, k2)
            if r is not None:
                s, k = r
                v = out.get(k, 0) + s * c1 * c2
                if v:
                    out[k] = v
                else:
                    out.pop(k)
    return out


def d_x(a, x):
    out = {}
    for (typ, e, S), c in a.items():
        if e[x] == 0:
            continue
        e2 = list(e)
        e2[x] -= 1
        vadd(out, {(typ, tuple(e2), S): Fraction(e[x])}, c)
    return out


def d_xi(a, x):
    out = {}
    for (typ, e, S), c in a.items():
        if x not in S:
            continue
        p = S.index(x)
        vadd(out, {(typ, e, tuple(i for i in S if i != x)): Fraction((-1) ** p)}, c)
    return out


def delta(a):
    out = {}
    for x in range(5):
        vadd(out, d_x(d_xi(a, x), x))
    return out


def split_parity(a):
    ev, od = {}, {}
    for k, c in a.items():
        (od if parity(k) else ev)[k] = c
    return ev, od


def bv_bracket(a, b):
    """{a,b} = del(ab) - del(a) b - (-1)^{|a|} a del(b), extended bilinearly over parity parts of a."""
    out = {}
    for sgn, part in ((1, split_parity(a)[0]), (-1, split_parity(a)[1])):
        if not part:
            continue
        vadd(out, delta(mul(part, b)))
        vadd(out, mul(delta(part), b), -1)
        vadd(out, mul(part, delta(b)), -sgn)
    return out


def schouten(a, b):
    """sum_x (d_xi_x a)(d_x b) + (-1)^{|a|} (d_x a)(d_xi_x b)."""
    out = {}
    for sgn, part in ((1, split_parity(a)[0]), (-1, split_parity(a)[1])):
        if not part:
            continue
        for x in range(5):
            vadd(out, mul(d_xi(part, x), d_x(b, x)))
            vadd(out, mul(d_x(part, x), d_xi(b, x)), sgn)
    return out


# ----------------------------------------------------------------- bases
SUBSETS = [S for r in range(6) for S in combinations(range(5), r)]


def basis(typ, W, pv=None):
    out = []
    for S in SUBSETS:
        if pv is not None and len(S) != pv:
            continue
        e = tuple(W[i] + (1 if i in S else 0) for i in range(5))
        if valid(typ, e):
            out.append((typ, e, S))
    return out


# ----------------------------------------------------------------- full complex T
def D(v, N=Fraction(1)):
    """D = t del + N {F,-} on T; v = {(k, typ, e, S): c}."""
    out = {}
    for (k, typ, e, S), c in v.items():
        m = {(typ, e, S): c}
        for kk, cc in delta(m).items():
            vadd(out, {(k + 1,) + kk: cc})
        if N:
            for kk, cc in bv_bracket(F, m).items():
                vadd(out, {(k,) + kk: N * cc})
    return out


def what(Tkey):
    k, typ, e, S = Tkey
    W = weight((typ, e, S))
    return W if typ == 2 else tuple(x - 1 for x in W)


def hdeg(Tkey):
    return Tkey[0] if Tkey[1] == 2 else Tkey[0] - 1


def T_basis(Wh, h):
    out = []
    if h >= 0:
        out += [(h,) + key for key in basis(2, Wh)]
    if h + 1 >= 0:
        out += [(h + 1,) + key for key in basis(0, tuple(x + 1 for x in Wh))]
    return out


def class_cohomology(Wh, h, N=Fraction(1)):
    """Exact data of H^h(T_What) and of the LES 0 -> coker phi_* -> H(T) -> ker phi_* -> 0.

    Returns dict with dim_H, coker (Laurent part = image of H(A2[[t]]) in H(T)),
    ker (regular part = image of H(T) in H(A0[[t]])), and the LES numbers computed from phi_*.
    """
    Bm, B0, Bp = T_basis(Wh, h - 1), T_basis(Wh, h), T_basis(Wh, h + 1)
    im_m = [D({b: 1}, N) for b in Bm]          # B_T^h
    im_0 = [D({b: 1}, N) for b in B0]
    for v in im_m + im_0:
        for kk in v:
            assert what(kk) == Wh and hdeg(kk) in (h, h + 1), (kk, Wh, h)
    r_in, r_out = rank(im_m), rank(im_0)
    dimH = len(B0) - r_out - r_in
    # cocycles of T^h
    Z = [{B0[i]: c for i, c in comb.items()} for comb in kernel(im_0)]
    # subcomplex S = A2 part, quotient Q = A0 part
    S0 = [b for b in B0 if b[1] == 2]
    Sm = [b for b in Bm if b[1] == 2]
    Sp = [b for b in Bp if b[1] == 2]
    Q0 = [b for b in B0 if b[1] == 0]
    Qm = [b for b in Bm if b[1] == 0]
    ZS = [{S0[i]: c for i, c in comb.items()} for comb in kernel([D({b: 1}, N) for b in S0])]
    coker_part = rel_rank(ZS, im_m)
    # regular part: projection of Z to Q modulo B_Q
    proj = lambda v: {k: c for k, c in v.items() if k[1] == 0}
    BQ = [proj(D({b: 1}, N)) for b in Qm]
    ker_part = rel_rank([proj(z) for z in Z], BQ)
    # LES through phi_*:  H^h(Q) -> H^{h+1}(S),  H^{h-1}(Q) -> H^h(S)
    def phi_data(Qdeg, Qprev, Sdeg, Sprev):
        imQ = [proj(D({b: 1}, N)) for b in Qdeg]
        ZQ = [{Qdeg[i]: c for i, c in comb.items()} for comb in kernel(imQ)]
        BQprev = [proj(D({b: 1}, N)) for b in Qprev]
        dimHQ = len(ZQ) - rank(BQprev)
        ZS_ = kernel([D({b: 1}, N) for b in Sdeg])
        BS = [D({b: 1}, N) for b in Sprev]
        dimHS = len(ZS_) - rank(BS)
        phi = [D(z, N) for z in ZQ]
        for v in phi:
            assert all(k[1] == 2 for k in v)
        rk = rel_rank(phi, BS)
        return dimHQ, dimHS, rk
    Qmm = [b for b in T_basis(Wh, h - 2) if b[1] == 0]
    # phi: H^{h-1}(Q) -> H^h(S)
    _, dimHS_h, rk_prev = phi_data(Qm, Qmm, S0, Sm)
    dimHQ_h, _, rk_next = phi_data(Q0, Qm, Sp, S0)
    les_coker = dimHS_h - rk_prev
    les_ker = dimHQ_h - rk_next
    return dict(dimH=dimH, coker=coker_part, ker=ker_part, les_coker=les_coker, les_ker=les_ker,
                dimHS=dimHS_h, dimHQ=dimHQ_h, Z=Z, B0=B0)


# ----------------------------------------------------------------- t^0 pieces
def ker_del(typ, W):
    B = basis(typ, W)
    return B, [{B[i]: c for i, c in comb.items()} for comb in kernel([delta({b: 1}) for b in B])]


def ker_phi_t0(Wr):
    """Regular t^0 classes at weight Wr: {alpha in A0_Wr : del alpha = 0, {F,alpha} = 0}."""
    B = basis(0, Wr)
    imgs = []
    for b in B:
        v = {}
        for kk, c in delta({b: 1}).items():
            vadd(v, {("del",) + kk: c})
        for kk, c in bv_bracket(F, {b: 1}).items():
            vadd(v, {("F",) + kk: c})
        imgs.append(v)
    return [{B[i]: c for i, c in comb.items()} for comb in kernel(imgs)]


def coker_phi_t0(Wn):
    """(ker del|A2_Wn, im phi_* = {F, ker del|A0_(Wn+1)}) as lists of vectors."""
    _, Z2 = ker_del(2, Wn)
    _, Z0 = ker_del(0, tuple(x + 1 for x in Wn))
    im = [bv_bracket(F, z) for z in Z0]
    im = [v for v in im if v]
    return Z2, im


def beta_of(Wn):
    """Boundary weight of a Laurent weight: beta = (W_z + 1, -W_w - 1)."""
    return (Wn[0] + 1, Wn[1] + 1, -Wn[2] - 1, -Wn[3] - 1, -Wn[4] - 1)


def ker_phi_formula(Wr):
    """Derived closed form of dim ker phi_*|t^0 at Wr (all entries >= -1)."""
    if min(Wr) < -1:
        return 0
    s = sum(1 for x in Wr if x >= 0)
    if s == 0:
        return 1
    rem = 1 if (Wr[0] >= 0 and Wr[1] >= 0 and all(x in (0, -1) for x in Wr[2:]) and Wr != ZERO5) else 0
    return 2 ** (s - 1) - rem


# ----------------------------------------------------------------- pairing candidates
def fischer_res(alpha_key, beta_key):
    """P(z^m w^b xi^S, z^m' w^(-a-1) xi^S') = m! [m=m'] [b=a] [S_z = S'_z] [S_w, S'_w complementary] * sign,
    sign defined by xi_w^{S_w} xi_w^{S'_w} = sign * xi_w1 xi_w2 xi_w3.  Weight rule: W_r = (W_n,z, -2 - W_n,w)."""
    t1, e1, S1 = alpha_key
    t2, e2, S2 = beta_key
    if t1 != 0 or t2 != 2:
        return 0
    if e1[:2] != e2[:2]:
        return 0
    if any(e1[2 + I] != -e2[2 + I] - 1 for I in range(3)):
        return 0
    Sz1 = tuple(i for i in S1 if i < 2)
    Sz2 = tuple(i for i in S2 if i < 2)
    if Sz1 != Sz2:
        return 0
    Sw1 = tuple(i for i in S1 if i >= 2)
    Sw2 = tuple(i for i in S2 if i >= 2)
    r = xi_mul(Sw1, Sw2)
    if r is None or r[1] != (2, 3, 4):
        return 0
    return r[0] * factorial(e1[0]) * factorial(e1[1])


def top_product(alpha_key, beta_key):
    """Contraction-type variant: Fischer in z-exponents, residue in w, complementary xi (all five):
    coefficient of z^m (x) w^(-1) xi_top in the product, times m!."""
    t1, e1, S1 = alpha_key
    t2, e2, S2 = beta_key
    if t1 != 0 or t2 != 2 or e1[:2] != e2[:2]:
        return 0
    if any(e1[2 + I] != -e2[2 + I] - 1 for I in range(3)):
        return 0
    r = xi_mul(S1, S2)
    if r is None or r[1] != (0, 1, 2, 3, 4):
        return 0
    return r[0] * factorial(e1[0]) * factorial(e1[1])


def pair(P, a, b):
    return sum((ca * cb * P(ka, kb) for ka, ca in a.items() for kb, cb in b.items()), Fraction(0))


def regular_partner(Wn):
    """Weight rule of fischer_res: W_r = (W_n,z, -2 - W_n,w) = beta(W_n) - (1,1,1,1,1)."""
    return (Wn[0], Wn[1]) + tuple(-2 - x for x in Wn[2:])


def preimage(target, keys, op):
    """x in span(keys) with op(x) = target, or None."""
    E = Echelon(track=True)
    for i, b in enumerate(keys):
        E.add(op({b: 1}), label=i)
    row, comb = E.reduce(target)
    if row:
        return None
    x = {}
    for i, c in comb.items():
        vadd(x, {keys[i]: Fraction(1)}, -c)
    assert not vadd(dict(op(x)), target, -1)
    return x


def pairing_data(Wn, P=fischer_res):
    """Regular modes K = ker phi_*|t^0 at W_r = regular_partner(Wn) against Laurent modes at Wn.

    direct:  <alpha, beta> = P(alpha, beta), beta in ker del|A2_Wn; descends to coker iff P(K, im phi_*) = 0.
    KS:      <alpha, beta> = P(alpha, del^{-1} beta) on exact beta; independent of the primitive iff
             P(K, ker del|A2_Wn) = 0; descends iff also P(K, F.ker del|A0_(Wn+1)) = 0.
    """
    Wr = regular_partner(Wn)
    K = ker_phi_t0(Wr)
    B2, Z2 = ker_del(2, Wn)
    _, Z0 = ker_del(0, tuple(x + 1 for x in Wn))
    im = [v for v in (bv_bracket(F, z) for z in Z0) if v]
    dimC = len(Z2) - rank(im)
    dir_desc = all(pair(P, k, v) == 0 for k in K for v in im)
    # rank on K x coker: rank of [P(k, z)] restricted to a complement of im (valid when dir_desc)
    E = Echelon()
    for v in im:
        E.add(v)
    comp = []
    for z in Z2:
        if E.add(z)[0]:
            comp.append(z)
    dir_rank = rank([{j: pair(P, k, z) for j, z in enumerate(comp)} for k in K])
    ks_indep = all(pair(P, k, z) == 0 for k in K for z in Z2)
    Fg = [mul(F, z) for z in Z0]
    ks_desc = ks_indep and all(pair(P, k, g) == 0 for k in K for g in Fg)
    gam = [preimage(z, B2, delta) for z in comp]
    ks_rank = rank([{j: pair(P, k, g) for j, g in enumerate(gam) if g is not None} for k in K])
    # closed regular x closed Laurent: Z(A0)_Wr x Z(A2)_Wn
    _, Z0r = ker_del(0, Wr)
    closed_rank = rank([{j: pair(P, a, z) for j, z in enumerate(Z2)} for a in Z0r])
    # K^perp(B^1): dim of {k in K : P(k, im phi_*) = 0}
    orth_dim = len(K) - rank([{j: pair(P, k, v) for j, v in enumerate(im)} for k in K]) if im else len(K)
    return dict(Wr=Wr, beta=beta_of(Wn), dimK=len(K), dimC=dimC, dir_desc=dir_desc, dir_rank=dir_rank,
                ks_indep=ks_indep, ks_desc=ks_desc, ks_rank=ks_rank, n_exact=sum(g is not None for g in gam),
                dimZ0=len(Z0r), dimZ2=len(Z2), closed_rank=closed_rank, dimB1=rank(im), orth_dim=orth_dim)
