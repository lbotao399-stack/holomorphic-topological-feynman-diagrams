"""Exact checks of holography.symmetry (bulk symmetry algebra g = ker phi_*|t^0, C^{2|3} dictionary, psl(3|3)).

Run:  python3 code/holography/tests/test_symmetry.py      (also collected by pytest)
"""
import itertools
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from fractions import Fraction
from holography import bulk_x as bx, derham as dr, linalg as la, symmetry as sy

BOX = list(itertools.product(range(-1, 2), repeat=5))
A2BOX = list(itertools.product(range(-1, 2), range(-1, 2), range(-3, 0), range(-3, 0), range(-3, 0)))


def _keys():
    rng = random.Random(5)
    A0 = [k for W in BOX for k in bx.basis(0, W)]
    A2 = [z for W in A2BOX for z in bx.basis(2, W)]
    return rng, A0, A2


def test_contraction_and_bracket_dictionary():
    """Phi(alpha.zeta) = iota_{Psi alpha} Phi(zeta);  Phi{alpha,zeta} = (L_{Psi alpha} - iota_{Psi del alpha}) Phi(zeta)."""
    rng, A0, A2 = _keys()
    for k in rng.sample(A0, 300):
        a = {k: Fraction(1)}
        for z in rng.sample(A2, 10):
            pr = bx.mul(a, {z: 1})
            assert (dr.Phi(pr) if pr else {}) == la.clean(sy.iota(a, dr.Phi({z: 1})))
            om = dr.Phi({z: 1})
            rhs = sy.L_pv(a, om)
            la.vadd(rhs, sy.iota(bx.delta(a), om), -1)
            assert la.clean(sy.bulk_on_forms(a, om)) == la.clean(rhs)
    for z in rng.sample(A2, 300):
        assert sy.phi_inv(dr.Phi({z: 1})) == {z: 1}


def test_flux_map_is_function_part():
    """Phi(F.alpha) = (Psi alpha)_{k=0}: {F, alpha} = del(F alpha) <-> d of the function part."""
    for W in itertools.product(range(-1, 3), repeat=5):
        for k in bx.basis(0, W):
            pr = bx.mul(bx.F, {k: 1})
            lhs = dr.Phi(pr) if pr else {}
            s, (m, K, c, T) = sy.psi_key(k)
            rhs = {(m, (), K, (0, 0, 0)): Fraction(s)} if sy.pv_degree(k) == 0 else {}
            assert lhs == rhs


def test_vector_fields_lie_derivative():
    rng, A0, A2 = _keys()
    for k in [k for k in A0 if sy.pv_degree(k) == 1]:
        a = {k: Fraction(1)}
        X, p = sy.vector_field_of(a), sy.spar(a)
        for z in rng.sample(A2, 3):
            om = dr.Phi({z: 1})
            assert la.clean(sy.L_pv(a, om)) == la.clean(la.vscale(sy.lie_derivative(X, om, p), (-1) ** p))


def test_lie_superalgebra():
    rng = random.Random(9)
    def relem(W):
        p = rng.randint(0, 1)
        B = [b for b in bx.basis(0, W) if len(b[2]) % 2 == p]
        return la.clean({rng.choice(B): Fraction(rng.randint(-3, 3)) for _ in range(3)}) if B else {}
    n = 0
    while n < 300:
        a, b, c = (relem(tuple(rng.randint(-1, 1) for _ in range(5))) for _ in range(3))
        if not (a and b and c):
            continue
        pa, pb = sy.spar(a), sy.spar(b)
        assert sy.lie(a, b) == la.clean(la.vscale(sy.lie(b, a), -(-1) ** (pa * pb)))
        assert sy.lie(a, sy.lie(b, c)) == la.clean(la.vadd(dict(sy.lie(sy.lie(a, b), c)), sy.lie(b, sy.lie(a, c)), (-1) ** (pa * pb)))
        n += 1


def test_psl33():
    basis = sy.psl33_basis()
    vecs = [v for (_, _, v) in basis]
    assert sy.psi_mat({(P, P): Fraction(1) for P in range(6)}) == {}
    assert la.rank(vecs) == 34
    assert all(not bx.delta(v) and sy.alpha0_is_constant(v) for v in vecs)
    for (na, A, va) in basis:
        for (nb, B, vb) in basis:
            assert sy.lie(va, vb) == sy.psi_mat(sy.supercomm(A, B))
    f = sy.structure_constants(vecs)
    assert not any(sy.killing_form(vecs, f).values())
    assert [sy.spar(v) for v in vecs].count(0) == 16
    assert all(not sy.lie(sy.UNIT, v) for v in vecs)


def test_g_structure():
    """dim g(W) = sum_{k>=1} dim ker Delta_k + [W=0] = dim Z(A0)(W) - dim B^1(gamma(W))."""
    for W in BOX:
        dims = sy.g_degree_dims(W)
        g = len(sy.g_basis(W))
        assert g == sum(v for k, v in dims.items() if k >= 1) + (W == sy.ZERO5)
        assert g == len(sy.Z_A0(W)) - sy.B1_dim(sy.c23_weight(W))


def test_excess_radicals_u1_sources():
    def inR(b):
        return all(x in (0, 1) for x in b[2:]) and min(b[:2]) >= 1 and b != (1,) * 5
    for b in itertools.product(range(3), repeat=5):
        if list(b[:2]) != sorted(b[:2], reverse=True) or list(b[2:]) != sorted(b[2:], reverse=True):
            continue
        Wr = tuple(x - 1 for x in b)
        Wn = (b[0] - 1, b[1] - 1) + tuple(-x - 1 for x in b[2:])
        Z2, im = bx.coker_phi_t0(Wn)
        cok = len(Z2) - la.rank(im)
        g = len(sy.g_basis(Wr))
        assert g - cok == sy.excess_formula(b) == (1 if (sy.in_M1(b) or b == (1,) * 5) else 0)
        d = sy.projected_pairing_data(b)
        assert d['left_rad'] == sy.B1_dim(b) and d['right_rad'] == (1 if inR(b) else 0) and d['dim_coker'] == cok
    for bz in itertools.product(range(3), repeat=2):
        for bt in itertools.product(range(2), repeat=3):
            b = bz + bt
            if b == (0,) * 5:
                continue
            us = sy.u1_source(b) if bz != (0, 0) else sy.ks_source_z0(b)
            assert len(us) == 1 and sy.proportional(us[0], sy.u1_formula(b))


def test_trivial_action():
    Wn_list = [(b[0] - 1, b[1] - 1) + tuple(-x - 1 for x in b[2:]) for b in itertools.product(range(3), repeat=5)]
    Wn_list.sort(key=lambda W: sum(abs(x) for x in W))
    cur, _ = sy.residual_trivial(sy.ZERO5, Wn_list, stop_at=1)
    assert len(cur) == 1 and sy.proportional(cur[0], sy.UNIT)
    for W in [(-1, 0, 0, 0, 0), (0, 0, 1, 0, 0), (-1, -1, 0, -1, -1), (1, 0, 1, 1, 1)]:
        cur, _ = sy.residual_trivial(W, Wn_list)
        assert cur == []


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
