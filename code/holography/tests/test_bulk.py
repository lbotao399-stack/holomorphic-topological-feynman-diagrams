"""Exact checks of holography.bulk_x and holography.derham.

Run:  python3 code/holography/tests/test_bulk.py      (also collected by pytest)
"""
import importlib.util
import itertools
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from fractions import Fraction
from holography import bulk_x as bx, derham as dr, linalg as la

P = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SINGLE = os.path.join(P, "explore", "twisted-holography-N4", "code", "single_particle.py")


def rand_elem(rng, typ, W, n=3):
    B = bx.basis(typ, W)
    if not B:
        return {}
    return la.clean({rng.choice(B): Fraction(rng.randint(-3, 3)) for _ in range(n)})


def test_bv_equals_schouten():
    rng = random.Random(1)
    for _ in range(300):
        Wa = tuple(rng.randint(-1, 1) for _ in range(5))
        Wb = tuple(rng.randint(-1, 1) for _ in range(5))
        a, b = rand_elem(rng, 0, Wa), rand_elem(rng, 0, Wb)
        assert bx.bv_bracket(a, b) == bx.schouten(a, b)
        Wc = (rng.randint(-1, 1), rng.randint(-1, 1)) + tuple(rng.randint(-3, -1) for _ in range(3))
        c = rand_elem(rng, 2, Wc)
        assert bx.bv_bracket(c, a) == bx.schouten(c, a)
        assert bx.bv_bracket(a, c) == bx.schouten(a, c)
    for W in itertools.product(range(-1, 2), repeat=5):
        for k in bx.basis(0, W):
            assert bx.bv_bracket(bx.F, {k: 1}) == bx.schouten(bx.F, {k: 1})


def test_bv_identities():
    """{a,b} = (-1)^{|a||b|} {b,a};  del{a,b} = -{del a, b} - (-1)^{|a|}{a, del b};  {F,F} = 0, del F = 0.
    ((-1)^{|a|}{a,b} is the Gerstenhaber-signed bracket, antisymmetric with -(-1)^{(|a|+1)(|b|+1)}.)"""
    rng = random.Random(2)
    for _ in range(200):
        Wa = tuple(rng.randint(-1, 1) for _ in range(5))
        Wb = tuple(rng.randint(-1, 1) for _ in range(5))
        for pa in range(6):
            a = {k: Fraction(rng.randint(1, 3)) for k in bx.basis(0, Wa, pa)[:2]}
            if not a:
                continue
            for pb in range(6):
                b = {k: Fraction(rng.randint(1, 3)) for k in bx.basis(0, Wb, pb)[:2]}
                if not b:
                    continue
                lhs = bx.bv_bracket(a, b)
                rhs = la.vscale(bx.bv_bracket(b, a), (-1) ** (pa * pb))
                assert lhs == la.clean(rhs)
                l2 = bx.delta(lhs)
                r2 = la.vadd(la.vscale(bx.bv_bracket(bx.delta(a), b), -1), bx.bv_bracket(a, bx.delta(b)), -(-1) ** pa)
                assert l2 == r2
    assert bx.bv_bracket(bx.F, bx.F) == {}
    assert bx.delta(bx.F) == {}


def test_D_squared():
    for Wh in itertools.product(range(-2, 2), range(-2, 2), range(-3, 1), range(-3, 1), range(-3, 1)):
        for h in range(-1, 3):
            for b in bx.T_basis(Wh, h):
                assert bx.D(bx.D({b: 1})) == {}
                assert bx.D(bx.D({b: 1}, N=Fraction(7)), N=Fraction(7)) == {}


def test_phi_chain_map():
    for W in itertools.product(range(-1, 2), range(-1, 2), range(-3, 0), range(-3, 0), range(-3, 0)):
        for k in bx.basis(2, W):
            assert dr.Phi(bx.delta({k: 1})) == dr.d_form(dr.Phi({k: 1}))
            assert dr.form_weight(next(iter(dr.Phi({k: 1})))) == bx.beta_of(W)
    assert dr.Phi(bx.F) == {((0, 0), (), (), (0, 0, 0)): 1}


def t0_coker(W):
    """dim ker(del|A2_W) - rank del(F . ker(del|A0_{W+1})) via the BV bracket."""
    Z2, im = bx.coker_phi_t0(W)
    return len(Z2) - la.rank(im)


def test_t0_sixty_slots():
    """Rerun of single_particle.py: coker(W) = dim B^{>=2}(beta) + [beta=0] on the 60 slots beta <= (2,2;2,2,2)."""
    spec = importlib.util.spec_from_file_location("single_particle", SINGLE)
    sp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sp)
    slots = 0
    for n1 in range(3):
        for n2 in range(n1, 3):
            for th in itertools.product(range(3), repeat=3):
                if list(th) != sorted(th, reverse=True):
                    continue
                beta = (n1, n2) + th
                W = (n1 - 1, n2 - 1) + tuple(-t - 1 for t in th)
                B = dr.B_list(beta)
                pred = sum(B[2:]) + (1 if beta == (0,) * 5 else 0)
                mine = t0_coker(W)
                old = sp.bulk_coker(W)[0]
                oldB = [sp.bdy_B(beta, p) for p in range(len(B))]
                assert mine == pred == old, (beta, mine, pred, old)
                assert oldB == B, (beta, oldB, B)
                r = bx.class_cohomology(W, 0)
                assert r['coker'] == pred and r['les_coker'] == pred
                slots += 1
    assert slots == 60


def test_les_small_window():
    for Wh in itertools.product(range(-2, 2), range(-2, 2), range(-3, 0), range(-3, 0), range(-3, 0)):
        for h in range(-1, 3):
            r = bx.class_cohomology(Wh, h)
            assert r['dimH'] == r['coker'] + r['ker']
            assert r['les_coker'] == r['coker'] and r['les_ker'] == r['ker']
            if h == -1:
                assert r['ker'] == bx.ker_phi_formula(tuple(x + 1 for x in Wh))


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
