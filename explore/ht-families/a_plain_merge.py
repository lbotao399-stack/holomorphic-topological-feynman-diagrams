"""(a) Plain k-merge: insertion I_k(J, Theta) of k parallel plain branches (exact check at rational t).
Conjecture:  I_k = s_k * t^{k-1} prod_i w_i^{-1} dw_1^...^dw_{k-1} * (2t)^{-(k-1)/2} He_{k-1}(sqrt(2t) hatJ),
hatJ = J + (dt/2t) Theta,  s_k = (-1)^{binom(k-1,2)},  t = t_red, w_i = t / t_i,  C_rec loops: branch a - branch k."""
import sys, random
from fractions import Fraction as Fr
import sympy as sp
import hf_recursion as hr
from hf_brute import wedge, fadd

def He_coeffs(n):
    x = sp.Symbol('x')
    return sp.Poly(sp.hermite_prob(n, x), x).all_coeffs()[::-1]   # coefficient of x^p

def predicted(k, t):
    tt = 1 / sum(1 / Fr(x) for x in t)
    w = [tt / Fr(x) for x in t]
    # forms in dt_0..dt_{k-1}
    tau = {(i,): w[i] / Fr(t[i]) for i in range(k)}            # dt/t = sum_i w_i dt_i / t_i
    dt = {(i,): w[i] ** 2 for i in range(k)}
    dw = [fadd({key: w[i] * v for key, v in tau.items()}, {(i,): -w[i] / Fr(t[i])}) for i in range(k)]
    vol = {(): Fr(1)}
    for a in range(k - 1):
        vol = wedge(vol, dw[a])
    pref = (-1) ** ((k - 1) * (k - 2) // 2) * tt ** (k - 1)
    for x in w:
        pref /= x
    vol = {key: pref * v for key, v in vol.items()}
    # (2t)^{-(k-1)/2} He_{k-1}(sqrt(2t) X) = sum_p h_p (2t)^{(p-(k-1))/2} X^p ; only p = k-1 mod 2 terms
    h = He_coeffs(k - 1)
    out = {}
    eps = {key: v / (2 * tt) for key, v in dt.items()}
    for p, hp in enumerate(h):
        if hp == 0:
            continue
        e2 = (p - (k - 1))
        assert e2 % 2 == 0
        c = Fr(int(hp)) * (2 * tt) ** (e2 // 2)
        # c * X^p, X = J + eps Theta:  J^p + p J^{p-1} eps Theta
        for key, v in vol.items():
            out[(key, (p,))] = out.get((key, (p,)), 0) + c * v
        if p >= 1:
            for key, v in wedge(vol, eps).items():
                kk = (tuple(key) + (hr.THETA,), (p - 1,))
                out[kk] = out.get(kk, 0) + c * p * v
    return {kk: v for kk, v in out.items() if v != 0}

rng = random.Random(11)
for k in range(2, 7):
    edges = [(0, 1)] * k
    plan = [(0, 1, [[e] for e in range(k)], [[1] for e in range(k)])]
    for trial in range(2):
        t = [Fr(rng.randint(1, 20), rng.randint(1, 7)) for _ in range(k)]
        rec = hr.Recursion(2, edges, plan)
        rec.run(t)
        I = rec.ins[k]
        P = predicted(k, t)
        print("k = %d: I_k == conjecture: %s  (#terms %d)" % (k, I == P, len(I)))
