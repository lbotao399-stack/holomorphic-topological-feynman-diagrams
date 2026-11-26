import random, time
from fractions import Fraction as Fr
from vce import *
from vnum import ev, d_explicit, bracket_explicit, basis_gl, basis_sl, inv_gram, rand_point

rng = random.Random(20261121)
t0 = time.time()


def rand_cochain(k, a, b, sl, nterms=6, tries=20):
    ms = monos(k, a, b)
    ps = [p for p in ALLPERMS[k] if not (sl and has_fixed(p))]
    for _ in range(tries):
        w = {}
        for _ in range(nterms):
            add(w, (rng.choice(ps), rng.choice(ms), rng.choice([0, 0, 1])), Fr(rng.randint(-3, 3)))
        w = alt(w, k)
        if w:
            return w
    return {}


def is_alt(w):
    k = arity(w)
    return all(lin((1, act(pi, w)), (-perm_sign(pi), w)) == {} for pi in ALLPERMS[k])


s0 = from_poly_tau({tuple([0] * 6): 1}, [(1, (1, 2, 0)), (-1, (2, 0, 1))])
print("s0 = tr(012) - tr(021) alternating:", is_alt(s0))
for sl in (False, True):
    for (k, a, b) in ((2, 1, 0), (3, 1, 1), (3, 2, 2), (4, 1, 1)):
        w = rand_cochain(k, a, b, sl)
        if not w:
            continue
        print(f"sl={sl} k={k} ({a},{b}) |w|={len(w)}: alternating {is_alt(w)}, d^2=0 {d(d(w, sl), sl) == {}},"
              f" d w alternating {is_alt(d(w, sl)) if d(w, sl) else 'dw=0'}")
    for (k, a, b) in ((2, 1, 0), (3, 1, 1), (3, 2, 2), (4, 2, 2)):
        w = rand_cochain(k, a, b, sl)
        if not w:
            continue
        L = bracket(s0 if not sl else restrict_sl(s0), w, sl)
        R = bracket(w, s0 if not sl else restrict_sl(s0), sl)
        dw = d(w, sl)
        e1 = lin((1, L), (1, dw)) == {}
        e2 = lin((1, R), (-1, L)) == {}
        e3 = lin((1, R), (1, L)) == {}
        print(f"sl={sl} k={k} ({a},{b}): {{s0,w}} = -dw {e1};  {{w,s0}} = {{s0,w}} {e2};  {{w,s0}} = -{{s0,w}} {e3}")
    print(f"sl={sl}: {{s0,s0}} = 0:", bracket(s0 if not sl else restrict_sl(s0), s0 if not sl else restrict_sl(s0), sl) == {})
    a3 = rand_cochain(3, 1, 1, sl)
    b3 = rand_cochain(3, 1, 1, sl)
    print(f"sl={sl}: {{a,b}} = {{b,a}} (3,3):", lin((1, bracket(a3, b3, sl)), (-1, bracket(b3, a3, sl))) == {})
    dab = d(bracket(a3, b3, sl), sl)
    x = bracket(d(a3, sl), b3, sl)
    y = bracket(a3, d(b3, sl), sl)
    for s1 in (1, -1):
        for s2 in (1, -1):
            if lin((1, dab), (-s1, x), (-s2, y)) == {}:
                print(f"sl={sl}: d{{a,b}} = {s1}{{da,b}} + {s2}{{a,db}}  (a,b arity 3)")

# explicit checks
for (N, sl) in ((2, False), (3, False), (2, True), (3, True)):
    B = basis_sl(N) if sl else basis_gl(N)
    Gi = inv_gram(B)
    okd = okb = True
    for (k, a, b) in ((3, 1, 1), (4, 2, 2)):
        w = rand_cochain(k, a, b, sl, nterms=4)
        dw = d(w, sl)
        for _ in range(2):
            X, p = rand_point(N, k + 1, sl, rng)
            okd &= (ev(dw, X, p, N) == d_explicit(w, X, p, N))
    for (m, n) in ((3, 3), (3, 4), (4, 3)):
        A = rand_cochain(m, 1, 1, sl, nterms=4)
        Bc = rand_cochain(n, 1, 1 if n == 3 else 0, sl, nterms=4)
        br = bracket(A, Bc, sl)
        for _ in range(2):
            X, p = rand_point(N, m + n - 2, sl, rng)
            v1 = ev(br, X, p, N)
            v2 = bracket_explicit(A, m, Bc, n, X, p, N, B, Gi)
            okb &= (v1 == v2)
    print(f"{'sl' if sl else 'gl'}{N}: formal d == explicit CE sum: {okd};  formal bracket == explicit dual-basis shuffle sum: {okb}")
print(f"time {time.time() - t0:.1f}s")
