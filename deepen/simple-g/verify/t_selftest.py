"""Self-tests: W_full vs W_point; alt/beta tensors vs pointwise evaluation; transpose convention."""
import sys, time, itertools
import numpy as np
import valg, vcoch as vc
rng = np.random.default_rng(7)
# transpose convention: V = transpose(W, argsort(s)) satisfies V[a0..a3] = W[a_s0..a_s3]
W = rng.integers(-9, 9, (2, 3, 4, 5))
ok = True
for s in itertools.permutations(range(4)):
    V = np.transpose(W, list(np.argsort(s)))
    for _ in range(5):
        shp = V.shape
        a = [rng.integers(0, shp[k]) for k in range(4)]
        ok &= V[tuple(a)] == W[a[s[0]], a[s[1]], a[s[2]], a[s[3]]]
print("transpose convention:", bool(ok))
for name in sys.argv[1:]:
    t0 = time.time()
    A = valg.build(name)
    S = vc.Setup(A)
    W = {G: S.W_full(G) for G in ('133', '223')}
    ok1 = True
    for _ in range(3):
        X = [rng.integers(-3, 4, A.n) for _ in range(4)]
        for G in ('133', '223'):
            full = int(np.einsum('abcd,a,b,c,d->', W[G].astype(object), *[x.astype(object) for x in X]))
            pt = S.W_point(G, X)
            ok1 &= (full == pt)
    ok2 = True
    for _ in range(2):
        X = [rng.integers(-3, 4, A.n) for _ in range(4)]
        P = vc.config4(rng)
        for G in ('133', '223'):
            T = vc.alt_tensor(W[G], P)
            full = int(np.einsum('abcd,a,b,c,d->', T.astype(object), *[x.astype(object) for x in X]))
            pt = vc.alt_point(lambda Y: S.W_point(G, Y), X, P)
            ok2 &= full == pt
        B = vc.beta_tensor(S.TT, P)
        full = int(np.einsum('abcd,a,b,c,d->', B.astype(object), *[x.astype(object) for x in X]))
        ok2 &= full == vc.beta_point(S, X, P)
    print("%-5s g = %d (kappa = K/g), s = %d ; W_full = W_point (3 random 4-tuples, both graphs): %s ; "
          "alt/beta tensors = pointwise: %s ; (%.1fs)" % (name, S.g, S.s, bool(ok1), bool(ok2), time.time() - t0), flush=True)
print("arith stats:", vc.STATS)
