"""Pointwise exact check (random integer X_0..X_3 in [-3,3]^dim, momenta in [-9,9]^2) of
F_133 + F_223 = c beta (c from the b_theta evaluation) and F_223 != c' beta for every c' (two points with
non-proportional (beta, F_223) values certify F_223 not in C beta independently of the injectivity argument)."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import numpy as np, sympy as sp
import algebras as al, core, closedness as cl
for key in sys.argv[1:]:
    t0 = time.time()
    L = al.build(key); T = core.Tensors(L); th = L.theta
    F133 = core.F_eval(T, '133', th); F223 = core.F_eval(T, '223', th); be = core.beta_eval(L, th)
    ok, sol = core.in_span(core.to_rows([be, F133 + F223]), 1)
    cint = sol[0]
    rng = np.random.default_rng(2026)
    pts = [cl.pointwise(L, T, rng) for _ in range(4)]
    rel = all(a + b == cint * c for a, b, c in pts)
    minors = [pts[0][2] * q[1] - pts[0][1] * q[2] for q in pts[1:]]
    print("%-14s c (s^7-scaled) = %s ; F133+F223 = c beta at 4 random points: %s ; beta values %s ; "
          "2x2 minors beta(P0)F223(Pk) - F223(P0)beta(Pk), k=1..3: %s  (%.1fs)"
          % (L.name, cint, rel, [q[2] for q in pts], minors, time.time() - t0), flush=True)
