"""S3: action of g on normalizable modes coker phi_*, kernel of the action, and comparison with the boundary
Lie-derivative action of psl(3|3) (projective vector fields on C^{2|3}) on single traces B^p.  Exact.

python3 run_S3.py > S3.log
"""
import sys, itertools, random, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from fractions import Fraction
from holography import bulk_x as bx, derham as dr, linalg as la, symmetry as sy

t0 = time.time()
rng = random.Random(11)
# ---------------------------------------------------------------- kernel of rho: g -> End(coker)
Wn_list = [(b[0] - 1, b[1] - 1) + tuple(-x - 1 for x in b[2:]) for b in itertools.product(range(3), repeat=5)]
Wn_list.sort(key=lambda W: sum(abs(x) for x in W))
nW = 0; resid = {}
for W in itertools.product(range(-1, 2), repeat=5):
    if list(W[:2]) != sorted(W[:2], reverse=True) or list(W[2:]) != sorted(W[2:], reverse=True):
        continue
    g = sy.g_basis(W)
    if not g:
        continue
    cur, used = sy.residual_trivial(W, Wn_list, stop_at=1 if W == sy.ZERO5 else 0)
    nW += 1
    resid[W] = [sy.pretty_pv(c) for c in cur]
nz = {W: r for W, r in resid.items() if r}
print("[chk] elements of g(W), W in [-1,1]^5 mod S2xS3 (%d weights, all of g there), acting trivially on coker(Wn)" % nW)
print("      for all Laurent Wn with beta(Wn) in [0,2]^5:", nz)
print("[drv] 1 is central for the Schouten bracket ({1,-} = 0), so C.1 acts trivially on all of coker.")
print("[drv] a nonzero vector field X on C^{2|3} acts nontrivially on B^2: L_X(dx_i dx_j) = 0 for all i,j forces dX^{theta_I} = 0")
print("      (i = j = theta_I, dtheta_I^2 != 0) and then dX^{z_a} = 0 (i = z_a, j = theta_I); constant X are detected by d(z_a^2 dtheta_1),")
print("      d(theta_1 theta_2 dtheta_3).  B^2 is in coker and L_X B^2 is in B^2, disjoint from B^1.")

# ---------------------------------------------------------------- boundary psl(3|3) = bulk psl(3|3) on normalizable modes
def y0(P):
    if P == 0:
        return {((0, 0), (), (), (0, 0, 0)): Fraction(1)}
    if P in (1, 2):
        m = [0, 0]; m[P - 1] = 1
        return {(tuple(m), (), (), (0, 0, 0)): Fraction(1)}
    return {((0, 0), (), (P - 3,), (0, 0, 0)): Fraction(1)}
def V_vf(P, Q):
    """y_P d/dy_Q (Q != 0) or -y_P E on C^{2|3}, coordinates 0,1 = z1,z2; 2,3,4 = theta1..3."""
    if Q != 0:
        return {Q - 1: y0(P)}
    X = {}
    for i in range(1, 6):
        X[i - 1] = la.vscale(sy.form_mul(y0(P), y0(i)), -1)
    return {k: la.clean(v) for k, v in X.items() if la.clean(v)}
par = sy.PAR6
forms = []
for b in itertools.product(range(3), repeat=5):
    for p in range(0, 5):
        for f in dr.forms(b, p):
            om = dr.d_form({f: Fraction(1)})
            if om:
                forms.append(om)
print("\nexact test forms d(monomial), weight in [0,2]^5:", len(forms))
bad = n = 0
for (name, M, x) in sy.psl33_basis():
    if name[0] == 'E':
        P, Q = int(name[1]), int(name[2])
        X = V_vf(P, Q); parity = (par[P] + par[Q]) % 2
    else:
        P, Q = int(name[1]), int(name[2])
        s = M[(Q, Q)]
        X = {}
        for (PP, cc) in ((P, Fraction(1)), (Q, s)):
            for k, v in V_vf(PP, PP).items():
                X.setdefault(k, {}); la.vadd(X[k], v, cc)
        X = {k: la.clean(v) for k, v in X.items() if la.clean(v)}
        parity = 0
    for om in rng.sample(forms, 150):
        lhs = la.clean(dr.Phi(sy.lie(x, sy.phi_inv(om))))
        rhs = la.clean(sy.lie_derivative(X, om, parity))
        n += 1
        bad += lhs != rhs
print("[chk] Phi [psi(A), Phi^{-1} omega] = L_{V_A} omega  (V_A = projective vector field of A, standard Lie derivative)")
print("      for all 34 basis elements of psl(3|3) on 150 exact forms each: %d pairs, %d mismatches" % (n, bad))
bad = n = 0
for (name, M, x) in sy.psl33_basis()[:30]:
    P, Q = int(name[1]), int(name[2])
    X = V_vf(P, Q); parity = (par[P] + par[Q]) % 2
    for f in rng.sample([k for b in itertools.product(range(3), repeat=5) for k in dr.forms(b)], 80):
        om = {f: Fraction(1)}
        l1 = dr.d_form(sy.lie_derivative(X, om, parity))
        l2 = la.vscale(sy.lie_derivative(X, dr.d_form(om), parity), (-1) ** parity)
        n += 1
        bad += la.clean(l1) != la.clean(l2)
print("[chk] d L_V = (-1)^{|V|} L_V d on %d (generator, form monomial) pairs: %d failures  (so L_V preserves B^p)" % (n, bad))
print("time %.1fs" % (time.time() - t0))
