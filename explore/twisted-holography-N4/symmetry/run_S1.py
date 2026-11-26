"""S1: bulk symmetry algebra g = ker phi_*|t^0: dictionary with PV(C^{2|3}), Lie structure, psl(3|3), U(1) sector.  Exact.

python3 run_S1.py > S1.log
"""
import sys, itertools, random, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from fractions import Fraction
from holography import bulk_x as bx, derham as dr, linalg as la, symmetry as sy

t0 = time.time()
rng = random.Random(2026)
box = list(itertools.product(range(-1, 2), repeat=5))
A0keys = [k for W in box for k in bx.basis(0, W)]
A2keys = [z for W in itertools.product(range(-1, 2), range(-1, 2), range(-3, 0), range(-3, 0), range(-3, 0))
          for z in bx.basis(2, W)]
print("A0 monomials with W in [-1,1]^5:", len(A0keys), "; A2 monomials in box:", len(A2keys))

# ---------------------------------------------------------------- C1-C3 dictionary
bad = n = 0
for k in A0keys:
    for z in rng.sample(A2keys, 40):
        pr = bx.mul({k: 1}, {z: 1})
        lhs = dr.Phi(pr) if pr else {}
        n += 1
        bad += lhs != la.clean(sy.iota({k: Fraction(1)}, dr.Phi({z: 1})))
print("[chk] C1  Phi(alpha.zeta) = iota_{Psi alpha} Phi(zeta): %d pairs, %d mismatches" % (n, bad))
bad = n = 0
for z in A2keys:
    bad += sy.phi_inv(dr.Phi({z: 1})) != {z: 1}
print("[chk]     Phi^{-1} o Phi = 1 on %d A2 monomials: %d mismatches" % (len(A2keys), bad))
bad = n = 0
for k in rng.sample(A0keys, 1200):
    a = {k: Fraction(1)}
    for z in rng.sample(A2keys, 15):
        om = dr.Phi({z: 1})
        rhs = sy.L_pv(a, om)
        la.vadd(rhs, sy.iota(bx.delta(a), om), -1)
        n += 1
        bad += la.clean(sy.bulk_on_forms(a, om)) != la.clean(rhs)
print("[chk] C2  Phi{alpha,zeta} = (L_{Psi alpha} - iota_{Psi del alpha}) Phi(zeta): %d pairs, %d mismatches" % (n, bad))
bad = n = 0
vf = [k for k in A0keys if sy.pv_degree(k) == 1]
for k in vf:
    a = {k: Fraction(1)}
    X = sy.vector_field_of(a)
    p = sy.spar(a)
    for z in rng.sample(A2keys, 8):
        om = dr.Phi({z: 1})
        lhs = la.clean(sy.L_pv(a, om))
        rhs = la.clean(la.vscale(sy.lie_derivative(X, om, p), (-1) ** p))
        n += 1
        bad += lhs != rhs
print("[chk] C3  C^{2|3}-vector fields (%d monomials): L_{Psi alpha} = (-1)^{p(alpha)} (standard Lie derivative): %d pairs, %d mismatches"
      % (len(vf), n, bad))
print("          hence rho(alpha) := [alpha,-] = (-1)^{|alpha|}{alpha,-} acts on Laurent modes as -L_{X(alpha)}, X(alpha) = vector field of Psi(alpha)")

# ---------------------------------------------------------------- Lie superalgebra axioms
def relem(W, n=3):
    p = rng.randint(0, 1)
    B = [b for b in bx.basis(0, W) if len(b[2]) % 2 == p]
    if not B:
        return {}
    return la.clean({rng.choice(B): Fraction(rng.randint(-3, 3)) for _ in range(n)})
bad = n = 0
while n < 1500:
    a, b, c = (relem(tuple(rng.randint(-1, 1) for _ in range(5))) for _ in range(3))
    if not (a and b and c):
        continue
    pa, pb = sy.spar(a), sy.spar(b)
    bad += sy.lie(a, b) != la.clean(la.vscale(sy.lie(b, a), -(-1) ** (pa * pb)))
    l = sy.lie(a, sy.lie(b, c))
    r = la.clean(la.vadd(dict(sy.lie(sy.lie(a, b), c)), sy.lie(b, sy.lie(a, c)), (-1) ** (pa * pb)))
    bad += l != r
    n += 1
print("[chk] [a,b] = (-1)^{|a|}{a,b}, parity |a|+1: super-antisymmetry and super-Jacobi on %d random triples: %d failures" % (n, bad))

# ---------------------------------------------------------------- g = C.1 x g',  g' = J cap ker del,  J = (w, xi_z)
def in_J(v):
    return all(any(x < 2 for x in S) or any(e[2:]) for (t, e, S) in v)
def rJ(W):
    v = relem(W)
    return {k: c for k, c in v.items() if any(x < 2 for x in k[2]) or any(k[1][2:])}
bad = n = 0
while n < 1500:
    a, b = rJ(tuple(rng.randint(-1, 1) for _ in range(5))), rJ(tuple(rng.randint(-1, 1) for _ in range(5)))
    if not (a and b):
        continue
    n += 1
    bad += not in_J(bx.bv_bracket(a, b))
print("[chk] {J,J} in J (J = ideal (w_I, xi_za) = ideal of the conormal Lagrangian {w=0, xi_z=0}): %d pairs, %d failures" % (n, bad))
bad = sum(1 for k in A0keys if sy.lie(sy.UNIT, {k: Fraction(1)}))
print("[chk] [1, A0] = 0 on %d monomials: %d failures;  1 has shifted parity %d (odd)" % (len(A0keys), bad, sy.spar(sy.UNIT)))
print("[drv] g = C.1 x g' (direct product), g' = J cap ker del = Psi^{-1}(PV^{>=1}_Delta(C^{2|3})); 1 is not in [g,g].")
nW = bad = 0
for W in itertools.product(range(-1, 3), repeat=5):
    if list(W[:2]) != sorted(W[:2], reverse=True) or list(W[2:]) != sorted(W[2:], reverse=True):
        continue
    dims = sy.g_degree_dims(W)
    pred = sum(v for k, v in dims.items() if k >= 1) + (1 if W == sy.ZERO5 else 0)
    zA0 = len(sy.Z_A0(W))
    g = len(sy.g_basis(W))
    gam = sy.c23_weight(W)
    b1g = sy.B1_dim(gam)
    nW += 1
    bad += (g != pred) or (g != zA0 - b1g)
print("[chk] dim g(W) = sum_{k>=1} dim ker(Delta|PV^k)(W) + [W=0] = dim Z(A0)(W) - dim B^1(gamma(W)),  gamma(W) = (W_z, -W_w):")
print("      %d weights W in [-1,2]^5 mod S2xS3, %d failures" % (nW, bad))

# polyvector-degree grading of the bracket
bad = n = 0
gb = {}
for W in box:
    gW = sy.g_basis(W)
    if gW:
        gb[W] = gW
keys = list(gb)
while n < 800:
    Wa, Wb = rng.choice(keys), rng.choice(keys)
    a, b = rng.choice(gb[Wa]), rng.choice(gb[Wb])
    da, db = sy.split_by_degree(a), sy.split_by_degree(b)
    for ka, va in da.items():
        for kb, vb in db.items():
            for p_ in (0, 1):
                va_p = {k: c for k, c in va.items() if len(k[2]) % 2 == p_}
                if not va_p:
                    continue
                vb_ = {k: c for k, c in vb.items()}
                r = sy.lie(va_p, vb_)
                bad += any(sy.pv_degree(k) != ka + kb - 1 for k in r)
    n += 1
print("[chk] [g_k, g_l] in g_{k+l-1} (k = C^{2|3} polyvector degree) on %d random pairs of g-basis vectors: %d failures" % (n, bad))
print("      g' = (+)_{k>=1} g_k is Z>=0-graded by k-1; g_1 = S(2|3) (divergence-free vector fields), g_{>=2} an ideal.")

print("dim g(W) by k for W in [-1,1]^5 mod S2xS3  (W: total; {k: dim}):")
line = []
for W in box:
    if list(W[:2]) != sorted(W[:2], reverse=True) or list(W[2:]) != sorted(W[2:], reverse=True):
        continue
    d = sy.g_degree_dims(W)
    g = len(sy.g_basis(W))
    line.append("%s:%d %s" % ("".join("%+d" % x for x in W), g, {k: v for k, v in sorted(d.items()) if v and k >= 1}))
for i in range(0, len(line), 3):
    print("   " + "  |  ".join(line[i:i + 3]))

# ---------------------------------------------------------------- psl(3|3)
basis = sy.psl33_basis()
vecs = [b[2] for b in basis]
print("\n[chk] psl(3|3): psi(E_PQ) = -V_PQ, V_PQ = y_P d/dy_Q (Q != 0), V_P0 = -y_P E;  y = (1 | z1, z2; theta1, theta2, theta3)")
print("      psi(identity) =", sy.psi_mat({(P, P): Fraction(1) for P in range(6)}), "; del psi(E_PQ) = str-weight:",
      [bx.delta(sy.V_elem(P, P)).get((0, sy.ZERO5, ()), 0) for P in range(6)], "(off-diagonal: 0)")
print("      34 elements, rank", la.rank(vecs), "; all in g:", all(not bx.delta(v) and sy.alpha0_is_constant(v) for v in vecs))
bad = 0
for (na, A, va) in basis:
    for (nb, B, vb) in basis:
        bad += sy.lie(va, vb) != sy.psi_mat(sy.supercomm(A, B))
print("      [psi A, psi B] = psi [A,B] on all %d basis pairs: %d mismatches" % (len(basis) ** 2, bad))
f = sy.structure_constants(vecs)
kap = sy.killing_form(vecs, f)
print("      Killing form str(ad ad): nonzero entries %d (identically zero, as for psl(n|n))" % sum(1 for v in kap.values() if v))
print("      dim [psl,psl] =", la.rank([sy.lie(a, b) for a in vecs for b in vecs]),
      "; center dim =", len(la.kernel([{(j,) + k: c for j in range(34) for k, c in sy.lie(vecs[i], vecs[j]).items()} for i in range(34)])))
print("      even|odd =", sum(1 for v in vecs if sy.spar(v) == 0), "|", sum(1 for v in vecs if sy.spar(v) == 1))
eps = {0: (0, 0, 0, 0, 0), 1: (1, 0, 0, 0, 0), 2: (0, 1, 0, 0, 0), 3: (0, 0, -1, 0, 0), 4: (0, 0, 0, -1, 0), 5: (0, 0, 0, 0, -1)}
bad = 0
for (name, M, v) in basis[:30]:
    P, Q = int(name[1]), int(name[2])
    W = bx.weight(next(iter(v)))
    bad += W != tuple(a - b for a, b in zip(eps[P], eps[Q])) or len({bx.weight(k) for k in v}) != 1
print("      T^5 weight of psi(E_PQ) = eps_P - eps_Q, eps = (0 | e_z1, e_z2; -e_w1, -e_w2, -e_w3): %d mismatches;" % bad)
print("      30 distinct root weights; Cartan = 4 divergence-free combinations of z_a d_za, theta_I d_thetaI;")
print("      the 5th torus direction (total weight, generated by the non-divergence-free Euler field) is an outer derivation.")
print("      generators (C^{2|3} polyvector of psi(E_PQ); dzA = d/dz_a, dtI = d/dtheta_I):")
for (name, M, v) in basis:
    W = bx.weight(next(iter(v))) if name[0] == 'E' else sy.ZERO5
    print("        %-4s W=%-18s label=%-16s %s" % (name, W, sy.label_of_regular(W), sy.pretty_pv(v)))

# weight-allowed components of a 2-cocycle psl x psl -> coker (from l_4(F,F,-,-): weight W_a + W_b - 2(1^5))
allowed = []
for i, (na, A, va) in enumerate(basis):
    for j, (nb, B, vb) in enumerate(basis):
        if j < i:
            continue
        Wa = bx.weight(next(iter(va))) if na[0] == 'E' else sy.ZERO5
        Wb = bx.weight(next(iter(vb))) if nb[0] == 'E' else sy.ZERO5
        Wn = tuple(a + b - 2 for a, b in zip(Wa, Wb))
        if not bx.basis(2, Wn):
            continue
        Z2, im = bx.coker_phi_t0(Wn)
        c = len(Z2) - la.rank(im)
        if c:
            allowed.append((na, nb, bx.beta_of(Wn), c))
print("[chk] pairs (x,y) in psl(3|3) with coker(W_x + W_y - 2(1^5)) != 0 (weight-allowed targets of an l_4(F,F,x,y) cocycle):")
for a in allowed:
    print("        [%s,%s] -> coker at beta=%s (dim %d)" % a)


# ---------------------------------------------------------------- gradings compatible with the flux
bad = 0
Gvals = set()
for (name, M, v) in basis:
    degs = {len(k[2]) - 1 for k in v}
    Ws = {sum(bx.weight(k)[2:]) for k in v}
    tot = {sum(bx.weight(k)) for k in v}
    bad += not (len(degs) == 1 and len(Ws) == 1 and degs.pop() == -Ws.pop())
    Gvals.add(Fraction(len(next(iter(v))[2]) - 1) + Fraction(2, 5) * tot.pop())
print("[chk] on psl(3|3): bulk degree |S|-1 = -W_w,tot for all 34 generators: %d failures" % bad)
print("[drv] gradings deg + L(W) preserved by D = t del + N{F,-} and the bracket need L(1^5) = 2 (deg F = 3, W_F = -(1^5));")
print("      on psl(3|3) deg + L(W) = (L - (0,0;1,1,1))(W); vanishing on all 30 roots (which span Z^5) forces L = (0,0;1,1,1), L(1^5) = 3.")
print("      Hence no flux-compatible grading puts psl(3|3) in degree 0; G = deg + (2/5)W_tot takes the values", sorted(Gvals))

# ---------------------------------------------------------------- U(1) sector U
Ks = [K for r in range(4) for K in itertools.combinations(range(3), r)]
def pb(m, mp):
    return m[0] * mp[1] - m[1] * mp[0], (m[0] + mp[0] - 1, m[1] + mp[1] - 1)
bad = n = 0
for m in itertools.product(range(4), repeat=2):
    for mp in itertools.product(range(4), repeat=2):
        if m == (0, 0) or mp == (0, 0):
            continue
        for K in Ks:
            for Kp in Ks:
                l = sy.lie(sy.hamiltonian_vf(m, K), sy.hamiltonian_vf(mp, Kp))
                c, mm = pb(m, mp)
                if set(K) & set(Kp) or c == 0 or mm == (0, 0):
                    ok = not l
                else:
                    ok = sy.proportional(l, sy.hamiltonian_vf(mm, tuple(sorted(K + Kp))))
                n += 1
                bad += not ok
print("\n[chk] U = span{theta^K H_h} + span{theta^K pi_0}, H_h = d2h d/dz1 - d1h d/dz2, pi_0 = d/dz1 ^ d/dz2:")
print("      [theta^K H_h, theta^K' H_h'] proportional to theta^{K u K'} H_{h,h'} ({,} Poisson for dz1^dz2): %d pairs, %d failures" % (n, bad))
bad = sum(1 for K in Ks for Kp in Ks for m in itertools.product(range(4), repeat=2)
          if sy.lie(sy.hamiltonian_vf((0, 0), K), sy.hamiltonian_vf(m, Kp)))
print("      theta^K pi_0 central in U: %d failures" % bad)
Us = [sy.hamiltonian_vf(m, K) for m in itertools.product(range(3), repeat=2) for K in Ks if (m, K) != ((0, 0), (0, 1, 2))]
def in_U(v):
    if not v:
        return True
    Ws = {bx.weight(k) for k in v}
    if len(Ws) != 1:
        return False
    beta = sy.label_of_regular(Ws.pop())
    return bool(sy.B1_dim(beta)) and sy.proportional(v, sy.u1_formula(beta))
stab = [name for (name, M, x) in basis if all(in_U(sy.lie(x, u)) for u in Us)]
print("      psl(3|3) elements preserving U:", stab)
print("      not preserving U:", [name for (name, M, x) in basis if name not in stab])
print("      psl(3|3) cap U:", [name for (name, M, x) in basis if in_U(x)])
print("time %.1fs" % (time.time() - t0))
