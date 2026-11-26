"""Part (c) for gl_N, general N (trace basis, formal_gl.py), and class structure for sl3, gl3.

omega = alt(K'_133 W_133) + alt(K'_223 W_223),  K' = (60/pi) K  (|Aut| = 4 common).
 1. combined colour tensor W_tot = -3 W_133 + 2 W_223 multiplying [p1 p2][p3 p4] (trace basis);
    its restriction to sl_N (terms without tr(X_i) factors); for sl2 in the kappa kappa basis.
 2. formal closedness d omega = 0 (identity in the trace basis, all N).
 3. exactness: ranks of d(C^3) (spanning set alt(tau_3 (x) mu)) and of [d(C^3) | omega(N)] at
    N = 2..10 (formal trace basis: faithful on gl_N for N >= 4).
 4. ranks of [d(C3) | omega_133 | omega_223] (explicit components) for sl3, gl3 and formally.
 5. dimension of the closed invariant 4-cochains and of H^4 at weight (2,2) in the formal basis
    (rank of d: C^4 -> C^5 modulo two primes; faithful for N >= 5).
"""
import sys, time, itertools
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import htk
from htk import theta
import lie, cochain as cc, formal_gl as fg
from c_obstruction import kernel_poly, trace_tensors

H = 2
mom = cc.Momenta(6, H)
P = mom.p


def br(a, b):
    return a[0] * b[1] - a[1] * b[0]


def vecs(cochains, keys):
    return np.array([[c.get(k, 0) for c in cochains] for k in keys], dtype=np.int64)


def rank_mod_p(A, p):
    M = np.array(A, dtype=np.int64) % p
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        if r == rows:
            break
        piv = np.nonzero(M[r:, c])[0]
        if piv.size == 0:
            continue
        i = r + piv[0]
        M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        nz = np.nonzero(M[:, c])[0]
        for j in nz:
            if j != r:
                M[j] = (M[j] - M[j, c] * M[r]) % p
        r += 1
    return r


t00 = time.time()
Wf, Kp = {}, {}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, ext, K = kernel_poly(pqr)
    Kp[pqr] = K
    Wf[pqr] = fg.colour_weight_formal(G.nv, G.edges, ext)
    print("theta%s: (60/pi) K =" % (pqr,), sp.factor(K))
same = sp.expand(Kp[(1, 3, 3)] / -3 - Kp[(2, 2, 3)] / 2) == 0
print("(60/pi) K_133 = -3 [p1p2][p3p4], (60/pi) K_223 = 2 [p1p2][p3p4]:", same,
      sp.expand(Kp[(1, 3, 3)] + 3 * br(P[0], P[1]) * br(P[2], P[3])) == 0)

print("\n=== 1. combined colour tensor W_tot = -3 W_133 + 2 W_223 (trace basis, {N power: coeff}) ===")
Wtot = {}
for pqr, c in (((1, 3, 3), -3), ((2, 2, 3), 2)):
    for tau, d in Wf[pqr].items():
        fg.add(Wtot, tau, d, c)
Wtot = fg.clean(Wtot)
for tau in sorted(Wtot, key=lambda t: (len(t), t)):
    print("  ", " ".join("tr(%s)" % "".join(str(i + 1) for i in c) for c in tau), Wtot[tau])
print("restriction to sl_N (drop terms containing tr(X_i)):")
Wsl = {tau: d for tau, d in Wtot.items() if all(len(c) > 1 for c in tau)}
for tau in sorted(Wsl, key=lambda t: (len(t), t)):
    print("  ", " ".join("tr(%s)" % "".join(str(i + 1) for i in c) for c in tau), Wsl[tau])
# sl2 in the kappa kappa basis
g2 = lie.sl(2)
T2 = fg.tensor_from_formal(Wtot, g2, 4)
bas = [fg.tau_tensor(((0, 1), (2, 3)), g2, 4), fg.tau_tensor(((0, 2), (1, 3)), g2, 4), fg.tau_tensor(((0, 3), (1, 2)), g2, 4)]
x = sp.symbols('x0:3')
eqs = sp.Matrix([sum(x[i] * int(b.reshape(-1)[j]) for i, b in enumerate(bas)) - int(T2.reshape(-1)[j]) for j in range(81)])
sol = sp.solve(list(eqs), x, dict=True)
print("sl2: W_tot = c12 k12 k34 + c13 k13 k24 + c14 k14 k23, k_ij = tr(X_i X_j):", sol)

print("\n=== 2. formal closedness ===")
om = {}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    om[pqr] = fg.alternate_formal(mom, Kp[pqr], Wf[pqr], 4)
omega = {}
for pqr in om:
    for k_, d in om[pqr].items():
        fg.add(omega, k_, d)
omega = fg.clean(omega)
d133 = fg.d_formal(mom, om[(1, 3, 3)], 4)
d223 = fg.d_formal(mom, om[(2, 2, 3)], 4)
print("d omega_133 = 0 formally:", d133 == {}, "; d omega_223 = 0 formally:", d223 == {})
print("N powers in omega:", sorted({p for d in omega.values() for p in d}))

print("\n=== 3. exactness for gl_N (formal trace basis) ===")
monos3 = cc.all_monos(mom, 3, (2, 2))
monos4 = cc.all_monos(mom, 4, (2, 2))
taus3 = fg.all_taus(3)
cols3 = []
for tau in taus3:
    for mu in monos3:
        eta = fg.alternate_formal(mom, mom.mono_expr(mu, 3, [P[0], P[1]]), {tau: {0: 1}}, 3)
        cols3.append(fg.at_N(fg.d_formal(mom, eta, 3), 0) if eta else {})
keys4 = sorted({k for c in cols3 for k in c} | set(omega), key=str)
A = vecs(cols3, keys4)
rA = cc.exact_rank(A)
print("spanning set of C^3: %d vectors; rank d(C^3) = %d" % (len(cols3), rA))
for N in range(2, 11):
    v = vecs([fg.at_N(omega, N)], keys4)
    r = cc.exact_rank(np.concatenate([A, v], axis=1))
    v1 = vecs([fg.at_N(om[(1, 3, 3)], N)], keys4)
    v2 = vecs([fg.at_N(om[(2, 2, 3)], N)], keys4)
    r12 = cc.exact_rank(np.concatenate([A, v1, v2], axis=1))
    print("  N = %2d: rank [d(C3) | omega] = %d  -> exact: %s;  rank [d(C3) | omega_133 | omega_223] = %d%s" % (
        N, r, r == rA, r12, "" if N >= 4 else "   (formal basis not faithful for N < 4)"), flush=True)

print("\n=== 4. explicit components: class structure of omega_133, omega_223 ===")
for g in [lie.sl(3), lie.gl(3)]:
    ws = {}
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        G, ext, K = kernel_poly(pqr)
        W, sc = cc.colour_weight(g, G.nv, G.edges, ext)
        ws[pqr] = cc.alternate(mom, K, W, 4)
    cols = []
    for T in [T for s, T in trace_tensors(g, 3)] + [g.fl]:
        for mu in monos3:
            eta = cc.alternate(mom, mom.mono_expr(mu, 3, [P[0], P[1]]), T, 3)
            cols.append(cc.flatten(cc.d_CE(mom, g, eta, 3), monos4, g.dim, 4) if eta else
                        np.zeros(len(monos4) * g.dim ** 4, dtype=np.int64))
    Ae = np.stack(cols, axis=1)
    v1 = cc.flatten(ws[(1, 3, 3)], monos4, g.dim, 4)[:, None]
    v2 = cc.flatten(ws[(2, 2, 3)], monos4, g.dim, 4)[:, None]
    r0 = cc.exact_rank(Ae)
    print("%s: rank d(C3) = %d, [.|w133] = %d, [.|w223] = %d, [.|w133|w223] = %d, [.|w133+w223] = %d, [.|w133-w223] = %d" % (
        g.name, r0, cc.exact_rank(np.concatenate([Ae, v1], 1)), cc.exact_rank(np.concatenate([Ae, v2], 1)),
        cc.exact_rank(np.concatenate([Ae, v1, v2], 1)), cc.exact_rank(np.concatenate([Ae, v1 + v2], 1)),
        cc.exact_rank(np.concatenate([Ae, v1 - v2], 1))), flush=True)

print("\n=== 5. invariant cohomology at weight (2,2), formal trace basis (N >= 5) ===")
taus4 = fg.all_taus(4)
c4 = []
for tau in taus4:
    for mu in monos4:
        w = fg.alternate_formal(mom, mom.mono_expr(mu, 4, [P[0], P[1], P[2]]), {tau: {0: 1}}, 4)
        c4.append(fg.at_N(w, 0))
keysC4 = sorted({k for c in c4 for k in c}, key=str)
C4 = vecs(c4, keysC4)
p1, p2 = 2147483629, 2147483587
dimC4 = (rank_mod_p(C4, p1), rank_mod_p(C4, p2))
print("dim C^4 (invariant, weight (2,2)): rank of %d spanning vectors mod p1, p2 = %s" % (len(c4), dimC4), flush=True)
t0 = time.time()
d4 = [fg.at_N(fg.d_formal(mom, {k: {0: v} for k, v in c.items()}, 4), 0) for c in c4]
keys5 = sorted({k for c in d4 for k in c}, key=str)
D4 = vecs(d4, keys5)
rd4 = (rank_mod_p(D4, p1), rank_mod_p(D4, p2))
print("rank d: C^4 -> C^5 mod p1, p2 = %s  (%.0fs)" % (rd4, time.time() - t0))
B = vecs(cols3, keysC4)
rB = (rank_mod_p(B, p1), rank_mod_p(B, p2))
print("dim Z^4 = %d, dim B^4 = %s, dim H^4_(2,2) = %d  (mod-p ranks; B^4 exact rank %d)" % (
    dimC4[0] - rd4[0], rB, dimC4[0] - rd4[0] - rA, rA))
print("total %.0fs" % (time.time() - t00))
