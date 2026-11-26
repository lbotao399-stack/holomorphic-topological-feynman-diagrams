"""gl_N (kappa = tr): order-hbar^3 classes.  N = 7 (trace basis faithful for k <= 5 < N; formal coefficients
evaluated at N), cross-checked at N = 11 for the dimensions.
 (1) H^3_(2,2) = Z^3_(2,2) (C^2_(2,2) = 0): dimension, C-parity, SL2-invariant part.
 (2) H^4_(3,3): dimension; C-odd part (subcomplex) ; SL2-invariant part.
 (3) phi_T(h) = [{T,h}] in H^4_(3,3) for T in {D, S1, S0}, h in Z^3_(2,2): rank of the span of all phi_T(h_i),
     kernel of T -> phi_T, rank of phi_J for J = j_M (Moyal), J = pi^*D_sl, random J; explorer's relations
     [{D,h1}] = -3e1, [{S1,h1}] = -(3N/2)e1 - 2v, [{S0,h1}] = -Nv, [{D,h2}] = 2u, [{S1,h2}] = Nu + 2w, [{S0,h2}] = Nw
     imply (basis-independent):  N(2 phi_S1 - N phi_D) = 4 phi_S0 on H^3, phi_{pi^*D_sl} = 0, rank{phi_D(h_i), phi_S0(h_i)} = 4.
Exact ranks over Q (sympy sparse rref)."""
import pickle, time, random, sys
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
from vlib import *
import vla

t0 = time.time()
SL = False
NN = int(sys.argv[1]) if len(sys.argv) > 1 else 7


def rows(cs, N):
    return vla.to_rows(cs, N)


def basis(cs, N):
    idx = vla.indep_rows(rows(cs, N))
    return [cs[i] for i in idx]


def nullspace_combos(B, images, N):
    """kernel of the linear map B[i] -> images[i] (coordinates at N): list of cochains sum c_i B[i]."""
    r = rows(images, N)
    cols = {}
    for x in r:
        for kk in x:
            cols.setdefault(kk, len(cols))
    M = sp.zeros(max(1, len(cols)), len(B))
    for i, x in enumerate(r):
        for kk, v in x.items():
            M[cols[kk], i] = sp.Rational(v.numerator, v.denominator)
    ns = DomainMatrix.from_Matrix(M).convert_to(QQ).nullspace().to_Matrix()
    out = []
    for j in range(ns.rows):
        z = {}
        for i in range(len(B)):
            if ns[j, i] != 0:
                z = lin((1, z), (Fr(int(ns[j, i].p), int(ns[j, i].q)), at_N(B[i], N)))
        out.append(z)
    return out


N = NN
print("gl_N at N = %d (trace basis)" % N, flush=True)
# (1) H^3_(2,2)
C2_22 = span_set(2, (2, 2), SL)
print("C^2_(2,2) spanning set nonzero elements: %d -> dim %d" % (len(C2_22), vla.rankQ(rows(C2_22, N)) if C2_22 else 0))
B3_22 = basis(span_set(3, (2, 2), SL), N)
dB3_22 = [d(w, 3, SL) for w in B3_22]
r3 = vla.rankQ(rows(dB3_22, N))
print("C^3_(2,2): dim %d, rank d3 %d -> Z^3 = H^3 = %d" % (len(B3_22), r3, len(B3_22) - r3), flush=True)
Z3 = nullspace_combos(B3_22, dB3_22, N)
print("  kernel vectors exactly closed: %s" % all(is_zero(d(z, 3, SL)) for z in Z3))
print("  C on Z^3_(2,2): C h = -h for all kernel vectors: %s" % all(is_zero(lin((1, Cop(z, 3)), (1, z))) for z in Z3))
for z in Z3:
    print("  h: %d terms; E h = 0 (SL2-invariant momentum polynomials): %s; perms: %s" %
          (len(z), is_zero(sl2_E(z, 3)), sorted(set(perm_str(p) for (p, m) in z))))
C3inv22 = basis(span_set(3, (2, 2), SL, sl2inv=True), N)
print("C^3_inv(2,2): dim %d; rank d on it: %d -> (H^3_(2,2))^{SL2} = %d" %
      (len(C3inv22), vla.rankQ(rows([d(w, 3, SL) for w in C3inv22], N)), len(C3inv22) - vla.rankQ(rows([d(w, 3, SL) for w in C3inv22], N))))
u2f = cubic('f', 2, SL)
print("  C^3_inv(2,2) = C u2f: rank[C3inv|u2f] = %d" % vla.rankQ(rows(C3inv22 + [u2f], N)), flush=True)

# (2) H^4_(3,3)
B3_33 = basis(span_set(3, (3, 3), SL), N)
dB3_33 = [d(w, 3, SL) for w in B3_33]
rd3 = vla.rankQ(rows(dB3_33, N))
print("C^3_(3,3): dim %d, rank d3 %d" % (len(B3_33), rd3), flush=True)
t1 = time.time()
B4_33 = basis(span_set(4, (3, 3), SL), N)
print("C^4_(3,3): dim %d (%.0fs)" % (len(B4_33), time.time() - t1), flush=True)
t1 = time.time()
dB4_33 = [d(w, 4, SL) for w in B4_33]
print("  d of the C^4 basis computed (%.0fs)" % (time.time() - t1), flush=True)
rd4 = vla.rankQ(rows(dB4_33, N))
rd4p = vla.rankp(rows(dB4_33, N))
print("rank d4 = %d (mod p %d);  H^4_(3,3) = %d - %d - %d = %d" % (rd4, rd4p, len(B4_33), rd4, rd3, len(B4_33) - rd4 - rd3), flush=True)
# C-odd part of H^4: C-odd projections
Cm3 = [scal(Fr(1, 2), lin((1, w), (-1, Cop(w, 3)))) for w in B3_33]
Cm4 = [scal(Fr(1, 2), lin((1, w), (-1, Cop(w, 4)))) for w in B4_33]
dCm4 = [scal(Fr(1, 2), lin((1, x), (-1, Cop(x, 5)))) for x in dB4_33]     # d commutes with C
dCm3 = [scal(Fr(1, 2), lin((1, x), (-1, Cop(x, 4)))) for x in dB3_33]
print("  check d C = C d on 3 basis elements: %s" % all(is_zero(lin((1, d(Cop(w, 4), 4, SL)), (-1, Cop(d(w, 4, SL), 5)))) for w in B4_33[:3]))
c4m = vla.rankQ(rows(Cm4, N)); r4m = vla.rankQ(rows(dCm4, N)); r3m = vla.rankQ(rows(dCm3, N))
print("C-odd subcomplex (3,3): dim C^4_- %d, rank d4|_- %d, rank d3|_- %d -> H^4_- = %d" % (c4m, r4m, r3m, c4m - r4m - r3m), flush=True)
# SL2-invariant part
C3i = basis(span_set(3, (3, 3), SL, sl2inv=True), N)
C4i = basis(span_set(4, (3, 3), SL, sl2inv=True), N)
r3i = vla.rankQ(rows([d(w, 3, SL) for w in C3i], N))
r4i = vla.rankQ(rows([d(w, 4, SL) for w in C4i], N))
print("SL2-invariant subcomplex (3,3): dim C^3 %d, C^4 %d, rank d3 %d, rank d4 %d -> H^4_inv = %d" %
      (len(C3i), len(C4i), r3i, r4i, len(C4i) - r4i - r3i), flush=True)

# (3) phi maps
D = cubic('d', 1, SL); S1 = cubic('s1', 1, SL); S0 = cubic('s0', 1, SL)
Bd = rows(dB3_33, N)
rB = vla.rankQ(Bd)
phi = {}
for nm, T in [('D', D), ('S1', S1), ('S0', S0)]:
    for i, h in enumerate(Z3):
        c = bracket(T, 3, h, 3, SL)
        phi[(nm, i)] = c
        print("  {%s,h%d}: closed %s, C-even %s, class nonzero: %s" %
              (nm, i + 1, is_zero(d(c, 4, SL)), is_zero(lin((1, Cop(c, 4)), (-1, c))),
               vla.rankQ(Bd + rows([c], N)) > rB), flush=True)


def cls_rank(cs):
    return vla.rankQ(Bd + rows(cs, N)) - rB


allphi = [phi[k] for k in sorted(phi)]
print("dim span of all classes [{T,h_i}] = %d" % cls_rank(allphi))
print("rank{[{D,h1}], [{S0,h1}], [{D,h2}], [{S0,h2}]} = %d" % cls_rank([phi[('D', 0)], phi[('S0', 0)], phi[('D', 1)], phi[('S0', 1)]]))
Nn = Fr(N)
for i in range(2):
    rel = lin((Nn * 2, phi[('S1', i)]), (-Nn * Nn, phi[('D', i)]), (-4, phi[('S0', i)]))
    print("  N(2 phi_S1 - N phi_D) - 4 phi_S0 on h%d: class zero: %s" % (i + 1, cls_rank([rel]) == 0 if rel else True))
pistar = lin((1, D), (-Fr(2) / Nn, S1), (Fr(4) / (Nn * Nn), S0))
jM = scal(Fr(1, 2), D)
fam = {'pi^*D_sl': pistar, 'j_M = D/2': jM, 'S1': S1, 'S0': S0, 'D': D,
       's j_M + t pi^*D_sl (s=1,t=3)': lin((1, jM), (3, pistar)),
       'branch beta=0: pi^*D_sl + 5 S0': lin((1, pistar), (5, S0))}
random.seed(3)
fam['random (a,b,c)'] = lin((random.randint(1, 9), D), (random.randint(1, 9), S1), (random.randint(1, 9), S0))
for nm, J in fam.items():
    imgs = [bracket(J, 3, h, 3, SL) for h in Z3]
    print("  J = %s: rank of h -> [{J,h}] on H^3_(2,2) = %d" % (nm, cls_rank(imgs)), flush=True)
# SL2-invariant I_2: the theta piece and the kappa_2 piece are exact (formal identities, v3); here at N:
A = pickle.load(open('data/assembled_gl.pkl', 'rb'))
for nm, c in [("theta'(d)", A['th']['d']), ("{D,u2f}", bracket(D, 3, u2f, 3, SL)), ("{S1,u2f}", bracket(S1, 3, u2f, 3, SL)), ("{S0,u2f}", bracket(S0, 3, u2f, 3, SL))]:
    print("  %s: exact at N=%d: %s" % (nm, N, cls_rank([c]) == 0 if not is_zero(c) else "zero cochain"))
print("time %.1fs" % (time.time() - t0))
