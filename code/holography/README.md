# holography

Exact (fractions.Fraction) tools for the twisted-holography dictionary: holomorphic twist of N=4 U(N) SYM vs BCOV on X = C^5 \ C^2 with flux N F (Costello–Li 1606.00365 §9).

Tests: `python3 code/holography/tests/test_bulk.py`, `test_cyclic.py`, `test_charges.py` (also pytest). `test_bulk.py::test_t0_sixty_slots` reruns the t^0 60-slot check against `explore/twisted-holography-N4/code/single_particle.py`.

This README covers linalg, bulk_x, derham, cyclic, charges (used by `explore/twisted-holography-N4/dictionary/run_A*.py`).

## linalg.py
Sparse vectors `{key: Fraction}`; `Echelon` (incremental row echelon, optional combination tracking), `rank`, `kernel(images)`, `rel_rank(U, V) = dim(U+V) - dim V`.

## bulk_x.py
- Coordinates (z1,z2,w1,w2,w3) = 0..4; xi_x = d/dx odd.
- Monomial key `(typ, e, S)`: typ 0 = A0 = C[z,w] (x) Lambda[xi] (e >= 0); typ 2 = A2 = C[z] (x) L (x) Lambda[xi], L = span{w^(-a-1)} (w-entries of e <= -1). S sorted; monomial z^m w^e xi_{S[0]} xi_{S[1]} ...
- Weight W = e - ind(S).
- Products: A0.A0 -> A0; A0.A2 -> A2, zero once a w-exponent is >= 0; A2.A2 = 0. Parity |S| mod 2.
- `delta` = sum_x d/dx d/dxi_x, with d/dxi from the left and d/dw w^(-a-1) = -(a+1) w^(-a-2).
- `bv_bracket(a,b)` = del(ab) - del(a) b - (-1)^{|a|} a del(b), equal to `schouten(a,b)` = sum_x (d_xi a)(d_x b) + (-1)^{|a|}(d_x a)(d_xi b).
- Symmetry {a,b} = (-1)^{|a||b|}{b,a}; del{a,b} = -{del a,b} - (-1)^{|a|}{a,del b}.
- `F` = (w1w2w3)^(-1) xi_z1 xi_z2; `XI_TOP` = xi_z1...xi_w3.
- T keys `(k, typ, e, S)` = t^k * monomial; `D(v, N)` = t del + N{F,-}.
- Grading: What = W - (1,...,1) on A0, What = W on A2; h = k - 1 on A0, h = k on A2. D maps T(What,h) -> T(What,h+1), and T(What,h) = t^h A2_What (+) t^(h+1) A0_(What+1) is finite-dimensional.
- `class_cohomology(What, h)` returns dim H, the Laurent part (image of H(A2[[t]])), the regular part (image in H(A0[[t]])), and the same two numbers from phi_* (LES).
- `ker_phi_t0`, `coker_phi_t0`, `ker_phi_formula`, `beta_of`.
- Pairing candidates `fischer_res`, `top_product`; `pairing_data(Wn)` (A2 of the dictionary run).

## derham.py
- Polynomial forms on C^{2|3}, key `(m, Sz, J, a)` = z^m dz^Sz theta^J dtheta^a. dz is odd and dtheta even.
- `d_form`; `B_dim(beta, p)` = rank d: Omega^{p-1} -> Omega^p at weight beta.
- `Phi`: A2 -> Omega, z^m w^(-a-1) xi_z^Sz xi_w^J -> eps(Sz)(-1)^|a|/a! z^m dz^{Sz^c} theta^J dtheta^a, with eps({0,1}) = eps({1}) = 1 and eps({0}) = eps({}) = -1.
- Phi o del = d o Phi and Phi(F) = 1.

## cyclic.py
- A = C[z1,z2] (x) Lambda[theta1..3], monomial `(m1, m2, mask)`.
- Hochschild `b`, cyclic sign s(x) = (-1)^{n + |a_n|(|a_0|+...+|a_{n-1}|)}, Connes complex C^lambda = C/(1-t).
- `HC(beta, n, reduced)`: reduced uses the augmentation ideal.
- Cyclic-word model: letters C_{m,J} of parity 1+|J|. Q C_{M,K} = (-1)^{|K|} sum binom(M,m) eps(J,J') (-1)^{|J'|(1+|J|)} C_{m,J} C_{m',J'}, from Q C = C.C.
- The relative complex drops c = C_{0,{}}.
- `CE_single_trace(beta, Y, relative)`.

## charges.py
- `llc_charges((theta), (n+,n-), Y)` returns (E; R_I; J_L, J_R) of LLC_main.tex eq-E-charges.
- `kmmr_monomial` gives (-1)^{Y+Theta} t^{2Theta+3n} y^{n+-n-} v^{theta1-theta2} w^{theta2-theta3}, KMMR (4.2) conventions.
- `kmmr_f` is KMMR (4.2); `kmmr_Isp` is KMMR (4.14).

## Finite N (Budzik-Gaiotto meson complex; used by explore/twisted-holography-N4/finiteN*/)
- `bg_words.py`: letter algebras (cobar word algebras A^!), Q on letters and words, graded-cyclic classes.
- `bg_complex.py`, `bg_complex_k.py`: meson complex C(gl'_{k|k}[A^!] | b_{k|k}; Ber^N), k = 1 and general k (Python, exact or mod p).
- `bg_fast.py` + `native/bgk.cpp`: same complex, C++ kernel, ranks mod p; reduced b-set {P0_{x,x+1}, MM0_{x,x+1}, R0_{1,k}}; block elimination by (T-monomial, multiset of words). `homology(alg, lam, k, N, p, ghs, mode)`.
- `bg_sizes.py`: dim V_{charge,g} from generating functions (no enumeration); gauge-theory cochain counts.
- `uN_q0.py`: finite-N H(Q_0) of the relative U(N)/SU(N) complex per sector.
- `rank_mod.py`: sparse ranks over F_p and Q.
- `giant.py`: 1/16-BPS U(N)/SU(N) index and the giant-graviton expansion.
- `census_io.py`: read-only access to the LLC census data.

## Bulk symmetry algebra (used by explore/twisted-holography-N4/symmetry/)
- `symmetry.py`: g = ker phi_*|t^0; Psi: A0 -> PV(C^{2|3}); iota, L_pi on Omega(C^{2|3}); psl(3|3) generators V_PQ. Tests: `tests/test_symmetry.py`.
