# B5: multi-loop obstruction classes of holomorphic Chern–Simons theory on C^3

Tags: [chk] means an exact check, with script and log in `/Users/libotao/Desktop/holomorphic-twist-program/explore/physics/anomalies/` (file names below are relative to it). [drv] means derived here. [src] means taken from a source (location given). [open] means not done.
No existing file was edited. Nothing was committed or downloaded. All arithmetic is exact (Fractions, sympy, int64 contractions).

## 0. Theory, sources, expectation

| statement | tag |
|---|---|
| twist of the 6d N=(1,0) vector multiplet = holomorphic BF theory on C^3 | [src arXiv:2110.15257 §3.1.2, eq. (12), ar5iv] |
| holomorphic BF + βγ hypers: each vertex has one B/β leg and each propagator uses one ⇒ \|I\| ≤ \|V\| ⇒ QME graphs have ≤ 1 loop ⇒ multi-loop classes of the (1,0) theory vanish identically. hCS on C^3 is **not** this twist (the brief's identification is incorrect). | [drv] |
| hCS on C^{3\|2} (6d N=(1,1)) = hCS on C^3 with 𝔤 = g⊗Λ[θ_1,θ_2]; every positive-loop bulk collision vanishes | [src 6d note S1; six-d-Q3/README §"Three-loop Ward equation"] |
| hCS on C^3 = gauge sector of the open B-model / type I topological string; one-loop anomaly ∝ ∫Tr_ad(c(∂A)^3) | [src arXiv:1905.09269 §II.E, eqs. (40)–(41)] |
| GS mechanism with type I KS: Tr_ad X^4 = (n−8)Tr X^4 + 3(Tr X^2)^2 ⇒ SO(8) on C^3; osp(n+8\|n) also works | [src 1905.09269 §III.B, eqs. (80)–(81)] |
| expectation: anomalies cancel at all loops (dimension argument, "no dimensionful counter-terms or potential gauge anomalies in the BRST cohomology") | [src 1905.09269 §II.J] |
| open–closed couplings Tr(μ⊢𝒜∧∂𝒜∧…∧∂𝒜): bilinear only for PV^{1,•} | [src 1905.09269 eqs. (16)–(20)] |

Action [src 6d note S1, odd directions dropped]:
$$S=\int_{\mathbb C^3}\Omega\wedge\Big[\tfrac12\delta_{AB}C^A\bar\partial C^B+\tfrac{i}{3!}C_{ABC}C^AC^BC^C\Big]_{(0,3)},\quad e^{S/\hbar};\qquad Q_0C^A=\tfrac i2C^A{}_{BC}C^BC^C\ [\mathrm{src\ S3}]\ \Rightarrow\ f_{ABC}:=\kappa(e_A,[e_B,e_C])=iC_{ABC}\ [\mathrm{drv}]$$
Matter [drv]:
- With a self-paired A, ½{S,S} ⊃ ⟨μ(γ,β),μ(γ,β)⟩ unless [odd,odd] ⊂ g satisfies Jacobi. So matter must be the odd part of a Lie superalgebra 𝔤 with an even invariant form.
- The colour weights of gl(N|M) and osp(m|2n) are the formal gl(t) and O(t) weights at t = sdim.

## 1. QME at order ħ^n: normalization and signs

$$|E|-1=3n,\ n=|E|-k+1\Rightarrow 2|E|=3k-4\Rightarrow k=2n+2,\ |E|=3n+1,\ 4\text{ legs}\quad[\mathrm{src\ S4}]$$
$$\sum_{n}\frac{\hbar^n}{(2n+2)!}\{\mathcal I^{2n+2}\}=\sum_n\hbar^n\frac{Q_n\mathcal I}{2n+2}=:\sum_n\hbar^nO_n=0,\qquad Q_n:=\frac{\{\mathcal I^{2n+1},\cdot\}}{(2n+1)!},\qquad \sum_{k}Q_kO_{n-k}=0$$
[drv: C^3 analogue of LLC_main.tex Lemma b-00452 (eq-2-40, eq-2-46); the quadratic identities are dimension-independent, six-d-Q3/README §1]

Ingredients:
$$Q_nO=\frac{i}{(8\pi^3)^n}\sum_{(\Gamma,r)}\frac{s_Q}{|\mathrm{Aut}(\Gamma,r)|}[\text{ordered contraction}]\,I^{S2}_\Gamma\ [\mathrm{src\ six\text{-}d\text{-}Q3\ "Ordered\ Wick\ weights"}],\quad I^{S2}=(-1)^{\frac{n(n+1)}2}K\ [\mathrm{src\ htk\ CONVENTIONS\ §3}]$$
$$\delta_{A_1}\!\cdots\delta_{A_k}\mathcal I=-iC_{A_1..A_k(E)}(C^E),\qquad \prod_{2n+2}C=(-1)^{n+1}\prod f\qquad[\mathrm{drv}]$$
Combining them:
$$Q_n\mathcal I=\frac{(-1)^{n+1+\frac{n(n+1)}2}}{(8\pi^3)^n}\sum_{(\Gamma,r)}\frac{F_Q(\Gamma,r)}{|\mathrm{Aut}(\Gamma,r)|},\quad F_Q=\rho_nF_{\rm rule},\quad \sum_r|\mathrm{Aut}(\Gamma,r)|^{-1}=\frac{2n+2}{|\mathrm{Aut}\,\Gamma|}$$
- ρ_n = −1 for every graph and every root at n = 1, 2, 3 [chk norm_check.log, norm_theta_insert.log].
- θ(2,2,3) needs current insertions for this check, because its bare kernel vanishes on Σλ = 0.

Define ω_n as the coefficient of θ_1θ_2θ_3θ_4 in O_n[Σθ_jx_je^{p_j·z}], with Σp = 0 and Δ = det(p_1,p_2,p_3):
$$\boxed{\omega_n=N_n\sum_\Gamma\frac1{|\mathrm{Aut}\,\Gamma|}\sum_{\sigma\in S_4}\mathrm{sgn}\sigma\,K_\Gamma(p_\sigma)W_\Gamma(x_\sigma),\quad N_n=\frac{(-1)^{n+1+\frac{n(n+1)}2}\rho_n}{(8\pi^3)^n}:\ N_1=\tfrac1{8\pi^3},\ N_2=-\tfrac1{(8\pi^3)^2},\ N_3=-\tfrac1{(8\pi^3)^3}}\ [\mathrm{drv};\rho\ \mathrm{chk}]$$
Here W_Γ is the orientation rule of anomaly-structure (c3), with f on the vertex words, κ^{-1} on the edges, and the external order equal to the kernel's leg order. K_Γ is the engine kernel.

| check | result |
|---|---|
| K·W(x^4) under random relabellings (root moved), edge reversals and edge orders; box + 9 n=3 graphs, kernels recomputed | invariant [chk relabel_n3.log] |
| Grassmann model sign = (−1)^{3m(m−1)/2} × rule sign, m = 2n+2, n ≤ 4, 61 graphs × 6 labellings | True [chk toy_sign.log] |
| ribbon (formal gl(t)) = explicit contraction, gl2 and gl3, 12 graphs | True [chk test_colour.log] |
| independent route from S6: Q_1𝓘 = (96π^3)^{-1}∫ε^{abc}tr(M(∂_aC)M(∂_bC)M(∂_cC)M(C)), M = −i ad ⇒ ω_1 = ΔStr(ad_x^4)/(16π^3) ⇒ N_1 = 1/(8π^3) | agrees [drv from src S6] |

Reduction to tensors:
$$K_\Gamma=c_\Gamma\Delta^n\ [\mathrm{src\ assembly/REPORT\ §2;\ chk\ kernels\_n3.log}]\ \Rightarrow\ \omega_n=24N_n\Delta^nT_n,\quad T_n=\sum_\Gamma\frac{c_\Gamma}{|\mathrm{Aut}\,\Gamma|}\{\mathrm{Sym}\ (n\ \mathrm{odd}),\ \mathrm{Alt}\ (n\ \mathrm{even})\}W_\Gamma$$
$$T_1=\tfrac1{48}\mathrm{Str}_{\rm ad}x^4\ (c_{\rm box}=-\tfrac16,\ |\mathrm{Aut}|=8),\quad \omega_1=\tfrac3{\pi^3}\Delta T_1;\qquad \omega_3=-\tfrac{3}{64\pi^9}\Delta^3T_3\quad[\mathrm{chk\ assemble\_n3.log}]$$

## 2. Three-loop data

Kernels are exact with sources e_i on the legs [chk kernels_n3.log; SP graphs = htk spred; K4-core values = nonSP/VERIFY §3.2 up to the labelling sign]. W is evaluated at x^4 in the formal gl(t) basis, p_k = tr x^k [chk assemble_n3.log].

| # | core | \|Aut\| | c_Γ | W_Γ(x^4) |
|---|---|---|---|---|
| 0 | banana | 4 | 1/155520 | −8p_1^4+32p_3p_1−24p_2^2+16tp_2p_1^2−8tp_4+8t^2p_3p_1−14t^2p_2^2−2t^3p_4 |
| 1 | banana | 4 | 1/155520 | −8p_1^4+32p_3p_1−24p_2^2+16tp_2p_1^2−8tp_4−8t^2p_2^2 |
| 2 | banana | 16 | 1/155520 | = #1 |
| 3 | K4 | 2 | −1/77760 | 4p_1^4−16p_3p_1+12p_2^2−8tp_2p_1^2+4tp_4+8t^2p_3p_1−2t^2p_2^2−2t^3p_4 |
| 4 | K4 | 2 | −1597/5443200+ζ3/4200 | w := −4p_1^4+16p_3p_1−12p_2^2+8tp_2p_1^2−4tp_4−4t^2p_2^2 |
| 5 | K4 | 2 | 0 | 0 |
| 6 | K4 | 1 | 17/226800−ζ3/16800 | w |
| 7 | K4 | 8 | 13/22400−ζ3/2100 | 2w |
| 8 | K4 | 8 | −13/22400+ζ3/2100 | w |

ζ(3) cancels, so T_3 is rational:
$$\mathrm{Sym}W_{6}=\mathrm{Sym}W_8=\mathrm{Sym}W_4=\tfrac12\mathrm{Sym}W_7;\quad \tfrac1{8400}-\tfrac1{16800}-\tfrac2{16800}+\tfrac1{16800}=0\ \Rightarrow\ T_3\in\mathbb Q$$
[chk zeta3_identity.log; the ζ3 part is also 0 by explicit contraction for sl2, sl3, so5–so9, sp4, sp6; rerun by the orchestrator 2026-11-23: zeta3_identity, check_sl23, assemble_n3 reproduce]

## 3. T_3 results

gl(t), all gl(N|M), str [chk assemble_n3.log, table_n3.log]:
$$T_1=\tfrac1{48}(2tp_4-8p_3p_1+6p_2^2),\quad T_3=\tfrac{t(t^2-6)}{103680}p_4-\big(\tfrac1{5760}+\tfrac{t^2}{34560}\big)p_2^2+\tfrac{6-t^2}{25920}p_3p_1+\tfrac{t}{8640}p_2p_1^2-\tfrac{1}{17280}p_1^4$$
Checks on the gl(t) formula:
- gl(1) (abelian): T_1 = T_3 = 0 [chk].
- so6 = sl4: the explicit so6 contraction equals the gl(t) formula at t = 4, for n = 1 and n = 3 [chk table_n3.log].

so(N) (explicit contraction, N = 5..9, degree ≤ 3 in N) [chk assemble_soN.log, table_n3.log]:
$$T_1=\tfrac{N-8}{48}\mathrm{tr}x^4+\tfrac1{16}(\mathrm{tr}x^2)^2,\qquad T_3=\tfrac{(N-8)(N-5)(N+4)}{829440}\mathrm{tr}x^4-\tfrac{N^2-11N+36}{276480}(\mathrm{tr}x^2)^2$$
- sp(2n): substitute N → −2n and tr x^4 → −tr x^4. This reproduces sp4 = (0, −1/2880) and sp6 = (77/207360, −23/46080) [chk assemble_classical.log, assemble_sp6.log].
- so8 fitted with Pf x included: the tr x^4 and Pf coefficients are both 0 [chk assemble_classical.log].

| 𝔤 | T_1 | T_3 | T_3/T_1 |
|---|---|---|---|
| sl2 | p_2^2/6 | −p_2^2/3240 | −1/540 [chk check_sl23.log] |
| sl3 | 3p_2^2/16 | −p_2^2/2560 | −1/480 [chk] |
| so8, osp(8+2m\|2m) | (tr x^2)^2/16 | −(tr x^2)^2/23040 | −1/1440 [chk; osp drv] |
| psl(N\|N) | (str x^2)^2/8 | −(str x^2)^2/5760 | −1/720 [chk] |
| gl(N\|N) | (str x^2)^2/8 − str x·str x^3/6 | −T_1/720 − (str x)^4/17280 | — [chk] |
| so(N), irreducible part | (N−8)/48 | (N−8)(N−5)(N+4)/829440 | (N−5)(N+4)/17280 |
| sl(t), irreducible part | t/24 | t(t^2−6)/103680 | (t^2−6)/4320 |

## 4. SL(3)-invariant local cohomology at ghost number 1

Setting: the translation-reduced CE complex of 𝔤[[z]], with cochains polynomial in the p's modulo Σp [src anomaly-structure VERIFY; freeness over C[∂] drv].

No coboundaries:
$$C^3_{SL_3,\deg>0}=0\ (\text{invariants of two vectors are constants})\Rightarrow B^4_{SL_3}=0;\qquad C^4_{SL_3}=\textstyle\bigoplus_n\Delta^n\otimes\{T(x_s)=\mathrm{sgn}(s)^{n+1}T\}\quad[\mathrm{drv}]$$
Consequences [drv]:
- There is no SL(3)-invariant counterterm at any loop order n ≥ 1, so the quantization, if it exists, is unique.
- [O_n] = 0 ⇔ ω_n = 0 as a cochain.

Closedness lemma:
$$d(\Delta^nT)=\pm\sum_{i<j}c^{(n)}_{ij}T([x_i,x_j],x_k,x_l,x_m),\quad c^{(n)}_{ij}=(-1)^{i+j}\det(p_k,p_l,p_m)^n$$
- For invariant T: n odd ⇒ T_{ij} lies in Z ≅ Λ^2(std) (dim 6, S_5-irreducible); n even ⇒ (−1)^{i+j}T_{ij} lies in Y ≅ [3,2] (dim 5, irreducible).
- c^{(1)} is a gradient; c^{(n)} for n = 2..6 is not; δ^n ≠ φ_i + φ_j for n = 2, 4, 6 [chk lemma_gradient.log].

Result:
$$\Rightarrow Z^4(\Delta^1)=\mathrm{Sym}^4(\mathfrak g^*)^{\mathfrak g};\qquad Z^4(\Delta^{n\ge2})\cap\mathrm{inv}=\mathrm{Sym}^4/\Lambda^4\big((\mathfrak g/[\mathfrak g,\mathfrak g])^*\big)\quad[\mathrm{drv}]$$
Direct check on sl2, sl3, gl2, gl3 with T ∈ {Str4, (tr xx)^2, Π tr x_i, alt(tr x_1 tr x_2[x_3,x_4])}, n = 0..5: closed exactly in the cases predicted [chk ce_check.log].

hCS weights vanish on central legs [src anomaly-structure (b2)], so reachable classes live on 𝔤/([𝔤,𝔤]+𝔷) [drv]:

| 𝔤 | classes at n ≥ 2 |
|---|---|
| perfect: simple, g⊗Λ[θ_1,θ_2], psl(N\|N), osp | none |
| reductive | none |
| gl(N\|N) | Δ^n(str x)^4 at odd n only |

## 5. Classes after one-loop cancellation

**(a) Cancellation inside the gauge complex (T_1 = 0).** By induction with Bianchi, dO_n = 0 ⇒ O_n = 0 for all n when 𝔤 is perfect [drv].
- The only realization found is g⊗Λ[θ_1,θ_2]: the handle element is 0, so W_Γ = 0 for n ≥ 1 [drv; src six-d-Q3 §1].
- Str_ad x^4 ≢ 0 for gl(m|n), psl, osp and D(2,1;α) [drv].

**(b) Pure hCS (T_1 ≠ 0).** [O_3] is undefined. dω_3 ≠ 0 whenever T_3 ≠ 0 (lemma), consistent with dO_3 = −Q_2O_1. Q_2O_1 itself was not computed [open].

**(c) Green–Schwarz mechanism (type I KS).** In the open-string sector:
$$T_1^{\rm irr}=0\Rightarrow T_3^{\rm irr}=0:\ \ so(N)\propto(N-8)(N-5)(N+4),\ \ sl(t)\propto t(t^2-6),\ \ so8:\ \mathrm{tr}x^4,\mathrm{Pf}\to0\quad[\mathrm{chk}]$$
So the three-loop hCS obstruction has GS-factorized colour exactly when the one-loop obstruction does (so8, psl(N|N), and trivially sl2, sl3 and the exceptional algebras).

Sub-blob colours for so8 (Sym parts) [chk theta_sym.log]:
- box: −3(tr x^2)^2; θ(1,3,3): 9(tr x^2)^2; θ(2,2,3): 0; tree: 0.
- In formal gl(t): W_θ(1,3,3) = −t·W_box.

Consequences:
- KS-exchange graphs with only bilinear couplings (eq. 17) give only (tr x^2)^2 for so8 [drv].
- Couplings with three or more open fields (eqs. 16, 19) can give single traces [src]. Their three-loop contribution, and the cancellation of the remainder −T_1/1440, are [open].

**(d) gl(N|N).** The genuine class receives
$$\omega_3\supset\frac{\Delta^3\prod_i\mathrm{str}\,x_i}{368640\,\pi^9}\quad[\mathrm{chk\ table\_n3.log}]$$
- It is absent at one loop: the box has 2 faces and T_1 has no p_1^4 term [chk].
- A two-annulus KS exchange carries the same colour [drv]; the total coefficient is [open].

## 6. Comparison with the expectation

| expectation | result | tag |
|---|---|---|
| 6d gauge anomalies are one-loop exact | (1,0) twist is BF-type: there are no multi-loop QME graphs | [drv, src] |
| one-loop anomaly cancelled ⇒ higher classes vanish | perfect 𝔤: H^1_{SL3}(n ≥ 2) = 0, and O_1 = 0 ⇒ O_n = 0 for all n | [drv, chk] |
| GS remains sufficient beyond one loop (1905.09269 §II.J) | T_3^irr vanishes whenever T_1^irr does | [chk] |
| no new classes | exception: gl(N\|N) Δ^3(str x)^4 | [drv, chk] |

## 7. Four loops

$$\omega_4=24N_4\Delta^4\sum_\Gamma c_\Gamma|\mathrm{Aut}|^{-1}\mathrm{Alt}\,W_\Gamma$$
Alt W_Γ = 0 graph by graph [chk even_orders.log]:
- formal gl(t), both thetas and all 50 four-loop graphs (37 with c_Γ ≠ 0);
- explicit so5 and sp4 on the thetas and the 11 graphs with banana or K4 core.

So ω_4 ≡ 0 at the cochain level for gl(t) and for every reductive g. The mechanism is H^4(𝔰) = 0 plus decoupling of central legs [src (b2)], now confirmed on the actual four-loop graph list [chk].

## 8. Side result

The orientation rule of anomaly-structure (c3) agrees with the nilpotency-fixed six-d-Q3 rule on the relative sign of θ(1,3,3) and θ(2,2,3) at H = 3 (ρ = −1 at all roots, with insertions) [chk norm_theta_insert.log]. This was listed as untested in (c3). Carrying it over to the 5d CS case (D = 3) is [drv], not run.

## 9. Not done / failed

- Q_2O_1 and the identity dO_3 = −Q_2O_1 [open].
- Three-loop KS-exchange graphs; the full class of the GS-cancelled theories [open].
- Secondary targets (KS on C^3, twistor type I) were not attempted.
- so(N) with N ≥ 10 and the so8, prism-core and K33-core four-loop alternations were stopped: the int64 einsum takes about 4 s per contraction. The so(N) polynomials rest on N = 5..9 plus the sp4, sp6 and so6 = sl4 cross-checks.
- One WebFetch of arXiv:2110.15257 returned a PDF binary, saved automatically by the tool and not read; the ar5iv HTML was used instead.
- arXiv:1810.06534 was read at abstract level only.

## Files

All under `/Users/libotao/Desktop/holomorphic-twist-program/explore/physics/anomalies/`:
- Kernels and colour: kernels_n3.py/.json/.log, colour.py, test_colour.py/.log.
- Signs and normalization: relabel_n3.py/.log, toy_sign.py/.log (+ htk_free_relabel.py), norm_check.py/.log, norm_theta_insert.py/.log.
- Assembly: assemble_n3.py/.json/.log, assemble_classical.py/.log, assemble_soN.py/.log (so5–so9; run stopped at so10), assemble_spN.py + assemble_sp6.log, table_n3.py/.log, check_sl23.py/.log, zeta3_identity.py/.log, theta_sym.py/.log.
- Cohomology: lemma_gradient.py/.log, ce_check.py/.log.
- Even orders: even_orders.py/.log.
