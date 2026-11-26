# VERIFY: twoloop-universal

Scripts and logs: `verify/` (vcore.py, v1_pointwise.py, v1b_unsorted.py, v2_values.py, v3_sources.py, v4_numerics.py, v5_formscheck.py, v6_expansion.py, v7_oracle_rerun.py, v8_polarization.py, *.log; S3_pointwise_run.log). No explorer script was imported or run. `vcore.py` is written from the definitions: graphs, paths and C3 built afresh; currents from the full KCL + KVL system with a BFS cycle basis (not C3, not the grounded Laplacian); forms evaluated on vectors,

$$D(t)=\frac{\Omega(v_1,\dots,v_{N-1})}{\det[t;v_1;\dots;v_{N-1}]},\qquad \Omega=U^{-H}\rho_2^T\wedge\prod_{i}\psi^i_1\psi^i_2,\quad \rho_2=U^{-1/2}A_{12},\ A_{12}=-\tfrac12\textstyle\sum_{ab}M^{-1}_{ab}\,dM_{1a}\wedge dM_{b2},$$

with two independent random rational vector sets required to agree (basicness of $\Omega$). `htk` is used only as an oracle (`kernel_exact`, `graph.theta`, `cycle_matrix`, one re-run of `kernel_numeric`).

REPORT.md does not exist in the folder; the proof text examined is the claim list.

## Verdicts

| claim | verdict |
|---|---|
| C0-setting | confirmed |
| C1-derivation | confirmed |
| C2-densities | confirmed |
| C3-integral-representation | confirmed |
| C4-displacements | confirmed |
| C5-T1-closed | confirmed |
| C6-T0-closed | confirmed |
| C7-basis-sign | confirmed |
| C8-interior-root | confirmed |
| C9-determinant-expansion | confirmed |
| C10-S1 | confirmed |
| C11-S2 | confirmed |
| C12-numerics | confirmed (own values; 2 of 6 oracle runs re-run) |
| C13-formscheck-value-sign | confirmed |
| C14-formscheck-cause | confirmed |

## Checks (printed outputs)

### C0, C1, C2, C4, C7, C8 pointwise (v1_pointwise.log, v1b_unsorted.log), exact rationals, all 32 graphs (T=0: H=1..6; T=1: H=0..5)

Per graph: 3 random rational $t$, sources at all non-root vertices (root $N'$); 2 random transforms (vertex relabelling, edge permutation, random edge reversals, root at a random vertex):

$$D^{\rm bf}(t)=\mathrm{sgn}(\pi_{\rm std})\times\begin{cases}(-1)^{H(H+1)/2+1}R_AR_BR_C\,Q(L(t))/U^{H+2}&T=0\\(-1)^{H(H+1)/2}D_\theta/(2U^{H+3/2})&T=1\end{cases}$$

`C2 density (root N'): True`, `C2 relabel/reverse/permute/any root: True` on all 32 lines; 53 graphs with unsorted path orders (B any of the three paths): `brute = closed density and sgn det = (-1)^(q+1): True`.

C0: $s_e\,y_{P_j}=w_PL-\hat m_P+Q^P_j$ exactly: `C0 currents: True` (32/32).

C1: (1) $\varphi^i_P(v)=V^i\,dR_P(v)/R_P+R_P\,d\hat m^i_P(v)$, $V=L/\sum_P R_P^{-1}$, on random vectors: True; (2) $E(L)=(-1)^{H+q}Q(L)$: True (16/16 T=0); (3) $\iota_E\mathrm{vol}=(-1)^{q-1}\prod_PR_P^{n_P-1}\,\Omega_R\wedge dx$ evaluated on random vectors: True; (4) $A_{12}=\Omega_R/(2U)$, i.e. $\rho_2=\Omega_R/(2U^{3/2})$: True; (5) $\det[t;\partial t/\partial\xi]=(-1)^{q-1}\prod_Pw_P^{-(n_P+1)}$, $\xi=(w_A,w_B,x^A_{2..},x^B_{2..},x^C_{2..})$: True (32/32). The S3 eq-rho2 numerator with $M_{11}=R_A+R_B$, $M_{12}=-R_B$, $M_{22}=R_B+R_C$ reduces by hand to $-\Omega_R$, giving the same $\rho_2$.

C4: $-\sum_ey_e\cdot z_e=\sum_{v\ne\rho}\lambda_v\cdot(X_v-X_\rho)$ with the claimed $X_v$, random rational $z$, roots $N'$ and the first interior vertex of each path: `C4 X_v: True` (31/31 graphs with $H\ge1$).

C7: $\mathrm{sgn}\det[B^{\mathsf T};C_3]=(-1)^{q+1}$: True (32/32), and on all 343 graphs theta(p,q,r), $1\le p,q,r\le7$: True. Density with htk's `cycle_matrix()` $=(-1)^{(q+1)D}D[C_3]$ pointwise: True (32/32).

C8: relabelling $\rho\leftrightarrow N'$, $s=\mathrm{sgn}\det[B^{\mathsf T};C_3]$, $s'$ for the relabelled grounding: $ss'=-1$ for every interior $\rho$ of all 32 graphs: True. Hence $K_{\rm htk}[G']=(ss')^DK_{\rm htk}[G]=(-1)^DK_{\rm htk}[G]$.

### C3 chart density (v1_pointwise.log)

$$g=D(t(\xi))\,|\det[t;\partial_\xi t]|=\begin{cases}(-1)^{H(H+1)/2+1}\prod_Pw_P^{H-n_P}\,Q(L(x))&T=0\\(-1)^{H(H+1)/2}\tfrac{D_\theta}2\prod_Pw_P^{H+1/2-n_P}&T=1\ (\text{sign and square})\end{cases},\qquad U=\frac1{w_Aw_Bw_C}:$$

`C3 chart density: True` (32/32). z ≠ 0: see C12 block (direct Schwinger integral = representation).

### C5, C6 values (v2_values.log)

Own integration: Cheng–Wu gauge $R_C=1$, $d^{n}t_P=R_P^{n-1}dR_P\,d^{n-1}x$, $\alpha_P=1/R_P$,

$$\int_{(0,\infty)^2}\alpha^{a-1}\beta^{b-1}(1+\alpha+\beta)^{-(a+b+c)}=\frac{\Gamma(a)\Gamma(b)\Gamma(c)}{\Gamma(a+b+c)},$$

$T=1$: $(a,b,c)=(H+\tfrac32-p,H+\tfrac32-q,H+\tfrac32-r)$; $T=0$: $(H+1-p,H+1-q,H+1-r)$. This R-integral computed independently with mpmath (dps 30) for $D\le4$: 11.4–25.4 agreeing digits. $E_x[Q(L)]$ from exact Dirichlet moments $E[x_j]=1/n$, $E[x_jx_l]=(1+\delta_{jl})/(n(n+1))$ including all covariances (moments also checked against sympy iterated integrals, n=2,3,4: True).

Printed: `c1 = claim: True` (16/16), `c0 = claim: True` (16/16), `sign_htk = claim: True` (32/32), `E[Q(L)] = Q(E[L]): True` (16/16), `Lbar = sum mu_v lam_v (claimed mu): True`, `K_htk = (-1)^((q+1)D) K[C3] (htk.kernel_exact): True` (32/32, 2 random source sets each, $\lambda_N\ne0$), `C5/C6 examples (htk convention) equal own closed values: True` (all seven examples, 3 source sets).

Lemma (v8_polarization.log): $B(u,v)=\det[(u;0),(0;v),\dots]$ vanishes on $\mathrm{span}\Lambda_A$, $\mathrm{span}\Lambda_C$; on $\mathrm{span}\Lambda_B$ it is antisymmetric and nonzero in general (e.g. H=4 (3,3,3): max $|B|=49726565/1728$). Only its symmetric part (the polarization) vanishes, which is what $E[Q(L)]=Q(E[L])$ needs.

### C8 interior root (v2_values.log)

38 cases (T=0: H=2,3,3,4; T=1: H=1,2,2,3,3; every interior $\rho$): `K_htk[relabelled] = (-1)^D K_htk, Lbar = sum (mu_v - mu_rho) lam_v: True` (htk.kernel_exact on the relabelled graph against own closed value with full sources).

### C9 (v6_expansion.log)

`eps1 / eps0 Laplace formulas = direct 2H x 2H determinants: True (168 graph/T cases, D <= 7, all path orders)`; `eps1 reproduces S3 signs (-1)^d (D_355): True, (-1)^{d+1} (D_445): True`; own term counts T=0 `[1, 4, 9, 9, 16, 16, 30, 25, 25, 25, 69, 36, 36, 36, 132, 102]` = claim; T=1 = claim; `results.json tables (32): terms = own expansion, K_htk_terms = own K_htk at random sources, c, sign_htk, Lbar weights: True`.

### C10 S1 (v3_sources.log)

S1 eq-1-40 with eq-1-38 at $n=2$, integrated by sympy over $(0,1)^2\times\Delta^2$, $\langle\cdot,\cdot\rangle$ any antisymmetric pairing:

$$I^{\rm S1}(\lambda;0)=\int d^4r\;r_2\,\langle\lambda_1,\lambda_2-r_3(r_1\lambda_1+\lambda_2)+r_4\lambda_3\rangle\,\langle r_1\lambda_1+\lambda_2,\lambda_3\rangle=\frac{(2\langle12\rangle+\langle13\rangle)(\langle13\rangle+2\langle23\rangle)}{24}$$

`equals claim -(1/24)[l1, l3+2 l2][l3, l1+2 l2]: True`, `at l2 = 0 equals eq-1-42 det(l1,l3)^2/24: True`, `equals minus RHS of eq-spider-two-arm-value with (v11,v21,A)=(1,3,2): True`. Pointwise: eq-operatope integrand $\bigwedge_{e\notin T}(-dy^1_e\wedge dy^2_e)$ (own implicit derivatives of the Kirchhoff solution) tree-independent over 8 trees, $D_{\rm S1}/D[C_3]=-1$ (fan edge order); $I^{\rm S1}$ = eq-1-40 value and $=-K_{\rm htk}$: True.

### C11 S2 (v3_sources.log)

From S2 eq-chain, eq-bubble and the box (eq-A-currents, eq-B-currents) with $\lambda_3\ne0$, $\lambda_q\ne0$:

$$\mathcal J_A=r_1^2(r_4\langle123\rangle+r_5\langle124\rangle+r_6\langle125\rangle)(\langle345\rangle+(1-r_2)\langle145\rangle+(1-r_3)\langle245\rangle),$$
$$\mathcal J_B=r_1(1-r_1)(r_4\langle uwq\rangle+r_5\langle uwa\rangle+r_6\langle uwb\rangle)(\langle qab\rangle+(1-r_2)\langle uab\rangle+(1-r_3)\langle wab\rangle),$$

integrated over $D_A$, $D_B$: `(3*b123 + 2*b124 + b125)*(b145 + 2*b245 + 3*b345)/432`, `-(2*bqab + buab + bwab)*(2*buwa + buwb + 3*buwq)/288`, both equal to the claim and to eq-IA, eq-IB at $\lambda_3=0$, $\lambda_q=0$. Pointwise from eq-kernel-def ($\varepsilon_T\det B_T\bigwedge_{e\notin T}d^3y_e$, $(B_T)_{ve}=+1$ at $e(0)$): tree-independent over 15 (A) and 16 (B) trees; ratio to $D[C_3]$ equal to 1; $I^{\rm S2}$ = general-source formula and $=-K_{\rm htk}$: True for both.

### C8/S4 ladder(2) (v3_sources.log)

check_ladder edges, cotree [0,4], root $b_2$ (interior vertex of C): $\bigwedge_{c}d^3y_c/D[C_3]=1$ at 3 points; corner sources: `I^S4 / Delta^2 at 3 corner-source sets: ['-1/432', '-1/432', '-1/432']`.

### C5/S3 (v3_sources.log, S3_pointwise_run.log)

S3 explicit $D_{166},\dots,D_{445}=D_\theta$: True. R-c values = Gamma formula, and with $\varsigma_\theta=-1$ (S3 proof of thm-theta): $\mathcal K^{\rm S3}=\varsigma_\theta K[C_3]=-K[C_3]=+c_{pqr}\pi D_\theta$: True. S3 pointwise.py (run read-only) prints ratios $D|{\rm jac}|/(\tfrac12D_\theta\prod w^{H+1/2-n})$: $+1$ (H=0), $-1$ (H=1), $-1$ (H=5, five graphs) $=(-1)^{H(H+1)/2}$.

### C12 numerics (v4_numerics_z0.log, v4_numerics_z.log, v7_oracle_rerun.log)

Own direct integral (no closed form, no explorer chart): pullback of $\Omega$ to the section $R_C=1$, chart $w_C=v^2$, $c=\sin^2(\pi\theta/2)$, $\alpha_A=\rho c$, $\alpha_B=\rho(1-c)$, $\rho=(1-v^2)/v^2$, stick-breaking fibres, $K=\int d\xi\,\Omega(\partial_\xi t)\,\mathrm{sgn}\det[t;\partial_\xi t]$, Gauss–Legendre. Explorer's sources (data); own exact $K[C_3]$ equals the explorer's exact value in all six cases.

| (H,T) pqr | exact K[C3] | own numeric (m=12,12,3) | (20,20,3) | (30,30,4) | digits |
|---|---|---|---|---|---|
| (1,1) (1,2,2) | 2π/9 | 0.698131700797724 | 0.698131700797743 | 0.698131700797667 | 14.0, 13.8, 13.0 |
| (2,1) (1,3,3) | −13π/60 | −0.680678408277784 | −0.680678408277804 | −0.680678408277738 | 14.2, 13.7, 13.1 |
| (2,1) (2,2,3) | −13π/90 | −0.453785605518526 | −0.453785605518526 | −0.453785605518525 | 14.9, 15.9, 14.7 |
| (2,0) (1,2,2) | 3/2 | 1.500000000000002 | 1.499999999999997 | 1.499999999999991 | 15.0, 14.7, 14.2 |
| (3,0) (1,3,3) | −25/192 | −0.130208333333333 | −0.130208333333333 | −0.130208333333333 | 15.0, 14.7, 14.2 |
| (3,0) (2,2,3) | −25/72 | −0.347222222222222 | −0.347222222222222 | −0.347222222222221 | 15.0, 15.1, 14.3 |

htk oracle re-run (default parameters): (1,1) (1,2,2) `0.6981311880 +- 2.30e-06`, |diff| 5.13e-07; (2,0) (1,2,2) `1.4999993109 +- 4.50e-06`, |diff| 6.89e-07 — identical to the explorer's log. The other four oracle values were not re-run (18–131 s each on a loaded machine); from the logged numbers, |oracle − exact K_htk| = 1.11e-4, 7.3e-5, 9.5e-6, 2.9e-5 against errors 7.1e-4, 4.7e-4, 1.8e-4, 3.7e-4, signs consistent with $K_{\rm htk}=(-1)^{(q+1)D}K[C_3]$.

z ≠ 0 (closes the explorer's open issue on numerics at z ≠ 0): random rational sources, $z_e\in\{-4..4\}/40$; direct Schwinger integral with $e^{-\sum y\cdot z}$ (own currents) and the claimed representation C3 ($e^{\sum\lambda_v\cdot X_v}$, own code from the claimed $X_v$) on the same grids:

| (H,T) pqr | m | direct | representation | digits |
|---|---|---|---|---|
| (1,1) (1,2,2) | (36,36,7) | 8.286679243087496 | 8.286679243088082 | 13.2 |
| (2,0) (1,2,2) | (36,36,7) | 5.755189168301206 | 5.755189168301199 | 15.0 |
| (2,1) (2,2,3) | (36,36,7) | −0.661325791539754 | −0.661325791539754 | 15.0 |
| (3,0) (2,2,3) | (36,36,7) | 10.840197556885590 | 10.840197556885572 | 14.8 |
| (3,0) (1,3,3) | (36,36,7) | 8.531232077676814 | 8.531232077676783 | 14.4 |

Grid convergence of the direct value: (2,0): stable to 15 digits from m=(30,30,5); (2,1) (2,2,3): 15 digits; (1,1): relative change ≤ 1.7e-13 among the last three grids; (3,0) (2,2,3) and (1,3,3): relative change 2.3e-11 and 4.7e-12 from fibre order 6 to 7.

### C13 (v5_formscheck.log)

`S3 mc_theta(1,2,2,1,200000) [seed 1]: mean -0.2539835710  err 0.0041449808  closed 0.2843069064`; seed-1 sources `[0.34558419 0.82161814 0.33043708]`;

$$K[C_3]/\pi=-\tfrac13\lambda_{B1}\lambda_{C1}=-\frac{22026136797979191689044246579867}{243388915243820045087367015432192}\ (\text{= claimed fraction: True}),\quad K=-0.2843069064265951,$$

`mc_theta 'closed' / K[C3] = -0.9999999999999994`. Own replay of the same stream with the own brute force: `mean -0.2539835710  err 0.0041449808`, relative differences 3.66e-15, 3.26e-14; $f$ min −10644.1034, max −2.1161; $(K-\text{mean})/\text{err}=-7.32$. mc_theta single-sample summand = own $D(t)/4!$ at 5 seeds (rel < 1e-12): True.

### C14 (v5_formscheck.log), own stream, $f=-D_\theta/(2U^{5/2})$ (pointwise-verified)

At the corners $\{t_{A1},t_{C1},t_{C2}\}\to0$ and $\{t_{A1},t_{B1},t_{B2}\}\to0$: $U=(1-\rho)\rho+R_AR_{C}$ resp. $(1-\rho)\rho+R_AR_B$, so $U/\rho\to1$; measured max $|U/\rho-1|$ = 8.3e-3, 8.8e-3 on $\rho<0.01$. $P(\rho<r)=4r^3-3r^4$ (Beta(3,2)); hence $P(|f|>x)\sim8(c/x)^{6/5}$, $c=|D_\theta|/2$, and $\int f^2\propto\int\rho^{-5}\rho^2d\rho=\infty$.

`x=1000: emp 1.900e-04 pred 1.830e-04  x=10000: emp 1.150e-05 pred 1.154e-05  x=100000: emp 7.500e-07 pred 7.284e-07`, fitted slope −1.202. Sample $\langle f^2\rangle$ at n = 1e4, 1e5, 1e6, 4e6: 816.9, 1.981e5, 3.855e4, 5.517e4 (no convergence). Runs:

| n | runs | median(mean/exact) | fraction below | fraction \|z\|>3 |
|---|---|---|---|---|
| 2000 | 2000 | 0.8634 | 0.844 | 0.378 |
| 20000 | 600 | 0.9054 | 0.817 | 0.405 |
| 200000 | 120 | 0.9467 | 0.767 | 0.317 |

slope of $\log(1-\text{median ratio})$ vs $\log n$: −0.204 (generalized CLT with $\alpha=6/5$: $n^{1/\alpha-1}=n^{-1/6}$).

## Logical gaps in the derivation as stated (REPORT.md absent)

1. C1 (2), T=1: dropping the "terms without $a\wedge b$" needs $\Omega_R\wedge a=\Omega_R\wedge b=0$ and $d\hat m_P$ free of $dR$ ($a,b,\Omega_R$ basic in $(R_A:R_B:R_C)$, $\Omega_R=\iota_{E_R}dR_AdR_BdR_C$); not stated. True; the end density is verified pointwise.
2. C1 (2): the signs $(-1)^{H(H-1)/2}$ and $(-1)^{H(H-1)/2+H+q-1}$ are asserted without a written computation. Only $E=(-1)^{H+q}Q$ and the final densities (which, with (3), (4), fix both signs) were checked.
3. C1 (5): the sign $(-1)^{q-1}$ of the chart Jacobian is stated without derivation; checked at rational points only.
4. C6 lemma: "polarization of Q vanishes on span(Λ_P)×span(Λ_P)" holds for the symmetric part only; the bilinear form $B(u,v)$ is antisymmetric and nonzero on $\mathrm{span}\Lambda_B$ (v8_polarization.log). The argument goes through because the covariance matrix is symmetric.
5. C7: the factors $\det B_T=(-1)^{p-1}$, $\varepsilon_T=(-1)^{(N-2)+(r-1)}$, $\sigma_T=(-1)^q$ were not checked individually; the resulting sign was checked on 343 + 53 graphs.
6. C7 source normalizations, C10: $I^{\rm S1}=-K_{\rm htk}$, $I^{\rm S2}=-K_{\rm htk}$ rest on CONVENTIONS in the explorer's evidence; here they were re-derived from S1 eq-1-40/eq-operatope and S2 eq-kernel-def/moves, and $\mathcal K^{\rm S3}=-K[C_3]$ from S3's $\varsigma_\theta=-1$.
7. C11: S2 states (A), (B) only at $\lambda_3=0$, $\lambda_q=0$; the general-source formulas are new and their derivation from S2's moves is not written in the claim (done here).
8. C14: "typically below the true magnitude" and the $O(n^{-1/6})$ deviation are stated from the stable-law limit without proof; supported numerically above.

## Not done

- Four of the six htk oracle runs of C12 were not re-run.
- The explorer's own deterministic digit counts (12.6–15.4) and their Monte Carlo statistics were not reproduced with their streams; own streams give the numbers above.
