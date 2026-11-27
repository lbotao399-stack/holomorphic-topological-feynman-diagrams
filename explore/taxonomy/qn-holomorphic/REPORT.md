# $Q_n=c^n\sum_\Gamma\kappa_\Gamma\hat Q_\Gamma$ for the holomorphic ($T=0$) twists, any $H$, any quadratic Lie superalgebra

Tags: [chk] exact (script, `.log` here); [drv]; [src]; [open]. Arithmetic mod $P=2^{31}-1$ with rational reconstruction; colour sums float64 with exactness check $|x|<2^{52}$; kernels from htk spred (Fractions) or mod $P$.

## open (newest first)
- 6d (1,0), $\mathfrak{sl}_2$ + 16 fund, order 2 on $O=\xi^2\varphi$: tree composition not finished in 24 min (U block 64); unicyclic compositions vanish on all tests, necessity untested.
- 8d order 3: 12 graphs, 9 slots ($9!$ rows), 3 block-triangular densities; $\sigma_4(3)$ undetermined.
- Closed formula for $\sigma_H(n)$: $H=2,3$ all $+$; $H=4$: $(+,-)$.
- 6d (1,0) $Q_3$: 14/18 unicyclic $n=3$ graphs non-SP, no exact evaluator without python-flint.
- 10d: no consistent $Q_1$ (hexagon anomaly); $Q_2$ only as kernels.
- Scans cut for time: 4d N=2 $N_f=3,4$ at $k\le4$, $|d|\le2$; 2d (4,4) $k=5$; 2d (8,8).

## 1. Engine

Theory $=(H,\mathfrak g_{\rm tot})$, $\mathfrak g_{\rm tot}$ quadratic Lie superalgebra (`superalg.SuperLie`): basis $T_A$, parity $p_A$, $f[A,B,C]$, pairing $w[A,B]$ of parity $p_\omega\equiv H+1$. Axioms [chk `test_alg.log`, 17 theory/algebra pairs]: graded antisymmetry, super Jacobi, $\omega(X,Y)=(-1)^{|X||Y|}\omega(Y,X)$, $\omega([X,Y],Z)=\omega(X,[Y,Z])$, $\hat Q_0^2=0$ (the other index convention for $G$ fails in 7/17).

Letters $\varphi^A$, $\pi_A=p_A+1$, jets $\partial^\alpha\varphi^A$, $\alpha\in\mathbb N^H$;
$$I=\tfrac16\omega(\Phi,[\Phi,\Phi])=\tfrac16V_{ABC}\varphi^A\varphi^B\varphi^C,\quad G=(w^{-1})^T,\quad \hat Q_0\varphi^D=\tfrac12(-1)^{\pi_DH}G^{AD}V_{ABC}\varphi^B\varphi^C$$

| $d$ | $\mathcal N$ | ESW | $H$ | $\mathfrak g_{\rm tot}$ |
|---|---|---|---|---|
| 10 | (1,0) | 4.4 | 5 | $\mathfrak g$ |
| 8 | 1 | 6.1 | 4 | $T^*_1\mathfrak g=\mathfrak g\otimes\Lambda[\theta]$ |
| 6 | (1,1) | 8.3 | 3 | $\mathfrak g\otimes\Lambda[\theta_1,\theta_2]$ |
| 6 | (1,0) | 8.1, 8.2 | 3 | $L_U=\mathfrak g\oplus\Pi U\oplus\mathfrak g^\vee$, $[u,u']=-\sum_b(\Omega\rho_b)_{uu'}X^b$ |
| 4 | 1 | 10.1 | 2 | $T^*_1(\mathfrak g\ltimes\Pi R)$ |
| 4 | 2 | 10.2 | 2 | $T^*_1L_U$ |
| 4 | 4 | 10.9 | 2 | $\mathfrak g\otimes\Lambda[\theta_1,\theta_2,\theta_3]$ |
| 2 | (2,2) | 12.1 | 1 | $T^*_0((\mathfrak g\ltimes\Pi R)\otimes\Lambda[\varepsilon])$ |
| 2 | (0,2) | 12.16 | 1 | $T^*_0(\mathfrak g\ltimes\Pi R)$ |
| 2 | (4,4) | 12.7 | 1 | $T^*_0(L_U\otimes\Lambda[\varepsilon])$ |
| 2 | (8,8) | §12 | 1 | $\mathfrak g\otimes\Lambda[\theta_1..\theta_4]$ |

Graphs (`rgraphs.py`) [drv]: $|I|=Hn+1$, $|V|=(H-1)n+2$, $(H-1)|I_S|\le H|V_S|-H-1$ $\Rightarrow$ $\Gamma\smallsetminus O$ connected, $s=n+1-b$ spokes; trees ($b=0$): $(H-1)n+1$ bulk vertices; unicyclic ($b=1$): $n$ spokes, cycle supertrace. Counts [chk `rgraphs.log`]: $H=3$ trees $1,2,5,15$ ($\sum1/|\mathrm{Aut}|=8$ at $n=4$); $H=4$: $1,3,12$; $H=5$: $1,4$; unicyclic $H=3$: $0,3,18$; $H=4$: $0,7$; $H=5$: $0,13$.

$$\hat Q_\Gamma O=\sum\Big[\prod_vV_{{\rm legs}(v)}\prod_eG^{A_eB_e}/G_{\rm den}^{|I|}\Big]\varepsilon_{\rm Koszul}\,[\lambda^\beta]I_\Gamma(\lambda;\alpha_{\rm in})\prod_v\partial^{\beta_v}\varphi^{B_v}\,O_{\rm rest},\qquad I_\Gamma=(-1)^{n(n+1)/2}K_\Gamma$$
Nilpotency: $\sum_{a+b=n}Q_aQ_b=0\iff\sum\kappa\kappa'\hat Q_{\Gamma_a}\hat Q_{\Gamma_b}=0$, $\kappa_{Q_0}=1$, $\kappa_1=1/|\mathrm{Aut}|$. Drivers: `explicit.GraphOp.apply`, `fmode.evaluate`, `nilp_explicit.py`, `run_nilp_f.py`, `lifts.py`, `run_lifts.py`.

## 2. Validation

| target | source | result | log |
|---|---|---|---|
| fan kernels $n=1..6$ | LLC eq-1-42 | $I=(-1)^{n(n-1)/2}[\lambda_a\lambda_b]^n/((n!)^3(n+1))$ | validate_kernels (a) |
| 4d N=1 $\kappa_n$ | LLC Thm b-00088 | $\kappa_n=\tfrac12$: $n=2,3$ (sl2 pure, sl2+fund, sl3+fund+afund); $n=4$ (4d N=2 sl2+fund) | nilp3_4dN1_*, a_nilp_*, a_nilp4_* |
| 4d N=4 $Q_1,Q_2$ | LLC eq-2-26–2-28 | $a_2a_0/a_1^2=-\tfrac1{24}$, $b_1/a_1=-1$ | n4_llc_check |
| pure N=2 SU(2), lowest $d_2$ | LLC table local-ch9-higher-lifts | sector $\{b,\beta_\Phi^3,\gamma_\Phi\}$, $d=0$: $\mathrm{rk}\,d_2=1$ | lifts_4dN2_sl2_pure_d0 |
| N=2 SU(2), $N_f=4$ | LLC local-ch9-sqcd-results | $d_1=0$, $d_2=0$ | a_lifts_4dN2_sl2_nf4_* |
| 6d kernels $n=4$ | 6d-nilpotency §2 | 15/15 | validate_kernels (b) |
| 6d (1,1) order 2, 3 | 6d-nilpotency §4 | so5 residual 0, rank 2/2, 5/5 | nilp2_6d11_*, nilp3_6d11_* |
| 8d $Q_3$ kernels | atlas e1_cone_kernels | 12/12 | b_validate_cone |
| relabelling $K\to\mathrm{sgn}(\pi_v)^H(-1)^{Hk_{\rm rev}}\mathrm{sgn}(\pi_e)^{H+1}K$ | CONVENTIONS §2.2 | 12/12, $H=2,3,4$ | b_validate_cone (d) |

Wick weights [chk]: $\kappa_\Gamma=\sigma_H(n)/|\mathrm{Aut}\,\Gamma|$; $H=2$: $\sigma=+1$, $n\le4$; $H=3$: $+1$, $n=2,3$; $H=4$: $\sigma(2)=-1$ (so5 rank 3/3); $H=1$: $\hat Q_{n\ge2}=0$ [drv].

## 3. Kernels and $Q_n$

$[v_1..v_H]=\det(\lambda_{v_1},..,\lambda_{v_H})$ [chk kernels_detforms, kernels_closed, b_theta_forms]:
- $H=1$ banana: $I(\lambda;a,b)=-\tfrac{a!b!}{(a+b+1)!}\lambda^{1+a+b}$.
- $H=2$ fan: $I_n=(-1)^{n(n-1)/2}[\lambda_0\lambda_n]^n/((n!)^3(n+1))$.
- $H=3$: box $-[012]/6$; tree #1 $[\lambda_0,\lambda_1,2\lambda_3+\lambda_4][\lambda_0+2\lambda_1,\lambda_3,\lambda_4]/432$; tree #2 $-[\lambda_0,\lambda_1,\lambda_3+\lambda_4][\lambda_0+2\lambda_1,\lambda_3,\lambda_4]/288$; unicyclic #3 $-[034]^2/432$, #4, #5 $=0$ without jets.
- $H=4$: pentagon $[0123]/24$; $\theta(1,4,4)$: $[\lambda_0,\lambda_1,\lambda_2,3\lambda_4+2\lambda_5+\lambda_6][\lambda_0+2\lambda_1+3\lambda_2,\lambda_4,\lambda_5,\lambda_6]/11520$; $\theta(2,3,4)$, $\theta(3,3,3)$: 4968, 6624 monomials.
- $H=5$: hexagon $-[01234]/120$; $n=2$: 4 trees, 13 unicyclic, exact [chk b_10d_q2].

2d: $\hat Q_1(\varphi^A\varphi^B)=2(-1)^{p_B}\partial[\omega(\Phi,[T^A,T^B])]$, second-order derivation [chk c_q1formula, c_q1jets, c_q1second].

4d N=1,2,4: $Q_n=c^n\tfrac12\hat Q_{\rm fan}$. 6d (1,1): $Q_1=c\tfrac12\hat Q_{\rm box}$, $Q_2=c^2\tfrac12(\hat Q_{2.1}+\hat Q_{2.2})$, $Q_3=c^3\sum\hat Q_\Gamma/|\mathrm{Aut}|$ ($|\mathrm{Aut}|=2,2,1,8,2$). 8d: $Q_1=c\tfrac12\hat Q_{\rm pent}$, $Q_2=-c^2\sum_\theta\hat Q/|\mathrm{Aut}|$ ($2,1,6$).

## 4. Per theory

**4d N=1.** Nilpotency through order 3, rank 1, $\kappa=\tfrac12$ (sl2 pure 1438, 5781 equations; sl2+fund 3398, 9577; sl3+fund+afund 16884). Anomalous sl3+fund: $\{Q_0,\hat Q_1\}\ne0$, order-2 system inconsistent (rank 1, augmented 2) [chk a_nilp_4dN1_sl3_fund_anomalous]; obstruction $\mathrm{str}_{\mathfrak g_{\rm tot}}X^3=2(\mathrm{tr}_{\rm adj}-\mathrm{tr}_R)X^3$ [drv].

**4d N=2, hypers $U$.** Letters $c$, $\varphi$, $\xi$ (odd), cotangents $b,\varsigma,a$ ($a$ adjoint scalar). Nilpotency, $\kappa=\tfrac12$, residual 0: sl2 + $N_f=1$ order 4 (6229 equations); $N_f=2,3,4$ order 3; pure order 3; sp4 + half-hyper in $\mathbf4$ order 3 (31067).

$T=\partial_2c^1\partial_1c^2-\partial_1c^1\partial_2c^2+2\partial_2c^3\partial_1c^3\propto\mathrm{tr}\,\varepsilon^{\dot a\dot b}\partial_{\dot a}c\partial_{\dot b}c$ [chk a_qn_examples]:
$$Q_1\mathrm{tr}(\xi a)=4T,\quad Q_1(\varphi_R\varsigma_R)=Q_1(\varphi_{R^*}\varsigma_{R^*})=T,\quad Q_1\mathrm{tr}(aa)=0,\quad Q_1\mathrm{tr}(ba)=4T'$$
Lowest lifts ($k=2$, $d=0$) [chk a_lifts_4dN2_*_k2d0*]:
$$O_J=\mathrm{tr}(\xi a)-\sum_f\varphi^f_R\varsigma_{f,R},\qquad Q_1O_J=(4-N_f)T;\qquad O'_J=-\mathrm{tr}(ba)+\sum_f\varsigma_{f,R}\varsigma_{f,R^*},\ Q_1O'_J\propto(4-N_f)$$
$\mathrm{rk}\,d_1=1$ per class for $N_f=1,2,3$, $0$ for $N_f=4$; sp4 + half $\mathbf4$: coefficient $11=2(2T({\rm adj})-T(U))$; $Q_1O_J\propto b_0=2T({\rm adj})-T(U)$, the one-loop N=2 $\beta$-function coefficient. $d_1\ne0$ also at $k=4$, $d=0$ and $k\le4$, $|d|\le2$; $d_2=0$ on every tested $\ker d_1$ ($N_f=1,2,3$); pure sl2: $d_1\ne0$ in 6 sectors ($k=4,5$), $\mathrm{rk}\,d_2=1$ on $\{b,\xi^3,a\}$ and $\{b^2,\xi^2,a\}$.

**4d N=4.** Nilpotency through order 3, $\kappa=\tfrac12$ [chk a_nilp_4dN4_su2]; su2: no $d_1\ne0$ for $k\le3$, $|d|\le2$ [chk lifts_4dN4_sl2_k3d2] (LLC first lift at level 36).

**6d (1,1).** Order 2 residual 0 (rank 2/2), order 3 residual 0 (rank 5/5), so5.

**6d (1,0), hypers.** $\{Q_0,\hat Q_{\rm box}\}=0\iff2\mathrm{tr}_{\rm adj}X^4=\mathrm{tr}_UX^4$ [drv: box colour $=\mathrm{str}_{\mathfrak g_{\rm tot}}(\mathrm{ad}X)^4$]: passes sl2 + 16 fund (with jets), sl2 + adj; fails sl2 pure, + 1 fund, + 15 fund [chk nilp1_6d10_anomaly]. Order 2, sl2 + 16 fund, $O=\xi^3$: trees residual 0, rank 2/2, $\kappa=\tfrac12,\tfrac12$; unicyclic columns zero on this test.

**8d N=1.** Order 1 (su2, so5); order 2 so5 rank 3/3, residual 0, $\sigma(2)=-1$ [chk b_nilp2_8d_so5*]. Lowest lift $[\mathrm{tr}\gamma\gamma]$ ($k=2$, $d=0$): $Q_1\mathrm{tr}(\gamma\gamma)=-\tfrac{16}3\varepsilon^{abcd}\mathrm{tr}_{\rm fund}(\partial_d\gamma\,\partial_ac\,\partial_bc\,\partial_cc)$ [chk b_8d_examples]; $d_1\ne0$ on 17 more sectors ($k\le4$, $|d|\le1$).

**10d (1,0).** $\{Q_0,\hat Q_{\rm hex}\}\ne0$ on $\varepsilon_{ABC}\partial_1c^A\partial_2c^B\partial_3c^C$ (927 monomials, su2); $\kappa_1=0$ forced (su2, so4, so5) $\Rightarrow$ no nilpotent $Q$ at one loop [chk b_10d_anomaly, b_nilp1_10d_*]; cause $O_1\propto\mathrm{tr}_{\rm adj}X^6$ [drv].

**2d.** $Q=Q_0+\hbar Q_1$, $\hat Q_1^2=0$ [chk c_nilp2d]. (0,2): $\{Q_0,\hat Q_1\}=0\iff T(R)=T({\rm adj})$ (passes sl2 + 4 fund, sl2 + adj, sl3 + 6 fund, sl3 + adj, sl3 + 3 fund + 3 afund; fails sl2 + 1, 2, 3, 5 fund, sl3 + 5, 7 fund). (2,2), (4,4), (8,8) pass (supertrace 0). Lifts: (2,2) sl2 + fund lowest $d_1\ne0$ at $k=5$, $d=0$ (4 sectors), $z_0=[-(\tilde r_0)^2\tilde\varepsilon_1+\tilde r_0\tilde r_1\tilde\varepsilon_2+(\tilde r_1)^2\tilde\varepsilon_0](\tilde\psi_0\tilde\psi_1)\mapsto-6[\tilde r_0\tilde\psi_0\tilde\psi_1\partial\tilde\psi_1+\tilde r_1\tilde\psi_0\partial\tilde\psi_0\tilde\psi_1]$ [chk c_example_2d22_q114_d0_k5]; (0,2) sl2 + adj ($k\le6$), sl2 + 4 fund ($k\le5$): $d_1=0$; (4,4) sl2 + fund: $d_1=0$ ($k\le4$); (8,8) sl2: $d_1=0$ in 651 sectors (partial).

## 5. Files
Engine: `superalg.py`, `theories.py`, `wick.py`, `polyops.py`, `rgraphs.py`, `kern_h.py`, `blocktri_h.py`, `rooted_trees.py`, `wick_weights.py`, `explicit.py`, `fmode.py`, `blocksutil.py`, `lifts.py`. Drivers: `test_alg.py`, `validate_kernels.py`, `kernels_closed.py`, `nilp_explicit.py`, `run_nilp_f.py`, `run_lifts.py`, `n4_llc_check.py`, `a_*.py` (4d), `b_*.py` (8d, 10d), `c_*.py` (2d). Note: `validate_kernels.log` part (c) ends in a TypeError; the 8d cone check is redone in `b_validate_cone.log`.

## 6. Replay (main session, scratch copy, `PYTHONPATH=code`)
- `n4_llc_check.py`: $a_2a_0/a_1^2=-\tfrac1{24}$, $b_1/a_1=-1$ (= LLC).
- `nilp_explicit.py 6d10 sl2 - 1 3 12 1`: $\{Q_0,\hat Q_{\rm box}\}\ne0$ (402 rows); `6d10 sl2 fund*16`: $=0$; `10d su2 -`: $\ne0$ (915 rows).
- `a_run_lifts.py 4dN2`: sl2 $N_f=1$, $N_f=3$: $d_1\ne0$ on the two $k=2$, $d=0$ classes; $N_f=4$: none; sp4 half $\mathbf4$: $Q_1$ image coefficient $\pm11$, $\pm\tfrac{11}2$.
- `c_nilp2d.py 2d02 sl2 fund 30 7 2`: $\{Q_0,\hat Q_1\}=0$ on 10/30; `fund*4`: 30/30; `2d22 sl2 fund`: 30/30; $\hat Q_1^2=0$ on 30/30 in all three.
- `test_alg.py`: same output as the agent's log.
- `b_validate_cone.py`: rerunning.
