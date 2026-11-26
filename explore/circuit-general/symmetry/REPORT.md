# Symmetry engine for weighted core periods

Tags: [chk] exact check (script, log in this folder), [drv] derived here, [src] source (location), [open]. Paths relative to `explore/circuit-general/symmetry/`; `cores/` = `explore/circuit-general/cores/`; skeleton names and edge orders = `fiveloop/skel5.json` (V8_3 = cube, V8_4 = Möbius–Wagner V8). K5 edge order (01,02,03,04,12,13,14,23,24,34).

| file | content | log |
|---|---|---|
| `engine.py` | completion/decompletion, canonical form of weighted completed graphs, moves M1 M2 FT TW PR, two-point terminal, search | — |
| `check_completion.py` | completion invariance; cut criterion vs `gutil.convergent` | `check_completion.log` |
| `check_wheel_hub.py` | wheel hub ↔ ∞ | `check_wheel_hub.log` |
| `check_moves.py` | each move, both sides tabulated | `check_moves.log` |
| `check_product_2sum.py`, `check_2sum_formula.py` | product identity vs msys; decompleted 2-sum formula | same-name `.log` |
| `check_two_point.py` | numerator / two-point formula | `check_two_point.log` |
| `check_kite_group.py` | orbit of a generic K4 vector | `check_kite_group.log` |
| `check_fiveloop.py` | dimension shift on five-loop cores; q=2 limits; Z5 vs msys | `check_fiveloop.log` |
| `check_routes.py`, `check_k5_examples.py` | independent chains for the same input | `check_routes.log`, `check_k5_examples.log` |
| `check_z5_closure.py`, `check_z5_joint.py` | msys + symmetry rows, Z5 | `check_z5_closure_q4.log`, `check_z5_joint.log` |
| `coverage.py`, `run_coverage.sh`, `check_large_q.py` | coverage | `coverage_q{7,8}_k{3,0}.log/.json`, `check_large_q.log` |
| `d0_planar_duality.py` | planar duality (orchestrator) | `d0_planar_duality.log` |

## 1. Position space and completion

### 1.1 Position-space period
$$|x|^{-2\beta}=\tfrac1{\Gamma(\beta)}\!\int_0^\infty\!dt\,t^{\beta-1}e^{-tx^2},\qquad \int\!\prod_{v\ne v_0}\!d^dx_v\,e^{-\sum_et_ex_e^2}=\pi^{d(V-1)/2}\Psi_G^{-d/2},\qquad \Psi_G=\sum_T\prod_{e\in T}t_e=\Big(\prod_et_e\Big)U_G(1/t)$$
$$d=2q,\ \ \sum_e\beta_e=q(V-1):\qquad \mathrm{Per}_x(G;\beta):=\int\prod_{v\ne0,1}\frac{d^dx_v}{\pi^q}\prod_e|x_e|^{-2\beta_e}\Big|_{x_0=0,\ |x_1|=1}$$
$J(x_1)=\pi^{q(V-2)}\mathrm{Per}_x\,|x_1|^{-d}$. In $\int d^dx_1J$: radial side $\frac{2\pi^q}{\Gamma(q)}\int\frac{dr}{r}$; Schwinger side $\frac{\pi^{q(V-1)}}{\prod\Gamma(\beta)}\int\frac{d\lambda}{\lambda}P_{\rm pos}$; $\frac{dr}{r}=\frac12\frac{d\lambda}{\lambda}$. Hence
$$\boxed{\mathrm{Per}_x(G;\beta)=\frac{\Gamma(q)}{\prod_e\Gamma(\beta_e)}\int_{\mathbb{RP}_{>0}}\prod_et_e^{\beta_e-1}\Psi_G^{-q}\,\Omega\overset{t\to1/t}{=}\frac{\Gamma(q)\,I_q(G;q-\beta)}{\prod_e\Gamma(\beta_e)}=\Gamma(q)\,p_\beta(G)},\qquad p_\beta:=I_q\Big/\prod_e\Gamma(\beta_e),\ \ \beta=q-\alpha$$
[drv; chk banana: $\mathrm{Per}_x=1\iff I_q=\prod\Gamma(\beta_i)/\Gamma(q)$, cores §0]

| | $U_G$, momentum, $\alpha$ | $\Psi_G$, position, $\beta=q-\alpha$ |
|---|---|---|
| normalization | $p_\alpha=I/\prod\Gamma(\alpha_e)$ | $p_\beta=I/\prod\Gamma(\beta_e)$ |
| series | $p_\alpha$ invariant, $\alpha_1+\alpha_2$ | $p_\beta=G_{2q}(\beta_1,\beta_2)\,p_\beta'$, $\beta'=\beta_1+\beta_2-q$ |
| parallel | $p_\alpha=G_{2q}(\alpha_1,\alpha_2)\,p_\alpha'$, $\alpha'=\alpha_1+\alpha_2-q$ | $p_\beta$ invariant, $\beta_1+\beta_2$ |
| weight 0 | $\alpha_e=0$: $p_\alpha(G)=p_\alpha(G/e)$ | $\beta_e=0$ ($\alpha_e=q$): $p_\beta(G)=p_\beta(G\setminus e)$ |
| unique | triangle $\sum\alpha=2q$ | vertex $\sum\beta=2q$ |

[drv; $G_{2q}$: cores §0]

### 1.2 Completion
$$\hat G=G\cup\{\infty\}:\quad \beta_{v\infty}=2q-\sum_{e\ni v}\beta_e\iff\hat\alpha_{v\infty}=(\deg v-1)q-\sum_{e\ni v}\alpha_e,\qquad \sum_v\beta_{v\infty}=2q,\qquad \sum_{e\ni v}\beta_e=d\ \ \forall v\in\hat G$$
Inversion: $|x_{ij}|^2\to|x_{ij}|^2/(x_i^2x_j^2)$, $d^dx\to d^dx/(x^2)^d$. A balanced integrand is invariant, edges at ∞ are constants, and the conformal group is transitive on triples. Therefore
$$\Pi(\hat G):=\mathrm{Per}_x(\hat G-c)\ \text{independent of }c$$
$$\boxed{\frac{I_q(G;\alpha)}{\prod_{e\in G}\Gamma(q-\alpha_e)}=\frac{I_q(\hat G-c;\hat\alpha)}{\prod_{e\in\hat G-c}\Gamma(q-\hat\alpha_e)}\quad\forall c,\qquad \hat\beta|_{\hat G-c}\in[1,q-1]\ (\hat\beta=0\text{ deleted, parallel }\beta\text{ added}),\ \hat G\text{ convergent}}$$
$$n(\hat G-c)=n(G)+\#\{v:\beta_{v\infty}\ne0\}-\deg_{\hat G}c,\qquad |V(\hat G)|-3=\#\text{integrated points (invariant)}$$
[drv]

Wheels: decompleting at the hub makes ∞ the new hub.
$$\frac{I_q(W_n;s,r)}{\prod_i\Gamma(q-s_i)}=\frac{I_q(W_n;s',r)}{\prod_i\Gamma(q-s'_i)},\qquad s'_i=2q-s_i-r_{i-1}-r_i\in[1,q-1]$$
[drv; chk `check_wheel_hub.log`: W3 q=3..7, W4 q=3..6, W5 q=3..5, W6 q=3,4: 52502 pairs, 0 mismatches]

### 1.3 Convergence
$$\hat G\text{ convergent}\iff \mathrm{cut}_\beta(W)>2q\quad\forall W\subset V(\hat G),\ 2\le|W|\le|V(\hat G)|-2$$
[drv: $W$ collapses iff $\mathrm{in}(W)<q(|W|-1)$; balance gives $\mathrm{in}(W)=q|W|-\frac12\mathrm{cut}(W)$; the condition is symmetric under $W\to W^c$; $\beta\equiv1,d=4$: internally 6-edge-connected]
Agrees with `gutil.convergent` on 23288 decompletions, 0 disagreements [chk `check_completion.log`]. Consequence: a 2-vertex cut of $\hat G$ means divergence [drv].

### 1.4 Completion checks
| core | q | orbits | pairs of tabulated decompletions | mismatches |
|---|---|---|---|---|
| K4→K5 | 3..7 | 2, 9, 25, 66, 145 | 8, 33, 94, 243, 529 | 0 |
| W4 | 3..6 | 9, 89, 541, 2377 | 28, 239, 1248, 4867 | 0 |
| W5 | 3, 4 | 16, 342 | 20, 399 | 0 |
| prism | 3, 4, 5 | 6, 106, 951 | 13, 172, 1142 | 0 |
| K33 | 3, 4, 5 | 3, 35, 279 | 4, 40, 196 | 0 |

Pair types include (W3,W3), (W4,W4), (W4,prism), (W4,K33), (W3,W5), (W_n,SP) [chk `check_completion.log`; rerun by the orchestrator 2026-11-23, identical].

## 2. Identities (completed, balanced $\hat G$; $\Pi$-normalization)

**(U) uniqueness.** $\deg_{\hat G}v=3$, neighbours $a,b,c$, $\beta_a+\beta_b+\beta_c=2q$ (automatic), $\beta_i\in[1,q-1]$:
$$\Pi(\hat G)=\prod_{i}\frac{\Gamma(q-\beta_i)}{\Gamma(\beta_i)}\;\Pi(\hat G-v+\triangle),\qquad \beta'_{bc}=q-\beta_a,\ \beta'_{ac}=q-\beta_b,\ \beta'_{ab}=q-\beta_c$$
[drv: $\int\frac{d^dx}{\pi^q}\prod_i|x-x_i|^{-2\beta_i}=\prod_i\frac{\Gamma(q-\beta_i)}{\Gamma(\beta_i)}\prod|x_{jk}|^{-2(q-\beta_i)}$]
Decompleted: $v\not\sim c$ gives the unique vertex $\sum_{e\ni v}\alpha_e=q$ and $I_q(H)=I_q(H')$ [src cores `startri.py`, chk 805 there]. $v\sim c$ gives the series rule, factor $G_{2q}(\beta_1,\beta_2)$.
[chk `check_moves.log` M1: 1347 pairs, 0 mismatches]

**(U\*) dual uniqueness.** Triangle with parts $\gamma_e\le\beta_e$, $\gamma\ge1$, $\sum_\triangle\gamma=q$; star vertex $s$ with $\beta_{sa}=q-\gamma_{bc}$, …:
$$\Pi(\hat G)=\prod_\triangle\frac{\Gamma(q-\gamma_e)}{\Gamma(\gamma_e)}\;\Pi(\hat G-\gamma_\triangle+\mathrm{star}_s)$$
Decompleted ($\gamma=\beta$, $c\notin\triangle$): momentum-space unique triangle $\sum_\triangle\alpha=2q$, $I_q(H)=I_q(H')$.
[drv; chk M2: 51471 pairs incl. edge splitting, 0 mismatches]

**(FT) Fourier / planar duality.** $H=\hat G-c$ planar, 2-connected, $\beta|_H\in[1,q-1]$; $H^*$ with $\alpha^*_{\pi(e)}=\beta_e$:
$$\Pi(\hat G)=\prod_{e\in H}\frac{\Gamma(q-\beta_e)}{\Gamma(\beta_e)}\;\Pi(\widehat{H^*}),\qquad I_q(H^*;\alpha^*)=I_q(H;q-\alpha^*\circ\pi),\qquad p_\alpha(H^*)=p_\beta(H)$$
[drv: $\Psi_{H^*}=U_H$; per edge $\int\frac{d^dx}{\pi^q}e^{ipx}|x|^{-2\beta}=4^{q-\beta}\frac{\Gamma(q-\beta)}{\Gamma(\beta)}|p|^{-2(q-\beta)}$; chk FT: 8035 pairs, 7036 non-isomorphic, 0 mismatches; `d0_planar_duality.log` 56 pairs]

**(PR) product at a 3-vertex cut.** $S=\{a,b,c\}$, $\hat G-S=C_1\sqcup C_2$, $w^{(i)}_s=\sum_{s\text{–}C_i}\beta_e$:
$$\Pi(\hat G)=\Pi(\hat G_1)\,\Pi(\hat G_2),\qquad \hat G_i=\hat G[C_i\cup S]-E(S)+\triangle^{(i)},\quad t^{(i)}_{ab}=q-\tfrac12\big(w^{(i)}_a+w^{(i)}_b-w^{(i)}_c\big)\in\mathbb Z$$
[drv: $a,b,c\mapsto0,1,\infty$ factorizes the integrand; $w_a+w_b+w_c=2q|C_i|-2\,\mathrm{in}(C_i)$ is even; $|C_1|=1$ reduces to (U)]
Decompleted at $c\in S$: 2-vertex cut $\{a,b\}$, $H=H_1\cup H_2\,(\cup\,ab)$:
$$\frac{I_q(H;\alpha)}{\Gamma(q-\alpha_{ab})}=\frac{\Gamma(q)\,I_q(H_1+e_{\omega_1})\,I_q(H_2+e_{\omega_2})}{\Gamma(q-\omega_1)\Gamma(q-\omega_2)},\qquad \omega_i=q(n_i+1)-\sum_{H_i}\alpha,\quad \omega_1+\omega_2=\alpha_{ab}$$
With no edge $ab$: $\omega_1+\omega_2=q$ and the left factor is 1.
[chk: K4+K4 vs msys (independent IBP system) q=3,4: 6+121 orbits; explicit formula = engine q=3..6: 11, 274, 2729, 18623 orbits; `check_moves` PR (both sides ≥2 vertices): 67; all 0 mismatches]

**(TW) weighted twist.** $S=\{a,b,c,d\}$ a 4-vertex cut, $C$ a union of components of $\hat G-S$, $w_s=\sum_{s\text{–}C}\beta_e$:
$$w_a=w_b,\ w_c=w_d\ \Longrightarrow\ \Pi(\hat G)=\Pi(\hat G^\sigma),\qquad \sigma=(ab)(cd)\text{ applied to the }C\text{–}S\text{ edges}$$
[drv: $f_C(x_a,\dots,x_d)$ is conformally covariant with weights $w$; a σ-invariant prefactor $\prod|x_{ij}|^{-2c_{ij}}$ with $c_{ac}=c_{bd}$, $c_{ad}=c_{bc}$ exists iff $w_a=w_b,w_c=w_d$; cross-ratios $u,v$ are σ-invariant. $\beta\equiv1$: Schnetz twist, src arXiv:0801.2856, not reread]
[chk TW: 5992 pairs, 509 non-isomorphic, 0 mismatches]

**(DS) dimension shift.**
$$I_{q-1}(G;\alpha)=\sum_{C\ \rm cotree}I_q(G;\alpha+1_C)\qquad[U^{1-q}=U\cdot U^{-q}],\qquad P^{(q-1)}_{\rm pos}(\beta)=\sum_TP^{(q)}_{\rm pos}(\beta+1_T)$$
Valid when every term converges. Not used by the engine.
[src cores §1.3]
[chk `check_fiveloop.log`, engine values on both sides, 0 mismatches]
- five-loop orbits: K5−e 112, W5 359, K4+K4 286, V7_1 503, V8_0 87 (total 1347)
- q=2: W5 → 70ζ7; K4+K4 → 36ζ3²

**(NT) two-point / numerator.** Place $a,b$ at $0,1$ and $c$ at ∞; the edge $ab$ drops out. $H'=\hat G-c-ab$, $\omega=\beta_{ab}$:
$$\Pi(\hat G)=\frac{\Gamma(q-\omega)}{\prod_{H'}\Gamma(\beta_e)}\int\prod_{H'}t_e^{\beta_e-1}\,\Psi_{H'}^{\omega-q}\,(\Psi^{ab}_{H'})^{-\omega}\,\Omega$$
[drv: $\int dt\,t^{\omega-1}(\Psi'+t\Psi^{ab})^{-q}=B(\omega,q-\omega)\Psi'^{\omega-q}(\Psi^{ab})^{-\omega}$, analytic in ω; $\Psi^{ab}$ = 2-forests separating $a,b$]
For $\omega=-k\le0$:
$$\Pi(\hat G)=\frac{\Gamma(q+k)}{\prod_{H'}\Gamma(\beta_e)}\sum_mc_m\,I_{q+k}(H';q+k-\beta-m),\qquad (\Psi^{ab}_{H'})^k=\sum_mc_mt^m,\ c_m>0$$
All terms are positive, so each term converges. A numerator edge becomes one loop fewer at level $q+k$.
[chk `check_two_point.log`: 1941 comparisons (k=1: 1408, 2: 418, 3: 79, 4: 32, 5: 4), 0 mismatches]
Several negative edges outside star(c) ∪ {ab}: Taylor coefficient with mixed signs [drv]; not implemented [open].

**Group on K4.** Labelled orbit of a generic K4 vector at fixed q (relabelling + completion + FT): 720 = |S6|; without FT 120 = |S5| (relabellings of $\hat G=K_5$); without completion 48 [chk `check_kite_group.log`]. Literature Z2×S6 of order 1440 for the two-loop master [src Broadhurst 1986, Barfoot–Broadhurst 1988, not reread]; the extra Z2 is not realized at fixed q [open].

## 3. Canonicalizer (`engine.period`)

### 3.1 Algorithm
- Input map: $I_q(G;\alpha)=\frac{\prod_{e\in G}\Gamma(q-\alpha_e)}{\Gamma(q)}\Pi(\hat G)$; edges with $\alpha_e=0$ are contracted first.
- State: canonical form of the weighted $\hat G$ (individualization–refinement).
- Search: best-first by $|V|$ over M1, M2 (exact triangles), FT (all c), TW (all 4-cuts and pairings).
- PR at a 3-cut with both sides ≥2 vertices: AND node, recursive, memoized.
- Terminal (a): some c with $\hat\beta|_{\hat G-c}\in[1,q-1]$ and `cores.period` evaluable (SP → banana, $W_n$ any q; prism and K33 only for q ≤ 7).
- Terminal (b): (NT) with k ≤ 3.
- Limits: budget 60 states, $|V|\le|V(\hat G_0)|+1$.
- Budget saturation (not logged): budget 400 or split-M2 changed reach on 30-sample subsets only for K33+e (11 → 12 of 30).

### 3.2 Coverage
Uniform samples of convergent $\alpha\in[1,q-1]^E$, $\sum\alpha=nq$. Columns "+NT": 500 samples with 95% Wilson interval. Columns "moves": 250 samples (the first 250 of the same sequence), NT off. q=7: D=2 bare; q=8: D=3 bare (five loops) and D=2 bare (K5). [chk `coverage_q{7,8}_k{3,0}.log`]

| core | V,E; planar; conn | q=7 +NT | q=7 moves | q=8 +NT | q=8 moves |
|---|---|---|---|---|---|
| K5−e | 5,9; P; 3 | 1.000 [.992,1] | 1.000 | 0.916 [.888,.937] | 0.600 |
| W5 | 6,10; P; 3 | 1.000 | 1.000 | 1.000 | 1.000 |
| K33+e | 6,10; N; 3 | 0.406 [.364,.450] | 0.340 | 0.246 [.210,.286] | 0.128 |
| Z5 | 6,10; P; 3 | 0.526 [.482,.569] | 0.412 | 0.488 [.444,.532] | 0.308 |
| K4+K4 | 6,10; P; 2 | 1.000 | 1.000 | 1.000 | 1.000 |
| V7_0 | 7,11; P; 3 | 0.412 [.370,.456] | 0.376 | 0.392 [.350,.435] | 0.196 |
| V7_1 | 7,11; P; 2 | 1.000 | 1.000 | 1.000 | 1.000 |
| V7_2 | 7,11; P; 3 | 0.288 [.250,.329] | 0.200 | 0.260 [.223,.300] | 0.148 |
| V7_3 | 7,11; N; 3 | 0.164 [.134,.199] | 0.144 | 0.106 [.082,.136] | 0.040 |
| V8_0 | 8,12; P; 2 | 1.000 | 1.000 | 1.000 | 1.000 |
| V8_1 | 8,12; P; 3 | 0.156 [.127,.190] | 0.092 | 0.156 [.127,.190] | 0.056 |
| V8_2 | 8,12; N; 3 | 0.056 [.039,.080] | 0.028 | 0.034 [.021,.054] | 0.008 |
| cube V8_3 | 8,12; P; 3 | 0.074 [.054,.100] | 0.048 | 0.062 [.044,.087] | 0.028 |
| Möbius–Wagner V8_4 | 8,12; N; 3 | 0.026 [.015,.044] | 0.024 | 0.014 [.007,.029] | 0.016 |
| K5 (6 loops) | 5,10; N; 4 | 0.426 [.383,.470] | 0.304 | 0.358 [.317,.401] | 0.148 |

Route statistics (q=7, +NT):
- K5−e: FT→prism 148, direct 352
- 2-connected cores: PR 323–368 of 500
- K5 (q=8): FT 86, M2 57, direct 34, FT+TW 2
- "direct" = completion and terminal only, no further move

Large q, 60 samples, q=20, +NT: K5−e 39, W5 60, K33+e 4, Z5 10, K4+K4 60, V7_0 12, V7_2 3, cube 2, Möbius 0, K5 6 [chk `check_large_q.log`; q=40 not run].

### 3.3 Exact values reached
| core | value | chain |
|---|---|---|
| K5−e | $I_7(4,6,3,2,4,4,6,4,2)=-\frac{269}{1108800}-\frac{2539}{2806650}\zeta_3+\frac{8}{6237}\zeta_5$ | FT → prism |
| K33+e | $I_7(3,3,4,6,4,2,5,3,4,1)=\frac{2550463}{2463436800}+\frac{328829}{119750400}\zeta_3-\frac{1391}{332640}\zeta_5$ | FT M1 FT FT → prism |
| Z5 | $I_7(6,6,2,3,6,1,4,1,5,1)=\frac{577499}{143700480}+\frac{10459}{935550}\zeta_3-\frac{5}{297}\zeta_5$ | FT M2 FT TW FT TW FT FT → prism |
| V7_0 | $I_7(1,6,5,5,1,3,2,2,6,2,2)=\frac{465089699}{1551965184000}+\frac{2578951}{3233260800}\zeta_3-\frac{545}{449064}\zeta_5$ | FT FT TW FT TW FT → NT(k=2) W4 |
| V7_2 | $I_7(6,6,1,5,1,5,2,2,4,2,1)=\frac{17455469}{387991296000}+\frac{159409}{1231718400}\zeta_3-\frac{463}{2395008}\zeta_5$ | … → NT(k=2) W4 |
| V7_3 | $I_7(5,5,3,4,4,1,1,2,3,2,5)=-\frac{2522141}{2874009600}-\frac{140131}{59875200}\zeta_3+\frac{37}{10395}\zeta_5$ | M2 M1 M1 M2 M1 FT M1 TW FT FT → prism |
| V8_1 | $I_7(5,2,4,1,1,3,6,1,1,3,5,3)=-\frac{139697}{8083152000}-\frac{1757}{32076000}\zeta_3+\frac{1}{12474}\zeta_5$ | … → prism |
| V8_2 | $I_7(2,3,2,6,2,2,4,5,1,5,2,1)=-\frac{152827}{1847577600}-\frac{3919}{17962560}\zeta_3+\frac{83}{249480}\zeta_5$ | … → prism |
| cube | $I_7(5,1,1,3,2,2,3,1,5,3,4,5)=-\frac{146563}{82771476480}-\frac{32891}{7185024000}\zeta_3+\frac{1}{142560}\zeta_5$ | … → NT(k=1) W4 |
| Möbius V8 | $I_7(1,2,3,3,5,4,4,3,2,5,2,1)=-\frac{110137}{16166304000}-\frac{2081}{112266000}\zeta_3+\frac{1}{35640}\zeta_5$ | … TW … FT → prism |
| K4+K4 | $I_4(1,1,1,1,3,1,3,3,3,3)=\frac{7}{960}+\frac{1}{40}\zeta_3$; $I_2(1^{10})=36\zeta_3^2$ | PR → K4·K4 |
| V7_1 | $I_4(1,1,1,1,3,1,3,3,1,3,2)=\frac{49}{86400}+\frac{7}{1800}\zeta_3+\frac{1}{150}\zeta_3^2$ | PR |
| V8_0 | $I_4(1,1,1,3,3,1,1,1,1,1,3,3)=\frac{7}{1920}+\frac{1}{80}\zeta_3$ | PR |
| W5 | $I_4(1^5;3^5)=-\frac{25255}{128}-490\zeta_3+350\zeta_5+420\zeta_7$ | tabulated |
| K5, q=8 | $I_8(5,5,6,1,6,4,5,6,5,5)=\frac{105943}{1569209241600}-\frac{17}{302702400}\zeta_3$ | 8 routes (FT→W3, FT→W4, M2→W4), equal |
| K5, q=8 | $I_8(4,4,6,7,6,1,5,5,5,5)=-\frac{62311}{224172748800}+\frac{1}{4324320}\zeta_3$ | 8 routes, equal |
| K5, q=7 | $I_7(4,2,5,5,5,6,5,2,5,3)=\frac{3251321}{64665216000}+\frac{52643}{449064000}\zeta_3-\frac{23}{124740}\zeta_5$ | 6 routes (FT→prism ×4, M2→K33 ×2), equal |

[chk `coverage_q7_k3.log`, `coverage_q8_k3.log`, `check_fiveloop.log`, `check_k5_examples.log`]

### 3.4 Consistency of engine values (identities not used by the engine)
- Independent routes: 900 samples (15 cores, q=5,7), 127 with ≥2 chains (different terminal tables), 0 disagreements [chk `check_routes.log`]. K5 q=8: 120 samples, 32 multi-route, 0 disagreements [chk interactive run; the logged examples are in `check_k5_examples.log`].
- Dimension shift: 1347 orbits, 0 mismatches (§2 DS) [chk].
- Z5 at q=4: 940/1956 orbits evaluated; 249 of these overlap the 260 orbits msys determines, all 249 equal [chk `check_fiveloop.log`].
- q=2 limits: W5 70ζ7 ✓; K4+K4 36ζ3² ✓; K33+e 36ζ3² directly. K33+e with $\beta\equiv1$ is a decompletion of the complement of $C_3\sqcup C_4$, which by PR has period $P(K_4)^2$ [drv]. P(Z5) not reached (§4.3).

### 3.5 Symmetry rows added to the IBP system (Z5)
| system | rank / unknowns | determined |
|---|---|---|
| msys q=4 | 1493/1956 | 260 |
| + engine values (940) | 1864/1956 | 972 |
| + 16093 symmetry equalities $I(Z_5;\alpha)=f\,I(Z_5;\alpha')$ | 1948/1956 | 986 |
| joint q=3,4 + 74 dimension links | 2022/2030 | 1016 (30/74, 986/1956) |

[chk `check_z5_closure_q4.log`, `check_z5_joint.log`]. P(Z5)=441/8 ζ7 not determined [open].

## 4. Structure

### 4.1 Which symmetry exchanges which move
$$\text{FT}:\ U_G\leftrightarrow\Psi_{G^*},\ \alpha\leftrightarrow\beta=q-\alpha,\ \text{S}\leftrightarrow\text{P},\ \text{Y (unique vertex)}\leftrightarrow\Delta\ (\text{unique triangle}),\ \text{contraction }(\alpha_e=0)\leftrightarrow\text{deletion }(\beta_e=0),\ p_\alpha\leftrightarrow p_\beta$$
$$\text{completion}:\ \text{S at }v\ (v\sim\infty)\ \text{and Y at a unique }v\ (v\not\sim\infty)\ \text{are one move M1 on }\hat G;\ \text{moving }\infty\ \text{exchanges them; P is trivial on }\hat G$$
Atlas circuit moves S, P, Y on core periods [drv]:
- S = series push-forward ($p_\alpha$ invariant)
- P = parallel push-forward ($G_{2q}$ factor)
- Y = star–triangle only at a unique vertex
- 2-vertex cut of H = 3-vertex cut of $\hat G$ = product (PR)

### 4.2 Loop number
$$|V(G^*)|=n(G)+1,\quad n(G^*)=|V(G)|-1\qquad(\text{Euler}),\qquad n(\hat G-c)=n(G)+\deg_{\hat G}\infty-\deg_{\hat G}c$$
The period is an integral over $(n-1)$ loop momenta or over $(|V|-2)$ positions, modulo scale. FT exchanges the two counts. Completion fixes $|V(\hat G)|-3$ and moves the loop number [drv]. Examples:
- prism (4 loops) ↔ K5−e (5 loops) under FT
- W4 ↔ W4 under hub ↔ ∞
- K33 → K5−e under (U)
- K5−e completion with one $\beta_{v\infty}=0$ → W4 (K6 minus two disjoint edges, decompleted at a degree-5 vertex)

### 4.3 Reachability
Trigger conditions [drv]:
- (T1) $\sum_{e\ni v}\alpha_e=(\deg v-2)q$, i.e. $\beta_{v\infty}=0$; for deg 3 this is the unique vertex.
- (T2) a triangle of $\hat G$ with $\sum\beta=q$.
- (T3) $w_a=w_b,\ w_c=w_d$ at a 4-vertex cut.
- (T4) a planar decompletion with $\beta\in[1,q-1]$.
- (T5) negative completion edges ⊂ star(c) ∪ {ab}, $k\le k_{\max}$.
- (T6) a 3-vertex cut of $\hat G$; always present when H is 2-connected but not 3-connected.

Reached for all convergent α [drv; chk tables]:
- $W_n$, K4, SP; prism and K33 for q ≤ 7.
- K5−e for q ≤ 7 (FT → prism).
- 2-connected skeletons K4+K4, V7_1, V8_0 (PR → completed pieces ≤ 5 vertices → K4).

Generic K5: $\hat G=K_6$ with all $\beta_{v\infty}\ne0$ and no triangle with $\sum\beta=q$.
- Every decompletion is K5 (non-planar).
- K6 has no 3- or 4-vertex cuts and no cubic vertex.
- (NT) gives K5−e at level q+k ≥ 8, not tabulated.
- Hence generic K5 is unreached [drv]. K5 reach ⊂ {∃v: $\sum_{e\ni v}\alpha_e=2q$} ∪ {∃ triangle of $\hat G$: $\sum\beta=q$} [drv; route statistics §3.2 consistent].

3-connected five-loop cores with |V| ≥ 6 other than W5: reach is a union of conditions of codimension ≥ 1 (T1–T3) plus bands $1\le k\le k_{\max}$ (T5). The fraction decreases with q: Z5 0.53 (q=7), 0.49 (q=8), 0.17 (q=20); Möbius V8 0.026, 0.014, 0 [chk `check_large_q.log`].

Zigzag, $q=2$, $\beta\equiv1$: $\hat G=C_7(1,2)$; every decompletion is Z5; Z5 is self-dual; twists map $C_7(1,2)$ to itself. Hence P(Z5)=441/8 ζ7 is a fixed point and is not reached [drv; chk engine None, `check_fiveloop.log` (B)].

Value spaces reached [drv from terminal tables]:
- 3-connected five-loop cores: $\mathbb Q+\mathbb Q\zeta_3+\mathbb Q\zeta_5$ through prism/K33/W4 terminals; $+\mathbb Q\zeta_7$ through W5 terminals (cube, q=8 example).
- 2-connected: products, e.g. $\zeta_3^2$.

## 5. Not done / failed
- Prism and K33 tables at q=8 not built (shared machine, load ≈ 16): D=3 coverage of K5−e is 0.916 instead of 1 [open].
- (NT) with ≥2 negative edges outside star(c) ∪ {ab} (mixed-sign Taylor coefficients) not implemented [open].
- P(Z5) at q=2, and Z5 q=4 (8-dimensional nullspace after symmetry rows), unresolved [open].
- The dimension-shift check needs every cotree term; partially covered cores were checked by routes and msys instead.
- Coverage is sampled (500 / 250 per core and q), not exhaustive; search budget 60.
- Twist only for double transpositions; FT only for planar 2-connected decompletions.
- Extra Z2 of the Z2×S6 kite group not identified [open].
