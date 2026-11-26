# Cores: exact periods of circuit-irreducible graphs

Tags: [chk] exact check run (script/log in this folder), [drv] derived here, [src] source (location), [open]. Paths relative to `explore/circuit-general/cores/`.

## 0. Kernels as core periods

$$K_\Gamma(\lambda;0)\Big[\prod_{\rm ins}y\Big]=\int_{\mathbb{RP}^{N-1}_{>0}}\frac{P(t)}{U^{q}}\,\Omega,\qquad P=U^{n+H}D\prod_{\rm ins}(Uy),\qquad q=n+H+m\quad(T=0)$$
[src `code/htk/CONVENTIONS.md` §2, `explore/nonSP/verify/k4_periods/README.md`; chk `numer.py`, `fastker.py`: interpolant = `htk.brute` at an independent point]

$$T=1,\ \kappa=0,\ (H,T)=(1,1),\ n=4:\qquad D=P\,U^{-9/2},\ \deg P=9\ \ (\text{prism});\ \text{smaller }q\text{ for degenerate densities}$$
[chk `numer.numerator_T1`, `t1_prism.log`, `t1_forests_rerun.log`]

$$I_q(G;\alpha)=\int\prod_e t_e^{\alpha_e-1}U_G^{-q}\,\Omega,\quad \sum_e\alpha_e=nq;\qquad \text{convergent}\iff\sum_{e\in S}\alpha_e>q\,h(S)\ \ \forall\,\emptyset\ne S\subsetneq E$$
[src k4_periods README; chk `gutil.convergent`]

$$p_q=\frac{I_q}{\prod_e\Gamma(\alpha_e)}:\quad\text{series: }p\text{ unchanged};\quad\text{parallel: }p=G_{2q}(s_1,s_2)\,p\big|_{s_1+s_2-q};\quad\alpha_e=0:\ p(G)=p(G/e);\quad \text{banana: }I_q=\frac{\prod_i\Gamma(q-s_i)}{\Gamma(q)}$$
$$G_{2q}(a,b)=\frac{\Gamma(q-a)\Gamma(q-b)\Gamma(a+b-q)}{\Gamma(a)\Gamma(b)\Gamma(2q-a-b)}$$
[src k4_periods README (Dirichlet push-forward); drv normalization; chk `check_c2_three_loop.log`, `check_c3_three_loop.log`, `check_fastker.log`]

Cores $n\le4$: $K_4=W_3$, $W_4$, prism, $K_{3,3}$ [src `explore/nonSP/verify/prism_forcer/README.md` §"Four-loop catalog"]. `atlas/census.json`: $(D,n)=(3,4)$: 39 graphs with prism or $K_{3,3}$ series core; $(4,3)$: 20 with $K_4$ [chk `d3_four_loop.py`, `d4_three_loop.log`].

## 1. Exact method

### 1.1 Wheels $W_n$, $n\ge3$
Radial–angular residue algorithm, integer powers, integer $q$ [src `explore/nonSP/verify/full_ward/README.md` §"Arbitrary integer powers on an n-spoke wheel"]; Fraction port `wheel.py` [chk `wheel_selftest.log`: radial masters $n\le10$; $(5,5,4,5,5;3,2,2,2,2;q=7)$: $-\frac{45001}{43545600}-\frac{116849}{29937600}\zeta_3+\frac{50837}{2993760}\zeta_5-\frac{7}{594}\zeta_7$].

### 1.2 Momentum-space triangle rule, $d=2q-2\epsilon$
$$(d-2a_1-a_2-a_3)I=a_2\big[I(1^-2^+)-I(2^+4^-)\big]+a_3\big[I(1^-3^+)-I(3^+5^-)\big]$$
zero power $=$ contraction; bubbles $G(a,b)$; $\Gamma(x+k\epsilon)=$ rational series $\times\,\Gamma(1+k\epsilon)$ or $\Gamma(\tfrac12+k\epsilon)$; truncated Laurent series with precision tracking [`prop.py`, `laurent.py`].
- $K_4$ ($q=2..7$, 31 vectors), $W_4$ spoke cut ($q=2..5$, 16): $=$ 1.1 [chk `check_momentum_vs_wheel.log`].
- Prism, rung cut (LA) and edge cut (BE): kite $K(a,b,c,d;\,n+\epsilon)$ (FORCER master t1star05), integer limit IR-divergent; reduction does not close [chk `check_t1star_obstruction.log`; src arXiv:1704.06650 §3.3].
- Same rule at $\epsilon=0$ on convergent periods only: no admissible rule for all 1481 bank monomials [chk `vac.py`].

### 1.3 Regulator-free relation systems (prism, $K_{3,3}$)
Relation used iff every period in it converges.
$$\text{(T)}\quad(2q-2\alpha_1-\alpha_2-\alpha_3)\,p(\alpha)=\alpha_2\big[p(\alpha-e_1+e_2)-p(\alpha+e_2-e_4)\big]+\alpha_3\big[p(\alpha-e_1+e_3)-p(\alpha+e_3-e_5)\big]$$
$A,B$ cubic; $L_1=AB$, $L_2=BC$, $L_3=AC$; $L_4,L_5$ third edges at $B,A$ [drv: cut integral with cut edge not at $A,B$; $\epsilon\to0$; $I_{\rm cut}=\Gamma(q)\,p$].
$$\text{(D)}\quad E\,p(\alpha)=\sum_i\alpha_{U_i}\big[p(\alpha+e_{U_i}-e_{D_i})-p(\alpha+e_{U_i}-e_{E_i})\big],\qquad E=(L+S)\,2q-\sum_i(\alpha_{U_i}+2\alpha_{D_i})-2\sum_j\alpha_{S_j}$$
$B$ with cubic neighbours $X_i$ ($D_i=BX_i$), common neighbour $T$ ($U_i=X_iT$), third edges $E_i$, $S_j=BT$ [src arXiv:1504.08258 eqs. 6–8; drv period level]; $K_{3,3}$: $B=a_1$, $X_i=b_i$, $T=a_2$, $E_i=b_ia_3$.
$$\text{(I)}\quad\sum_{C\not\ni i}I_q(b+1_C)=\frac{q-1-b_i}{b_i}\sum_{C\ni i}I_q(b+1_C),\qquad b\ge1,\ \textstyle\sum b=n(q-1),\ b_i<q-1$$
[src full_ward README §"Wheel integration by parts"; drv]
$$\text{(B)}\quad I_{q-1}(G/i;\,b_{\hat\imath})=(q-1)\sum_{C\ni i}I_q(b+1_C),\qquad b_i=0$$
[drv: $t_i=0$ boundary of $\int\partial_{t_i}\prod_{e\ne i}t_e^{b_e-1}U^{1-q}$]
$$\text{(S)}\quad\sum_{e\ni v}\alpha_e=q:\qquad I_q(G;\alpha)=I_q\big(G-v+\triangle;\ q-\alpha_1,q-\alpha_2,q-\alpha_3\big)$$
[src full_ward README ($x=S/r_1$, $U_\triangle=(x+y+z)U_{\rm star}$); drv $S^{\sum\alpha_\triangle-2q}=1$; chk `check_startri_bank.log` 805/805]

$C$: cotrees; boundary data: $W_n$ (1.1; half-integer $q$: 1.2), $K_4$ + series–parallel.

| core | $q$ | convergent orbits | relations | determined |
|---|---|---|---|---|
| prism | 3, 4, 5 | 6, 106, 951 | 28, 814, 9030 | all |
| prism | 6, 7 | 5416, 22789 | 59052, 273833 | all |
| prism | 7/2, 9/2 | 54, 615 | 297, 4622 | all |
| $K_{3,3}$ | 3, 4, 5 | 3, 35, 279 | 7, 210, 2226 | all |
| $K_{3,3}$ | 6, 7 | 1476, 5988 | 14158, 64967 | all |

[chk `check_closure_and_bank.log`, `tables_prism_q7.log`, `tables_K33_q7.log`, `t1_prism.log`; `tables/*.json`]

Checks:
- prism, $q=6$: 1481/1481 FORCER values [chk `check_closure_and_bank.log`; src `explore/nonSP/verify/prism_forcer/periods.json`]
- $I_{q-1}(\alpha)=\sum_CI_q(\alpha+1_C)$, not among the solving relations: prism 5416, $K_{3,3}$ 1476 orbits ($q=6\to7$); prism 54 ($\tfrac72\to\tfrac92$); $W_4$: 16 (I) and 2 dimension relations at $q=\tfrac72$ [chk `check_dimension_relation.log`, `t1_prism.log`, `check_halfq_wheel_ibp.log`]
- $\mathbb C^2$, $n=4$: prism $-\frac{577}{40320}-\frac{2}{35}\zeta_3+\frac{5}{63}\zeta_5$; prism (four marks) $\frac{169}{151200}-\frac{3}{2800}\zeta_3$; $K_{3,3}$ $-\frac{789}{640}+\frac{59}{80}\zeta_3$ [chk `check_c2_four_loop.log`, `check_fastker.log`; src nonSP README, k33_forcer README]
- $\mathbb C^2$, $n=3$: $-\frac{71}{12}$; $\frac{222317}{25200}-\frac{3021}{700}\zeta_3$; $-\frac{4003781}{302400}+\frac{3772}{525}\zeta_3$ [chk `check_c2_three_loop.log`; src nonSP README]
- $\mathbb C^3$, $n=3$: six $K_4$ subdivisions [chk `check_c3_three_loop.log`; src k4_periods README]
- reversal $(-1)^{Dk}$, edge order $\mathrm{sgn}^{D+1}$ [chk `check_sign_rules.log`; src CONVENTIONS §2.2]

### 1.4 Choice among (a), (b), (c)
- (a) hyperlogarithms: fibration basis over Dodgson-polynomial letters, $\le8$ variables, numerator degree $4q-9$ [drv]; not implemented [open].
- (b) Gegenbauer: $W_n$: one harmonic sum (1.1); prism, $K_{3,3}$: theta angular graph for every origin [drv]:
$$\sum_{m_1,m_2,m_3}\Delta_\lambda(m)\prod_i\Big(\frac{\lambda}{m_i+\lambda}\Big)^{L_i}\rho_i(m_i)\,R(m_1+m_2+m_3),\qquad \Delta_\lambda=\frac{\lambda}{s+\lambda}\,\frac{(2\lambda)_s}{(\lambda)_s}\prod_{x\in\{a,b,c\}}\frac{(\lambda)_x}{x!}$$
$s=\tfrac{m_1+m_2+m_3}{2}$, $(a,b,c)=(s-m_1,s-m_2,s-m_3)\in\mathbb Z^3_{\ge0}$ [chk `check_dougall.log`, $\lambda\le4$, $m_i\le5$]; exact triple summation [open].
- (c) at the level of scalar periods: $\mathbb Q$-linear algebra, wheel boundary data, closure 1.3 [chk].

### 1.5 Period spaces
$$I_q(\text{prism}),\ I_q(K_{3,3})\in\mathbb Q+\mathbb Q\zeta_3+\mathbb Q\zeta_5\quad(q=3,\dots,7);\qquad I_q(\text{prism})\in\pi^2\mathbb Q+\pi^4\mathbb Q\quad(q=\tfrac72,\tfrac92)$$
[drv: rational relations, boundary values in these spaces; chk tables]; uniform in $q$ [open].

## 2. Wheels and zigzags
$$I_q\big(W_n;\ q-1,\dots,q-1;\ 1,\dots,1\big)=\binom{2n-2}{n-1}\sum_{s\ge1}\frac{M_q(s)}{s^{2n-1}}=\binom{2n-2}{n-1}\frac{2}{(2q-2)!}\sum_{k=1}^{q-1}c_{q,k}\,\zeta(2n-1-2k)$$
$$M_q(s)=\frac{2s^2\prod_{j=1}^{q-2}(s^2-j^2)}{(2q-2)!},\qquad\prod_{j=1}^{q-2}(x-j^2)=\sum_{k}c_{q,k}x^{k-1}$$
[src full_ward README ($M_q$); drv closed form; chk `check_radial_closed_form.log`, 55 pairs, $3\le n\le12$]; $q=2$: $\binom{2n-2}{n-1}\zeta(2n-3)=P(WS_n)$; $n=5$: $70\zeta_7$, $\frac{35}{6}(\zeta_5-\zeta_7)$, $\frac{7}{36}(\zeta_3-5\zeta_5+4\zeta_7)$ [chk].
- Wheel-cored holomorphic graphs: $W_n$ with $n(H-2)+1$ subdivision vertices; $H=2$: one subdivided spoke or rim [drv: $|E|=nH+1$, loop number $n$]; every source, insertion: $K\in\mathbb Q+\sum_{j=2}^{n-1}\mathbb Q\zeta_{2j-1}$ [src full_ward README; drv path push-forward].
- Parity: summand $M_q(s)\times$(odd in $s$); $\lambda=q-1\in\mathbb Z$: $M_q$ even $\Rightarrow$ odd $\zeta$ [src full_ward README]; $\lambda\in\tfrac12+\mathbb Z$ ($T=1$): multiplicity odd, $\sum_{s\in\frac12+\mathbb N}s^{-w}=(2^w-1)\zeta(w)$, $w$ even $\Rightarrow\pi^{2k}$ [drv]; $W_4$, $q=\tfrac52$, $(1,1,1,1,1,2,1,2)$: $\frac{56}{3}\pi^2-\frac{16}{9}\pi^4$ [chk `prop.cut_period`, `check_halfq_wheel_ibp.log`].
- Zigzags: $Z_3=K_4$, $Z_4=W_4$ ($C_{n+2}(1,2)$ minus a vertex) [chk `cores.zigzag`]; $P(Z_n)=\frac{4(2n-2)!}{n!(n-1)!}\big(1-\frac{1-(-1)^n}{2^{2n-3}}\big)\zeta(2n-3)$ [src arXiv:1208.1890]; $Z_5$: rank 25/74 ($q=3$), 1493/1956 ($q=4$) [chk `zigzag_z5_q3.log`, `zigzag_z5_q4.log`]; joint $q=3..5$ with dimension links unsolved [open, `zigzag_multilevel.log`]; holomorphic zigzag-cored graphs at $H=2$: one subdivision [drv]; values [open].

## 3. Structure
- Series–parallel: Dirichlet integrals; $T=0$: $K\in\mathbb Q[\lambda]$; $T=1$, $n$ even: $K\in\pi^{n/2+2E}\mathbb Q[\lambda]$ [src `explore/ht-families/README.md` §2].
- Cores: $I_q(G;\alpha)$ = period of $\mathbb P^{N-1}\setminus\{U_G=0\}$ relative to the coordinate simplex, polynomial numerator; integer powers and numerators add no divisor [src `prism_forcer/README.md`, arXiv:1305.2161]; linear reducibility $\Rightarrow$ MZV [src arXiv:0910.0114]; prism cuts MZV, $K_{3,3}$ cut alternating sums admitted [src arXiv:1305.2161 Thm 1.1]; computed $n\le4$: 1.5 [chk].
- Non-MZV: $\phi^4$ counterexamples to Kontsevich's conjecture, point counts by a modular form of a singular K3 [src arXiv:1006.4064]; 8 loops, weight 3, level 7; first quasi-constant $c_2$ at 7 loops [src arXiv:1304.5342]; $\phi^4$ data complete to 7 loops [src arXiv:1603.04289]. Holomorphic cores: same hypersurface $U_G=0$ [drv]; non-MZV expected first for 8-loop cores [drv conditional on src]; weighted, numerated periods [open].
- $T=1$, prism at $(H,T)=(1,1)$, $\kappa=0$:
$$K(1,2,-1,3,-2,-3)=-\tfrac{808}{63}\pi^2+\tfrac{2041}{1680}\pi^4,\qquad K(1,1,0,-1,0,-1)=-\tfrac{33}{70}\pi^2+\tfrac{17}{360}\pi^4,\qquad K(1,0,-1,0,0,0)=0$$
[chk `t1_prism.log`; float MC $-8.50\pm0.34$ vs $-8.2411$, `mc_t1_prism.log`]; field $\mathbb Q(\pi)$; non-series–parallel: $\pi^2\mathbb Q+\pi^4\mathbb Q$; series–parallel: one power [chk].
- $(1,1)$ product forests, $n=4$, sources of `ht-d2/nonsp_numeric.py` (roots $2,-1$; $\lambda_A=-2$) [chk `t1_forests.log`, `t1_forests_rerun.log`, `t1_tree_check.log`]:

| forest | $K$ | $K/\pi^2$ | ht-d2 MC |
|---|---|---|---|
| $(A,(B,(A,B)))(A,B)$ | $-\frac{19}{105}\pi^2+\frac{1}{35}\pi^4$ | $0.101036$ | $0.10195(42)$, $0.10111(49)$ |
| $(B,(A,(A,B)))(A,B)$ | $\frac{23}{105}\pi^2-\frac{1}{35}\pi^4$ | $-0.062941$ | $-0.06289(34)$, $-0.06053(218)$ |
| $(B,(A,(A,(A,B))))$ | $-\frac{8}{315}\pi^2$ | $-0.025397$ | $-0.02544(8)$ |
| $(B,(A,(B,(A,B))))$ | $-\frac{4}{315}\pi^2$ | $-0.012698$ | $-0.01414(81)$ |
| $(B,(B,(A,(A,B))))$ | $\frac{52}{315}\pi^2$ | $0.165079$ | $0.16800(244)$, $0.16449(57)$ |
| $((B,(A,B)),(A,B))$ | $-\frac{64}{315}\pi^2$ | $-0.203175$ | $-0.20343(419)$ |
| $((A,(A,B)),(A,B))$ | $-\frac{8}{315}\pi^2$ | $-0.025397$ | $-0.02514(9)$ |
| $((A,B),(A,B))(A,B)$, $(A,((A,B),(A,B)))$, $(B,((A,B),(A,B)))$ | $0$ | $0$ | $\approx0$ |

Series–parallel rows: $=$ ht-d2 exact table (10 forests) [chk]. $((A,(A,B)),(A,B))$, $(l_0,l_A)=(1,0)$: $-\frac{4}{315}\pi^2$ (MC $c_4=-0.01263(9)$); $(1,1)$: $-\frac{67}{630}\pi^2+\frac{3}{560}\pi^4=-0.053476\,\pi^2$ (MC $\sum c_k=-0.05315$) [chk `t1_tree_check.log`].

## 4. Evaluator
- `evaluate.py kernel NV EDGES SOURCES [INSERTIONS] [T]`, `evaluate.py period EDGES POWERS Q`; API `cores.kernel`, `cores.period`, `fastker.kernel` (bare, $T=0$, modular Newton interpolation, CRT).
- Period: any graph whose series–parallel reduction ends in a banana, $W_n$ ($n\ge3$), prism or $K_{3,3}$; integer $q$; half-integer $q$ for wheels and prism; prism/$K_{3,3}$ tables built on demand.

| class | route | time |
|---|---|---|
| $W_n$ periods | `wheel.py` | ms |
| prism tables | `msys` | $q=6$: 21 s; $q=7$: 110 s; $q=\tfrac72$: 5 s; $q=\tfrac92$: 131 s |
| $K_{3,3}$ tables | `msys` | $q=6$: 330 s; $q=7$: 201 s |
| $n\le3$, $D=2$, insertions | `cores.kernel` | 0.3–2 s |
| $n=3$, $D=3$ | `fastker` / `cores.kernel` | 0.1 s / 8 s |
| $n=3$, $D=4$ | `fastker` | 2 s (census[4,3][22–24]) |
| $n=4$, $D=2$ | `fastker` / `cores.kernel` | 2 s / 45–70 s |
| $n=4$, $D=3$, bare | `fastker` | 154–156 s, 4 primes, 1.5 GB |
| $n=4$, $T=1$, $(1,1)$ | `cores.kernel(T=1)` | 80–120 s |

$D=3$, $n=4$ [chk `d3_four_loop.log`]: census[3,4][191] (prism core): $\frac{539994681409}{634894848000}+\frac{1638600697}{734832000}\zeta_3-\frac{2779639}{816480}\zeta_5$; census[3,4][219] ($K_{3,3}$ core): $\frac{14383751986301}{317447424000}+\frac{89673345011}{734832000}\zeta_3-\frac{302364397}{1632960}\zeta_5$; sources $(1,0,0),(0,1,0),(0,0,1),(1,1,-1),(1,-1,1),(-1,1,1),(2,0,-1),(0,-2,1),(-1,2,0),(-3,-2,-2)$.

## 5. Not done
- $D=3$, $n=4$ with insertions: $q\ge8$ tables and modular tree-current product in `fastker` [open].
- $n=3$, $D\ge5$: dense-cube interpolation ($3^{N-1}$) infeasible; core periods exact [open].
- Exact Gegenbauer triple sums; hyperlogarithm integration [open].
- Five-loop cores other than $W_n$; $P(Z_5)$ via relations [open].
- $K_{3,3}$ at half-integer $q$: not built [open].
- $T=1$: individual numerator monomials can diverge; orbit sums cancel in the computed case; general statement [open].
- `ibp_system.py`, `check_prism_bank.py`, `vac.py`: superseded or negative results.
