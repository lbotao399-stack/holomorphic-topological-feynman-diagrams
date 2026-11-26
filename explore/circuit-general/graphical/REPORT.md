# Graphical functions in d = 2q: exact core periods beyond circuit moves

Tags: [chk] exact check (script/log in `explore/circuit-general/graphical/`), [drv] derived here, [src] source, [open]. Paths are relative to `explore/circuit-general/graphical/`. Tables are in `tables/<core>_q<q>.json`, as `[alpha, {label: value}]` with labels `'1','z3','z3*z3',…` (same format as `cores/tables`). Edge orders are those of `fiveloop/skel5.json` and `cores`. K5 edge order: 01,02,03,04,12,13,14,23,24,34.

## 0. Dictionary

$$A_G:=\int\prod_{v\ne h_0,h_1}\frac{d^{2q}x_v}{\pi^q}\prod_e|x_e|^{-2\beta_e},\qquad x_{h_0}=0,\ |x_{h_1}|=1,\qquad \beta_e=q-\alpha_e,\quad \lambda=q-1$$
$$|x|^{-2\beta}=\tfrac1{\Gamma(\beta)}\!\int_0^\infty\!a^{\beta-1}e^{-ax^2}da,\qquad \int\prod_{v\ne h_0}\tfrac{d^dx_v}{\pi^{q}}e^{-\sum a_ex_e^2}=\Phi_G^{-q},\quad \Phi_G=\sum_T\prod_{e\in T}a_e$$
$$\tfrac{d\lambda}{\lambda}=2\tfrac{d\rho}{\rho}:\ \ \int\prod a_e^{\beta_e-1}\Phi_G^{-q}\Omega=\tfrac{\prod\Gamma(\beta_e)}{\Gamma(q)}A_G;\qquad a_e=t_e^{-1},\ \Phi_G(t^{-1})=U_G(t)\textstyle\prod t_e^{-1}$$
$$\boxed{I_q(G;\alpha)=\frac{\prod_e\Gamma(q-\alpha_e)}{\Gamma(q)}\,A_G}$$
[drv; chk banana, §4; equals symmetry/REPORT §1.1]. Borinsky–Schnetz: $\nu_e=\beta_e/\lambda$, $f^{(\lambda)}_G(z)=A_G(0,e,z)$, $P^{\rm BS}_G=A_G$ [src arXiv:2105.05015 eqs. (3),(12)].

Chains: by the series rule ($p_\alpha$ invariant, cores §0), an edge with weight $\alpha_e$ equals a chain of $\alpha_e$ edges of weight $\beta=\lambda$ ($\nu=1$):
$$A_G=A_{G_{\rm sub}}\,\Gamma(\lambda)^{\alpha_e}\Gamma(\alpha_e)/\Gamma(\beta_e)$$
[drv; chk $\alpha_e=2$: $G_{2q}(\lambda,\lambda)=\Gamma(\lambda-1)/\Gamma(\lambda)^2$].

## 1. Single-valued polylogarithms, MZVs (`mzv.py`, `rat.py`, `sv.py`, `expand.py`)

$$\partial_zL_{wa}=\tfrac{L_w}{z-a},\quad L_{0^k}=\tfrac{\log^kz}{k!},\quad \partial L=L\big(\tfrac{x_0}{z}+\tfrac{x_1}{z-1}\big),\quad Z=\sum\zeta_ww=L(1)_{\rm reg},\quad \zeta_{1\,0^{m_1}\cdots1\,0^{m_d}}=(-1)^d\zeta(m_d{+}1,\dots,m_1{+}1)$$
$$P(z)=\tilde L_{x_0,x_1'}(\bar z)\,L_{x_0,x_1}(z),\qquad Z_{-x_0,-x_1'}\,x_1'\,Z^{-1}_{-x_0,-x_1'}=Z\,x_1Z^{-1}$$
[src arXiv:1302.6445 §2.6 (2.30); drv the monodromy derivation: $M_0L=e^{2\pi ix_0}L$, $M_1L=Ze^{2\pi ix_1}Z^{-1}L$, $\tilde Z_{x_0,x_1}=Z^{-1}_{-x_0,-x_1}$]

$$\partial_zP_{wa}=\tfrac{P_w}{z-a},\qquad \partial_{\bar z}P=\big(\tfrac{x_0}{\bar z}+\tfrac{x_1'}{\bar z-1}\big)P,\qquad P_uP_v=P_{u\,ш\,v},\qquad P_w(1)_{\rm reg}=[\tilde Z'Z]_w$$
$$x_1'-x_1=\zeta_3(-2x_0^2x_1^2+4x_0x_1x_0x_1+2x_0x_1^3-4x_1x_0x_1x_0-6x_1x_0x_1^2+2x_1^2x_0^2+6x_1^2x_0x_1-2x_1^3x_0)+\dots$$
[chk `test_sv.log`]

Operations:
- $\int_{\rm sv}dz$: partial fractions + IBP, $\int P_w/(z-a)=P_{wa}$.
- $\int_{\rm sv}d\bar z$: $\int P_w/(\bar z-a)=P_{aw}-\sum_v c_v\int P_{(aw)\setminus v}/(\bar z-1)$.
- $\sigma P_w=\sum_s[x_1'\text{-subst}(\tilde s)]_wP_s$.
- $P_w(1-z)$ and $P_w(1/z)$: re-integration, constants from $P(1)_{\rm reg}$.
- Expansions at $0,1,\infty$ (log-Laurent, total-degree truncation).

Rational coefficients: canonical $N/(z^a\bar z^b(z-1)^c(\bar z-1)^d(z-\bar z)^e)$. Exact Fractions throughout.

MZV algebra: double shuffle + Hoffman relations, weights ≤ 11. The free dimensions are 1,0,1,1,1,2,2,3,4,5,7,9 (asserted). Basis: $\zeta_2,\zeta_3,\zeta_5,\zeta_7,\zeta_9,\zeta_{11},\zeta(5,3),\zeta(7,3),\zeta(5,3,3)$.

| check | result | log |
|---|---|---|
| 14 reduced MZVs (w ≤ 11) vs 40-digit numerics | $\le1.2\cdot10^{-39}$ | `mzv.log` |
| $P_{01}-P_{10}=4iD$ | $9\cdot10^{-32}$ | `test_sv.log` |
| $\sigma P_w(z)=P_w(\bar z)$, $\|w\|\le6$ | $5\cdot10^{-40}$ | same |
| $P_w(1-z)$ expansion, $\|w\|\le6$ (tests single-valuedness and $x_1'$) | $2.5\cdot10^{-35}$ | same |
| $P_w(1/z)$, $\|w\|\le5$ | $4\cdot10^{-20}$ | same |
| $\partial,\bar\partial$ numerical; $\partial\!\int=$ id | $8\cdot10^{-24}$; exact | same |

## 2. Appending an edge in d = 2λ+2

### 2.1 Equation
$$z=x_1+i\rho,\ w=z-\bar z:\qquad \tfrac14\Delta_d=M_\lambda:=\partial\bar\partial-\lambda T,\qquad T:=w^{-1}(\partial-\bar\partial)$$
[drv]
$$\Delta_d|x|^{-2\lambda}=-\tfrac{4\pi^{d/2}}{\Gamma(\lambda)}\delta^d\ \Rightarrow\ M_\lambda f_{G_1}=-\tfrac1{\Gamma(\lambda)}f_G$$
[src BS (7)–(9)]

### 2.2 Identities
$$[\partial\bar\partial,T]=T^2,\qquad M_{\mu+1}T=TM_\mu,\qquad M_1=w^{-1}\partial\bar\partial\,w,\qquad T^n\big(h(z)/w\big)=w^{-(n+1)}d_nh,\quad d_n=w^{n+1}\partial^nw^{-(n+1)}$$
$d_n$ is the BS kernel operator [chk `check_identities_2.log`; src BS Thm 33].
$$M_\lambda(\partial\bar\partial)^{\lambda-1}w^{-1}=w^{-1}(\partial\bar\partial)^\lambda,\qquad (\partial\bar\partial)^mQ^{m-1}\log Q=(m-1)!\,m!\,\rho^{2m}Q^{-m-1},\qquad \int\tfrac{d^dy}{\pi^{d/2}}=\tfrac1{\Gamma(\lambda)}\int\tfrac{d^4y}{\pi^2}\rho^{2\lambda-2}$$
[chk `check_identities_1.log`, λ ≤ 4]. Hence
$$f_{G_1}=-\tfrac1{\Gamma(\lambda)}(\partial\bar\partial)^{\lambda-1}w^{-1}(\partial\bar\partial)^{-\lambda}(wf_G)+\ker M_\lambda$$
This form needs $wf_G$ pole-free in $w$. After one append, $f$ has poles up to $w^{-(2\lambda-1)}$, so the next step meets $\int dz\,w^{-1}\Phi$, i.e. the letter $\bar z$ [drv]; BS: GSVHs in letters $0,1,\bar z$ [src BS Thm 45].

### 2.3 Pole-free representation $V_d$
$$f=\sum_{k=0}^{\lambda-1}T^k\big(w^{-1}\Phi_k\big),\qquad \Phi_k\ \text{SV in letters }0,1,\ \text{no }(z-\bar z)\text{ poles}$$
$$M_\lambda\sum_kT^kX_k=\sum_kT^k\big[M_1X_k-(\lambda-k)X_{k-1}\big]\ \Rightarrow\ X_k=M_1^{-1}\Big[(\lambda-k)X_{k-1}-\tfrac{A_k}{\Gamma(\lambda)}\Big],\qquad M_1^{-1}(w^{-1}\Psi)=w^{-1}(\partial\bar\partial)^{-1}_{\rm sv}\Psi$$
[drv]. Only pole-free SV double primitives occur, so letters stay in $\{0,1\}$.

Edge factors ($r$ symmetric rational; $T^mr$ has no $w$-pole):
$$rT^kA=\sum_j\tbinom kj(-1)^{k-j}T^j\big((T^{k-j}r)A\big)$$
[drv]. Numerator edges ($\beta\le0$) are allowed.

Swap 0↔1: $\Phi_k(z)\to-\Phi_k(1-z)$ [drv].

Inversion: Kelvin gives $M_\lambda K=(z\bar z)^{-2}KM_\lambda$, but $TK\ne cKT$ [chk `check_kelvin.log`]. So relabelings that move $z$ are not available in $V_d$.

### 2.4 Kernel
$$\Phi_k\to\tfrac12(\Phi_k-\sigma\Phi_k);\qquad f\to f+\sum_\varphi c_\varphi T^{\lambda-1}\big(w^{-1}(\varphi(z)-\varphi(\bar z))\big),\quad \varphi\in\{z^{-j},(z-1)^{-j},z^j\}$$
These are the symmetric part of BS $d_n\varphi+\bar d_n\bar\varphi$ [drv]. Conditions [src BS Thm 5 (G3), eqs. (17)–(20); drv power counting]:
- Regularity on $z=\bar z$: the homogeneous parts of $w^Ef$ at $0,1,\infty$ are exactly divisible by $w^E$.
- $m,\bar m\ge M_s$ ($s=0,1$) and $m,\bar m\le M_\infty$, with

$$M_s=-\max_V\Big[\!\sum_{e\subseteq V\cup\{s,z\}}\!\beta_e-(\lambda+1)|V|\Big],\qquad M_\infty=-\min_V\Big[\!\sum_{e\cap(V\cup\{z\})\ne\emptyset}\!\beta_e-(\lambda+1)|V|\Big]$$

Any undetermined direction must be the zero function (runtime check). The DP for the bounds equals brute force on 300 random cases [chk].

d=6 star (2,2,2) gives exactly $1/(z\bar z(z-1)(\bar z-1))$ [chk `check_basic.log`; src BS (13)].

## 3. Periods, completion, plans (`gf.period`, `seq.py`)
$$\int\tfrac{d^dz}{\pi^{d/2}}f=(-1)^\lambda\tfrac{\lambda!}{(2\lambda)!}\int_{\mathbb C}w^{2\lambda}f\,\tfrac{d^2z}{\pi}$$
[drv; = BS (35)]
$$\int_{\mathbb C}\Phi\tfrac{d^2z}{\pi}=c_\infty-c_0-c_1,\quad H=\int_{\rm sv}\Phi\,dz,\quad c_s=[(z-s)^0(\bar z-\bar s)^{-1}]H,\ c_\infty=[z^0\bar z^{-1}]H$$
[drv Stokes; $w^{2\lambda}f$ is pole-free]
$$\hat G:\ \beta_{v\infty}=2q-\textstyle\sum_{e\ni v}\beta_e,\qquad A_{\hat G-c}=A_G\ \ \forall c$$
[src BS §2.3, symmetry REPORT §1.2; chk §4]

Sequential plan: hubs $h_0,h_1\to0,1$; $(\hat G-c)-\{h_0,h_1\}$ is a path; path edges have $1\le\beta\le\lambda$ (chain of $q-\beta$ appends); hub edges may have any integer $\beta$. Convergence conditions: $-M_\infty>1$ before each append; at the end $-M_\infty>\lambda+1$ and $M_s+\lambda+1>0$ [drv].

Vertex elimination with function-valued response: each append integrates one vertex whose star has legs $0,1$ (factors), the previous vertex (response $f(z)$) and the new $z$ [drv].

## 4. Validation (exact agreement with independent values)

| object | q | orbits agreeing | reference | log |
|---|---|---|---|---|
| W3, W4, W5, W6, Z5, unit weights | 2 | 5/5: $6\zeta_3$, $20\zeta_5$, $70\zeta_7$, $252\zeta_9$, $\tfrac{441}8\zeta_7$ | [src 1208.1890, cores §2] | `check_basic.log` (rerun by the orchestrator 2026-11-23, identical) |
| K4 | 3, 4, 5 | 2, 9, 25 | `cores.period` (wheel.py) | `logs_validate_K4q3.log`, `logs_validate1.log` |
| W4 | 3, 4 | 9, 89 | wheel.py | `logs_validate1.log` |
| W5 | 3, 4 | 16, 342 | wheel.py | `logs_validate_W5q3.log`, `logs_validate_W5q4.log` |
| prism | 3; 6; 7 | 6; 3 samples; 3 samples | msys tables | `logs_validate1.log`, `logs_highq.log` |
| completion route: K33, K4, W4, prism | 3, 4 | K33 q=3: 1, q=4: 10; K4 q=4: 9; W4 q=3: 9; prism q=3: 6 | cores tables | `logs_validate_completion.log` |
| Z5; V7_0; V7_2; V8_1; cube V8_3 | 3 | 74; 64; 52; 4; 1 | fiveloop relation tables | `logs_tables_q3.log`, `logs_tables_V_q3.log` |
| Z5 | 4 | 16 (run stopped at 16; fiveloop table has all 1956) | fiveloop | `logs_Z5_q4_part1.log` |
| K33+e | 3, 4 | 6, 111 | fiveloop | `logs_tables_q3.log`, `logs_K33e_q4.log` |
| K5 (three engine values, §7) | 8, 8, 7 | 3/3 | symmetry engine | `logs_k5_ref.log` |

Dimension relation $I_{q-1}(a)=\sum_C I_q(a+1_C)$ [chk]:
- Z5, q=3→2: $\tfrac{441}{8}\zeta_7$
- W5, q=3→2: $70\zeta_7$
- K5−e, q=4→3: 6/6

Plan independence (structurally different operation sequences) [chk]:
- K5, q=3: 6 plans
- six-loop q=3: 259/259
- K5, q=8: 2 values × 3 plans
- Z5, q=7: 2 plans

0 mismatches overall.

## 5. Coverage (sequential directly, or after completion + decompletion)

Five loops, plus four-loop cores. Entries are plan / convergent orbits [chk `coverage_q34.log`]:

| core | sequential as a graph | q=3 | q=4 |
|---|---|---|---|
| K4, W4, prism | yes | 2/2, 9/9, 6/6 | 9/9, 89/89, 106/106 |
| K33 | no | 1/3 | 10/35 |
| K5−e | yes | 6/6 | 106/106 |
| W5 | yes | 16/16 | 342/342 |
| Z5 | yes | 74/74 | 1956/1956 |
| V7_0, V7_2 | yes | 64/64, 52/52 | 3147/3147, 2301/2301 |
| K33+e | no | 6/18 | 111/429 |
| K4+K4 | no | 3/11 | 56/274 |
| V7_1, V7_3 | no | 2/10, 19/65 | 56/493, 516/2862 |
| V8_0, V8_1, V8_2 | no | 0/0, 4/22, 3/13 | 4/87, 230/2161, 135/1157 |
| cube V8_3, Möbius V8_4 | no | 1/8, 7/19 | 20/465, 183/1398 |

Six loops: 32 class-I graphs, 10 series-reduced skeletons; none is sequential as a graph [chk `coverage6_q3.log`, `coverage6_q4_*.log`]:

| census indices | V, E | q=3 | q=4 |
|---|---|---|---|
| 146,153,306,347,357,383,434,444,463,467 | 6, 11 | 20/29 | 675/1229 |
| 309,350,382 (K5) | 5, 10 | 1/2 | 15/30 |
| 338,345,431,453,497,558,566 | 7, 12 | 92/257 | — |
| 358,384,441,515,539,560 | 7, 12 | 22/122 | — |
| 583; 588; 591; 594; 601; 603 | 8, 13 | 27/353; 36/356; 38/430; 0/34; 22/348; 1/95 | — |

## 6. New exact values (`tables/`)
- K5−e: q=3, 6/6; q=4, 106/106. Value space $\mathbb Q+\mathbb Q\zeta_3+\mathbb Q\zeta_5$ [chk].
- V7_1 2, V7_3 19, V8_2 3, Möbius V8_4 7, K4+K4 3 orbits at q=3; K4+K4 56 orbits at q=4 [chk]. Example: $I_3(V8_4;1^9,2,2,2)=\tfrac53\zeta_3-\tfrac53\zeta_5$.
- Six-loop K5-cored, q=3: 259 orbits, each confirmed by a second plan (`tables/K5core_<idx>_q3.json`). Value space up to weight 9: $1,\zeta_3,\zeta_5,\zeta_7,\zeta_9,\zeta_3^2,\zeta_3\zeta_5,\zeta_3^3$. Example: census-338 skeleton, $\alpha=(1,1,2,1,2,1,2,2,2,1,1,2)$:
$$\zeta_3-\tfrac{13}6\zeta_5+2\zeta_3^2+\tfrac{28}3\zeta_7-\tfrac23\zeta_3^3-\tfrac{1063}{108}\zeta_9$$
[chk `logs_table6_q3.log`]
- $I_3(K_5;1,1,2^8)=\tfrac53\zeta_3-\tfrac53\zeta_5$ [chk `logs_K5_q3_plans.log`, 6 plans].
- Z5 at q=7, not reached by the symmetry engine [chk `logs_Z5q7_engine.log`, `logs_Z5q7_plans.log`]:
$$I_7(Z_5;4,3,3,4,3,4,3,3,4,4)=\tfrac{109885461859}{1385683200}-\tfrac{687732763}{997920}\zeta_3-\tfrac{11385931}{18480}\zeta_5+\tfrac{550599}{400}\zeta_7$$
- The Z5 q=3 table (74 orbits) is independent of the fiveloop joint relations, which used $P(Z_5)|_{q=2}$ as input; it agrees with them on all 74 orbits.

## 7. K5 by vertex elimination with function-valued response (coordinator)

Reference values, $\alpha$ in the K5 edge order, all confirmed [chk `logs_k5_ref.log`]:

| value | $\beta_{v\infty}$ | time | engine (symmetry §3.3) |
|---|---|---|---|
| $I_8(5,5,6,1,6,4,5,6,5,5)=\frac{105943}{1569209241600}-\frac{17}{302702400}\zeta_3$ | (1,4,6,5,0) | 7 appends, 143 s | equal |
| $I_8(4,4,6,7,6,1,5,5,5,5)=-\frac{62311}{224172748800}+\frac1{4324320}\zeta_3$ | (5,0,4,1,6) | 80 s | equal |
| $I_7(4,2,5,5,5,6,5,2,5,3)=\frac{3251321}{64665216000}+\frac{52643}{449064000}\zeta_3-\frac{23}{124740}\zeta_5$ | (2,6,0,2,4) | 61 s | equal |

Reach: a sequential plan exists ⟺ T1 ($\exists v:\beta_{v\infty}=0$, i.e. $\sum_{e\ni v}\alpha_e=2q$). This was checked on every convergent orbit at q=3..6 [chk `k5_reach_fraction.log`]:

| q | orbits | T1 (= GF reach) | T2 (triangle $\sum\beta=q$) | T1 ∪ T2 |
|---|---|---|---|---|
| 3 | 2 | 1 | 2 | 2 |
| 4 | 30 | 15 | 26 | 28 |
| 5 | 312 | 133 | 277 | 289 |
| 6 | 2107 | 792 | 1749 | 1876 |

Values at q=8 that the symmetry engine (budget 60) does not reach [chk `logs_k5_missed_q8.log`; the first two confirmed by 3 plans each, `logs_k5_missed_plans.log`]:
$$I_8(5,6,4,1,6,5,6,7,2,6)=\tfrac{263371}{7472424960000}-\tfrac{1}{34594560}\zeta_3\quad(\beta_{4\infty}=-1:\ \text{numerator edge})$$
$$I_8(3,3,7,5,5,4,7,3,5,6)=-\tfrac{8711086259}{3923023104000}-\tfrac{26946407}{4540536000}\zeta_3+\tfrac{5689}{630630}\zeta_5$$
$$I_8(2,5,7,6,6,5,4,7,4,2)=\tfrac{378433667}{4707627724800}+\tfrac{21719}{100900800}\zeta_3-\tfrac{5}{15288}\zeta_5$$
$$I_8(5,5,7,1,7,2,7,6,4,4)=-\tfrac{735555239}{201180672000}-\tfrac{17797783}{1816214400}\zeta_3+\tfrac{3755}{252252}\zeta_5$$

Generic K5 (no T1): the one-variable class is insufficient [drv]. $\hat G=K_6$ with all 15 $\beta\ne0$. For any frame $\{0,1,\infty\}$, the remaining three vertices $a,b,c$ form a triangle. The first elimination is
$$F_c(x_a,x_b)=\int\frac{d^dx_c}{\pi^q}|x_c|^{-2\beta_{c0}}|x_c-e|^{-2\beta_{c1}}|x_c-x_a|^{-2\beta_{ca}}|x_c-x_b|^{-2\beta_{cb}},\qquad \text{leg sum}=2q-\beta_{c\infty}\ne2q$$
- This is the one-loop non-conformal four-point star, i.e. a conformal 5-point function with weight $\beta_{c\infty}$ at ∞.
- It depends on $|x_a|^2,|x_a-e|^2,|x_b|^2,|x_b-e|^2,|x_a-x_b|^2$: two complex variables $z_a,z_b$ plus the dihedral angle between their half-planes.
- For $\beta_{c\infty}=0$ it reduces to one cross-ratio, which is the T1 case above.
- This is BS's "$\bar K_6$ not constructible" [src BS §7].

Not pursued further.

## 8. Timings (machine load 12–31)

| item | time |
|---|---|
| MZV w ≤ 11 | 13 s |
| SV setup at SV_WMAX = 9 | ~10 s |
| $x_1'$ at weight 11 | >90 s, stopped |
| d=4 periods | <0.1 s |
| q=3 five loops | 0.1–11 s/orbit (Z5: 44 s for 74) |
| q=4 | 1–10 s/orbit (K5−e: 70 s for 106) |
| prism q=6 | 13–32 s |
| prism q=7 | 37–66 s |
| Z5 q=7 | 130–210 s |
| K5 q=8 | 35–143 s |
| six loops q=3 | 259 orbits × 2 plans ≈ 20 min |

Memory < 200 MB.

## 9. Files
- Code: `mzv.py`, `rat.py`, `sv.py`, `expand.py`, `gf.py`, `seq.py`, `validate.py`, `table.py`, `table6.py`, `coverage.py`, `coverage6.py`, `check_plans.py`, `check_highq.py`, `k5_elim.py`, `single.py`, `test_sv.py`, `check_basic.py`, `check_identities_{1,2}.py`, `check_kelvin.py`
- Data: `data/mzv_w*.json`, `tables/*.json` (plus checkpoint `*.partial.jsonl` files)
- Logs: `*.log`
- API: `seq.period(edges, alphas, q) -> (mzv dict, plan)`

## 10. Not done / failed
- Products of two non-rational graphical functions, and relabelings $1\leftrightarrow z$: they leave $V_d$. Products need GSVHs (letter $\bar z$) even at d=4 [src BS Thm 45]. Not implemented [open].
- The direct $(\partial\bar\partial)^{-\lambda}$ form is obstructed at the second append for λ ≥ 2 [drv]. BS's $\mathcal I_n$ algorithm (their Section 25) was not accessible; $V_d$ is an independent derivation [drv].
- Non-sequential orbits: K33, K33+e, Möbius, cube etc. (§5); K33 at q=6 samples have no plan [chk `logs_highq.log`]. Generic K5 needs the two-variable function of §7 [open].
- Z5 q=4: only 16 orbits recomputed (the fiveloop table has all 1956).
- Half-integer q not implemented [open].
- SV_WMAX ≤ 10 used (no overflow occurred; overflow raises an error). MZV basis: the double shuffle + Hoffman reduction is exact, with numerical spot checks only.
