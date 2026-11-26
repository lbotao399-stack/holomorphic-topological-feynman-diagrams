# Circuit cores, census of admissible graphs, atlas of diagrams by theory

Tags: [chk] exact check run (script, log in this folder); [drv] derived, not independently checked; [src] source with location; [open].
P = project root. Paths without prefix are in `explore/circuit-general/atlas/`. Engine `code/htk` (brute, spred, catalog, linalg, conventions) imported, not edited.

| file | content | log |
|---|---|---|
| `graphs_enum.py` | admissible graphs, skeleton enumeration | — |
| `reduce_search.py` | moves as vertex eliminations; searches A, B, R$_k$, I; cores; isomorphism classes | — |
| `polyreduce.py` | exact move maps (times, sources, currents), forward AD, polynomial densities, Dirichlet integration, partial reductions | — |
| `test_polyreduce.py` | full reductions vs brute and known kernels | `test_polyreduce.log` |
| `check_factorization.py` | factorization identities $T=0,1$ | `check_factorization.log` |
| `check_yd_counts.py`, `check_llc_forests.py`, `check_closure.py`, `check_kirr.py`, `check_tight_composites.py` | cross-checks | `.log`, same name |
| `census.py` | census | `census.log`, `census.json` |
| `atlas.py` | theory graph sets | `atlas.log`, `atlas.json` |

## 1. Moves and cores

### 1.1 Moves as vertex eliminations

Network $(V,I)$, $D=H+T$, sources $\lambda_v$, $\sum_v\lambda_v=0$, $S=\{v:\lambda_v\not\equiv0\}$.

| move | condition | new coordinates | sources | reverse currents | transfer |
|---|---|---|---|---|---|
| S | $\deg v=2$, edges $vu,vw$ | $(a_1,a_2)\in\Delta^1$: $t_{vu}=a_1R,\ t_{vw}=a_2R$ | $\lambda_u\mathrel{+}=a_2\lambda_v,\ \lambda_w\mathrel{+}=a_1\lambda_v$ | $y_{v\to u}=a_2\lambda_v-Y,\ y_{v\to w}=a_1\lambda_v+Y$ | polynomial |
| P | $k$ edges $p$–$q$ | $w\in\Delta^{k-1}$: $t_i=T/w_i$ | — | $y_i=w_iY$ | polynomial |
| Y0 | $\deg v=3$, $\lambda_v\equiv0$ | none (merges are P) | — | $y_{v\to u}=-Y_{uw}-Y_{ux},\ y_{v\to w}=Y_{uw}-Y_{wx},\ y_{v\to x}=Y_{ux}+Y_{wx}$ | polynomial |
| Y1 | $\deg v=3$, $\lambda_v\not\equiv0$ | none | $\lambda_a\mathrel{+}=g_a\lambda_v$, $(g_u,g_w,g_x)=\frac{(R_{wx},R_{ux},R_{uw})}{R_{uw}+R_{ux}+R_{wx}}$ | Y0 row $+\,g_a\lambda_v$ | rational |

$$t_{vu}=\frac{R_{uw}R_{ux}}{P_\Delta},\quad t_{vw}=\frac{R_{uw}R_{wx}}{P_\Delta},\quad t_{vx}=\frac{R_{ux}R_{wx}}{P_\Delta},\quad P_\Delta=R_{uw}+R_{ux}+R_{wx}$$
[src P/sources/LLC_main.tex eq-1-64 l.2380, eq-1-63 l.2405, eq-sum-star-triangle l.999, thm-sourced-star l.5176; chk `test_polyreduce.log`, `check_factorization.log`: reconstructed currents $=$ `brute.currents` at every point]

Elimination graph after eliminating $X\subset V$:
$$G_X=G[V\setminus X]\ \cup\bigcup_{C\in\pi_0(G[X])}K(\partial C),\qquad u\in S_X\iff u\in S\ \ \text{or}\ \ \exists C:\ C\cap S\ne\emptyset,\ u\in\partial C$$
[drv: Kron reduction (Schur complement) independent of order; harmonic measure of $C$ on $\partial C$ positive]
$\Rightarrow$ state $=X$; elimination of $v$ valid iff $\deg_{G_X}v\in\{2,3\}$; class (S, Y0, Y1) fixed by $X$ [drv].
Search: 0-1 BFS over $X$, cost $=\#$Y1 [chk `check_closure.log`: closure search $=$ exhaustive search on every (graph, pattern) pair, $D=2$ $n\le6$, $D=3$ $n\le4$, $D=4,5,6$ $n=3$, 0 differences; chk `check_kirr.log`: independent Kron DFS, same irreducible set at $D=2$ $n=6$].

### 1.2 Tightness

$$(D-1)|I|=D|V|-D-1\ \text{preserved by the chain move}\iff\sum_{i=1}^ka_i=D(k-1)+1$$ [src P/PLAN.md §3.3 E1]
Eliminated connected set $C$, $G_X$ tight, configuration $C\cup\partial C$ admissible (LLC Lemma b-00194 l.2033):
$$(D-1)\binom{|\partial C|}2\le D|\partial C|-D-1\ \Rightarrow\ |\partial C|=2\ (D\ge3),\qquad |\partial C|\in\{2,3\}\ (D=2)$$
[drv; chk `check_tight_composites.log`: all components of $G[X]$ in tight states reached only through stars have $|\partial C|=2$ ($D=3$ $n\le4$, $D=4$ $n=3$); $|\partial C|\in\{2,3\}$ ($D=2$ $n\le5$)]
$\Rightarrow$ $D\ge3$: tight composite moves $=$ two-terminal reductions; chain moves $=$ series–parallel case; stars only inside two-terminal bridge subnetworks ($D=3$ $n=3$: 4 graphs $\to$ box) [chk].

### 1.3 Definitions

| object | definition |
|---|---|
| A-core$(\Gamma,S)$ | stuck state of tight polynomial moves: chain moves; $D=2$ also LLC b-01159 (l.3678), b-01173 (l.3745), source-free centres; minimal $|V|$ |
| B-core$(\Gamma,S)$ | minimal stuck $G_X$ under S, P, Y0 |
| $k(\Gamma,S)$ | minimal number of Y1 in a reduction to one edge |
| I-core | stuck $G_X$ under S, P, Y0, Y1 (minimum degree $\ge4$) |

| class | condition | kernel |
|---|---|---|
| SP | no $K_4$ minor | polynomial integral [src htk spred] |
| A | A-core $=$ edge | polynomial, factorizes through tight kernels |
| B | $k=0$ | polynomial integral |
| R$_k$ | $k\ge1$ | rational transfer |
| I | no reduction | no circuit chart |

$\mathrm{SP}\subset\mathrm B$, $\mathrm A\subset\mathrm B$; Y1 allowed $\Rightarrow$ reducibility independent of $S$ [drv; chk I-sets of patterns gen and all coincide, `census.log`].

### 1.4 Factorization, $T=0$

Tight polynomial partial reduction $\Gamma\to\Gamma_c$, coordinates $\xi$ ($k=\dim\xi$), affine core times $T$, $\ell$ = Kruskal chord currents, $\sigma$ = `Graph.sigma()`, $\epsilon(n,H)=(-1)^{nH+\binom n2\binom H2}$ [src P/code/htk/CONVENTIONS.md §2.3]:
$$D_\Gamma\big(t(\xi,T)\big)\,J_{\rm aff}\,(-1)^k=\sigma_\Gamma^H\epsilon(n,H)\,\sigma_c^H\epsilon(n_c,H)\;J(\xi;\lambda,\ell_c)\;D_{\Gamma_c}\big(T;\lambda_c(\xi)\big),\qquad J=\det\frac{\partial\ell_\Gamma}{\partial(\xi,\ell_c)},\quad J_{\rm aff}=\det\frac{\partial t}{\partial(\xi,T)}$$
[chk `check_factorization.log`: 132/132 points exact, all nonzero; 44 reductions: $D=2$ $n=4$ (8 single, 2 paired, 6 mixed), $D=2$ $n=5$ (6 chain, 6 paired, 6 mixed), $D=3$ $n=4$ (10 chain); cores $n_c=1,2,3,4$: 10, 2, 31, 1]
$J$ polynomial in $(\xi,\lambda,\ell_c)$ [drv: each reverse map polynomial in its coordinates, linear in currents and sources]
$$K_\Gamma(\lambda;z)=\pm\int d\xi\;e^{\sum_v\lambda_v\cdot X_v(\xi;z)}\;J\big(\xi;\lambda,-\partial_{z_c}\big)\;K_{\Gamma_c}\big(\lambda_c(\xi);z_c(\xi)\big),\qquad z_{c,f}(\xi)=\sum_eA_{ef}(\xi)z_e,\quad X_v=-\sum_eB_{ev}z_e$$
($y_e=\sum_fA_{ef}(\xi)y_{c,f}+\sum_vB_{ev}(\xi)\lambda_v$) [drv; $z=0$: chk above]
LLC local factors:
$$J_{\rm single}=r_1r_2\det(y^c_{ux},y^c_{wx}),\qquad J_{\rm paired}=(r_1-r_2)(1-r_3-r_4)\langle y^c_{uA},y^c_{AB}\rangle\langle y^c_{uB},y^c_{AB}\rangle\,g_m$$
($g_m$ = Jacobian of the nested merge coordinates onto $(r_3,r_4)$) [chk `check_factorization.log`: ratio $=1$ at every point; src LLC eq-kernel-star-triangle l.3712, eq-kernel-two-leaves l.3788]

Full reduction to one edge, non-tight intermediates allowed:
$$K_\Gamma(\lambda;0)=o\,\sigma_\Gamma^H\epsilon(n,H)\int_{\prod\Delta}\det\frac{\partial\ell}{\partial\xi}\,d\xi,\qquad o=\mathrm{sgn}\det[t;\partial_\xi t]$$
[chk `test_polyreduce.log`: $D(t(\xi))\,|\det[t;\partial_\xi t]|=$ integrand at every point, Y1 reductions included]
B: $\ell(\xi)$ polynomial $\Rightarrow$ $K\in\mathbb Q[\lambda]$, every current moment rational [drv; chk symbolic densities `test_polyreduce.log`]
Y1: $\ell(\xi)$ rational [chk `test_polyreduce.log`: polynomial pass stops at the sourced centre; pointwise identity holds]

### 1.5 Factorization, $T=1$, chain moves

$k=2$ tight chain move, merged edge $e_*\in\Gamma_c$, $x$ = series fractions, $w$ = parallel fractions:
$$D_\Gamma J_{\rm aff}(-1)^k=s\,(w_1w_2)^{-1/2}\,J_x(x,\ell_c)\,D^{[-\eta_*]}_{\Gamma_c},\qquad \eta_*=\sqrt{T_{e_*}}\,w_{e_*},\qquad J_x=\det\frac{\partial\ell_\Gamma}{\partial(x,\ell_c)}$$
$$D^{[F]}=U^{-H-\frac12}\Big\langle F\prod_r\psi^\tau_r\Big\rangle\wedge\prod_{i,r}\psi^i_r\Big/\iota_E{\rm vol},\qquad \langle w_ew_f\rangle=\tfrac12\big(\delta_{ef}/t_e-(C^{\mathsf T}M^{-1}C)_{ef}\big)$$
[src P/explore/circuit-general/topological/REPORT.md §2.1 boxed formula; chk `check_factorization.log`: Wick density without insertion $=$ `brute.density(T=1)` on 15 graphs; 27/27 points exact (sign and square), 9 reductions ($D=2$ $n=4,6$; $D=3$ $n=4$; cores not SP), $s=\pm1$ constant per reduction]
[open] $s$ as a closed rule; $k\ge3$ and nested merges at $T=1$; star moves at $T=1$ [src P/explore/circuit-general/ht-d2/REPORT.md §1].

### 1.6 Exact values by polynomial reduction

| graph | sources | plan | $K(\lambda;0)$ | reference |
|---|---|---|---|---|
| fan(3), $D=2$ | 3 sourced | S, S, S | $-343/864$ | $=$ spred [chk `test_polyreduce.log`] |
| $\theta(1,3,3)$, $D=3$ | all | S$^4$ | $-7/6$ | $=$ spred [chk] |
| $\mathbb C^2$ $K_4$ subdivision $(01,02,03,14,42,13,23)$, $\lambda_3=0$ | 2 integer sets | S4, Y0 at 3, S0 | $\frac{[AD]}{2304}(-3[AB]^2+6[BC]^2-[BD]^2)$: $125/288$, $469/144$ | [chk; src P/explore/nonSP/README.md] |
| $\mathbb C^3$ $K_4[1,1,1,2,2,3]$a, $\lambda_{0,3,6}=e_{1,2,3}$, root 7 | four marks | S$^4$, Y0 at 1, S | $-1/77760$ | $=$ [src P/explore/nonSP/VERIFY.md §3.2] [chk] |
| other five $\mathbb C^3$ $K_4$ subdivisions, four marks | | $k=1$ | $\mathbb Q+\mathbb Q\zeta(3)$ | [chk $k$; src nonSP/VERIFY.md §3.2] |

Four legs at $D=3$: $K(\lambda;0)=c_\Gamma\det(\lambda_1,\lambda_2,\lambda_3)^n$ [drv: $K(g\lambda)=\det(g)^nK(\lambda)$, src P/explore/anomaly-structure/VERIFY.md l.15–36; three independent sources].

## 2. Census

Patterns: gen = bivalent vertices sourced; all = every vertex sourced; ins = gen $\cup\{v\}$, $\deg v\ge3$.
Counts $=$ P/results/graph_counts.txt ($D=2$ $n\le6$; $D=3$ $n\le3$; $D=4$ $n\le2$; $D=5,6$ $n=1$), $=$ topological/REPORT.md §3 ($D=3$ $n=4$: 230; $D=4$ $n=3$: 42) [chk `graphs_enum.py`, `census.log`].
LLC (a), (b), no source condition: $3,11,57,412$ ($n=3..6$) $=$ YD column; with (c): $3,11,60,454$; gen with sourced single star (LLC App. E): $3,11,60,447$ [chk `check_yd_counts.log`, `census.log`].

| $D$ | $n$ | graphs | SP | A | B | B∖SP | R$_1$ | R$_2$ | R$_3$ | I | all: $k=0$ | all: $k$ distribution |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | — |
| 2 | 2 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | — |
| 2 | 3 | 3 | 2 | 3 | 3 | 1 | 0 | 0 | 0 | 0 | 2 | 1: 1 |
| 2 | 4 | 13 | 5 | 11 | 13 | 8 | 0 | 0 | 0 | 0 | 5 | 1: 4, 2: 2, 3: 2 |
| 2 | 5 | 70 | 12 | 53 | 62 | 50 | 6 | 2 | 0 | 0 | 12 | 1: 20, 2: 21, 3: 14, 4: 3 |
| 2 | 6 | 608 | 39 | 339 | 426 | 387 | 84 | 62 | 4 | 32 | 39 | 1: 95, 2: 158, 3: 174, 4: 92, 5: 18, I: 32 |
| 3 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | — |
| 3 | 2 | 2 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 2 | — |
| 3 | 3 | 16 | 10 | 10 | 11 | 1 | 5 | 0 | 0 | 0 | 10 | 1: 6 |
| 3 | 4 | 230 | 66 | 64 | 77 | 11 | 89 | 59 | 5 | 0 | 66 | 1: 83, 2: 42, 3: 39 |
| 4 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | — |
| 4 | 2 | 3 | 3 | 3 | 3 | 0 | 0 | 0 | 0 | 0 | 3 | — |
| 4 | 3 | 42 | 22 | 19 | 23 | 1 | 19 | 0 | 0 | 0 | 22 | 1: 20 |
| 5 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | — |
| 5 | 2 | 4 | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 | 4 | — |
| 5 | 3 | 116 | 60 | 50 | 61 | 1 | 55 | 0 | 0 | 0 | 60 | 1: 56 |
| 6 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | — |
| 6 | 2 | 5 | 5 | 5 | 5 | 0 | 0 | 0 | 0 | 0 | 5 | — |
| 6 | 3 | 242 | 111 | 78 | 112 | 1 | 130 | 0 | 0 | 0 | 111 | 1: 131 |

[chk `census.log`]. A at $D\ge3$: chain moves only (two-terminal bridge reductions of §1.2 not included).
gen with $<2$ sources ($K\equiv0$): $D=2$: $n=3$: 1, $n=4$: 5, $n=5$: 25, $n=6$: 189; $D\ge3$: 0 [chk].
ins, polynomial for every/some/no choice of $v$: $D=2$ $n=6$: 161/426/182; $D=3$ $n=4$: 66/77/153; $D=4,5,6$ $n=3$: 22/23/19, 60/61/55, 111/112/130 [chk].

Core classes, gen (number of sourced core vertices in parentheses) [chk `census.log`, `check_kirr.log`]:

| $D,n$ | A-cores | B-cores | I-cores |
|---|---|---|---|
| 2,4 | prism (0), $K_{3,3}$ (0) | — | — |
| 2,5 | 15 classes, 17 graphs | $K_4$ (4) ×6, $K_5-e$ (2) ×1, $K_5-e$ (3) ×1 | — |
| 2,6 | 173 classes, 269 graphs: $(|V|,|E|)=(5,7)$ ×44, $(6,9)$ ×41 + prism ×24 + $K_{3,3}$ ×5, $(7,11)$ ×77, graph itself ×78 | $K_4$ (4) ×84, $K_5-e$ ×50, $W_4$ ×12, $K_5$ ×32, 6-vertex ×4 | $K_5$ ×32 |
| 3,3 | the six $K_4$ subdivisions | $K_4$ (4) ×5 | — |
| 3,4 | 127 classes, 166 graphs | $K_4$ (4) ×89, $W_4$ (5) ×42, $W_4$ (4) ×5, $K_5-e$ (5) ×8, $K_5-e$ (4) ×4, $K_5$ (5) ×3, prism (6) ×2 | — |
| 4,3; 5,3; 6,3 | 23; 66; 164 graphs | $K_4$ (4) ×19; ×55; ×130 | — |

- pattern all: B $=$ SP in every row [chk `census.log`].
- pattern gen: B∖SP $\ne\emptyset$ for $D=3$ ($n=3$: 1, $n=4$: 11) and $D=4,5,6$ ($n=3$: 1 each); "reducible $\iff$ SP at generic sources" (P/PLAN.md §10) fails [chk]. $n=3$ member: $K_4$ subdivision with chain lengths $[1,1,1,D-1,D-1,D]$ [chk `atlas.log`].
- A $\subsetneq$ SP at $D\ge3$: $(3,4)$ 64/66, $(4,3)$ 19/22, $(5,3)$ 50/60, $(6,3)$ 78/111 [chk]; cf. graph 11 at $n=4$ [src P/explore/six-d-Q3/README.md l.120].
- $D=2$, $n=6$: 32 graphs without reduction by S, P, Y0, Y1 for any sources; I-core $K_5$; all 32 non-planar; non-planar graphs 99 ($n=6$), 5 ($n=5$), 1 ($n=4$), the others reducible [chk `check_kirr.log`].
- $D\ge3$ in range: $k\le3$ ($D=3$ $n=4$), $k\le1$ ($D=4,5,6$ $n=3$) [chk].

## 3. Atlas

Graph sets [chk `atlas.log`]:
$$\text{bulk3}(D,n):\ \max\deg\le3,\ L=(D-3)n+4\ \text{legs at bivalent vertices},\ K\equiv0\ \text{unless}\ L-1\ge H$$
$$\text{cone}(D,n):\ \Gamma\setminus O\ \text{a forest of }p\text{ trees},\ \deg O=n+p,\ \text{others}\ \deg\le3;\qquad \text{twistor}(n):\ \Gamma\setminus O\ \text{connected},\ b_1\le1,\ \deg\le3$$
($L=3|V|-2|I|$ [drv]; span condition [drv: currents are real combinations of sources]). Sources: legs/outputs at bivalent vertices, and $O$. Counts: unlabeled graphs; rooted orbits for cone, twistor. Wick weights, colour factors, derivative-vertex insertions not included.

| set | $n=1$ | $n=2$ | $n=3$ | $n=4$ |
|---|---|---|---|---|
| bulk3 $D=2$ | 1 | 1 ($L=2$: 0) | 1 ($L=1$: 0) | 2 ($L=0$: 0) |
| bulk3 $D=3$ | 1 | 2 | 9: 3 SP, 1 B, 5 R$_1$ | 50: 3 SP, 2 B, 9 R$_1$, 31 R$_2$, 5 R$_3$ |
| bulk3 $D=4$ | 1 | 3 | 27: 7 SP, 1 B, 19 R$_1$ | — |
| bulk3 $D=5$ | 1 | 4 | 83: 27 SP, 1 B, 55 R$_1$ | — |
| bulk3 $D=6$ | 1 | 5 | 185: 54 SP, 1 B, 130 R$_1$ | — |
| cone $D=2$ (fans) | 1 | 1 | 1 | 1 ($n=5,6$: 1, 1) |
| cone $D=3$ | 1 | 2 | 5 | 15 (all SP; 14 A) |
| cone $D=4$ | 1 | 3 | 12 (all SP; 9 A) | — |
| cone $D=5$ | 1 | 4 | 28 (all SP; 21 A) | — |
| cone $D=6$ | 1 | 5 | 51 (all SP; 30 A) | — |
| twistor tree | 1 | 2 | 5 (SP) | 15 (SP) |
| twistor unicyclic | 0 | 3 (SP) | 18: 4 SP, 2 B, 12 R$_1$ | 95: 10 SP, 7 B, 50 R$_1$, 28 R$_2$ |
| product forests $D=2$ (ordered $A,B$) | 1 | 3 | 8 | 23 ($n=5$: 65, $n=6$: 194), all A |

[chk `atlas.log`, `check_llc_forests.log`]; $=$ [src: cone $D=3$: six-d-Q3/README.md l.79–118; twistor 1/0, 2/3, 5/18, 15/95: P/explore/twistor-bulk/operators/README.md l.5–31; bulk3 $D=3$ 9, 50: P/explore/kodaira-spencer/README.md, P/deepen/o3/higher_orders/four_loop_reduction/README.md; product forests 8 ($n=3$): LLC Table tab-ch7-three-loop-forests, 23 ($n=4$): nonSP/README.md; 1, 3, 8, 23, 65, 194: ht-d2/REPORT.md §2.4].

Methods: CF closed all-loop family; SPI circuit polynomial integral (spred: SP; `polyreduce`: B); CT core table / exact evaluator of a non-polynomial core; WB Ward bootstrap; V vanishing theorem; MC Monte Carlo only.

### 3.1 $(1,0)$ 2d chiral; $(0,1)$ topological QM

| $(H,T)$ | computation | $n$ | graphs | method | status |
|---|---|---|---|---|---|
| $(1,0)$ | all | all | banana ($|V|=2$) | CF (parallel move) | computed [src P/PLAN.md §10] |
| $(0,1)$ | $\star$, bracket | all | banana $N=n+1$ | CF: $K_N=\pi^{n/2}$ ($n$ even), $I^\star_N=2^{-N}$ | computed [src topological/REPORT.md §5] |

### 3.2 $(2,0)$ 4d $\mathcal N=1$ gauge theories; $\mathcal N=4$ boundary

| computation | $n$ | graphs | class | method | status |
|---|---|---|---|---|---|
| $Q_n$, $\deg\bar W\le3$ | all | fan, 1 | SP | CF | computed [src LLC thm-sum-all-kernels (i) l.1018] |
| $Q_n$, $\deg\bar W\ge4$ | all | fans, spiders | SP | CF | computed [src LLC (ii)] |
| $A\star B$, $\deg\bar W\le3$ | 1–6 | 1, 3, 8, 23, 65, 194 | A (all) [chk] | SPI (LLC moves) | computed through 7 loops [src LLC l.382, tab-sum-complexity l.1064] |
| $A\star B$, $\deg\bar W\ge4$ | — | trees with several sources | R (LLC App. E) | rational transfer | [open] |
| bulk | 1 | triangle | SP | polygon $|\det|/2$ | computed [src P/PLAN.md §3.1] |
| bulk | $\ge2$ | bulk3 $D=2$, $L\le2$ | — | V (span) | 0 [drv; chk `atlas.log`] |
| $\mathcal N=4$ boundary $Q_n$ | all | fan | SP | CF | computed [src P/explore/twisted-holography-N4/README.md l.1–19, 404–458] |
| $\mathcal N=4$ three-graviton product | $\le2$ | triangle, two-loop | SP | SPI | computed [src same, l.2033–2103] |

### 3.3 $(2,0)$ Poisson BF on $\mathbb C^2$ (T1)

| computation | $n$ | graphs | method | status |
|---|---|---|---|---|
| $Q_n$ | all | fan, 1 | CF | computed [src P/explore/poisson-bf/README.md l.22–42] |
| anomaly | 1 | directed triangle, 2 orientations | SPI | computed; alternating cochain 0 [src l.1–20] |
| $\{\mathbf S\mathbf S\mathbf S\}$ | 1 | 2 Wick triangles | SPI | computed [src l.44–54] |
| $\{\mathbf S\mathbf S\mathbf S\}$ | $\ge2$ | none (ghost count) | V | 0 [src l.44–54] |

### 3.4 $(2,0)$ arbitrary sources (E6: non-BF bulk, several insertions)

| $n$ | graphs | gen: SP / B∖SP / R / I | all: SP / R / I | method | status |
|---|---|---|---|---|---|
| 3 | 3 | 2 / 1 / 0 / 0 | 2 / 1 / 0 | SPI | computed [src nonSP/README.md] |
| 4 | 13 | 5 / 8 / 0 / 0 | 5 / 8 / 0 | WB, CT | computed, all 13, arbitrary sources, $z=0$; $\mathbb Q(\zeta_3,\zeta_5)$ [src nonSP/README.md, P/explore/nonSP/verify/full_ward/README.md] |
| 5 | 70 | 12 / 50 / 8 / 0 | 12 / 58 / 0 | SPI (B), rational (R) | [open] |
| 6 | 608 | 39 / 387 / 150 / 32 | 39 / 537 / 32 | SPI, rational; I: none | [open]; I-core $K_5$ [chk] |

### 3.5 $(3,0)$ 6d $\mathcal N=(1,1)$ SYM: hCS on $\mathbb C^{3|2}$ (T2)

| computation | $n$ | graphs | class | method | status |
|---|---|---|---|---|---|
| $Q_n$ kernels | 1, 2, 3, 4 | 1, 2, 5, 15 (cone $D=3$) | SP | SPI | computed [src six-d-Q3/README.md l.100–124]; weights $\kappa$, $n\le3$ [src l.3–18] |
| $Q_n$ kernels | 5 | 49 | SP | SPI | computed [src l.122] |
| $Q_n$, ladder channel | all | box ladder | SP | CF | computed [src l.126–144] |
| signed $Q_4$, $Q_5$ | 4, 5 | 15, 49 | — | assembly | [open] [src l.116, 122] |
| bulk | $\ge1$ | — | — | V (odd $\delta$) | 0 [src l.5–7] |

### 3.6 $(3,0)$ hCS on $\mathbb C^3$, Kodaira–Spencer on $\mathbb C^3$ (T4): bulk, 4 legs

| $n$ | bulk3 $D=3$ | method | hCS $\mathbb C^3$ status | KS status |
|---|---|---|---|---|
| 1 | box | CF | $K=-\Delta/6$ [src anomaly-structure/VERIFY.md l.66–101] | fixed labels $-1/120$ [src kodaira-spencer/README.md l.38–51] |
| 2 | $\theta(1,3,3)$, $\theta(2,2,3)$ | CF | $K\equiv-\Delta^2/432$, $0$ mod conservation [src AS l.48–90] | $1/3024000$, $1/2721600$ [src KS l.38–51] |
| 3 | 9: 3 SP, $K_4[1,1,1,2,2,3]$a (B), 5 $K_4$ (R$_1$) | SPI, CT (`nonSP/verify/k4_periods`) | all 9 kernels ($c_\Gamma$, nonSP/VERIFY.md §3.2; §1.6); obstruction [open] [src AS l.199–210] | $R=H^\bullet(K3),H^\bullet(T^4)$: 0, $n\ge2$ [src KS l.3–7]; $R=\mathbb C$ polarized: 0 pointwise, all $n$ [src KS l.13–36]; fixed labels [open] |
| 4 | 50: 3 SP, 2 B, 9 R$_1$, 31 R$_2$, 5 R$_3$ | SPI (5), rational (45) | [open] | polarized 0 (50 replays) [src KS l.13–36]; one-descendant $D_4\in\mathbb Q+\mathbb Q\zeta_3+\mathbb Q\zeta_5$ [src KS l.53–68]; individual kernels [open] |

Classes [chk `atlas.log`]. hCS + type I KS on twistor space (T4.5): [open] [src P/PLAN.md §4 T4].

### 3.7 $(3,0)$ twistor BF and Poisson BF (T3)

| computation | $n$ | graphs | class | method | status |
|---|---|---|---|---|---|
| bulk | 1 | box | SP | CF | computed [src P/explore/twistor-bulk/README.md] |
| bulk | $\ge2$ | — | — | V (BF) | 0 [src P/PLAN.md §4 T3] |
| $Q_n$ | 1 | 1 tree | SP | SPI | computed [src twistor-bulk/operators/README.md l.5–31] |
| $Q_n$ | 2 | 2 trees + 3 unicyclic | SP | SPI | computed [src l.54–76] |
| $Q_n$ | 3 | 5 trees + 18 unicyclic | 9 SP, 2 B, 12 R$_1$ | SPI, rational | [open] [src l.22–27] |
| $Q_n$ | 4 | 15 trees + 95 unicyclic | 25 SP, 7 B, 50 R$_1$, 28 R$_2$ | SPI, rational | [open] |
| integrated bilinears | all | — | — | V | 0 [src l.31–52] |
| line defect $f^{(m)}$ | all | spine fan | SP | CF | computed [src P/explore/line-defects/higher_f/README.md] |

### 3.8 $(4,0)$ 8d SYM; $(5,0)$ 10d SYM, BCOV, type I (T5)

| theory | computation | $n$ | graphs | class | method | status |
|---|---|---|---|---|---|---|
| 8d ($\mathbb C^{4|1}$) | $Q_n$ | 1, 2 | 1, 3 (cone $D=4$) | SP | CF, SPI | computed [src P/explore/msym-D8/VERIFY.md l.33–130] |
| 8d | $Q_n$ | 3 | 12 | SP | SPI | computed, 12/12 at random integer sources, $\le5$ s each [chk `e1_cone_kernels.log`] |
| 8d | bulk | $\ge1$ | — | — | V (odd $\delta$, one odd direction) | 0 [drv from six-d-Q3/README.md l.5–7] |
| $D=4$ kernels (engine) | bulk3 | 1, 2, 3 | 1, 3, 27 (7 SP, 1 B, 19 R$_1$) | | CF; SPI, rational | $n\le2$ computed [src P/PLAN.md §3.1; P/explore/twoloop-universal/VERIFY.md]; $n=3$ selected ($\zeta_3$) [src P/PLAN.md §10], rest [open] |
| 10d / type I | bulk | 1 | hexagon | SP | CF | computed [src P/PLAN.md §3.1; P/explore/type-i/README.md l.5–21] |
| type I | bulk, string loop 2 | 2 | 4 thetas | SP | CF | computed [src type-i/README.md l.36–49] |
| BCOV | bulk, genus 2, 3, 4 | 2, 3, 4 | 4; 83 (27 SP, 1 B, 55 R$_1$); — | | SPI, CT, WB, IBP | potentials computed; genus 3 $\zeta_3$ (tetrahedron core), genus 4 in $\{1,\zeta_3,\zeta_5\}$ [src P/explore/iib-bcov/README.md l.964–1106, 1584–1708]; general-source kernels [open] |
| BCOV | $Q_n\varphi(0)$, $\varphi^2/2$, $\varphi^{n+1}$ | all | marked graphs | | WS, CF, SPI | computed [src iib-bcov/README.md l.27–29, 605–930] |
| BCOV | D3 pole $C_4$ | 4 | one-star, prism chambers | | CT | partial; 5286 + 5929 chambers [open] [src l.2317–2611] |

### 3.9 $(1,1)$ 3d $\mathcal N=2$ HT

| computation | $n$ | graphs | method | status |
|---|---|---|---|---|
| all | odd | — | V ($\rho_{\rm odd}=0$) | 0 [src P/explore/ht-families/README.md] |
| $Q_n$, $\deg W\le3$ | 2 | fan, tadpole$(3,0)$ | CF, SPI | computed; tadpole$(3,0)=\frac\pi3l^2$ [src ht-d2/REPORT.md §3.3] |
| $Q_n$, $\deg W\le3$ | 4 | fan, tadpole$(5,0)$, tadpole$(4,1)$, tadpole$(3,2)$ | CF; V; MC | fan computed [src P/explore/fans-ht/VERIFY.md]; tadpole$(4,1)$, $(3,2)$: 0 [src ht-d2 §3.3]; tadpole$(5,0)$ MC [open] |
| $Q_n$, $\deg W\ge4$ | 4 | spiders | SPI | computed [src ht-d2/REPORT.md §2.3] |
| products | 2, 4, 6 | 3; 23 (10 SP, 13 non-SP); 194 | SPI; MC | $n=2$: 3/3; $n=4$: 10 exact, 13 MC [open]; $n=6$ [open] [src ht-d2/REPORT.md §2.4, §2.7] |

### 3.10 $(2,1)$ 5d CS (T6a), $(H,1)$ MSYM, $(5,1)$ 11d supergravity (T6b), $(H,T\ge2)$

| theory | computation | $n$ | graphs | method | status |
|---|---|---|---|---|---|
| 5d CS | bulk QME | 2 | 2 thetas | CF | computed [src P/notes/flagship/flagship.tex l.243–279] |
| 5d CS | bulk | 3 | 9 | V ($\rho_3=0$) | 0 [src flagship l.849–1012] |
| 5d CS | bulk | 4 | 50: 3 SP, 47 non-SP | SPI, CT (cut, IBP, FORCER) | computed, all 50 [src P/deepen/o3/higher_orders/four_loop_reduction/README.md] |
| 5d CS | marked, $\hbar^3$ | 2 | 12 ($J_1$ placements) | SPI | computed [src flagship l.849–1012] |
| 5d CS | $\hbar^{\ge5}$ | $\ge5$ | — | — | [open] [src P/deepen/o3/higher_orders/README.md l.3–46] |
| 5d CS | defects (Oh–Zhou §5.3), $f_2(\kappa)$ | — | — | — | [open] [src P/PLAN.md §4 T6a] |
| MSYM $(1,1),(2,1),(3,1)$ | $Q_2$ | 2 | 1, 2, 3 thetas (cone $D=2,3,4$) | CF | computed [src P/explore/msym-odd/VERIFY.md l.72–183] |
| MSYM $(2,1)$ | $Q_4$ | 4 | 15 (cone $D=3$) | SPI ($T=1$) | computed, 15/15 at random integer sources, values in $\mathbb Q\,\pi^2$ (4 zero), $\le11$ s each [chk `e1_cone_kernels.log`] |
| $(4,1)$ 9d | $Q_2$ | 2 | 4 thetas | CF | computed [src P/PLAN.md §10] |
| 11d sugra | bulk | 2 | 5 thetas | CF | computed $c_{pqr}$ [src flagship l.268; P/explore/sugra-11d-qme/VERIFY.md l.91–114]; class $[X_1]$ [open] |
| 11d sugra | bulk | 3 | 185 (bulk3 $D=6$) | V ($\rho_3=0$) | 0 [src P/PLAN.md §10] |
| $(H,T\ge2)$ | bulk | all | — | V ($\rho\wedge\rho=0$) | 0 [src topological/REPORT.md §4] |
| line defects, $T\ge1$ | — | — | — | — | [open] [src P/PLAN.md §10] |

### 3.11 Summary: theory × computation × $n$

$N$ = diagrams (topologies), $C$ = computed exactly (value or vanishing per diagram).

| theory | computation | $n$: $N$/$C$ | open |
|---|---|---|---|
| $(1,0)$ | all | all: 1/1 | — |
| $(0,1)$ | $\star$ | all: 1/1 | — |
| $(2,0)$ 4d $\mathcal N=1$ | $Q_n$ | 1–6: 1/1 each | — |
| | $A\star B$ | 1: 1/1, 2: 3/3, 3: 8/8, 4: 23/23, 5: 65/65, 6: 194/194 | multi-source forests ($\deg\bar W\ge4$) |
| | bulk | 1: 1/1; $\ge2$: 0 diagrams | — |
| $(2,0)$ Poisson BF | $Q_n$; anomaly; $\{SSS\}$ | all: 1/1; $n=1$: 2/2; $n=1$: 2/2 | — |
| $(2,0)$ E6 | kernels | 3: 3/3, 4: 13/13, 5: 70/58, 6: 608/0 (all K5-core periods $q\le8$ [k5]) | $n=5$: 12 graphs with cores $Z_5$ (6), K33+e (3), V7_0, V7_2, V7_3 (need $\sim10^5$–$10^6$ distinct $q=7$ periods each); $n=6$: 39 SP, 387 B, 150 R, 32 I ($K_5$) |
| $(3,0)$ 6d (1,1) | $Q_n$ kernels | 1: 1/1, 2: 2/2, 3: 5/5, 4: 15/15, 5: 49/49 | — |
| | $(Q^2)_{\hbar^n}=0$ ($\mathfrak{so}(5)$) | 2, 3, 4: all compositions (4: 44); $\kappa_\Gamma=\sigma_\Gamma/|\mathrm{Aut}|$ fixed at $n=4$ (rank 15/15) [physics/6d-nilpotency] | $n=5$ |
| | bulk | 0 diagrams | — |
| $(3,0)$ hCS $\mathbb C^3$ / KS | bulk kernels | 1: 1/1, 2: 2/2, 3: 9/9, 4: 50/50 (with first and second current insertions) | KS fixed labels $n\ge3$; $n=5$: 407 non-SP |
| | obstruction $\omega_n$ | 1, 3: exact; 4: $\equiv0$ (50 graphs) [physics/anomalies] | $Q_2O_1$; KS exchange |
| $(3,0)$ twistor BF | bulk | 1: 1/1 | — |
| | $Q_n$ | 1: 1/1, 2: 5/5, 3: 23/23, 4: 110/110 (with first and second insertions); 5: 555/77 | $n=5$: 470 non-SP, 8 SP over time limit |
| $(4,0)$ 8d | $Q_n$ | 1: 1/1, 2: 3/3, 3: 12/12 | signed $Q_3$ assembly |
| $(4,0)$ engine | bulk3 kernels | 1: 1/1, 2: 3/3, 3: 27/selected | $n=3$: 7 SP, 1 B, 19 R$_1$ |
| $(5,0)$ 10d / BCOV / type I | bulk | 1: 1/1, 2: 4/4, 3: 83/83 (general sources, first and second insertions); cone$(5,3)$ 28/28 | BCOV genus-4 global trace; D3 pole chambers |
| $(1,1)$ 3d $\mathcal N=2$ | $Q_n$ ($\deg W\le3$) | 2: 2/2, 4: 4/4 | — |
| | $\{Q_0,Q_2\}=0$; $(Q^2)_{\hbar^4}$ at $z=0$ | 8 theory families; tadpole$(5,0)$ fixed [physics/3dN2-supercharge] | $\epsilon_2$; derivative letters at $n=4$ |
| | products | 2: 3/3, 4: 23/23, 6: 194/0 | $n=6$ |
| $(2,1)$ 5d CS | bulk | 2: 2/2, 3: 9/9 (0), 4: 50/50 | $\hbar^{\ge5}$; defects |
| $(H,1)$ MSYM | $Q_n$ | 2: 1, 2, 3, 4 / all; 4: 15/15 ($D=3$) | signed $Q_4$ assembly |
| $(5,1)$ 11d | bulk | 2: 5/5, 3: 185/185 (0) | $n\ge4$ |
| | order-$g$ obstruction | 2: $[X_1]\ne0$, $-\tfrac{63g}{256\pi^{10}}$ on the three-plus-six cycle [sugra-11d-qme/collaboration/current_complex] | one-loop counterterms beyond the detector |
| $(H,T\ge2)$ | bulk | all: 0 | — |
| defects | $T=0$ line | all $m$: 1/1 | $T\ge1$ defects |

## 4. Not done; open

- Exact B-class kernels with $n\ge4$ or $D\ge4$ (12×12 polynomial determinant in 12 variables; run stopped at 1.4 GB) [open]; route: $\det\partial\ell/\partial\xi=J_A\cdot\det\partial\ell_c/\partial\xi_c$, $J_A$ the $3\times3$ local chain-move factor, linear in $y_{e_*}$ [drv].
- $T=1$: polynomiality of Y0 moves (contact term, ht-d2/REPORT.md §1.5–1.6) [open]; sign $s$ of §1.5 [open].
- R$_k$ classes: exact evaluators only for three loops (`nonSP/verify/k4_periods`), $\mathbb C^2$ four loops (Ward), 5d CS four loops (FORCER) [src]; others [open].
- I-core $K_5$ ($D=2$, $n=6$, 32 graphs): no S/P/star–triangle chart; $\Delta$–Y moves (rational current map) not implemented [open].
- Not enumerated: $D=4$ $n=4$; $D=5,6$ $n\ge4$. Theory weights (Wick, colour, derivative vertices) not combined with the counts [open].
- $z\ne0$ factorization only [drv]. P/explore/circuit-general/cores/ empty at writing.

## 5. Addendum (orchestrator run, 2026-11-23)

`e1_cone_kernels.py` → `e1_cone_kernels.log`: `htk.kernel_exact` (series–parallel reduction, self-check against `brute` at one random point), random integer sources at bivalent vertices, root $O$.
- cone$(4,3)$, $(H,T)=(4,0)$, 8d $Q_3$: 12/12 exact rationals, $\le4.3$ s each [chk].
- cone$(3,4)$, $(H,T)=(2,1)$, 5d MSYM $Q_4$: 15/15 exact, $K\in\mathbb Q\,\pi^2$, four zero, $\le10.5$ s each [chk].

`e2_twistor_kernels.py` → `e2_twistor_kernels.log`: twistor BF $Q_n$, $D=3$, $T=0$, random integer sources in $\mathbb C^3$.
- series–parallel: $n=3$: 9/9, $n=4$: 25/25 exact rationals by `htk.kernel_exact` (self-checked) [chk].
- $n=3$ non-SP (14): `nonSP/verify/k4_periods/kernel.py` imports python-flint, absent from every interpreter on this machine; not evaluated [open].

`e3_open_kernels.py` (cores evaluator `explore/circuit-general/cores/fastker.py`, `cores.py`):
- twistor BF $Q_3$, $D=3$: 23/23 exact, 3.8 s total [chk `e3_twistor3.log`]. Series–parallel 9: fastker $=$ `htk.kernel_exact` exactly. Non-SP 14: B-class $\#10,\#11$ rational; R$_1$-class 12 in $\mathbb Q+\mathbb Q\zeta_3$ with nonzero $\zeta_3$ coefficient.
- $(1,1)$ tadpole$(5,0)$, $T=1$: $K=\frac{\pi^2}{630}\lambda^4$ ($\lambda=1,2$; $79$ s each) [chk `e3_tadpole50.log`]; ht-d2 MC $0.00159(1)\pi^2$.
- bulk3$(3,4)$, hCS / KS on $\mathbb C^3$: 50/50 exact, 0 failures, $\le203$ s each; 15 zero, 2 rational nonzero (= the 2 B-class graphs), 33 with $\zeta_3$, 28 with $\zeta_5$; SP 3: fastker $=$ `htk.kernel_exact` [chk `e3_bulk3_s0.log`, `e3_bulk3_s1.log`].
- twistor BF $Q_4$ non-SP: 85/85 exact, 0 failures, $\le170$ s each; 4 zero, 7 rational nonzero (= the 7 B-class graphs), 74 with $\zeta_3$, 27 with $\zeta_5$ [chk `e3_twistor4_s0.log`, `e3_twistor4_s1.log`].
- All values in $\mathbb Q+\mathbb Q\zeta_3+\mathbb Q\zeta_5$ [chk].

Five-loop rows (`../fiveloop/REPORT.md` §4) [chk]:
- $(2,0)$, $n=5$, generic sources: 58/70 exact; values in $\mathbb Q$, $+\zeta_3$, $+\zeta_5$, $+\zeta_7$ (W5 cores), $+\zeta_3^2$ (2-sums); B-class patterns rational.
- $(3,0)$ bulk3$(3,5)$: 410 graphs; SP 3/3 exact; non-SP 407 open (whole-graph numerator $4.9\cdot10^9$ points; V8_1–V8_4 periods at $q=8$).
- $(3,0)$ twistor $Q_5$: 555 graphs (49 trees + 506 unicyclic); SP 77/85 exact (8 over 300 s); non-SP 470 open (periods available at $q=8$; numerator assembly open).

Insertions (`../assembly/REPORT.md`) [chk]: bulk3$(3,4)$ 50, twistor $Q_4$ 110, twistor $Q_3$ 23, bulk3$(5,3)$ 83, cone$(5,3)$ 28 graphs: bare kernel, all first and all second current moments exact; prism/K33 tables at $q=8,9$.
