# Holomorphic twist program: plan

2026-11-18

Status tags on every statement:

- **[LLC]** Li–Liu–Chang, `main.tex`, version 20261118-140936.
- **[6d]** `Desktop/6d-holomorphic-twist-note` (hCS on $\mathbb C^{3|2}$ through two loops).
- **[chk]** exact check run on 2026-11-18 (`code/`, `Desktop/check_cd`, `results/`).
- **[lit]** literature, read on 2026-11-18 through a web extraction; formulas not yet compared with TeX source.
- **[drv]** derived here by counting; no independent check yet.
- **[open]** not established.

## 1. Admission rule

A problem is admitted when

1. the theory lives on $\mathbb C^d$ (or on $\mathbb R^T\times\mathbb C^H$ once E8 is done);
2. the object is $\sum_\Gamma \mathcal I_\Gamma(\lambda;z)\,W_\Gamma/|\mathrm{Aut}\,\Gamma|$ with finitely many graphs at each order;
3. the answer is an operator or a number with coefficients in $\mathbb Q$ times a fixed power of $\pi$, or the number field is identified;
4. there is a check in exact arithmetic that does not use the derivation;
5. the action and its normalization are imported from a TeX source, not from memory.

## 2. Pipeline

Same steps as [LLC §3–§5] and [6d §2–§8]; each step has its own check.

| step | content | check |
|---|---|---|
| P1 | nilpotent $Q$; $S=Q\Psi+S_{\rm hol}+{}$total derivative; BV; contractible pairs; Dolbeault superfield action | explicit matrices (`twist10d_check.py` pattern) |
| P2 | letters, $Q_0$ | $Q_0^2=0$ |
| P3 | $[P]_\Omega=\int_\varepsilon^L dt\,\bar\partial^*[K_t]_\Omega$; $Q_{\rm BV}e^{\hbar\partial_P}(\mathcal O e^{I/\hbar})\to e^{\hbar\partial_P}((\mathbf Q\mathcal O)e^{I/\hbar})$, $\mathbf Q=\sum_n\hbar^nQ_n$ | $\int\Omega\wedge[K_t]_\Omega=1$ |
| P4 | $y_e=(\bar x_{e(0)}-\bar x_{e(1)})/2t_e$; $\sum_{e(0)=v}y_e-\sum_{e(1)=v}y_e=\lambda_v$, $\sum_{e\in C}\pm t_ey_e=0$; $\mathcal I_\Gamma=\varepsilon_T\det B_T\int_{\mathbb{RP}^{m-1}_{>0}}e^{-\sum_ey_e\cdot z_e}\bigwedge_{e\notin T}d^dy_e$ | normalization on the one-loop polygon |
| P5 | counting; selection rules of the theory; list of graphs | `enumerate_graphs.py` with the rule |
| P6 | kernels by network reduction: polynomial integrals over simplices | $\det[\partial y/\partial t]$ at random rational $t$ against the reduced density |
| P7 | $Q_n=c_d^{\,n}\sum_\Gamma\mathcal I_\Gamma(\partial;z)W_\Gamma/|\mathrm{Aut}\,\Gamma|$ | — |
| P8 | $\sum_{k+l=n}Q_kQ_l=0$ on explicit operators; literature value | exact linear algebra |
| P9 | list of assumptions | — |

## 3. Engine

### 3.1 Established

Counting on $\mathbb C^d$ [LLC §4.1; 6d (5.11), (5.12); Wang 2401.08113, Def. 16 (4.3)–(4.4) [lit]]:

$$(d-1)|I|=d|V|-d-1\ \Longleftrightarrow\ |I|=dn+1,\ |V|=(d-1)n+2;\qquad (d-1)|I_S|\le d|V_S|-d-1 .$$

$$\text{cycle length}\ \ge d+1,\qquad \sum_v(\deg v-2)=2n-2 .$$

One-loop $(d+1)$-gon, $1>\sigma_1>\dots>\sigma_d>0$ [chk, $d=2,\dots,5$]:

$$y_1=-\sum_{j=1}^d\sigma_j\lambda_j,\qquad \text{density}=(-1)^d\det(\lambda_1,\dots,\lambda_d),\qquad |\mathcal I(\lambda;0)|=\frac{|\det(\lambda_1,\dots,\lambda_d)|}{d!} .$$

Chain move: $d-1$ vertices of degree two on $u\to v_1\to\dots\to v_{d-1}\to w$, parallel to an edge $u\to w$ [chk, pointwise with sign, $d=2,\dots,5$]:

$$\sigma_j=\frac1T\sum_{k>j}t_k,\quad\rho=\frac{T}{T+t_0},\quad y_0=\rho\,y^{\rm red},\quad y_1=(1-\rho)y^{\rm red}-\sum_j\sigma_j\lambda_j,$$

$$\text{local density}=\rho^{\,d-1}\det\big(y^{\rm red},\lambda_1,\dots,\lambda_{d-1}\big),\qquad \lambda_u^{\rm red}=\lambda_u+\sum_j\sigma_j\lambda_j,\quad\lambda_w^{\rm red}=\lambda_w+\sum_j(1-\sigma_j)\lambda_j .$$

Bubble move on $\mathbb C^3$ [6d (7.4)–(7.6)].

Closed forms:

| kernel | value | status |
|---|---|---|
| fan, $\mathbb C^2$, $n$ loops, no input derivatives | $c^{(n)}=1/[(n!)^3(n+1)]$ | [LLC] |
| box, $\mathbb C^3$ | $\frac16\det(\lambda_1,\lambda_2,\lambda_3)$ | [6d][chk] |
| double box (A), spider (B), $\mathbb C^3$ | [6d (1.5)–(1.9)] | [6d][chk] |
| $n$-box ladder, $\mathbb C^3$, sources at the four corners, $\Delta=\det(A,B,C)$ | $\mathcal I/\Delta^n=\dfrac{(-1)^{n-1}2^{n+1}}{(n+2)!\,n!^2\,(n+1)!^2}$ | [All-order density, rung moments and eight-loop Kirchhoff replay](explore/six-d-Q3/README.md#four-corner-box-ladder) |

Operators: $Q_0^2=0$, $\{Q_0,Q_1\}=0$, $Q_1^2+\{Q_0,Q_2\}=0$ for hCS on $\mathbb C^{3|2}$ [6d §8].

Files:

| path | content |
|---|---|
| `Desktop/check_cd/check_cd.py`, `check_ladder.py` | polygon, chain move, ladder |
| `Desktop/6d-holomorphic-twist-note/code/kirchhoff.py`, `check_twoloop_kernels.py`, `check_z_dependence.py` | direct Kirchhoff solver, two-loop kernels |
| `Desktop/6d-holomorphic-twist-note/code/salg.py`, `ops.py`, `ops2.py`, `test_q0q1.py`, `test_q2.py` | superfield operators, nilpotency |
| `Desktop/6d-holomorphic-twist-note/code/twist10d.py`, `twist10d_check.py` | twist with explicit $\Gamma$ matrices |
| `code/enumerate_graphs.py`, `results/graph_counts.txt` | graphs, reach of the moves |

### 3.2 Reach of the moves [chk]

`results/graph_counts.txt`. SP = series–parallel (no $K_4$ minor). YD = reducible to an edge by [LLC] moves (a), (b) inside Laman graphs. The $d=2$ totals reproduce the Laman numbers $1,1,3,13,70,608$.

| $d$ | $n$ | $|V|$ | $|I|$ | graphs | SP | not SP | YD ($d=2$) |
|---|---|---|---|---|---|---|---|
| 2 | 1 | 3 | 3 | 1 | 1 | 0 | 1 |
| 2 | 2 | 4 | 5 | 1 | 1 | 0 | 1 |
| 2 | 3 | 5 | 7 | 3 | 2 | 1 | 3 |
| 2 | 4 | 6 | 9 | 13 | 5 | 8 | 11 |
| 2 | 5 | 7 | 11 | 70 | 12 | 58 | 57 |
| 2 | 6 | 8 | 13 | 608 | 39 | 569 | 412 |
| 3 | 1 | 4 | 4 | 1 | 1 | 0 | |
| 3 | 2 | 6 | 7 | 2 | 2 | 0 | |
| 3 | 3 | 8 | 10 | 16 | 10 | 6 | |
| 4 | 1 | 5 | 5 | 1 | 1 | 0 | |
| 4 | 2 | 8 | 9 | 3 | 3 | 0 | |

- Every SP graph of the table reduces to one edge by the parallel-chain moves of E1 alone (column CH of the output equals SP in all rows).
- $d=2$, $n=4$: the two graphs not reduced are the prism and $K_{3,3}$.
- $d=3$, $n=3$: the six graphs not SP are subdivisions of $K_4$ with chain lengths $[1,1,1,1,3,3]$, $[1,1,1,2,2,3]$ (three), $[1,1,2,2,2,2]$ (two); all have maximal degree three.
- $n=2$, any $d$: theta graphs with path lengths $a\le b\le c\le d$, $a+b+c=2d+1$:

| $d$ | $(a,b,c)$ |
|---|---|
| 2 | $(1,2,2)$ |
| 3 | $(1,3,3),(2,2,3)$ |
| 4 | $(1,4,4),(2,3,4),(3,3,3)$ |
| 5 | $(1,5,5),(2,4,5),(3,3,5),(3,4,4)$ |
| 6 | $(1,6,6),(2,5,6),(3,4,6),(3,5,5),(4,4,5)$ |

- If $\Gamma\setminus\mathcal O$ is a forest, $\Gamma$ has treewidth $\le2$ and is SP. This is the case for $Q_n$ in 4d $\mathcal N=1$ gauge theories [LLC] and in hCS on $\mathbb C^{3|2}$ [6d (5.13)].
- hCS on $\mathbb C^{3|2}$, rule [6d (5.13)] $\deg_\Gamma\mathcal O=n+p$ with $\Gamma\setminus\mathcal O$ a forest of $p$ trees: $n=1$: 1 graph; $n=2$: 2 graphs; $n=3$: 5 of the 16 graphs, all SP, $\deg\mathcal O=4$, $p=1$, skeleton chain lengths $[1,1,2,3,3]$ (two), $[1,2,2,2,3]$ (two), $[2,2,2,2,2]$.

### 3.3 Engine problems

**E1. Series–parallel reduction on $\mathbb C^d$.** $k$ chains of lengths $a_1,\dots,a_k$ between $p$ and $q$, interior vertices of degree two, replaced by one edge:

$$\underbrace{\sum_i(a_i-1)}_{\text{series ratios}}+\underbrace{(k-1)}_{\text{parallel ratios}}=d\,(k-1)\ \Longleftrightarrow\ \sum_ia_i=d(k-1)+1 .$$

$k=2$: chain $(1,d)$ [chk], bubble $(2,2)$ at $d=3$ [6d]. First cases with no pair $a_i+a_j=d+1$: $(3,3,3)$ at $d=4$, $(3,4,4)$ at $d=5$.
Deliverable: the local density for general $(a_1,\dots,a_k)$ with sources on all interior vertices; `sp_reduce(G, d)` returning ratio coordinates and polynomial density.
Check: pointwise against the direct Kirchhoff density for every SP graph of §3.2 ($d=3$: 13 graphs; $d=4$: 4; $d=5$, $n\le2$: 5).

**E2. Displacements.** $-\sum_ey_e\cdot z_e=\sum_v\lambda_v\cdot X_v(r;z)$ with $X_v$ polynomial in the ratios, for every SP graph. Known: box [6d (6.4)], (A), (B) [6d (1.8)]; fans [LLC §2.1, theorem on the all-loop supercharge].
Check: random rational points (`check_z_dependence.py` pattern).

**E3. Derivative vertices and derivative propagators.** $-\partial/\partial z_e$ inserts $y_e$ [LLC §2.4, circuit dictionary]. Kernels with insertions $\prod_ey_e^{m_e}$ are polynomial integrals after E1.
Needed by: Poisson brackets (two $\partial$ per vertex), $\alpha\mu_i\partial_i\gamma$, kinetic operators $\bar\partial\partial^{-1}$.
Check: moments of the direct Kirchhoff solution.

**E4. Operator algebra.** `salg.py`, `ops.py`, `ops2.py` generalized to arbitrary $d$, letters, and cubic vertices with derivatives; relations $\sum_{k+l=n}Q_kQ_l=0$ as linear systems over $\mathbb Q$.

**E5. Selection rules.** For each theory, the rule that cuts the list of §3.2 (odd-coordinate delta functions [6d (5.13)], antifield count [LLC §4.3], orientation of propagators in BF-type theories). Output: graph list per $(\text{theory},n)$.

**E6. Graphs that are not series–parallel.** $d=2$: prism, $K_{3,3}$. $d=3$: the six $K_4$ subdivisions.

1. Is $\mathcal I_\Gamma(\lambda;0)\in\mathbb Q[\lambda]$? The four-loop $K_{3,3}$ admits integer momenta with $K=-789/640+59\zeta(3)/80$. The prism has an integer-source period $-43/20160+\zeta(3)/560$. The three-loop $K_4$ subdivision $[1,1,1,1,3,3]$ admits a rational-source period $(39-32\zeta(3))/67200$. These periods are irrational by Apéry’s theorem. [Six-source $K_{3,3}$ calculation](explore/nonSP/verify/k33_forcer/README.md) · [$K_4$ calculations](explore/nonSP/README.md)
2. For the eight target graphs the independent Kirchhoff currents give constant signed volume density on two explicitly described semialgebraic regions. At three loops a rational linear current insertion has an irrational period, excluding a finite signed rational-polytope decomposition of those regions. [Current coordinates](explore/nonSP/README.md#signed-current-coordinates)
3. At four loops on $\mathbb C^2$ all sixty zero-displacement quadratic identities are verified at arbitrary sources. The edge-contraction map has rank $2344$ and nullity $15$. Local numerator zeros, wheel integration by parts and radial integrals $20\zeta(5)$ and $5(\zeta(3)-\zeta(5))/3$ determine every coefficient of the thirteen kernels. [Exact Ward equations and analytic reconstruction](explore/nonSP/verify/full_ward/README.md)

Tools: sourced star–triangle [LLC App. E]; quadratic relations [operatope 2207.14321 (51); Wang 2401.08113, third main theorem] [lit]; direct integration of the Kirchhoff density.

**Computational update, 26 November 2026 [drv, chk]:** [Exact series fibres and graph-rank support](explore/circuit-general/series-fibres/README.md). At arbitrary original sources, bivalent-chain column differences are constant source columns. Their exact elimination leaves an $(m-1)$-square determinant for $m$ maximal chains, polynomial in the chain fractions with total degree at most $m-1+k$ for $k$ current factors. Exact rational product-simplex cubature integrates those fractions before core reconstruction. Separately, the all-loop polynomial numerator obeys $\deg_S P\le n\operatorname{rank}C_S-1$ for every nonempty edge set; the assembly and five-loop evaluators now use these bounds. Complete coefficient comparisons, the independent complete five-loop fan period, and nonzero five/six-loop fibre certificates pass. Pointwise full-current controls extend through twelve loops. Remaining at that stage: core integration/period tables for the retained resistance variables, and the $T=1$ contact-term extension. The previous 4,855,374,080-point whole-graph route is not a claim that a general five-loop $H=3$ period is now solved.

**Computational update, 27 November 2026 [drv, chk]:** [Certified reciprocal core reconstruction](explore/circuit-general/series-fibres/CORE_RESULTS.md) completes the remaining core integration for three explicitly sourced six-loop $K_5$ graphs, covering every admissible $K_5$ subdivision topology at $H=2$. Core resistance inversion lowers degree 35 to 15; graphic-flat inequalities give a 683,834-point interpolation grid, and exact three-node interval rules require 27 contracted determinants per point. A permanent/Mahler bound certifies the entire integer-scaled coefficient dictionary after seven primes. The three complete kernels lie in $\mathbb Q+\mathbb Q\zeta(3)+\mathbb Q\zeta(5)$, with rational positivity certificates and six independent rational fibre replays. The period bridge exposes all 54,391 existing $K_5$ orbits through level 8. Scope: bare $T=0$ at supplied exact sources, with all sources retained; compiled current insertions, longer-chain rules, additional core period tables and $T=1$ remain extensions.
Benchmarks at $d=2$ from the bootstrap of 2207.14321 §4 [lit]: bitriangle, square-triangle, tritriangle (three loops); square-triangle at leading order $\frac1{24}(\lambda_1\wedge(\lambda_2+2\lambda_3))(\lambda_2\wedge(\lambda_1+2\lambda_3))$, the fraction bar reconstructed.
**Computational update, 27 November 2026 [drv, chk]:** [All-loop cubic-forest current messages](explore/circuit-general/forest-moments/README.md) implement the manuscript's single/paired-star theorem without global ratio charts or precomputed tree tables. Subtrees retain three boundary currents, an exact outside marked-degree argument, and at most four local ratios; component assembly performs the complete insertion-source translation. This gives every finite original-edge current word at arbitrary supplied rational root sources. All 100 topologies through five loops pass 2,292 original-chart comparisons. Complete certificates include a nonzero twenty-loop degree-three moment, two independently integrated nine-loop branched/multi-root examples, and a thirty-six-loop known fork regression. The fan specialization batches a whole finite displacement expansion with at most $\binom{k+3}{3}$ coefficient states through order $k$, independently of loop count: all nine coefficients through degree eight at 128 loops are recorded and replayable. Higher-degree product forests are recognized separately; the series-parallel portion retains all sources through the existing generic determinant integrator. At that stage, the remaining non-series-parallel sourced-star rational messages were open; the following continuation supplies a source-complete core route.

**Computational continuation [drv, chk]:** [Source-complete product forests](explore/circuit-general/sourced-forests/README.md). Passive single/paired stars require zero source only at the removed centres, so the planner retains arbitrary higher-valence parent sources. Every antisymmetric current bracket has one Kirchhoff denominator; exact factor cancellation preserves the original denominator level under passive elimination. The residual core engine combines the existing exact linear current-minor determinant, chain beta integrals, graph symmetries, and read-only finite period tables. Complete certificates include irreducible physical five-loop and nonplanar physical six-loop quartic forests in the rational/zeta(3)/zeta(5) channels, independently replayed through the planar face dual, and rational current moments through nineteen loops independently checked from an all-order boundary formula. When no finite period reduction is available, the public evaluator returns the full explicit convergent rational series with certified rational tails. The support theorem proves separate convergence of every scalar monomial at every finite current order. Exact computability at rational sources is established with explicit resource limits; a universal finite zeta classification and scalable completion of all large cores remain open.

**Computational continuation [drv, chk]:** [Forest periods and quantum Leibniz](explore/circuit-general/forest-leibniz/README.md). The existing cochain identity couples the new product coefficient under $Q_0$ to all lower product coefficients under the higher supercharges. Rational forcing leaves rational homogeneous directions whose coefficients can be actual periods; closedness does not imply exactness. The new three-loop Ward evaluator retains all original currents, including their affine limits on contracted edges, and integrates the proper diamond–triangle boundary independently. All 153 words through degree two and 25 additional complete controls replay exactly. A factorial paired-star transfer gives full-spinor lifts at five, seven, nine, thirteen and nineteen loops, with all 42 quotient moments explicitly evaluated. Five- and seven-loop face censuses are exhaustive; larger examples are inductive lifts of the base all-word identity. One period-valued outer forest reconstructs its partner from the rational middle term and boundary. The new weighted conformal-combination reducer preserves every unresolved period at search bounds; its complete five-loop control resolves 1,007 input orbits, cancels five unresolved components, and retains 5,694 residual terms. Full Wick/coupling assembly for a specified higher superpotential and finite closure of all period cores remain open.

**Nonrenormalizable continuation [drv, chk]:** [Complete quartic products](explore/circuit-general/nonrenormalizable-product/README.md) assemble all Wick, derivative and representation coefficients for specified abelian and SU(2) quartic theories. The entire three-loop product is rational in both examples, with a proof that higher coefficients vanish for these inputs. Keeping a neutral spectator gives a stronger full-source calculation: the nonabelian matter–matter and matter–gauge sums have opposite nonzero zeta(3) polynomials. Independent literal delta-operator tests check every candidate graph; separate color contractions verify every output color entry. These are color-singlet cochains, not asserted quantum cocycles. Universal rationality of all assembled products remains open.

**Further computational reduction [drv, chk]:** [Sourced degree-two ears](explore/circuit-general/forest-general-reduction/README.md) retain the sources at removed vertices and all original current words on arbitrary tight residual cores. Any number of ears sharing an edge is integrated using exact simplex factorials and at most two aggregate source coordinates. Complete Schwinger polynomials are summed before chain integration and period lookup. Physical higher-superpotential fixtures at four, five, six, eight and eleven loops have complete bare and first-current spinor answers; the interaction degree grows from five to twelve. Independent original-graph controls and sourced four/five-loop Ward identities verify the reduction. Individual kernels retain zeta(3), while the full signed Ward sums are rational. Irreducible cores still require the finite period tables, further relations, or the complete convergent-series fallback.

**Whole-numerator period reduction [drv, chk]:** [Targeted family relations](explore/circuit-general/forest-target-reduction/README.md) reconstruct the complete missing five-loop numerator and close only its reached scalar-period family. All 10,577 original coefficients and 7,368 scalar coefficients are checked against immutable input digests. The local exact system has 47,092 variables, 163,907 rows, and 46,879 pivots; its entire target has nineteen remaining scalar-period terms. Adding 687 lower-dimension identities and closing their reached family gives 61,130 variables, 215,623 rows, 60,844 pivots, and seventeen residual terms, with no budget stop. Every remaining power vector and coefficient is recorded with the finite contribution. This improves the previous 5,694-term residual representation without asserting independence of the remaining terms or nonvanishing of their combination. All 894 lower-level local identities and eighteen double-dimension controls agree with independent full tables. The following continuation completes this target.

**Finite five-loop completion [drv, chk]:** [Polynomial annihilators and exact proof weights](explore/circuit-general/forest-finite-target/README.md) evaluate the complete source-preserving kernel as $(1187-816\zeta(3))/435456000$, strictly positive by rational bounds on $\zeta(3)$. Three degree-one and fourteen degree-two annihilators couple the power-sum blocks preserved by the earlier local identities. The new system reaches 14,553 variables, including 4,096 outside the previous system; 2,580 pivots solve the complete target. A certificate of 1,039 rationally weighted original identities includes every coordinate face and directly reproduces all seventeen target coefficients. All residual periods and the zeta(5) coefficient cancel exactly. Default verification reconstructs the complete original numerator and sums the witness, without repeating elimination. Independent tests cover 917 lower-level identities, 310 face/convergence guards, and all seventeen generators. The existing generator cache is transported and verified; the reusable bounded engine accepts other supplied polynomial archives and complete targets.

**Universal product cancellation [drv, chk]:** [All color indices and arbitrary-degree pure matter superpotentials](explore/circuit-general/universal-quartic-product/README.md). The full quartic product of $\beta^{1i_1}(\partial_{+\dot a}\beta^{2i_2})\gamma_5$ with $b_A\beta^{3i_3}$, for a neutral free $\gamma_5$, has twelve formal color words before the representation relation. Its entire zeta(3) coefficient factors into $[t_A,t_B]-iC^C{}_{AB}t_C$ acting on the first or second superpotential slot. No gauge-invariance identity for the superpotential tensor is needed for this cancellation. The remaining ten rational coefficient spinors are complete. A finite classification by whether an additional superpotential leg participates extends cancellation to every pure-matter degree at least three; all extra surviving sources remain present. Lower-loop graphs are series-parallel and fourth and higher orders vanish by the antifield/even-field counts. This proves rational normalized coefficients at every order for the fixed input family, without claiming an explicit evaluation of each lower-order coefficient or rationality for arbitrary inputs and field-strength-dependent superpotentials.

**Next reduction targets [open]:** Extend the verified polynomial-annihilator archives to further sourced cores, and combine the complete numerator with symbolic representation identities before choosing color entries. The [sourced-star chart calculation](explore/circuit-general/forest-general-reduction/SOURCED_STAR_OBSTRUCTION.md) exhibits a nonzero local logarithm despite a rational complete kernel, so a general sourced-star integration step must retain rational functions and their exact primitives. These are distinct remaining computational requirements; increasing the former local-search budget alone does not supply them.

Needed by: theories whose graphs are not cones over forests (bulk graphs of non-BF theories, brackets with several insertions).

**E7. Line defects.** Vertices restricted to a curve $C\subset\mathbb C^3$ (bulk–defect diagrams of Koszul duality) [lit]. Deliverable: current representation with defect vertices. [open]
Benchmarks [lit]: $\mathcal M_{0000}=\frac{1}{96\pi^2}\big(\frac{3z_0}{z_{12}}+\frac12\big)d\bar z_0$, $\alpha=\frac1{32\pi^2}$, $\beta=\frac{h^\vee}{48\pi^2}$ (Fernández 2302.14292 (36), (43)–(44)); $\mathcal I_{0,0,0,0}=-\frac{3}{16\pi^2(w-z)}$, $\alpha=\frac{2}{5\pi^2}$ (Bittleston 2211.06417 (C.15), (3.21)).

**E8. $\mathbb R^T\times\mathbb C^H$.** Counting with $n=H+T$ (Gaiotto–Kulp–Wu 2403.13049 (3.10)–(3.12); Wang–Williams 2407.08667 Def. A.10) [lit]:

$$n|V|=(n-1)|I|+n+1,\qquad n|V_S|\ge(n-1)|I_S|+n+1 .$$

Representation (2403.13049 (3.15), (3.18), (3.52)) [lit], with $y_e=\bar x^{\mathbb C}_e/t_e$, $s_e=x^{\mathbb R}_e/t_e^{1/2}$:

$$I_\Gamma(\lambda;z)=\int_{\Delta_\Gamma}\prod_{v\ne v_*}\delta\Big(\lambda_v-\!\!\sum_{e(0)=v}\!\!y_e+\!\!\sum_{e(1)=v}\!\!y_e\Big)\prod_e\pi^{-T/2}e^{-s_e^2-y_e\cdot z_e}\,d^Hy_e\,d^Ts_e,\qquad\omega_\Gamma=\alpha_\Gamma^{\,T}\prod_{i=1}^H\rho_\Gamma[\lambda^i;z^i] .$$

Vanishing [lit]: $\alpha_\Gamma=0$ at odd loop number; $\alpha_\Gamma\wedge\alpha_\Gamma=0$ for every graph with a loop (Balduf–Gaiotto 2408.03192 Thm 2; 2407.08667 Props. 3.23, 3.26). For $T\ge1$, loop corrections exist only at $T=1$ and even loop number.
Known integrands, not integrated [lit]: $H=T=1$ bitriangle, $\alpha\propto U^{-3/2}$, $\omega\propto U^{-5/2}$ with $U$ the spanning-tree polynomial (2403.13049 (106), (109)). No two-loop graph with $T\ge1$, $H\ge1$ integrated to a closed value was found.
Deliverable: $\alpha_\Gamma$ in the ratio coordinates of the E1 moves; closed values for the two-loop theta graphs at $T=1$: $(1,2,2)$ for $H=1$; $(1,3,3),(2,2,3)$ for $H=2$; the five graphs of the $d=6$ row of §3.2 for $H=5$. [open]

## 4. Targets

| | theory | space | vertices | bulk type | first use of |
|---|---|---|---|---|---|
| T1 | 4d $\mathcal N=1$ gauge theory with the twisted gravity multiplet (Poisson BF) | $\mathbb C^2$ | cubic, two $\partial$ | BF | E3 |
| T2 | hCS on $\mathbb C^{3|2}$ (6d $\mathcal N=(1,1)$ SYM), $Q_3$ and all $n$ | $\mathbb C^3$ | cubic | CS | E1 |
| T3 | twistor BF and Poisson-BF with the field $\eta$ (SDYM, SDGR) | $\mathbb{PT}$, locally $\mathbb C^3$ | cubic, up to two $\partial$ | BF | E3, E7 |
| T4 | Kodaira–Spencer on $\mathbb C^3$ (IIB on $T^4$, K3); hCS with type I KS on twistor space | $\mathbb C^3$ | cubic, $\partial^{-1}$ kinetic term | not BF | E3, E6 |
| T5 | hCS on $\mathbb C^5$ (10d SYM), BCOV on $\mathbb C^5$ (IIB), type I | $\mathbb C^5$ | cubic | not BF | E1 at $d=5$ |
| T6 | 5d Chern–Simons on $\mathbb R\times\mathbb C^2$; twisted 11d supergravity on $\mathbb R\times\mathbb C^5$; Poisson BF on $\mathbb R\times\mathbb C^2$ (5d $\mathcal N=1$) | $\mathbb R\times\mathbb C^H$ | cubic, up to two $\partial$ | not BF; BF for the last | E8 |

Legs of a graph with all vertices cubic and no insertion [drv]: $3|V|-2|I|=(d-3)n+4$. At $d=3$: four legs at every loop order. A nonzero kernel needs sources spanning $\mathbb C^d$.

### T1. Poisson BF on $\mathbb C^2$ with matter

Source [lit]: Elliott–Williams 2008.02302, Conjecture 4.3, eq. (1), §4.

$$S=\int_{\mathbb C^2}B\wedge\big(\bar\partial A+\tfrac12\{A,A\}_\Pi\big)+\int\beta\wedge\bar\partial\gamma+\int\beta\wedge\{A,\gamma\}_\Pi,\qquad\{f,g\}_\Pi=\partial_zf\,\partial_wg-\partial_wf\,\partial_zg .$$

Benchmarks [lit]:

- pure Poisson BF on $\mathbb C^2$ has no one-loop anomaly (2008.02302 Thm 3.4, Rem. 4.5; by Gelfand–Fuchs cohomology, no graph evaluated);
- Bomans–Wu 2311.04304 (3.14), (4.15), (4.24): $\{\mathbf S\,\mathbf S\,\mathbf S\}=\tfrac43(3c-2a)\Lambda_{\mathcal N=4}+16(a-c)\Lambda_{\rm vector}$, with $(a,c)=(\tfrac3{16},\tfrac18)$ for a free vector and $(\tfrac{N_f}{48},\tfrac{N_f}{24})$ for chirals at $\zeta=\tfrac13$.

Problems:

1. The directed triangle integrates to the fourth power of the external momentum determinant times $-1/120$. Its alternating local cochain vanishes. [Calculation](explore/poisson-bf/README.md)
2. With the displayed cubic vertices, only the fan survives. Its ordered Poisson kernel and arbitrary input jets are computed at every loop order. [Calculation](explore/poisson-bf/README.md#local-operator-with-cubic-poisson-vertices)
3. For the improved bilinear stress tensor, the connected central bracket receives only the two Wick triangles. Their coefficients reproduce the superconformal anomalies at fixed R-charge. In $\mathcal N=4$ the bracket equals $(\dim\mathfrak g)\Lambda_{\mathcal N=4}/3$. [Derivation and replay](explore/poisson-bf/README.md#ternary-stress-tensor-bracket)

### T2. hCS on $\mathbb C^{3|2}$ beyond two loops

1. $Q_3$ has weights $(-1/2,1/2,-1,1/8,-1/2)$ in the five-graph basis. The cycle cancellation and BV insertion identity give the universal Ward equation; the exact rank minor is $6\bmod101$. [Derivation and replays](explore/six-d-Q3/README.md)
2. The oriented tree generator and Schwinger reduction give signed graph kernels at every order. The full $Q_n$ is nonzero on conjugation-invariant cochains, even at fixed $\mathfrak{sl}_2$ modulo total holomorphic derivatives, for all $n\ge1$. [Graph generator](explore/six-d-Q3/README.md#rooted-graph-count-at-every-loop-order) · [Exact matrix element](explore/six-d-Q3/README.md#a-nonzero-matrix-element-of-the-full-differential-at-every-loop) · [Fixed-rank proof](explore/six-d-Q3/README.md#nonvanishing-modulo-total-derivatives-at-fixed-rank)
3. Every invariant on reduced commuting $\mathfrak{sl}_N$ scalar pairs lifts by finite marked products. The lift of a scalar of bidegree $(p,q)$ stops by $\hbar^{\min(p,q)}$. On the regular semisimple locus this extends to all finite eigenvalue jets. [Quantum representatives](explore/six-d-Q3/lifting/scalar_products/README.md) · [Regular jets and translation charge](explore/six-d-Q3/lifting/scalar_products/regular_jets.md) At field length $p$ and rank parameter $N\ge p$, every local classical class survives in the absolute matrix current complex. [Stable-range calculation](explore/six-d-Q3/lifting/README.md#stable-rank-at-fixed-field-length) · [Descent and integrated cohomology in the same range](explore/six-d-Q3/lifting/README.md#integrated-cohomology-in-the-same-rank-range)

### T3. Twistor space

Sources [lit]: Bittleston–Sharma–Skinner 2208.12701 (2.5), (2.15)–(2.18), (5.27); Costello 2111.08879.

$$S[b,a]=\frac1{2\pi i}\int_{\mathbb{PT}}b\wedge\big(\bar\partial a+\tfrac12[a,a]\big),\qquad S[g,h]=\frac1{2\pi i}\int_{\mathbb{PT}}g\wedge\big(\bar\partial h+\tfrac12\{h,h\}\big) .$$

Bulk graphs stop at one loop (2208.12701 §4.2). Benchmarks for the box with the heat-kernel regulator of 2208.12701 §4.2–4.3 [lit]:

$$\frac1{3!}\Big(\frac i{2\pi}\Big)^3\int_{\mathbb{PT}}\mathrm{tr}_{\rm ad}\big(c\,(\partial a)^3\big)\ \ (3.25),\qquad\frac1{5!}\Big(\frac i{2\pi}\Big)^3\int\mathrm{tr}\big(\psi\,(\partial s)^3\big)\ \ (4.52),\qquad\text{mixed}\ (3.29) .$$

Problems:

1. The six cyclic orders give $1/3!$, $-1/24$, and $(1+\dim\mathfrak g)/5!$ in the common heat normalization. [Bulk box replay](explore/twistor-bulk/README.md)
2. $Q_n$ on local operators of the BF theories: tree and unicyclic lists, including three two-loop unicyclic operators. On flat space, all positive-loop contributions on integrated bilinears vanish in both BF theories. [Calculation and replays](explore/twistor-bulk/operators/README.md)
3. The ordered defect coefficient $f^{(m)}$ is an ordered interval integral at every loop order and every input and output jet. It obeys the associativity recursion of Fernández–Paquette 2412.17168 (74), (80), with seed $f^{(1)}_{(0,1)(1,0)}(0,0)=1/(32\pi^2)$. [Integral, proof and exact replays](explore/line-defects/higher_f/README.md)
4. The coefficient of $\epsilon^{ijkl}f_{abc}[12]$ in 2604.21015 §5.2.2 vanishes for the completed insertion. [Two-component Ward identity](explore/line-defects/README.md#the-one-loop-three-point-form-factor-in-self-dual-supersymmetric-yangmills-theory). In the commuting-current convention of 2201.02595 (7.3.5), $C=-1/[8(2\pi i)^3]$. [Normalized contour residues](explore/line-defects/README.md#the-commuting-current-coefficient-in-735). The two-loop divergence for $G\ne SO(8)$ in 2111.08879 §4.3 is tracked [here](explore/twistor-type-i/README.md).

### T4. Kodaira–Spencer on $\mathbb C^3$

Sources [lit]: Costello–Paquette 2001.02177 (3.3.10); Fernández–Paquette–Williams 2404.14318 (2.3.16).

$$S=\frac12\int\epsilon^{ijk}\,\bar\partial\mu_i\,(\partial^{-1}\mu)_{jk}+\int\alpha\,\bar\partial\gamma+\frac16\int\epsilon_{ijk}\,\mu_i\mu_j\mu_k+\int\alpha\,\mu_i\,\partial_{z_i}\gamma .$$

The $\mu$ sector is not of BF type: bulk graphs exist at every loop order, with four legs [drv].

Problems:

1. propagator of $\mu$ in the current representation (E3);
2. one loop: the box with the two vertices; benchmark: the loop computations of 2404.14318 App. A (not read);
3. two loops: $(1,3,3)$ and $(2,2,3)$ with these vertices;
4. For $R=H^\bullet(K3)$ or $H^\bullet(T^4)$ all 16 vanish already on internal labels, as do both two-loop graphs and all higher connected bulk graphs. [Degree count and check](explore/kodaira-spencer/README.md). In the full minimal theory with $R=\mathbb C$, polarizing all four slots cancels every connected bare local contact at fixed Schwinger times. [All-order proof and 59 graph checks](explore/kodaira-spencer/README.md#minimal-theory-with-scalar-coefficients).
5. same steps for hCS with type I Kodaira–Spencer on twistor space (2111.08879 (2.2.1)–(2.2.5); 2306.00940 (2.18)–(2.35)); the cancellation of higher-loop anomalies is known only for $\mathfrak{so}(8)$ (2111.08879 §2.2).

### T5. $\mathbb C^5$

Sources: Baulieu 1009.3893 and `twist10d.py` for hCS on $\mathbb C^5$; Costello–Li 1606.00365, 1505.06703, 1905.09269 (ids checked, bodies not read).

1. hexagon: $|\mathcal I(\lambda;0)|=|\det(\lambda_1,\dots,\lambda_5)|/5!$ [chk]; anomaly coefficient against 1905.09269;
2. two loops: $(1,5,5),(2,4,5),(3,3,5),(3,4,4)$, $|V|=10$, $|I|=11$, eight legs for cubic vertices [drv].

3. Minimal BCOV on the fivefold in the pure closed sector: the primary marked scalar square has nonzero bare differential at every loop order. The [type IIB calculation](explore/iib-bcov/README.md) gives an exact all-loop matrix element. The eight-primary two-loop theta sum requires a seven-field action term which also cancels the integrated cubic defect. [Genus-two symbols](explore/iib-bcov/README.md#the-genus-two-bulk-action-and-integrated-cubic). Every constant top-polyvector scalar background lifts simultaneously through genus two in the local model. [Scalar deformations](explore/iib-bcov/README.md#the-integrated-bilinear-and-simultaneous-scalar-deformations) The maximal scalar input derivative order at loop number $n$ is $n+1$. The scalar square at three loops has a nonzero $\zeta(3)$ coefficient after summing every graph on the repeated-momentum test. [Three-loop scalar square](explore/iib-bcov/README.md#the-three-loop-scalar-square). A nonzero transverse bulk coefficient containing $\zeta(3)$ requires a nine-primary genus-three action term on a ten-wave family. [Three-loop bulk calculation](explore/iib-bcov/README.md#a-transverse-three-loop-bulk-term)

### T6. Holomorphic–topological ($T=1$)

**T6a. 5d Chern–Simons on $\mathbb R\times\mathbb C^2$.** Sources [lit]: Costello 1610.04144 §5, 1705.02500 §5, §9; Oh–Zhou 2002.07343 (23), (87)–(89).

$$S=\frac1{\epsilon_1}\int_{\mathbb R\times\mathbb C^2}dz_1dz_2\,\mathrm{Tr}\Big(\tfrac12A\star_{\epsilon_2}dA+\tfrac13A\star_{\epsilon_2}A\star_{\epsilon_2}A\Big),\qquad P=\frac{\bar z\,d\bar w\,dt-\bar w\,d\bar z\,dt+t\,d\bar z\,d\bar w}{(t^2+|z|^2+|w|^2)^{5/2}} .$$

Benchmarks [lit]: $2^5\pi^2\{X_{\alpha\beta},X_{\gamma\delta}\}=\hbar(\partial_{z_2}X_{\gamma\beta}\partial_{z_1}X_{\alpha\delta}-\partial_{z_1}X_{\gamma\beta}\partial_{z_2}X_{\alpha\delta})+\dots$ (1705.02500 Prop. 9.3.1); $[t[2,1],t[1,2]]_{\epsilon_1}=\epsilon_1\epsilon_2^2\,t[0,0]t[0,0]$ (2002.07343 (107)); $c_{m,n}=(-1)^{m+n}(\tfrac{8\pi}{15})^2(m+n+1)!$ (2103.00963 (63)).

Problems:

1. the diagram of 2002.07343 §5.3 (four internal propagators, three internal vertices), bounded but not evaluated there (eq. (137));
2. the two-loop graphs $(1,3,3)$, $(2,2,3)$ without defects; a two-loop anomaly is expected for the commutative theory by 2407.08667 (introduction), while 1610.04144 §9 states no anomaly for non-commutative $\mathfrak{gl}_N$;
   [chk, runs 1–2, verified] $K_{133}=-\tfrac{\pi}{20}[p_1p_2][p_3p_4]$, $K_{223}=+\tfrac{\pi}{30}[p_1p_2][p_3p_4]$, $G_\Gamma=-K_\Gamma/(64\pi^5)$, $[O_2]=\tfrac{\pi N_{\rm norm}}{48}[F_{223}]$; $[F_{223}]\neq0$ for $\mathfrak{sl}_N,\mathfrak{gl}_N$ ($N\ge3$), $\mathfrak{so}_{5,6,7,8}$, $\mathfrak{sp}_6$, $G_2$, $F_4$; $=0$ for $\mathfrak{sl}_2,\mathfrak{gl}_{1,2}$; cancelled by $J_1=a[p_1p_2]d_{abc}$ for $\mathfrak{sl}_N,\mathfrak{gl}_N$ ($N\ge3$); conditional on H1 (BRST ≅ CE) [open]. Order $\hbar^3$ [open, run 3]. Draft `notes/flagship/flagship.tex`.
   [chk, run 3, verified] H1 in $\eta$-invariant form: $H^g(V^\eta)=H^{g+3}(C^0_{\mathfrak g})$; translation invariance alone gives 8 summands. $[F_{223}]=0\iff\mathfrak g=\mathfrak{sl}_2$ for all simple $\mathfrak g$, $\mu(F_{223})=\tfrac{3d^2(d-3)\alpha\beta\gamma}{2t^3}$; no closed $J_1$ of any bidegree for simple $\mathfrak g\neq\mathfrak{sl}_2$ with $\dim(\mathfrak g^{*\otimes3})^{\mathfrak g}=1$. Order $\hbar^3$: $[O_3]=0$ for $\mathfrak{sl}_N$ ($o_3=a\pi\nu\,d(-\tfrac{2N^2+5}{30}[p_1p_2]^3d)$, forced by $C$-parity) and for $\mathfrak{gl}_N$ with $SL_2$-invariant $I_2$. $\varepsilon U=+1$ [derived, unverified]. Order $\hbar^4$ [chk]: all 50 bare kernels, the decorated diagrams, and the full traceless QME have been evaluated for arbitrary conserved momenta. [Cochain replay](deepen/o3/higher_orders/four_loop_reduction/README.md#all-external-momenta). Draft v2 `notes/flagship/flagship.tex`.
3. $f_2(\kappa)$ in $\hbar_{\rm QFT}\mapsto2^5\pi^2\hbar_{\rm comb}+f_2(\kappa)\hbar_{\rm comb}^2+\dots$ (1705.02500 Thm 12.0.1), left undetermined there.

**T6b. Twisted 11d supergravity on $\mathbb R\times\mathbb C^5$.** Source [lit]: Raghavendran–Saberi–Williams 2111.03049 (7), (9), (13), Thms 1.1, 1.4; the identification with the twist is their Conjecture 0.1.

$$S=\int^\Omega\Big[\beta\wedge(\bar\partial+d_S)\nu+\gamma\wedge(\bar\partial+d_S)\mu+\beta\wedge\partial_\Omega\mu+\tfrac12\tfrac1{1-\nu}\mu^2\vee\partial\gamma\Big]+\frac g6\int\gamma\wedge\partial\gamma\wedge\partial\gamma .$$

One loop: finite and anomaly-free (Thm 1.4, no graph evaluated). For $g\ne0$ graphs exist at every loop order; existence of the quantization is left open there (§1.6).

Problems:

1. propagators of the $(\beta,\nu)$, $(\gamma,\mu)$ system with the mixing $\beta\,\partial_\Omega\mu$ in the representation of E8;
2. two loops: the five theta graphs $(1,6,6),(2,5,6),(3,4,6),(3,5,5),(4,4,5)$, $|V|=12$, $|I|=13$, with the vertices of $S$.

**T6c. Poisson BF on $\mathbb R\times\mathbb C^2$** (5d $\mathcal N=1$ supergravity, Elliott–Williams 2008.02302 Conjecture 4.1) [lit]: BF type, anomaly-free by their Thm 3.4; $Q_n$ on local operators as in T1.

**Tree-level exact** [lit]: the non-minimal twist of 11d supergravity, Poisson–Chern–Simons on $\mathbb C^2\times\mathbb R^7$ (2111.03049 Thm 4.1), has $T=7$.

### Not admitted (no action in the sources opened)

4d $\mathcal N=2$ Poincaré supergravity; 6d $\mathcal N=(1,0)$, $(1,1)$ Poincaré supergravity; conformal supergravities, for which only the symmetry algebras are given ($\mathrm{Vect}^{\rm hol}(\mathbb C^2)$, vector fields on $\mathbb C^{2|1}$, $E(3|6)$: 1910.04120, 2410.08176) [lit].

## 5. Order

| stage | work | needs | output |
|---|---|---|---|
| S0 | E1, E2, E3, E4 | — | `code/` engine with checks |
| S1 | T2.1, T2.2 | E1, E2, E4 | $Q_3$ on $\mathbb C^{3|2}$ |
| S2 | T3.1 | E3 | box anomalies on twistor space |
| S3 | T1 | E3, E4, E5 | Poisson BF on $\mathbb C^2$ |
| S4 | T4.1–T4.3 | E3 | Kodaira–Spencer through two loops |
| S5 | E6 | — | $K_4$, prism, $K_{3,3}$ |
| S6 | E7, T3.3–T3.4 | E3 | defect diagrams |
| S7 | T5 | E1 at $d=5$ | hexagon, two loops |
| S8 | E8, T6a, T6b | E1, E3, E7 | two-loop graphs at $T=1$ |

Each stage ends with a note in the format of [6d], its scripts and checks in `code/`, and an entry in `LOG.md`.

## 6. Background to be supplied

1. T1: the twist of 4d $\mathcal N=1$ supergravity as Poisson BF is a conjecture (2008.02302); the dictionary between its fields and the supergravity multiplet; the relation between $\mathbf S_{\dot\alpha}$ of 2311.04304 and the coupling to $A$; the sign in 2311.04304 (1.14) against (3.14); $\Lambda_{\rm vector}$ (3.15).
2. T3: the form of the propagator and of $\bar\partial^\dagger$ in 2208.12701 (4.18)–(4.28); the prefactor $\frac1{4\pi}$ or $\frac1{4\pi i}$ of $S_\eta$; the propagator of $\eta$.
3. T4: the printed form of 2001.02177 (3.3.10) and 2404.14318 (2.3.16); the constraint $\partial\mu=0$ and the propagator for $\bar\partial\partial^{-1}$; the ring $R=H^\bullet(K3)$.
4. T5: the BCOV action, propagator and anomaly coefficients of Costello–Li.
5. T6, E8: the closed formula for $\alpha_\Gamma$ (2408.03192 Thm 1); the gauge and normalization of the 5d propagator (2408.15712 (2.9)–(2.12) against 2002.07343 (87)–(89)); the propagators of the 11d system; the parities and the $\mathbb Z/2$ grading of 2111.03049 §1.3.
6. E7: the bulk–defect couplings, e.g. 2201.02595 (7.0.1), 2211.06417 (5.1).

## 7. Conventions

- Sources: TeX first; web extractions are used only to locate equations.
- Numbers: `fractions.Fraction` or `sympy.Rational`; no floating point in a check.
- Every kernel is compared with the direct Kirchhoff density at random rational points before it is used.
- `LOG.md`: newest unresolved difficulty first.

## 8. Phase 2: loop corrections in twisted supergravity (added 2026-11-21)

Admission: a bulk loop order beyond one exists (not BF type, or $T=1$ with a non-BF vertex), action and propagator in a source. Tools carried over: theta kernels at $T\in\{0,1\}$, $H+T\le6$ (flagship Thms 2.1–2.3); derivative $\to$ current insertion; local cohomology $H^g(V^\eta)=H^{g+3}(C^0)$ (flagship Thm 3.2) and its proof scheme; certified-rank CE code; Lie (super)algebra builders.

Sources [src, run-3 check]: RSW 2111.03049 Thm 1.4 (one loop finite, anomaly-free, via GRW21), §1.6 ($J=\tfrac16\gamma\partial\gamma\partial\gamma$ gives graphs at every loop order; higher orders deferred), Thm 2.1 (flat-space symmetry: trivial central extension of $E(5,10)$, §2.2.1: even part $\mathrm{Vect}_0(\mathbb C^5)$, odd part $\Omega^2_{cl}(\mathbb C^5)$). Costello–Li 1905.09269: type I KS + hCS consistent only for $SO(32)$ ($\mathbb C^5$), $SO(8)$ ($\mathbb C^3$); one-loop selection, higher loops then cancelled by their method. Costello 1610.04144 §1.2–1.5: twisted M-theory on $TN_k\times\mathbb R^3\times\mathbb C^2$ $\simeq$ non-commutative 5d CS, $\mathfrak{gl}_{k+1}$, existence for $|\delta|\ll|e|$.

### Stage M (5d CS inside twisted M-theory; uses the flagship directly)

M1. $\hbar^2$ analysis in the full translation-invariant complex $V=\bigoplus_S\psi_SH^{\bullet+|S|}(C^0)$ (drop $\eta$, $V_{ht}$): closed ghost-0 $J_1$ with $S\neq U$ components ($G_\epsilon$, $T_\lambda$ types), brackets into the $S=U$ summand, positive-weight $O_n^{(w)}$ (referee-2 R1). Decides whether the no-go and (1.6) of `notes/interpretation/J1_and_5d_CS_level.md` survive.
M2. Non-commutative $\mathfrak{gl}_N$ at order $\delta^2$: Moyal vertices to all orders in $e$, counterterms $e^{-r}\delta^2$, $r=1,2$; coefficient of $e^{-1}\delta^2\pi^*D_{\mathfrak{sl}}$ (leading value (3.1) of the note); compare Costello Thms 9.0.2–9.0.3 (uniqueness).
M3. Two-loop defect diagrams: Oh–Zhou 2002.07343 §5.3 (bound (137)), $f_2(\kappa)$ of 1705.02500 Thm 12.0.1; theta kernels with defect insertions (E7); compare with the M2-brane algebra side (Gaiotto–Oh 1907.06495; 2103.00963 (63)).
M4. $\mathfrak{so}_N,\mathfrak{sp}_N$ (O6 planes): orientifold-projected Moyal product; whether $e\neq0$ cancels $[O_2]$ for class (B) through $e^{-r}\delta^2$ terms.

### Stage S (twisted 11d supergravity on $\mathbb R\times\mathbb C^5$, $(H,T)=(5,1)$)

S1. Fields, $\mathbb Z/2$ grading, propagators of $(\beta,\nu)$, $(\gamma,\mu)$ with the mixing $\beta\,\partial_\Omega\mu$, heat-kernel gauge of GRW21, in the form of E8.
S2. Vertices of $S_{BF,\infty}$ ($\tfrac12\tfrac1{1-\nu}\mu^2\vee\partial\gamma$, all $\nu^k$) and $gJ$; selection rules of run 1 (all 13 edges $\mu\to\gamma$; sectors $g^m\mu^{11-2m}\gamma^{2m-1}\nu^K$).
S3. Odd loops: $\rho_{\rm odd}=0$ reproduces RSW Thm 1.4 at one loop; three-loop vanishing.
S4. Local cohomology: analogue of flagship Thms 3.1–3.3 with $C^0$ replaced by reduced CE cochains of the formal $L_\infty$ algebra of RSW (central extension of $E(5,10)$) modulo translations; $\eta$-invariant form; proof by the h1 scheme.
S5. Obstruction space at the weights of $O_2^{11d}$ by $SU(5)\times U(1)$ characters and certified ranks, replacing the $6.3\times10^{16}$-monomial jet basis of run 1.
S6. Two-loop class: five thetas $(1,6,6),(2,5,6),(3,4,6),(3,5,5),(4,4,5)$ with current insertions (up to 23 derivatives); $[O_2^{11d}]$; if nonzero, order-$\hbar$ cancelling terms in $H^0$ at that weight. Hypothesis to test: the cancelling term is the twisted image of a one-loop topological coupling of M-theory ($C\wedge X_8$), as $J_1$ is the twisted image of the 5d CS level (note §1).

### Stage K (type I / Kodaira–Spencer, $T=0$)

K1. Type I KS + hCS on $\mathbb C^5$: hexagon coefficient vs 1905.09269; explicit two-loop counterterms for $SO(32)$ (theta kernels at $(5,0)$, eight legs).
K2. Twistor space, hCS + type I KS: the two-loop divergence for $G\neq SO(8)$ of 2111.08879 §4.3; run-1 result (even-$n$ hCS obstruction on $\mathbb C^3$ vanishes for reductive $\mathfrak g$) as input.
K3. KS on $\mathbb C^3$ (T4 of §4): bulk two loops, $(1,3,3)$, $(2,2,3)$ at $(3,0)$ with $\partial^{-1}$ kinetic term.

### Stage I (physical interpretation of the 5d results)

I1. The pure minimal 5d SYM twist is the cotangent BF theory $T^*[-1]\operatorname{Map}(\mathbb C^2\times\mathbb R_{\rm dR},B\mathfrak g)$. Its $A$-only level insertion $-3F_D$ occurs only at tree level. Every Lie algebra section of the cotangent target pulls its kinetic pairing back to zero on a simple factor. Matching a supersymmetrically completed level to the flagship and fixing $c_{\rm tw}$ and $\hbar/g^2$ remain open. [ESW, §9.1.1](https://arxiv.org/pdf/2002.10517#page=67) · [Counting, trace identity and exact replay](notes/interpretation/J1_and_5d_CS_level.md#cotangent-field-counting-and-the-kinetic-pairing)
I2. With $\kappa=\mathrm{tr}$ and $\hbar=\lambda/N$, the single-trace contribution is a CE boundary. The double-trace coefficient is $U\lambda^2/(512\pi^4N)$ and the extra $\mathfrak{gl}_N$ triple-trace coefficient is $-U\lambda^2/(256\pi^4N^2)$, both planar. [Exact formal-$N$ calculation](notes/interpretation/README.md)
I3. The four-point probe equals $-2d$ times the $X_2$-channel trace of $F_{223}^{uv}$ [derived; 22 contraction tables checked]. Identification of an untwisted observable with this trace remains open. See `notes/interpretation/x2_channel_check.py`.

### Order

| stage | needs | output |
|---|---|---|
| M1, I1 | flagship code, h1 jet code, ESW TeX | physical status of the no-go and of $J_1$ |
| I2, I3 | simple-all code | $N$-scaling; $X_2$ reading |
| M2 | M1 | $e^{-r}\delta^2$ counterterms of twisted M-theory |
| S1–S3 | RSW, GRW21 | 11d propagators, vertices, odd-loop checks |
| S4–S5 | h1 scheme | 11d obstruction space |
| S6 | S1–S5, theta kernels at $(5,1)$ | two-loop status of twisted 11d supergravity |
| M3, M4 | M2, E7 | two-loop M2-brane algebra; orientifolds |
| K1–K3 | E1 at $d=3,5$ | type I / KS two loops |

## 9. Twisted holography of the holomorphic twist of $\mathcal N=4$ (added 2026-11-21)

Sources: LLC arXiv:2611.24160 (all-loop $Q_n$, multi-graviton protection, census); Costello–Li 1606.00365 §9 (BCOV on $\mathbb C^5\setminus\mathbb C^2$ with $NF$); Budzik–Gaiotto 2611.18539 (finite-$N$ BPS Hilbert space, Koszul dual of the D-brane system).

| stage | content | status |
|---|---|---|
| H1 | planar single-trace spectrum: bulk $\operatorname{coker}[N\{F,\cdot\}]$ vs boundary $\bigoplus_pB^p$ | [chk] $t^0$, 60 slots; differs by $B^1$ ($U(1)$) |
| H2 | $t^{\ge1}$ towers; regular modes vs sources | [chk] towers through $t^3$; KMMR index through $t^{30}$; no perfect pairing (excess on $M_1\cup\{1^5\}$). Regular modes $\mathfrak g=\mathbb C1\times PV^{\ge1}_\Delta(\mathbb C^{2|3})\supset\mathfrak{psl}(3|3)$; excess $=$ $U(1)$-letter sources $u_\beta$ and the central $1$; $U(1)$ shifts Higgsed against $B^1$. [open] $\ell_4(F,F,\cdot,\cdot)$. [Dictionary](explore/twisted-holography-N4/dictionary/REPORT.txt) · [Symmetry](explore/twisted-holography-N4/symmetry/REPORT.txt) |
| H3 | P1 on the stable relative single-trace complex | [drv] Filtered linear conjugacy to $Q_0$; explicit nonzero planar $Q_2$ on four letters. A local product-compatible conjugacy is open. [Homotopy and color calculation](explore/twisted-holography-N4/NOTE.md#4-loop-orders-drv) [chk] The product-compatible version (M) fails at $\hbar^2$, $O(N^{-2})$: $\Phi^3_{Q_2}\ne0$ [Verify](explore/twisted-holography-N4/verify/VERIFY.txt) |
| H4 | One-loop BCOV on punctured $X$ vs nonplanar boundary $Q_n$ | [open] The boundary ratios begin at $6/N^2$. On bulk primary scalar probes the genus-two and genus-three tree contributions are fixed, and ordinary ambient loops at that arity vanish. A puncture-sensitive theta edge survives both angular tensors. Its physical integrated amplitude and comparison to the boundary remain open. [Sewing and puncture calculation](explore/twisted-holography-N4/README.md#sewing-the-local-primary-vertices) |
| H5 | finite $N$ vs giant-graviton quantization | [chk] BG $k=1$ complex $=E_1^{U(N)}$ ($N\le4$, low levels); $k=2$ removes the $k=1$ excess at levels 16–20 (computed sectors, mod $p,p'$); onset rule at level $2(k+1)(N+k+1)$. $\mathcal O_0$ (level 24) not reached: $\dim V_T=1.06\cdot10^{12}$ at $k=2$. [Finite N](explore/twisted-holography-N4/finiteN/REPORT.txt) · [$k=2$](explore/twisted-holography-N4/finiteN_k2/REPORT.txt) |

## 10. Circuit-reduction formalism by $(H,T)$ (added 2026-11-22)

$D=H+T$; admissible graphs: $(D-1)|I|=D|V|-D-1$, $(D-1)|I_S|\le D|V_S|-D-1$, cycle length $\ge D+1$. [chk] At $D\ge3$ every tight composite move is a two-terminal reduction ($|\partial C|=2$); source-free stars still occur inside two-terminal bridge subnetworks, so polynomial reducibility $\supsetneq$ SP at generic sources ($D=3$: 1 graph at $n=3$, 11 at $n=4$; $D=4,5,6$: 1 at $n=3$). With every vertex sourced, polynomial $=$ SP in all computed ranges. Census, cores and the theory × computation atlas: [atlas](explore/circuit-general/atlas/REPORT.md).

| $(H,T)$ | theories | established | open |
|---|---|---|---|
| $(1,0)$ | 2d chiral | $|V|=2$, bananas; parallel move; all moments rational [chk] | — |
| $(2,0)$ | 4d $\mathcal N=1,2,4$ holomorphic twist; Poisson BF | fans, spiders, cubic product forests, all loops (series–parallel + single/paired/sourced $Y$–$\Delta$) [LLC]; non-SP from $n=4$, periods $\mathbb Q(\zeta_3,\zeta_5)$; four-loop Ward reconstruction of all 13 kernels [chk] [nonSP](explore/nonSP/README.md) | — |
| $(3,0)$ | 6d $(1,1)$ hCS on $\mathbb C^{3|2}$; KS on $\mathbb C^3$; twistor BF, WZW$_4$ | chain, bubble moves; box ladder all $n$; $Q_3$; graph generator at every order; non-SP from $n=3$, periods $\mathbb Q(\zeta_3)$ [chk] [six-d-Q3](explore/six-d-Q3/README.md); four-loop bulk (hCS/KS, 50/50) and twistor $Q_4$ (110/110) exact with first and second current insertions, values in $\mathbb Q+\mathbb Q\zeta_3+\mathbb Q\zeta_5$ [chk] [atlas](explore/circuit-general/atlas/REPORT.md) · [assembly](explore/circuit-general/assembly/REPORT.md); five-loop cores all tabulated at $q=3,4$ [fiveloop](explore/circuit-general/fiveloop/REPORT.md); 6d $Q_4$ assembled (15 kernels; non-tight plan by block-triangular density), $(Q^2)_{\hbar^4}=0$ in $\mathfrak{so}(5)$ on $b^5$ (16 functionals, rank 15/15) fixes $\kappa_\Gamma=\sigma_\Gamma/|\mathrm{Aut}\,\Gamma|$, $Q_5$ kernels 49/49 [chk] [6d-nilpotency](explore/physics/6d-nilpotency/REPORT.md); hCS on $\mathbb C^3$: $\omega_3=-\tfrac{3}{64\pi^9}\Delta^3T_3$ ($\zeta_3$ cancels), $T_3\supset\tfrac{(N-8)(N-5)(N+4)}{829440}\mathrm{tr}x^4$ for $\mathfrak{so}(N)$, $\omega_4\equiv0$ (50 graphs) [chk] [anomalies](explore/physics/anomalies/REPORT.md) | $n=5$ non-SP kernels (assembly size); $n\ge6$; order-5 nilpotency; 10d operator layer; $Q_2O_1$, KS exchange |
| $(4,0)$, $(5,0)$ | 8d, 10d SYM; BCOV; type I | polygons, two-loop thetas; three-loop bulk on $\mathbb C^5$ (83 graphs) and cone$(5,3)$ with first and second insertions [chk] [assembly](explore/circuit-general/assembly/REPORT.md); 8d $Q_3$ (12) [chk] | four loops at $D\ge4$ |
| $(H,0)$ | — | SP with sources on all vertices (`code/htk/spred.py`); loop- and dimension-raising identities; signed current coordinates for every graph; first irrational loop order [chk] [ht-families](explore/ht-families/README.md); cores: $W_n$ radial masters in closed form, prism and $K_{3,3}$ periods at every $q=3..7$ and $q=\tfrac72,\tfrac92$ (prism) by regulator-free relations (FORCER bank 1481/1481); exact evaluator: $n\le4$ at $D=2$ with insertions, $n=3$ at $D=3,4$, $n=4$ at $D=3$ bare [chk] [cores](explore/circuit-general/cores/REPORT.md) | $n=3$, $D\ge5$ kernels; five-loop cores beyond $W_n$; non-MZV expected from 8-loop cores |
| $(0,1)$ | topological QM | bananas, all $n$: $K_N=\pi^{n/2}$ ($n$ even), $I^\star_N=2^{-N}$, $f\star g=m\circ e^{\frac\hbar2\eta}(f\otimes g)$ (Moyal) [chk $N\le7$] [topological](explore/circuit-general/topological/REPORT.md) | — |
| $(1,1)$ | 3d $\mathcal N=2$ HT twist | fans, all even $n$, closed Taylor coefficients [chk] [fans-ht](explore/fans-ht/VERIFY.md); $E>0$ kernels vanish [chk]; single, paired, sourced $Y$–$\Delta$ at $T=1$: $U_\Gamma=U_{\rm red}t_{vx}/(r_1r_2)$, local factor $(r_1r_2T_{ux}T_{wx}/Q)^{1/2}$ (`code/htk/star.py`, 155 sites) [chk]; spiders, SP product forests at $n=4$ exact; edge-dipole vanishing (122 cases) [chk]; $Q_n$ graphs: fans, tadpole$(n+1,0)$ ($K=\tfrac\pi3l^2$ at $n=2$), spiders for $\deg W\ge4$; Wess–Zumino: $Q_{n\ge1}=0$ [chk] [ht-d2](explore/circuit-general/ht-d2/REPORT.md); non-SP four-loop forests and prism exact, $K\in\mathbb Q\pi^2+\mathbb Q\pi^4$; tadpole$(5,0)=\tfrac{\pi^2}{630}\lambda^4$ [chk] [cores](explore/circuit-general/cores/REPORT.md); $Q_n=N_n\hat Q_n$, $N_n=\epsilon_n(\pm\hbar/4\pi i)^n$; $\{Q_0,Q_2\}=0$ in 8 theory families; $(Q^2)_{\hbar^4}=0$ at $z=0$ reproduces tadpole$(5,0)$; tadpoles $\propto d^V_{abc}$; two-loop lifts in $U(1)_0$ + chiral = induced CS level [chk] [3dN2](explore/physics/3dN2-supercharge/REPORT.md) | $\epsilon_2$; $Q_4$ with derivative letters; non-abelian $d_2$ with $d^V\ne0$; $n=6$ products |
| $(2,1)$ | 5d CS, 5d $\mathcal N=1$ | thetas; $[O_2]$, $\hbar^3$, $\hbar^4$ (50 kernels) [chk] [flagship](notes/flagship/flagship.tex) | — |
| $(3,1)$, $(4,1)$ | 7d, 9d SYM | thetas; $Q_{\rm odd}=0$ [chk] | — |
| $(5,1)$ | twisted 11d supergravity | thetas $c_{pqr}$ [chk]; two-loop order-$g$ defect: five theta shapes sum to $1008\pi$ on the three-plus-six cycle, unrooted $-\tfrac{63g}{256\pi^{10}}$, rooted $Q_2\mathcal I$ $-\tfrac{189g}{64\pi^{10}}$; not exact among local two-loop counterterms; order-$g$ self-bracket of one-loop counterterms $=0$ on the detector (ranks 3, 80, 83 mod 1009, 251) [chk] [current complex](explore/sugra-11d-qme/collaboration/current_complex/README.md) | one-loop finite counterterms beyond the detector; comparison with HPR 2506.02246 §4.3.2; $C\wedge X_8$ reading |
| $(H,1)$ | — | SP: $K\in\pi^{n/2+2E}\mathbb Q[\lambda]$; thetas $\tfrac{D_\theta}2\prod_P\Gamma(H+\tfrac32-n_P)/(n_P-1)!/\Gamma(H+\tfrac32)$; two-loop current polygon; fans, general $H$ [chk] | $\kappa$ transfer (run) |
| $(H,T\ge2)$, $(0,T)$ | 4d CS, 4d $\mathcal N=2$ HT $(1,2)$; 2d topological $(0,2)$; 3d CS, RW, 3d $\mathcal N=4$ topological $(0,3)$; DW $(0,4)$ | $\alpha_\Gamma=(-1)^n(4\pi)^{-n/2}\rho_n$; $\rho_n\wedge\rho_n=(-\tfrac12)^n\det(X\wedge X)=0$, $X=M^{-1}dM$, $\operatorname{tr}X^{2j}=0$ $\Rightarrow$ bulk kernels vanish at $T\ge2$ [chk $D\le4$, $n\le3$]; $\kappa$ transfer under chain and sourced $Y$–$\Delta$ moves [chk]; half-plane Poisson sigma model: Kontsevich order-$\hbar^2$ weights $1,\tfrac12,\tfrac13,-\tfrac13,-\tfrac16$ [chk-num] [topological](explore/circuit-general/topological/REPORT.md) | 4d CS with Wilson lines; exact boundary weights; $\kappa\ne0$ slice dependence |
| defects | line defects at $T=0$ | ordered interval integrals, all loops [chk] [line-defects](explore/line-defects/README.md) | $T\ge1$ defects |

Beyond circuit moves (2026-11-23): symmetry group of weighted core periods $I_q(G;\alpha)$ [chk, 0 mismatches in $\sim1.5\cdot10^5$ checks] [symmetry](explore/circuit-general/symmetry/REPORT.md): planar duality $I_q(G^*;\alpha)=I_q(G;q-\alpha\circ\pi)$ (S$\leftrightarrow$P, Y$\leftrightarrow\Delta$, $n\leftrightarrow|V|-1$); conformal completion in $d=2q$ ($\beta=q-\alpha$, balanced $\hat G$, invariance under the choice of $0,1,\infty$); uniqueness and its dual; product at 3-vertex cuts; weighted twist; dimension shift; two-point numerator identity. Reach: K5$-$e (FT$\to$prism), W5, 2-connected five-loop skeletons for all weights; other five-loop cores and K5 on codimension-$\ge1$ conditions (fractions 0.03–0.53 at $q=7$); P($Z_5$) a fixed point of the group.

Graphical functions in $d=2q$ (2026-11-23) [chk] [graphical](explore/circuit-general/graphical/REPORT.md): single-valued polylogarithms (letters $0,1$, MZV $w\le11$), edge appending via $M_\lambda=\partial\bar\partial-\lambda T$ in the pole-free representation $f=\sum_kT^k(w^{-1}\Phi_k)$; vertex elimination with function-valued response. Reach: all weights of K5$-$e, W5, $Z_5$, V7_0, V7_2 at $q=3,4$ (and $Z_5$ at $q=7$, $P(Z_5)=\tfrac{441}8\zeta_7$); K5 $\iff\exists v:\sum_{e\ni v}\alpha_e=2q$; 259 six-loop K5-cored periods at $q=3$ (values to weight 9: $\zeta_3^3,\zeta_9,\dots$). Frontier: generic K5 requires the non-conformal one-loop four-point star, a two-variable function (non-constructible $\bar K_6$).

Generic K5 (2026-11-23) [chk] [k5](explore/circuit-general/k5/REPORT.md): $J_q(\beta)=\Pi(K_6;\beta)$ on $\mathcal C_q$; level $q$ from level $q-1$ by vertex IBP (X) and dimension shift (DS); inputs $J_2=20\zeta_5$, $1$, $\zeta_3-\tfrac13$, $-4\zeta_3+5\zeta_5$; all convergent weights $q\le8$ (generic included) in $\mathbb Q+\mathbb Q\zeta_3+\mathbb Q\zeta_5$; $1.07\cdot10^7$ rows exact. K5 family = complete three-loop massless two-point family. Open: all-$q$ proof; origin of one input relation.

Locality of Ohm's law at kernel level (2026-11-23) [chk] [locality](explore/circuit-general/locality/REPORT.md): for any subgraph $S$ with $k$ terminals, $\Omega_\Gamma=\pm\,\Omega_S\wedge\pi^*\Omega_{\rm out}$, $\pi:t\mapsto(t_O,Z=L_r^+,h)$, $L_r$ = Schur complement on $\partial S$, $h$ = harmonic transfer of internal sources; $T=0$ and $T=1$ (Gaussian Dirichlet source $\xi$ = the $T=1$ contact term). $k=2$ dressed edge: $K_\Gamma=\kappa\int D_S\,K_{\Gamma'}[J\to-\partial_{z'}]$, moments $=\frac{\Gamma(A+B)}{\Gamma(A)\Gamma(B)}I_{A+B}(S+e_*)$ (power $-B$ of $R_{\rm eff}$ = weight $A$ on the closing edge); reproduces all 23 bridge and 18 bubble $D=3$, $n=4$ kernels, 5–10× faster. Degree-$k$ vertex: response = toric subvariety $c_{ab}c_{cd}=c_{ac}c_{bd}$ (codim $k(k-3)/2$), no product chart; generic K5 stays two-variable.

Runs (2026-11-22): `explore/circuit-general/ht-d2/` ($(1,1)$: $Y$–$\Delta$ at $T=1$, spiders, forests, 3d $\mathcal N=2$ $Q_n$); `explore/circuit-general/topological/` ($T\ge2$ vanishing in circuit form, $\kappa$ transfer, $(0,1)$ Moyal, one defect/boundary case).
