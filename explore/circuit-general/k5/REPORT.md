# K5 (six-loop core): all weighted periods, generic weights included

Tags: [chk] means an exact check, with script and log in `P/explore/circuit-general/k5/`. [drv] means derived here. [src] gives a source location. [open] means not settled. Paths below are relative to `k5/`. K5 edge order is 01,02,03,04,12,13,14,23,24,34. K6 has vertices 0..5 and edges (i<j) in lexicographic order. Every value is exact (Fractions, or modular arithmetic followed by CRT and rational reconstruction, then an exact Fraction check of every row). Numerics are used only as cross-checks.

## 0. Family

$$J_q(\beta):=\Pi(K_6;\beta)=\frac{\Gamma(q)\,I_q(K_6-c;\,q-\beta)}{\prod_{e\in K_6-c}\Gamma(\beta_e)}\quad\forall c,\qquad I_q(K_5;\alpha)=\frac{\prod_{e\in K_5}\Gamma(q-\alpha_e)}{\Gamma(q)}\,J_q(\hat\beta)$$
[src symmetry/REPORT.md §1.1–1.2]

$$\mathcal C_q=\{\beta\in\mathbb Z^{15}:\ \textstyle\sum_{e\ni v}\beta_e=2q,\ \beta_e\le q-1,\ \sum_\triangle\beta\le2q-1\}/S_6$$
- Convergence is the cut condition for |W|=2 and |W|=3; |W|=4 is the complement of |W|=2 [drv from symmetry §1.3].
- Complementary triangles: $\mathrm{cut}(W)=6q-2\sum_W\beta\Rightarrow\sum_W\beta=\sum_{W^c}\beta$ [drv].
- $\beta_e\ge4-2q$ [drv].
- Kinds: T1 means some β=0. T2 means some triangle has Σβ=q. G means neither.
- d: a K5 decompletion exists, i.e. β∈[1,q−1] on K6−c. a: auxiliary, i.e. nonpositive edges that are not in one star.

Position space with $x_0=0$, $x_1=e$, integrating $x_2,x_3,x_4$:
- The integrand has 9 nonconstant propagators $x_v^2,(x_v-e)^2,x_{vw}^2$.
- These are in bijection with the 9 scalar products $x_v\!\cdot\!x_w$, $x_v\!\cdot\!e$.
- So the family is a complete 3-loop massless two-point family [drv].

## 1. Relations

(X): $\int d^dx\,\partial_\mu[(x-x_1)^\mu\prod_v|x-x_v|^{-2\beta_{xv}}]=0$ and $d=\sum_{e\ni x}\beta_e$ give
$$(\beta_{xc}-\beta_{xx_1})J(\beta)=\sum_{i\notin\{x,x_1,c\}}\beta_{xi}\big[J(\beta-e_{xx_1}+e_{xi}+e_{x_1c}-e_{ic})-J(\beta+e_{xi}-e_{x_1i}+e_{x_1c}-e_{xc})\big]$$
- This holds for every ordered triple (x, x₁, c) of distinct vertices.
- Each relation is used only if every term with a nonzero coefficient converges.
- It is rel5 (X) in Π normalization: $2a_1+\sum a_i-3q=\beta_{xc}-\beta_{xx_1}$ [drv; src fiveloop/rel5.py].

(DS): $\Psi^{1-q}=\Psi\cdot\Psi^{-q}$ on K6−c gives
$$J_{q-1}(\beta)=\frac1{q-1}\sum_{T\subset K_6-c}\prod_{e\in T}\beta_e\;J_q\big(\beta+1_T\big)^{\rm completed\ at\ c}$$
[drv from symmetry §2 (DS)]

Check: all 166 rows at level 3 (X, plus DS from 2 to 3) hold on independently computed level-3 values [chk `test_q3.log`, `tower.py`; rerun by the orchestrator 2026-11-23]. Rows are generated in C (`xrows.c`) [chk C rows = Python rows, q=4,5].

## 2. Counts [chk `enum6.c` = `k6.enumerate_convergent` for q=2..4; `k5tables.log`]

| q | K6 classes | T1 (d/a) | T2 (d/a) | G (d/a) | K5 S5-orbits: T1 / T2 / G |
|---|---|---|---|---|---|
| 2 | 1 | 0/1 | – | – | – |
| 3 | 15 | 1/13 | 1/0 | 0/0 | 1 / 1 / 0 |
| 4 | 159 | 11/127 | 8/8 | 2/3 | 15 / 13 / 2 |
| 5 | 1259 | 88/931 | 71/116 | 16/37 | 133 / 156 / 23 |
| 6 | 7195 | 496/4951 | 427/848 | 123/350 | 792 / 1084 / 231 |
| 7 | 32187 | 2133/20527 | 2038/4456 | 735/2298 | 3469 / 5480 / 1563 |
| 8 | 117744 | 7444/69773 | 7668/17759 | 3536/11564 | 12201 / 21152 / 8075 |

The K5 orbit totals 2, 30, 312, 2107 and T1 = 15, 133, 792 equal `graphical/k5_reach_fraction.log` [chk].

## 3. Closure (Task 1) [chk `rank_c_q45.log`, `rank_c_q6.log`, `rank_c_q7.log`, `solve_c_q*.log`, `masters.log`]

Modular rank (p = 2³¹−1). Modes:
- A: (X) only, every class unknown.
- B: (X) only, T1 and T2 known.
- B1: only T1 known.
- B2: only T2 known.
- Bd: only the K5-decompletable T1/T2 classes known.
- C: (X) plus (DS) from level q−1, with every level-q class unknown.

| q | A: rank/classes (nullity) | B, B1, B2, Bd: generic d-classes determined | C: determined |
|---|---|---|---|
| 3 | 11/15 (4) | – | 12/15 with J₂ known (3 free) |
| 4 | 156/159 (3) | 2/2 in all four modes | 159/159 |
| 5 | 1255/1259 (4) | 16/16 in all four; 1 auxiliary class free | 1258/1259 |
| 6 | 7192/7195 (3) | 123/123 in all four | 7195/7195 |
| 7 | 32183/32187 (4) | Bd: 735/735; 1 auxiliary class free | 32186/32187 |
| 8 | not run | not run | 117744/117744 (every 4th X row + all DS rows) |

- With T1∪T2 as boundary data, every generic orbit is determined for q=4..7. T1 alone already suffices (q=4..6), and so do the K5-decompletable T1/T2 alone (q=4..7) [chk].
- Stronger result: (X)+(DS) determine level q from level q−1 alone, so T1/T2 values at q ≥ 4 are not needed as input [chk q=4..8].
- Classes left free by C are only the odd-q class $\kappa_q=(-k,-k,\lambda,\lambda,\lambda,-k,\lambda^6,-k,-k,-k)$, with $k=\frac{q-3}2$ and $\lambda=q-1$. It is a harmonic K₃,₃ with two numerator triangles; at q=3 it is a level-3 input. These are fixed in §4 [chk].
- Tower inputs: J₂ plus 3 level-3 classes. 79 of the 455 triples are master sets, and all of them contain κ₃.
  - Example basis: J₂(octahedron) = 20ζ₅; κ₃ = (0,0,2,2,2,0,2,2,2,2,2,2,0,0,0), J = 1; (0,0,2,2,2,1,1,2,2,2,1,2,1,0,0), J = ζ₃ − 1/3; (1,1,1,1,2,1,1,2,1,2,1,1,1,1,1), J = −4ζ₃ + 5ζ₅ = 2 I₃(K5;1,2,2,2,2,2,2,1,2,2) [chk `masters.log`].
  - These satisfy $J_2=4J^{(14)}+16J^{(6)}+\tfrac{16}3J^{(5)}$, a relation that (X)+(DS) do not produce [drv; origin open].
- No genuinely new masters remain for q ≥ 4.

## 4. Values (`tables/K5_q{3..8}.json`, `data/J_q{2..8}_tower.json`)

- Level 2: J = 20ζ₅ (engine = graphical functions) [chk `test_q3.log`].
- Level 3: 15 classes, values in {1, ζ₃−1/3, (10/3)(ζ₃−ζ₅), −4ζ₃+5ζ₅, …}; engine, graphical functions and per5 agree on all 15 [chk `test_q3.log`].

Generic K5 examples [chk]:
$$I_4(1,1,2,2,3,3,3,3,3,3)=-\tfrac5{12}-\tfrac{13}{12}\zeta_3+\tfrac53\zeta_5,\qquad I_4(2,2,2,3,3,3,2,3,2,2)=\tfrac1{360}-\tfrac{13}{30}\zeta_3+\tfrac12\zeta_5$$
$$I_5(K_5;3^{10})=-\tfrac{13}{15120}+\tfrac{19}{720}\zeta_3-\tfrac5{168}\zeta_5$$
$$I_6(1,1,3,5,5,5,3,5,5,3)=-\tfrac{719897}{24494400}-\tfrac{1003}{12600}\zeta_3+\tfrac{38}{315}\zeta_5$$
$$I_7(2,2,3,4,6,3,6,6,4,6)=-\tfrac{229357639}{2309472000}-\tfrac{2840093}{10692000}\zeta_3+\tfrac{109}{270}\zeta_5,\qquad I_7(1,2,5,5,6,3,3,6,6,5)=-\tfrac{60023}{17244057600}+\tfrac1{332640}\zeta_3$$
$$I_8(1,4,6,7,5,7,6,5,3,4)=\tfrac{1897778101}{4903778880000}+\tfrac{3965987}{3783780000}\zeta_3-\tfrac{2671}{1681680}\zeta_5$$

Isolated classes:
$$J_q(\kappa_q)=\int d\mu^3\,(|x_{01}|^2|x_{02}|^2|x_{12}|^2)^k,\quad d\mu=\tfrac{d^{2q}x}{\pi^q}|x|^{-2\lambda}|x-e|^{-2\lambda}$$
- Method: direction average in $\mathbb R^{2q-1}$, then one-point moments as sums of $G_{2q-2\epsilon}(\lambda-a,\lambda-b)$ with ε→0 [drv].
- Results: κ₃ = 1, κ₅ = 23/69984, κ₇ = 75409/9331200000000 [chk `check_aux_odd.log`, `check_aux_q5_exact.py`].
- κ₅ equals the value from its 6 DS rows to level 6, which all agree [chk `fix_u5.log`].
- κ₇ is consistent with all 8,243,378 level-8 rows, including its 6 DS rows [chk `solve_c_q8.log`].
- (−5,1,6,6,6,1,6,6,6,4,4,4,−1,−1,−1) at q=7 = 79/31104000000. It was left free only by the row subsample and was fixed by the full rows [chk `fix_rows_q7.log`].

## 5. Checks

Exact checks of every row on the solution, all with 0 failures [chk `solve_c_q*.log`, `set_aux_q7.log`]:

| q | rows checked |
|---|---|
| 4 | 4,562 |
| 5 | 54,428 |
| 6 | 398,172 |
| 7 | 2,041,216 |
| 8 | 8,243,378 |

Tower values compared with independent direct evaluations: E = symmetry engine, S = graphical functions, P = per5 via K33+e tables or dual prism. There were 0 differences in total [chk `check_direct_q{4,5,6,7}.log`].

| q | sample | T1 agreements | T2 agreements |
|---|---|---|---|
| 4 | all classes | 138/138 (E, S); P 10 | E 8, P 8 |
| 5 | 25 per kind | E 50, S 50, P 25 | E 26, P 25 |
| 6 | 40 per kind | E 79, S 80, P 16 | E 45, P 21 |
| 7 | 8 per kind | E 16, S 16 | E 4, P 2 |

Earlier project values: 9 of 9 equal [chk `check_ref.log`; rerun by the orchestrator 2026-11-23]. Sources are graphical §6/§7, symmetry §3.3 and fiveloop §5:
- q=3: two values.
- q=7: one value.
- q=8: six values, including the four that the engine does not reach.

Tropical sampling with 10⁷ samples per vector [chk `check_tropmc_q4.log`, `check_tropmc_q567.log`, `check_tropmc_q8.log`, `check_tropmc_sym.log`; method src fiveloop/tropmc.py]:
- 20 generic runs at q=4..8 (16 vectors), all within 0.0–1.7σ.
- q=4, both vectors, two seeds: 0.1–1.5σ.
- K5(3¹⁰) at q=5: 1.2σ and 1.0σ.

## 6. Structure (Task 3)

Every $J_q$ for q ≤ 8 is a ℚ-linear combination of the 4 tower inputs (20ζ₅, 1, ζ₃−1/3, −4ζ₃+5ζ₅). Hence
$$I_q(K_5;\alpha)\in\mathbb Q+\mathbb Q\zeta_3+\mathbb Q\zeta_5,\qquad q\le8,\ \text{all convergent }\alpha$$
[drv from chk §3–5]. No ζ₇, no ζ₃², and no non-MZV occurs. The weight is ≤ 5.

Constant content per orbit, as orbit counts [chk `label_stats.log`]:

| q | G | T1 | T2 |
|---|---|---|---|
| 4 | 1+ζ₃+ζ₅: 2 | 1+ζ₃: 5; 1+ζ₃+ζ₅: 8; ζ₃+ζ₅: 2 | 1+ζ₃+ζ₅: 13 |
| 5 | 1+ζ₃+ζ₅: 23 | 1: 1; 1+ζ₃: 31; 1+ζ₃+ζ₅: 101 | 1+ζ₃: 5; 1+ζ₃+ζ₅: 151 |
| 6 | 1+ζ₃+ζ₅: 231 | 1: 5; 1+ζ₃: 202; 1+ζ₃+ζ₅: 585 | 1+ζ₃: 35; 1+ζ₃+ζ₅: 1049 |
| 7 | 1+ζ₃: 10; 1+ζ₃+ζ₅: 1553 | 1: 25; 1+ζ₃: 950; 1+ζ₃+ζ₅: 2494 | 1: 1; 1+ζ₃: 217; 1+ζ₃+ζ₅: 5262 |
| 8 | 1+ζ₃: 83; 1+ζ₃+ζ₅: 7992 | 1: 74; 1+ζ₃: 3399; 1+ζ₃+ζ₅: 8728 | 1: 5; 1+ζ₃: 942; 1+ζ₃+ζ₅: 20205 |

- The non-constructible response $F_c$ (graphical §7) introduces no new constants. By (X), generic values are ℚ-linear in the T1 values of the same level (mode B1) [chk q=4..6].
- The constructible (T1∪T2) and generic subsets have the same value space [chk].
- All q: [open]. Closure of level q from q−1 is a computed fact for q=4..8, not a proof.

## 7. Task 2

No masters remain for q ≥ 4, so no two-variable method was needed. Gegenbauer, Mellin–Barnes and the $F_c$ integration were not implemented. The cross-checks are tropical sampling (§5) and exact G-functions for κ_q.

## 8. Task 4: D=2, n=6 kernel (census[2,6] #309, #350, #382, all with a K5 core)

- Period side: done. `tables/K5_q8.json` contains all 41,428 orbits at q=8 [chk].
- Numerator side, for #309 (8 vertices, 13 edges, chains [1⁷,2³]), not run [drv]:
  - fk5-type whole-graph interpolation needs the lower set {e∈[0,5]¹²: Σe ≤ 35} = 1,788,632,625 points per prime, against 8,209,410 at n=5 (52 s per prime there).
  - The dense cube has 6¹² = 2.18·10⁹ int64 entries, about 17.4 GB, which is over the 5 GB limit. Sampling alone would take ≥ 3 h per prime.
  - The contraction afterwards involves ≤ 8,177,785 labelled core vectors (compositions of 48 into 10 parts in [1,7]), each a table lookup.

## 9. Files (`P/explore/circuit-general/k5/`)

- Code:
  - `k6.py`, `linsys.py`, `bval.py`
  - `enum6.c`, `xrows.c`, `modsolve.c`, `modsolve2.c`, `modrank.c`
  - `solve_level.py`, `solve_c.py`, `rank_c.py`, `rank_x.py`, `tower.py`
  - `fix_u5.py`, `fix_u.py`, `fix_rows.py`, `fix_rows_lib.py`, `set_aux.py`
  - `check_aux_odd.py`, `check_aux_q5_exact.py`, `check_direct.py`, `check_ref.py`, `check_tropmc.py`
  - `k5tables.py`, `label_stats.py`, `masters.py`, `test_q3.py`, `count_classes.py`
- Data: `data/classes_q{2..8}.txt`, `data/J_q{2..8}_tower.json`, `data/J_q{2,3}_direct.json`, `tables/K5_q{3..8}.json`.
- Logs: `*.log`.
- The C binaries and the row files (up to 640 MB) are in the session scratchpad, which the system can clear. They are regenerated from `*.c`.

## 10. Not done / failed

- The first q=7 run with 3 primes gave 5,284 reconstruction failures. It was rerun with 6 primes plus a check prime, giving 0 failures [chk].
- Every-4th-row subsampling at q=7 lost one rank; the full rows closed it (`fix_rows.py`) [chk].
- The Python sparse elimination at q=6 was too slow and was killed; it was replaced by C. The mpmath quadrature for κ₅ was too slow and was replaced by exact G-functions.
- Modes A and B at q=8 were not run (about 1.3 h per mode) [open].
- The origin of the extra input relation in §3 is [open], as is closure for q ≥ 9.
- (L) syzygies, planar duality and (Y) to K33+e were not used; they were not needed.
- The D=2, n=6 kernel was not evaluated (§8).
- Nothing was committed or installed, and no existing files were edited.
