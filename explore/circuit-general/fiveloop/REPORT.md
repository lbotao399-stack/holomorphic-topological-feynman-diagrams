# Five-loop cores: periods and kernels

Tags: [chk] = exact check run (script, log in this folder); [drv] = derived here; [src] = source (location); [open]. Paths are relative to `P/explore/circuit-general/fiveloop/`, with P = project root. No existing files were edited: `msys5.py` is a copy of `cores/msys.py` plus two changes (a planar-duality hook and a degree-signature prefilter). Nothing was committed or installed.

| file | content | log |
|---|---|---|
| `skel5.py` | five-loop skeletons; connectivity, planarity | `skel5.log`, `skel5.json` |
| `census_cores.py`, `core_occurrence.py` | series–parallel core of each census graph; core occurrence at D=2, 3 | `census_cores.log`, `core_occurrence.log` |
| `enum_d3.py` | bulk3(3,5) and twistor Q₅: graphs, cores, classes | `enum_d3.log`, `enum_d3.json`, `check_enum_d3.log` |
| `per5.py` | period evaluator: reductions, 2-sums, duality, tables, (Y) fallback | `check_lower.log` |
| `rel5.py`, `syz5.py` | relations (Y), (X), (P), (C), (L) | `check_rel5.log`, `check_completion.log`, `check_syz5.log` |
| `build5.py`, `build5L.py`, `multi5.py`, `joint5.py`, `buildup.py` | relation systems: single-level, multilevel, joint, build-up | `build*.log`, `multi_*.log`, `joint_*.log`, `buildup_*.log` |
| `tropmc.py`, `numcheck.py` | float cross-checks: tropical sampling, and naive Gamma sampling (fails) | `check_tropmc.log`, `check_k5_tropmc.log` |
| `compare_graphical.py` | exact comparison with `../graphical/tables` | `compare_graphical.log` |
| `fk5.py`, `run_d2.py`, `summarize_d2.py` | D=2 kernels | `check_fk5.log`, `run_d2*.jsonl`, `d2_kernels.md`, `cnum/*.json` (144 MB) |
| `check_gen_rational.py` | rationality test of B-class kernels | `check_B_rational.log` |
| `run_d3_sp.py` | D=3, n=5 series–parallel kernels | `d3_sp_bulk3.jsonl`, `d3_sp_twistor.jsonl` |
| `k5_deltay.py` | n=6 graphs without a circuit chart | `k5_deltay.log` |
| `rec5.py` | recursive single-period evaluation (negative result) | `check_z5_q7_graphical.log` |
| `tables/*.json` | period tables | `table_summary.log` |

## 1. Cores

### 1.1 Skeletons
Constraints: simple, 2-connected, deg ≥ 3, h = 5 ⇒ |E| = |V|+4 and 5 ≤ |V| ≤ 8 [drv]. There are 14 skeletons [chk `skel5.log`].

| name | \|V\|,\|E\| | degrees | κ | planar | triangles | \|Aut\| | remark |
|---|---|---|---|---|---|---|---|
| K5−e | 5, 9 | 4³3² | 3 | yes | 7 | 12 | = prism* |
| W5 | 6, 10 | 5·3⁵ | 3 | yes | 5 | 10 | self-dual |
| K33+e | 6, 10 | 4²3⁴ | 3 | no | 3 | 12 | = K5 after Δ–Y |
| Z5 | 6, 10 | 4²3⁴ | 3 | yes | 4 | 2 | C₇(1,2)−v, self-dual |
| K4⊕K4 | 6, 10 | 4²3⁴ | 2 | yes | 4 | 16 | 2-sum |
| V7_0 | 7, 11 | 4·3⁶ | 3 | yes | 3 | 2 | dual: 6 vertices, h=6 |
| V7_1 | 7, 11 | 4·3⁶ | 2 | yes | 4 | 8 | 2-sums of K4 |
| V7_2 | 7, 11 | 4·3⁶ | 3 | yes | 2 | 4 | dual: 6 vertices, h=6 |
| V7_3 | 7, 11 | 4·3⁶ | 3 | no | 1 | 4 | |
| V8_0 | 8, 12 | 3⁸ | 2 | yes | 4 | 16 | 2-sums of K4 |
| V8_1 | 8, 12 | 3⁸ | 3 | yes | 2 | 4 | dual: 6 vertices, h=7 |
| V8_2 | 8, 12 | 3⁸ | 3 | no | 1 | 12 | |
| V8_3 | 8, 12 | 3⁸ | 3 | yes | 0 | 48 | cube |
| V8_4 | 8, 12 | 3⁸ | 3 | no | 0 | 16 | Wagner (Möbius ladder M₈) |

- Edge lists: `skel5.json`.
- |V|=6 with max deg ≤ 4: these are decompletions of C₇(1,2) (giving Z5) and of K₃,₄+M (giving K33+e and K4⊕K4) [drv].

### 1.2 Occurrence
The tables below count graphs whose series–parallel core is the given five-loop core [chk `core_occurrence.log`, `census_cores.log`, `enum_d3.log`].

**D=2, n=5** (|I|=11): 70 graphs.

| core | SP | K4 | W4 | prism | K33 | K5−e | W5 | Z5 | K33+e | K4⊕K4 | V7_0 | V7_1 | V7_2 | V7_3 | V8_* |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| graphs | 12 | 20 | 13 | 2 | 1 | 5 | 2 | 6 | 3 | 2 | 1 | 1 | 1 | 1 | 0 |

- V8_* cores cannot occur: the graph has |V(Γ)|=7 [drv].

**D=3, n=5** (|I|=16): all 14 cores occur.

| | K5−e | W5 | K33+e | Z5 | K4⊕K4 | V7_0 | V7_1 | V7_2 | V7_3 | V8_0 | V8_1 | V8_2 | cube | Wagner |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| graphs | 111 | 120 | 157 | 681 | 104 | 456 | 85 | 325 | 416 | 11 | 136 | 74 | 38 | 93 |

**bulk3(3,5):** 410 graphs on 16 cubic skeletons. Classes use bivalent sources.

| core | count (class) |
|---|---|
| SP | 3 |
| K4 | 11 (2 B, 9 R₁) |
| prism | 32 (6 R₁, 26 R₂) |
| K33 | 12 (10 R₂, 2 R₃) |
| V8_0 | 11 (1 B, 3 R₁, 7 R₂) |
| V8_1 | 136 (7 R₂, 109 R₃, 20 R₄) |
| V8_2 | 74 (28 R₃, 42 R₄, 4 R₅) |
| cube | 38 (21 R₃, 15 R₄, 2 R₅) |
| Wagner | 93 (I: no S, P, Y0, Y1 move) |

**Twistor Q₅:** b₁(Γ∖O) = n+1−deg O ≤ 1 ⇒ deg O ∈ {5,6} [drv]. Classes use bivalent ∪ {O} sources.

| set | count (class) |
|---|---|
| tree | 49 (all SP) |
| unicyclic | 506: SP 36, K4 159 (29 B, 130 R₁), W4 191 (74 R₁, 117 R₂), W5 120 (3 R₁, 55 R₂, 62 R₃) |

**Enumeration checks** [chk `check_enum_d3.log`]:
- Connected loopless cubic multigraphs on 2, 4, 6, 8 vertices: 1, 2, 6, 20 [src OEIS A000421]; 2-connected on 8 vertices: 16.
- All 965 graphs pass the brute-force admissibility test.
- The same code at n=4 gives bulk3 50 and twistor 15+95, equal to the atlas.
- Twistor tree count 49 = cone(3,5) [src P/explore/six-d-Q3/README.md l.122].
- The V8 counts and W5 = 120 agree between `core_occurrence.log` and `enum_d3.log`.

## 2. Reductions and relations

Reductions (`per5.lower`, `msys5` hook):

$$\text{(2-sum)}\quad I_q(G)=\frac{\Gamma(q)}{\Gamma(s_A)\Gamma(s_B)}\,I_q(A+e_{uv};s_B)\,I_q(B+e_{uv};s_A),\qquad s_X=\textstyle\sum_Xa-q\,h(X),\ \ s_A+s_B=q$$

- [drv: two-terminal push-forward f_A = Γ(q) p(A + e(q−s_A)), which reproduces G₂q for a parallel pair and 1 for a series pair.]
- [chk `check_lower.log`: 16 vectors agree with the series–parallel rules; K4⊕K4 at q=2 gives 36ζ₃²; dimension relations hold for K4⊕K4 (q=2→3) and V7_1 (64 vectors, q=3→4).]

$$\text{(P)}\quad I_q(G^*;a)=I_q(G;q-a\circ\pi)$$

- [src: orchestrator, `../symmetry/d0_planar_duality.log`]
- K5−e = prism* [chk: dimension relations 45 vectors (q=3→4), 200 vectors (q=4→5), 0 mismatches].
- Hook: a reduced non-family graph with 2|V|−2 < |E| is replaced by its dual.
- W_n: `cores/wheel.py`; W5 at q=2 gives 70ζ₇ [chk].

New relations. Each is used only if every term in it is convergent [chk `check_rel5.log`, `check_completion.log`, `check_syz5.log`; tested on K4, W4, prism, K33, W5, K5−e, K4⊕K4 for q ≤ 7: (Y) 845, (X) 4713, (P) 29296, (C) 15638, (L) 3850 relations, 0 failures].

- **(Y)** Δ→Y at a triangle with c₁+c₂+c₃ = 2q:
  $$I_q(G;c)=I_q(G-\triangle+Y;a),\quad a_{sx_k}=q-c_{\rm opp}$$
  [drv: inverse of (S)]
- **(X)** position-space vertex IBP at vertex x, with x₁ adjacent to all other neighbours x_i (e₁ = xx₁, e_i = xx_i, f_i = x₁x_i):
  $$X_q(G;q-a)=\frac{I_q(G;a)}{\prod_e\Gamma(q-a_e)}$$
  $$(2a_1+\textstyle\sum_ia_i-(k-1)q)\,v(a)I(a)=\sum_i(q-a_i)\big[v(a+e_1-e_i)I(a+e_1-e_i)-v(a-e_i+f_i)I(a-e_i+f_i)\big],\qquad v(a)=\prod_e\Gamma(q-a_e)^{-1}$$
  - [drv: per-edge Fourier factor; 4^{qh−Σa} = 1; the π-powers cancel]
  - Used only if every term has 1 ≤ a_e ≤ q−1. Terms with a_e = 0 carry an IBP boundary term [chk: dropping them gives failures].
  - A first normalization Π Γ(q−a)/Γ(a)² was wrong at q ≥ 5 [chk, corrected].
- **(C)** conformal completion. Ĝ = G + ∞ with weights w_v = 2q − Σ_{e∋v}(q−a_e) ∈ [0, q−1]:
  $$\frac{I_q(G;a)}{\prod_e\Gamma(q-a_e)}=\frac{I_q(\hat G-u;\hat a)}{\prod_e\Gamma(q-\hat a_e)},\quad \hat a_{\infty v}=q-w_v$$
  [drv]
- **(L)** syzygies Σ_i a_i∂_iU = bU, with deg a = k. From Σ_i∂_i(a_i t^{β−1}U^{−q}):
  $$\sum_\gamma c_\gamma I_q(G;\gamma)+\sum_i\sum_\delta d_\delta\,I_q(G/i;\delta)=0$$
  - Face terms t_i = 0 contribute −∫ F_i. Used only if every β+μ−e_i (μ ∈ supp a_i) is convergent [chk: without this condition some relations fail].
  - Minimal generators [chk `syz_cache/`]:

| core | generators at k=1, 2 (3) |
|---|---|
| Z5 | 3, 17 (28) |
| V7_0 | 7, 14 |
| V7_2 | 3, 39 |
| V7_3 | 3, 29 |
| V8_1 | 11, 11 |
| V8_2 | 9, 18 |
| cube | 1, 74 |
| Wagner | 1, 68 |

## 3. Systems and tables

| system | unknowns | rank | result | log |
|---|---|---|---|---|
| Z5 q=3, old relations (T, I, S) | 74 | 25 | [src `cores/zigzag_z5_q3.log`] | |
| Z5 q=3, + K5−e boundary | 74 | 44 | | `build_q3.log` |
| Z5 q=4: + K5−e / + Y, X, P, C / + L (k≤2); k=3 adds nothing | 1956 | 1813 / 1927 / 1952 | | `build_q4.log`, `build5L_q4*.log` |
| Z5 q=3..4 + dimension links + L | 2030 | 2026 | all free columns at q=3 | `multi_34_L.log` |
| Z5 q=3..5 + L | 24698 | 24694 | q=5 determined given q=3 | `multi_345_L.log` |
| q=3: Z5, V7_0, V7_2, V7_3, V8_1..4; q=4: Z5 | 2274 | 2272 | | `joint_3all_4Z5.log` |
| same + input I₂(Z5;1¹⁰) = 441/8 ζ₇ | 2274 | 2273 | Z5 q3 74/74, q4 1956/1956; V7_0 64, V7_2 52, V8_1 22, V8_3 8 | `joint_3all_4Z5_z5q2.log` |
| V7_3, V8_2, V8_4, q=3..4 | 5579 | 5579 | complete | `joint_np_34.log` |
| no input; q=4 adds V7_3, V8_2, V8_4 | 7830 | 7828 | deficiency 2: the input is necessary | `joint_noinput_34.log` |
| K33+e q=3..4 (joint with Z5, no input) | | | 18/18 | `multi_34.log` |
| K33+e q=3..5 | | | 429/429 at q=4 | `multi_345_ZK.log` |
| K33+e q=5, build-up | 4662 | 4662 | | `buildup_K33e_q5.log` |
| Z5 q=5, build-up from the q=4 table | 22668 | 22668 | | `buildup_Z5_q5.log` |
| V7_0, V7_2, V8_1, V8_3 q=4, build-up | 8190 | 8074 | all 8074 family orbits | `buildup_V_q4.log` |

- Source of the input: [src arXiv:1208.1890 (zigzag theorem); `cores/REPORT.md` §2]. P(Z₅) = 4·8!/(5!4!)·(63/64)·ζ₇ = 441/8 ζ₇.

Tables in `tables/`, with the labels that appear [chk `table_summary.log`]:

| core | q (orbits) | labels |
|---|---|---|
| Z5 | 3 (74), 4 (1956), 5 (22668) | 1, ζ₃, ζ₅, ζ₇, ζ₃² |
| K33+e | 3 (18), 4 (429), 5 (4662) | 1, ζ₃, ζ₅, ζ₃² (no ζ₇) |
| V7_0, V7_2 | 3 (64, 52), 4 (3147, 2301) | 1, ζ₃, ζ₅, ζ₇, ζ₃² |
| V7_3 | 3 (65), 4 (2862) | same |
| V8_1, V8_3 | 3 (22, 8), 4 (2161, 465) | same |
| V8_2, V8_4 | 3 (13, 19), 4 (1157, 1398) | same |
| K5−e | every q with a prism table (q ≤ 8 with `assembly/tables`) | 1, ζ₃, ζ₅ |
| W5 | any q (`cores/wheel.py`) | 1, ζ₃, ζ₅, ζ₇ |
| K4⊕K4, V7_1, V8_0 | any q (products of K4 values) | 1, ζ₃, ζ₃² |

Validation:
- **q=2** (φ⁴ five-loop primitives): W5 = 70ζ₇ [chk; src Broadhurst wheel formula, `cores/REPORT.md` §2]; K4⊕K4 = 36ζ₃² [chk]; K33+e = 36ζ₃², obtained from the q=3 table by the dimension relation, not used as input [chk; src product identity P₃² = (6ζ₃)², Schnetz census arXiv:0801.2856]; Z5 = 441/8 ζ₇ (input). K5−e, V7, V8 have no convergent q=2 vector [drv].
- **Independent exact** [chk `compare_graphical.log`]: the graphical-function tables agree on every entry. Z5 q=3 74/74, K33+e q=3/4 6/111, K5−e q=3/4 6/106, K4⊕K4 3/56, V7_0 64, V7_1 2, V7_2 52, V7_3 19, V8_1 4, V8_2 3, cube 1, Wagner 7, W5 16; 0 differences. This confirms the Z5 tables, which used the input.
- **Internal:** the 65 duplicate level-3 V7_3 unknowns inside the q=4 system equal the q=3 table [chk `check_np_dup.log`].
- **Float:** tropical sampling [method src arXiv:2008.12310]. About 60 values of K33+e, Z5, V7_*, V8_* at q=3, 4, 5, and K5 q=3, agree within 2.7σ (≤1σ when rerun with 3·10⁶ samples) [chk `check_tropmc.log`, `check_k5_tropmc.log`]. Naive Gamma-sampling Monte Carlo fails because of heavy tails [chk, run only].
- **Negative result:** `rec5.py`, recursion without a table (T rules with Y/S/C as terminals), reaches 194/400 Z5 q=4 orbits (0 wrong). It does not reach I₇(Z5;4,3,3,4,3,4,3,3,4,4) [chk `check_z5_q7_graphical.log`].

## 4. Kernels

### 4.1 D=2, n=5 (E6 row of the atlas)

Method: `fk5`.
- Whole-graph modular Newton interpolation of P = U^{n+H}D on 8,209,410 points, primes p < 2³¹, about 52 s per prime, 3 primes.
- The exact series-skeleton numerator is reconstructed by CRT and saved in `cnum/<i>.json`; it is then contracted with `per5` periods at q=7.

Sources: generic; λ_v ∈ [−3,3]² at every vertex, `random.Random(1000+i)`, root = last vertex.

Checks [chk `check_fk5.log`, `check_B_rational.log`]:
- C² four-loop FORCER fixtures (prism, K33): match.
- SP graph #5 against `htk.kernel_exact`: 1432537/86400, match.
- B-class source patterns (bivalent vertices plus extra vertices chosen so that S, P, Y0 reduce the graph; 3 sources) must give rational kernels. #32 (K5−e) = 7533/128000; #45 (W5) = −9/1600; #42, #64, #69, #56, #63 = 0. All rational, with every ζ part cancelling.

Coverage: 58/70.

| # | core | gen class | K(λ;0) |
|---|---|---|---|
| 0 | SP | SP | −12399/1600 |
| 1 | SP | SP | 166391/800 |
| 2 | SP | SP | 77557/43200 |
| 3 | SP | SP | −259/86400 |
| 4 | SP | SP | 110079/1600 |
| 5 | SP | SP | −19/19200 |
| 6 | SP | SP | 16251323/172800 |
| 11 | SP | SP | 0 |
| 13 | SP | SP | −1164367/19200 |
| 18 | SP | SP | −2930479/108000 |
| 24 | SP | SP | −55203/8000 |
| 38 | SP | SP | 3055783/1152000 |
| 7 | K4 | B | −56737453/403200 + 8819/840 ζ₃ |
| 8 | K4 | B | −70126741/8064000 + 531983/84000 ζ₃ |
| 9 | K4 | B | −110237/2800 + 5478/175 ζ₃ |
| 10 | K4 | R₁ | 335738861/16128000 − 220149/14000 ζ₃ |
| 12 | K4 | B | 1716551/604800 − 841/280 ζ₃ |
| 14 | K4 | B | −72515791/2268000 − 105961/14000 ζ₃ |
| 15 | K4 | R₁ | −94030751/1512000 + 8651557/168000 ζ₃ |
| 16 | K4 | R₁ | −65396131/72576000 + 5981/252000 ζ₃ |
| 17 | K4 | R₁ | −11813/10125 + 352/375 ζ₃ |
| 20 | K4 | B | 67396021/6912000 − 125659/24000 ζ₃ |
| 22 | K4 | B | −78398591/6804000 + 8280719/756000 ζ₃ |
| 25 | K4 | B | −8007232327/435456000 + 20678633/1512000 ζ₃ |
| 27 | K4 | B | −16929779/403200 + 22151/700 ζ₃ |
| 28 | K4 | B | 230171/201600 − 683/700 ζ₃ |
| 30 | K4 | B | 44003101/107520 − 102863/280 ζ₃ |
| 34 | K4 | B | −23572161/89600 + 1386423/28000 ζ₃ |
| 40 | K4 | R₁ | 635878391/1814400 + 906113/12600 ζ₃ |
| 43 | K4 | B | 16080847/10368000 − 1633/1500 ζ₃ |
| 46 | K4 | B | 4608673/3628800 − 4327/9450 ζ₃ |
| 60 | K4 | B | −541308529/7257600 − 61861/420 ζ₃ |
| 19 | W4 | B | 131168919487/82944000 + 7023740683/1814400 ζ₃ − 6066713/1008 ζ₅ |
| 21 | W4 | B | 147927523/7838208000 + 98/10125 ζ₃ − 2491/75600 ζ₅ |
| 23 | W4 | B | 0 (numerator ≡ 0) |
| 26 | W4 | B | 11952969809/2612736000 + 1429277/144000 ζ₃ − 8509/560 ζ₅ |
| 29 | W4 | B | 0 |
| 31 | W4 | B | 17646715589/163296000 + 4724612341/18144000 ζ₃ − 683479/1680 ζ₅ |
| 33 | W4 | B | −852042053/2551500 − 13947287/141750 ζ₃ + 215842/525 ζ₅ |
| 35 | W4 | R₁ | 8257933456397/3919104000 + 239016887/43200 ζ₃ − 425664751/50400 ζ₅ |
| 36 | W4 | B | 6164289827/8164800 + 283564661/226800 ζ₃ − 666068/315 ζ₅ |
| 41 | W4 | B | 96592769953/9331200 + 4139570081/151200 ζ₃ − 584165/14 ζ₅ |
| 47 | W4 | B | −287199343/1306368000 + 86417/504000 ζ₃ |
| 52 | W4 | B | −13076087/537600 + 66447/4480 ζ₃ |
| 59 | W4 | B | 1851583523/9676800 + 59471/1680 ζ₃ − 1569/16 ζ₅ |
| 56 | prism | B | 45148139719/139968000 + 758724559/850500 ζ₃ − 25442803/18900 ζ₅ |
| 62 | prism | B | 485080099013/559872000 + 891868367/504000 ζ₃ − 29159203/10080 ζ₅ |
| 63 | K33 | B | −212562375377/1567641600 − 2113051901/5443200 ζ₃ + 731669/1260 ζ₅ |
| 32 | K5−e | B | −53846263/4032000 − 468140747/13608000 ζ₃ + 1002851/18900 ζ₅ |
| 37 | K5−e | R₂ | 1109389163/4838400 + 417563047/907200 ζ₃ − 468331/630 ζ₅ |
| 39 | K5−e | B | −6905295080041/326592000 − 763703887199/13608000 ζ₃ + 107399903/1260 ζ₅ |
| 42 | K5−e | B | 2244623/126000 − 130273/15750 ζ₃ |
| 44 | K5−e | R₂ | −377873561/3024000 − 74912161/252000 ζ₃ + 45862/105 ζ₅ |
| 45 | W5 | B | 4726885375291/7838208000 + 8595128581/3024000 ζ₃ − 56442787/5040 ζ₅ + 151389/20 ζ₇ |
| 48 | W5 | B | 17068519090957/391910400 + 48429455113/453600 ζ₃ − 435653/7 ζ₅ − 319415/3 ζ₇ |
| 55 | K4⊕K4 | B | −9313/1400 + 717/175 ζ₃ |
| 64 | K4⊕K4 | B | 629779393/31752000 − 14180749/2646000 ζ₃ − 76297/6125 ζ₃² |
| 69 | V7_1 | B | 180037715/48771072 − 850867/88200 ζ₃ + 8907/4900 ζ₃² |

- Full table with core-vector counts and timings: `d2_kernels.md` (110–360 s per graph).
- Period field [chk]: K4 → Q+Qζ₃; W4, prism, K33, K5−e → Q+Qζ₃+Qζ₅; W5 → +ζ₇; K4⊕K4, V7_1 → Q+Qζ₃+Qζ₃².

Still lacking a route [open], 12 graphs:

| graphs | core |
|---|---|
| 49, 50, 51, 53, 54, 58 | Z5 |
| 57, 61, 65 | K33+e |
| 66 | V7_0 |
| 68 | V7_2 |
| 67 | V7_3 |

- Skeleton numerators are saved for 49, 57, 66: 175,671, 208,098 and 921,867 core vectors at q=7.
- Full tables at q=7 are too large: Z5 has about 2.2·10⁶ orbits (4,395,456 compositions), V7 cores about 2.1·10⁷ compositions.
- The graphical engine takes 130–210 s per Z5 orbit at q=7, i.e. about 9·10⁴ orbits × 170 s ≈ 4·10³ h for graph 49.

### 4.2 D=3, n=5
- **SP graphs** (`htk.kernel_exact`; sources in [−3,3]³, seed 1000·idx+7): bulk3 3/3 (e.g. 1/75600000, −1419857/2419200000, 1/774144); twistor 77/85 exact (4 of them zero), 0.1–18 s each. 8 hit the 300 s limit: tree #31, 33, 44, 46, 47, 48 and unicyclic #120, 121 [chk `d3_sp_*.jsonl`]. Period field Q [src htk spred].
- **Non-SP:** [open]. The whole-graph lower set has 4,855,374,080 points (compare 8,209,410 at D=2, n=5 and 4,384,432 at D=3, n=4) [chk]. The skeleton-level numerator would need 1,829,850 points for 10 core edges and 34,219,650 for 12. The assembly route (`assembly/numerator.py`) also interpolates the whole graph.
- **q=8 periods available:** K4, W4, W5 (wheel), prism/K33 (`assembly/tables`), K5−e (dual prism), V8_0 (2-sum). This covers all 470 non-SP twistor graphs and 66 non-SP bulk3 graphs. V8_1–V8_4 at q=8 [open]; 341 bulk3 graphs need them.

## 5. Stretch: D=2, n=6, the 32 graphs without a circuit chart [chk `k5_deltay.log`]
- Series–parallel cores are six-loop:

| (\|V\|,\|E\|) | graphs | remark |
|---|---|---|
| (5,10) | 3 | K5 |
| (6,11) | 10 | |
| (7,12) | 13 | |
| (8,13) | 6 | |

- Δ–Y on the triangles of Γ gives the five-loop cores K33+e, V7_3, V8_2, V8_4 (all non-planar). Two graphs have no triangle.
- Period level: K5 → K33+e through (Y) at a unique triangle.
  $$I_3(K_5;1,1,2^8)=\tfrac53\zeta_3-\tfrac53\zeta_5\ \ (\text{1's adjacent}),\qquad I_3(K_5;\ 1\text{'s disjoint})=-2\zeta_3+\tfrac52\zeta_5$$
  The value is the same from all 5 (resp. 4) unique triangles, 4 (C) relations hold, and tropical sampling agrees at 0.7σ and 0.0σ [chk].
- n=6 kernels need q=8 K5-type periods [open].

## 6. Open / failed
- Five-loop tables at q ≥ 6 (q ≥ 5 for V7/V8): too large in exact Python.
- Relations alone leave 2 free q=3 parameters for Z5; the q=2 literature value fixes one, and the q=4 V7_3, V8_2, V8_4 levels fix the other [chk].
- `rec5` recursion without a table is insufficient.
- D=2: the 12 graphs listed in §4.1. D=3: all non-SP kernels (numerator size) and V8_1–V8_4 at q=8.
- Timings: Z5 q=3..4 joint 494 s; non-planar q=3..4 765 s; K33+e q=5 273 s; Z5 q=5 308 s; V q=4 951 s; peak memory ≤ 1.6 GB.
