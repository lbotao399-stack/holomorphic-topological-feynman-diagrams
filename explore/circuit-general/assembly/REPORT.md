# Assembly: kernels with current insertions through the series reduction

Tags: [chk] exact check run (script → log in this folder), [drv] derived here, [src] source (location), [open].
Paths are relative to `P/explore/circuit-general/assembly/` (P = project root). Imported unchanged: `P/code/htk` and `../cores/{msys,gutil,vac,wheel,halfq,startri}.py`. Prism/K33 tables for q ≤ 7 are read from `../cores/tables/`. No existing file was edited. Nothing was committed.

## 0. Objects

$$K_\Gamma(\lambda;z)=\int_{\mathbb{RP}^{N-1}_{>0}}e^{-\sum_ey_e\cdot z_e}\,\Omega_\Gamma=K-\sum_{e,i}z^i_e\,m_{ei}+\sum_{(e,i)\le(f,j)}c\,z^i_ez^j_f\,m_{ei,fj}+O(z^3),\qquad c=\tfrac12\ (=),\ 1\ (<)$$
$$m_{ei}=\int D\,y^i_e,\qquad m_{ei,fj}=\int D\,y^i_e\,y^j_f$$
[src P/code/htk/CONVENTIONS.md §2]

$$P=U^{n+H}D=\sum_\alpha c_\alpha t^{\alpha-1},\quad \deg P=n^2-1,\quad \deg_{t_e}P\le n-1,\qquad K=\sum_\alpha c_\alpha I_{q_0}(\Gamma;\alpha),\quad q_0=n+H$$
[src ../cores/fastker.py; chk every prime: interpolant = direct sample at a random point]

## 1. Series reduction Γ → Γ_s

Chain c = skeleton edge p_c → q_c. Its edges e_{c,1..L_c} start at p_c. s_{c,j} = +1 iff e_{c,j} is oriented p_c → q_c. Interior vertices are u_{c,1..L_c−1}.

$$t_{c,j}=x_{c,j}R_c,\quad\textstyle\sum_jx_{c,j}=1,\quad\sigma_{c,k}=\sum_{j>k}x_{c,j},\qquad\prod_edt_e=\prod_cR_c^{L_c-1}dR_c\,d^{L_c-1}x_c,\qquad U_\Gamma(t)=U_s(R)$$
[drv: a spanning tree of Γ contains all of c (c ∈ T_s) or all but one of its edges (c ∉ T_s)]

$$y_{e_{c,j}}=s_{c,j}\big[Y_c(R;\lambda^s(\sigma))-\hat m_c+Q_{c,j}\big],\qquad \hat m_c=\sum_k\sigma_{c,k}\lambda_{u_{c,k}},\quad Q_{c,j}=\sum_{k<j}\lambda_{u_{c,k}}$$
$$\lambda^s_w(\sigma)=\lambda_w+\sum_{c,k}\lambda_{u_{c,k}}\big([p_c=w]\,\sigma_{c,k}+[q_c=w](1-\sigma_{c,k})\big)$$
[src CONVENTIONS §5, one-chain move; chk `Network.check_point`: exact against brute.currents for 20 graphs × 3 rational points × all (e,i), `check_pointwise.log`]

$$U_sY_c=\sum_{w\ne z}\Phi_{c,w}(R)\,\lambda^s_w,\qquad\Phi_{c,w}=\sum_{T\in\mathcal T(\Gamma_s)}f_T(c;w)\,R^{\kappa_T},\qquad\kappa_T=1_{E_s\setminus T}$$
Here f_T(c;w) is the flow on c in tree T from w to the sink z [src ../cores/numer.py tree_current_polys].

$$\Rightarrow\ y^i_e=s_e\Big[\sum_{T,\beta}A^{ei}_{T\beta}\frac{x^\beta R^{\kappa_T}}{U_s}+\sum_\beta B^{ei}_\beta x^\beta\Big],\quad x^\beta\in\{1\}\cup\{x_{c,j}\}_{j\ge2},\quad A^{ei}_{T\beta}=\sum_wf_T(c;w)[x^\beta]\lambda^{s,i}_w,\quad B^{ei}_\beta=[x^\beta](Q^i_{c,j}-\hat m^i_c)$$
This is the atlas REPORT §1.4 map $y_e=\sum_fA_{ef}(\xi)y_{{\rm core},f}+\sum_vB_{ev}(\xi)\lambda_v$ with ξ = x and core = Γ_s. Γ_s is not tight: at D=3, n=4, |E_s|−1 = 8 < nH = 12 [drv].

## 2. Push-forward and assembly

$$\int_{\Delta_c}\prod_{e\in c}x_e^{\alpha_e+\beta_e-1}dx=\frac{\prod_{e\in c}\Gamma(\alpha_e+\beta_e)}{\Gamma(\alpha_c+|\beta_c|)},\qquad\alpha_c=\sum_{e\in c}\alpha_e$$
[src fastker for β=0; drv]

$$PF_\beta(\gamma)=\sum_{\alpha:(\alpha_c)=\gamma}c_\alpha\prod_c\frac{\prod_{e\in c}\Gamma(\alpha_e+\beta_e)}{\Gamma(\alpha_c+|\beta_c|)},\qquad S^{(k)}_\beta(\kappa)=\sum_\gamma PF_\beta(\gamma)\,I_{q_0+k}(\Gamma_s;\gamma+\kappa)$$

$$K=S^{(0)}_1(0),\qquad m_{ei}=s_e\Big[\sum_{T,\beta}A^{ei}_{T\beta}S^{(1)}_\beta(\kappa_T)+\sum_\beta B^{ei}_\beta S^{(0)}_\beta(0)\Big]$$

$$m_{ei,fj}=s_es_f\Big[\sum_{\beta\beta'}B_\beta B'_{\beta'}S^{(0)}_{\beta\beta'}(0)+\sum_{T\beta\beta'}(A_{T\beta}B'_{\beta'}+B_\beta A'_{T\beta'})S^{(1)}_{\beta\beta'}(\kappa_T)+\sum_{TT'\beta\beta'}A_{T\beta}A'_{T'\beta'}S^{(2)}_{\beta\beta'}(\kappa_T+\kappa_{T'})\Big]$$
[drv; chk §5]

- Only the bare P is interpolated.
- Insertions and the levels q_0+1, q_0+2 are handled on Γ_s only (|E_s| ≤ 9).
- Dot products mod p use exact float64 BLAS on 16-bit limbs; one entry per call is re-checked in int64.

GL(H) covariance [drv: y → gy and D → det(g)^n D pointwise]:
$$K(g\lambda)=\det(g)^nK(\lambda),\quad m_e(g\lambda)=\det(g)^n\,g\,m_e(\lambda),\quad m_{ef}(g\lambda)=\det(g)^n(g\otimes g)m_{ef}(\lambda)$$
For bulk3(3,4): $K=c_\Gamma\det(\lambda_1,\lambda_2,\lambda_3)^4$ and $m_e=\det(\lambda_1,\lambda_2,\lambda_3)^4\sum_am_{e,a}\lambda_a$ [drv; chk check_gl.log]. These are computed at λ_a = e_a, with the legs ordered as in `run_batch.std_sources`.

## 3. Core periods

$$P(G;a)=\frac{\prod_e\Gamma(a_e)\prod_{\rm P}G_{2q}(s_a,s_b)}{\prod_{e'}\Gamma(a'_{e'})}P(G_f;a'),\quad G_{2q}(a,b)=\frac{\Gamma(q-a)\Gamma(q-b)\Gamma(a+b-q)}{\Gamma(a)\Gamma(b)\Gamma(2q-a-b)}$$

- Series step: s_a+s_b. Parallel step: s_a+s_b−q, convergent iff 0 < s_a, s_b < q and s_a+s_b > q [src ../cores/REPORT.md §0, msys.classify; drv vectorized].
- Final cores G_f:
  - banana: ∏Γ(q−s_i)/Γ(q);
  - W_n: wheel.py, any q;
  - prism, K33: orbit tables expanded by automorphisms, sorted-key lookup.
- [chk test_corelookup.log: 11 types (q = 6–9) × 300 vectors equal msys.classify, 0 mismatches]

| table | q | orbits | relations | rows | time | RSS | check |
|---|---|---|---|---|---|---|---|
| prism | 8 | 77322 | T, I, B, S | 998879 | 522 s | 2.60 GB | I_7 = Σ_C I_8(a+1_C) on all 22789 q=7 orbits: 0 failures |
| K33 | 8 | 19807 | D, I, B | 235343 | 748 s | 0.95 GB | same relation on 5988 orbits: 0 failures |
| prism, lean | 8 | 77322 | T, B, S | 793778 | 594 s | 1.12 GB | 0 differences from the full table |
| K33, lean | 8 | 19807 | D, B | 181451 | 701 s | 0.54 GB | 0 differences |
| prism, K33, lean | 6 | 5416, 1476 | | | 15 s, 40 s | | 0 differences |
| K33, lean | 9 | 56407 | D, B | 538934 | 2399 s | 1.47 GB | 8→9 on all 19807 orbits: 0 failures |
| prism, lean | 9 | 223549 | T, B, S | 2358333 | 3588 s | 3.17 GB | 8→9 on all 77322 orbits: 0 failures |

- Logs [chk]: tables_*_q8.log, check_tables_q8.log, lean_*_q8_check.log, tables_*_q9.log, check_tables_q9.log.
- "Lean" means the parametric IBP rows (I, up to |cotrees| terms each) are dropped. Full rank is verified at q = 6 and 8.
- W4 and K4 values come from wheel.py on demand (any q), memoized in-process. No persistent W4/K4 table file is written.

## 4. Bare numerator; cost of core-level interpolation

- Lower set L = {e ∈ [0,n−1]^{N−1} : |e| ≤ n²−1}.
- Samples P(1, 1+e) mod p via numba (`fastsample.py`): solve M ℓ = −C t y_0, then take det[C_re Y^i_e]_{e≠0}.
- Newton differences and conversion to monomials are done on L directly (L is downward closed) [src fastker, dense cube; drv lower-set form].

|L| = 4384432 (D=3, n=4) and 272034 (H=5, n=3). One prime takes 6–7 s for D3n4 (3–4 threads) and 1–2 s for H=5 [chk results logs].

Core-level interpolation would cost more than the whole-graph lower set [drv]:
- D=3, n=4: Γ_s has 153061 candidate monomials (bulk3 #12). Each value needs a 4-fold σ integral of degree ≤ 6 per σ, i.e. ≥ 4⁴ Gauss nodes. Total ≥ 3.9·10⁷ evaluations, versus |L| = 4.4·10⁶.
- H=5, n=3: 974 core monomials, 10 σ's, degree ≤ 4 per σ. Total ≥ 3¹⁰ evaluations per R, i.e. ≥ 5.7·10⁷, versus 2.7·10⁵.
- Hence the bare P is still interpolated on the whole graph. The reduction carries the insertions and the levels q > q_0.

Direct whole-graph route for one insertion (P·Uy, per-variable degree ≤ n, |e| ≤ n²−1+n):
- |L| = 44227105 points, about 95 s per (e,i) per prime [chk check_direct_n4.log].
- The assembly does the bare kernel plus all 39 (e,i), over 3–5 primes, in 35–85 s per graph.

## 5. Validation

| check | scope | result |
|---|---|---|
| decomposed y vs brute.currents (exact); P vs brute.density·U^{n+H} mod p; P·y vs brute | 20 graphs (bulk3 ×8, twistor Q4 ×6, bulk3(5,3) ×4, cone(5,3) ×2) × 3 rational points × all (e,i) | 0 failures [chk check_pointwise.log] |
| existing route `cores.kernel(…, insertions)`, n=3, D=3 | 6 K4 subdivisions: 6 bare, 30 first, 18 second moments | 0 mismatches [chk check_n3.log, test_small.log, test_small2.log; rerun by the orchestrator 2026-11-23: 54 compared, 0 mismatches] |
| independent whole-graph route at n=4 (interpolate P·Uy, q=8 prism table) | bulk3 #12: m_{0,0} and m_{12,2} mod p | equal [chk check_direct_n4.log] |
| bare kernels vs e3 (fastker, same seeds) | bulk3(3,4) 50/50 (via c·det⁴), twistor Q4 85/85, twistor Q3 23/23 | equal [chk compare_e3_bulk3.log, compare_e3_twistor4.log, compare_e3_twistor3.log] |
| htk spred (independent SP layer): bare + 3 first + 1 second moment per graph | bulk3(3,4) 3, twistor Q4 24, twistor Q3 9, bulk3(5,3) 24, cone(5,3) 21 | all equal [chk results/*.jsonl, field spred_match] |
| GL(3) | bulk3 #6, #12, #40 (order 1); #4 (order 2, 1718 values); #40 (order 2, K33 at q=9, 2577 values); twistor Q4 #30, #60 | 0 mismatches [chk check_gl.log, check_glcov_bulk3_o2.log, check_glcov_bulk3_40_o2.log, check_glcov_twistor4.log] |
| GL(5) | bulk3(5,3) #30, #45 (order 2, 6962 values each); cone(5,3) #3 (3481 values) | 0 mismatches [chk check_glcov_h5.log, check_glcov_cone5.log] |
| htk.oracle, lattice | θ(1,3,3), θ(2,2,3): bare, 2 first, 2 second moments each | \|Δ\| ≤ 0.13σ, σ/\|K\| ≈ 10⁻³ [chk check_oracle_theta.log] |
| htk.oracle, QMC (N−1 = 9) | K4[1,1,1,2,2,3]a: 5 values | \|Δ\| ≤ 1.33σ, σ/\|K\| ≈ 0.15–0.5 (weak) [chk check_oracle.log] |
| order-1 run vs order-2 run | bulk3 1799 first moments; twistor Q4 4056 first moments | 0 differences [chk check_o1_o2_consistency.log] |

## 6. Results

Files `results/*.jsonl` hold one record per graph: edges, sources, core, `bare`, and `moments` keyed "(e, i)" or "(e, i, f, j)", with labels '1', 'z3', 'z5' and timing in seconds.

| set | graphs | Γ_s → G_f | order | bare | nonzero m_{ei} | nonzero m_{ei,fj} | sources | files |
|---|---|---|---|---|---|---|---|---|
| bulk3(3,4) (hCS / KS on C³) | 50 | prism 28, K33 11, K4 8, banana 3 | 2 | 13 zero; 2 in Q; 5 in Q+Qζ3; 30 in Q+Qζ3+Qζ5 | 1799/1950 | 37883/39000 | λ_a = e_a | bulk3_o1_std, bulk3_o2_std{,_prismA,_prismB} |
| twistor Q4 | 110 | W4/K4 85, banana 25 | 2 | 6 zero; 30 in Q; 47 in Q+Qζ3; 27 in Q+Qζ3+Qζ5 | 4056/4290 | 81120/85800 | e3 seeds | twistor4_o1_s2_{0,1}, twistor4_o2{,_restA,_restB} |
| twistor Q3 | 23 | K4 14, banana 9 | 2 | 11 in Q; 12 in Q+Qζ3 | 690 | 10683 | e3 seeds | twistor3_o2 |
| bulk3(5,3) (BCOV / 10d on C⁵) | 83 | K4 56, banana 27 | 2 | 31 in Q; 52 in Q+Qζ3 | 6640/6640 | 268920/268920 | e3 seed rule, H=5 | bulk3_5_o2 |
| cone(5,3) | 28 | banana 28 | 2 | 28 in Q | 2240/2240 | 90720/90720 | e3 seed rule, H=5 | cone5_o2 |

- All values lie in Q+Qζ3+Qζ5 at n=4 and in Q+Qζ3 at n=3 [chk].
- bulk3(3,4): c_Γ = 0 for #1, 2, 3, 7, 8, 9, 10, 13, 24, 34, 42, 48, 49 [chk].
  - 12 of these have m_{ei} ≠ 0, so K = O(z).
  - #2 has all m_{ei} = 0 but 664/780 nonzero m_{ei,fj}, so K = O(z²).
- The e3 zeros #15 and #43 come from det(λ_1,λ_2,λ_3) = 0 at the e3 seeds (#10, #24 too); c_15, c_43 ≠ 0 [chk compare_e3_bulk3.log].
- twistor Q4 #17, #24, #43, #57, #88, #95: K, all m_{ei} and all m_{ei,fj} vanish [chk].

c_Γ (K = c_Γ det(λ_1,λ_2,λ_3)⁴), nonzero cases [chk results/bulk3_o1_std.jsonl]:

| # | core | 1 | ζ3 | ζ5 |
|---|---|---|---|---|
| 0 | banana | 1/186624000 | 0 | 0 |
| 4 | K4 | −90379/42326323200 | 29/16329600 | 0 |
| 5 | K4 | −1/62985600 | 0 | 0 |
| 6 | K4 | −69743/105815808000 | 11/20412000 | 0 |
| 11 | prism | −136903/9405849600 | −23941/587865600 | 5/81648 |
| 12 | prism | 516283/169305292800 | 763/83980800 | −11/816480 |
| 14 | prism | −90283/67722117120 | 1/20995200 | 1/816480 |
| 15 | prism | 1969714721/253957939200 | 1742063/83980800 | −25741/816480 |
| 16 | prism | −1771133/16124313600 | −37957/130636800 | 241/544320 |
| 17 | prism | −5474243/5079158784000 | 1777/1959552000 | 0 |
| 18 | prism | 966996061/677221171200 | 51269/13436928 | −18943/3265920 |
| 19 | prism | −1376521747/677221171200 | −304513/55987200 | 4499/544320 |
| 20 | prism | −17876711/1693052928000 | −8003/326592000 | 1/25920 |
| 21 | prism | 6626281/1058158080 | 58609/3499200 | −433/17010 |
| 22 | prism | 11629/16930529280 | −19/32659200 | 0 |
| 23 | prism | −5387821/1504935936 | −750487/78382080 | 7921/544320 |
| 25 | prism | 706474481/225740390400 | 2185063/261273600 | −173/13608 |
| 26 | prism | −765137/314928000 | −379/58320 | 4/405 |
| 27 | prism | 4642972919/1693052928000 | 14374987/1959552000 | −12137/1088640 |
| 28 | prism | 66991/17635968000 | −43/13608000 | 0 |
| 29 | prism | 230775907/60466176000 | 738469/72576000 | −13/840 |
| 30 | prism | −193542799/52907904000 | −7181059/734832000 | 1213/81648 |
| 31 | prism | −237259/12093235200 | −3083/48988800 | 5/54432 |
| 32 | prism | −197320301/253957939200 | −606637/293932800 | 10261/3265920 |
| 33 | prism | 3717317/56435097600 | 30901/167961600 | −113/408240 |
| 35 | prism | 4721359/3224862720 | 1533761/391910400 | −1619/272160 |
| 36 | prism | 2076679/587865600 | 76933/8164800 | −65/4536 |
| 37 | prism | 3476941/1045094400 | 7441/839808 | −5503/408240 |
| 38 | prism | 101978921/25395793920 | 2526137/235146240 | −1333/81648 |
| 39 | K33 | −339933143/94058496000 | −4742693/489888000 | 4003/272160 |
| 40 | K33 | −12694111093/3386105856000 | −57479/5740875 | 3551/233280 |
| 41 | K33 | −304937/734832000 | −44911/40824000 | 19/11340 |
| 43 | K33 | 17809481/22574039040 | 1239569/587865600 | −2617/816480 |
| 44 | K33 | 4417188833/1693052928000 | 13743059/1959552000 | −1159/108864 |
| 45 | K33 | −209968763/846526464000 | −7554973/11757312000 | 1607/1632960 |
| 46 | K33 | −153527/116121600 | −460847/130636800 | 73/13608 |
| 47 | K33 | 37777763/9405849600 | 311461/29393280 | −1649/102060 |

Example, bulk3 #12 with standard sources: m_{0,0} = −1624793/4389396480000 − (1988659/2263282560000) ζ3 + (173/125737920) ζ5. This value was checked against the direct whole-graph route [chk].

## 7. Timings

Wall time per graph on a shared machine, load average 7–48.

| run | per graph |
|---|---|
| bulk3(3,4), order 1 | 35–85 s (prism 54–85, K33 58–83, K4 41–50, banana 35–42); e3 fastker, bare only: 79–203 s |
| twistor Q4, order 1 | 34–140 s |
| bulk3(3,4), order 2 | banana 158–208 s, K4 197–427 s, K33 296–391 s, prism 322–490 s (q = 9 tables) |
| twistor Q4, order 2 | banana 46–115 s, W4/K4 52–256 s |
| bulk3(5,3), order 2 | K4 8–16 s, banana 4–7 s |
| cone(5,3), order 2 | 5–7 s |
| twistor Q3 and C³ three-loop, order 2 | 0.2–2 s (`cores.kernel`: 8 s per single moment) |
| direct whole-graph route, D3n4, one (e,i), one prime | ≈ 95 s, 3.4 GB |

## 8. Evaluator and API

Command line:
- `evaluate.py moments NV 'EDGES' 'SOURCES' [ORDER]`
- `evaluate.py set SET IDX [ORDER] [--std]`
- `evaluate.py period 'EDGES' 'POWERS' Q`
- Batches: `run_batch.py SET ORDER [K S] [--std] [--only=…] [--tag=…]`
- Tables: `lean_tables.py NAME Q --noI`

Python API:
- `assembly.Network(nv, edges, lam)` — chains, skeleton, the A/B map `insertion_coefficients(e,i)`, and `check_point`.
- `assembly.Assembly(net, order).run()` and `assembly.kernel_moments(...)`.
- `numerator.BareNumerator(G, lam_full).coefficients_mod(p)` and `numerator.LowerSet` (dense index map up to 3·10⁸ points, else searchsorted).
- `corelookup.CoreLookup(edges, q).prepare / evaluate` — final cores banana, W_n, and any family in `periods.FAMILIES` that has tables.
- `periods.table` and `periods.period`.

Coverage:
- T = 0, integer sources, ≥ 2 vertices of degree ≠ 2, and no chain that returns to its own start.
- Γ_s must reduce series-parallel to a banana, W_n, prism or K33.
- Level q_0 + order ≤ 9 for prism/K33; any q for W_n.
- |L| up to about 5·10⁷ points.

## 9. Not done / failed

- The bare P is still interpolated on the whole graph (12 or 15 variables). Core-level interpolation was costed in §4 and is more expensive, so it was not implemented [drv].
- The QMC oracle at N−1 = 9 gives only 15–50 % precision.
- htk spred references timed out (60–240 s cap) for 3 bulk3(5,3) graphs (#18, #19, #26) and 7 cone(5,3) graphs (#15, #17, #18, #19, #25, #26, #27). Assembled values exist for these, but there is no independent reference.
- [open] T = 1; Y1 (rational transfer) reductions; non-series final cores beyond prism/K33/W_n (e.g. five-loop cores); third order.
