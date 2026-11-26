# 3d N=2 HT twist: all-loop supercharge Q = Σ_n ħ^n Q_n (item B3)

Tags: [chk] exact check (script and log in this folder), [drv] derived, [src] source, [open].
Folder: `explore/physics/3dN2-supercharge/`.
Code:
- `ht3d.py`: letters, Q_0, and the general Q_2 graph engine.
- `kernels.py`: htk Taylor tables of the six labelled bitriangles, giving `kernels_4.pkl` and `kernels.log`.
- `fanengine.py`: fan(n) and tadpole(n+1,0) at z = 0.
- `cohom.py` and `lift_general.py`: H(Q_0) and d_2.
- `theories.py`: the concrete theories.

Run everything with `PYTHONPYCACHEPREFIX=<scratchpad>/pycache_3d`.

## 0. Letters, interaction, Q_0

Letters (0-form jets):
- ∂^d c^a odd
- ∂^d B_a even
- ∂^d φ^i even
- ∂^d ψ_i odd

Gradings:
- J = N_B + N_ψ + D + Σ_k (r_k/2)(N_φk − N_ψk)
- R = N_c + Σ_k[r_k N_φk + (1−r_k)N_ψk]
- af = N_B + N_ψ

Conventions: A = c^aT_a, t_a = −iT_a, [t_a,t_b] = f_ab^c t_c, ρ_a = −iT_a. Then
$$I=\tfrac12f_{bc}{}^aB_ac^bc^c+\psi_i\,c^a\rho_a{}^i{}_j\phi^j-W(\phi)+\tfrac\kappa2c^a\partial c^a$$
[drv from src CDG (3.32), (3.34)]

$$\eta(c^a,B_b)=\delta^a_b,\quad \eta(\phi^i,\psi_j)=\delta^i_j\quad(\eta\ \text{symmetric})$$
[chk `test_q0.log`: with these signs in I, Q_0² = 0 requires η(c,B)η(φ,ψ) = +1; the choice η(c,B) = −1 fails on 482 of 575 monomials]

$$Q_0c^a=\tfrac12f_{bc}{}^ac^bc^c,\quad Q_0\phi=c^a\rho_a\phi,\quad Q_0\psi_i=-\partial_iW-c^a\psi_j\rho_a{}^j{}_i,\quad Q_0B_a=\kappa\partial c^a+f_{ab}{}^cc^bB_c-\psi_j\rho_a{}^j{}_i\phi^i,\quad[Q_0,\partial]=0$$
[chk `test_q0.log`, from the 1-edge graph with K = e^{−λz}; this is CDG (3.34) restricted to 0-forms up to ψ → −ψ and the normalizations of μ and k]

Q_0² = 0 holds on all monomials of degree ≤ 3 (degree ≤ 2 for three adjoints) [chk `test_q0.log`] in:
- SU(2)+adj with W = ½φ·φ, κ = 1/3
- U(1) with charges ±1, W = φ₁φ₂, κ = 2
- SU(2)+3 adj with W = εXYZ

## 1. Normalization N_n

**1.1 Propagator.** Higher.tex (3.15): 𝒫 = π^{−1/2}e^{−s²−zY}dY ds, with Y = z̄/τ and s = u/τ^{1/2}.

$$d\mathcal P=0\ (\mathrm{mod}\ dz),\qquad P_{\rm GKW}=-\!\int_0^\infty\!\iota_{\partial_\tau}\mathcal P=\frac{2\bar z\,du-u\,d\bar z}{4r^3},\quad r^2=u^2+|z|^2$$
[chk `normalization.log`]

$$d'P_{\rm GKW}=0\ (x\ne0),\qquad\int_{B^3\ni0}\mathrm{dVol}\wedge d'P_{\rm GKW}=1\quad(\text{orientation of (3.2)})$$
[chk]

$$P_{\rm GKW}\wedge\mathrm{dVol}=\frac1{8\pi i}\frac{(2\bar z\,dt-t\,d\bar z)\,dz}{|x|^3}$$
[chk]

CDG (5.7) as printed is $\frac{3}{16\pi i}\frac{(\bar z\,dt-t\,d\bar z)dz}{|x|^3}$ [src CDG (5.7), ar5iv]. It is not d'-closed:
$$d'\Big[\frac{\bar z\,dt-t\,d\bar z}{|x|^3}\Big]=\frac{2t^2-|z|^2}{2|x|^5}\,dt\,d\bar z\neq0$$
This also fails for |x|² = t² + 2|z|² and t² + 4|z|². Closure forces the numerator 2z̄ dt − t dz̄ [chk]. I used P_GKW.

**1.2 Kernel conversion.**
$$I^{\rm GKW}_\Gamma=\epsilon_n(4\pi)^{-n/2}K^{\rm htk}_\Gamma$$
- Topological factor: α_Γ = (−1)^n(4π)^{−n/2}ρ_n [src `circuit-general/topological/REPORT.md` §1, checked there].
- Holomorphic factor: ∫_{iℝ}(dx/2πi)e^{xw} = δ(w), so ρ^GKW = ±U^{−1}∏ψ, with no factors of π [drv].
- ε_n is a single sign for each n [drv]; its value is [open].

GKW (3.59), edge order (t12,t23,t24,t31,t43):
$$\omega^{(3.59)}=+2\pi^3\cdot(4\pi)^{-1}\Omega^{\rm htk}=\tfrac{\pi^2}{2}\Omega^{\rm htk}$$
[chk `normalization.log`, 3 points]. So (3.59) disagrees with (3.2), (3.15) and (3.52) by a factor 2π³.

**1.3 N_n.** With path-integral weight e^{S_CDG/ħ}, S_CDG = ±2πi S_GKW, so the loop parameter is g = ±ħ/(2πi) [drv]. The interacting bracket is {O}_I = Σ(1/n!){I^n,O} [src Higher.tex §4.5]. Then

$$Q_n=N_n\hat Q_n,\quad\hat Q_n=\frac1{(n+1)!}\sum_{\Gamma\ \text{labelled 2-Laman}}\frac{K^{\rm htk}_\Gamma}{\pi^{n/2}}W_\Gamma\ (\text{rational}),\quad N_n=\epsilon_n\Big(\frac{\pm\hbar}{4\pi i}\Big)^n$$

$$N_2=-\epsilon_2\frac{\hbar^2}{16\pi^2},\qquad N_4=\epsilon_4\frac{\hbar^4}{256\pi^4},\qquad Q_{\rm odd}=0$$
[drv; src ht-d2 §2.1]

Checks:
- n = 0: the 1-edge graph gives Q_0 = {I,·} (§0) [chk].
- ε_4 = +1 is fixed by nilpotency at ħ⁴ (§3.3) [chk], so N_4/N_2² = +1.
- With weight e^{iS/ħ} instead, N_2 changes sign [drv].
- When W = κ = 0, rescaling (B,ψ) → α(B,ψ) maps Q_2 → α²Q_2 [drv]. Only the size of ħ² relative to W and κ is physical.

## 2. Graph rule and colour/Wick weights

**2.1 Rule** (Higher.tex §4.1):
1. Move the contracted pairs (x_e, y_e) to the front of I_0⋯I_nO in edge order; record the Koszul sign.
2. Replace each pair by η(x_e,y_e). Survivors stay in their original order.
3. Apply K^htk with the same edge order and orientation:
   - a letter with d derivatives at e(0) gives (−y_e)^d;
   - a letter with d derivatives at e(1) gives (+y_e)^d;
   - λ_v → ∂ acting on the survivors at vertex v.
4. Moments: M_p = (−1)^{|p|}∏_e p_e![z^p]K.

The result does not depend on edge order or orientation [chk `check_conventions.log` (a)]. Each bitriangle was randomly permuted and reversed and its kernel recomputed by htk:
- SU(2)+adj: 1102 monomials (372 with Q̂_2 ≠ 0), 0 mismatches.
- U(1) with charges 1,2,−3: 343 monomials (66 nonzero), 0 mismatches.

**2.2 The n = 2 graph.** The bitriangle K₄ − e is the only 2-Laman graph on 4 vertices [src Higher, bitriangle figure]. It appears with two rootings:
- O at a degree-3 vertex: fan(2).
- O at a degree-2 vertex: tadpole(3,0).

Kernels [chk `kernels.log`]:
- fan(2), edges (01,03,12,13,23): K_0 = (π/3)l_0l_2 for any l_1;
  K_1 = −(π/30)l_0l_2[l_0(4z01+6z03+z12+3z13+z23) + 2l_1(−z01+z03+z12+3z13+z23) + l_2(−z01+z03−4z12+3z13+6z23)].
- tadpole(3,0), edges (01,02,12,13,23): K_0 = (π/3)l_0(l_0+l_1+l_2).

All orders ≤ 4 are in `kernels_4.pkl`.

Graphs containing a CS vertex contribute 0 [chk (b): Q̂_2 is identical with and without κ for U(1) with charges 1,2,−3 (220 monomials) and SU(2)+adj (1140)].

**2.3 Weights.** The sum over labelled graphs equals ½ × the patterns on the canonical fan(n) or tadpole(n+1,0) [drv; chk: the pattern engine equals the general engine on all undifferentiated monomials of degree ≤ 4].

$$\hat Q_2\big|_{z=0}=\tfrac16\sum\varepsilon_KI_{DA_1R}\eta^{RR'}I_{R'A_2S}\eta^{SS'}I_{S'A_3E}\,\partial x^D\partial x^E\,\partial_{A_3}\partial_{A_2}\partial_{A_1}\;+\;\tfrac13d^V_{abc}\,\partial^2c^c\,\partial_{B_a}\partial_{B_b}$$

Tadpole term, with d^V_abc = tr_V ρ_(aρ_bρ_c):
- The two cycle orientations add: Q̂_2(B_aB_b)|_tad = ⅓[tr(ρ_aρ_bρ_c) + tr(ρ_bρ_aρ_c)]∂²c^c [chk (d): GL(2) with 1F, 1F+1AF and 2F, 10/10 each].
- Gauge loops and mixed loops give 0 [chk (d)].
- The term is supported only on B_aB_b [chk `check_tadpole_support.log`].
- It vanishes for real V, e.g. adjoint matter or V⊕V* [chk].
- For U(1), per unit Σq³: B∂B → ⅓∂³c; ∂B∂B → (13/105)∂⁴c; B∂²B → (22/105)∂⁴c [chk (c)].
- In GL(2)+1F the fan part and the tadpole part are each not Q_0-closed (70 failures each), but their sum is [chk `check_split.log`].

U(1) with one chiral of charge q, complete at z = 0 [chk `tables.log`]:
$$\hat Q_2|_{z=0}=\frac{q^3}{3}\big[\partial^2c+\partial c\,\partial\phi\,\partial_\phi-\partial c\,\partial\psi\,\partial_\psi-\partial\phi\,\partial\psi\,\partial_B\big]\partial_B^2$$

Examples:
- Q̂_2(B²∂ψ) = −(1/15)∂c∂²ψ + (4/15)∂²c∂ψ.
- SU(2)+adj: Q̂_2(B_0B_0c^0) = −2∂c^1∂c^2 [chk].

## 3. Nilpotency

**3.1 Order ħ²: {Q_0,Q_2} = 0** [chk `test_q2_*.log`, `check_conventions.log` (d); U(1) 1,2,−3 rerun by the orchestrator 2026-11-23: 1352 monomials, 0 failures]

| theory | W | κ | monomials (deg, D) | Q̂_2 ≠ 0 | both Q_0Q̂_2 and Q̂_2Q_0 ≠ 0 | tadpole ≠ 0 | failures |
|---|---|---|---|---|---|---|---|
| SU(2)+adj | 0 | 0 | 1288 (≤4, 0) | 543 | 540 | 0 | 0 |
| SU(2)+adj | ½φ·φ | 1/3 | 1396 (≤3, ≤1) | 309 | 309 | 0 | 0 |
| SU(2)+2adj | φ₁·φ₂ | 2/5 | 1158 (≤3, 0) | 126 | 126 | 0 | 0 |
| SU(2)+3adj | εXYZ | 0 | 2624 (≤3, 0) | 261 | 261 | 0 | 0 |
| SU(2) pure | – | 3/7 | 506 (≤4, ≤1) | 345 | 342 | 0 | 0 |
| U(1) ±1 | φ₁φ₂ | 2 | 506 (≤4, ≤1) | 97 | 85 | 0 | 0 |
| U(1) 1,2,−3 | φ₁φ₂φ₃ | 1/2 | 1352 (≤4, ≤1) | 276 | 246 | 145 | 0 |
| GL(2) 1F / 1F+1AF / 2F | 0 | 0 | 364 / 816 / 816 (≤3, 0) | 109 / 136 / 159 | – | ≠0 | 0 |

Negative control: doubling the order-1 moments breaks the identity on 27 of 824 monomials in U(1) 1,2,−3 [chk `check_split.log`].

**3.2 Q_4** (cubic theories, z = 0, κ = 0): Q̂_4 = fan(4) + tadpole(5,0) [src ht-d2 §3.3–3.4].
- fan(4): K = +(π²/180)(l_0l_4)² [chk `fan4_kernel.log`].

**3.3 Order ħ⁴.** Solve Q̂_2² + u{Q_0,Q̂_4^fan} + v{Q_0,T̂_{5,0}} = 0, with u = ε_4 and v = ε_4t_4:

| theory | monomials | Q̂_2² ≠ 0 | equations | solution |
|---|---|---|---|---|
| U(1), q = 1 | 60 | 11 | 19 | u = 1 ({Q_0,T̂} ≡ 0) |
| SU(2)+adj | 150 | 26 | 158 | u = 1 ({Q_0,T̂} ≡ 0) |
| GL(2)+1F | 40 (B⁴X) | 27 | 270 | u = 1, v = 1/630 |

[chk `test_q4_u1.log` (rerun by the orchestrator, identical), `test_q4_su2adj.log`, `test_q4_t4.log`]

$$K_{\rm tadpole(5,0)}(\lambda;0)=\frac{\pi^2}{630}\lambda^4$$
Here v/u = t_4 = 1/630 [chk: determined by nilpotency]. It agrees exactly with the direct evaluation [src `circuit-general/atlas/e3_tadpole50.log`] and with the earlier Monte Carlo value 0.00159(1) [src ht-d2 `mc_summary.log`].

## 4. Q_2 on H(Q_0)

d_2 = [Q̂_2]: H^R → H^{R+1} at fixed J and flavour. It is well defined because Q_1 = 0 and {Q_0,Q_2} = 0 [drv]. It does not change the index: R shifts by 1 while J and flavour are fixed [drv]. Ranks are computed exactly [chk `liftg_*.log`].

**U(1)_0 + one chiral** (r = ½; q = 1 and q = 2 give identical ranks). Dimensions are listed by R = 1,2,3,4 [chk `liftg_U1_24.log`, `liftg_U1k_24.log`, `liftg_U1q2_24.log`, `liftg_U1q2k_24.log`]:

| J | H(Q_0), κ = 0 | rank d_2 (1→2, 2→3, 3→4) | E_3 | H(Q_0), κ ≠ 0 |
|---|---|---|---|---|
| 1 | 1,1 | 0 | 1,1 | 1,1 |
| 2 | 2,2 | 0 | 2,2 | 2,2 |
| 3 | 3,4,1 | 1,1 | 2,2,0 | 2,2 |
| 4 | 4,6,2 | 2,2 | 2,2,0 | 2,2 |
| 5 | 5,10,5 | 3,3 | 2,4,2 | 2,4,2 |
| 6 | 6,14,9,1 | 4,5,1 | 2,5,3,0 | 2,5,3 |

$$Q\big(-\tfrac12B^2\partial c+B\,\partial\phi\,\psi\big)=\tfrac12N_2\,\partial c\,\partial^2c+O(\hbar^4)$$
At q = 1 the lifted classes are [chk]:
- (J,R) = (3,1) → (3,2): the line above;
- (3,2) → (3,3): c × the same;
- (4,1) → (4,2): ½B²∂²c + B∂B∂c − B∂φ∂ψ − B∂²φψ + ∂Bφ∂ψ ↦ −½∂c∂³c.

Prediction [chk J ≤ 6]: E_3 at κ = 0 equals the classical H(Q_0) at κ ≠ 0. Interpretation as a radiatively induced CS-type term [open]. q = 2 is consistent at k = 0 with the parity anomaly; q = 1 at k = 0 is not.

**No two-loop lifts** [chk] in:
- 3d N=4 SU(2) SYM (= SU(2)+adj, W = 0): J ≤ 3, about 7·10³ monomials per sector.
- N=4 SQED with N_f = 1: J ≤ 3.
- U(1)_{0,1} with charges (1,1) or (1,−1): J ≤ 3. Here H(κ=0) = H(κ=1).
- 3d N=8 SU(2) SYM: J ≤ 2, 74 sectors.
- U(1) with charges 1,2,−3: J ≤ 2.

Literature:
- CDG §1.4.1: whether a chiral analogue of formality holds for the bulk is unknown [src]. The table above is a non-formality instance.
- GKW: no non-renormalization theorem is expected at T = 1 [src Higher.tex l.199].
- Wess–Zumino models have no loop corrections [src CDG §5.1.1; ht-d2 §3.6].

## 5. Open / failed
- The overall sign ε_2 of Q_2 is [open]; it is also convention-dependent (e^{S/ħ} vs e^{iS/ħ}).
- CDG (5.7) as printed is not d'-closed. GKW (3.59) disagrees with GKW's own definitions by 2π³.
- Q_4 is only built at z = 0 and κ = 0 for cubic theories; derivative letters need the fan(4) and tadpole(5,0) z-moments [open]. Q_6 not attempted.
- d_2 not computed for non-abelian theories with d^V ≠ 0 (e.g. GL(2)+F) [open].
- The 2-loop anomaly (1/4!){I⁴} was not computed directly. {Q_0,Q_2} = 0 on every tested operator is consistent with it vanishing [drv].
