# VERIFY: anomaly-structure

Scripts and logs: `explore/anomaly-structure/verify/`. None of them imports the explorer's code; `htk` is imported only for test-graph lists (`catalog.tight_graphs`) and for side-by-side engine values.

| file | content |
|---|---|
| `vk.py` | Kirchhoff solver (exact), positive cycle matrix from a chosen tree, $D(t)$ from $\Omega=U^{-H}\rho^T\wedge\prod_{i}\prod_r\psi^i_r$, $\rho_2=U^{-1/2}A_{12}$, $A=-\tfrac12 dM\,M^{-1}dM$; forward-mode derivatives $\partial y/\partial t_j$; Grundmann–Möller rules (rational nodes) |
| `sp_coords.py` | coordinates $t^P_j=x^P_j/w_P$ for $k$ parallel paths, $g(\xi)=D(t(\xi))\lvert\det[t;\partial_\xi t]\rvert$, exact GM integration |
| `lie_eval.py` | own bases ($\mathfrak{sl}_N$: $E_{ij}$, $E_{kk}-E_{NN}$), $f_{abc}=\kappa([T_a,T_b],T_c)$, $\kappa^{-1}$, colour weights from the index forms of (c3) and from the orientation rule as stated |
| `ce_eval.py` | cochains as functions of $(c_i,p_i)$ with $\sum p_i=0$; $d\omega(x_0..x_k)=\sum_{i<j}(-1)^{i+j}\omega([x_i,x_j],\dots)$; exact $\mathbb Q$-ranks |
| `v_a_polynomiality.py`, `v_kernels.py`, `v_S2_conservation.py`, `v_b1.py`, `v_lie_ce.py`, `v_ce_large.py`, `v_formalN.py`, `v_orientation.py` | checks below |

Exactness criterion used for cochains: value matrix $V$ (rows = random samples, columns = $d\eta_j$), $\omega$ a further column. $\operatorname{rank}_{\mathbb Q}[V|\omega]>\operatorname{rank}_{\mathbb Q}V$ is a certificate that $\omega\notin d(C^3_{\rm inv})$. Equality: $\eta$ solved exactly, $d\eta=\omega$ re-evaluated at fresh samples with coordinates and momenta in $[-30,30]$ (`v_ce_large.log`). $\eta$ spanned by $\mathrm{alt}(\mu\,\mathrm{tr}_\sigma)$, $\sigma\in S_3$ (first fundamental theorem for $GL_N$, Reynolds).

## (a1) confirmed

$$D_\psi(t)=\sigma_T^H\,\epsilon(n,H)\,D_{dy}(t),\qquad \epsilon(n,H)=(-1)^{nH+\binom n2\binom H2},\qquad \sigma_T=\varepsilon_T\det B_T .$$

Pointwise, exact, 3 random rational $t$, sources at all non-root vertices, 45 graphs (bananas $k=2..5$; all tight graphs of $(H,n)=(2,1..4),(3,1..3),(4,1..2)$ incl. non-SP), $C$ from a random tree, $D_{dy}$ for a second random tree with forward-mode derivatives:

```
A1 (basis independence) A2 (identity above) A3 (comp. scaling s^n) A4 (det(g)^n): True for 45/45 graphs
A5 (D_dy = sum_S prod_r [lam_S_r] dens(wedge_r dR_{c_r,S_r})): True for the 16 graphs with <= 3000 products
all pointwise checks: True
```

Integrated, exact GM quadrature in own coordinates (index $s=1$ and $s=2$ agree), random sources:

```
polygon4 H=3: K(s=1) = 81/2, K(s=2) = 81/2, K(g lam) = det(g)^1 K: True, engine 81/2
paths(1, 3, 3) H=3 n=2: K = 4880/27, stable True; K(g lam) = det(g)^n K: True; K(scaled comp 0 by 4) = s^n K: True; engine 4880/27
paths(2, 3, 4) H=4 n=2: K = 134801/120, stable True; K(g lam) = det(g)^n K: True; K(scaled comp 2 by 4) = s^n K: True; engine 134801/120
(14 cases: polygons m=3..6, bananas k=2..5, thetas (1,2,2), (1,3,3), (2,2,3), (1,4,4), (2,3,4), (3,3,3): all True)
```

Gaps: the hypothesis "the kernel integral converges absolutely for every $\lambda$" is assumed, not shown. Non-series-parallel graphs: pointwise identities only, no integrated value (explorer and here).

## (b1) confirmed

$$2(|I|-d)-(3(|V|-1)-4)=3-2d\le0\ \Rightarrow\ d\ge2,\qquad 3|V|-2|I|=4,$$
$$[124]=-\Delta,\quad[134]=\Delta,\quad[234]=-\Delta,\qquad \Delta(p_{s1},p_{s2},p_{s3})=\operatorname{sgn}(s)\Delta\pmod{\textstyle\sum p},$$
$$\mathrm{alt}(\Delta^nK)=\Delta^n\sum_s\operatorname{sgn}(s)^{n+1}K^s\pmod{\textstyle\sum p}\quad(n=1,2,3,\ \text{generic }K\text{ on }2^4\text{ index values}).$$

`v_b1.log`: all True. Tight $D=3$ graphs with all degrees in $\{2,3\}$: $n=1$: 1 (box), $n=2$: 2 (both thetas), each with four vertices of degree 2.

Gaps: $O_n=c_n\int\Omega\,\Delta^nK\,A^4$ with $K=\sum_\Gamma\kappa_\Gamma W^\Gamma/|\mathrm{Aut}\,\Gamma|$ is not taken from a source equation for hCS on $\mathbb C^3$ (S2 eq-Qn is the assembly of $Q_n$ with an operator insertion); $c_n$ and the Koszul sign between $W^\Gamma$ and $K_\Gamma$ at $H=3$ are not fixed.

## (b2) confirmed

Steps (1)–(4) re-derived: $f(z,\cdot,\cdot)=0$; $\kappa(z,[x,y])=\kappa([z,x],y)=0$; external legs sit on $f$-slots, internal central lines meet $f$; $W^\Gamma|_{\mathfrak s}$ is $\mathfrak s$-invariant. Step (5): $H^*(\mathfrak s)=\Lambda(x_{2m_i+1})$, $2m_i+1\ge3$ odd, so $H^4(\mathfrak s)=0$; $(\Lambda^4\mathfrak g^*)^{\mathfrak g}=\Lambda^4\mathfrak z^*\oplus\mathfrak z^*\otimes H^3(\mathfrak s)$.

`v_lie_ce.log` part 1 (exact nullity via DomainMatrix over $\mathbb Q$, and mod $p$):

```
sl2                                   Alt(W_133)=0 True, Alt(W_223)=0 True; dim(Lambda^4 g*)^g = (0, 0), predicted 0
gl2                                   True, True; (1, 1), predicted 1
sl3                                   True, True; (0, 0), predicted 0
gl3 (twisted basis, tr XY + 2trXtrY)  True, True; (1, 1), predicted 1
so4 = sl2+sl2 (levels 1, 3)           True, True; (0, 0), predicted 0
so5                                   True, True; (0, 0), predicted 0
gl2+gl1 (non-diagonal invariant form) True, True; (2, 2), predicted 2
```

Gaps: step (5) for general semisimple $\mathfrak s$ is the cited theorem (checked on the seven algebras above). The vanishing assumes only the classical cubic vertex; by (b4) the one-loop obstruction is not exact, so the statement concerns the two-loop graphs of the classical vertex only.

## (b3) confirmed (sign remark on the S2 comparison)

$$K_{\theta(1,3,3)}=-\tfrac1{432}(2[123]+[124])(2[134]+[234]),\qquad K_{\theta(2,2,3)}=\tfrac1{288}(2[123]+[124])([134]+[234]),\qquad K_{\rm box}=-\Delta/6 .$$

`v_kernels.log` (exact GM, own solver; $s=1$ and $s=2$ agree):

```
theta(1, 3, 3): K equals the claimed polynomial at 4 random points: True; mod conservation: True
theta(2, 2, 3): K equals the claimed polynomial at 4 random points: True; mod conservation: True
  explorer's numeric point (Delta = 80): K_133 = -400/27, K_223 = 0, K_box = -40/3
box (polygon(4), edges (k,k+1 mod 4), root 3): K = -Delta/6: True
graph A: -K = S2 eq-IA(lambda;0) at 3 random source sets: True
graph B: -K = S2 eq-IB(lambda;0) at 3 random source sets: True
```

`v_b1.log`: $|\mathrm{Aut}|$ = 8 (box), 4 ($\theta(1,3,3)$), 4 ($\theta(2,2,3)$), 2 (S2 graph A with $\mathcal O$ marked). Cochains: $\mathrm{alt}(K_\theta W_\theta)=\kappa_\theta\Delta^2\,\mathrm{Alt}(W_\theta)=0$ by (b2) part 1; box: nonzero, closed (b4).

Sign remark. On $\sum\lambda=0$, S2 eq-IA in S2's own labelling gives (`v_S2_conservation.log`)

```
A: Delta' = -248, I_A(eq-IA) = -3844/27, -Delta'^2/432 = -3844/27
A: Delta' = 50,   I_A(eq-IA) = -625/108, -Delta'^2/432 = -625/108
```

i.e. $I_A|_{\sum\lambda=0}=-\Delta'^2/432$, $\Delta'=\det(\lambda_1,\lambda_2,\lambda_4)$, whereas the report states $I_{133}=+\Delta^2/432$. The two labellings differ by the vertex transposition $(0\,2)$ and four edge reversals, total sign $(-1)^3(-1)^{3\cdot4}=-1$, so the two values agree after this sign; "$I_{133}=\Delta^2/432$, consistent with eq-IA" holds only for $I_{133}:=-K$ in engine labelling.

## (b4) confirmed

$\omega_{\rm box}=\mathrm{alt}(\det(p_1,p_2,p_3)\,W_{\rm box})$, weight $(1,1,1)$; $\eta$: 6 trace types $\times$ 8 monomials.

```
sl2: closed True; rank D 0, rank [D|w] 1      gl2: closed True; rank D 2, rank [D|w] 3
sl3: closed True; rank D 1, rank [D|w] 2      gl3: closed True; rank D 3, rank [D|w] 4
sl4: closed True; rank D 1, rank [D|w] 2
large-range closedness (R = 30, 6 samples): sl2, gl2, sl3, gl3 True
```

## (c1) confirmed

Pointwise, exact, own $\rho_2$, sources at all five non-root vertices ($\lambda_N$ included), 5 random points each:

$$U=\frac1{w_Aw_Bw_C},\qquad q\,\lvert\det J\rvert=s\,\frac{D_\theta}2\prod_Pw_P^{H-n_P},\qquad g=qU^{-1/2}\lvert\det J\rvert=s\,\frac{D_\theta}2\prod_Pw_P^{H+\frac12-n_P},$$
$$K=s\frac{D_\theta}2\frac{\prod_P\Gamma(H+\frac32-n_P)}{\Gamma(H+\frac32)\prod_P(n_P-1)!}:\quad \theta(1,3,3):\ s=-1,\ -\tfrac{\pi}{20}D_\theta;\qquad \theta(2,2,3):\ s=+1,\ \tfrac{\pi}{30}D_\theta .$$

```
theta(1, 3, 3) T=1: ... with s = -1 at 5 random points/sources: True ; K/D_theta = -0.15707963267949 (= -pi/20)
theta(2, 2, 3) T=1: ... with s = 1 at 5 random points/sources: True ;  K/D_theta =  0.10471975511966 (= pi/30)
numerical integral of the density from the definition (float, Gauss-Legendre after smoothing maps):
theta(1,3,3), D_theta = -8:   exact 1.256637061435917, m=12: 1.256637061435911, rel.dev 4.6e-15 (14 digits)
theta(2,2,3), D_theta = -115: exact -12.042771838760872, m=12: -12.042771838760734, rel.dev 1.2e-14 (13 digits)
```

## (c2) unverified

Status derived-unchecked; not re-derived. Checked only: $(1/(8\pi^{5/2}))^2=1/(64\pi^5)$; $\pi^{-1}\cdot\pi/20=1/20$; $\int 4\cdot\tfrac14 e^{-r^2/4t}(4\pi t)^{-5/2}d^5x=1$.

## (c3) unverified

Status derived-unchecked. Partial checks: the index forms of $W_{133}$, $W_{223}$ equal the orientation rule applied to the engine graphs, sign $\operatorname{sgn}(L\to L')=-1$ for both (sl2, gl2, sl3, gl3). `v_orientation.log`: $\mathrm{alt}(K_\Gamma W^\Gamma)$ ($(H,T)=(2,1)$, own kernel and own rule) unchanged under 4 random vertex relabellings (root moved), edge reversals and edge-order permutations: True for both thetas, sl2 and gl2. The Wick-sign identity (relative sign between $\theta(1,3,3)$ and $\theta(2,2,3)$) was not tested.

## (c4) confirmed

Evaluation at 6 random integer matrix 4-tuples, $W$ from the index forms (exact tensors, own basis), $\kappa=\mathrm{tr}$:

```
sl2: W_tot gl_N formula True; W_133 sl_N formula True; W_223 sl_N formula True; sl2 -104 k12 k34 + 16 k13 k24 + 40 k14 k23: True
gl2, gl3, gl4: W_tot = -3W_133 + 2W_223 equals the gl_N trace formula: True
sl3, sl4: W_tot True; W_133 True; W_223 True
```

$O_2=\frac{N_{\rm norm}}4\int\Omega\,(K_{133}W_{133}+K_{223}W_{223})A^4=\frac{\pi N_{\rm norm}}{240}\int[d_1d_2][d_3d_4]W_{\rm tot}$, $-\frac{\pi}{20}=-\frac{3\pi}{60}$, $\frac{\pi}{30}=\frac{2\pi}{60}$.

## (c5) confirmed

```
d omega_133 = 0, d omega_223 = 0: sl2, gl2, sl3, gl3 (4 samples), sl4, gl4, gl5, sl5 (2 samples): True
large-range (R = 30, 6 samples, Schwartz-Zippel bound 1.0e-05): sl2, gl2, sl3, gl3: True
omega_133, omega_223 nonzero: True for all eight algebras
```

## (c6) confirmed

```
             rank D   rank [D|omega]   exact
sl2            2          2            True (eta re-verified, 6 large samples)
gl2            5          5            True (eta re-verified, 6 large samples)
sl3            3          4            False
gl3            6          7            False
sl4            3          4            False
gl4            6          7            False
gl5            6          7            False
sl5            3          4            False
```

All $N\ge4$ (`v_formalN.log`). Formal $d$ is $N$-independent; $\omega(N)=\sum_kN^k\omega_k$ from (c4) (degree $\le2$: at most 3 faces of the cubic ribbon graph, at least one carries legs; formula checked at $N=2,3,4$); at $N_0=5$ the 24 (gl) / 9 (sl) four-argument and 6 / 2 three-argument trace monomials are independent functions, so evaluation is injective and $\omega(N)$ is exact iff $(1,N,N^2)\in R$:

```
sl5: rank D = 3; omega: rank [D|omega_0,omega_1,omega_2] = 4; R = span{(1,0,0), (0,0,1)}  ->  [omega(N)] = N[omega_1] != 0
gl5: rank D = 6; omega: rank [D|omega_0,omega_1,omega_2] = 8; R = span{(0,0,1)}           ->  (1,N,N^2) not in R for any N
sl6: same as sl5
```

gl$_N$ versus sl$_N$: $W^\Gamma$ vanishes on the centre, so the gl$_N$ cochain is the pull-back along the homomorphism $X\mapsto X-\tfrac{\mathrm{tr}X}N$; restriction and pull-back commute with $d$, hence exact on gl$_N$ iff exact on sl$_N$.

## (c7) confirmed

$$\omega_{133}=-3F_{133},\quad\omega_{223}=2F_{223},\quad\omega_{133}-\tfrac32\omega_{223}=-3(F_{133}+F_{223}),\quad O_2=\tfrac{N_{\rm norm}}4\big(\tfrac{\pi}{20}+\tfrac{\pi}{30}\big)F_{223}=\tfrac{\pi N_{\rm norm}}{48}F_{223}\ \text{mod exact}.$$

```
omega_133-(3/2)omega_223 exact: sl2, gl2, sl3, gl3, sl4, gl4, gl5, sl5 True (re-verified at large samples: sl3, gl3, gl4)
F_133+F_223 exact: all eight True
F_223 exact: sl2 True, gl2 True; sl3, gl3, sl4, gl4, gl5, sl5 False
formal sl_N (N0 = 5, 6): F_133+F_223 = sum_k N^k (...)_k with R = span{(1,0,0),(0,1,0),(0,0,1)}: exact for every N
formal sl_N: F_223 = F_0 + N F_1, R = span{(1,0)}: [F_223(N)] = N [F_1] != 0 for N >= 4
```

Remark: with one kernel sign flipped the class is $\pm\frac{\pi N_{\rm norm}}{240}[F_{223}]$; it is $-$ when the sign of $K_{133}$ is flipped and $+$ when the sign of $K_{223}$ is flipped.

## (c8) confirmed (sample ranks)

gl$_5$ and gl$_6$, explicit matrices, 864 spanning cochains $\mathrm{alt}(\mu\,\mathrm{tr}_\sigma)$, $\sigma\in S_4$, $\mu$ of weight $(2,2)$ in $p_1,p_2,p_3$:

```
gl5: dim C^4_inv(2,2) = 36 (70 samples), 36 (90 samples); rank d: C^4 -> C^5 = 26 (60), 26 (75); dim Z^4 = 10
     dim B^4 = rank d: C^3 -> C^4 = 6 ; dim H^4 = 4
gl6: dim C^4 = 36, 36; rank d = 26, 26; dim Z^4 = 10
gl5: the 10 kernel vectors of d are cocycles at 4 fresh samples (range 12) and independent as cochains (rank 10)
```

Ranks of sample value matrices are exact over $\mathbb Q$ and are rigorous lower bounds; the upper bounds rest on stability under additional random samples.

## Logical gaps (REPORT.md was not written by the explorer; gaps taken from the claims and scripts)

1. (a1) Absolute convergence for every $\lambda$ is a hypothesis; no convergence proof for the tight graphs is given.
2. (b1) The obstruction formula, its constant $c_n$ and the $H=3$ Koszul sign between $W^\Gamma$ and $K_\Gamma$ are not derived from a source equation for hCS on $\mathbb C^3$.
3. (b2) Step (5) for arbitrary semisimple $\mathfrak s$ is the cited theorem. The vanishing assumes no counterterm or Green–Schwarz vertices. Since the one-loop obstruction is not exact (b4), the order-2 obstruction is defined only after something cancels the one-loop class, and that is not treated.
4. (b3) The S2 comparison holds only up to the labelling sign $-1$ (numbers above).
5. (b4, c5–c8) Two steps are argued in the report, not taken from the sources: the identification of local BRST cohomology with the CE cohomology of $\mathfrak g[[z]]$ modulo $\sum p$ (freeness over $\mathbb C[\partial]$), and the claim that the $\tau$ direction gives only a shift.
6. (c4, c6, c7) The relative sign of $K_{133}W_{133}$ and $K_{223}W_{223}$ comes from (c3), which is derivation only. Non-exactness of $O_2$ does not depend on it (both $\omega_{133}$, $\omega_{223}$ separately non-exact for $N\ge3$), but the class coefficient $\frac{\pi N_{\rm norm}}{48}$ does.
7. (c2) Overall constant and sign of $O_2$: derivation only.
8. (c8) Upper bounds on dimensions are probabilistic (sample stability), both here and in the explorer's mod-$p$ computation.

## Not done

Integrated values for non-series-parallel graphs. A float direct numerical integral for the $H=3$ thetas in the original Schwinger coordinates (exact GM quadrature in own coordinates was used instead). Re-derivation of (c2) and of the Wick-sign part of (c3).
