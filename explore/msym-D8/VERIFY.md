# msym-D8: independent verification

Scripts and logs: `P/explore/msym-D8/verify/`. No script imports `htk` or any file of `P/explore/msym-D8/code/`; the explorer's scripts were read for conventions only and not run. `P/explore/msym-D8/REPORT.md` does not exist (the folder contains only `code/`); logical gaps are read from the structured report and the code.

## 0. Independent layer (`verify/common.py`)

Graphs built from the claim texts; edge lists coincide with the explorer's logs:
pentagon `[(4,0),(0,1),(1,2),(2,3),(3,4)]`, O = 4; theta(p,q,r): O = 7, N = 0, paths P1,P2,P3 oriented away from O, edges in the order P1,P2,P3.

Mesh-current solver (BFS spanning tree, not Kruskal; no node Laplacian):

$$y=y^0+C^{\mathsf T}\ell,\quad (CtC^{\mathsf T})\ell=-Cty^0,\quad \frac{\partial\ell_r}{\partial t_j}=-(M^{-1}C_{:,j})_r\,y_j,\quad D(t)=\frac{(-1)^{e_0}}{t_{e_0}}\det\Big[\frac{\partial\ell^a_r}{\partial t_j}\Big]_{(r,a);\,j\ne e_0},\ e_0=N-1,$$

$$K_\Gamma(\lambda;z)=\int_{t>0,\ \sum t=1}D(t)\,e^{-\sum_ey_e\cdot z_e}\,dt_0\cdots dt_{N-2}.$$

$H=4$: $\ell'=A\ell+\text{const}$ with $\det A=\pm1$ gives $(\det A)^4=1$, $\sigma_T^4=1$, $\epsilon(n,4)=(-1)^{4n+6\binom n2}=1$, so $K_\Gamma$ does not depend on the tree. The orientation is $\iota_E\,dt_0\wedge\cdots\wedge dt_8$ in the edge order above.

## C2. Edge and vertex factors (`v1_edge_vertex.py`)

$$[K_t]_\Omega=\frac{e^{-x\cdot\bar x/2t}}{2^d(2\pi t)^d}d\bar x^{1..d},\qquad \bar\partial^*=-2\sum_a\iota_{\partial/\partial\bar x^a}\partial_a,\qquad y=\bar x/2t .$$

```
(a) d = 4: dx^{1..d} ^ dxbar^{1..d} = 16 d^{8}x_E          (d = 3: -8*I, d = 2: 4)
(b) d = 4: int Omega ^ [K_t] = 1                          (d = 3: 1)
(c) d = 3: [K_t] - dt ^ dbar^*[K_t] == (2 pi)^{-3} e^{-y.x} d^3y = 1/(8 pi^3) e^{-y.x} d^3y: True   [S2 eq-y]
(c) d = 4: [K_t] - dt ^ dbar^*[K_t] == (2 pi)^{-4} e^{-y.x} d^4y = 1/(16 pi^4) e^{-y.x} d^4y: True
    d = 4: dbar^*[K_t] has form degree [3] -> propagator parity odd   (d = 3: degree [2], even, as S2 Sec. 5.1)
(d) c_Omega i^d (2 pi)^d = 8 pi^3 (d = 3, c_Omega = i, S2 eq-vertex), 16 pi^4 (d = 4, c_Omega = 1)
```

Verdict: confirmed.

## C5. Pentagon (`v2_pentagon.py`)

Derivation: $y_e=y_{e_0}+Q_e$, $Q_e=\sum_{v<e}\lambda_v$; on $\sum t=1$, $\partial y_{e_4}/\partial t_j=-y_j$;

$$D=\frac{\det(y_0,y_1,y_2,y_3)}{t_4}=\frac{\det(y_0,\lambda_0,\lambda_1,\lambda_2)}{t_4}=\frac{-t_4\det(\lambda_3,\lambda_0,\lambda_1,\lambda_2)}{t_4}=\langle0123\rangle,$$

$$y_{e_0}=-\sum_jr_j\lambda_j,\quad y_{e_4}=\sum_j(1-r_j)\lambda_j,\quad -y_{e_0}\cdot z_1+y_{e_4}\cdot z_2=\sum_j\lambda_j\cdot x_j,\quad \int_{1>r_1>\dots>r_4>0}e^{\sum_ja_jr_j}=\exp[0,s_1,s_2,s_3,s_4],\ s_k=\textstyle\sum_{j\le k}a_j .$$

```
(a) point 0..5: D(t) T^5 = <0123> (1812, -857, 96, -3606, 648, 310): True; y_e0 = -sum r_j lam_j: True; -sum y.z = sum lam.x: True
(d) moment formula == iterated sympy integral for 210 exponent tuples (a+b+c+d <= 6): True
(b) order 0..3 in z, t-simplex Dirichlet integral == <0123> int_r (sum lam.x)^k/k!: -35/4, 931/4, -363839/108, 29693887/864: True
(c) trial 0: quadrature (14^4 Gauss nodes, fresh float solver) = -0.154059498367233, <0123> e^b exp[0,s1..s4] = -0.1540594983672329, 15 digits
(c) trial 1: -2.37148942701332 vs -2.37148942701332, 15 digits;  trial 2: -0.616358333290734 vs -0.6163583332907343, 15 digits
v6: pentagon at the explorer's QMC sources: <0123>/24 = 9/2
```

(c) is the Schwinger integral over the t-simplex at $z\ne0$, $|z^a|\le0.4$. Verdict: confirmed.

## C6. Theta displacement maps (`v3_thetas.py` (i), (ii))

Derivation (KCL at interior vertices, $\lambda_N=0$, equal potential drop $\Phi$ on the three paths):

$$y_{j,k}=\ell_j+\sum_{i<k}\lambda_{v_{j,i}},\quad \sum_j(\ell_j+\Lambda_j)=0,\quad \ell_jR_j+R_ja_j=\Phi,\quad a_j=\sum_{v\in P_j}\sigma_v\lambda_v,$$

$$\ell_j=-c_j\sum_v(1-\sigma_v)\lambda_v-\sum_{v\in P_j}\sigma_v\lambda_v,\qquad -\sum_j\ell_j\cdot z_j=\sum_v\lambda_v\cdot\big(\sigma_vz_{j(v)}+(1-\sigma_v)Z_N\big).$$

```
theta(1,4,4), (2,3,4), (3,3,3): (i) -sum_e y_e.z_e == sum_v lam_v.(sigma_v z_j(v) + (1-sigma_v) Z_N), 5 rational points: True
(ii) spoke currents of the mesh solver = own formula for l_j: True at 4 points per graph
```

Verdict: confirmed.

## C7, C8, C9. Theta kernels (`v3_thetas.py`, `v6_values.py`)

Own coordinates for all three thetas (three-chain move; the explorer used two moves for (1,4,4), (2,3,4)): $\xi=(c_2,c_3,x_{j,k\ge2})$, section $t_{j,k}=x_{j,k}/c_j$;

$$g(\xi)=D(t(\xi))\,|\det[t;\partial_\xi t]|=\mathrm{sgn}\cdot\det\Big[\frac{\partial(\ell_2,\ell_3)}{\partial\xi}\Big],\qquad K_\Gamma(\lambda;0)=\int_{\Delta^2\times\prod_j\Delta^{L_j-1}}g .$$

```
(ii) D(t) det[t;dt/dxi] == det[d(l2,l3)/dxi], 4 rational points per graph: True; sgn = -1 (144), +1 (234), +1 (333)
(iii) C7 claimed density -om1^3 u^4 det(l1,l2,l3,mu3) det(mu2,l4,l5,l6) == D|det[t;dt]| in the explorer's two-move section, 4 points: True
(iii) C8 claimed density -om1^2 om2 u^4 det(l1,l2,l3,mu3) det(mu12,l4,l5,l6) == D|det[t;dt]|, 4 points: True
(iii) C9 g = det L == D|det[t;dt]| (sign +1), 4 points: True
(iv) exact integral of sgn*det[dl/dxi] (sympy DomainMatrix det, numeric lam) == claimed closed form:
     (1,4,4): -543433/480, -433927/192, 16263/64; (2,3,4): -8979/40, 279911/720, -312137/720; (3,3,3): -272207/120, -5920901/2160, -27533/135: True
v6: claimed closed forms expanded in monomials == own K(lam;0) polynomial: (1,4,4) 4968 = 4968 True; (2,3,4) 4968 = 4968 True; (3,3,3) 6624 = 6624 True
v6: at the explorer's QMC sources: -287899/1440, -16549/80, -2243657/8640 (own = closed form = explorer exact): True
(vi) tensor Gauss (5 nodes per direction, 390625 points) of the mesh-solver float density D|det[t;dt]| (no reduction formula):
     (1,4,4): K = -2260.0364583333 vs -433927/192, rel 1.6e-15; first order (solver spoke currents) 686.04307870371 vs 29637061/43200, rel 3.5e-15
     (2,3,4): K = 388.76527777778 vs 279911/720, rel 1.0e-15; first order 4127.4428935185 vs 178305533/43200, rel 1.8e-15
     (3,3,3): K = -2741.1578703704 vs -5920901/2160, rel 1.5e-15; first order 96237.902263374 vs 93543241/972, rel 1.2e-15
(vii) (1,4,4), P2 <-> P3: K(lam o pi; z o pi)/K(lam; z) = 1 at z^0 and z^1, 3 points
(vii) (3,3,3), P1<->P2, P2<->P3, P1<->P3: ratio -1 at z^0 and z^1, 3 points each
```

The integrand in these coordinates is a polynomial of degree at most 9 in each Duffy variable, so the 5-point rule is exact up to rounding (14-15 agreeing digits). Verdicts: C7, C8, C9 confirmed.

## C11. $\{\hat Q_0,\hat Q_1\}=0$ (`v4_ops.py`, `v4_c11.py`)

Own component form ($C=c+\theta\gamma$, $\delta_j=\theta\partial_{c_j}-\partial_{\gamma_j}$, $F=F_0+\theta F_1$):

$$\delta_1\delta_2=-\theta(\partial_{c_1}\partial_{\gamma_2}+\partial_{\gamma_1}\partial_{c_2})+\partial_{\gamma_1}\partial_{\gamma_2},\qquad \int d\theta\,F\,\delta_1\delta_2=-F_0(\partial_{c_1}\partial_{\gamma_2}+\partial_{\gamma_1}\partial_{c_2})+F_1\partial_{\gamma_1}\partial_{\gamma_2},$$

$$\hat Q_0c^A=\tfrac12f_{ABC}c^Bc^C,\qquad \hat Q_0\gamma^A=f_{ABC}c^B\gamma^C,$$

$W_1(k_1,k_2;m)=[u^m]\int_{1>r_1>\dots>r_4>0}\prod_a(\sum_qr_qu_{qa})^{k_{1a}}(\sum_q(1-r_q)u_{qa})^{k_{2a}}$. Test operators and the algebra so(4) are not the explorer's.

```
su2: hatQ0^2 = 0: True; C10 closed forms of hatQ1(g g), hatQ1(g c), hatQ1(c c) = 0 for (A,B) = (0,1), (1,1), (2,0): True
{hatQ0,hatQ1} = ZERO on: g0, d3g1, d2d4g2, g1 g2, g2 g2, c1 g2, d3g0 d3g1, d4g2 c0, g1 d2d3c2, d1d4g0 g2, d2g1 d4c1,
                         d3d3g2 g0, g1 g2 d4c0, g0 g0 g2, d2g0 g1 c1   (15 operators, up to two derivatives)
so4: hatQ0^2 = 0: True; closed forms for (0,5), (2,2), (1,3): True
{hatQ0,hatQ1} = ZERO on: g0, d2g3, g0 g5, g1 g1, g2 c3, d3g0 g4, d2g1 d4c4, g0 g3 g5
CONTROL (moments of 0 < r_1 < .. < r_4 < 1): NONZERO on d3g0 d3g1 (380 terms), d1d4g0 g2 (870), d2g1 d4c1 (142)
```

Verdict: confirmed.

## C12. $\hat Q_1^2+\sum_\Gamma k_\Gamma\{\hat Q_0,\hat Q_\Gamma\}=0$, su(2) (`v4_c12.py`)

$$\int d\theta\,F\,\delta_1\delta_2\delta_3=F_0(\partial_{c_1}\partial_{\gamma_2}\partial_{\gamma_3}+\partial_{\gamma_1}\partial_{c_2}\partial_{\gamma_3}+\partial_{\gamma_1}\partial_{\gamma_2}\partial_{c_3})-F_1\,\partial_{\gamma_1}\partial_{\gamma_2}\partial_{\gamma_3},$$

$T_\Gamma$ and the output order as stated in C10; theta weights = own $K_\Gamma(\lambda;0)$ and own first-order coefficients $\sum_v\lambda^b_v\int g\,c_{v,j}$ from C7-C9.

```
operator d2c0 g1 g2: 2720 equations, rank 3, consistent True, solution (2, -4, -2/3)
operator d3g2 g2 c0: 2762 equations, rank 3, consistent True, solution (2, -4, -2/3)
operator g1 d4g0 c2: 2720 equations, rank 3, consistent True, solution (2, -4, -2/3)
operator g0 g0 c1: 320 equations, rank 2, consistent True;  g2 g1 c1: 320, rank 2, True;  g0 g1 g1: 2240, rank 2, True
explorer's operators recomputed: g0 g1 c2 (288 eq.), g0 g1 g2 (2016 eq.): rank 1, a = (1, -23/3, 18), b = -62/3, a.(2,-4,-2/3) + b = 0
ALL 8 operators: 13386 equations, rank 3, consistent True, (k144, k234, k333) = (2, -4, -2/3), |k| |Aut| / 4 = (1, 1, 1), residual zero
zero-derivative operators together (5 operators, 5184 equations): rank 2, consistent True
CONTROL (z_1 <-> z_2 in the first-order weights of theta(2,3,4)), d2c0 g1 g2: inconsistent
```

Verdict: confirmed. Correction to a side statement: "on the operators without derivatives alone the rank is 1" holds for the explorer's two operators (distinct colours); zero-derivative operators with a repeated colour give rank 2.

## C4. Graphs (`v5_count.py`)

```
n = 1: 222 connected graphs (5 vertices, 5 edges); sparse (3|I_S| <= 4|V_S| - 5 on all edge subsets): only the 5-cycle
n = 2 cores with 9 edges, 8 vertices: theta (1,4,4), (2,3,4), (3,3,3) sparse; theta (1,2,6), (1,3,5), (2,2,5), all dumbbells, all figure-eights violate
admissible O (Gamma\O forest, deg O = n + p): the two trivalent vertices; |Aut fixing O| = 2, 1, 6; outputs 3(3n+1) - 2(4n+1-deg O) - deg O = 6
pentagon: deg O = 2, p = 1, |Aut fixing O| = 2, outputs 4
```

Trees attached to the 2-core $C$ are excluded by $3|C|\le4(|C|-1)-5\Rightarrow|C|\ge9$ (n = 2), $|C|\ge5$ (n = 1). Theta rule with one $\theta$ per edge: $4n+1-\deg O=3n+1-p\Leftrightarrow\deg O=n+p$. Each $\ell_r$ contributes one factor per direction: $n$ derivatives per direction. Verdict: confirmed.

## C1, C3, C10, C13

C1. Checked by hand: $\int d\theta\frac1{3!}f_{ABC}C^AC^BC^C=\frac12f_{ABC}\gamma^Ac^Bc^C$ for $C=c+\theta\gamma$; $(\theta\partial_c-\partial_\gamma)(c+\theta'\gamma)=\theta-\theta'$; $\frac12(-\frac i4)^2\cdot16=-\frac12\Rightarrow\hbar=g^2/2$; $\Omega\wedge d\bar x^{1234}=16\,d^8x_E$ for $\Omega=dx^{1234}$ (v1 (a)). Not checked: the 8d reduction (contractible pairs, the 8d analogue of S2 eq-Lhol-form with the $\chi_{a5}$ components). Verdict: unverified.

C3. Magnitude: $(16\pi^4)^{-m}(16\pi^4)^{k-1}=(16\pi^4)^{-n}$ from C2. Sign chain (i)-(v) not re-derived. Verdict: unverified.

C10. Confirmed parts: closed forms of $\hat Q_1$ on two letters (C11 run); $|k_\Gamma|=4/|\mathrm{Aut}\,\Gamma|$ reproduced by the independent fit (C12). Not tested: $Q_n\to\beta^nQ_n$ and $Q_1\to-Q_1$ leave the C12 system invariant, so the fit fixes only $k_\Gamma\propto|\mathrm{Aut}\,\mathrm{pent}|^2/|\mathrm{Aut}\,\Gamma|$, not $(16\pi^4)^{-n}$, $\hbar$, $s_4(n)$ or $\varepsilon_W$. Verdict: unverified.

C13. Every term of $\hat Q_1$ changes $\#\gamma$ by $-1$, of $\hat Q_\Gamma$ by $-2$ (component form above), so $\gamma\to-\gamma$ ($\theta\to-\theta$) gives $Q_n\to(-1)^nQ_n$. $\varepsilon_W$ itself is not derived. Verdict: unverified (open).

## Logical gaps

1. C1: the reduction 10d → 8d is quoted from S2 Sec. 3, written there for 10d → 6d with two $\theta$'s, without an 8d computation of $\chi\wedge\bar\partial\chi$. $c_\Omega=1$ is a convention (S2 gives no rule on $\mathbb C^4$), and $\hbar$ depends on it.
2. C3: the sign chain reproduces S1 eq-1-7 at d = 2 but differs from S2 eq-B-norm by $(-1)^n$ at d = 3 (explorer's open issue, unresolved). $s_4(n)=+1$ rests on the same chain. For D = 4, $\varsigma^{\rm ext}=+1$ follows from $\binom42$ being even, so it does not test the chain.
3. C10/C12: the signs of $k_\Gamma$ (= $\varepsilon_W$ of each theta relative to the pentagon) are fitted, not derived. The statement "$Q_1^2+\{Q_0,Q_2\}=0$ holds" means: a sign assignment with $|k_\Gamma|=4/|\mathrm{Aut}\,\Gamma|$ exists. The fit does not test $(16\pi^4)^{-n}$ or $\hbar$ (scaling invariance above).
4. C12: only su(2), and only operators of total derivative order $\le1$ (first-order theta weights). The identity is third order, so 3-letter operators suffice in number of letters, but not in derivative order.
5. C12 side statement on the rank of the zero-derivative operators: see the correction under C12.
6. C9: "three-chain move (no two path lengths sum to 5)" is a criterion for the engine's local formula. The kernel value does not depend on it: the generic Jacobian in the same coordinates is used above for all three thetas.
7. C7-C9 QMC evidence (4-5 digits) is weak. Replaced above by an exact-polynomial Gauss rule on the mesh-solver density (14-15 digits).

## Not done

so(5) or other algebras for C12; second-order theta weights; the C3 sign chain; the 8d reduction of C1; $\varepsilon_W$.
