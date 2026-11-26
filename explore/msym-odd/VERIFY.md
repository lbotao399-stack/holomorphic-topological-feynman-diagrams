# VERIFY: msym-odd

Scripts: `verify/` (no import of `htk` or of the explorer's `code/`). Logs: `verify/logs/`.
REPORT.md does not exist in `explore/msym-odd/`; the gaps below refer to the structured report and to the explorer's code.

| script | content |
|---|---|
| `vkern.py` | exact Kirchhoff solver ($\Lambda\varphi=\lambda$, $y=t^{-1}B\varphi$, both Kirchhoff laws asserted), BFS-tree cycle basis with $\det[B^{\mathsf T};C]>0$, $\rho_2=U^{-1/2}\bigl(-\tfrac12\sum_{ab}dM_{1a}M^{-1}_{ab}\wedge dM_{b2}\bigr)$, $\psi^i_r$, wedge on $dt$-space, $D(t)=[\Omega]_{0..N-2}/((-1)^{N-1}t_{N-1})$ with the extraction repeated for every omitted edge; chart $t^P_j=x^P_j/w_P$, exact $\lvert\det[t;\partial_\xi t]\rvert$; Dirichlet integrals |
| `kernels_verify.py` | claims theta-kernels, Q1-zero |
| `numeric_verify.py`, `numeric_wx.py` | direct numerical integrals of the brute-force density |
| `twist_verify.py` | claims twist-5d, action-hbar, propagator-T1 |
| `graphs_verify.py` | claim graphs-Q2 |
| `valg.py`, `vops.py`, `nilp_verify.py` | Grassmann algebra, $\hat Q_0$, $\delta/\delta C(\theta)$, $T_\Gamma$, $\hat Q_\Gamma$, exact rank/null space; claims nilpotency-5d, 7d-relative-coefficients, 3d-closure |

## twist-5d — confirmed

$x^{H+1}=\tau+i\sigma$, $\partial_\sigma=0$: $\partial_{\bar x^{H+1}}=\partial_{x^{H+1}}=\tfrac12\partial_\tau$.

$$\Phi(d\bar x^{H+1}\,\tfrac12\partial_\tau\omega)=d\tau\,\partial_\tau\Phi(\omega)\ \Rightarrow\ \Phi\bar\partial_{C^{H+1}}=(\bar\partial+d_\tau)\Phi,$$
$$-2\iota_{\overline{H+1}}\,\tfrac12\partial_\tau(d\bar x^{H+1}\alpha)=-\partial_\tau\alpha,\qquad -\tfrac12\iota_\tau\partial_\tau(2d\tau\,\Phi\alpha)=-\partial_\tau\Phi\alpha\ \Rightarrow\ \Phi\bar\partial^*=D^*\Phi,$$
$$[d\tau\,\partial_\tau,-\tfrac12\iota_\tau\partial_\tau]=-\tfrac12\partial_\tau^2,\quad[\bar\partial,-2\iota_{\bar a}\partial_a]=-2\partial_a\partial_{\bar a}=-\tfrac12(\partial_{2a-1}^2+\partial_{2a}^2)\ \Rightarrow\ [D,D^*]=-\tfrac12\Delta_{\mathbb R^{2H+1}} .$$

$C\supset-i\,d\bar x^{3}A_{\bar 3}$, $A_{\bar3}=\tfrac12(A_5+iA_6)$ (S2 eq-cx, eq-C-dict) $\Rightarrow\Phi(C)\supset-i\,d\tau(A_5+iA_6)$. $\Phi$ rescales one odd generator, so it is an algebra map and $\Phi(\bar\partial C+\tfrac i2fCC)=D\Phi C+\tfrac i2f\Phi C\Phi C$.

`logs/twist_verify.log` (1):
```
H=1: Phi dbar = (dbar+d_tau) Phi: True; Phi dbar* = D* Phi: True; [D,D*] = -Lap/2: True; pull-back: Phi' dbar = (dbar + d_tau/2) Phi': True
H=2: ... True; True; True; True
H=3: ... True; True; True; True
```
Not checked: "the edge form is then not basic under the scaling" (pull-back convention).

## action-hbar — confirmed

$dx\wedge d\bar x=-2i\,dx_E^{1}dx_E^{2}$:
$$dx^1dx^2\,d\bar x^1d\bar x^2\,d\tau=-(-2i)^2d^5x_E=4\,d^5x_E,\quad 2i\,dx\,d\bar x\,d\tau=4\,d^3x_E,\quad \tfrac i2dx^{1..3}d\bar x^{1..3}d\tau=\tfrac i2(-8i)\,d^7x_E=4\,d^7x_E .$$
$\chi\wedge\bar\partial\chi=-\tfrac14 e^je^ke^ie^le^m\chi_{jk}\partial_i\chi_{lm}=e^1\cdots e^5\mathcal L_{\rm hol}$; with $C^{(2)}=-\tfrac i4\Phi(\chi)$, $\Phi(e^{H+1})=2d\tau$:
$$\int\Omega\wedge\Bigl[\int d^k\theta\,\tfrac12C^{(2)}DC^{(2)}\Bigr]=\tfrac12\bigl(-\tfrac i4\bigr)^2\cdot2\cdot4\int d^Dx_E\,\mathcal L^{\rm red}_{\rm hol}=-\tfrac14\int d^Dx_E\,\mathcal L^{\rm red}_{\rm hol},\qquad \hbar=\tfrac{g^2}4 .$$
The $\theta$-ordering sign $(-1)^{k(H+T)}$ is $+1$ for $(H,T,k)=(1,1,3),(2,1,2),(3,1,1)$.

`logs/twist_verify.log` (2), (3):
```
H=2 T=1 Omega = 1 dx^1..dx^H : Omega ^ dxbar^1..dxbar^H ^ dtau = 4 d^5x_E      (also H=1: 4, H=3: 4, (3,0): 8)
10d: chi ^ dbar chi = e^1..e^5 L_hol: True
(H,T,k)=(3,0,2): ... = {-1/4} (same support: True)
(H,T,k)=(1,1,3): ... = {-1/4} (same support: True)
(H,T,k)=(2,1,2): ... = {-1/4} (same support: True)
(H,T,k)=(3,1,1): ... = {-1/4} (same support: True)
```
Input not re-derived: S2 eq-decomp and eq-Stw-phys reduced to $\partial_6=\dots=\partial_{10}=0$ (taken from S2 §3, §8.1 item 1).

## propagator-T1 — confirmed

$$\int\Omega\wedge K_t=4\,\frac{(2\pi t)^{H+1/2}}{4(2\pi t)^{H+1/2}}=1,\qquad D^*K_t=\Bigl(\tfrac{\bar x^a}t\iota_{\bar a}+\tfrac{\tau}{2t}\iota_\tau\Bigr)K_t=t^{-1}\iota_EK_t,$$
$$d^Hy\wedge d\eta=(2t)^{-H-\frac12}\bigl[d^H\bar x\,d\tau-dt\wedge t^{-1}\iota_E(d^H\bar x\,d\tau)\bigr]\ \Rightarrow\ \frac{e^{-\eta^2-y\cdot x}}{4\pi^{H+1/2}}d^Hy\wedge d\eta=K_t-dt\wedge D^*K_t,$$
$$(d_t+D)\Theta=dt\wedge(\partial_t+[D,D^*])K_t=0,\qquad [P]_\Omega=\iota_E(d^H\bar x\,d\tau)\int_0^\infty\frac{dt\,e^{-r^2/2t}}{4(2\pi)^{H+\frac12}t^{H+\frac32}}=\frac{\Gamma(H+\frac12)}{4\pi^{H+\frac12}r^{2H+1}}\iota_E(d^H\bar x\,d\tau),$$
$D[P]=\int dt\,DD^*K_t=-\int dt\,\partial_tK_t=0$ ($r>0$); $\deg P=H$.

`logs/twist_verify.log` (4):
```
H=1: int Omega^K = 1; (d_t + [D,D*])K = 0: True; D*K = iota_E K/t: True; K - dt^D*K = ... d^Hy^deta: True; (d_t + D)Theta = 0: True; P coefficient 1/(8*pi) (...: True); D P = 0: True; P degree 1
H=2: ... P coefficient 3/(16*pi**2) (...: True); D P = 0: True; P degree 2
H=3: ... P coefficient 15/(32*pi**3) (...: True); D P = 0: True; P degree 3
```

## raw-normalization — unverified

Re-derived ingredients: with $e^{-\frac12\tau^{\mathsf T}\Lambda\tau}$, $\langle\tau\tau^{\mathsf T}\rangle=\Lambda^{-1}$, $w_e=\tau_e/2t_e$; $B\Lambda^{-1}B^{\mathsf T}t^{-1}+tC^{\mathsf T}M^{-1}C=1$ (complementary projectors, $CB=0$) gives
$$\langle w_ew_f\rangle=\tfrac14\bigl[\delta_{ef}/t_e-(C^{\mathsf T}M^{-1}C)_{ef}\bigr],\qquad\langle\psi^\tau_r\psi^\tau_{r'}\rangle=-\tfrac14(dM\,M^{-1}dM)_{rr'},\qquad\Bigl\langle\prod_r\psi^\tau_r\Bigr\rangle=2^{-n/2}U^{1/2}\rho_n .$$
Not re-derived: the constant $(8\pi^{H+1/2})^{-n}$ and the sign $s_n$ (vertex continuation, extraction, reordering). No end-to-end position-space integral was done (as in the explorer's report).

## Q1-zero — confirmed

$n=1$: $-\tfrac12dM\,M^{-1}dM=-\tfrac1{2M}dM\wedge dM=0$, $\mathrm{Pf}$ of a $1\times1$ antisymmetric matrix $=0$; $\rho_n=0$ for odd $n$ (S3 prop-rho, cor-odd). `logs/kernels_verify.log`: `H=1,2,3 polygon: dM ^ dM = {}`.

## graphs-Q2 — confirmed

Fresh enumeration: subdivisions of theta, dumbbell, figure-eight cores with $|I|=2D+1$, $|V|=2D$, every edge subset tested against $(D-1)|I_S|\le D|V_S|-D-1$; graphs with a pendant vertex tested separately (none pass). Admissibility: $\Gamma\setminus\mathcal O$ a forest, $\deg\mathcal O=2+p$ (S2 eq-theta-rule; the derivation uses only $\delta^k(0)=0$, valid for $k\ge1$). Outputs $=3(|V|-1)-2(|I|-\deg\mathcal O)-\deg\mathcal O$.

`logs/graphs_verify.log`:
```
D=2: ... satisfying the counting: [('theta', (1, 2, 2))]; with a pendant vertex: 0
   theta(1, 2, 2): admissible O (vertex, deg O, p, |Aut fixing O|, outputs): [(0, 3, 1, 2, 2), (1, 3, 1, 2, 2)]
D=3: ... [('theta', (1, 3, 3)), ('theta', (2, 2, 3))]; with a pendant vertex: 0
   theta(1, 3, 3): [(0, 3, 1, 2, 4), (1, 3, 1, 2, 4)]
   theta(2, 2, 3): [(0, 3, 1, 2, 4), (1, 3, 1, 2, 4)]
D=4: ... [('theta', (1, 4, 4)), ('theta', (2, 3, 4)), ('theta', (3, 3, 3))]; with a pendant vertex: 0
   theta(1, 4, 4): [(0, 3, 1, 2, 6), (1, 3, 1, 2, 6)]
   theta(2, 3, 4): [(0, 3, 1, 1, 6), (1, 3, 1, 1, 6)]
   theta(3, 3, 3): [(0, 3, 1, 6, 6), (1, 3, 1, 6, 6)]
```

## theta-kernels — confirmed

(1) Pointwise, exact, 8 random rational points per graph, sources at all vertices including $N'$:
$$D(t(w,x))\,\bigl\lvert\det[t;\partial_\xi t]\bigr\rvert=s_\Gamma\,\frac{D_\Gamma}2\prod_Pw_P^{H+\frac12-n_P}\quad(\text{squares equal, one sign}).$$
Negative control: $w_B^{+1}$ inserted changes the squared ratio (1089, 2304/625, …).
(2) $D_\Gamma$ equal to the explicit determinants of the claim at 5 random integer $\lambda$.
(3) $K(\lambda;0)/(\pi D_\Gamma)=\tfrac{s}{2}\,\frac{\prod_P\Gamma(H+\frac32-n_P)}{\Gamma(H+\frac32)}\prod_P\frac1{(n_P-1)!}$:
$$\tfrac12\cdot\tfrac{2}{5}\cdot\tfrac14=\tfrac1{20},\quad-\tfrac12\cdot\tfrac{2}{15}\cdot\tfrac12=-\tfrac1{30},\quad\tfrac12\cdot\tfrac27\cdot\tfrac1{36}=\tfrac1{252},\quad\tfrac12\cdot\tfrac2{35}\cdot\tfrac1{12}=\tfrac1{420},\quad\tfrac12\cdot\tfrac{2}{105}\cdot\tfrac18=\tfrac1{840},\quad-\tfrac12\cdot\tfrac23=-\tfrac13 .$$
(4) Displacements: $-\sum_ey_e\cdot z_e=\sum_v\lambda_v\cdot X_v$ with fresh currents; with $y^P_1=\ell_P=w_P(\lambda_N+\sum_Q\hat m_Q)-\hat m_P$, $\hat m_P=\sum_k\sigma^P_k\lambda^P_k$:
$$-\sum_Py^P_1\cdot z_P=\sum_{P,k}\lambda^P_k\cdot\bigl[\sigma^P_kz_P+(1-\sigma^P_k)\bar Z\bigr]+\lambda_{N'}\cdot\bar Z,\qquad\bar Z=\sum_Pw_Pz_P .$$

`logs/kernels_verify.log`:
```
5d (A) S2 data   n_P=[1, 3, 3]: ... single sign s=+1 (claim +1: True); D_Gamma explicit: True; K(lam;0)/(pi D) = 1/20 (claim 1/20: True); X rule (4 points): True
5d (B) S2 data   n_P=[2, 2, 3]: ... s=-1 (claim -1: True); D_Gamma explicit: True; K(lam;0)/(pi D) = -1/30 (claim -1/30: True); X rule: True
5d theta_rootN(1, 3, 3): s=+1; 1/20;  theta_rootN(2, 2, 3): s=-1; -1/30
7d theta(1,4,4): s=+1 (True); True; 1/252 (True); True
7d theta(2,3,4): s=+1 (True); True; 1/420 (True); True
7d theta(3,3,3): s=+1 (True); True; 1/840 (True); True
3d theta(1,2,2): s=-1 (True); True; -1/3 (True); True
explicit X maps of the claim (5d):  (A): True   (B): True     (6 rational points each)
```
Direct numerical integrals of the brute-force density (float Kirchhoff solve and forms; agreement with the exact density at a rational point to 14–15 digits):

(a) chart-free, $K=\int_{\mathbb R^{N-1}}D(e^s)\prod_et_e\,d^{N-1}s$, scrambled Sobol, $s=5\,\mathrm{logit}(u)$, $|s_e|\le30$ (`logs/numeric_verify.log`):
```
case 122: K(lam;0)        QMC -3.1415529 +- 6.9e-05 (2^22 x 8)  exact -3.141592654  rel. diff 1.3e-05 (4.9 digits)
          int D y^0_{e=1} QMC 5.9689229 +- 0.00011             exact 5.969026042   rel. diff 1.7e-05 (4.8 digits)
```
For $N=7,9$ this chart gives no usable result: points with $t$-ratios up to $e^{60}$ give non-finite or garbage values (`case A: QMC -1.37e+14 +- 1.4e+14`; 7d: singular Laplacian).

(b) chart $t^P_j=x^P_j/w_P$, integrand $D(t)\lvert\det[t;\partial_\xi t]\rvert$ evaluated by brute force, $w\sim\mathrm{Dir}(\frac12,\frac12,\frac12)$ importance sampling, $2^{20}\times8$ points (`logs/numeric_wx.log`):
```
122  K(lam;0)  -3.141592649  exact -3.141592654  (8.9 digits);  int D y^0_e  5.969026037  exact 5.969026042  (9.1 digits)
A    K(lam;0)  10.05309647   exact 10.05309649   (8.6 digits);  int D y^0_e -26.32953818  exact -26.32953843 (8.0 digits)
B    K(lam;0)  -6.702064321  exact -6.702064328  (9.0 digits);  int D y^0_e  10.05309763  exact 10.05309649  (6.9 digits)
144  K(lam;0)  0.8975978981  exact 0.897597901   (8.5 digits);  int D y^0_e -4.463056224  exact -4.46305623  (8.9 digits)
234  K(lam;0)  0.5385587401  exact 0.5385587406  (9.0 digits);  int D y^0_e -2.617994322  exact -2.617993878 (6.8 digits)
333  K(lam;0)  0.1009797639  exact 0.1009797639  (10.0 digits); int D y^0_e -0.2804993379  exact -0.2804993441 (7.7 digits)
```
Exact insertion values: $\int D\,y^0_{e}=-\mathrm{sg}_e\int s\frac{D_\Gamma}2\prod w^{\alpha}\sum_v\lambda^0_vX^0_v(z_B=e_0)$ by Dirichlet integrals. Taylor orders in $z$: implied by (1) and (4); not separately compared with an engine.

## Q2-5d-operator — unverified

Consistent with the confirmed pieces: $J_\Gamma=-s_\Gamma\frac{D_\Gamma}2\prod w^\alpha$ ($I=-K$), $J_A=\tfrac12\langle12\rangle\langle45\rangle w_A^{3/2}w_B^{-1/2}w_C^{-1/2}$, $J_B=-\tfrac12\langle uw\rangle\langle ab\rangle(w_Aw_B)^{1/2}w_C^{-1/2}$, $\int J_A=\frac\pi{20}\langle12\rangle\langle45\rangle$, $\int J_B=-\frac\pi{30}\langle uw\rangle\langle ab\rangle$, the X maps, and $\frac{i}{128\pi^5}\cdot\frac{\pi}{20}=\frac{i}{128\pi^4}\cdot\frac1{20}$. $\deg_\lambda D_\Gamma=4$ with one $\lambda^1$ and one $\lambda^2$ in each $2\times2$ determinant, and $\sum\lambda\cdot X$ moves derivatives without changing their number: two $\partial_1$, two $\partial_2$. The prefactor follows from $Q_2=i^5\sum_\Gamma G_\Gamma W_\Gamma/|\mathrm{Aut}\,\Gamma|$ (S2 eq-Qn) with $G=-K/64\pi^5=I/64\pi^5$, $|\mathrm{Aut}|=2$: $Q_2=\frac{i}{128\pi^5}\sum_\Gamma I_\Gamma W_\Gamma$. It rests on raw-normalization (unverified) and on S2's $W_\Gamma$ Koszul sign (S2 §8.2 item 5, fixed there by nilpotency, not derived).

## nilpotency-5d — confirmed

Independent implementation: superfield $C^A=\sum_S\theta_Sc^A_S$; $\hat Q_0c^A_S=(-1)^{|S|}[\theta_S]\tfrac12f_{ABC}C^BC^C$, extended by $[\hat Q_0,\partial_a]=0$; $\delta/\delta C^A(\theta)=\sum_SF_S(\theta)\partial_L/\partial c^A_S$ with $F_S$ obtained by expanding $\prod_i(\theta_i-\theta'_i)$; $\hat Q_\Gamma\mathcal O=\int d^k\theta\sum_{l_1l_2l_3}G_\Gamma\,\delta_{l_1}\delta_{l_2}\delta_{l_3}\mathcal O$, $G_\Gamma=\sum_m\partial_z^{j}K_\Gamma/\pi\big|_{\lambda^m}T_\Gamma(m)$, weights from the pointwise-verified density and X maps by my own Dirichlet code, $T_A,T_B$ of S2 eq-res-TA, eq-res-TB, $\mathfrak{so}(N)$ constants from matrix commutators (antisymmetry and Jacobi asserted). Only the $\theta$-degree $\sum_j|S_j|-2k$ of $G$ is built.

Weights at $z=0$ (`logs/nilp_5d_so3.log`): A: $\pm\frac1{20}$, B: $\pm\frac1{30}$ on the four monomials $\lambda_1^{a}\lambda_2^{b}\lambda_4^{c}\lambda_5^{d}$, equal to the explorer's.

```
5d so3, 40 random operators (seed 11, 1-3 input derivatives): 14 with Q_G O != 0, 12 with nonzero anticommutators; 798 equations, rank 1, null space [1, 1]
5d so3 control (z-derivative weights doubled), same operators: 798 equations, rank 2, null space []
   (single operators give [-37/40, 1], as in the explorer's control_mutate_5d.log)
5d so3, 7 operators without input derivatives: all anticommutators vanish (5 with Q_G O != 0)
5d so5, 5 operators without input derivatives: all anticommutators vanish (4 with Q_G O != 0)
5d so5 fixed operators:
   d2 g1_0 g2_1 b_2       248 eqs  rank 1  null [1, 1]
   g1_0 d1 g2_1 b_2       248 eqs  rank 1  null [1, 1]
   d1 b_0 b_1 g1_2       1240 eqs  rank 1  null [1, 1]
   b_0 d2 b_1 c_2         248 eqs  rank 1  null [1, 1]
   g1_0 g2_1 d1d2 b_2    2344 eqs  rank 1  null [1, 1]
   b_0 d1 g1_1 d2 b_2   11680 eqs  rank 1  null [1, 1]
   d1 g1_0 d2 g2_1 b_2   2336 eqs  rank 1  null [1, 1]    (explorer: 2336)
   total 18344 equations, rank 1, null space [1, 1]
5d so5, 12 random operators (seed 21), 4 with nonzero anticommutators:
   0_4^01 1_5^02 01_6^00            5458 eqs  rank 1  null [1, 1]
   01_0^10 1_3^01 01_6^01 0_2^00   60325 eqs  rank 1  null [1, 1]
   1_6^00 01_9^10 0_2^10 01_3^00   10851 eqs  rank 1  null [1, 1]
   01_9^12 0_3^00 1_9^00 1_0^00    17155 eqs  rank 1  null [1, 1]
   total 93789 equations, rank 1, null space [1, 1]
5dN (theta_rootN data, arm tensor): so3 40 random 798 eqs null [1, 1]; so5 3 operators 3832 eqs null [1, 1]
Q0hat^2 = 0 on all random operators tested: True
```
(Tokens in the logs: `<S>_<colour>^<jet>`, $S$: `e` = c, `0` = $\gamma^1$, `1` = $\gamma^2$, `01` = b.) The two operators the explorer stopped ($\gamma^1\gamma^2\partial_1\partial_2\mathsf b$, $\mathsf b\,\partial_1\gamma^1\partial_2\mathsf b$) are included.

## 7d-relative-coefficients — confirmed

Same implementation, $k=1$, $H=3$, theta_rootN data, arm tensor, outputs ordered $v^A_1..,v^B_1..,v^C_1..$, spoke $j$ on the first edge of arm $j$.
```
7d so3, 6 random operators (seed 3): 2 with nonzero anticommutators
   e_2^010 0_0^000 0_0^110 e_1^000   404 eqs, rank 2, null [-3, 6, 1]
   e_2^011 0_2^000 0_2^000            46 eqs, rank 2, null [-3, 6, 1]
   total 450 equations, rank 2, null space [-3, 6, 1]
control (z-derivative weights doubled), same operators: 450 equations, rank 3, null space []
7d so3, 8 random operators (seed 4): 4 with nonzero anticommutators (14, 18, 14, 34 eqs); 80 equations, rank 2, null space [-3, 6, 1]
```
(Tokens for $k=1$: `e` = $c$, `0` = $\theta$-component.)
$(-3,6,1)=6\,(-\tfrac12,1,\tfrac16)$; $(-1)^{m_Bm_C}$: $m=(0,3,3)\to-1$, $(1,2,3)\to+1$, $(2,2,2)\to+1$; 5d theta_rootN $(0,2,2),(1,1,2)\to(+1,+1)$, consistent with the 5dN null space $[1,1]$ and $|\mathrm{Aut}|=(2,2)$. $k_\Gamma K_\Gamma/\pi$: $-\tfrac12\cdot\tfrac1{252}=-\tfrac1{504}$, $\tfrac1{420}$, $\tfrac16\cdot\tfrac1{840}=\tfrac1{5040}$. The sign rule remains an observation on 4 graphs; so(5) not run.

## 3d-closure — confirmed

```
3d so3, 1500 random operators: 269 with Q_G O != 0: all anticommutators vanish
3d so5,  600 random operators:  92 with Q_G O != 0: all anticommutators vanish
controls (z-derivative weights doubled): so3 43 operators fail, 257 equations, rank 1, null space []
                                         so5 64 operators fail, 1585 equations, rank 1, null space []
```
$K=-\frac\pi3\lambda_{B_1}\lambda_{C_1}$ confirmed above. Magnitude $\frac{1}{64\pi^3}\cdot\frac\pi3\cdot\frac12=\frac1{384\pi^2}$: arithmetic only; it rests on raw-normalization and on the assembly rule for odd $P$ (unverified).

## S2-Bnorm-sign — unverified

Requires the position-space sign bookkeeping that was not re-derived (see raw-normalization). S2 eq-B-norm attributes $(-1)^n=(-1)^{3(k-1)(m-1)}$ to the ordering of the Schwinger differentials; the explorer's extra $(-1)^{\binom{K'}2}$ was not independently computed.

## Logical gaps

1. nilpotency-5d: $Q_1=0$ makes $\{Q_0,Q_2\}=0$ homogeneous; it fixes $k_A:k_B$ only, not the normalization nor the sign of $Q_2$ relative to $Q_0$.
2. nilpotency-5d: "(1,1) equals the kernel prediction $I_AW_A+I_BW_B$" holds with S2's $T_A,T_B$ and S2's relative Koszul sign of $W_A,W_B$; in S2 that sign is fixed by $Q_1^2+\{Q_0,Q_2\}=0$ (S2 §8.2 item 5), not derived, and in 5d $\Omega$ is even. The agreement is a consistency statement, not a derivation of the 5d $W_\Gamma$.
3. nilpotency-5d, 7d, 3d: closure on finitely many operators; $Q_0$-exactness of $\sum k_\Gamma\hat Q_\Gamma$ (triviality in cohomology) not addressed.
4. 7d-relative-coefficients: the coefficients refer to the arm convention (vertex labels of theta_rootN, output order, spoke placement); the physical relative signs need the Koszul signs of $W_\Gamma$ for odd $P$, not derived. The rule $(-1)^{m_Bm_C}$ is fitted to 4 graphs. The statement "$|k_\Gamma|=1/|\mathrm{Aut}\,\Gamma|$" is true up to the common factor 6 of the null vector.
5. 3d-closure: closure of a single graph is checked; the normalization $1/(384\pi^2)$ and the overall sign are not.
6. Q2-5d-operator: built on raw-normalization (no end-to-end position-space integral) and on S2's assembly; overall sign relative to $Q_0$ not checked (stated by the explorer).
7. raw-normalization / S2-Bnorm-sign: the sign $s_n$ and the $(-1)^n$ discrepancy with S2 eq-B-norm come from the same unverified bookkeeping; at $n=2$ the two agree, so no number in the 5d result depends on the discrepancy.
8. action-hbar: for odd $k$ the placement of $\int d^k\theta$ (odd) relative to $\Omega$ is a convention; for the three cases computed the reordering sign is $+1$, so $-\tfrac14$ does not depend on it.
9. theta-kernels: the general formula is checked on the six graphs used (plus two relabellings); $s_\Gamma$ is determined per graph by the brute force, not by a formula.
10. twist-5d: the remark on the non-basic edge form in the pull-back convention has no check in the explorer's logs or here.
