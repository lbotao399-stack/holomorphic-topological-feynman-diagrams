# nonSP: adversarial verification

Scripts and logs: `explore/nonSP/verify/`. No explorer script is run; the explorer's `graphs.py` is imported for its edge lists only (`K4SUB`, `PRISM`, `K33`; identical to `results/graph_counts.txt` lines 5, 6, 26–31). The C^2 $n=3$ graph is the edge list of claim b5. `htk` is imported only in `t_conv.py` (convention check of the independent layer against `htk.brute.density`).

`explore/nonSP/REPORT.md` does not exist (directory listing of `explore/nonSP/`). Section 5 is based on the structured report and the explorer scripts.

## 1. Independent layer (`vk.py`)

$$B_{ev}=\delta_{v,e(0)}-\delta_{v,e(1)}\ (v<\text{root}=nv-1),\qquad \Lambda=B^{\mathsf T}t^{-1}B,\quad \varphi=\Lambda^{-1}\lambda,\quad y=t^{-1}B\varphi ,$$

$C$: fundamental cycles of a BFS tree, first row flipped until $\det\binom{B^{\mathsf T}}{C}>0$;

$$M=CtC^{\mathsf T},\quad U=\det M,\quad \Psi_{(i,r),e}=C_{re}y^i_e\ (i\text{ outer}),\qquad D=\frac{\det[\mathbf 1;\Psi]}{U^H\sum_et_e},\qquad K=\int_{t>0,\ \sum t=1}D\,dt_0\cdots dt_{N-2}.$$

($\det[a;\Psi]=\sum_e(-1)^ea_e\,\mathrm{minor}_e$ and $U^{-H}\prod\psi=D\,\iota_E\mathrm{vol}$ give $\mathrm{minor}_e=(-1)^eU^HDt_e$, hence $\det[a;\Psi]=U^HD\,(a\cdot t)$.)

`t_conv.py` (`out_t_conv.log`): `density_exact == htk.brute.density` (exact rationals) at 3 random rational $t$ with random rational sources for each of the 9 graphs (6 K4 subdivisions, prism, K33, C^2 n=3 graph): True in all 27 cases; numpy `NG.D` relative deviation $\le4.3\cdot10^{-14}$.

## 2. Exact pointwise identities

### 2.1 a0: currents (`check_a0_a1.py`, `out_check_a0_a1.log`)

Formulas of a0 re-implemented from the claim text; chart $\xi=(s_{E,k};g_A,g_B,r_A,r_B,r_C)$, $g_C=1-g_A-g_B$, $R_{D_cv}=1/g_v$, $R_{wx}=(1-r_v)/(r_vg_wg_x)$, $t_{E,j}=R_E(s_{E,j-1}-s_{E,j})$. Compared with `currents_exact` (full graph) at 3 random rational $\xi$, every edge and component:

```
6 graphs x 4 centres, sources e1,e2,e3,-(sum) at the degree-two vertices: currents(formula) == Kirchhoff: True (24/24)
K4[1,1,1,1,3,3] Dc=1,2 and K4[1,1,1,2,2,3]a Dc=1,2, generic integer sources at all 8 vertices: True (4/4)
C^2 n=3 graph, Dc = 0,1,2,4, sources lam3=e1, lam0=e2, lam2=(1/3,2/5), lam1=0: True (4/4)
```

### 2.2 a1: density

$g_{\rm ref}=D(t(\xi))\,|\det[t;\partial_\xi t]|$ (sympy Jacobian of the map, `density_exact`), $g_f=\det\begin{pmatrix}Q&\partial Q\\N&\partial N\end{pmatrix}/Q^{10}$, $N=Q\,y$ on the first edges of the three triangle chains (chain direction), rows ($i$ outer, chain inner). Schur complement: $\det\begin{pmatrix}Q&\partial Q\\N&\partial N\end{pmatrix}=Q\det(\partial N-N\partial Q/Q)=Q^{10}\det\partial_\xi(N/Q)$, i.e. $g_f=\det\partial_\xi\ell$, $\ell=N/Q$ the cotree currents.

```
g_ref/g_f at 3 points per (graph, centre):
K4[1,1,1,1,3,3]  Dc=1: -1,-1,-1  Dc=2: -1,-1,-1  Dc=4: 1,1,1  Dc=5: 1,1,1
K4[1,1,1,2,2,3]a Dc=1: 1,1,1     Dc=2: 1,1,1     Dc=4: 1,1,1  Dc=5: 1,1,1
K4[1,1,1,2,2,3]b Dc=1: -1,-1,-1  Dc=2: -1,-1,-1  Dc=4: -1,-1,-1  Dc=6: -1,-1,-1
K4[1,1,1,2,2,3]c Dc=2: 1,1,1     Dc=3: 1,1,1     Dc=4: 1,1,1  Dc=5: -1,-1,-1
K4[1,1,2,2,2,2]a Dc=0: 1,1,1     Dc=2: -1,-1,-1  Dc=3: 1,1,1  Dc=5: 1,1,1
K4[1,1,2,2,2,2]b Dc=0: 1,1,1     Dc=1: -1,-1,-1  Dc=3: 1,1,1  Dc=4: 1,1,1
```

(The sign depends on the row order of $N$ and the order of $\xi$; it is constant per (graph, centre), as claimed.)

### 2.3 Automorphisms (`check_sym.py`, `out_check_sym.log`)

$(Pt)_{\pi(e)}=t_e$; Lebesgue measure on the simplex is $P$-invariant, so $D(Pt,\lambda')=\sigma D(t,\lambda)\ \forall t\Rightarrow K(\lambda')=\sigma K(\lambda)$. Exact ratios at 4 random rational $t$:

```
a7 K4[1,1,2,2,2,2]a pi=(3 5)(4 7): D(Pt,lam)/D(t,lam) = ['-1', '-1', '-1', '-1']
a7 generic sources, lam' = lam o pi^-1: D(Pt,lam')/D(t,lam) = ['1', '1', '1', '1']
b3 K33 src {0,1,3,5}, pi=(2 4): ['-1', '-1', '-1', '-1']
b3 K33 src {0,1,3,5} (other lam), pi=(2 4): ['-1', '-1', '-1', '-1']
b3 K33 src {0,1,2}, pi=(3 5): ['-1', '-1', '-1', '-1']
b3 K33 src {1,2,3,5} (+root 5 implied), pi=(0 4): ['-1', '-1', '-1', '-1']
b5 src {3,2,4}, pi=(0 1): ['-1', '-1', '-1', '-1']
b2 prism: D(Pt, src{0,1,3})/D(t, src{0,1,2}), pi=(2 3)(4 5): ['-1', '-1', '-1', '-1']
```

a7: with $\lambda=(e_1,e_2,e_3)$ at $(1,4,6)$ the ratio is $-1$ with the same sources, hence $K=-K=0$. (Equivalently: $\sigma=+1$ for $\lambda'=\lambda\circ\pi^{-1}$, and $\lambda'=g\lambda$ with $\det g=-1$, $D(t,g\lambda)=\det(g)^3D(t,\lambda)$.)
b2: $K_{\{0,1,2\}}=-K_{\{0,1,3\}}=-K_{\rm b1}(\lambda_2=0)=-c_{400}$.

### 2.4 Covariance (`check_cov.py`, `out_check_cov.log`)

$D(t,g\lambda)=\det(g)^nD(t,\lambda)$, exact, 2 random $(t,\lambda,g)$ each for K4[1,1,1,1,3,3], K4[1,1,1,2,2,3]c, prism, K33: True (8/8). Hence $K(g\lambda)=\det(g)^nK(\lambda)$: on C^3 with sources at three independent vertices $K=c\,\Delta^3$ (degree 9, $SL_3$-invariant); on C^2 $K$ is a polynomial in the brackets.

### 2.5 General formulas of a2 and the weight statement of c2 (`check_a2_c2.py`, `out_check_a2_c2.log`, mpmath, 40 digits)

$$\int_{\Delta^2}\frac{g^\beta\,dg_Adg_B}{(\rho\cdot g)^{10}}=\frac{\beta!}{9!}\prod_v\rho_v^{-\beta_v-1}\ (|\beta|=7),\qquad \frac1{Q^{10}G^b}=\frac{\Gamma(10+b)}{\Gamma(10)\Gamma(b)}\int_0^1\frac{\xi^9(1-\xi)^{b-1}d\xi}{(\sum_v(1-\xi r_v)g_v)^{10+b}}$$

```
Feynman formula beta=(3, 2, 2): lhs/rhs - 1 = -1.15e-41
Feynman formula beta=(7, 0, 0): lhs/rhs - 1 = 0.0
Feynman formula beta=(1, 4, 2): lhs/rhs - 1 = 0.0
Feynman parametrization b=1..4: rel dev -4.59e-41
int xi^-1 log^1(1-xi): [v,1,z2,z3,z4] relation [1, 0, 1, 0, 0]
int xi^-1 log^2(1-xi): [1, 0, 0, -2, 0]
int xi^-2 log^2(1-xi): [1, 0, -2, 0, 0]
int xi^-1 log^3(1-xi): [-1, 0, 0, 0, -6]
int xi^-2 log^3(1-xi): [-1, 0, 0, -6, 0]
int xi^-3 log^3(1-xi): [1, 0, 3, 3, 0]
```

and $\int_0^1\xi^m\log^k(1-\xi)d\xi=\sum_{j\le m}\binom mj(-1)^{j+k}k!/(j+1)^{k+1}\in\mathbb Q$ for $m\ge0$.

## 3. Direct numerical integration of the Schwinger integral

### 3.1 Method (`hepp.py`, `lattice2.py`)

K4 subdivision, chain $E$ with $a_E$ edges: $t_{E,j}=R_Ex_{E,j}$, $\sum_jx_{E,j}=1$, $d^Nt=\prod_ER_E^{a_E-1}dR_E\,d^{a_E-1}x_E$,

$$K=\int_{\mathbb R^5}F(s)\,d^5s,\qquad F(s)=\prod_ER_E^{a_E}\int_{\prod\Delta}D(t(R,x))\,dx,\qquad R=(e^{s_1},\dots,e^{s_5},1).$$

$D$: `vk.NG.D` on the full graph (no reduction formula). $x$-integral: Gauss–Legendre product rule ($2$-simplices by $x=(u,(1-u)v,(1-u)(1-v))$), $m=2$; `xdeg.py`: $m=2,3,4$ equal $m=6$ to $\le6.6\cdot10^{-14}$ relative at 4 random $R$ per graph, i.e. exact. $s$-integral: trapezoidal rule on $h\mathbb Z^5$ and $h\mathbb Z^5+h/2$, flood fill from $0$ through $\max_S|F_S|>\mathrm{tol}\cdot\max|F|$, every evaluated point summed.

Second quadrature (`hepp.py`, `sector_integral`): the 720 Hepp sectors $R_{\sigma_0}=1\ge R_{\sigma_1}\ge\dots$, $R_{\sigma_k}=u_0\cdots u_{k-1}$, $d^5R=u_0^4u_1^3u_2^2u_3\,d^5u$, tensor Gauss–Legendre $n^5$.

Calibration on the C^2 n=3 graph, sources $\lambda_3=e_1$, $\lambda_0=e_2$, $\lambda_2=-(e_1+e_2)$ (both quadratures; logs `out_bench1_hepp_c2n3.log`, `out_bench2_hepp_c2n3.log`, `out_bench5_lattice_c2n3.log`):

```
Hepp n=3,4,5,6,8:  -1.714303533e-3, -1.732651750e-3, -1.735532413e-3, -1.735983121e-3, -1.7361029671e-3
lattice h=0.80 tol=1e-9:  K0=-1.7361109148812e-03 K1=-1.7361109018561e-03  (pts 6419199)
lattice h=0.70 tol=1e-10: K0=-1.7361110794182e-03 K1=-1.7361110800775e-03  (pts 17539671)
-1/576 = -1.7361111111111e-03:  rel dev (lattice) 1.1e-7 (h=0.8), 1.8e-8 (h=0.7)
```

At $(h,\mathrm{tol})=(0.8,10^{-9})$ the two shifts differ by $7.5\cdot10^{-9}$ relative, the error is $1.1\cdot10^{-7}$ relative: the shift difference is not an error bound; the error is read off from the calibration and from the $h=0.7$ run in §3.2.

### 3.2 C^3, $n=3$: $c_\Gamma=K(\lambda;0)$ at $\lambda_{d_1,d_2,d_3}=e_1,e_2,e_3$ (`run_c3.py`, `out_c3_A.log`, `out_c3_B.log`)

$h=0.8$, $\mathrm{tol}=10^{-9}$:

| graph | lattice, shift 0 | lattice, shift $h/2$ | claim | rel. dev. |
|---|---|---|---|---|
| K4[1,1,1,1,3,3] | 7.949092625444e-06 | 7.949092918987e-06 | $(39-32\zeta(3))/67200$ = 7.949093733527e-06 | $-1.4\cdot10^{-7}$, $-1.0\cdot10^{-7}$ |
| K4[1,1,1,2,2,3]a | -1.286007948236e-05 | -1.286007973335e-05 | $-1/77760$ = -1.286008230453e-05 | $-2.2\cdot10^{-7}$, $-2.0\cdot10^{-7}$ |
| K4[1,1,1,2,2,3]b | -7.189565615267e-06 | -7.189565738186e-06 | $(1296\zeta(3)-1597)/5443200$ = -7.189567442895e-06 | $-2.5\cdot10^{-7}$, $-2.4\cdot10^{-7}$ |
| K4[1,1,1,2,2,3]c | 3.404901849036e-06 | 3.404901767675e-06 | $(34-27\zeta(3))/453600$ = 3.404902148790e-06 | $-8.8\cdot10^{-8}$, $-1.1\cdot10^{-7}$ |
| K4[1,1,2,2,2,2]a | -1.467421640293e-20 | -1.655046543488e-20 | 0 | (scale $10^{-5}$) |
| K4[1,1,2,2,2,2]b | 7.949094500577e-06 | 7.949093622488e-06 | $(39-32\zeta(3))/67200$ | $+9.7\cdot10^{-8}$, $-1.4\cdot10^{-8}$ |

$h=0.7$, $\mathrm{tol}=10^{-10}$ (`out_c3_fine2.log`):

```
K4[1,1,1,2,2,3]b h=0.70 tol=1e-10  K(shift0)=-7.189567307242e-06  K(shift h/2)=-7.189567313254e-06  claim=-7.189567442895e-06  rel dev -1.89e-08 -1.80e-08  pts 16981218
```

(deviation $2.5\cdot10^{-7}\to1.9\cdot10^{-8}$, as in the calibration $1.1\cdot10^{-7}\to1.8\cdot10^{-8}$.)

Hepp-sector Gauss (`run_hepp_c3.py`, `out_hepp_c3.log`):

```
K4[1,1,1,2,2,3]a Hepp n=4: K = -1.2798260927e-05   sectors >0: 0 <0: 720
K4[1,1,1,2,2,3]a Hepp n=5: K = -1.2847807031e-05   sectors >0: 0 <0: 720
K4[1,1,1,1,3,3]  Hepp n=5: K =  7.9620495901e-06   sectors >0: 720 <0: 0
```

(relative deviation from the claim at $n=5$: $9.5\cdot10^{-4}$ (a4), $1.6\cdot10^{-3}$ (a3).)

### 3.3 b5: C^2, $n=3$ (`run_c2n3.py`, `out_c2n3.log`)

$\lambda_3=e_1$, $\lambda_0=e_2$, $\lambda_2=(a,b)$, $\lambda_1=0$: $p=1$, $q=b$, $s=-a$; 16 configurations $(a,b)\in\{-1,0,1,2\}^2$ in one lattice run ($h=0.8$, tol $=10^{-9}$, 11183157 points); cubic form in $(p,q,s)$ by least squares (max residual $9.7\cdot10^{-17}$):

```
                 576*coef (shift 0)   576*coef (shift h/2)   claim
p^3              -0.999999815         -0.999999809           -1
p^2 q            -0.999999857         -0.999999854           -1
p^2 s            +2.999999473         +2.999999445           +3
p q^2            +4.999999362         +4.999999304           +5
p q s            +11.999998616        +11.999998491          +12
p s^2            +2.999999869         +2.999999838           +3
q^3, q^2 s, q s^2, s^3   +7.0e-8, +1.8e-7, +1.2e-8, -9.6e-15                 0
K(a,b)=(0,0):  -1.736110790405e-03   (claim -1/576 = -1.736111111111e-03)
```

i.e. $K=-(p^3+p^2q-3p^2s-5pq^2-12pqs-3ps^2)/576$ to $\le1.5\cdot10^{-6}$ in $576\cdot$coefficient. The configuration $(a,b)=(-1,-1)$ is sources $\{3,0,2\}$: $K=-p^3/576$ from the same fit.

### 3.4 b1, b2: prism, C^2, $n=4$ (`prism.py`, `out_prism.log`)

Chart built here (bijection of positive orthants): K4 $\{0,1,2,3\}$ resistances $R_{01},R_{02},R_{03}=1,R_{12},R_{13},R_{23}$ and fractions $f_{13},f_{23},f_{15}\in(0,1)$;
$t_{23}=R_{23}/f_{23}$, $r''_{23}=R_{23}/(1-f_{23})$, $r_{13}=R_{13}/f_{13}$, $r''_{13}=R_{13}/(1-f_{13})$, $r''_{12}=R_{12}$; $\Delta\to Y$ on $\{1,2,3\}$ with centre 5: $\rho_1=r''_{12}r''_{13}/S_2$, $\rho_2=r''_{12}r''_{23}/S_2$, $\rho_3=r''_{13}r''_{23}/S_2$; $t_{25}=\rho_2$, $t_{15}=\rho_1/f_{15}$, $r'_{15}=\rho_1/(1-f_{15})$, $r_{35}=\rho_3$; $\Delta\to Y$ on $\{1,3,5\}$ with centre 4: $t_{14}=r_{13}r'_{15}/S_1$, $t_{34}=r_{13}r_{35}/S_1$, $t_{45}=r'_{15}r_{35}/S_1$.
$K=\int D(t(\xi))|\det[t;\partial_\xi t]|d\xi$, log-lattice in the five $\log R$, Gauss–Legendre in the fractions; $m=2$ equals $m=3,4,6$ to $\le10^{-15}$ (`out_t_prism.log`, exact). Sources $\lambda_0=e_1$, $\lambda_1=e_2$, $\lambda_2=(S_A,S_B)$, $\lambda_3=-(\lambda_0+\lambda_1+\lambda_2)$, $\lambda_4=\lambda_5=0$; 19 configurations; quartic fit in $(S_A,S_B)$ (residual $\le10^{-18}$); $c_{4-j-k,j,k}=(-1)^k[S_A^kS_B^j]K$. $h=0.8$, tol $=10^{-9}$, 5285733 points:

```
  c_400  shift0 -2.8935180898e-05  shift h/2 -2.8935182581e-05  claim -2.8935185185e-05  diff 4.3e-12
  c_310  shift0 -5.6161426550e-05  shift h/2 -5.6161429611e-05  claim -5.6161436216e-05  diff 9.7e-12
  c_301  shift0 +9.0223417310e-05  shift h/2 +9.0223418129e-05  claim +9.0223423863e-05  diff -6.6e-12
  c_220  shift0 +5.6161438680e-05  shift h/2 +5.6161438787e-05  claim +5.6161436216e-05  diff 2.5e-12
  c_211  shift0 +2.7234029153e-04  shift h/2 +2.7234029227e-04  claim +2.7234030727e-04  diff -1.6e-11
  c_202  shift0 -4.9364592210e-05  shift h/2 -4.9364593226e-05  claim -4.9364598078e-05  diff 5.9e-12
  c_130  shift0 -8.6805547529e-05  shift h/2 -8.6805548864e-05  claim -8.6805555556e-05  diff 8.0e-12
  c_121  shift0 -3.7273951990e-04  shift h/2 -3.7273951747e-04  claim -3.7273953910e-04  diff 1.9e-11
  c_112  shift0 -6.4507980644e-04  shift h/2 -6.4507980018e-04  claim -6.4507984637e-04  diff 4.0e-11
  c_103  shift0 -1.8211686959e-04  shift h/2 -1.8211686636e-04  claim -1.8211688340e-04  diff 1.4e-11
  c_040, c_031, c_022, c_013, c_004 (shift 0): 1.0e-19, 3.9e-13, 1.5e-12, 1.3e-11, 3.4e-12        claim 0
```

Relative deviations $\le1.6\cdot10^{-7}$ on all ten nonzero coefficients. At $(S_A,S_B)=(1/3,2/5)$ the fitted polynomial gives $-1.21886646\cdot10^{-4}$ (claim $-1.2188666\cdot10^{-4}$). $K(\lambda_2=0)=c_{400}=-2.8935181\cdot10^{-5}$; with §2.3, $K_{\{0,1,2\}}=+2.8935181\cdot10^{-5}$ (claim $1/34560=2.8935185\cdot10^{-5}$).

### 3.5 b3, b4: K33 (`k33.py`, `run_k33.py`)

Chart: K4 $\{1,2,3,5\}$ ($R_{35}=1$), fractions $f_{13},f_{15},f_{35}$: $r_{ab}=R_{ab}/f_{ab}$ (triangle of the star at 4), $r'_{ab}=R_{ab}/(1-f_{ab})$ (triangle of the star at 0), $\Delta\to Y$ on both; $t_{12},t_{23},t_{25}=R$. The integrand is not polynomial in the fractions; quadrature $f=\phi(u)$, $\phi(u)=u^p/(u^p+(1-u)^p)$, Gauss–Legendre in $u$. Pointwise convergence (`out_t_k33b.log`, max over 6 random $R$, relative to $\max|F|$, reference $m=32,p=3$ vs $m=28,p=4$: $6.5\cdot10^{-8}$):

```
p=1 (plain Gauss, out_t_k33.log, 4 points): m=8 vs converged: +1.1%, +1.7%, +0.6%, +0.8% in magnitude
p=2 m=8 1.3e-03   m=10 2.8e-04   m=12 7.4e-05   m=16 8.2e-06
p=3 m=8 1.6e-02   m=10 2.3e-03   m=12 5.3e-04   m=16 2.0e-05
```

Configurations: [b4] $\lambda_0=e_1,\lambda_1=e_2,\lambda_2=(1/3,2/5),\lambda_3=-(\mathrm{sum}),\lambda_4=\lambda_5=0$; [b3a] sources $\{0,1,3,5\}$; [b3b] sources $\{0,1,2\}$.

```
out_k33_m10.log:        h=1.2 tol=1e-5  m=10 p=2  (168237 points)
  b4 src{0,1,2,3}: shift0 6.32467390e-05  shift h/2 6.32326280e-05
  b3 src{0,1,3,5}: shift0 3.20254421e-09  shift h/2 1.40582661e-08
  b3 src{0,1,2}:   shift0 -3.40326088e-10 shift h/2 -2.98302396e-10
out_k33_m10_coarse.log: h=1.4 tol=3e-5  m=10 p=2  (65571 points)
  b4 src{0,1,2,3}: shift0 6.32889617e-05  shift h/2 6.31721031e-05
  b3 src{0,1,3,5}: shift0 4.33730733e-08  shift h/2 1.21330422e-07
  b3 src{0,1,2}:   shift0 -5.08507082e-10 shift h/2 -3.98007390e-10
out_k33_m12_coarse.log: h=1.4 tol=3e-5  m=12 p=2  (65577 points)
  b4 src{0,1,2,3}: shift0 6.32918699e-05  shift h/2 6.31669087e-05
  b3 src{0,1,3,5}: shift0 4.42817511e-08  shift h/2 1.25900044e-07
  b3 src{0,1,2}:   shift0 -5.07900166e-10 shift h/2 -3.99046804e-10
```

$m=10\to12$ on the same lattice: $4.6\cdot10^{-5}$, $8\cdot10^{-5}$ relative; lattice: shift spread $2.2\cdot10^{-4}$ ($h=1.2$), $1.8\cdot10^{-3}$ ($h=1.4$); shift means $6.3240\cdot10^{-5}$ ($h=1.2$), $6.3229\cdot10^{-5}$ ($h=1.4$); boundary-layer mass $1.3\cdot10^{-7}$ ($h=1.2$), $3.0\cdot10^{-7}$ ($h=1.4$). Hence

$$K_{K33}(\lambda_0=e_1,\lambda_1=e_2,\lambda_2=(1/3,2/5),\lambda_3=-\textstyle\sum)=(6.324\pm0.005)\cdot10^{-5}.$$

The b3 configurations give $|K|\le1.3\cdot10^{-7}$, i.e. zero within the lattice error at the scale $6\cdot10^{-5}$ (exact zero: §2.3).

## 4. Verdicts

| id | verdict | evidence |
|---|---|---|
| a0 | confirmed | §2.1 |
| a1 | confirmed | §2.2 |
| a2 | unverified | general formulas confirmed (§2.5); the case statements (pole-freeness of $R_k$ at $\xi=1$, which centres give $b=0$, the exported $R_k$) were not re-derived |
| a3 | confirmed | §3.2: 6–7 digits; Hepp $n=5$: $1.6\cdot10^{-3}$ |
| a4 | confirmed | §3.2: 6–7 digits; Hepp $n=4,5$ |
| a5 | confirmed | §3.2: 6 digits ($h=0.8$), 7–8 digits ($h=0.7$, rel. dev. $1.9\cdot10^{-8}$) |
| a6 | confirmed | §3.2: 7 digits |
| a7 | confirmed | §2.3 exact antisymmetry; §3.2 $|K|<2\cdot10^{-20}$ |
| a8 | confirmed | §3.2: 7 digits; equals a3 within $2.4\cdot10^{-7}$ |
| b1 | confirmed | §3.4: all 15 coefficients, $\le1.6\cdot10^{-7}$ relative |
| b2 | confirmed | §2.3 exact sign relation + §3.4 $c_{400}$ |
| b3 | confirmed | §2.3 exact antisymmetry (4 configurations); §3.5 |
| b4 | refuted | §3.5: $K=(6.324\pm0.005)\cdot10^{-5}$; nonzero as claimed, but the claimed two digits "6.4e-5" are wrong (two digits: $6.3\cdot10^{-5}$); the explorer's $m=8$ value $6.3876\cdot10^{-5}$ is $1.0\%$ high, $m=6$ ($6.4188\cdot10^{-5}$) $1.5\%$ high |
| b5 | confirmed | §3.3: all 10 coefficients; §2.3 for $\{3,2,4\}$ |
| b6 | unverified | not computed; logic in §5 |
| c1 | confirmed | (i)–(iii) are a3–a8, b1–b5 at the stated precision; irrationality rests on the closed forms, which the independent numerics confirm to 6–7 digits, not on an independent exact derivation |
| c2 | unverified | the lemma (§2.5) holds; the form $\sum_kR_k\log^k(1-\xi)$ with $R_k\in\mathbb Q[\xi,1/\xi]$ was not re-derived |
| c3 | unverified | status open |

## 5. Logical gaps in the explorer's evidence

1. a3–a8, b1, b5: the exact values come from one Mathematica pipeline (P→H step: Feynman formula, $r$-integration, export). Its internal checks are the agreement of two centres and xi_cert.py, which certifies only the $\xi$-quadrature of the exported 1D integrand. The explorer's reduction-free checks reach 4–5 digits. (§3 closes this to 6–7 digits.)
2. pslq_test.py: basis $\{1,\zeta(3),\pi^2,\log2\}$ omits $\zeta(4)$ (= $\pi^4/90$), which c2 lists in the a priori span. Not consequential, since the values are exact closed forms from the pipeline.
3. a1: "$s_\Gamma$ constant" is checked at 3 points; constancy follows because the chart is a bijection of positive orthants ($r_v=h_v/(h_v+g_wg_x)$ inverts it), so $\det[t;\partial_\xi t]\ne0$ on the connected domain.
4. a7, b3, b5 $\{3,2,4\}$, b2: the sign bookkeeping (vertex sign, reversals, $\mathrm{sgn}^{D+1}$) is replaced in §2.3 by the exact pointwise identity $D(Pt)=\pm D(t)$; all signs agree with the explorer's.
5. b3 "every 3-source coefficient vanishes": needs that among the three source-free vertices two lie on one side (pigeonhole, 3 vertices, 2 sides); not stated.
6. b4: the fraction integrand in the chart (K4 $\{1,2,3,5\}$, three merge fractions) is not polynomial; plain Gauss–Legendre converges algebraically and from above (pointwise errors $+0.6\%$ to $+1.7\%$ in magnitude at $m=8$, `t_k33.py`, `out_t_k33.log`), so the $m=6,8$ sequence does not bracket the value and "about 2 digits" is not supported. Endpoint clustering $\phi(u)=u^2/(u^2+(1-u)^2)$ restores fast convergence (§3.5).
7. b6: the conclusion "cannot determine" needs only that the relations form a linear system over $\mathbb Q$ with rational right-hand side (a uniquely determined unknown is then rational). The hypothesis that all $n\le3$ inputs with insertions are rational is untested (only the $z=0$ four-source value of the C^2 $n=3$ graph is computed).
8. c1 (iii) "rational even through a sourced star center": the rationality of b5 is established by the exact pipeline only; §3.3 confirms the rational coefficients to $\le1.5\cdot10^{-6}$ in $576\cdot$coefficient.
9. c2: the b = 0 argument (a Laurent polynomial in $\rho$ integrable on $(0,1)^3$ has no negative powers) is correct; the general $b>0$ span is conditional on the form of the reduced integrand.

## 6. Not done

- No independent exact (symbolic) derivation of any $\zeta(3)$-valued coefficient; the precision of §3 (6–7 digits) does not allow an independent PSLQ identification with denominators $\sim10^5$–$10^6$.
- a2 case statements and the exported $R_k$ not re-derived.
- b6 relations not implemented.
- K33 sources $\{0,1,2,3\}$: no exact value.
