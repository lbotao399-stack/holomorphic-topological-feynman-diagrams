# VERIFY: fans-ht

Scripts: `explore/fans-ht/verify/*.py`; printed outputs: `explore/fans-ht/verify/logs/*.log`.
`verify/vlib.py` is independent of `code/fanlib.py` and of `code/htk` (neither is imported). `explore/fans-ht/REPORT.md` does not exist; the logical-gap review below is made on the structured report.

## 0. Brute force used here

Fan of S1 Theorem thm-fan-kernel: path vertices $1,\dots,n+1$, root $O=n+2$, $a_i=(i,i+1)$, $b_i=(i,O)$, edge order $a_1,b_1,\dots,a_n,b_n,b_{n+1}$.

Cycle basis: fundamental cycles of the path tree $T=\{a_1,\dots,a_n,b_{n+1}\}$,
$$c_r=b_r-b_{n+1}-\sum_{k=r}^na_k,\qquad r=1..n,$$
first row flipped if $\det[B^{\mathsf T};C]<0$ (the triangle basis is used only for claim a2).

Currents (mesh analysis, not nodal): $y=y^0+C_T^{\mathsf T}\ell$, $B_T^{\mathsf T}y^0_T=\lambda$, $(C_TtC_T^{\mathsf T})\ell=-C_Tty^0$; current and voltage laws asserted exactly.

$$\mathrm{d}M_{ra}=\sum_eC_{re}C_{ae}\,\mathrm{d}t_e,\quad A_{rs}=-\tfrac12\sum_{a,b}\mathrm{d}M_{ra}(M^{-1})_{ab}\wedge\mathrm{d}M_{bs},\quad \mathrm{Pf}(A)=\sum_{\pi}\mathrm{sgn}(\pi)\bigwedge_{(i,j)\in\pi}A_{ij},$$
$$\Omega=U^{-H}\bigl(U^{-1/2}\mathrm{Pf}A\bigr)^T\wedge\prod_{i,r}\psi^i_r,\qquad D=\frac{(-1)^{e_0}}{t_{e_0}}[\Omega]_{\widehat{\mathrm dt_{e_0}}}\ \ (e_0=N-1\text{ and }e_0=0\text{ agree, asserted}),$$
$$g_{\rm brute}(r)=D(t(r))\,\bigl|\det[t;\partial t/\partial r]\bigr|,$$
$t(r)$ = own inversion of S1 eq-fan-local-ratios with the final triangle of eq-1-25, $\partial t/\partial r$ by sympy differentiation, all arithmetic in `Fraction`. Contraction of S1: $\lambda_{+\dot a}\mu_+{}^{\dot a}=\det(\lambda,\mu)$ from $\epsilon^{\dot+\dot-}=1$ (S1 eq-sum-spinor-conventions).

## 1. Verdicts

| id | verdict |
|---|---|
| a0-chart-U | confirmed |
| a1-green | confirmed |
| a2-rho | confirmed |
| a3-gaussian-sources | confirmed |
| a4-potential-form-S1 | confirmed |
| a5-density-general | confirmed |
| a6-monomial-times-polynomial | confirmed |
| a7-endpoint-closed-form | confirmed |
| a8-general-closed-form | confirmed |
| a9-K2-K4 | confirmed |
| a10-K6-K8 (derived-unchecked; not required) | n=6: confirmed at 3 source vectors; n=8: unverified |
| a11-structure | confirmed |
| a12-reflection | confirmed |
| b1-taylor | confirmed |
| c1-n2-numeric | confirmed |
| c2-n4-numeric | confirmed (in the S1 chart only) |
| d1-control-20 | confirmed |

## 2. Checks

### 2.0 a0 (`v1_chart_green.py`, `jac_sym.py`)

$$t(r)\ (n=2)=\Bigl[\tfrac{u_2(1-u_1)}{1-v_1},\ \tfrac{u_1u_2}{1-v_1},\ 1-u_2-v_2,\ \tfrac{u_2}{v_1},\ v_2\Bigr]$$
(every $t_e(r)$ for $n=2,4,6,8$ is a product of powers of $u_k,1-u_k,v_k,1-v_k,1-u_n-v_n$, positive on the domain). Symbolic, $n=2,4$: $r(t(r))=r$ and $t(r(t))=c(t)\,t$.

$$\det[t;\partial_{u_1}t;\partial_{v_1}t;\dots;\partial_{v_n}t]=(-1)^{n-1}\prod_{k=1}^{n-1}\frac{u_{k+1}^{2k}}{v_k^2(1-v_k)^{2k+1}}=(-1)^{n-1}\,U^2\prod_{k=1}^{n-1}(1-v_k)^{-1}.$$
Printed: `n=2: ... = u2**2/(v1**2*(v1 - 1)**3)`; `n=4: ... = u2**2*u3**4*u4**6/(v1**2*v2**2*v3**2*(v1 - 1)**3*(v2 - 1)**5*(v3 - 1)**7)` (symbolic); $n=6,8$: equal to the product formula exactly at 3 points each. Hence sign $-1$ on the whole domain for $n=2,4$; at the tested points for $n=6,8$.

$U=\det M$ (path-tree basis) $=\prod_k u_{k+1}^k(1-v_k)^{-k}v_k^{-1}=\prod_{i<n}R_i/v_i$, exact, 3 points each, $n=2,4,6,8$.

### 2.1 a1 (`v1_chart_green.py`)

$G=\Lambda^{-1}$ from own Laplacian (exact inverse): $G_{uv}=g_uh_v$ for all $u\le v$, $g_ih_{i+1}-g_{i+1}h_i=-a_i$, $G_{nn}=u_n(1-u_n)$, $G_{n,n+1}=u_nv_n$, $G_{n+1,n+1}=v_n(1-v_n)$; exact at 3 points each, $n=2,4,6,8$.

### 2.2 a2 (`v3_rho.py`)

Plain triangle basis $c_r=a_r+b_{r+1}-b_r$: $\det[B^{\mathsf T};C_{\triangle}]=-8,\ 55,\ -377$ ($n=2,4,6$). Printed for $n=2,4,6$ (2 points each): `Pf - RHS = 0: True;  Pf + RHS = 0: False`, identity of $n$-forms on $\mathbb R^N$ in $\mathrm dt_e$ ($8$, $55$, $377$ terms), and the ratio identities $a_r/b_r$, $a_r/b_{r+1}$ in $r$ exact.

Derivation (own): S3 eq-currents, eq-psi in the triangle basis,
$$\psi^\tau_r=w_{a_r}\mathrm da_r+w_{b_{r+1}}\mathrm db_{r+1}-w_{b_r}\mathrm db_r=\tfrac12(\tau_r\alpha_r-\tau_{r+1}\beta_r);$$
S3 eq-ww, $\langle\tau_u\tau_v\rangle=2G_{uv}=2g_{\min}h_{\max}$; for $i<j$ (so $i+1\le j$)
$$\langle\psi^\tau_i\psi^\tau_j\rangle=\tfrac12\bigl(g_i\alpha_i-g_{i+1}\beta_i\bigr)\wedge\bigl(h_j\alpha_j-h_{j+1}\beta_j\bigr);$$
S3 prop-rho: $\mathrm{Pf}=\sum_\pi\mathrm{sgn}(\pi)\bigwedge_{(i,j)\in\pi}\langle\psi_i\psi_j\rangle=2^{-n/2}\sum_\pi\mathrm{sgn}(\pi)^2\bigwedge_rC^\pi_r$.

### 2.3 a3 (`v3_rho.py`)

$\mathrm{Pf}(-\tfrac12\mathrm dM\,M^{-1}\mathrm dM)=E_\kappa\bigl[\prod_r\psi_r[\kappa]\bigr]$, $\kappa\sim N(0,\Lambda/2)$, by exact expansion of $\prod_r\sum_v\kappa_v\omega_{r,v}$ with $\kappa$-monomials kept and Isserlis moments; and $y[\Lambda\tau/2]=t^{-1}B\tau/2$ (mesh solver). Printed `True` for: fan2 (path-tree and triangle bases), fan4, doubled triangle ($n=2$), wheel $W_4$ ($n=4$), $K_4$ with a doubled edge ($n=4$); 2 random rational $t$ each. Then
$$U^{-1}\rho_n\wedge\prod_r\psi_r[\lambda]=U^{-3/2}E_\kappa\Bigl[\prod_r\psi_r[\kappa]\Bigr]\wedge\prod_r\psi_r[\lambda]=U^{1/2}E_\kappa\bigl[\Omega^{(2,0)}[\kappa,\lambda]\bigr].$$

### 2.4 a4 (`v2_density.py`, `v7_d1_all_n.py`)

$$g^{(2,0)}_{\rm brute}=(-1)^{\frac{n(n+1)}2}\,\mathcal I^{\rm S1}_{\text{eq-1-40}}(r)=(-1)^{\frac{n(n+1)}2}J(r)\prod_{i=1}^n\det(\Phi_i,\Phi_{i+1}),$$
$\mathcal I^{\rm S1}_{\text{eq-1-40}}$ computed from S1 eq-1-40 with currents eq-1-38 and $\lambda^{\rm red}$ of eq-B-fan-first-cycle, random integer 2-vector sources at all path vertices: exact, $n=1..7$ (2 points), $n=2,4,6$ (3,3,2 points), $n=8$ (1 point). Intermediates exact ($n=1..7$): eq-1-38 currents $=$ mesh currents; $\lambda^{\rm red}_i=(\Phi_i-\Phi_{i+1})/a_i+\Phi_i/B_i$; $\lambda^{\rm red}_{i+1}-y_{i+1,i+2}=\Phi_{i+1}/B_{i+1}$.

### 2.5 a5 (`v2_density.py`, `v8_chain_a5.py`)

$g^{(1,1)}_{\rm brute}=2^{-n/2}\Pi\,\mathrm{Haf}_n[L,M]$ (sign and square), random integer sources at all path vertices: $n=2$ (3 points), $4$ (3), $6$ (2), $8$ (1). At the same points: $\Pi^2=U\,(J\prod a_i)^2$ with equal signs, $M_j=\sum_{w>j}h_w\lambda_w$, $L_jP_j=\lambda^{\rm red}_j$. Chain steps, $n=2,4,6,8$ (3 points): $\hat\alpha_i=-a_iL_i$, $\hat\beta_j=a_jM_j$, $\langle l_il_j\rangle=2\hat\alpha_i\hat\beta_j$ ($i<j$). Sign: $(-1)^{\frac{n(n+1)}2+\frac n2}=(-1)^{\frac n2(n+2)}=+1$.

### 2.6 a6 (`v4_integrate.py`)

$p_n=g_n/m_0=2^{-n/2}(\Pi/m_0)\,\mathrm{Haf}_n$,
$$\frac{\Pi}{m_0}=v_n^{-1}\prod_{k<n,\ k\ \rm odd}(1-v_k)^{\frac{k-1}2}\prod_{k<n,\ k\ \rm even}(1-v_k)^{\frac k2-1}\prod_{k\ \rm odd}u_k^{-\frac{k+1}2}\prod_{k\ \rm even}u_k^{-\frac k2}.$$ `sympy.cancel`: denominator free of $r$ for $n=2,4$ (symbolic $\lambda$; 1 and 27 monomials in $r$) and $n=6$ (3 integer source vectors). Exponents re-derived from $L_i=\sum_{w\le i}\lambda_w/P_w$, $P_k^2c_k/P_w$, $c_k=u_k(1-u_k)R_k$:
$$e(u_l)=-\tfrac{l+1}2-\#\{w_i>l\}+\sum_{\chi_j=(w,k)}\bigl(1+[w\le l<k]\bigr)+p=\tfrac{n-1-l}2-\#\{w_i>l\}+\#\{w\le l<k\},$$
$e(u_n)=-\tfrac12$, $e(1-u_l)=\#\{k=l\}$, $e(v_k)=\tfrac12$, $e(v_n)=p-1$, $e(1-v_m)=\tfrac m2-1-c_m$; equal to the statement.

### 2.7 a7, d1 integrals (`v5_closed_structure_taylor.py`)

Own Beta products: $\int g_{\rm end}=(4\pi)^{n/2}/((n+1)(n!)^2)=\pi^{n/2+1}/(2^{n+1}\Gamma(\tfrac n2+1)^2\Gamma(\tfrac{n+1}2)\Gamma(\tfrac{n+3}2))$, $n=2..30$ even; printed `pi/3, pi**2/180, pi**3/56700, pi**4/57153600, pi**5/141455160000` ($n=2..10$). $\int$(eq-1-41 density) $=1/((n!)^3(n+1))$, $n=1..30$. Pointwise endpoint (1,1) density: $n=2,4,6,8$ (`v2_density.py`).

### 2.8 a8 (`v4_integrate.py`, `v5_closed_structure_taylor.py`)

Own implementation of the configuration sum with $N(O)=\prod_{j\notin O}\bigl(\#\{i\in O:i<j\}-\#\{j'\notin O:j'<j\}\bigr)$ equals termwise integration of $m_0p_n$: $n=2,4$ symbolic $\lambda$, $n=6$ at $\lambda=(3,-4,4,4,3,-4,-3)$, $(-3,-1,3,1,2,3,-3)$, $(-1,3,0,3,-2,-2,-3)$: $\tfrac{237}{350}\pi^3$, $\tfrac{3086}{525}\pi^3$, $\tfrac{197}{1050}\pi^3$. Ensemble: $m_jP_nu_nv_n=M_j$ and $g_n=w(r)\,\mathrm{Haf}_n(l,m)/(n-1)!!$, $w=g_{\rm end}|_{\lambda_1=\lambda_{n+1}=1}$, exact at 3 points, $n=2,4,6,8$; $w\propto\prod_{m<n}u_m^{\frac{n-1-m}2}\prod_{k<n}v_k^{\frac12}(1-v_k)^{\frac k2-1}\cdot u_n^{-\frac12}v_n^{\frac n2-1}$, i.e. Beta$(\tfrac{n+1-m}2,1)$, Beta$(\tfrac32,\tfrac k2)$, Dir$(\tfrac12,\tfrac n2,1)$, $\int w$ from 2.7.

### 2.9 a9 (`v4_integrate.py`, `v6_numeric.py`)

Termwise exact: `K_2 = pi*l1*l3/3`; `K_4 = pi**2*l1*l5*(2*l1*l3 + 2*l1*l4 + l1*l5 + 4*l2*l3 + 4*l2*l4 + 2*l2*l5 + 12*l3**2 + 4*l3*l4 + 2*l3*l5)/180`; equal to the claim and to `general_kernels.txt`. Numerics: 2.13, 2.14.

### 2.10 a10 (not required)

$n=6$: 2.8 (three source vectors, `general_kernels.txt` equal). $n=8$: `v4b_n8.py` (termwise integration at fixed integer sources) killed after 37 min wall time without result; $K_8$ unverified here.

### 2.11 a11 (`v5_closed_structure_taylor.py`)

(i) $1$ is the smallest element of every pair containing it, $n$ the largest: $L_1=\lambda_1$, $M_n=P_nu_nv_n\lambda_{n+1}$ in every term, so $\lambda_1\lambda_{n+1}\mid g_n$ for every $n$. (ii) For $i\ge2$, $L_i=u_1^{-1}\bigl[(\lambda_2+u_1\lambda_1)+\sum_{w\ge3}\lambda_w\,u_1/P_w\bigr]$; for $j\ge2$, $M_j=u_1\times(\text{free of }u_1,\lambda_1,\lambda_2)$; $M_1$ never occurs; $\Pi\propto u_1^{-1}$; powers of $u_1$: $-1-(\tfrac n2-1)+\tfrac n2=0$. Symbolic check $n=2,4$: `(g/m_0)/lam_1 with lam_2 = s - u_1 lam_1 is free of u_1 and lam_1: True`. (iii) $K=\Psi(\lambda_1+\lambda_2)-\Psi(\lambda_2)$, $\Psi(x)-\Psi(0)=C_nx^{n/2}\lambda_{n+1}^{n/2}$ by a7. Printed on the polynomials, $n=2,4,6,8$: `lam_1 lam_{n+1} divides K: True; K(lam_1, lam_2, 0.., lam_{n+1}) = C_n lam_{n+1}^{n/2}[(lam_1+lam_2)^{n/2} - lam_2^{n/2}]: True`.

### 2.12 a12 (`v5_closed_structure_taylor.py`)

Reflection $i\mapsto n+2-i$, $t'_{a_j}=t_{a_{n+1-j}}$, $t'_{b_j}=t_{b_{n+2-j}}$, $\lambda'_j=\lambda_{n+2-j}$. Brute force: $D(t';\lambda')/D(t;\lambda)=1$ exactly, $n=2,4,6$ (2 points), $n=8$ (1 point). Positions $\pi$: evens reversed ($(-1)^{\lfloor n/2\rfloor}$), $\{3,5,\dots,2n-1\}$ reversed ($(-1)^{\lfloor (n-1)/2\rfloor}$), transposition $(1,2n)$:
$$\mathrm{sgn}\,\pi_E=(-1)^{\lfloor n/2\rfloor+\lfloor(n-1)/2\rfloor+1}=(-1)^n ;$$
computed $+1$ for even $n=2..40$. For $D=2$ the cycle-basis normalization, edge reversals and vertex relabelling enter $\Omega$ through even powers; only $\iota_E\prod\mathrm dt_e\to\mathrm{sgn}(\pi_E)\,\iota_E\prod\mathrm dt_e$ remains, so $K(\lambda')=(-1)^nK(\lambda)=K(\lambda)$ for even $n$. Polynomials $n=2,4,6,8$ (`general_kernels.txt`): difference $0$.

### 2.13 b1 (`v5_closed_structure_taylor.py`, `v6_numeric.py`)

Mesh currents with endpoint sources: $y_{i,i+1}=\sum_{k\ge i}\mu_{ik}\lambda_1-\nu_i\lambda_{n+1}$ and $\sum_ey_ez_e=\lambda_1z_{1,O}+\lambda_{n+1}z_{n+1,O}+\sum_iy_{i,i+1}Z_i$, exact at 3 points, $n=2,4,6$. $c^{(n)}$ (closed form) $=\int g_{\rm end}(\lambda_1\lambda_{n+1})^{-n/2}\prod\mu_{ik}^{r_{ik}}\prod\nu_i^{s_i}$ (own Beta integration of the collected monomial exponents) for 60 multi-indices each (entries $\le2$, including $0$), $n=2,4,6$; $c_0=\pi/3,\ \pi^2/180,\ \pi^3/56700$. Order 1, $-\sum_ez_e\int g_{\rm end}y_e$ exact, equal to the stated $n=2$ and $n=4$ expressions (`equals claim: True`). Chart-free lattice, $n=2$, $\lambda=(1,0,1)$, $h=0.1$, $X=3.6$: coefficients of $z_{a_1},z_{b_1},z_{a_2},z_{b_2},z_{b_3}$: $-0.3141596920,\ -0.7330389024,\ 0.3141593963,\ -0.6283190883,\ -0.7330391981$ vs $-\tfrac\pi{30}(3,7,-3,6,7)$, rel.diff $1.36\cdot10^{-6},8.41\cdot10^{-7},4.17\cdot10^{-7},8.87\cdot10^{-7},1.24\cdot10^{-6}$ (5–6 digits).

### 2.14 c1 (`v6_numeric.py`)

Float brute force (own vectorized code, path-tree basis, nodal solve, $A_{rs}[e,f]=-\tfrac12Q_{ef}(C_{re}C_{sf}-C_{rf}C_{se})$, $Q=C^{\mathsf T}M^{-1}C$); self-test against the exact brute force: rel.diff $\le1.8\cdot10^{-15}$ ($n=2,4$).

Chart-free: $K=\int_{\mathbb R^{4}}D(e^s)\prod_ee^{s_e}\,\mathrm d^4s$ ($s_{b_3}=0$), $s=2\sinh x$, trapezoid in $x\in[-X,X]^4$, points with $\max t/\min t>e^{32}$ skipped (counted):

| $h$ | $X$ | $K(1,0,1)$ | rel.diff | $K(3,-2,5)$ | rel.diff |
|---|---|---|---|---|---|
| 0.20 | 3.2 | 1.047489938373 | 2.79e-04 | 15.712349075599 | 2.79e-04 |
| 0.15 | 3.2 | 1.047229838658 | 3.08e-05 | 15.708447579876 | 3.08e-05 |
| 0.10 | 3.2 | 1.047192793907 | 4.54e-06 | 15.707891908598 | 4.54e-06 |
| 0.10 | 3.6 | 1.047198594342 | 9.96e-07 | 15.707978915134 | 9.96e-07 |

(6 digits; exact $\pi/3=1.047197551197$, $5\pi=15.707963267949$.) In the S1 chart, Gauss–Legendre in $\theta$ with $r=\sin^2\theta$, $v_n=(1-u_n)\sin^2\phi$, $m=8$: $1.0471975511970$, $15.7079632679550$, rel.diff $3.86\cdot10^{-13}$.

### 2.15 c2 (`v6_numeric.py`)

S1 chart, same substitution, $|\det|$ from 2.0 (symbolic at $n=4$), brute-force float density:

| $m$ (nodes) | $K(1,0,0,0,1)$ | rel.diff | $K(2,-1,3,1,-2)$ | rel.diff |
|---|---|---|---|---|
| 4 (65536) | 0.054804055873 | 4.94e-04 | -23.684044525825 | 1.27e-04 |
| 5 (390625) | 0.054831118009 | 3.20e-07 | -23.687061932382 | 4.80e-07 |
| 6 (1679616) | 0.054831157747 | 4.05e-07 | -23.687052043720 | 6.25e-08 |
| 7 (5764801) | 0.054831134790 | 1.41e-08 | -23.687050507889 | 2.31e-09 |

exact $\pi^2/180=0.054831135562$, $-\tfrac{12}{5}\pi^2=-23.687050562614$ (7–8 digits at $m=7$). No chart-free number at $n=4$.

### 2.16 d1 (`v5_closed_structure_taylor.py`, `v7_d1_all_n.py`)

$g^{(2,0)}_{\rm brute}|_{\rm endpoint}=(-1)^{\frac{n(n+1)}2}\det(\lambda_1,\lambda_{n+1})^nv_n^{n-1}\prod_{j<n}u_j^{n-j}v_j(1-v_j)^{j-1}$ exact, $n=1..7$ (2 points), $n=8$ (1 point); $\det(\Phi_i,\Phi_{i+1})=\det(\lambda_1,\lambda_{n+1})P_nu_nv_na_i$ and $J\prod a_i(P_nu_nv_n)^n=v_n^{n-1}\prod_{j<n}u_j^{n-j}v_j(1-v_j)^{j-1}$ exact, $n=1..7$; integral $1/((n!)^3(n+1))$, $n=1..30$.

## 3. Logical gaps in the explorer's report

1. REPORT.md is absent; a5 cites "(3.1)", a10 "(3.4)", c2 "(5.2)", a8 "(5.1)", equations of a document that does not exist.
2. a0: the sign of $\det[t;\partial t/\partial r]$ is asserted for all $n$ from rational points ($n=2,4,6,8$); no derivation. Here: closed form $(-1)^{n-1}U^2\prod(1-v_k)^{-1}$, proved for $n=2,4$, pointwise for $n=6,8$; general $n$ open. Bijectivity evidence in a0 is one-sided ($r(t(r))=r$ at points); $t(r(t))\propto t$ (surjectivity onto $\mathbb{RP}^{2n}_{>0}$) is not in the explorer's evidence (checked here symbolically for $n=2,4$; general $n$ from the stepwise inversion of eq-fan-local-ratios).
3. a2 holds for the plain triangle basis. $\det[B^{\mathsf T};C_\triangle]=-8,55,-377$ ($n=2,4,6$); in the basis normalized by $\det[B^{\mathsf T};C]>0$ (htk CONVENTIONS Section 1) $\rho_n$ is minus the stated expression for $n=2,6$. The statement does not say this; $\Omega$ at $D=2$ is unaffected.
4. a3 ("any graph") had no direct evidence in the report; checked here on three non-fan graphs.
5. a4: the explorer's evidence (P4) compares the potential form with the brute force only; equality with eq-1-40 as printed in S1 was not shown (shown here, 2.4).
6. a5: the step from a3 to the hafnian needs the $(2,0)$ potential form a4 for the pair of sources $(\kappa,\lambda)$ with $\Phi=(\tau/2,\varphi)$ and Isserlis for the commuting $l_i$ with only $i<j$ covariances; these steps are stated without proof in the report (checked here, 2.5). General $n$ rests on that derivation; the pointwise checks cover $n\le8$.
7. a6: the explorer's evidence (minimal monomials of the expansion) does not test polynomiality of $p_n$; the bounds argument in the statement is correct (re-derived, 2.6).
8. a8: $N(O)$ is not defined in the statement (definition used here in 2.8).
9. a11 (iii) and a12 for $n=6,8$ use the polynomials of `general_kernels.txt`; $K_6$ is checked here at three source vectors only, $K_8$ not at all. a11 (i)–(iii) hold for all even $n$ by 2.11.
10. a12: general-$n$ sign not derived in the report; derived here, $\mathrm{sgn}\,\pi_E=(-1)^n$ (2.12).
11. b1: the closed form of $c^{(n)}$ is stated without the exponent bookkeeping; here compared with direct Beta integration for random multi-indices ($n=2,4,6$) only.
12. c2: no chart-free determination at $n=4$ beyond 1 digit (explorer); none here.

## 4. Not done / open

- Chart-free numerical integration at $n=4$ (8-dimensional) not done.
- General-$n$ proof of the Jacobian formula of 2.0 not done.
- a10 at $n=8$: `v4b_n8.py` killed after 37 min (machine load average about 30) without result; $K_8$ of `general_kernels.txt` unverified here. a11 (iii) and a12 at $n=8$ on the polynomial level rest on that unverified $K_8$; a12 at $n=8$ is checked pointwise (1 point) and a11 holds for all even $n$ by the derivation in 2.11.
- $n=8$ pointwise checks of a5 and a12 use one rational point each.
