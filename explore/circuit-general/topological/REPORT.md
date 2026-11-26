# Topological directions in the circuit picture: $T\ge2$, $H=0$

Tags: [chk] exact check run (script, log); [chk-num] floating-point check; [drv] derived, not independently checked; [src] source; [open].
Sources:
- GKW = `研究项目_非Git/HT_Papers_LaTeX/2_Higher_Operations/Higher.tex`. Equation numbers are by section and were counted by script; (3.15), (3.18), (3.52) agree with PLAN.md E8.
- CONV = `code/htk/CONVENTIONS.md`.
- S3 = `Twisted_11d_sugra_loops_2026-11-18/main.tex`.
- LLC = `sources/LLC_main.tex`.
- BG = arXiv HTML 2408.03192.
- WW = arXiv HTML 2407.08667.
- BK = arXiv HTML 2209.14438v2.

Code (own; htk not imported):
- `tb.py`: circuit, forms, Wick, Pfaffian, hafnian, spanning-tree fiber integral.
- `graphs.py`: admissible graphs.
- `wick.py`: exact derivatives, form-valued Gaussian polynomials.
- `moves.py`, `psm.py`, `psm2.py`.
- Checks `t1_*.py`, `t2_*.py`, `t3_banana.py`, `t4_psm.py`, logs `*.log`.

## 0. One topological direction

$$w_e=\frac{\tau_e}{2t_e},\quad \eta_e=\frac{\tau_e}{2\sqrt{t_e}}=\sqrt{t_e}\,w_e,\quad \Theta_e=\pi^{-1/2}e^{-\eta_e^2}d\eta_e,\quad d\eta_e=\frac{d\tau_e-w_e\,dt_e}{2\sqrt{t_e}},\quad \sum_e\eta_e^2=\tfrac14\tau^{\mathsf T}\Lambda\tau$$
[src S3 eq-currents, eq-ww; GKW (3.15) with $x^{\mathbb R}=\tau/2$, $s_e=\eta_e$]

$$\mathcal T_\Gamma[F;\kappa]:=U^{-1/2}e^{P}\Big\langle F\prod_{r=1}^n\psi_r\Big\rangle,\quad \psi_r=\sum_eC_{re}w_e\,dt_e,\quad P=\sum_et_e\bar w_e^2=\kappa^{\mathsf T}\Lambda^{-1}\kappa,$$
$$\langle w_e\rangle=\bar w_e\ (\text{Kirchhoff, sources }\kappa),\qquad \langle\delta w_e\,\delta w_f\rangle=\tfrac12\big(\delta_{ef}/t_e-K_{ef}\big),\qquad K=C^{\mathsf T}M^{-1}C .$$
[src CONV §2, S3 eq-ww]. $\mathcal T_\Gamma[1;0]=\rho_n(M)$. With $T$ directions the factor is $\prod_a\mathcal T_\Gamma[F_a;\kappa^a]$: independent Gaussians, common ratio coordinates.

## 1. $\alpha_\Gamma$ and $\rho_n$

$$\int_{\mathbb R^{V-1}}\!e^{\kappa\cdot\tau}F\prod_e\Theta_e=(2\sqrt\pi)^{-n}U^{-1/2}e^{P}\,\Phi_{\rm GKW}[F],\qquad \Phi_{\rm GKW}[F]=(-1)^n\sum_{T}\epsilon_T\det B_T\Big\langle F\!\!\prod_{e\notin T}\!w_e\Big\rangle dt_{T^c}$$
($\prod_vd\tau_v$ in front; $\epsilon_T$ = sign of moving the tree factors to the front) [drv; chk t1_normalization.log]

$$\Phi_{\rm GKW}=(-1)^n(\det G)^{-1}\Big\langle\prod_r\psi_r\Big\rangle,\qquad \det G=1\ \text{for the positive }C\text{ of CONV §1}$$
[drv from CONV §4 with $c=1$; chk on 88 graphs ($D\le4$, $n\le3$; $D\le2$, $n=4$; banana $n=6$) and on non-fundamental bases]

$$\boxed{\alpha_\Gamma:=\int_{\mathbb R^{V-1}}\prod_e\pi^{-1/2}e^{-s_e^2}ds_e=(-1)^n(4\pi)^{-n/2}\rho_n(M_C)},\qquad n\ \text{even}:\ \alpha_\Gamma=(4\pi)^{-n/2}\rho_n,\qquad \omega_\Gamma=\alpha_\Gamma^{\,T}\prod_i\rho_\Gamma[\lambda^i;z^i]$$
[src GKW (3.52)–(3.54), l.1445–1458; chk]
- GKW bitriangle (3.56), l.1500–1508: edges $(t_{12},t_{23},t_{24},t_{31},t_{43})$, vertex order $v_1,v_2,v_3$, root $v_4$. $\alpha=\frac1{8\pi}U^{-3/2}{\rm Br}\iff 2U\Phi_{\rm GKW}={\rm Br}$, sign included [chk at 4 points].
- Banana: $I_N=(4\pi)^{-n/2}K_N$ (§5) [chk].

$$\rho_n=U^{-1/2}{\rm Pf}\big(-\tfrac12dMM^{-1}dM\big)=(-2)^{-n/2}U^{-1/2}\!\!\sum_{|S|=n}\!\det C_S\,{\rm Haf}(K_{SS})\,dt_S=(-1)^{n/2}\pi^{-n/2}\!\int\! d^n\ell\,e^{-\ell^{\mathsf T}M\ell}\prod_r(dM\,\ell)_r$$
$K_{ef}$ = current in $f$ per unit EMF in $e$. [src S3 prop-rho; Pf = Wick = Haf: chk on 88 graphs; last form drv]

## 2. Moves on the topological factor

**2.0 Decomposition.** Subnetwork $S$, internal vertices $I$: $\tau_I=H\tau_\partial+\delta$ (harmonic extension plus internal noise).
$$\tau^{\mathsf T}\Lambda\tau=\tau_{\rm red}^{\mathsf T}\Lambda'\tau_{\rm red}+\delta^{\mathsf T}\Lambda_I\delta,\quad \kappa^{\rm red}=\kappa_{\rm red}+H^{\mathsf T}\kappa_I,\quad P=P'+\kappa_I^{\mathsf T}\Lambda_I^{-1}\kappa_I,\quad U=U_S\,U'$$
- $\delta$ is independent of the $\Gamma'$ fluctuations.
- $\psi(\text{lift}_r)=\psi'_r+\sum_EC'_{rE}\mu_E$, i.e. the odd shift $\zeta_E\to\zeta_E+\mu_E$ with $\zeta_E=w'_Edt'_E$, plus local cycles $\psi_c$.
$$\textstyle\prod\psi'(\zeta+\mu)=\sum_B\mu_B\wedge\partial_B\prod\psi',\qquad \prod\psi=\zeta_B\wedge\partial_B\prod\psi+\dots$$
[drv; chk through 2.1–2.2]. Cut-edge form: $\mathcal T^{(B)}_{\Gamma'}[F;\kappa]:=U'^{-1/2}e^{P'}\langle F\,\partial_B\prod\psi'\rangle$.

**2.1 Chain move.** Notation of CONV §5: chains $i=1..k$ from $p$ to $q$, $R_i$, $x_{ij}$, $t^{\rm red}$, $r_i=t^{\rm red}/R_i$, $\sigma_{ik}$, $\hat m_i$, $Q_{ij}$. Series is $k=1$; parallel is $a_i=1$. Basis: lifts of the $\Gamma'$ basis through chain 1, then chain$_i-$chain$_1$.

$$\kappa^{\rm red}_p=\kappa_p+\sum_i\hat m_i,\quad \kappa^{\rm red}_q=\kappa_q+\sum_{i,k}(1-\sigma_{ik})\kappa_{ik},\quad \bar w_{ij}=r_i\bar w_*-\hat m_i+Q_{ij},$$
$$P=P'+\sum_iR_i\,{\rm Var}_{x_i}(Q_i),\quad {\rm Var}_{x_i}(Q_i)=\sum_jx_{ij}Q_{ij}^2-\hat m_i^2,\qquad U=\frac{\prod_iR_i}{t^{\rm red}}\,U'$$
[drv; chk t2_moves.log, 1078/1078 cases: every chain bundle of every admissible graph with $D\le4$, $n\le3$ (root not interior), single chains, bananas, multigraphs; random $t$ with random $\kappa$ and with $\kappa=0$]

$$w_{ij}=r_iw_*+\nu_{ij},\quad \langle\nu_{ij}\rangle=Q_{ij}-\hat m_i,\quad \langle\delta\nu_{ij}\delta\nu_{ij'}\rangle=\tfrac12\big(\delta_{jj'}/t_{ij}-1/R_i\big),\quad \nu_i\perp\nu_l,\ \nu_i\perp\Gamma',$$
$$\psi(\text{lift}_r)=\psi'_r+C'_{r*}\chi_1,\quad \psi(\text{chain}_i-\text{chain}_1)=\chi_i-\chi_1,\quad \chi_i=-t^{\rm red}w_*\,d\log r_i+R_i\sum_j\nu_{ij}\,dx_{ij}$$
[drv; chk t2_noise.log, 300/300 cases with insertions $F\in\{1,w_e,w_ew_f,w_e^2w_f\}$]

$$\Big\langle F\prod_\Gamma\psi\Big\rangle=\Big\langle F\prod\psi'\wedge\partial(\chi_1\cdots\chi_k)\Big\rangle+(-1)^{n'-1}\Big\langle F\,\partial_\zeta\Big(\prod\psi'\Big)\wedge\chi_1\cdots\chi_k\Big\rangle,\quad \partial(c_1\cdots c_k)=\sum_j(-1)^{j-1}c_1\cdots\widehat{c_j}\cdots c_k$$
For $F=1$: $\chi_i\to\bar\chi_i=R_i\,(d\hat m_i-w_*\,dr_i)$. [drv; chk t2_moves.log, and the three-chain moves of theta $(p,q,r)$ in t2_noise.log]

$\kappa_I=0$ (in particular $\kappa=0$), using $\partial(R_1dr_1\cdots R_kdr_k)=(t^{\rm red})^{k-1}\prod_ir_i^{-1}dr_2\cdots dr_k$ on the simplex:
$$\boxed{\mathcal T_\Gamma[\kappa]=\mathcal T_{\Gamma'}\big[(-\eta_*)^{k-1};\kappa^{\rm red}\big]\wedge\prod_{i=1}^kr_i^{-1/2}\,dr_2\wedge\cdots\wedge dr_k,\qquad \eta_*=\sqrt{t^{\rm red}}\,w_*}$$
- $k=1$: $\mathcal T_\Gamma=\mathcal T_{\Gamma'}|_{t_*=R}$ (pullback; no $dx_{ij}$).
- Theta, $\Gamma'$ = one edge: $\rho_2=\tfrac12\prod_Pr_P^{-1/2}dr_2dr_3$. [drv; chk]

Series move in LLC labels ($v$ between $u,w$, $r=t_{vw}/(t_{vu}+t_{vw})$, $R=t_{vu}+t_{vw}$):
$$\mathcal T_\Gamma=e^{r(1-r)R\kappa_v^2}\Big(\mathcal T_{\Gamma'}+R\kappa_v\,dr\wedge\mathcal T^{(e^*)}_{\Gamma'}\Big),\quad \kappa^{\rm red}_u=\kappa_u+r\kappa_v,\ \kappa^{\rm red}_w=\kappa_w+(1-r)\kappa_v,\quad \bar w_{vu}=r\kappa_v-\bar w^{\rm red},\ \bar w_{vw}=(1-r)\kappa_v+\bar w^{\rm red}$$
[drv; chk as the case $k=1$]

**2.2 Star–triangle with source** (LLC (eq-sum-star-triangle), App. E). Edges $a=vu$, $b=vw$, $c=vx$, $d=ux$, $f=wx$; no edge $uw$.
- $\Sigma=t_at_b+t_bt_c+t_ct_a$; $T_U=\Sigma/t_c$, $T_X=\Sigma/t_b$, $T_W=\Sigma/t_a$; $S=T_U+T_X+T_W$.
- $t'_X=(1/t_d+1/T_X)^{-1}$, $t'_W=(1/t_f+1/T_W)^{-1}$; $r_1=t'_X/t_d$, $r_2=t'_W/t_f$; $\Lambda_{vv}=1/t_a+1/t_b+1/t_c$.

$$\kappa^{\rm red}_u=\kappa_u+\tfrac{t_bt_c}{\Sigma}\kappa_v\ (\text{cyclic}),\quad P=P'+\frac{\kappa_v^2}{\Lambda_{vv}},\quad U=\frac{t_c}{r_1r_2}U',$$
$$w_a=-w'_U-(1-r_1)w'_X+\tfrac{\epsilon_v}{t_a},\ \ w_b=w'_U-(1-r_2)w'_W+\tfrac{\epsilon_v}{t_b},\ \ w_c=(1-r_1)w'_X+(1-r_2)w'_W+\tfrac{\epsilon_v}{t_c},\ \ w_d=r_1w'_X,\ w_f=r_2w'_W$$
- $\epsilon_v\perp\Gamma'$, $\langle\epsilon_v\rangle=\kappa_v/\Lambda_{vv}$, ${\rm var}\,\epsilon_v=1/(2\Lambda_{vv})$.
- Harmonic $Y$–$\Delta$ identity: $\sum_{e\in Y}J_ew^h_e\,dt_e=\sum_{E\in\Delta}\mathcal J_EW_E\,dT_E$ (polarized $dP|_{\rm injections}$).
- $\psi_\Delta=W_UdT_U+W_WdT_W-W_XdT_X$, with $W_X=(1-r_1)w'_X$, $W_W=(1-r_2)w'_W$.

$$\mu_U=-\tfrac{T_U}{S}\psi_\Delta+\epsilon_v\,d\log\tfrac{t_b}{t_a},\quad \mu_X=-\tfrac{t'_Xw'_X}{r_1}dr_1,\quad \mu_W=-\tfrac{t'_Ww'_W}{r_2}dr_2,$$
$$\psi_{c_1}=-\tfrac{T_X}{S}\psi_\Delta-\tfrac{t'_Xw'_X}{r_1(1-r_1)}dr_1+\epsilon_v\,d\log\tfrac{t_a}{t_c},\qquad \psi_{c_2}=\tfrac{T_W}{S}\psi_\Delta-\tfrac{t'_Ww'_W}{r_2(1-r_2)}dr_2+\epsilon_v\,d\log\tfrac{t_b}{t_c},$$
$$\Big\langle\prod_\Gamma\psi\Big\rangle=\Big\langle\prod_{r\ne\Delta}\big(\psi'_r+C'_{rU}\mu_U+C'_{rX}\mu_X+C'_{rW}\mu_W\big)\wedge\psi_{c_1}\wedge\psi_{c_2}\Big\rangle$$
Bases: the $\Gamma'$ basis ends with $\Delta=U+W-X$; lifts $U\mapsto-a+b$, $X\mapsto d$, $W\mapsto f$; $c_1=a+d-c$, $c_2=b+f-c$. [drv; chk t2_moves.log: 5 graphs with $n=2..5$, $\kappa\ne0$ and $\kappa=0$, 20/20 all identities true]

## 3. $\rho_n\wedge\rho_n=0$, $n\ge1$

**3.1 Series–parallel graphs.** At $\kappa=0$ the $k=1$ moves are pullbacks. At the first move with $k\ge2$, 2.1 gives $\rho_\Gamma=A\wedge\omega_k$, where $\omega_k=\prod r_i^{-1/2}dr_2\cdots dr_k$ is a top form on $\Delta^{k-1}$ and $A$ has no $dr$:
$$\rho_\Gamma\wedge\rho_\Gamma=\pm A\wedge A\wedge\omega_k\wedge\omega_k=0 .$$
[drv; local factor chk in 2.1]. Same mechanism: at $T\ge2$ any multiple edge violates the local count $Tk\le T+k-1$ [drv].

**3.2 All graphs, and every family $M(t)>0$.** With $X=M^{-1}dM$ (entries are 1-forms):
$$\rho_n\wedge\rho_n=U^{-1}\det\big(-\tfrac12dMM^{-1}dM\big)=(-\tfrac12)^n\det(X\wedge X),\qquad {\rm tr}(X\wedge X)^j={\rm tr}X^{2j}=-{\rm tr}X^{2j}=0\ \Rightarrow\ \det(X\wedge X)=0$$
The sign comes from a cyclic shift past $2j-1$ odd entries; the conclusion follows from the Newton identities over the commutative $\mathbb Q$-algebra of even forms. [drv; chk t1_rho_square.log]
- (P) ${\rm Pf}\wedge{\rm Pf}=0$ on $P_2$ and $P_4$; ${\rm tr}(X\wedge X)^j=0$ for $n=2,4,6$.
- (C) random non-graphic integer $C$, $n=2,4$.
- An earlier exact check of the same statement in the project: explore/ht-families/a_rho.out, A2.

**3.3 Circuit form of Balduf–Gaiotto** [src BG §4, Thm 2]. Cauchy–Binet gives $\det C_{S_1}\det C_{S_2}=U\det K_{S_1S_2}$ ($K$ = transfer admittance):
$$\rho_n\wedge\rho_n=2^{-n}\!\!\sum_{S_1\cap S_2=\emptyset}\!\!\det K_{S_1S_2}\,{\rm Haf}K_{S_1S_1}\,{\rm Haf}K_{S_2S_2}\,dt_{S_1}\wedge dt_{S_2},\quad Q_E=\!\!\sum_{E=S_1\sqcup S_2}\!\!{\rm sgn}(S_1,S_2)\det K_{S_1S_2}{\rm Haf}K_{S_1}{\rm Haf}K_{S_2}=0$$
Involution:
- Each term is a permutation word $(S_1,\beta S_1)$ together with matchings $\mu_1,\mu_2$.
- $\beta\cup\mu_1\cup\mu_2$ is a union of alternating cycles of length $4k$.
- Reflect the cycle through $\min E$, fixing it and its antipode: $2k-1$ transpositions, same monomial ($K$ symmetric), opposite sign.

[drv; chk (H): symbolic for $2L=4,8$; rational for $2L=12$; control with the sign dropped is nonzero]

**Checks on graphs** [chk t1_rho_square.log (G)]:
- Every admissible graph with $D=H+T\le4$, $n\le3$, at 2 rational points each: ${\rm Pf}\wedge{\rm Pf}=0$ as a form.
- Counts for $n=1,2,3$: $D=1$: 1,1,1; $D=2$: 1,1,3; $D=3$: 1,2,16; $D=4$: 1,3,42. These equal results/graph_counts.txt where listed.
- Odd $n$: ${\rm Pf}=0$, and the spanning-tree fiber integral is $0$.
- $D\le3$, $n=4$ (1+13+230 graphs), 1 point each.

## 4. Vanishing of bulk loop kernels

$$\kappa=0:\quad \Omega_\Gamma=U^{-H}\rho_n^{\wedge T}\wedge\prod_{i,r}\psi^i_r=0\quad (T\ge2,\ n\ge1)\qquad\text{and}\qquad (T=1,\ n\ \text{odd}:\ {\rm Pf}_{\rm odd}=0,\ \langle\textstyle\prod\psi\rangle\overset{\tau\to-\tau}=(-1)^n\langle\textstyle\prod\psi\rangle).$$
[drv from 3.2; src WW Props 3.23, 3.26 and GKW l.1509–1513]
[chk t1_vanishing.log: all admissible graphs with $D\le4$, $n\le3$, every splitting $H+T=D$. Controls: $T=1$, $n$ even, $\kappa=0$ is nonzero; $\kappa\ne0$ is nonzero for $T\ge2$ and for odd $n$]

## 5. $H=0$, $T=1$: bananas and Moyal–Weyl

$D=1$ forces $|V|=2$, $|I|=n+1$ [src GKW (3.10), l.1223]. Parallel move ($k=N=n+1$, $\Gamma'$ = one edge, source $\kappa$ at vertex 0):
$$\mathcal T_N[\kappa]=(-1)^ne^{t_{\rm red}\kappa^2}\,\mathbb E\big[(\sqrt{t_{\rm red}}\,\kappa+\xi)^n\big]\prod_{i=1}^Nr_i^{-1/2}\,dr_2\wedge\cdots\wedge dr_N,\quad \xi\sim\mathcal N(0,\tfrac12),\quad U=\frac{t_{\rm red}^n}{\prod r_i},\ P=t_{\rm red}\kappa^2$$
[chk t3_banana.log (a), $N=2..7$]

htk density at $\kappa=0$: $D(t)=\frac{(n-1)!!}{2^{n/2}}t_{\rm red}^{N/2}\prod_et_e^{-3/2}$ for $n$ even, $0$ for $n$ odd [chk (b): sign and square].

$$K_N=\frac{\pi^{N/2}}{\Gamma(N/2)}\frac{(n-1)!!}{2^{n/2}}=\pi^{n/2}\ (n\ \text{even}),\ 0\ (n\ \text{odd});\qquad I_N=(4\pi)^{-n/2}K_N=\begin{cases}2^{1-N}&N\ \text{odd}\\0&N\ \text{even}\end{cases}$$
[chk (c), $N\le7$; src GKW (3.28), l.1234–1240; chk-num (d): $N=3$ lattice gives $3.141553$, $N=5$ Monte Carlo gives $9.95\pm0.41$]

Sources and slices:
- On the slice $t_{\rm red}=s$: $K_N(\kappa;s)=\frac{\pi^{N/2}}{\Gamma(N/2)}e^{s\kappa^2}\mathbb E[(\sqrt s\,\kappa+\xi)^n]$. The orientation is $o_N=(-1)^n$ [chk (g)].
- For $\kappa\ne0$ the form is not basic, so the integrated value depends on the slice.
- On the UV faces of $[\epsilon,\infty)^N$, $t_{\rm red}\le\epsilon$, so $\lim_{\epsilon\to0}=K_N$, independent of $\kappa$. This equals $\int d\tau\,e^{\kappa\tau}\partial_\tau P^N$, the jump at $0$ [drv].

Normalization chain:
- $P=\tfrac12{\rm sign}$ [src GKW l.1227] → $\Theta$ (3.15) → $\alpha$ (3.53) $=(4\pi)^{-n/2}\rho_n$ (§1) → parallel move → $I_N$.
- Brackets: $\{f,g\}=\sum_N\hbar^{N-1}\frac{I_N}{N!}\eta^{\otimes N}$ [src GKW (4.23), l.1828–1835]. Coefficients $I_N/N!$: $1,0,\tfrac1{24},0,\tfrac1{1920},0,\tfrac1{322560}$ [chk (e)].
- Ordered product, with $\theta=\tfrac12+P$:
$$I^\star_N={\rm jump}(\theta P^N)=\tfrac12I_N+I_{N+1}=2^{-N}\quad (N\le6),\qquad f\star g=\sum_N\frac{I^\star_N}{N!}\hbar^N\eta^{\otimes N}(\partial^Nf,\partial^Ng)=m\circ e^{\frac\hbar2\eta}(f\otimes g)$$
[drv; chk (e); src GKW l.1842–1846]
- With $\hbar\eta^{ij}=i\hbar\Pi^{ij}$ (i.e. $[x^i,x^j]_\star=i\hbar\Pi^{ij}$) the coefficient is $\frac1{N!}(\frac{i\hbar}2)^N\Pi^{\otimes N}$ [drv].

## 6. Stretch (b): Poisson sigma model on the upper half plane, heat-kernel gauge with images

Option (a), 4d Chern–Simons, was not done.

$$\frac{d\phi(p,h)}{2\pi}=\int_0^\infty dt\,\big[\psi(h-p;t)-\psi(h-\bar p;t)\big],\quad \psi(Z;t)=\frac{e^{-|Z|^2/t}}{2\pi t^2}(Z^1dZ^2-Z^2dZ^1),\quad \phi(p,h)=\arg(h-p)-\arg(h-\bar p)$$
[drv; $\phi$: src BK Def. 1 / Kontsevich §6.2; images: src GKW l.2394–2397, sign $-$ at the tail]

From the half plane to the plane:
- $\phi(\bar p,h)=-\phi$ and $\phi(p,\bar h)=\phi$; every aerial vertex has two outgoing edges, so $\int_{\mathbb H}=\int_{\mathbb R^2}\frac{{\rm sgn}\,y_v}2$.
- $\frac{{\rm sgn}\,y}2=\int_0^\infty\frac{y\,e^{-y^2/\tau}}{2\sqrt\pi\,\tau^{3/2}}d\tau$: an auxiliary $T=1$ edge from $v$ to ground. [drv; chk-num t4_psm.log (2)]

Circuits per image pattern:
- Parallel direction: ordinary incidence with Dirichlet data $x_L=0$, $x_R=1$ (factor $e^{-1/R_{LR}}$).
- Perpendicular direction: signed incidence (an image edge has incidence $(+1,+1)$), ground at $y=0$.
- $\omega=\omega_\parallel\wedge\omega_\perp$. With $t=\lambda u$: $J=\lambda^{-a}e^{-E/\lambda}(P_0+\lambda P_1)$, $P_1=0$ [chk-num (1)], and $I_\Gamma=\Gamma(a-2n)\int_{\Delta}P_0E^{2n-a}$.

Wedge (tree): $J_\tau=\frac{(t_Lt_R)^{-1/2}}{2\pi^{3/2}R^{3/2}}e^{-1/R}$ with $R=t_L+t_R$, hence $I=\frac{1}{2\pi^{3/2}}B(\tfrac12,\tfrac12)\Gamma(\tfrac12)=\tfrac12$ [drv; chk-num]. Product graph = wedge$^2=\tfrac14$ [drv; chk-num].

| graph (BK encoding) | aerial part | $I_\Gamma$ (quadrature) | $c_\Gamma=\frac{\#\text{labellings}}{n!}I_\Gamma$ | BK App. A.2 |
|---|---|---|---|---|
| $h^1$ `0 1` | tree | 0.4999965 | 0.9999931 | 1 |
| $h^2$ `0 1 0 1` | forest | 0.2499978 | 0.4999955 | 1/2 |
| $h^2$ `0 3 0 1` | tree | 0.0833327 | 0.3333307 | 1/3 |
| $h^2$ `1 3 0 1` | tree | $-0.0833327$ | $-0.3333307$ | $-1/3$ |
| $h^2$ `1 3 0 2` | one loop (2-cycle) | $-0.0416664$ | $-0.1666658$ | $-1/6$ |

[chk-num t4_psm.log (3): Hepp sectors plus tanh-sinh, relative deviation $\le9\times10^{-6}$. The values $I=\tfrac12,\tfrac14,\tfrac1{12},-\tfrac1{12},-\tfrac1{24}$ are identified by agreement with BK]

## 7. Not done, open, suggested edits elsewhere

- [open] Option (a), 4d CS with two Wilson lines. In task 4, the weights other than the wedge and the product graph come from floating-point quadrature, not exact arithmetic.
- [open] The holomorphic factor $\rho_\Gamma[\lambda;z]$ of GKW (3.52) was not compared with $U^{-1}\prod_r\psi_r$; only $\alpha_\Gamma$ is fixed.
- [open] The star–triangle move at $\kappa=0$ has no top-form local factor, so there is no local proof of §3 for star–triangle-reducible non-SP graphs; 3.2 and 3.3 cover them.
- [open] Graph enumeration: $D=4$, $n=4$ and $n\ge5$ were not enumerated; the skeleton enumerator handles at most 6 skeleton vertices.
- [open] Integrated kernels with $\kappa\neq0$ depend on the slice (the forms are not basic); only the form identities, the banana slice formula and the UV-face argument [drv] are given.
- Suggested edits, not made:
  - CONV §2: add $\alpha_\Gamma=(-1)^n(4\pi)^{-n/2}\rho_n(M_C)$, $\det G=1$, and the hafnian form.
  - PLAN.md E8: the deliverable "$\alpha_\Gamma$ in ratio coordinates" is covered by §2.1 and §3.
  - Compendium: banana kernels $K_N=\pi^{n/2}$ and $I^\star_N=2^{-N}$.
