# Circuit reduction at $D=H+T=2$, $T=1$: star–triangle moves, families, 3d $\mathcal N=2$ HT classification

Tags: [chk] exact check run (script/log), [drv] derived, not independently checked, [src] source (location), [open]. Float Monte Carlo (MC) values are oracle estimates and are tagged [open]. Paths relative to `explore/circuit-general/ht-d2/` unless stated. Engine: `code/htk/star.py`; test `code/tests/test_star.py` (log `test_star.log`).

## 0. Conventions

$$\Omega_\Gamma=U^{-H}\rho_n(M)^T\wedge\prod_{i=1}^H\prod_{r=1}^n\psi^i_r=D_\Gamma\,\iota_E{\textstyle\prod_e}dt_e,\quad M=CtC^{\mathsf T},\ U=\det M,\ \rho_n=U^{-1/2}\mathrm{Pf}\big(-\tfrac12dM\,M^{-1}dM\big)$$ [src `code/htk/CONVENTIONS.md` §2]

$$t_ey_e=\Phi_e:=\Phi_{e(0)}-\Phi_{e(1)},\qquad \psi_r=\sum_eC_{re}\,y_e\,dt_e=\sum_eC_{re}\,\Phi_e\,d\log t_e$$ [drv]

$$\mathrm{Pf}\big(-\tfrac12dMM^{-1}dM\big)=\mathbb E_X\Big[\prod_r\psi^\tau_r\Big],\quad \psi^\tau_r=\sum_eC_{re}X_e\,d\log t_e,\quad X\sim Z^{-1}e^{-X^{\mathsf T}\Lambda X},\ \Lambda=B^{\mathsf T}t^{-1}B$$
[src `explore/fans-ht/VERIFY.md` §2.3, $\kappa=\Lambda X\sim N(0,\Lambda/2)$; chk `star.density_kron` = brute, `test_star.log`]

$$\Omega^{(1,1)}[\lambda]=U^{1/2}\,\mathbb E_\kappa\big[\Omega^{(2,0)}[\kappa,\lambda]\big]$$ [src `explore/fans-ht/VERIFY.md` §2.3]

Move coordinates: fractions $f$ ($k=2$ or $4$) and all reduced times $T$ (no gauge); $t(f,cT)=c\,t(f,T)$:
$$\iota_E{\textstyle\prod_e}dt_e=J_{\rm aff}\,df_1\wedge\dots\wedge df_k\wedge\iota_{E'}{\textstyle\prod_{e'}}dT_{e'},\quad J_{\rm aff}=\det\frac{\partial(t_0,\dots,t_{N-1})}{\partial(f,T)},\quad D':=D_\Gamma J_{\rm aff},\quad K_\Gamma=\mathrm{sgn}(J_{\rm aff})\!\int\!D'\,d\mu$$ [drv; chk `selfcheck`]

## 1. Star–triangle moves at $T=1$

### 1.1 Single move: maps (LLC Theorem b-01159 configuration)
$v\ne$ root, $\deg v=3$, neighbours $u,w,x$; $ux,wx\in I$, $uw\notin I$. $\Gamma^{\rm red}$: $V\setminus\{v\}$, edges $(uw,ux,wx,\text{rest})$, orientations $u\to w$, $u\to x$, $w\to x$.

$$R_{uw}=T_{uw},\ R_{ux}=\frac{T_{ux}}{1-r_1},\ R_{wx}=\frac{T_{wx}}{1-r_2},\ P=R_{uw}+R_{ux}+R_{wx},$$
$$t_{vu}=\frac{R_{uw}R_{ux}}P,\ t_{vw}=\frac{R_{uw}R_{wx}}P,\ t_{vx}=\frac{R_{ux}R_{wx}}P,\ t_{ux}=\frac{T_{ux}}{r_1},\ t_{wx}=\frac{T_{wx}}{r_2}$$
[src LLC eq-sum-star-triangle (inverted); chk forward $\circ$ inverse $=$ id]

$$\Sigma=t_{vu}t_{vw}+t_{vw}t_{vx}+t_{vx}t_{vu},\ (R_{uw},R_{ux},R_{wx})=\Sigma\,(t_{vx}^{-1},t_{vw}^{-1},t_{vu}^{-1}),\ T_{ux}^{-1}=t_{ux}^{-1}+R_{ux}^{-1},\ r_1=\frac{T_{ux}}{t_{ux}},\ T_{wx}^{-1}=t_{wx}^{-1}+R_{wx}^{-1},\ r_2=\frac{T_{wx}}{t_{wx}}$$ [src LLC]

$$\det\frac{\partial\log(t_{vu},t_{vw},t_{vx})}{\partial\log(R_{uw},R_{ux},R_{wx})}=\det\Big(\begin{smallmatrix}1&1&0\\1&0&1\\0&1&1\end{smallmatrix}-\mathbf 1p^{\mathsf T}\Big)=-2\big(1-\tfrac12{\textstyle\sum}p_i\big)=-1,\qquad p_i=\frac{R_i}P$$ [drv]

$$J_{\rm loc}=\det\frac{\partial(t_{ux},t_{wx},t_{vu},t_{vw},t_{vx})}{\partial(r_1,r_2,T_{uw},T_{ux},T_{wx})}=-\frac{T_{uw}T_{ux}^2T_{wx}^2}{r_1^2r_2^2(1-r_1)^3(1-r_2)^3P^3}=-\frac{T_{uw}T_{ux}^2T_{wx}^2}{r_1^2r_2^2Q^3},$$
$$Q:=(1-r_1)(1-r_2)P=(1-r_1)(1-r_2)T_{uw}+(1-r_2)T_{ux}+(1-r_1)T_{wx},\qquad J_{\rm aff}=\mathrm{sgn}(\pi_{\rm rows})J_{\rm loc}$$
[drv; chk `SingleMove.J_closed` = exact AD Jacobian, `test_star.log`]

### 1.2 Kron decomposition (any $H$; $T\in\{0,1\}$)
$$X^{\mathsf T}\Lambda X=G\,(X_v-\bar X_v)^2+X_r^{\mathsf T}\Lambda_{\rm red}X_r,\quad G=\sum_at_{va}^{-1},\quad \bar X_v=\sum_a\frac{X_a}{t_{va}G},\quad \Lambda_{\rm red}=\Lambda_{rr}-\frac{gg^{\mathsf T}}G\ \ (g_a=t_{va}^{-1})$$
$\Lambda_{\rm red}$ = Laplacian of $\Gamma^{\rm red}$ (triangle conductances $g_ag_b/G=R_{ab}^{-1}$, merged with $ux,wx$) [drv]
$$X_v=\bar X_v+\xi,\qquad \xi\sim N\Big(0,\frac1{2G}\Big)\ \perp\ X_{\rm red}\sim Z^{-1}e^{-X^{\mathsf T}\Lambda_{\rm red}X},\qquad \frac1{2G}=\frac{T_{uw}R_{ux}R_{wx}}{2P^2}$$ [drv; chk `density_kron` = brute]
$$\Phi_v=\bar\Phi_v+\frac{\lambda_v}G,\qquad \lambda^{\rm red}_a=\lambda_a+\frac{\lambda_v}{t_{va}G}\ \ (a=u,w,x)$$ [src LLC app-sourced-star; chk `SingleMove.currents` = brute]
$$\det\Lambda_\Gamma=G\det\Lambda_{\rm red}\ \Rightarrow\ U_\Gamma=U_{\rm red}\,G\,\frac{\prod_\Gamma t}{\prod_{\rm red}T}=U_{\rm red}\frac{t_{vx}}{r_1r_2}=U_{\rm red}\frac{T_{ux}T_{wx}}{r_1r_2\,Q}$$ [drv; chk]

### 1.3 $M$ in an adapted basis
$B'=(\gamma_1,\gamma_2,c'_1,\dots,c'_{n-2})$: $\gamma_1=ux-vx+vu$, $\gamma_2=wx-vx+vw$; $c'$ = cycles of $\Gamma\setminus v=\Gamma^{\rm red}\setminus uw$ ($ux^{\rm red}\to ux$, $wx^{\rm red}\to wx$); reduced basis $(\delta,c)$, $\delta=uw+wx-ux$ [drv]
$$M_{\gamma\gamma}=\begin{pmatrix}t_{ux}+t_{vx}+t_{vu}&t_{vx}\\t_{vx}&t_{wx}+t_{vx}+t_{vw}\end{pmatrix},\quad M_{\gamma_1c'}=C_{c,ux}t_{ux},\quad M_{\gamma_2c'}=C_{c,wx}t_{wx},$$
$$M_{c'c''}=M^{\rm red}_{cc''}+C_{c,ux}C_{c'',ux}\frac{(1-r_1)T_{ux}}{r_1}+C_{c,wx}C_{c'',wx}\frac{(1-r_2)T_{wx}}{r_2}$$ [drv]
$$\det M_{\gamma\gamma}=M^{\rm red}_{\delta\delta}\frac{t_{vx}}{r_1r_2},\qquad M_{cc}-M_{c\gamma}M_{\gamma\gamma}^{-1}M_{\gamma c}=M^{\rm red}_{cc}-\frac{M^{\rm red}_{c\delta}M^{\rm red}_{\delta c}}{M^{\rm red}_{\delta\delta}}$$ [chk `check_components`, all forests with a site, $n\le5$]

### 1.4 Current forms
$$t_{vu}=Pp_1p_2,\ t_{vw}=Pp_1p_3,\ t_{vx}=Pp_2p_3,\ t_{ux}=\frac{Pp_2(1-r_1)}{r_1},\ t_{wx}=\frac{Pp_3(1-r_2)}{r_2},\ T_{uw}=Pp_1,\ T_{ux}=Pp_2(1-r_1),\ T_{wx}=Pp_3(1-r_2)$$ [drv]
$a=\Phi_u-\Phi_x$, $b=\Phi_w-\Phi_x$ (reduced potentials); source-free part $\Phi_v-\Phi_x=ap_3+bp_2$ [drv]
$$\psi_{\gamma_1}=a\,d\log t_{ux}-(ap_3+bp_2)\,d\log t_{vx}+(ap_3+bp_2-a)\,d\log t_{vu}=-\frac{a\,dr_1}{r_1(1-r_1)}+p_2\,\omega,$$
$$\psi_{\gamma_2}=-\frac{b\,dr_2}{r_2(1-r_2)}-p_3\,\omega,\qquad \omega=a\,d\log\frac{p_2}{p_1}+b\,d\log\frac{p_1}{p_3}=-\psi^{\rm red}_\delta+\frac{a\,dr_1}{1-r_1}-\frac{b\,dr_2}{1-r_2},$$
$$\psi_{c'}=\psi^{\rm red}_c-C_{c,ux}\,a\,\frac{dr_1}{r_1}-C_{c,wx}\,b\,\frac{dr_2}{r_2}$$
$\Phi_v\to\Phi_v+s$:
$$\psi_{\gamma_1}\mathrel{+}=s\,\alpha,\quad \psi_{\gamma_2}\mathrel{+}=s\,\beta,\quad \psi_{c'}\mathrel{+}=0,\qquad \alpha=d\log\frac{t_{vu}}{t_{vx}}=d\log\frac{T_{uw}}{T_{wx}}-\frac{dr_2}{1-r_2},\quad \beta=d\log\frac{t_{vw}}{t_{vx}}=d\log\frac{T_{uw}}{T_{ux}}-\frac{dr_1}{1-r_1}$$
[drv; chk `check_components`: each component, $H=1,2,3$, sourced centre $s=\lambda_v/G$, against the pull-back of the brute-force $\psi$]
$H=1$:
$$\psi_{\gamma_1}\wedge\psi_{\gamma_2}=\frac{ab\,(1-p_2r_1-p_3r_2)}{r_1r_2(1-r_1)(1-r_2)}dr_1\wedge dr_2-\Big(\frac{p_3\,a\,dr_1}{r_1(1-r_1)}+\frac{p_2\,b\,dr_2}{r_2(1-r_2)}\Big)\wedge\psi^{\rm red}_\delta$$ [drv]

### 1.5 $\rho_n$
$$\mathrm{Pf}_\Gamma:=\mathrm{Pf}\big(-\tfrac12dMM^{-1}dM\big)_{B'}=\mathbb E_{X_{\rm red},\xi}\Big[\prod_{B'}\psi^\tau\Big]=\mathbb E_{X_{\rm red}}\Big[\prod_{B'}\psi^\tau\Big]_{\xi=0}+\frac1{2G}\,\alpha\wedge\beta\wedge\mathbb E_{X_{\rm red}}\Big[\prod_c\psi^\tau_{c'}\Big]$$
$$\rho_n(M_\Gamma)=U_\Gamma^{-1/2}\mathrm{Pf}_\Gamma=U_{\rm red}^{-1/2}\Big(\frac{r_1r_2Q}{T_{ux}T_{wx}}\Big)^{1/2}\mathrm{Pf}_\Gamma$$ [drv; chk `check_components` ($n$-form identity after pull-back)]

### 1.6 Density; local factor
$$D'=(\det A)^{H+1}\,U_\Gamma^{-H-\frac12}\Big[\mathrm{Pf}_\Gamma\wedge\prod_{i=1}^H\prod_{B'}\psi^i\Big]_{df\wedge\iota_{E'}\prod dT},\qquad B'=A\,C_{\rm pos}$$
[drv; chk `density_closed` = brute: $T=1$, $H=1$; $T=0$, $H=2$ (no $\rho$, $U^{-H}$); sourced; `test_star.log`]

$H=1$, $\lambda_v=0$, $\xi^0$ part; from $\Omega^{(1,1)}=U^{1/2}\mathbb E_\kappa\Omega^{(2,0)}$ and $d^2y_{ux}\wedge d^2y_{wx}\equiv r_1r_2\det(y^{\rm red}_{ux},y^{\rm red}_{wx})\,dr_1\wedge dr_2\wedge d^2y^{\rm red}_{uw}$:
$$D'\big|_{\xi^0}=(-1)^{n-1}(r_1r_2t_{vx})^{1/2}\,D'_{\rm red}\big[w_{ux}y_{wx}-y_{ux}w_{wx}\big],\qquad (r_1r_2t_{vx})^{1/2}=\Big(\frac{r_1r_2T_{ux}T_{wx}}Q\Big)^{1/2}$$
$$U_{\rm red}^{-3/2}\,\mathbb E_{X_{\rm red}}\Big[F\prod_r\psi^{{\rm red},\tau}_r\Big]\wedge\prod_r\psi^{\rm red}_r[\lambda]=D'_{\rm red}[F]\,\iota_{E'}{\textstyle\prod}dT,\qquad w_e=\frac{X_e}{T_e}$$
[src LLC eq-star-current-wedge; drv sign $\epsilon(n,2)\epsilon(n-1,2)=(-1)^{n-1}$; chk `SingleMove.xi0_T1_identity`]
Contact part:
$$D'-D'\big|_{\xi^0}=(\det A)^2\,U_\Gamma^{-3/2}\Big[\frac1{2G}\,\alpha\wedge\beta\wedge\mathbb E_{X_{\rm red}}\Big[\prod_c\psi^\tau_{c'}\Big]\wedge\prod_{B'}\psi[\lambda]\Big]_{\rm top}$$ [drv; chk as the difference of two checked quantities]

$T=0$, $H=2$, $\lambda_v=0$:
$$D'=(-1)^{n-1}\,r_1r_2\det(y^{\rm red}_{ux},y^{\rm red}_{wx})\,D_{\rm red}$$ [chk `density_llc_T0`].
With LLC's edge order $(vu,vw,vx,ux,wx,\dots)$: $\mathrm{sgn}J_{\rm aff}=-1$, so $K_\Gamma=(-1)^n\!\int r_1r_2\langle y_{ux},y_{wx}\rangle(\cdot)$; this is LLC eq-kernel-star-triangle with $I^{\rm S1}=(-1)^{n(n+1)/2}K$ [drv; src CONVENTIONS §3].

Example: fan(2), $v$ = middle vertex, $\Gamma^{\rm red}$ = triangle:
$$D_\Gamma=\frac{\lambda_0\lambda_2}{2U^{5/2}}\ \ (\text{any }\lambda_1),\qquad |D'|=\frac{\lambda_0\lambda_2}2\,U_{\rm red}^{-5/2}\,T_{uw}(T_{ux}T_{wx})^{-1/2}(r_1r_2)^{1/2}Q^{-1/2}$$
[chk `bitriangle_density.log` (symbolic, matches brute); drv second equality; cf. src Higher.tex (3.59), $\omega\propto U^{-5/2}$]
$$\omega_{\rm bitriangle}^{\rm Higher\,(3.59)}=\pm\frac{\pi^2}2\,\Omega^{\rm htk}_{\rm fan(2)}$$ [drv]

$Q$ is not a Laurent monomial in $(r,T)$, so LLC's Dirichlet scheme (polynomial density on intervals × simplices) does not apply to star–triangle charts at $T=1$ [drv]. Exact $T=1$ kernels of non-series–parallel graphs [open].

### 1.7 Sourced single move
$$\psi_{\gamma_1}[\lambda]=\psi_{\gamma_1}[\lambda^{\rm red}]+\frac{\lambda_v}G\alpha,\quad \psi_{\gamma_2}[\lambda]=\psi_{\gamma_2}[\lambda^{\rm red}]+\frac{\lambda_v}G\beta,\quad \psi_{c'}[\lambda]=\psi_{c'}[\lambda^{\rm red}],\quad \frac{\lambda_v}G=\lambda_v\frac{T_{uw}R_{ux}R_{wx}}{P^2}$$
$\mathrm{Pf}_\Gamma$ unchanged; $D'$ of §1.6 [chk `check_components`, `density_closed`, $T=1$ $H=1$ and $T=0$ $H=2$]

### 1.8 Paired move (LLC Theorem b-01173 configuration)
$$R^1_{uA}=\frac{T_{uA}}{r_1},\ R^2_{uA}=\frac{T_{uA}}{1-r_1},\ R^1_{uB}=\frac{T_{uB}}{r_2},\ R^2_{uB}=\frac{T_{uB}}{1-r_2},\ R^1_{AB}=\frac{T_{AB}}{r_3},\ R^2_{AB}=\frac{T_{AB}}{r_4},\ t_{AB}=\frac{T_{AB}}{1-r_3-r_4},$$
$$t_{v_iu}=\frac{R^i_{uA}R^i_{uB}}{P_i},\quad t_{v_iA}=\frac{R^i_{uA}R^i_{AB}}{P_i},\quad t_{v_iB}=\frac{R^i_{uB}R^i_{AB}}{P_i}$$ [src LLC b-01180; chk]
$$|J_{\rm aff}|=\frac{T_{uA}T_{uB}T_{AB}^2}{r_1^2(1-r_1)^2r_2^2(1-r_2)^2r_3^2r_4^2(1-r_3-r_4)^2}\prod_{i=1,2}\frac{R^i_{uA}R^i_{uB}R^i_{AB}}{P_i^3},\qquad \frac{U_\Gamma}{U_{\rm red}}=\frac{T_{uA}T_{uB}T_{AB}^2}{r_1(1-r_1)r_2(1-r_2)r_3r_4(1-r_3-r_4)P_1P_2}$$ [drv; chk `J_abs_closed`, `U_ratio_closed`]
Kron: $X_{v_i}=\bar X_{v_i}+\xi_i$, $\xi_i\sim N(0,\frac1{2G_i})$, independent ($v_1,v_2$ non-adjacent) [drv; chk `density_kron`, $T=1$ $H=1$ sourced, $T=0$ $H=2$]
$$D'\big|_{\xi^0}=-\Big(\frac{U_\Gamma}{U_{\rm red}}\Big)^{1/2}(r_1-r_2)(1-r_3-r_4)\,D'_{\rm red}\big[(w_{uA}y_{AB}-y_{uA}w_{AB})(w_{uB}y_{AB}-y_{uB}w_{AB})\big]$$ [chk `PairedMove.xi0_T1_identity`]
$T=0$, $H=2$:
$$D'=-(r_1-r_2)(1-r_3-r_4)\det(y_{uA},y_{AB})\det(y_{uB},y_{AB})\,D_{\rm red}$$ [chk `density_llc_T0`; src LLC eq-kernel-two-leaves].
Closed $\psi$ formulas for the paired move [open].

### 1.9 Odd $n$; $H\ge2$
$n=3$ forests: $D_\Gamma\equiv0$ at $T=1$ ($\rho_3=0$) [chk `test_star.log`; src].
§1.2–1.5 are $H$-independent [chk $H=1,2,3$].
$N^{\rm red}-1-n^{\rm red}D=D-2$, so $\Gamma^{\rm red}$ is tight only for $D=2$ ($H=1$ at $T=1$) [drv].
No tight graph with $(D,n)\in\{(3,2),(3,3),(4,2)\}$ has a single or paired site [chk `general_h_sites.log`].

### 1.10 Checks (exact Fractions)
`python3 code/tests/test_star.py`: 138 (graph, centre) single sites and 17 paired sites on fan(2), fan(4), all forests with $n=3,4,5$, and spiders $(1,1,2),(1,2,1,2)$. Covered: maps, currents, $U$, $J$, Kron densities, closed densities, $\psi$/Pf/$M$ components, LLC $T=0$ factors, $\xi^0$ identities. Run time 8 s [chk `test_star.log`].

## 2. Families at $(H,T)=(1,1)$

### 2.1 Odd $n$
$K_\Gamma\equiv0$ ($\rho_{2k+1}=0$) [src CONVENTIONS, ht-families]; spiders with $n=3,5$, all arm lists, have $K=0$ [chk `spiders_exact.log`].

### 2.2 Identifications
- spider$(\ell_1..\ell_a)$ = product forest of $a$ B-chains $(B,(B,\dots,(A,B)))$, hub $=A$ [drv; chk spider$(1,1,2)$ = $(B,(A,B))(A,B)(A,B)$, `spiders_exact.log`, `probe_forests.log`].
- $(B,(B,\dots(A,B)))$ and $(A,(A,\dots(A,B)))$ are fans: $K=\pm C_n(\lambda_0\lambda_A)^{n/2}$, resp. $\pm C_n(\lambda_0\lambda_B)^{n/2}$, with $C_n=\frac{(4\pi)^{n/2}}{(n+1)(n!)^2}$ [src fans-ht a7; chk $n=2,4$].
- Arms $(1^n)$ = forest $(A,B)^n$ = banana: $K=(-1)^{n/2}\frac{\pi^{n/2}}{n+1}\prod_{i=1}^n\lambda_i$, independent of $\lambda_A$ [chk $n=2,4$ (`probe_forests.log`, `spiders_exact.log`, `spiders_n6_numeric.log`); modulus for all even $n$: src ht-families VERIFY a5, $\frac{(n-1)!!}{2^{n/2}}\frac{\pi^{n/2}\Gamma(3/2)}{\Gamma((n+3)/2)}=\frac{\pi^{n/2}}{n+1}$ drv; sign for $n\ge6$ open].

### 2.3 Spiders, $n=4$
$Q_n$ sources $l_i$ at $v_{i,1}$, $l_A$ at $A$, $B$ = root:
$$K_{(1,1,2)}=-\frac{\pi^2}{90}l_1l_2l_3(3l_1+3l_2-2l_3+6l_A),\qquad K_{(1,1,1,1)}=\frac{\pi^2}5l_1l_2l_3l_4$$ [chk `spiders_exact.log`]
$$K_1^{(1,1,1,1)}=-\frac{\pi^2}{70}l_1l_2l_3l_4\Big[\sum_il_i\big(7z_i+{\textstyle\sum_j}z_j+3z_A\big)+2l_A\big({\textstyle\sum_j}z_j+3z_A\big)\Big]$$
$z_i$ is the displacement on spoke $(v_{i,1},B)$ and $z_A$ on $(A,B)$; $K_1,K_2$ for both spiders [chk `taylor_sp.log`].

### 2.4 Product forests: $n=2$ (3 forests), $n=4$ (23: 10 series–parallel, 13 not)
Counts $1,3,8,23,65,194$ ($n=1..6$) [chk `forests.py`; = LLC].
$n=2$: $(A,B)(A,B)$: $-\frac\pi3l_0l_1$; $(A,(A,B))$: $\frac\pi3l_0(l_0+l_A)$; $(B,(A,B))$: $\frac\pi3l_0l_A$ [chk `probe_forests.log`].

$n=4$, series–parallel: root sources $l_v$, $\lambda_A=l_A$, lexicographic edge order with $(A,B)$ last [chk `probe_forests.log`]:

| forest | $K/\pi^2$ |
|---|---|
| $(A,B)^4$ | $\frac15l_0l_1l_2l_3$ |
| $(A,(A,B))(A,B)^2$ | $-\frac1{90}l_0l_2l_3(8l_0+3l_2+3l_3+6l_A)$ |
| $(A,(A,B))^2$ | $\frac1{180}l_0l_2(10l_0^2+21l_0l_2+22l_0l_A+10l_2^2+22l_2l_A+12l_A^2)$ |
| $(B,(A,B))(A,B)^2$ | $\frac1{90}l_0l_2l_3(2l_0-3l_2-3l_3-6l_A)$ |
| $(B,(A,B))(A,(A,B))$ | $\frac1{180}l_0l_2(l_0l_2+2l_0l_A+10l_2^2+22l_2l_A+12l_A^2)$ |
| $(B,(A,B))^2$ | $\frac1{180}l_0l_2(l_0l_2+2l_0l_A+2l_2l_A+12l_A^2)$ |
| $(A,(A,(A,B)))(A,B)$ | $\frac1{180}l_0^2l_3(2l_0+l_3+2l_A)$ |
| $(B,(B,(A,B)))(A,B)$ | $-\frac1{180}l_0^2l_3(l_3+2l_A)$ |
| $(A,(A,(A,(A,B))))$ | $-\frac1{180}l_0^2(l_0+l_A)^2$ |
| $(B,(B,(B,(A,B))))$ | $\frac1{180}l_0^2l_A^2$ |

Symmetric forks: swapping two isomorphic, source-free sibling subtrees ($m$ vertices, $2m+1$ edges each) is an odd edge permutation. At $D=2$: $K\to\mathrm{sgn}(\pi)^{D+1}K$ and vertex relabelling gives $\mathrm{sgn}^D=+1$, so $K(\lambda;0)=0$ [drv; src CONVENTIONS §2.2; LLC eq-product-four-loop-evaluated at $(2,0)$]. At $n=4$: $((A,B),(A,B))(A,B)$, $(A,((A,B),(A,B)))$, $(B,((A,B),(A,B)))$; MC consistent with 0 [open, `mc_summary.log`].

Non-series–parallel, MC (`mc_chain.py`: chain of star–triangle moves then a series–parallel chart, float brute density, $2\cdot10^6$ points; calibration $(A,(A,B))(A,B)^2$: $0.66647\pm0.00204$ vs $2/3$) [open, `mc_summary.log`]:
- Two-tree forests, root sources $(2,-1)$, $\lambda_A=-2$:
  - $(A,(B,(A,B)))(A,B)$: $K/\pi^2=0.10195\pm0.00042$ (single-move chart: $0.10111\pm0.00049$)
  - $(B,(A,(A,B)))(A,B)$: $K/\pi^2=-0.06289\pm0.00034$ (single-move chart: $-0.06053\pm0.00218$)
- Single trees, source polynomial $K/\pi^2=\sum_kc_kl_0^kl_A^{4-k}$ (`run_mc_poly.py`; calibration $(A,(A,(A,(A,B))))$: $c_{2,3,4}=-0.00558,-0.01116,-0.00558$ vs $-\frac1{180},-\frac2{180},-\frac1{180}$):

| tree | $c_1$ | $c_2$ | $c_3$ | $c_4$ |
|---|---|---|---|---|
| $(A,(A,(B,(A,B))))$ | $-0.01076(5)$ | $-0.03339(17)$ | $-0.03290(17)$ | $-0.01026(6)$ |
| $(A,(B,(A,(A,B))))$ | $-0.00587(7)$ | $-0.00999(12)$ | $-0.00329(8)$ | $+0.00084(3)$ |
| $(A,(B,(B,(A,B))))$ | $+0.00070(0)$ | $+0.00471(2)$ | $+0.00558(3)$ | $+0.00157(1)$ |
| $(B,(A,(A,(A,B))))$ | $+0.00072(0)$ | $-0.00263(2)$ | $-0.00175(2)$ | $0$ |
| $(B,(A,(B,(A,B))))$ | $-0.00580(10)$ | $-0.00749(15)$ | $-0.00089(3)$ | $0$ |
| $(B,(B,(A,(A,B))))$ | $-0.01080(5)$ | $+0.00113(3)$ | $+0.00160(1)$ | $0$ |
| $((A,(A,B)),(A,B))$ | $+0.00247(2)$ | $-0.01473(10)$ | $-0.02826(18)$ | $-0.01263(9)$ |
| $((B,(A,B)),(A,B))$ | $-0.00244(5)$ | $-0.02237(68)$ | $-0.00878(20)$ | $-0.00158(1)$ |

$c_0=0$ in all rows. Exact non-series–parallel values [open].

### 2.5 Edge-dipole vanishing
$$\lambda_{e(0)}=-\lambda_{e(1)}=\lambda,\ \text{all other sources }0\ \Longrightarrow\ D_\Gamma(t)=0\quad(H=T=1)$$
[chk `dipole_check.log`: all tight $D=2$ graphs, $n=2$: 5/5 edges, $n=4$: 117/117 (graph, edge), exact at 2 random rational $t$ each; non-adjacent pairs: 0/1 and 52/78 vanish; general proof open]

Partial derivation, $\Gamma=e_0\parallel\Gamma'$ [drv]: $X=X''+V\varphi$, with $\varphi$ the unit harmonic potential of $\Gamma'$, $X''$ the Gaussian of $\Gamma''=\Gamma'/(p{=}q)$, $V\perp X''$. Then
$$\psi^\tau_r=\psi_r[X'']+\frac V{J'}\psi_r[\lambda]\ \Rightarrow\ \Omega_\Gamma=U_\Gamma^{-3/2}J'^n\,\mathrm{Pf}_{\Gamma''}\wedge{\textstyle\bigwedge_c}\psi'_c\wedge\frac{ds}s;$$
the vanishing of $\mathrm{Pf}_{\Gamma''}\wedge\bigwedge_c\psi'_c$ is [open].

Consequences [drv from the [chk] statement; chk at $n=4$ in `zero_density_check.log`]:
- Fans: $\lambda_1\lambda_{n+1}\mid g_n$ (consistent with src fans-ht a11).
- Single trees: $l_0\mid K$ (dipole on $AB$); $\lambda_B\mid K$ for root spoke $A$; $\lambda_A\mid K$ for root spoke $B$; this matches the MC table ($c_0=0$; $c_4=0$ for $(B,\cdot)$; $\sum_k(-1)^kc_k=0$ for $(A,\cdot)$).
- Tadpoles with $L\ge1$: $D\equiv0$ (§3.3).

### 2.6 Taylor coefficients
tadpole(3,0) ($n=2$), order 3 [chk `taylor_tadpole.log`]:
$$K_0=\frac\pi3l^2,\qquad K_1=-\frac\pi6l^3(z_{01}+z_{02}+z_{13}+z_{23})$$
Fans: all orders [src fans-ht b1]. Series–parallel spiders and forests: `spred.taylor` [chk `taylor_sp.log`].

### 2.7 $n=6$
Not done: exact `spred` runs at $T=1$, $n=6$ (symbolic and integer sources) did not finish within 1700 s [open].

## 3. Classification for 3d $\mathcal N=2$ HT theories

### 3.1 Action
Generic HT action $\int[(\Phi,\mathbb d\Phi)+\mathcal I(\Phi)]\,d^Hx^{\mathbb C}$ [src Higher.tex (1.1)]. Costello–Dimofte–Gaiotto, arXiv:2005.00083, ar5iv HTML:
$$S=\int\mathbf BF(\mathbf A)+\int\mathbf\Psi\,d_{\mathbf A}\mathbf\Phi-\int W(\mathbf\Phi)+\frac k{4\pi}\int\mathbf A\,\partial\mathbf A,\qquad d_{\mathbf A}=d'-i\mathbf A$$ [src CDG (3.32)]
- Superfields: $\mathbf A\in\mathbf\Omega^\bullet\otimes\mathfrak g[1]$, $\mathbf B\in\mathbf\Omega^{\bullet,(1)}\otimes\mathfrak g$, $\mathbf\Phi\in\oplus_r\mathbf\Omega^{\bullet,(r/2)}\otimes V^{(r)}$, $\mathbf\Psi\in\oplus_r\mathbf\Omega^{\bullet,(1-r/2)}\otimes(V^{(r)})^\vee[1]$ [src CDG (3.26)].
- BV brackets: $\{\mathbf\Phi(x),\mathbf\Psi(y)\}=\{\mathbf A(x),\mathbf B(y)\}=\delta(x-y)\,d\mathrm{Vol}$ [src (3.44)]; $Q\mathbf\Phi=\delta S/\delta\mathbf\Psi$, $Q\mathbf\Psi=-\delta S/\delta\mathbf\Phi$ [src (3.45)].
- $Q$ on superfields: $Q\mathbf\Phi=d_{\mathbf A}\mathbf\Phi$, $Q\mathbf\Psi=d_{\mathbf A}\mathbf\Psi+\partial W/\partial\mathbf\Phi$ [src (3.34)].
- Chirals only: $S=\int[\mathbf\Psi d'\mathbf\Phi-W(\mathbf\Phi)]$ [src (5.20)]; $R(W)=2$, $J(W)=1$ [src (3.3)]; propagator [src (5.7)].

### 3.2 Feynman rules
- Propagators $\langle\mathbf\Phi\mathbf\Psi\rangle$, $\langle\mathbf A\mathbf B\rangle$.
- Vertices: $\mathbf\Psi\mathbf A\mathbf\Phi$ and $\mathbf B\mathbf A\mathbf A$ (one antifield terminal each); $W$ (degree $m$, none); $\mathbf A\partial\mathbf A$ (degree 2, none) [drv from (3.32)].
- Dictionary to LLC eq-B-10: $(b,c,\beta,\gamma,\bar W)\leftrightarrow(\mathbf B,\mathbf A,\mathbf\Psi,\mathbf\Phi,W)$; the CS vertex is new [drv].
- Graphs: 2-Laman, $|V|=n+2$, $|I|=2n+1$ [src Higher.tex (3.10)–(3.12)].

### 3.3 Counting for $Q_n\mathcal O$ ($n+1$ interaction vertices)
$k=\deg\mathcal O\le n+1$, $m=2n+1-k$ internal edges, $w=\#(W,\mathrm{CS})$:
$$n\le m\le n+1-w\ \Rightarrow\ (m,k,w)\in\{(n,n+1,0),\,(n+1,n,0),\,(n,n+1,1)\}$$
$$\#(\mathbf\Psi,\mathbf B)_{\rm out}=\#(\mathbf\Psi,\mathbf B)_{\mathcal O}-n-w,\qquad Q_n\mathcal O=0\ \ (n>\#(\mathbf\Psi,\mathbf B)_{\mathcal O})$$ [drv, as LLC eq-1-15, eq-antifield-loss]

Laman on (component $\cup\,\mathcal O$): $e_j+k_j\le2c_j-1$, so the interaction graph is connected. Cubic sources: $s=3(n+1)-w-2m-k$ [drv]. Three cases:
- (a) $(n,n+1,0)$: $s=2$, tree, path: fans [drv; src LLC Prop b-00210].
- (b) $(n+1,n,0)$: $s=1$, unicyclic. Laman on $C\cup\mathcal O$ gives $k_C\le|C|-1$, so the unspoked vertex $z$ lies on $C$. One degree-two vertex $q$ carries the source; $q=z$, or $q$ ends a tail at $z$. These are tadpole$(c,L)$, $c\ge3$, $c+L=n+1$ [drv]. Abelian $\mathfrak g$: $L=0$ only [drv].
- (c) $(n,n+1,1)$: tree with one $W$ or CS vertex. $W$: spiders, $a\le\deg W-1$ [src LLC Prop b-00240]. CS: a source-free fan end, $K\equiv0$ for all $z$ ($\lambda_1\lambda_{n+1}\mid g_n$ pointwise) [src fans-ht a11(i); drv].

At $(2,0)$ case (b) vanishes ($y_e=c_e\lambda$ collinear) [src LLC Lemma b-00213]. At $(1,1)$:
$$K_{\text{tadpole}(3,0)}(l;0)=\frac\pi3l^2\ne0$$ [chk `taylor_tadpole.log`]
$n=4$ [chk `zero_density_check.log`; MC `mc_summary.log`]:
- tadpole$(3,2)$, tadpole$(4,1)$: $q$ adjacent to $\mathcal O$, $D\equiv0$ (§2.5).
- tadpole$(5,0)$: $K/\pi^2=0.00159\pm0.00001$ at $\lambda_q=1$ [open].

### 3.4 Contributing families
- $n$ odd: none [src].
- $n$ even, $\deg W\le3$: fans, plus tadpole$(n+1,0)$ ($n=2$ exact, $n=4$ MC); tadpoles with $L\ge1$ vanish at $n=4$ [chk], general $n$ [open].
- $\deg W\ge4$: add spiders.
- CS vertex: none.

"Only fans for cubic $W$" is false at kernel level [chk $n=2$]. Whether the Wick weights $W_\Gamma$ cancel over the two cycle orientations [open].

### 3.5 $Q_n$, $n$ even
$$Q_n\mathcal O=\mathcal N_n\sum_{\Gamma\in\{\text{fans},\ \text{tadpole}(n+1,0),\ \text{spiders}\}}\Big[K^{(1,1)}_\Gamma(\lambda;z)\,W_\Gamma[\mathcal I_1,\dots,\mathcal I_{n+1},\mathcal O]\Big],\qquad Q_{\rm odd}=0$$
[drv, structure of LLC eq-2-3; $\mathcal N_n$ from CDG (5.7) vs Higher (3.15) open]
$$K^{(1,1)}_{\rm fan}(\lambda;z)=(\lambda_1\lambda_{n+1})^{n/2}\!\int dr\,w(r)\,e^{\lambda_1X_1(r,z)+\lambda_{n+1}X_{n+1}(r,z)},\qquad \int w=C_n$$
$w$: (1,1) endpoint density [src fans-ht a6–a8]; $X_{1},X_{n+1}$: LLC eq-fan-output-points [src].
$\lambda_v^p\mapsto\partial_z^p$ on the surviving field at $v$; $z_e$ derivatives act on the letters of $\mathcal O$; $K(c\lambda;z/c)=c^nK(\lambda;z)$, so $Q_n$ adds $n$ derivatives $\partial_z$ [drv].
Tadpole: $K=c_n\lambda_z^n$ at $z=0$; $c_2=\pi/3$ [chk]; $c_4\approx0.00159\,\pi^2$ [open].

### 3.6 $\sum_{k+l=n}Q_kQ_l=0$, $W=\Phi^3/3$ and $W=XYZ$
$W$ vertices carry no $\mathbf\Psi$, so every edge ends at $\mathcal O$; the graph is a star, $|I|=n+1=2n+1$, hence $n=0$, and $Q_{n\ge1}=0$ [drv; chk all tight $D=2$ graphs $n\le4$ with every choice of $\mathcal O$: 0 admissible, `q_check_wz.log`].
$Q_0\psi^{(k)}=\partial^kW'(\phi)$, $Q_0\phi^{(k)}=0$. $Q_0^2=0$ on all 92 monomials ($\Phi^3/3$, weight $\le6$) and all 344 monomials ($XYZ$, weight $\le4$) [chk `q_check_wz.log`]. Hence $\sum_{k+l=n}Q_kQ_l=0$ for $n\le4$ [drv]. With gauge fields, $\{Q_0,Q_2\}=0$ [open].

## 4. Code and files
- `code/htk/star.py`:
  - classes and search: `SingleMove(G,v)`, `PairedMove(G,v1,v2)`, `find_single`, `find_paired`;
  - maps and Jacobians: `tmap`, `forward`, `jacobian`, `J_closed`, `J_abs_closed`, `U_ratio`, `U_ratio_closed`;
  - sources and currents: `red_sources`, `currents`;
  - densities: `density_kron(f,T,lam,Tdir,xi_var)`, `density_closed` (single move), `density_llc_T0`;
  - identities and checks: `xi0_T1_identity`, `check_components` (single move), `selfcheck(lam,H,Tdir)`;
  - helpers: `wick`, `reduced_insertion_density`, `Dual` (exact AD).
- Test: `code/tests/test_star.py`.
- Scripts in this folder:
  - `forests.py` (forests, spiders, tadpoles); `probe_forests.py`, `spiders_exact.py`, `spiders_n6_numeric.py`;
  - `taylor_sp.py`, `taylor_tadpole.py`; `bitriangle_density.py` with `polybrute.py`;
  - `dipole_check.py`, `zero_density_check.py`, `general_h_sites.py`;
  - `mc_chain.py`, `run_mc_chain.py`, `run_mc_poly.py`, `nonsp_numeric.py` (single-chart MC); `q_check_wz.py`.
- Logs: `*.log`; MC collected in `mc_summary.log`.
- Environment: reads of `code/htk/__pycache__/*.pyc` hang (probably iCloud-evicted files); runs from 23:38 on used `PYTHONPYCACHEPREFIX=<scratchpad>/pycache`.

## 5. Open / not done
- Exact $T=1$ kernels of non-series–parallel forests and of tadpole$(5,0)$; only MC, about 3 digits.
- Proof of §2.5 for general graphs.
- Closed $\psi$ formulas for the paired move.
- Spiders and forests at $n=6$.
- Closed forms for spiders other than arms $(1^n)$ and $(n)$.
- $\mathcal N_n$, the Wick/Koszul weights $W_\Gamma$ for gauge theories, and $\{Q_0,Q_2\}=0$ with gauge fields.
- Higher.tex numbering: a scratch compile (aux) gives (3.10) global constraint, (3.14), (3.15), (3.18) operatope, (3.52) factorization, (3.53) $\alpha_\Gamma$, (3.56) $\alpha_{\rm bitriangle}$, (3.59) $\omega_{\rm bitriangle}$, (3.60) $d\mathrm{Vol}_t$. Sections 1–2 end at (1.6) and (2.39), so (3.56) and (3.59) are the 101st and 104th equations. The cited (106), (109) are identified with them by content; the version mismatch is open.
