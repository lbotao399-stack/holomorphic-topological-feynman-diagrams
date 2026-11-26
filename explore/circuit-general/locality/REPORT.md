# Locality of Ohm's law at kernel level: Kron reduction of subnetworks

Tags: [chk] exact check (script, log in this folder), [drv] derived here, [src] source (location), [open]. P = project root; paths without prefix in `explore/circuit-general/locality/`. Engine `code/htk` and `explore/circuit-general/cores` imported, not edited.

| file | content | log |
|---|---|---|
| `locality.py` | `Split` (Kron data, derivatives, adapted cycle basis), `factor_T0`, `factor_T1` (form identities, modes R/P), `TwoTerminal` (chart identity k=2), `rho_S_T1`, `density_T1_insert`, Wick with linear Gaussian 1-forms | — |
| `check_locality.py` | form level T=0 (R, P), chart level k=2, currents, $\det M_{\rm out}P_\Delta=U_{\Gamma'}$ | `check_locality.log` |
| `check_locality_T1.py` | T=1 form level, even n; T=1 chart level with $\rho_S^{(j)}$ and $D'^{[w^j]}$ | `check_locality_T1.log` |
| `check_sourcefree.py` | vanishing for source-free tight S | `check_sourcefree.log` |
| `response_rank.py` | response dimensions; star rank structure; Kirchhoff star–mesh identity | `response_rank.log` |
| `enumerate_pieces.py` | two-terminal decompositions of the D=3, n=4 sets | `enumerate_pieces.log`, `pieces_d3n4.json` |
| `routeB_bubble.py` | dressed tight SP bubble + 3-loop core $\Gamma'$ | `routeB_bubble.log`, `routeB_bubble2.log` |
| `routeA_bridge.py` | dressed bare bridge (D=2), K4 periods | `routeA_bridge.log`, `routeA_bridge_n4.log` |
| `routeA_general.py` | dressed tight two-terminal S, any H, $\hat S$ periods | `routeA_general_n3.log`, `routeA_general_n4{,b,c,d}.log` |
| `check_k3.py` | three-terminal stars, composite stars, T=0,1 | `check_k3.log` |
| `check_k5_elim.py` | degree-4 vertex elimination in the 32 K5-cored graphs | `check_k5_elim.log` |

## 1. Kernel-level locality theorem

### 1.1 Data
$\Gamma=(V,I)$, $D=H+T$, $S\subset I$, $O=I\setminus S$, $\partial S=V(S)\cap V(O)$ ($k$ terminals), $S^\circ=V(S)\setminus\partial S$ ($m$ vertices, root $\notin S^\circ$), $N_S=|S|$, $n_S=b_1(S)$.
$$L_S=B_S^{\mathsf T}t^{-1}B_S=\begin{pmatrix}L_{\partial\partial}&L_{\partial\circ}\\L_{\circ\partial}&L_{\circ\circ}\end{pmatrix},\quad L_r=L_{\partial\partial}-L_{\partial\circ}L_{\circ\circ}^{-1}L_{\circ\partial},\quad E=-L_{\circ\circ}^{-1}L_{\circ\partial},\quad h=E^{\mathsf T},\quad Z=L_r^{+}$$
$$\lambda'_{\partial}=\lambda_\partial+h\lambda_\circ,\qquad \textstyle\sum_ah_{av}=1$$
$$y_S=F(j)+y^D,\qquad y^D=t^{-1}B_S\begin{pmatrix}0\\L_{\circ\circ}^{-1}\lambda_\circ\end{pmatrix},\qquad j_a=\mathrm{div}_a(y_S)+(h\lambda_\circ)_a$$
$F(d)$ = electrical flow of $S$ with terminal divergence $d$ ($\sum_ad_a=0$); $j$ = injections of $\Gamma'$ into the response element [drv].

### 1.2 Adapted cycle basis, block structure
Kruskal tree with $S$ first: cycles $C^S$ of $S$ ($n_S$); outer cycles $c_r$ ($n-n_S$), $d^{(r)}$ = terminal divergence of $c_r|_S$;
$$c'_r=c_r|_O+F(d^{(r)}),\qquad c'_r-c_r\in H_1(S;\mathbb R)\ \Rightarrow\ \text{unimodular change}$$
$$\sum_{e\in S}C^S_{re}F(d)_et_e=0\ (\text{KVL of }F)\ \Rightarrow\ M=\begin{pmatrix}M_S&0\\0&M_{\rm out}\end{pmatrix},\quad M_{\rm out}=M_O+D^{\mathsf T}ZD,\quad U_\Gamma=U_S\det M_{\rm out}$$
[drv; chk `check_locality.log` (all form checks use it)]

### 1.3 Two identities
$$\partial_{t_e}\big(d^{\mathsf T}Zj\big)=F(d)_eF(j)_e\ \Rightarrow\ \sum_{e\in S}F(d)_eF(j)_e\,dt_e=d^{\mathsf T}dZ\,j$$
$$\sum_{e\in S}F(d)_ey^D_et_e=\big(B_S^{\mathsf T}F(d)\big)\cdot\phi^D=0\ \Rightarrow\ \sum_{e\in S}F(d)_ey^D_e\,dt_e=-\sum_et_eF(d)_e\,dy^D_e=-\phi^{F}\cdot d(\mathrm{div}\,y^D)=d^{\mathsf T}Z\,dh\,\lambda_\circ$$
[drv: Rayleigh (envelope) for the first; $\phi^D|_\partial=0$, $\mathrm{div}_\partial y^D=-h\lambda_\circ$ for the second]
$$\psi^{{\rm out},i}_r:=\psi^i_{c'_r}=\sum_{e\in O}C_{re}y^i_e\,dt_e+d^{(r)\mathsf T}\big(dZ\,j^i+Z\,dh\,\lambda^i_\circ\big)$$
= pull-back of a 1-form on $(t_O,Z,h)$ [drv; chk: $dZ,dh$ computed from $\partial_{t_e}L_S$ independently of the flows, `factor_T0` mode R].

### 1.4 Theorem, $T=0$
$$\Omega_\Gamma=\varepsilon\;\Omega_S\wedge\pi^*\Omega_{\rm out},\qquad \Omega_S=U_S^{-H}\prod_{i,r}\psi^i_r(S),\quad \Omega_{\rm out}=(\det M_{\rm out})^{-H}\prod_{i,r}\psi^{{\rm out},i}_r,\quad \varepsilon=\pm1$$
$\pi:t\mapsto(t_O,Z(t_S),h(t_S))$; $\Omega_S$ = $T=0$ form of $S$ alone with sources $\lambda_\circ$ at $S^\circ$, $j-h\lambda_\circ$ at $\partial S$ [drv; chk `check_locality.log`: $H=2,3$, $k=2..4$, tight and non-tight $S$, $\varepsilon=\pm1$, 0 failures]
$$\boxed{K_\Gamma(\lambda;z)=\varepsilon\int_{\mathcal R_S\times\mathbb R_{>0}^{O}}\rho_S\big(Z,h;\,j\to-\partial_{z'},\lambda_\circ\big)\wedge\Omega_{\rm out}\big(t_O,Z,h;\lambda_O+h\lambda_\circ\big)\,e^{-\sum_{e\in O}y_e\cdot z_e},\qquad \rho_S=\pi_*\Big(\Omega_S\,e^{-\sum_{e\in S}y_e\cdot z_e}\Big)}$$
$$\sum_{e\in S}y_e\cdot z_e=\sum_aj_a\cdot z'_a+\sum_{v\in S^\circ}\lambda_v\cdot\Big(\sum_eg_{v,e}z_e\Big),\qquad z'_a=\sum_{e\in S}F(e_a-e_{a_0})_e\,z_e$$
($g_{v,\cdot}$ = Dirichlet flow of a unit source at $v$; $z'$ = mapped displacements on the response element; $j$-monomials = insertions on it) [drv]
$$\deg\rho_S=n_SH-f_S,\qquad f_S=N_S-\dim\mathcal R_S,\qquad \mathcal R_S=\mathrm{im}\,(t_S\mapsto(L_r,h))$$
$\rho_S$ polynomial in $(j,\lambda_\circ)$ of degree $n_SH$ [drv].
Response element of $\Gamma'=O+(\cdot)$: $k=2$ an edge of time $R_{\rm eff}$; $k=3$ a 2-port (triangle without its circulation, §1.7); $k\ge4$ a $(k-1)$-port $Z$ restricted to $\mathcal R_S$ [drv].

### 1.5 Theorem, $T=1$
$$\phi\sim N(0,\tfrac12\Lambda^{-1}),\quad w=t^{-1}B\phi,\quad \psi^\tau_r=\sum_eC_{re}w_e\,dt_e,\quad \Omega_\Gamma=U^{-H-\frac12}\Big\langle\prod_r\psi^\tau_r\Big\rangle\wedge\prod_{i,r}\psi^i_r$$
(equal to the brute covariance $\tfrac12(\delta/t-C^{\mathsf T}M^{-1}C)$: $t^{-1}B\Lambda^{-1}B^{\mathsf T}t^{-1}+C^{\mathsf T}M^{-1}C=t^{-1}$ [drv])
Markov: $\phi_\circ=E\phi_\partial+\phi^D$, $\phi^D\sim N(0,\tfrac12L_{\circ\circ}^{-1})$ independent of $\phi_{\Gamma'}\sim N(0,\tfrac12\Lambda'^{-1})$, $\Lambda'=\Lambda_O+L_r$ on $V(O)\cup\partial S$.
$$\psi^{\tau,{\rm out}}_r=\sum_{e\in O}C_{re}w_e\,dt_e+d^{(r)\mathsf T}\big(dZ\,j^\tau+Z\,dh\,\xi\big),\qquad j^\tau=L_r\phi_\partial,\qquad \xi=L_{\circ\circ}\phi^D\sim N(0,\tfrac12L_{\circ\circ})$$
$$\Omega_\Gamma=\varepsilon\,(U_S\det M_{\rm out})^{-H-\frac12}\Big\langle\prod_r\psi^\tau_r(S)\prod_r\psi^{\tau,{\rm out}}_r\Big\rangle_{\phi_{\Gamma'},\phi^D}\wedge\prod_{i,r}\psi^i_r(S)\wedge\prod_{i,r}\psi^{{\rm out},i}_r$$
[drv: second identity of §1.3 with $\lambda_\circ\to\xi$; chk `check_locality_T1.log`: $(H,T)=(1,1),(2,1)$, $n=2,4$, $k=2..5$, modes R; odd $n$: $D\equiv0$ ($\rho_{\rm odd}=0$, src P/explore/circuit-general/ht-d2/REPORT.md §1.9)]
$\xi$ = Dirichlet Gaussian source at $S^\circ$; its Wick pairs $\langle\xi\xi\rangle=\tfrac12L_{\circ\circ}$ between outer forms = contact term [drv; §3].

### 1.6 $k=2$: dressed edge
$\Gamma'=O+e_*$ ($e_*:p\to q$, $t'_{e_*}=R:=R_{\rm eff}(t_S)=U_{S/pq}/U_S$), $Z=R$, $h_v$ = harmonic measure of $p$ at $v$, $s=\sum_{v\in S^\circ}h_v\lambda_v$, $J=y'_{e_*}$, $\sigma=J-s$ = current from $p$ into $S$.
$$\psi^{{\rm out},i}_r=\psi'^i_r\big|_{t'_{e_*}=R}+C'_{r e_*}R\,ds^i$$
Tightness and complement:
$$\delta_S:=N_S-1-n_SD\ge0,\qquad \delta_S+\delta_O=D-1\quad(\Gamma=S\cup_{\{p,q\}}O)$$ [drv: $N-1=nD$, $n=n_S+n_O+1$; $\ge0$ = admissibility]
$\delta_S=0$: $\Omega_S$ is top degree on the fiber $\{R_{\rm eff}=R\}$, the $ds$ terms drop:
$$\boxed{K_\Gamma(\lambda;z)=\kappa\int_{\mathbb{RP}^{N_S-1}}D_S\big(t_S;\sigma,\lambda_\circ\big)\Big|_{\sigma\to J-s(t_S)}e^{\lambda_\circ\cdot X_\circ(t_S;z)}\;K_{\Gamma'}\big[J\to-\partial_{z'_{e_*}}\big]\big(\lambda_O,\lambda_p+s,\lambda_q+\textstyle\sum\lambda_\circ-s;\ z'(t_S)\big)}$$
$$z'_e=z_e\ (e\in O),\qquad z'_{e_*}=\sum_{e\in S}a_e(t_S)z_e,\qquad X_v=-\sum_{e\in S}g_{v,e}(t_S)z_e,\qquad \kappa=\frac{D_\Gamma|J_\Phi|}{D_SD_{\Gamma'}}=\pm1$$
$D_S$ = `htk.brute.density` of $S$ alone ($q$ root, source $\sigma$ at $p$, $\lambda_\circ$ inside); $\Phi:(t_S,t')\mapsto t$, $t_S\to t_S\,t'_{e_*}/R_{\rm eff}(t_S)$; $J_\Phi$ its Jacobian on the charts $t_{S,{\rm last}}=1$, $t'_{\rm last}=1$, $\partial_{t_e}R_{\rm eff}=a_e^2$ [drv; chk `check_locality.log`: pointwise $D_\Gamma J_\Phi=\varepsilon D_SD_{\Gamma'}$, $\varepsilon=1$ ($H=2,3$); currents $y_O=y'_O$, $y_S$ = currents of $S$ alone with $(\lambda_\circ,\sigma)$; mode P fails for $\delta_S>0$: 7/7 $(2,0)$, 9/10 $(3,0)$]
Moments as periods of $\hat S=S+e_*$, $U_{\hat S}=t_*U_S+U_{S/pq}$:
$$D_S=\sum_\beta\frac{P_\beta(t)}{U_S^{\,n_S+H}}\sigma^\beta,\qquad h_v=\frac{\mathcal H_v(t)}{U_{S/pq}},\qquad \int_{\mathbb{RP}^{N_S-1}}\frac{P(t)}{U_S^{A}U_{S/pq}^{B}}=\frac{\Gamma(A+B)}{\Gamma(A)\Gamma(B)}\,I_{A+B}\big(\hat S;\alpha_S,\alpha_*=A\big)$$
$\frac{P}{U_S^AU_{S/pq}^B}=\frac{P}{U_S^{A+B}}R_{\rm eff}^{-B}$: dressed edge = power $-B$ of the effective resistance = weight $A$ on the closing edge [drv: $\int_0^\infty dR\,R^{A-1}(RU+V)^{-A-B}=B(A,B)U^{-A}V^{-B}$; chk routes A, §2]
Source-free tight $S$, $H\ge2$, $n_S\ge1$: $y_S=aJ$, $\psi^i_r(S)=J^i\alpha_r$, $\prod_i\alpha_r\wedge\alpha_r=0\Rightarrow D_S\equiv0\Rightarrow D_\Gamma\equiv0$ [drv; chk `check_sourcefree.log`: $H=2,3$, 3 graphs, 0 nonzero without internal sources, 9/9 nonzero with]
$T=1$, $\delta_S=0$ ($N_S-1=n_SD$), point with $t'_{e_*}=R_{\rm eff}(t_S)$:
$$D_\Gamma J_\Phi=\varepsilon\sum_{j=0}^{n_S}\rho_S^{(j)}(t_S)\,D_{\Gamma'}^{[w_{e_*}^j]}(t'),\qquad \rho_S^{(j)}=[\omega^j]\;U_S^{-H-\frac12}\Big\langle\prod_r\sum_{e\in S}C^S_{re}\big(a_e\omega+w^D_e\big)dt_e\Big\rangle_{\phi^D}\wedge\prod_{i,r}\psi^i_r(S)$$
$$D^{[F]}_{\Gamma'}=U'^{-H-\frac12}\Big\langle F\prod_r\psi'^\tau_r\Big\rangle\wedge\prod\psi'^i_r\Big/\iota_E{\rm vol},\qquad \rho_S^{(j)}(ct)=c^{j/2}\rho_S^{(j)}(t)\ \Rightarrow\ \eta_*=\sqrt R\,w_{e_*}$$
$\rho_S^{(j)}$ from $S$ alone (Dirichlet field, $p,q$ grounded); $D'^{[w^j]}$ from $\Gamma'$ alone; radicand $U_\Gamma=U_SU'$ [drv; chk `check_locality_T1.log`: $(1,1)$, Wheatstone bridge + $O$, $n=4$, 4 points, ratio 1; $\rho^{(1)}=0$] — $=$ §1.5 of atlas/REPORT.md ($k=2$ chain moves, $\eta_*$) for $S$ series–parallel [src P/explore/circuit-general/atlas/REPORT.md §1.5].

### 1.7 $k=3$
$\mathcal R_S$: $L_r$ = triangle conductances $c_{ab}=-L_{r,ab}$, $Z$ = 2-port impedance (third terminal grounded), $h\in(\Delta^2)^m$.
$$\det M_{\rm out}\cdot P_\Delta=U_{O+\triangle}\big(\rho_{ab}=1/c_{ab}\big),\qquad P_\Delta=\rho_{12}+\rho_{13}+\rho_{23}$$
[drv: Schur complement of the triangle circulation; chk `check_locality.log`, 17/17]
$\Omega_{\rm out}$ = form of $O$ + 2-port: the triangle circulation is not a loop current; $\Omega_{O+\triangle}$ has $H$ more degrees than $\dim(t_O,c)$ ($N'-1-n'H=3-H-N_S+n_SH$) [drv].
Single star ($m=1$, $n_S=0$): $\mathcal R_S$ = cone, $f_S=0$, $\rho_S=1$, $h_a=\rho_{bc}/P_\Delta$ (rational transfer, Y1); $\Omega_{\rm out}$ contains $d^{\mathsf T}Z\,dh(c)\,\lambda_v$ = Jacobian of the rational transfer [drv; chk `check_k3.log`].
LLC configuration ($ux,wx\in O$): merge $c_{ux},c_{wx}$ with $O$-edges; push-forward of $\Omega_{\rm out}$ along the merge fractions $(r_1,r_2)$ = $J_{\rm single}=r_1r_2\det(y_{ux},y_{wx})$ [src atlas/REPORT.md §1.4, ht-d2/REPORT.md §1.6].

### 1.8 Response space
[chk `response_rank.log`, exact Jacobian ranks at a random rational point]

| $S$ | $k$ | $m$ | $N_S$ | $n_S$ | $\dim L_r$-image | $\dim(L_r,h)$-image | cone $\binom k2$ |
|---|---|---|---|---|---|---|---|
| star | 3 | 1 | 3 | 0 | 3 | 3 | 3 |
| star, one leg subdivided | 3 | 2 | 4 | 0 | 3 | 4 | 3 |
| star, three legs subdivided | 3 | 4 | 6 | 0 | 3 | 6 | 3 |
| paired stars $K_{2,3}$ | 3 | 2 | 6 | 2 | 3 | 6 | 3 |
| $K_4$, one internal vertex | 3 | 1 | 6 | 3 | 3 | 5 | 3 |
| star | 4 | 1 | 4 | 0 | 4 | 4 | 6 |
| star | 5 | 1 | 5 | 0 | 5 | 5 | 10 |
| two adjacent centres | 4 | 2 | 5 | 0 | 5 | 5 | 6 |
| two adjacent centres + edge | 4 | 2 | 6 | 1 | 5 | 6 | 6 |
| $K_4$ on terminals | 4 | 0 | 6 | 3 | 6 | 6 | 6 |
| $K_{2,4}$ | 4 | 2 | 8 | 3 | 6 | 8 | 6 |
| $K_{3,4}$ | 4 | 3 | 12 | 6 | 6 | 12 | 6 |
| bridge $K_4-pq$ | 2 | 2 | 5 | 2 | 1 | 3 | 1 |
| bridge, two subdivisions | 2 | 4 | 7 | 2 | 1 | 5 | 1 |

- $k=2,3$: $L_r$-image = full cone for every connected $S$ with all $c_{ab}>0$ [drv: $k=2$ scaling; $k=3$ Y–Δ bijective; chk table].
- Degree-$k$ vertex: $c_{ab}=g_ag_b/G$, $h_a=g_a/G$; image = toric subvariety $c_{ab}c_{cd}=c_{ac}c_{bd}$ of dimension $k$, codimension $k(k-3)/2$ ($k=4$: 2, $k=5$: 5) [drv; chk `response_rank.log`, `check_k5_elim.log`].
- Full cone at $k=4$ needs $\ge6$ independent response directions (e.g. $K_4$, $K_{2,4}$); two adjacent centres give codimension 1 [chk].

### 1.9 Checks (exact Fractions, all against `htk.brute`)

| check | cases | result | log |
|---|---|---|---|
| form, mode R, $T=0$ | $(2,0)$ $n=3,4$; $(3,0)$ $n=3$; $k=2,3,4$; random $S$ | all $\pm1$ | `check_locality.log` |
| form, mode P, $T=0$ | $k=2$ | tight: 1; non-tight: fails | `check_locality.log` |
| chart, $k=2$, $T=0$ | D2 bridge, D3 bridges | $\kappa=1$ | `check_locality.log` |
| form, mode R, $T=1$ | $(1,1)$ $n=2,4$; $(2,1)$ $n=2,4$; $k=2..5$ | all $\pm1$ | `check_locality_T1.log` |
| chart, $k=2$, $T=1$ | $(1,1)$ $n=4$ bridge | 1 | `check_locality_T1.log` |
| $k=3$ stars, composite stars | §3 | §3 | `check_k3.log` |
| degree-4 elimination | 8 K5-cored graphs | §4 | `check_k5_elim.log` |

## 2. Dressed two-terminal elements

### 2.1 Pieces in the $D=3$, $n=4$ sets
[chk `enumerate_pieces.log`, `pieces_d3n4.json`; sources as atlas/e3_open_kernels.py]
Non-SP graphs: 132 (bulk3 47, twistor4 85). Two-vertex cuts $\{p,q\}$, $S$ = edges at a union of components:

| piece (tight, $\delta_S=0$) | $N_S$ | $n_S$ | internal sources | $\Gamma'$ | bulk3 | twistor4 |
|---|---|---|---|---|---|---|
| bridge $K_4-pq$ + 2 subdivisions | 7 | 2 | 2 | 2 loops, SP | 8 | 15 |
| core minus a length-3 chain | 10 | 3 | 2, 3 or 4 | box | 24 | 56 |
| bubble (chains 2+2 or 1+3) | 4 | 1 | 2 | 3 loops, non-SP | 4 | 19 |
| none | — | — | — | — | 17 | 15 |

Graphs with a bridge or a bubble: 37; only core-minus-chain: 63; no tight piece: 32 [chk].

$\delta_S+\delta_O=D-1$: $D=2$: exactly one side of a two-vertex cut is tight; $D=3$: at most one [drv].
Bridge types (chain lengths $ab;\,pa,pb,qa,qb$): $(3;1,1,1,1)$ 5 graphs, $(2;2,1,1,1)$ 8, $(1;$ two sides of length 2$)$ 10 [chk].

### 2.2 Dressed series–parallel bubble, $\Gamma'$ = 3-loop core (route B)
$g\lambda_u=e_1$, $g\lambda_w=e_2$; $D_S=\sigma^3c(t_S)$ (rows of component 3 contain only $\sigma$); $s^3=0$:
$$K_\Gamma(\lambda)=\det(g)^{-n}\,\kappa\sum_{a,b\le3}m[a,b]\,[\alpha^a\beta^b]\,K_{\Gamma'}\big[y^3_{e_*}\big]\big(\lambda'_p=\lambda_p+\alpha e_1+\beta e_2\big),\qquad m[a,b]=\int c\,h_u^ah_w^b$$
$K_{\Gamma'}[y^3_{e_*}]$: degree $\le3$ in $\alpha$ and in $\beta$ (3 loops); 16 exact evaluations (`cores.kernel`, insertion $y^3_{e_*}$); $m$: exact (ratio coordinates, polynomial) [drv; chk `routeB_bubble.log`]

| graphs | equal to atlas logs | $\kappa$ | $\Gamma'$ data | time per graph | whole graph (`fastker`) |
|---|---|---|---|---|---|
| twistor4 #44–#56 | 13/13 | $-1$ | 16 exact kernels, $16\times8272=132\,352$ samples | 126–143 s | 4 384 432 points/prime, 120–167 s |
| twistor4 #58–#62 | 5/5 | $-1$ | same | 124–133 s | same |
| twistor4 #57, bulk3 #7–#10 ($K=0$) | via route A (§2.3) | — | — | — | — |

[chk `routeB_bubble.log` (#44–#56; ends at #57 with $D\equiv0$ at the $\kappa$ point, traceback in `routeB_bubble.err`), `routeB_bubble2.log`]; 18/18 nonzero bubble graphs exact (15 with $\zeta_3$); bubble moments $m[a,b]$: $<1$ s; sample count reduced 33× per prime (×2–4 primes), wall time equal (exact Python vs vectorized modular) [chk].

### 2.3 Dressed non-SP bridge (route A)
$$K_\Gamma=\kappa\sum_\beta\sum_\gamma\Big[\textstyle\prod_v h_v^{\gamma_v}\Big]\text{-coefficient of }K_{\Gamma'}\big[(J-s)^\beta\big](h)\ \times\ \frac{\Gamma(A+B)}{\Gamma(A)\Gamma(B)}\sum_{\rm mon}I_{A+B}\big(\hat S;\alpha,\alpha_*=A\big)$$
$A=n_S+H$ (reduced by exact division by $U_S$), $B=|\gamma|$ (reduced by $U_{S/pq}$); $P_\beta$, $\mathcal H_v$, $U_S$, $U_{S/pq}$ by exact interpolation from `brute` on $S$ alone (verified at independent points); $K_{\Gamma'}[J^{b}](h)$ by `htk.kernel_exact` with symbolic $h_v$; $I$ by `cores.period` on $\hat S$ (series reduction to $K_4=W_3$); no divergent period met (`cores.period` raises on divergent input) [drv; chk `routeA_bridge*.log`, `routeA_general_*.log`]

| case | graphs | equal to reference | moment monomials | $\hat S$ periods | time | reference |
|---|---|---|---|---|---|---|
| $D=2$, $n=3$, bare bridge, $\Gamma'$ = triangle | 1 graph × 3 source sets | 3/3 | 489 | 224 ($K_4$) | 37–96 s (sympy) | `cores.kernel` 0.3 s |
| $D=2$, $n=4$, bare bridge, $\Gamma'$ 2 loops | 1 (only one in the census) | 1/1, $\frac{132197}{43200}-\frac{1251}{200}\zeta_3$ | 1535 | 505 | 127 s | `fastker` 0.8 s |
| $D=3$, $n=3$, bridge + 2 subdivisions, $\Gamma'$ = box (sources at all vertices) | K4[3,2,1,1,2,1] ×2, K4[3,1,2,1,1,2] ×2 | 4/4 | 0.9–1.4·10^5 | 2438–4951 | 27 s | `fastker` 0.1 s |
| $D=3$, $n=4$, bridge + 2 subdivisions, $\Gamma'$ 2 loops SP | all 23 (bulk3 #3–#10, twistor4 #28, 29, 36–38, 47, 48, 55–57, 63–67) | 23/23 (17 nonzero, 13 with $\zeta_3$, 6 zero) | 63–49182 | ≤ 12724 (cumulative cache) | 15–18 s | atlas logs, 77–170 s |

[chk `routeA_bridge.log`, `routeA_bridge_n4.log` (symbolic $D_S$, sympy), `routeA_general_n3.log`, `routeA_general_n4.log`, `routeA_general_n4b.log` (ends at #57: $D\equiv0$ at the $\kappa$ point), `routeA_general_n4c.log`, `routeA_general_n4d.log` (interpolated $D_S$); $n=4$ sources = atlas logs (bulk3: bivalent vertices; twistor4: bivalent vertices and $O$)]
- $D=3$, $n=4$: $\kappa=\pm1$ per graph; root internal to $S$ (bulk3): relabel, $K\to\mathrm{sgn}(\pi)^3K$ (CONVENTIONS §2.2, checked with `brute` where $D\ne0$).
- $D=3$, $n=4$ dressed-edge data: $\sum_\beta|P_\beta|\le12$ monomials ($\deg P_\beta=3$ in 7 variables), $U_S$, $U_{S/pq}$, $\mathcal H_v$ (degree 2, 3, 3): 15 s per (bridge, $\lambda_\circ$).
- Whole-graph interpolation (`fastker`, $N=13$): $|\{e\in[0,3]^{12}:\sum e\le15\}|=4\,384\,432$ points per prime, 2–4 primes, 1.5 GB; route A: $\le5\cdot10^4$ monomials, $\le1.3\cdot10^4$ $K_4$-type periods [chk count].
- $n=3$: route A slower than whole-graph ($|\{e\in[0,2]^9:\sum e\le8\}|=8272$ points); the gain is at $n=4$.

## 3. Three-terminal elements, $D=2$: $(H,T)=(2,0),(1,1)$
[chk `check_k3.log`: all census graphs $D=2$, $n=4$; $S$ = edges at $W$]

| $S$ | interior | $(2,0)$ mode R | $(2,0)$ mode P | $(1,1)$ mode R | $(1,1)$ mode P |
|---|---|---|---|---|---|
| single star | sourced | 29 ±1, 1 zero | 30 FAIL | 27 ±1, 3 zero | 29 FAIL |
| single star | unsourced | 29 ±1, 1 zero | 30 ±1 | 28 ±1, 2 zero | 30 FAIL |
| star + subdivided leg | sourced | 4 ±1 | 4 FAIL | 4 ±1 | 4 FAIL |
| star + subdivided leg | unsourced | 4 zero | 4 zero | 4 zero | 4 zero |
| paired stars $K_{2,3}$ | sourced | 4 ±1 | 4 FAIL | 4 ±1 | 4 FAIL |
| paired stars $K_{2,3}$ | unsourced | 4 ±1 | 4 ±1 | 4 ±1 | 4 FAIL |

Mode P = outer forms without $d^{\mathsf T}Z\,dh\,\lambda_\circ$ and without $d^{\mathsf T}Z\,dh\,\xi$.
- $T=0$: P exact $\iff$ $dh\,\lambda_\circ$ absent (unsourced interior) = LLC Y0, paired moves; sourced interior: the rational-transfer Jacobian $d^{\mathsf T}Z\,dh\,\lambda_\circ$ is required (Y1) [chk].
- $T=1$: P fails also for unsourced interior; R exact. The missing term $d^{\mathsf T}Z\,dh\,\xi$, $\xi=L_{\circ\circ}\phi^D$, $\langle\xi\xi\rangle=\tfrac12L_{\circ\circ}$ ($=\tfrac G2$ for a star) is the contact term $\frac1{2G}\alpha\wedge\beta$ of P/explore/circuit-general/ht-d2/REPORT.md §1.5 [chk; identification drv].
- Composite stars: $\dim(L_r,h)$-image $=N_S$ for star + legs and $K_{2,3}$ (§1.8) $\Rightarrow f_S=0$, $\rho_S=\Omega_S$ itself, of degree $n_SH$ on $(c,h)$; base $c$ only ($K_{2,3}$): $f=3$, $\deg\rho_S=n_SH-3=1$ ($H=2$) [drv].
- $(1,1)$ non-SP forests (cores/REPORT.md §3): not recomputed by a 3-terminal reduction [open]; obstruction: $h_a=\rho_{bc}/P_\Delta$ and $\langle\xi\xi\rangle$ are not Laurent monomials in a product-of-simplices chart (as $Q$ in ht-d2 §1.6) [drv].

## 4. Degree $\ge4$ vertices; $K_5$
Star at $v$, conductances $g_a$ ($a=1..k$), $G=\sum g_a$:
$$L_r=\mathrm{diag}(g)-\frac{gg^{\mathsf T}}G,\qquad h_a=\frac{g_a}G,\qquad f_S=0,\quad \rho_S=1,\qquad \Omega_\Gamma=\varepsilon\,\Omega_{\rm out}$$
[chk `check_k5_elim.log`: 8 of the 32 $D=2$, $n=6$ K5-cored graphs, $v$ of degree 4: mode R 16/16 exact (sourced and unsourced $v$); mode P exact iff $\lambda_v=0$; $L_r$, $h$ formulas 8/8]
Kirchhoff polynomial (conductances $x$, $K_G=\sum_T\prod_{e\in T}x_e$):
$$K_{K_5}(x,g)=G\,K_{K_4}\Big(x_{ab}+\frac{g_ag_b}G\Big)$$ [chk `response_rank.log`, symbolic]
Period ($\beta=q-\alpha$, $I_q(G;\alpha)=\int\prod x_e^{\beta_e-1}K_G^{-q}\Omega$):
$$I_q(K_5;\alpha)=\int_{\Delta^3}d\hat g\prod_a\hat g_a^{\beta_a-1}\;\Phi\big(\hat g\hat g^{\mathsf T}\big),\qquad \Phi(c)=\int_{\mathbb R^6_{>0}}\prod_{ab}x_{ab}^{\beta_{ab}-1}K_{K_4}(x+c)^{-q}\,d^6x$$ [drv: gauge $G=1$]
Fiber over the merged $K_4$ conductances $x'=x+gg^{\mathsf T}/G$: dimension 4; merge fractions $r_{ab}=\frac{g_ag_b}{Gx'_{ab}}\in(0,1)^6$ subject to
$$r_{ab}r_{cd}x'_{ab}x'_{cd}=r_{ac}r_{bd}x'_{ac}x'_{bd}=r_{ad}r_{bc}x'_{ad}x'_{bc}$$
(2 equations; $k=3$: none, LLC $(r_1,r_2)\in(0,1)^2$) $\Rightarrow$ no product-of-simplices chart, no Dirichlet scheme [drv].
Star in Schwinger form:
$$\int\frac{d^{2q}X_v}{\pi^q}\prod_a|X_v-X_a|^{-2\beta_a}=\frac{\Gamma(\sum\beta-q)}{\prod\Gamma(\beta_a)}\int_{\mathbb{RP}^{3}}\prod_ag_a^{\beta_a-1}\,G^{-\beta_{v\infty}}\Big(\sum_{a<b}g_ag_b|X_{ab}|^2\Big)^{q-\sum\beta},\qquad \beta_{v\infty}=2q-\sum_a\beta_a$$ [drv]
- variables: 6 distances, 5 after scaling = 2 complex + 1 angle of graphical/REPORT.md §7 [drv; src].
- $\beta_{v\infty}=0$ (T1): $G$ drops out, conformal 4-point function, one cross-ratio after fixing three points = the graphical-function reach criterion [drv; src graphical/REPORT.md §7].
- $\beta_{v\infty}\ne0$ (generic $K_5$): Schwinger elimination = reparametrization ($f_S=0$); representation of dimension 3 (star) + 6 ($K_4$); no lower-dimensional exact formula; nothing new to test on the generic orbits of `graphical/k5_reach_fraction.log` [drv].

## 5. Not done; failed
- $T=1$ dressed elements integrated exactly (moments of $\rho_S^{(j)}$ with $\eta_*^j$ insertions in $\Gamma'$; half-integer $\hat S$ levels): not implemented; $(1,1)$ forests of cores/REPORT.md §3 not recomputed [open].
- $k=3$ composite stars at $D=2$: $f_S=0$ relative to $(c,h)$ for every composite star met (§1.8), nothing to precompute; exact R-class values not recomputed through $k=3$ [open].
- "Core minus a length-3 chain" pieces ($N_S=10$, $n_S=3$, $\Gamma'$ = box; present in 80 of the 132 graphs, the only tight piece in 63): route A needs $P_\beta$ of degree $(n_S+H)n_S-N_S=8$ in 10 variables ($\binom{17}{9}=24310$ coefficients, $64$ $\sigma$-points each) and $\hat S$ levels $q\ge8$ (prism/K33 tables only to $q=7$); not run [drv].
- First route B run: overall sign wrong (orientation of $\Phi$ omitted, $\kappa=\varepsilon\,\mathrm{sgn}J_\Phi$); fixed. First route A run at $D=3$: period cache keyed without the edge list (wrong value for the second bridge type); fixed, rerun equal.
- Generic $K_5$: no exact Schwinger-space formula (§4).
