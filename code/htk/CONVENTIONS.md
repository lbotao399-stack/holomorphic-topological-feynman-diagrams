# htk conventions

Sources: S1 = `sources/LLC_main.tex` (version 20261120-231110; Section "Positive Schwinger times and the operatope", Section "Fan kernels"), S2 = `Desktop/6d-holomorphic-twist-note` (sec5–sec7), S3 = `Twisted_11d_sugra_loops_2026-11-18/main.tex` (Sections 4–6, Appendix A), S4 = `Desktop/check_cd`, S6 = `code/enumerate_graphs.py`. Equation labels are those of the TeX sources.

Every statement below marked [t: name] is checked in exact arithmetic by the test `name` of `code/tests/test_known.py` (k) or `code/tests/test_engine.py` (e).

## 1. Graph

`Graph(nv, edges)`: vertices $0,\dots,nv-1$; root $=nv-1$ (position $0$, source $-\sum_{v<\mathrm{root}}\lambda_v$); `edges[e] = (e(0), e(1))`, ordered, oriented; $N=|I|$, $n=N-nv+1$.

$$B_{ev}=\delta_{v,e(0)}-\delta_{v,e(1)}\quad(v<\mathrm{root}),\qquad t\,y=B\varphi,\quad B^{\mathsf T}y=\lambda,\quad \varphi=\Lambda^{-1}\lambda,\quad\Lambda=B^{\mathsf T}t^{-1}B .$$

Spanning tree $T$ of `Graph.spanning_tree()`: Kruskal in edge order. Cotree $c_1<\dots<c_n$. Fundamental cycles $C_T$ ($n\times N$): $(C_T)_{r,c_r}=+1$, cycle $r$ runs along $c_r$ and back through $T$, $(C_T)_{re}=\pm1$ for $e$ along/against.

$$\sigma_T=\varepsilon_T\det B_T=\pm1,$$

$\varepsilon_T$ the sign of the permutation (tree edges, cotree edges), each in increasing order; $B_T$ the rows of $B$ on $T$ (S2 eq-kernel-def). With columns ordered $(T,\text{cotree})$, $C_T=(X,1)$, $B_{\rm cot}=-XB_T$:

$$\det\begin{pmatrix}B^{\mathsf T}\\ C_T\end{pmatrix}=\varepsilon_T\det\begin{pmatrix}B_T^{\mathsf T}&-B_T^{\mathsf T}X^{\mathsf T}\\X&1\end{pmatrix}=\sigma_T\det(1+XX^{\mathsf T})=\sigma_T\,\#\{\text{spanning trees}\}.$$

Cycle matrix of the kernel, `Graph.cycle_matrix()`: $C=C_T$ with its first row multiplied by $\sigma_T$, so that $\det\binom{B^{\mathsf T}}{C}>0$. [t: e test_tree_independence]

## 2. Reduced kernel

$$M=C\,t\,C^{\mathsf T},\qquad U=\det M,\qquad\psi^i_r=\sum_eC_{re}\,y^i_e\,dt_e,\qquad\iota_E\psi^i_r=0,\quad E=\sum_et_e\partial_{t_e},$$

$$\rho_n(M)=(\det M)^{-1/2}\,\mathrm{Pf}\bigl(-\tfrac12\,dM\,M^{-1}dM\bigr),\qquad\mathrm{Pf}(A)=\sum_{j=2}^n(-1)^jA_{1j}\,\mathrm{Pf}(A_{\hat1\hat\jmath}),$$

$$\Omega_\Gamma=U^{-H}\,\rho_n(M)^T\wedge\prod_{i=1}^H\prod_{r=1}^n\psi^i_r=D(t)\,\iota_E\,(dt_0\wedge\dots\wedge dt_{N-1}),\qquad D(ct)=c^{-N}D(t),$$

$$D(t)=\frac{(-1)^{e_0}}{t_{e_0}}\,\bigl[\Omega_\Gamma\bigr]_{dt_0\cdots\widehat{dt_{e_0}}\cdots dt_{N-1}}\ \ (\text{any }e_0,\ 0\text{-based}),\qquad D\,(a\cdot t)=\frac{\det[a;\Psi]}{U^{H}}\ \ (T=0,\ \Psi_{(i,r),e}=C_{re}y^i_e),$$

$$K_\Gamma(\lambda;z)=\int_{\mathbb{RP}^{N-1}_{>0}}e^{-\sum_ey_e\cdot z_e}\,\Omega_\Gamma=\int_{t_e>0,\ \sum t_e=1}D(t)\,e^{-\sum_ey_e\cdot z_e}\,dt_0\cdots dt_{N-2},$$

$\mathbb{RP}^{N-1}_{>0}$ oriented by $\iota_E\prod_edt_e$ (positive form). Order of factors: $\rho_n$, then $i=1..H$ (outer), $r=1..n$ (inner). Degree: $n(H+T)=N-1$. $y^i_e\,z_e^i$ summed over $i$; $-\partial/\partial z^i_e$ inserts $y^i_e$.

Topological sources $\kappa_v$ (S3 Section 6), brute-force layer only: $\bar w=$ Kirchhoff solution with sources $\kappa$,

$$D_\kappa(t)\,\iota_E\mathrm{vol}=U^{-H-\frac12}\,e^{\sum_et_e\bar w_e^2}\Bigl\langle\prod_r(\bar\psi^\tau_r+\delta\psi^\tau_r)\Bigr\rangle\wedge\prod_{i,r}\psi^i_r,\qquad\langle\delta\psi^\tau_r\,\delta\psi^\tau_{r'}\rangle=\bigl(-\tfrac12dM\,M^{-1}dM\bigr)_{rr'} .$$

`brute.Density(coeff, radicand, expo)` $=\mathrm{coeff}\cdot\mathrm{radicand}^{-1/2}\,e^{\mathrm{expo}}$; $T=0$: radicand $=1$; $T=1$: radicand $=U$, expo $=\sum_et_e\bar w_e^2$.

### 2.1 Change of cycle basis

$C'=AC$: $\psi'=A\psi$, $M'=AMA^{\mathsf T}$, $\rho_n(M')=\mathrm{sgn}(\det A)\,\rho_n(M)$, and for $|\det A|=1$

$$\Omega_\Gamma[C']=(\det A)^{H+T}\,\Omega_\Gamma[C],\qquad K_\Gamma=\sigma_T^{\,H+T}\,K_\Gamma[C_T]\ \text{ for every spanning tree }T .$$

[t: e test_tree_independence]

### 2.2 Edge reversal, edge order, vertex labels ($D=H+T$)

| operation | $K_\Gamma\to$ | reason |
|---|---|---|
| reverse $k$ edges, $z_e\to-z_e$ on them | $(-1)^{Dk}K_\Gamma$ | $C_{re}y_e$, $M$ unchanged; $\det\binom{B^{\mathsf T}}{C}\to(-1)^k$ |
| edge order permuted by $\pi$ | $\mathrm{sgn}(\pi)^{D+1}K_\Gamma$ | $\iota_E\prod dt_e\to\mathrm{sgn}\,\pi$; $\det\binom{B^{\mathsf T}}{C}\to\mathrm{sgn}\,\pi$ |
| vertices relabelled by $\pi$ (root = vertex labelled last) | $\mathrm{sgn}(\pi)^{D}K_\Gamma$ | rows of $B^{\mathsf T}$ permuted; deleting column $j$ instead of the last column of the full incidence matrix gives $(-1)^{nv-1-j}$ |

[t: e test_sign_rules]

### 2.3 Independent currents

$y_e=y^0_e+\sum_r(C_T)_{re}\ell_r$, $\ell_r=y_{c_r}$; $\sum_eC_{re}t_ey_e=0\Rightarrow\psi^i_r=-\sum_{r'}M_{rr'}d\ell^i_{r'}$:

$$\prod_{i=1}^H\prod_{r=1}^n\psi^i_r=\epsilon(n,H)\,U^H\bigwedge_{r=1}^nd^Hy_{c_r},\qquad\epsilon(n,H)=(-1)^{nH+\binom n2\binom H2},\qquad d^Hy=dy^1\wedge\dots\wedge dy^H,$$

$$K_\Gamma=\sigma_T^{\,H+T}\,\epsilon(n,H)\int e^{-\sum y\cdot z}\,\rho_n^T\wedge\bigwedge_{r}d^Hy_{c_r}.$$

`brute.jacobian_density` $=D_{dy}$, $\bigwedge_rd^Hy_{c_r}=D_{dy}\,\iota_E\mathrm{vol}$, $D_{dy}=(-1)^{N-1}\det[\partial y^i_{c_r}/\partial t_j]_{j<N-1}/t_{N-1}$. [t: e test_tree_independence]

## 3. Relation to the source notes

$$\boxed{\;I^{\rm S1}_\Gamma=(-1)^{\frac{n(n+1)}2}K_\Gamma\ (H=2,T=0),\qquad I^{\rm S2}_\Gamma=(-1)^{\frac{n(n+1)}2}K_\Gamma\ (H=3,T=0),\qquad\mathcal K^{\rm S3}_\Gamma=(-1)^{\frac{n(n+1)}2}K_\Gamma\ (H=5,T=1)\;}$$

$$I^{\rm S4}_{\Gamma,T}:=\int_{\omega_m}\bigwedge_{c\notin T}d^dy_c=\epsilon(n,d)\,\sigma_T^{\,d}\,K_\Gamma .$$

Factors printed in the sources for position-space integrals: (S1 eq-1-7) heat-edge sum $=(-8\pi^2)^{-n}I^{\rm S1}$; (S2 eq-B-norm) $=(-1/8\pi^3)^nI^{\rm S2}$; (S3 Theorem thm-master) $=(-i/64\pi^{11/2})^n\mathcal K^{\rm S3}$. Functions: `conventions.S1_kernel`, `S2_kernel`, `S3_kernel`, `S4_integral`, `raw_factor`.

In particular, `raw_factor("S2", n)` returns the factor printed in S2 eq-B-norm, $(-1/(8\pi^3))^n$. Direct evaluation with $\Omega\wedge d^3\bar x=8\,d^6x_E$ and all vertex volume forms already at the front gives $+(8\pi^3)^{-n}I^{\rm S2}$. Gathering the initially interleaved vertex volumes across the odd cubic words supplies $(-1)^n$ in the full Wick graph. [Real-volume derivation](../../explore/six-d-Q3/README.md#ordered-wick-weights) · [Marked Gaussian replay](../../explore/six-d-Q3/verify/marked_volume.py)

**S1** (eq-operatope): $I^{\rm S1}=\int e^{-\sum y\cdot z}\bigwedge_{e\notin T}(-dy_{e,+\dot+}\wedge dy_{e,+\dot-})$, orientation of (eq-UV-boundary) $=\iota_E\prod dt_e$, components $(\dot+,\dot-)=(1,2)$, $\lambda_{1,+\dot a}\lambda_{2,+}{}^{\dot a}=\det(\lambda_1,\lambda_2)$ ($\epsilon^{\dot+\dot-}=1$):
$I^{\rm S1}=(-1)^n\epsilon(n,2)K=(-1)^{n+\binom n2}K$. Edge order of the fan (Theorem b-01142, $(vu,vw)$ first): $(1,2),(1,n{+}2),(2,3),(2,n{+}2),\dots,(n,n{+}1),(n,n{+}2),(n{+}1,n{+}2)$ = `graph.fan(n)` with S1 vertex $i\mapsto i-1$.
Confirmed: $I^{\rm S1}_{\rm fan}=\det(\lambda_1,\lambda_{n+1})^n/((n!)^3(n+1))$, $n=1,2,3$ (eq-1-42), density (eq-1-41); density (eq-1-40) and currents (eq-1-38), (eq-1-39) with sources on all path vertices, $n=1..4$, as polynomial identities in $r_1..r_{2n}$; triangle expansion (eq-1-46) through order 3. [t: k test_S1_fans, test_S1_fan_general_sources, test_S1_triangle_displacements]

**S2** (eq-kernel-def): $I^{\rm S2}=\varepsilon_T\det(B_T)\int_{\omega_m}e^{-\sum y\cdot z}\bigwedge_{e\notin T}d^3y_e=\sigma_T\,\epsilon(n,3)\,K[C_T]=\epsilon(n,3)K$, $\omega_m=\iota_E\prod dt_e$.
Confirmed: $\sigma_T=-1,+1,-1$ for the trees of S2; $I_{\rm box}=\frac16\det(\lambda_1,\lambda_2,\lambda_3)$; $I_A$, $I_B$ of (eq-IA), (eq-IB) with symbolic sources; densities $\mathcal J_A,\mathcal J_B$ with sources at all vertices; displacement maps (eq-res-X). [t: k test_S2_box, test_S2_graph_A, test_S2_graph_B]

**S4**: kernel convention of `check_cd.py` (no factor $\varepsilon_T\det B_T$). Confirmed: polygon $d=2..5$, glued polygons $d=2..5$, ladder $n=1..4$: $1/6,-1/432,1/155520,-1/186624000$. [t: k test_S4_polygons, test_S4_chain_move, test_S4_ladder]

**S3** (eq-master): $\mathcal K^{\rm S3}=\varsigma_\Gamma\int e^{-\sum y\cdot z}U^{-11/2}\langle\prod_r\psi^\tau_r\rangle\wedge\prod_{i,r}\psi^i_r=\varsigma_\Gamma K_\Gamma$ ($D=6$: independent of the cycle basis), $\varsigma_\Gamma=(-1)^n(-1)^{(5n+1)5n/2}\varsigma^{\rm ext}_\Gamma$. By Section 4 with $N=6n+1$, $nv=5n+2$: $\varsigma^{\rm ext}_\Gamma=(-1)^n$ for every graph, hence

$$\varsigma_\Gamma=(-1)^{\frac{5n(5n+1)}2}=(-1)^{\frac{n(n+1)}2}.$$

Confirmed: $\varsigma^{\rm ext}=-1$ ($C_7$), $+1$ (five $\theta_{pqr}$), computed exactly from the Jacobian of (eq-extraction); $\mathcal K_7=-\det[q_1..q_6]/6!$; $\mathcal K_{\theta_{pqr}}=c_{pqr}\pi D_\theta$ with (R-c); density $-\frac12D_\theta\prod_Pw_P^{11/2-n_P}$. [t: k test_S3_extraction_signs, test_S3_heptagon, test_S3_thetas]

**Two-loop theta graphs at $T=1$**, `graph.theta(p,q,r)` (edge order $A,B,C$; $N=0$, $N'=$ root), positive cycle basis, one chain move with chains $(A,B,C)$, $H=0..5$ (16 graphs):

$$g=s\,\frac{D_\theta}2\prod_Pw_P^{H+\frac12-n_P},\qquad K_\theta(\lambda;0)=s\,\frac{D_\theta}{2}\,\frac{\prod_P\Gamma(H+\frac32-n_P)}{\Gamma(H+\frac32)\prod_P(n_P-1)!},\qquad s=(-1)^{\lceil H/2\rceil+(H+1)(q+1)} .$$

$H=0$: $K=\pi$; S3 remark after Corollary cor-odd: $\frac1{4\pi}K=\frac14$. [t: e test_theta_T1_all_H; k test_S3_H0_remark]

## 4. Position-space identity (S3 eq-extraction, general $D$)

$X^a_v$ ($a=1..D$, $v<$ root), $X^a_e=\sum_vB_{ev}X^a_v$, $\theta^a_e=dX^a_e-c_aJ^a_e\,dt_e$, $\psi^a_r=\sum_eC_{re}J^a_edt_e$. With $P$ such that $PB=1$ ($P=B_T^{-1}$ on $T$), $G=\binom PC$, $\det G=\det(A)\,\sigma_T$ for $C=AC_T$:

$$\bigwedge_e\theta^a_e=(\det G)^{-1}\bigwedge_v\bigl(dX^a_v-\dots\bigr)\wedge\bigwedge_r(-c_a\psi^a_r),$$

$$\bigwedge_e\bigl(\theta^1_e\wedge\dots\wedge\theta^D_e\bigr)=\varsigma^{\rm ext}\prod_ac_a^{\,n}\Bigl[\prod_vdX^1_v\cdots\prod_vdX^D_v\Bigr]\wedge\prod_{a}\prod_r\psi^a_r+(\text{lower degree in }dX),$$

$$\varsigma^{\rm ext}=(-1)^{\binom N2\binom D2+nD+(nv-1)\,n\binom D2}\quad\text{for }\det\tbinom{B^{\mathsf T}}{C}>0;\qquad\times\,\sigma_T^{\,D}\ \text{for }C=C_T .$$

($\binom N2\binom D2$: reordering $\bigwedge_e\bigwedge_a\to\bigwedge_a\bigwedge_e$; $nD$: the signs $-c_a$; $(nv-1)n\binom D2$: moving the $dX^a$ blocks in front of the $\psi^{a'}$, $a'<a$.) `conventions.extraction_sign` evaluates the ratio exactly at random rational points; `extraction_sign_formula` is the formula. [t: e test_extraction_sign_formula]

## 5. Series–parallel layer (`spred`)

Chain move $m$ ($m=1,2,\dots$): $k\ge1$ chains from $p$ to $q$; chain $i$ has edges $j=1..a_i$ (in the order from $p$), interior vertices $u_{i,1..a_i-1}$ of degree two with stage sources $\lambda_{i,k}$; replaced by the edge $p\to q$ with id $N+m-1$.

$$R_i=\sum_jt_{i,j},\quad x_{i,j}=\frac{t_{i,j}}{R_i},\qquad\frac1{t^{\rm red}}=\sum_i\frac1{R_i},\quad w_i=\frac{t^{\rm red}}{R_i},\qquad\sum_jx_{i,j}=1,\quad\sum_iw_i=1 .$$

$$\sigma_{i,k}=\sum_{j>k}x_{i,j},\quad\hat m_i=\sum_k\sigma_{i,k}\lambda_{i,k},\quad Q_{i,j}=\sum_{k<j}\lambda_{i,k},\qquad\lambda_p^{\rm red}=\lambda_p+\sum_i\hat m_i,\quad\lambda_q^{\rm red}=\lambda_q+\sum_{i,k}(1-\sigma_{i,k})\lambda_{i,k},$$

$$y_{i,j}=w_i\,y^{\rm red}-\hat m_i+Q_{i,j}\quad(\text{direction }p\to q;\ \text{sign }-1\text{ for an edge oriented against it}).$$

Last edge $p_f\to q_f$: $t=1$ (gauge), $y=\lambda^{\rm red}_{p_f}$. Hence $t_e=\prod(\text{ratios})^{\pm1}$ (`texp`), $y_e=\sum_v$ `resp[e][v]` $\lambda_v$ with polynomial `resp`.

Names: `w{m}_{i}`, `x{m}_{i}_{j}`. Blocks (simplices): `w{m}_1..w{m}_k` if $k\ge2$; `x{m}_{i}_1..x{m}_{i}_{a_i}` if $a_i\ge2$. In each block the first ratio is eliminated. Coordinates $\xi$ (`coords`): for $m=1,2,\dots$: `w{m}_2..`, then for $i=1..k$: `x{m}_{i}_2..`. Domain: product of the open simplices. `orientation` $=\mathrm{sgn}\det[t;\partial t/\partial\xi_1;\dots;\partial t/\partial\xi_{N-1}]$.

Displacement map: $-\sum_ey_e\cdot z_e=\sum_{v<\rm root}\lambda_v\cdot X_v$, $X_v=\sum_e$ `displacement()[v][e]` $z_e$, `displacement()[v][e]` $=-$`resp[e][v]`.

Density (`RatioDensity`): $K_\Gamma(\lambda;0)=\int_{\prod\Delta}g(\xi)\,d\xi$,

$$g(\xi)=D(t(\xi))\,\bigl|\det[t;\partial_\xi t]\bigr|=\prod_bb^{\,\alpha_b}\cdot P(\xi),$$

$b$ over all ratios (including the eliminated ones), $\alpha_b\in\frac12\mathbb Z$ (`exps`; $0$ for $T=0$), $P$ a polynomial in $\xi$ (`num`, product of `factors`).

$T=0$, generic: $g=\texttt{orientation}\cdot\sigma_T^{H}\epsilon(n,H)\,\det[\partial y^i_{c_r}/\partial\xi_j]$, rows $(r,i)$, $T$ the Kruskal tree.

$T=0$, local (plan of tight chain moves, $\sum_ia_i=H(k-1)+1$): with a reduced tree containing or not containing the merged edge, $d\ell_i-w_i\sum_jd\ell_j\equiv y^{\rm red}dw_i+\sum_{j,l}(w_i-\delta_{ij})Q_{j,l}\,dx_{j,l}$ modulo the reduced loop currents, $\ell_i=y_{i,1}$:

$$g=\pm\prod_m\det L_m,\qquad L_m=\Bigl[\;y^{{\rm red},a}\,\delta_{ii'}\ \Big|\ (w_i-\delta_{ij})\,Q^a_{j,l}\;\Bigr]_{(i,a);\ (i'\ge2),\ (j,\,l\ge2)},\quad i=2..k,\ a=1..H;$$

the sign is fixed by the generic formula at one rational point. [t: e test_sp_all_tight_graphs_T0: `num` identical for both methods]

$T=1$: $\det M$ is a Laurent monomial in the ratios; $\rho_n(M(\xi))=\sum_S\rho_S\,d\xi_S$ computed for $C_T$;

$$g=\texttt{orientation}\cdot\sigma_T^{H+1}\epsilon(n,H)\sum_{|S|=n}\mathrm{sgn}(S,S^c)\,\rho_S\,\det\bigl[\partial y^i_{c_r}/\partial\xi_j\bigr]_{j\in S^c}.$$

Integration: $\int_{\Delta^m}\prod_{j=0}^mb_j^{a_j}=\prod_j\Gamma(a_j+1)/\Gamma(m+1+\sum_ja_j)$; result: `Fraction` ($T=0$), `PiNumber` $\sum_kq_k\pi^k$ ($T=1$), sympy expression (symbolic sources).

Every call of `density` with `check>0` compares, at random rational $\xi$: `ratios_of_t(tmap(xi)) == xi`; polynomial currents $=$ Kirchhoff currents; $g(\xi)=D(t(\xi))|\det[t;\partial_\xi t]|$ including sign ($T=1$: sign and square).

Coordinates of the sources in terms of the ratios:

| source | relation |
|---|---|
| S1 fan, default plan, moves $i=1..n-1$ | $r_{2i-1}=$ `x{i}_1_2`, $r_{2i}=$ `w{i}_2`; last move: $(r_{2n-1},r_{2n})=((1-w)(1-x),(1-w)x)$, $w=$ `w{n}_2`, $x=$ `x{n}_2_2` |
| S2 chain move, chains $[(e_1,e_2,e_3),(e_0)]$ | $r_1=w_2$, $r_2=x_{1,2}+x_{1,3}$, $r_3=x_{1,3}$ |
| S2 bubble move, chains $[(pv,vq),(pv',v'q)]$ | $r_1=w_2$, $r_2=x_{1,2}$, $r_3=x_{2,2}$ |
| S3 theta, one move, chains $(A,B,C)$ | $w_P=$ `w1_P`, $x^P_j=$ `x1_P_j` |
| S4 chain $(1,d)$ | $\rho=w$ of the single edge, $\sigma_j=\sum_{l>j}x_l$ |

## 6. Oracle (`oracle`)

$S(t)=D(t)\prod_et_e$ (scale invariant); $K=\int_{\mathbb R^{N-1}}S(e^{s})\,d^{N-1}s$ ($s_{N-1}=0$); trapezoidal rule on $h\mathbb Z^{N-1}+\{0,h/2\}$, flood fill of $\{|S|>\mathrm{tol}\cdot\max|S|\}$ from $s=0$. Currents from $M\ell=-C\,t\,y^0$; $T=0$: $S=\det[1;V]/(NU^H)$, $V_{(i,r),e}=C_{re}t_ey^i_e$; $T=1$, $n=2$: $S=\det W\cdot u^{\mathsf T}\hat Av/(NU^{H+1/2})$, $W=[1;V;R]$, $u,v$ the last two columns of $W^{-1}$, $\hat A_{ef}=t_et_fA_{12}[e,f]$.

Measured relative errors with the default $(h,\mathrm{tol})$: $N-1\le4$: $\le10^{-6}$; $N-1=6$: $10^{-4}$ ($T=0$, $8\cdot10^6$ evaluations), $2\cdot10^{-4}$ ($T=1$, $2\cdot10^7$ evaluations). Error estimate returned: half the difference of the two shifted sums, plus the sub-threshold layer divided by $1-e^{-h/2}$, plus a floor $5\cdot10^{-8}|K|$ ($T=1$). [t: k test_oracle, e test_oracle_T1_slow]

## 7. Catalog (`catalog`)

Functions of S6 with $d\to D=H+T$: `enumerate_tight`, `is_sp`, `chain_reducible`, `yd_reducible`, `two_loop_thetas`, `cone_over_forest`; `tight_graphs(D, n)`, `theta_graphs(D)` return `Graph` objects (`from_networkx`: vertices in sorted order, edges sorted, oriented from the smaller to the larger label). Counts equal `results/graph_counts.txt`. [t: e test_catalog_counts]

## 8. API

```
graph      Graph(nv, edges); .ne .root .loops .incidence() .spanning_tree() .cotree(tree) .fundamental_cycles(tree)
           .sigma(tree) .cycle_matrix(tree, positive=True) .reverse_edges(which) .permute_edges(order) .relabel(new_label)
           polygon(m)  theta(p,q,r)->(G, interior)  banana(k)  fan(n)  box()  graph_A()  graph_B()
           ladder(n)->(G, cotree, a, b)  glued_polygons(d)->(G, cotree, chain, base)  from_networkx(G, root)->(Graph, labels)
brute      sources(G, lam)  currents(G, t, lam)->y[e][i]  response(G, t)->c[e][v]  potentials  laplacian
           period_matrix(G, t, C)  symanzik_U(G, t, C)  psi_forms(G, t, y, C)  topological_form(G, t, C, kappa)
           density(G, t, lam, T=0, kappa=None, C=None, tree=None, check_basic=False)->Density
           integrand(G, t, lam, T, powers, z)->Density  jacobian_density(G, t, lam, cotree)
           exponent(y, z)  monomial(y, powers)  theta_D(lam, interior)  check_kirchhoff(G, t, lam, y, C)
           Density: .coeff .radicand .expo .rational() .sign() .square() .scaled(c) .same_as(other) float()
oracle     kernel(G, lam, T=0, powers=None, z=None, method='lattice'|'gauss'|'qmc', h, tol, m, npoints, nrep)->(value, error)
           kernel_lattice(...)->(value, error, info)  taylor_coefficient(G, lam, alpha, T)  density_batch  scaled_density_batch
spred      sp_reduce(G, D=None, plan=None, prefer='small'|'large')->SPReduction
           plan_chain_moves(G, D, prefer)  plan_generic(G)  find_chains(cur, p, q)
           SPReduction: .moves .blocks .bary .coords .ring .texp .resp .orientation .tight_plan .describe()
             .random_point(rng) .complete_point(pt) .tmap(pt) .t_expr() .ratios_of_t(t) .jacobian_t(pt)
             .currents(lam) .currents_at(lam, pt) .displacement() .displacement_expr()
             .density(lam, T=0, method='auto'|'local'|'generic', check=2)->RatioDensity
             .selfcheck(lam, T, npoints, seed, dens)  .integrate(lam, T, powers, insertion, method, check, dens)
             .kernel(lam, T)  .taylor(lam, z, order, T)
           RatioDensity: .num .factors .exps .exponents() .as_expr(factored) .evaluate(vals)->(rational, radicand) .integrate(F)
simplex    dirichlet(exps)->(q, k): q pi^k   gamma_half(x)   PiNumber({k: q})
catalog    enumerate_tight(D, n)  tight_graphs(D, n)  theta_graphs(D)  two_loop_thetas(D)  is_sp  is_series_parallel
           chain_reducible  is_chain_reducible  yd_reducible  skeleton  describe  cone_over_forest  sparse_ok  tight_ok
conventions  sign_n(n)  sign_dy(n, H)  S1_kernel(K, n)  S2_kernel(K, n)  S3_kernel(K, n)  S4_integral(K, G, cotree, H)
           K_from_dy(I, G, cotree, H, T)  raw_factor(source, n)  extraction_sign(G, D, C, seed, c, k)
           extraction_sign_formula(G, D)  varsigma_S3(n, ext)  edge_reversal_sign  edge_permutation_sign  vertex_relabelling_sign
htk        kernel_exact(G, lam, T=0, powers=None, plan=None, check=2)   kernel_numeric(G, lam, T=0, powers=None, z=None)
```

`lam`: list over the non-root vertices (or over all vertices, summing to zero; or dict vertex → vector) of $H$-vectors; `None` for $H=0$. Indices $e$, $v$, $i$ are 0-based; `powers = {(e, i): p}`; `z`: list over edges of $H$-vectors or dict.

Tests: `python3 code/tests/test_known.py` (16 tests), `python3 code/tests/test_engine.py` (11 tests); `HTK_FAST=1` skips the $|I|=7$ oracle runs, the $d=5$ chain move, the enumeration at $(d,n)=(2,6)$, the reduction sweep at $(2,5)$, and restricts the comparison of the two $T=0$ density methods to $|I|\le7$.
