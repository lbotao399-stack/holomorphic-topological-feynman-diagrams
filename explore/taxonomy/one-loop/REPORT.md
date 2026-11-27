# One-loop obstructions of all holomorphic twists; vanishing at $T\ge2$ and odd loops at $T=1$

Tags: [chk] exact (script, log here); [drv]; [src]; [open]. ESW = arXiv:2002.10517 Part II, local text `~/Desktop/研究项目_Git/qcohomology_project/09_logs/taxonomy_notion_import/generated/` ("file l.line"). CONV = `code/htk/CONVENTIONS.md`.

## open (newest first)
- C4 UV step for $\kappa\ne0$ (§5.1 iv): weights $c^{j/2}$, $j\ge1$, exact [chk]; integrability at subgraph faces of $\mathbb{RP}^{N-1}_{>0}$ argued, not proven.
- Overall constant $N^{(H)}$ of $\omega_1$ fixed only at $H=3$ (project convention) and $H=5$ (Gwilliam–Williams convention).
- Twisted 6d Green–Schwarz exchange (tensor multiplet) not computed; only the factorization condition.
- Components of $O_1$ with matter or $B$ legs not computed; only $H+1$ gauge legs.
- Fermi multiplet not in ESW; its twist ($bc$ system, $F$ odd) [drv].
- $\rho_n\wedge\rho_n=0$ on graphs exhaustive only for skeletons with $n\le4$; $n=6$ (6 random cubic skeletons), $n=8$ (1) sampled.

## 1. Polygon kernel at $(H,0)$

$m=H+1$, $e_j=(j\to j+1\bmod m)$, root $H$, sources $\lambda_0..\lambda_{H-1}$, $\Delta=\det(\lambda_0..\lambda_{H-1})$ (`htk.polygon(m)`):
- KCL $y_{e_v}-y_{e_{v-1}}=\lambda_v\Rightarrow y_e=\ell+Q_e$, $Q_e=\sum_{u\le e}\lambda_u$ ($e<H$), $Q_H=0$; $C=(1,\dots,1)$, $\sigma_T=+1$.
- KVL $\sum_et_ey_e=0\Rightarrow\ell=-U^{-1}\sum_{e<H}t_eQ_e$.
- $[\Omega]_{dt_0\cdots dt_{H-1}}=U^{-H}\det(y_0..y_{H-1})=U^{-H}(t_H/U)\Delta$.
- $D=(-1)^Ht_H^{-1}[\Omega]_{\widehat{dt_H}}=(-1)^H\Delta U^{-H-1}$, $K=\int_{\Sigma t=1}D=(-1)^H\Delta/H!$.

$$c_H=\frac{(-1)^H}{H!}:\ -1,\ \tfrac12,\ -\tfrac16,\ \tfrac1{24},\ -\tfrac1{120}\quad(H=1..5)$$
[drv; chk polygon.log: $H\le8$ every omitted edge; `htk.kernel_exact` $H\le6$]

$\sum p_j=0\Rightarrow\mathrm{sgn}(s)\det(p_{s(0)},..,p_{s(H-1)})=\det(p_0..p_{H-1})$ $\forall s\in S_{H+1}$ [drv; chk $H\le5$].

Cross-checks [chk polygon.log (d)]: bigon $c_1=-1$ (banana(2) $+1$, one reversed edge); LLC triangle fan(1) $-\Delta/2$, $I^{S1}=+\Delta/2$ = LLC eq-1-42; Poisson-BF $|c_2|=\tfrac12$; box $-\tfrac16$ (anomalies); msym-D8 pentagon $+\tfrac1{24}$; type-i hexagon $-\tfrac1{120}$.

## 2. Polygon $\to O_1$

$$\omega_1:=[\theta_0\cdots\theta_H]\,O_1[A=\textstyle\sum\theta_jx_je^{p_j\cdot z}]=N^{(H)}\sum_\Gamma|\mathrm{Aut}\,\Gamma|^{-1}\sum_{s\in S_{H+1}}\mathrm{sgn}(s)K_\Gamma(p_s)W_\Gamma(x_s)$$
directed loop $|\mathrm{Aut}|=H+1$, $W=\mathrm{Str}_V(x_{s0}\cdots x_{sH})$; self-paired $|\mathrm{Aut}|=2(H+1)$.

$$\omega_1=(-1)^HN^{(H)}\Delta\,\mathcal S_H(x),\qquad \mathcal S_H=\sum_{\rm pairs}\mathrm{Str}_V\rho_V(x)^{H+1}+\tfrac12\sum_{\rm self}\mathrm{Str}_V\rho_V(x)^{H+1}$$
- pair: $\varphi\in\Omega^{0,*}\otimes V$, $\varphi^\vee\in\Omega^{H,*}\otimes V^\vee[H-1]$, vertex $\langle\varphi^\vee,\rho(A)\varphi\rangle$; self-paired: $\tfrac12\omega(u,\rho(A)u)$ [drv].
- one loop = $\log\mathrm{Ber}$ (pair), $\tfrac12\log\mathrm{Ber}$ (self) of $\bar\partial+\rho(A)$ $\Rightarrow$ colour factor $\mathrm{Str}_V$ [drv].
- $\mathrm{Str}_{V^\vee[H-1]}((-\rho^T)^{H+1})=\mathrm{Str}_V\rho^{H+1}$ [drv; chk supertrace.log "sides"].
- $T[1]Y$: $(1-1)\mathrm{Str}_L=0$; $\mathbb C^H_{\rm Dol}$: $(1-1)^H\mathrm{Str}_L=0$; $\mathbb C^{3|2}$: $\mathrm{sdim}\,\Lambda[\theta_1,\theta_2]=0$ [drv].

Normalization [chk assemble.log]: $H=3$: $\mathcal S_3=-\tfrac12\mathrm{Str}_{\rm ad}$, $N^{(3)}=1/(8\pi^3)$ [src physics/anomalies §1] $\Rightarrow\omega_1=\Delta\,\mathrm{Str}_{\rm ad}x^4/(16\pi^3)=3\Delta T_1/\pi^3$. $H=5$: $N^{(5)}c_5=1/(120(2\pi i)^5)\Rightarrow\omega_1=-\Delta\,\mathrm{Tr}_{\rm ad}x^6/(2(2\pi i)^5)$ [src type-i "The hexagon"]. Gwilliam–Williams: $P_{H+1}(X)=\mathcal S_H(X)/((H+1)!(2\pi i)^H)$ [src GW 1810.06534 Prop. 3.9].

Physical side: $I_{d+2}\supset\sum_{\Pi S_\chi\otimes E}\tfrac12\chi\,\mathrm{Tr}_EF^{H+1}$, $\chi(S_\pm)=\pm1$.

$$\mathcal S_H=-I_{\rm phys}\quad\text{in every row}\quad[\text{chk supertrace.log, 38/38; sl(3), sl(2), random rational }x]$$

## 3. Per twist

$S_{BF}=\int\langle B,\bar\partial A+\tfrac12[A,A]\rangle+\sum_{\rm pairs}\int\langle\beta,\bar\partial\gamma+\rho(A)\gamma\rangle+\sum_{\rm self}\tfrac12\int\omega(u,\bar\partial u+\rho(A)u)$, $S_{hCS}=\int\Omega\wedge(\tfrac12\langle A,\bar\partial A\rangle+\tfrac16\langle A,[A,A]\rangle)$. $\mathrm{Tr}_X:=\mathrm{Tr}_X\rho(x)^{H+1}$.

| theory | ESW | $V$ (parity) | $c_H$ | $\mathcal S_H$ | physical gauge term | agree |
|---|---|---|---|---|---|---|
| 2d (2,2), $R$ | 12.1 | $\mathfrak g[1](-)$, $\mathfrak g[2](+)$, $R(+)$, $R[1](-)$ | $-1$ | $0$ | $0$ | yes |
| 2d (0,2), $R$, Fermi $F$ | 12.16 | $\mathfrak g[1](-)$, $R(+)$, $F(-)$ | $-1$ | $-\mathrm{Tr}_{\rm ad}+\mathrm{Tr}_R-\mathrm{Tr}_F$ | $\mathrm{Tr}_{\rm ad}-\mathrm{Tr}_R+\mathrm{Tr}_F$ | yes |
| 2d (4,0), $U$ | 12.17 | $\mathfrak g[1](-)$, $\mathfrak g^\vee[-1](-)$, $U(+)$ | $-1$ | $-2\mathrm{Tr}_{\rm ad}+\mathrm{Tr}_U$ | $2\mathrm{Tr}_{\rm ad}-\mathrm{Tr}_U$ | yes |
| 2d (4,4), $U$ | 12.7 | (4,0) $\otimes\mathbb C[\varepsilon]$ | $-1$ | $0$ | $0$ | yes |
| 2d (8,8) | §12 | $(\mathfrak g[1]\oplus\mathfrak g)\otimes\mathbb C[\varepsilon]\otimes\Lambda[dz]$ | $-1$ | $0$ | $0$ | yes |
| 2d $(\mathcal N_+,0)$ | 12.18 | $\mathfrak g[1](-)$; self $\mathfrak g\otimes\mathbb C^{\mathcal N_+-2}(-)$ | $-1$ | $-\tfrac{\mathcal N_+}2\mathrm{Tr}_{\rm ad}$ | $\tfrac{\mathcal N_+}2\mathrm{Tr}_{\rm ad}$ | yes ($\mathcal N_+=2..7$) |
| 4d N=1, $R$ | 10.1 | $\mathfrak g[1](-)$, $R(+)$ | $\tfrac12$ | $\mathrm{Tr}_R$ | $-\mathrm{Tr}_R$ | yes |
| 4d N=2, $U$ | 10.2 | $\mathfrak g[1](-)$, $\mathfrak g[2-4n](+)$, $U[2n](+)$ | $\tfrac12$ | $0$ | $0$ | yes |
| 4d N=4 | 10.9 | $\mathfrak g[1]\otimes\Lambda[dz_1,dz_2]$ | $\tfrac12$ | $0$ | $0$ | yes |
| 6d (1,0), $U/\!/\mathfrak g$ | 8.1 | $\mathfrak g[1](-)$; self $U\otimes K^{1/2}(+)$ | $-\tfrac16$ | $-\mathrm{Tr}_{\rm ad}+\tfrac12\mathrm{Tr}_U$ | $\mathrm{Tr}_{\rm ad}-\tfrac12\mathrm{Tr}_U$ | yes |
| 6d (1,0), BF $R/\mathfrak g$ | 8.2 | $\mathfrak g[1](-)$, $R(+)$ | $-\tfrac16$ | $-\mathrm{Tr}_{\rm ad}+\mathrm{Tr}_R$ | $\mathrm{Tr}_{\rm ad}-\sum_{\rm hyp}\mathrm{Tr}_R$ | yes |
| 6d (1,1) | 8.3 | $\mathfrak g[1](-)$, $\mathfrak g(+)$ | $-\tfrac16$ | $0$ | $0$ | yes |
| 8d N=1 | 6.1 | $\mathfrak g[1](-)$ | $\tfrac1{24}$ | $-\mathrm{Tr}_{\rm ad}x^5\equiv0$ | $0$ | yes |
| 10d (1,0), hCS $\mathbb C^5$ | 4.4 | self $\mathfrak g[1](-)$ | $-\tfrac1{120}$ | $-\tfrac12\mathrm{Tr}_{\rm ad}$ | $\tfrac12\mathrm{Tr}_{\rm ad}$ | yes |

[chk supertrace.log, formal column, both test algebras]. Consistency [chk]: 6d $U/\!/\mathfrak g$ with $U=T^*R$ = BF $R/\mathfrak g$ (both $851929/50$ at the same $x$); formal (4,0) at $U=0$ = $(\mathcal N_+=4)$: $-2\mathrm{Tr}_{\rm ad}$.

## 4. Class, counterterm, Green–Schwarz

- $Q_0$ preserves $B$-degree and matter degree $\Rightarrow$ pure-gauge $dI$ from pure-gauge $I$ only [drv].
- Translation- and $SL(H)$-invariant $H$-input cochains of momentum degree $H$: polynomials in $H-1$ momenta of $\mathbb C^H$, none by FFT $\Rightarrow B^1=0$ [drv; $H=3$ src physics/anomalies §4] $\Rightarrow$ $[O_1]=0\iff\mathcal S_H\equiv0$.
- GS: tree exchange $\Rightarrow\sum_aP_aQ_a$; cancellation needs $\mathcal S_H\in\mathrm{span}\{P_aQ_a\}$ [drv].

Traces [chk casimirs.log; sl(2..5), so(3..8), sp(2,4,6), $k=2..6$], $p_j=\mathrm{tr}_{\rm fund}x^j$:
- $\mathfrak{sl}_N$: $\mathrm{Tr}_{\rm ad}x^2=2Np_2$; $x^4$: $2Np_4+6p_2^2$; $x^6$: $2Np_6+30p_2p_4-20p_3^2$.
- $\mathfrak{so}_N$: $(N-2)p_2$; $(N-8)p_4+3p_2^2$; $(N-32)p_6+15p_2p_4$.
- $\mathfrak{sp}_N$: $(N+2)p_2$; $(N+8)p_4+3p_2^2$; $(N+32)p_6+15p_2p_4$.
- $\mathfrak{sl}_N$: $\mathrm{tr}_{\Lambda^2}x^3=(N-4)p_3$, $\mathrm{tr}_{S^2}x^3=(N+4)p_3$; $\mathrm{tr}_{\Lambda^2}x^4=(N-8)p_4+3p_2^2$, $\mathrm{tr}_{S^2}x^4=(N+8)p_4+3p_2^2$.

| $H$ | $\mathcal S_H$ | $[O_1]$ | GS |
|---|---|---|---|
| 1 | $\mathfrak{sl}_N$, $n_f(N\oplus\bar N)$, $n_F$ Fermis: $-(2N+n_F-2n_f)p_2$; $(\mathcal N_+,0)$: $-\mathcal N_+Np_2$ | $\ne0$ unless quadratic indices cancel | none (semisimple) |
| 2 | $\mathrm{Tr}_Rx^3$, $A(\Lambda^2)=N-4$, $A(S^2)=N+4$, $A(\bar N)=-1$ | $\ne0$ unless $\sum A(R)=0$ | none |
| 3 | $\mathfrak{sl}_N$, $n_fN+n_A\Lambda^2+n_SS^2$: $(n_f+(N-8)n_A+(N+8)n_S-2N)p_4+3(n_A+n_S-2)p_2^2$; $\mathfrak{so}_N$, $n_v$: $(n_v-N+8)p_4-3p_2^2$; $\mathfrak{sp}_N$, $n_f$: $(n_f-N-8)p_4-3p_2^2$ | $\ne0$ unless both vanish | iff $p_4$ coefficient $=0$; $p_2^2$ by tensor exchange [open] |
| 4 | $\equiv0$ | $0$ | — |
| 5 | $\mathfrak{so}_N$: $-\tfrac12(N-32)p_6-\tfrac{15}2p_2p_4$ | $\ne0$ $\forall N$ | $N=32$: $P_6=-p_2p_4/(96(2\pi i)^5)$ cancelled by type I KS exchange $+1/(96(2\pi i)^5)$; residual $-(N-32)p_6/(1440(2\pi i)^5)$ [chk; = type-i "closed exchange"; Costello–Li 1905.09269 (82)–(85)] |

## 5. C4

### 5.1 Common
(i) Tight graph: $\deg\prod\Theta_e=D|I|=D(|V|-1)+(N-1)+\ell\Rightarrow\ell=0$; each topological direction $a$: $\mathcal T_\Gamma[\kappa^a]=U^{-1/2}e^P\langle\prod_{r=1}^n\psi^a_r\rangle$, $\psi_r=\sum_eC_{re}w_edt_e$ [drv].

(ii) $\langle w_e\rangle=\bar w_e$, $\langle\delta w_e\delta w_f\rangle=\tfrac12(\delta_{ef}/t_e-K_{ef})$, $K=C^TM^{-1}C$ [src CONV §2]; $G_S=\det C_{[:,S]}\,E[\prod_{f\in S}w_f]$; $\kappa=0$: $\langle\prod\psi\rangle=\mathrm{Pf}(A)$, $A=-\tfrac12dM\,M^{-1}dM$, $\mathcal T=\rho_n$ [chk c4_s4.log "pf_route"].

(iii) $\varphi_c:t\to ct$: $\varphi_c^*\mathcal T_\Gamma[\kappa]=\mathcal T_\Gamma[\sqrt c\,\kappa]$ $\Rightarrow$ $\mathcal T=\sum_j\mathcal T^{(j)}$, $\varphi_c^*\mathcal T^{(j)}=c^{j/2}\mathcal T^{(j)}$, $\mathcal T^{(0)}=\rho_n$, $j\equiv n\ (2)$ [drv; chk c4_s4.log].

(iv) Weight $s>0$ contributes $\varepsilon^s\to0$ on the UV face; only $\rho_n$ enters [drv; integrability open].

### 5.2 (S3) $T\ge2$, $n\ge1$
$$\rho_n\wedge\rho_n=U^{-1}\mathrm{Pf}(A)^2=U^{-1}\det A=(-\tfrac12)^n\det(X\wedge X),\quad X=M^{-1}dM$$
$Y=X\wedge X$: $\mathrm{tr}\,Y^j=\mathrm{tr}\,X^{2j}=-\mathrm{tr}\,X^{2j}=0$ $\Rightarrow$ (Newton) $e_k(Y)=0$ $\Rightarrow\det Y=0$ $\Rightarrow\Omega_\Gamma=U^{-H}\rho_n^{\wedge T}\wedge\prod\psi^i_r=0$; $\kappa\ne0$: nonzero multidegree has weight $>0$ [drv]. Every bulk and $Q_n$ kernel with $n\ge1$ vanishes at $T\ge2$, all $H$, $n$. Uses only $M(t)$ symmetric positive; on graphs $\rho_n(\Gamma)$ = pullback from the skeleton.

[chk c4_s3.log] generic $M$ mod $2^{61}-1$, $10^9+7$: $n=2,4,6$ ($\mathrm{Pf}\ne0$, $\mathrm{Pf}\wedge\mathrm{Pf}=0$, $\mathrm{tr}(X\wedge X)^j=0$), $n=8$; all skeletons $n=2,3,4$ (1, 4, 17), $\rho_4\ne0$; random cubic skeletons $n=6$, $n=8$; admissible graphs $D=5,6$, $n=2,3$; thetas; $n=4$ subdivisions $D=4,5,6$; with sources two directions $n=2$, $D=2..5$.

### 5.3 (S4) $T=1$, odd $n$, all $H$
$\kappa=0$: $E[\prod_{f\in S}\delta w_f]=0$ for $|S|$ odd $\Rightarrow\rho_n=0\Rightarrow\Omega_\Gamma=0$; $H$ enters only via $U^{-H}\prod\psi^i$; $\kappa\ne0$: $j$ odd $\Rightarrow$ weight $\ge\tfrac12$ [drv]. $\Rightarrow K_\Gamma=0$, $Q_{\rm odd}=O_{\rm odd}=0$ at $T=1$. Reflection $\tau\to-\tau$: $K=(-1)^nK$ [src S3 cor-odd].

[chk c4_s4.log] $n=1$ polygons $H=1..6$; $n=3$ all admissible $D=2..6$ (3, 16, 42, 116, 242); $n=5$ all 70 at $D=2$; $\kappa=0\Rightarrow G=0$, $\kappa\ne0\Rightarrow G\ne0$ with $\kappa$-degrees $\{1,3\}$, $\{1,3,5\}$; even controls (thetas $D=2..7$; $n=4$ at $D=2$ (13), $D=3$ (230)).

## 6. Files
`common.py`, `lie.py`, `c4graphs.py`; `polygon.py`, `supertrace.py`, `casimirs.py`, `assemble.py`, `c4_s3.py`, `c4_s4.py` (each $\to$ `.log`, < 2 min; `-W ignore` for `c4_*`).

## 7. Replay (main session)
`polygon`, `supertrace`, `casimirs`, `assemble`: reran, logs identical (0 diff lines); `supertrace.log` 40 True, 0 False.
`c4_s3` (57 s), `c4_s4` (96 s): reran, logs identical, 0 False, "all checks: True".
