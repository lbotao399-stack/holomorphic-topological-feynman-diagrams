# $Q=\sum\hbar^nQ_n$ on local operators, HT twists at $(H,1)$

Tags: [chk] exact (script, `.log` here); [drv]; [src]; [open]. $D=H+1$; $Q_{\rm odd}=0$ (S4), content $Q_2$, $Q_4$.

## open (newest first)
- $d_2$ beyond the §6 ranges (non-abelian 5d $\Sigma w\ge8$, 7d $\Sigma w\ge6$); limit: exact elimination in sectors of dimension $>4\cdot10^3$.
- 5d $T_{\rm vec}+\sum T_R$ as a quartic invariant (not the symmetric trace: $T_{\rm fund}\not\propto T_{\rm vec}$).
- $Q_4$ at $(2,1)$, $(3,1)$: only 15 cone kernels at random integer sources [src atlas §5]; BF with matter needs non-cone graphs and polynomial $\lambda$-dependence.
- $Q_4$ at $(1,1)$ with derivative letters.
- $N_2$ at $H\ge2$ (normalization, $\epsilon_2$); $d_2$ ranks and nilpotency independent of it [drv].

## 0. Engine (`hteng.py`, `thetak.py`), all $H$

Letters $\partial^\alpha X$, $\alpha\in\mathbb N^H$, $X\in\{c^a,B_a,\varphi^i,\psi_i,u^\alpha\}$; $|c|=1$, $|B|\equiv D$, $|\varphi|=|u|=0$, $|\psi|\equiv D+1$. $\eta(c^a,B_b)=\delta$, $\eta(\varphi^i,\psi_j)=\delta$, $\eta(u^\alpha,u^\beta)=\Omega^{\alpha\beta}$; $\eta(y,x)=(-1)^{D+|x||y|}\eta(x,y)$ [drv: $K\to(-1)^DK$ under edge reversal, CONVENTIONS §2.2] $\Rightarrow$ self-paired $u$ only at $D$ odd.

$$I=\tfrac12f_{bc}{}^aB_ac^bc^c-c^a\varphi\rho_a\psi-W(\varphi)-\tfrac12c^au\,\omega\rho_au,\qquad\omega=\Omega^{-1}$$
$Q_0^2=0$ only for $s_Bs_M=(-1)^D$, $s_Bs_U=-1$ [chk `test_q0.log`, $H=1,2,3$; other signs fail on 659–3275 monomials]. All $D$ [chk `display.log`]: $Q_0c=\tfrac12fcc$, $Q_0B_a=f_{ab}{}^cc^bB_c-\varphi\rho_a\psi$, $Q_0\varphi=c\rho\varphi$, $Q_0\psi=-c\psi\rho-\partial W$.

Graphs: thetas $\theta(n_0,n_1,n_2)$, $\sum n_P=2D+1$, $n_P+n_Q\ge D+1$ [src msym-odd VERIFY]; fan type ($O$ at $N$, $\deg O=3$, $2H$ outputs); tadpole type ($O$ interior, $\deg O=2$, $2H-1$ outputs, contains a one-loop $(D+1)$-cycle avoiding $O$). Rooted shapes, $|\mathrm{Aut}_O|$ [chk `check_thetak.log`]: $H=1$: (1,2,2)@N 2, @(1,1) 2; $H=2$: (1,3,3)@N 2, @(1,1) 1, (2,2,3)@N 2, @(0,1) 2, @(2,1) 2; $H=3$: (1,4,4)@N 2, @(1,1) 1, @(1,2) 2, (2,3,4)@N 1, @(0,1) 2, @(1,1) 1, @(2,1) 1, @(2,2) 2, (3,3,3)@N 6, @(0,1) 2.

Kernels [src msym-odd `theta_kernels.py`, S3 thm-theta]: $y^P_j=w_P\lambda_N+\sum_{Q,k}(w_P-\delta_{PQ})\sigma^Q_k\lambda^Q_k+\sum_{k<j}\lambda^P_k$, $\sigma^Q_k=\sum_{j>k}x^Q_j$; $g=s_\Gamma\tfrac{D_\Gamma}2\prod_Pw_P^{H+1/2-n_P}$, $D_\Gamma=\det\begin{pmatrix}\Lambda_0&\Lambda_1&0\\0&\Lambda_1&\Lambda_2\end{pmatrix}$; $M_p=\pi^{-1}\int g\prod_{e\ni O}y_e^{p_e}$; $\lambda_O=-\sum\lambda_v$, $\lambda_N=\lambda_{N'}=0$; $s_\Gamma$ fitted to `htk.kernel_exact`. $z=0$: $K/\pi=s\,c_\Gamma D_\Gamma$,
$$c_\Gamma=\frac{\prod_P\Gamma(H+\tfrac32-n_P)}{2\pi\,\Gamma(H+\tfrac32)\prod_P(n_P-1)!}:\ \tfrac13\ (H=1);\ \tfrac1{20},\tfrac1{30}\ (H=2);\ \tfrac1{252},\tfrac1{420},\tfrac1{840}\ (H=3)$$
tadpole example $H=2$, (1,3,3)@(1,1): $K/\pi=\tfrac1{20}\langle\lambda_0\lambda_1\rangle\langle\lambda_O\lambda_4\rangle$. [chk `check_thetak.log`: moments = htk, 17 rooted shapes $H\le3$; Taylor to order 2 = htk `taylor`, $H=1,2$]

Rule: sequential left multiplication over $v=n..0$; contraction $\eta(x,y)\partial^L_y\partial^L_x$, edges $(v,w)$ in order $(-v,w)$; no extra sign ($(-1)^{D(D+1)}=1$) [drv].
$$\hat Q_2O=\sum_{(\Gamma,O)}|\mathrm{Aut}_O\Gamma|^{-1}\sum_{S\subset O,|S|=\deg O}\varepsilon_SF_\Gamma(O_S)O_{S^c},\qquad Q_2=N_2\hat Q_2$$
$Q_4$ at $H=1$, $z=0$: fan(4) $K=\pi^2(\lambda_0\lambda_4)^2/180$, tadpole(5,0) $K=\pi^2\lambda_0^4/630$ [src 3dN2 §3].

Interface: `Theory(H, G, f, bf=[(dim,rho)], W, sym=[(dim,Omega,rho)], omlow)`; `Engine(th).Q0(P)`, `.Q2(P, which)`; `thetak.rooted_shapes(H)`; `lift.Lifter(spec).run_sector(w, letters, out)`.

## 1. Validation
- $H=1$ vs 3dN2 engine (`ref_ht3d.py`), fan and tadpole separately, 0 mismatches [chk `check_h1_vs_3dN2.log`]: U(1) charges $1,2,-3$, $W=\varphi_1\varphi_2\varphi_3$ (3696 monomials); SU(2)+adj (5812); GL(2)+1F (376); 3d N=4 SQED (1352).
- $(Q^2)_{\hbar^4}=0$, $H=1$: $\hat Q_2^2+u\{Q_0,\hat Q_4^{\rm fan}\}+v\{Q_0,\hat T_{5,0}\}=0$: GL(2)+1F $u=v=1$ (296 equations; $t_4=\tfrac1{630}$, $\epsilon_4=1$) [chk `test_q4_gl2f_B4X.log`]; SU(2)+adj $u=1$ (336) [chk]; U(1) $u=1$ [chk].
- $d_2$ control U(1)$_0$ + chiral ($r=\tfrac12$): $\hat Q_2(-\tfrac12B^2\partial c+B\partial\varphi\psi)=\tfrac12\partial c\partial^2c$; $J=4$: $\mapsto-\tfrac12\partial c\partial^3c$; ranks $\{1,1\}$ ($J=3$), $\{2,2\}$ ($J=4$) [chk `lift_3d_u1chiral_16.log`, `corrections_3d_u1chiral.log`] = 3dN2 §4.
- Labelling independence (random search order, sign refitted), 0 mismatches [chk `check_labelling.log`]: 5d SU(2)+fund (1389), 5d half-hyper (710), 7d (56), 3d N=4 SU(2) (5940).
- MSYM: tadpole part $\equiv0$ in 5d N=2 and 7d; fan weights $1/|\mathrm{Aut}_O|$: $(\tfrac12,\tfrac12)$ (5d), $(\tfrac12,1,\tfrac16)$ (7d) = msym-odd null spaces [src].

## 2. 3d N=4, $(1,1)$
Fields [src ESW 11.3]: $c,B$; $\Phi$ ($r=1$); $Q,\tilde Q$ ($r=\tfrac12$); $W=\tilde Q\Phi Q$ [drv]. SQED $N_f=1$, $z=0$, fan [chk `display.log`]:
$$\hat Q_2(B^3)=-2(\partial Q\partial\psi_Q-\partial\tilde Q\partial\psi_{\tilde Q}),\ \hat Q_2(B^2\psi_\Phi)=-2\partial Q\partial\tilde Q,\ \hat Q_2(B^2Q)=\tfrac23\partial c\partial Q,\ \hat Q_2(B^2\psi_Q)=-\tfrac23(\partial c\partial\psi_Q+\partial\Phi\partial\tilde Q),\ \hat Q_2(B\psi_\Phi\psi_Q)=\tfrac23\partial c\partial\tilde Q$$
Tadpole $\propto d^U_{abc}=0$ for $U=T^*V$ [drv; chk]; $\{Q_0,\hat T_{5,0}\}\equiv0$.

| theory | $\{Q_0,Q_2\}=0$ | $(Q^2)_{\hbar^4}$, $z=0$ |
|---|---|---|
| SQED $N_f=1$ | 1352 | $u=1$, 14 equations |
| SQED $N_f=2$ | 5788 ($\deg\le4$, $\partial\le1$) | $u=1$, 32 |
| SU(2) $N_f=1$ | 5940 | $u=1$, 354 |

Lowest classes with $\hat Q_2z\ne0$ (exact $\Rightarrow$ correction, no lift): SQED $N_f=1$, $J=3$: $z=\tfrac12B^2\Phi\psi_\Phi-\tfrac12B^2Q\psi_Q+B\psi_\Phi\psi_Q\psi_{\tilde Q}$, $\hat Q_2z=Q_0y$, $y=-\tfrac13B\partial^2c-\tfrac23\partial B\partial c+\tfrac13\Phi\partial^2\psi_\Phi+\tfrac13\partial\Phi\partial\psi_\Phi+\tfrac23\partial^2\Phi\psi_\Phi-\tfrac13Q\partial^2\psi_Q+\tfrac13\partial Q\partial\psi_Q-\tfrac13\tilde Q\partial^2\psi_{\tilde Q}$ [chk `corrections_3d_sqed1.log`]; SU(2) $N_f=1$, $J=2$: $\hat Q_2z=\tfrac{16}3(c^e\partial c^f\partial c^h-\partial c^ec^f\partial c^h)=Q_0(\tfrac83\partial c^e\partial c^f)$ [chk].

## 3. 3d N=8, $(1,1)$
Fields [src ESW 11.8]: $c,B$; $\Phi,X,Y$ adjoint; $W=\kappa(\Phi,[X,Y])$. $\hat Q_2(B_0^2c^0)=-2\partial c^1\partial c^2$, $\hat Q_2(B_0B_1c^0)=\tfrac23\partial c^0\partial c^2$, $\hat Q_2(B_1^2c^1)=2\partial c^0\partial c^2$ (= SU(2)+adj; matter arms cancel on $BBc$) [chk]. $\{Q_0,Q_2\}=0$ on 10088 [chk `test_nilp_3d_su2N8_3_1.log`]; $(Q^2)_{\hbar^4}$: $u=1$, 314 equations [chk]. $J=2$: $\hat Q_2z=\tfrac83(c^e\partial c^f\partial c^h-\partial c^ec^f\partial c^h)=Q_0(\tfrac43\partial c^e\partial c^f)$ [chk].

## 4. 5d, $(2,1)$
ESW 9.2 ($U/\!/\mathfrak g$) = 9.3 (BF $R/\mathfrak g$) for $U=T^*R$ under $u^i\mapsto\varphi^i$, $u^{\dim R+i}\mapsto\psi_i$: 0 mismatches (SU(2) fund 1419; U(1) $q=2$ 128); $u^{\dim R+i}\mapsto-\psi_i$ fails (862) [chk `check_forms.log`].

$X^a=\partial_1c^a\partial_2c^a$, $e_2(X)=X^0X^1+X^0X^2+X^1X^2$: $\hat Q_2^{\rm fan}(B_aB_bc^c)=\tfrac{11}{10}\varepsilon_{abc}e_2(X)$; N=2: $\hat Q_2(B_0\varphi^1\psi_2)=\tfrac{11}{10}e_2(X)$; pure N=1 tadpole: $\hat Q_2^{\rm tad}(B_aB_b)=\tfrac{11}{10}\varepsilon_{abd}\sum_eT(c^e;c^d)$, $T(x;y)=\partial_2x\partial_{12}x\partial_1y-\partial_2x\partial_1^2x\partial_2y+\partial_2^2x\partial_1x\partial_1y+\partial_1x\partial_{12}x\partial_2y$ [chk].

Tadpole part $T=T_{\rm vec}+\sum T_R$ (402 coefficients) [chk `tadpole_trace.log`]: $T_{T^*\mathfrak g}=-T_{\rm vec}$ (5d N=2 cancels); $T_{T^*R}=2T_{\frac12R}$; $T_{\rm abelian}=0$; $T_{T^*{\rm fund}}/T_{\rm vec}$ not constant (23 ratios in $[\tfrac{43}{48},1]$).

| theory | $\{Q_0,Q_2\}=0$ | fan / tadpole alone | control |
|---|---|---|---|
| N=1 pure | 1127 ($\partial\le2$) | 438 / 438 failures | all 5 shapes needed |
| N=1 + fund ($R/\mathfrak g$) | 5337 | 822 / 822 | all 5 |
| N=1 + fund ($U/\!/\mathfrak g$) | 1389 | 162 / 162 | all 5 |
| N=1 + half-hyper | 710 | 126 / 126 | all 5 |
| N=1 U(1), $q=1,2$ | 976 | tadpole $\equiv0$ | — |
| N=2 | 9148 ($\partial\le2$) | tadpole $\equiv0$, fan closed | both fan shapes |

Lowest classes with $\hat Q_2z\ne0$ (all exact): pure N=1 $z=B_eB_fB_h$ ($\Sigma w=6$), $y=\tfrac{22}5(B_e\partial_2^2c^e\partial_1^2c^h+\dots)$; half-hyper $z=(-2B_eu^0u^0+B_hu^0u^1)c^ec^fc^h$, $\hat Q_2z=\tfrac{27}{10}c^e\partial_2c^e\partial_1c^e c^f\partial_2c^f\partial_1c^f c^h+\dots$ [chk `corrections_5d_*.log`].

## 5. 7d N=1, $(3,1)$
Fields [src ESW 7.2]: $c$ odd, $B$ even ($=$ CS on $\mathbb C^{3|1}\times\mathbb R$, $C=c+\theta B$ [drv]). $Y^a=\partial_3c^a\partial_2c^a\partial_1c^a$: $\hat Q_2(B_0^2c^0)\ni\tfrac{88}{35}Y^1Y^2$; coefficients $\{\tfrac23,\tfrac29,\tfrac{124}{315},\tfrac{88}{35},\tfrac3{35}\}$ [chk]. $\{Q_0,Q_2\}=0$ on 488 ($\partial\le1$); tadpole $\equiv0$; (1,4,4)@N, (2,3,4)@N, (3,3,3)@N each needed [chk `test_nilp_7d_so3_3_1.log`]. $z=B_eB_fc^ec^fc^h$, $w=(2,2,2)$: exact [chk].

## 6. $d_2=[\hat Q_2]$ on $H(Q_0)$
$d_2$: af $\to$ af $-2$, af $=N_B+N_\psi$ [drv]. [chk `lift_*.log`, `lift_summary.log`]

| theory | range | $\sum\dim H(Q_0)$ | cocycles with $\hat Q_2z\ne0$ | $\mathrm{rk}\,d_2$ |
|---|---|---|---|---|
| 3d SQED $N_f=1$ | $J\le6$ | 2126 | 2487 | 0 |
| 3d SQED $N_f=2$ | $J\le4$ | 2108 | 962 | 0 |
| 3d N=4 SU(2) $N_f=1$ | $J\le13/4$ | 190 | 4005 | 0 |
| 3d N=4 SU(2) pure | $J\le9/2$ | 388 | 22470 | 0 |
| 3d N=8 SU(2) | $J\le9/4$ | 690 | 2929 | 0 |
| 5d U(1) + hyper | $\Sigma w\le8$ | 958 | 10 | 0 |
| 5d SU(2) pure | $\Sigma w\le7$, (4,4), (5,3) | 168 | 2152 | 0 |
| 5d SU(2) + fund | $\Sigma w\le13/2$ | 118 | 287 | 0 |
| 5d SU(2) + half-hyper | $\Sigma w\le7$ | 102 | 3247 | 0 |
| 5d N=2 SU(2) | $\Sigma w\le13/2$ | 766 | 1703 | 0 |
| 7d SU(2) | $\Sigma w\le5$ | 102 | 5 | 0 |

No two-loop lifts in these ranges: $\hat Q_2\ne0$ on cochains, every image $Q_0$-exact.

## 7. Literature
- CDG §1.4.1: chiral formality of the bulk algebra unknown [src `/Users/libotao/Documents/Codex/2026-10-01/xi/outputs/costello-twists-corpus/text/2005.00083.txt` l.450–453]; 3d N=2 U(1)$_0$ + chiral is a non-formal instance [src 3dN2 §4]; no $\mathcal N\ge4$ instance in the §6 ranges.
- CDG §1.4.2: 3d N=4 HT twist deforms to A/B [src l.479–488]; A perturbatively trivial, B = CS $U/\!/\mathfrak g$ on $\mathbb R^3$, $(0,3)$ $\Rightarrow Q_{n\ge1}=0$ (S3); whether $d_2=0$ follows [open].
- ESW §9.1: 5d tadpole type = theta containing the 6d box; cancels for N=2 [chk §4; analogy drv].

## Files
Engine: `hteng.py`, `thetak.py`, `theories.py`, `specs.py`, `lift.py`, `corrections.py`. Checks: `check_thetak.py`, `check_h1_vs_3dN2.py`, `check_labelling.py`, `check_forms.py`, `test_q0.py`, `test_nilp.py`, `test_q4.py`, `tadpole_trace.py`, `display.py`, `summarize_lifts.py`. References (copied): `ref_ht3d.py`, `ref_fanengine.py`, `ref_theories.py`, `kernels_4.pkl`.

## 8. Replay (main session, scratch copy)
`check_h1_vs_3dN2.py`, `check_forms.py`, `test_q4.py gl2f 5 40 B4X`, `test_nilp.py 3d_su2N8 3 1`, `test_nilp.py 5d_N2_so3 3 1`, `check_labelling.py`: logs identical to the agent's up to timings (0 diff lines).
